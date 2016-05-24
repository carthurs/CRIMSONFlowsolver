from CRIMSONPython import *
from math import pi, cos

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = 1
		self.periodicTime = 0.0; #\todo think about this for restarts!
		# self.proportionOfCycleTakenToReachMaxElastance = 0.4#0.31;
		# self.proportionOfCycleTakenToRelax = 0.2#0.16;
		self.minimumElastance = 8.10246e-3 # 4.10246e-3
		self.initialMaxElastance = 0.157443 #3.0827e-1
		self.maximumElastance = self.initialMaxElastance
		self.heartPeriod = 1.0#0.86;
		self.currentCycleIndex = 0

	def updateControl(self, currentParameterValue, delt, dictionaryOflPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):
		
		self.heartRate = self.getRecievedBroadcastValue('baroreceptor','heartRate')
		self.heartPeriod = 60/self.heartRate

		self.updatePeriodicTime(delt)

		elastanceControlSignal = self.getRecievedBroadcastValue('baroreceptor','maximumElastanceControlSignal')
		self.maximumElastance = elastanceControlSignal * self.initialMaxElastance
		
		elastance = self.getElastance(currentParameterValue)
		return elastance


	def updatePeriodicTime(self, delt):
		self.periodicTime = self.periodicTime + delt
		# Keep periodicTime in the range [0,heartPeriod):
		if self.periodicTime >= self.heartPeriod:
			self.periodicTime = self.periodicTime - self.heartPeriod
			self.currentCycleIndex += 1

	def getElastance(self, currentParameterValue):
		# See Fig. 2 of Fananapazir et al. Contribution of heart rate to QT interval shortening during
		# exercise. Eur Heart J, 1983 4:265-271 (Fananapazir83.pdf)
		#
		# POPULATION: PATIENTS WITH PACEMAKERS!
		# QtInterval = (515.0 - 0.85 * self.heartRate)/1000.0 # surrogate for period during which LV is contracting

		
		# Kligfield et al. QT interval-heart rate relation during exercise in normal men and women: definition 
		# by linear regression analysis. JACC, 1996 28:1547-55
		#
		# Kligfield69.pdf - table 1
		QtInterval = (490.0 - 1.45 * self.heartRate)/1000.0 # surrogate for period during which LV is contracting

		# *** analytical elastance function from:
		#     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
		#     estimation and identification of parameters in a lumped cerebrovascular model.
		#     math biosci eng, 2009, 6, 93-115

		# timeToMaxElastance = self.heartPeriod * self.proportionOfCycleTakenToReachMaxElastance
		# timeToRelax = self.heartPeriod * self.proportionOfCycleTakenToRelax

		timeToMaxElastance = 0.8 * QtInterval
		timeToRelax = 0.4 * QtInterval

		#never used this, just an idea: # timeToRelax = QtInterval - timeToMaxElastance

		# This is the elastance function. It's defined piecewise:
		if self.periodicTime <= timeToMaxElastance:
			elastance = self.minimumElastance \
	        + 0.5*(self.maximumElastance - self.minimumElastance) \
	        * (1.0 - cos((self.periodicTime*pi)/timeToMaxElastance))

		elif self.periodicTime <= (timeToMaxElastance + timeToRelax):
		 	elastance = self.minimumElastance \
		    + 0.5*(self.maximumElastance-self.minimumElastance) \
		    * (1.0 + cos((self.periodicTime-timeToMaxElastance)*(pi/timeToRelax)))

		elif self.periodicTime > (timeToMaxElastance + timeToRelax):
			elastance = self.minimumElastance

		return elastance;