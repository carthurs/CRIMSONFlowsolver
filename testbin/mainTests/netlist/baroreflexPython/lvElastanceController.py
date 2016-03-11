from CRIMSONPython import *
from math import pi, cos

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = 1
		self.periodicTime = 0.0; #\todo think about this for restarts!
		self.proportionOfCycleTakenToReachMaxElastance = 0.4#0.31;
		self.proportionOfCycleTakenToRelax = 0.2#0.16;
		self.minimumElastance = 4.10246e-3;
		self.initialMaxElastance = 3.0827e-1;
		self.maximumElastance = self.initialMaxElastance
		self.heartPeriod = 1.0#0.86;
		self.currentCycleIndex = 0

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):
		
		heartRate = self.getRecievedBroadcastValue('baroreceptor','heartRate')
		self.heartPeriod = heartRate/60

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
		# *** analytical elastance function from:
		#     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
		#     estimation and identification of parameters in a lumped cerebrovascular model.
		#     math biosci eng, 2009, 6, 93-115

		timeToMaxElastance = self.heartPeriod * self.proportionOfCycleTakenToReachMaxElastance
		timeToRelax = self.heartPeriod * self.proportionOfCycleTakenToRelax

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