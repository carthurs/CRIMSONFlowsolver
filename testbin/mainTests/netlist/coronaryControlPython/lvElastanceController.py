from CRIMSONPython import *
from math import pi, cos

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = -3
		self.periodicTime = 0.0; #\todo think about this for restarts!
		self.proportionOfCycleTakenToReachMaxElastance = 0.31;
		self.proportionOfCycleTakenToRelax = 0.16;
		self.minimumElastance = 4.10246e-3;
		self.maximumElastance = 3.0827e-1;
		self.heartPeriod = 0.86;
		self.currentCycleIndex = 0

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.updatePeriodicTime(delt)
		elastance = self.getElastance(currentParameterValue)

		# Step the heart rate to increase MVO2 so we can test the control response:
		if self.currentCycleIndex < 5:
			self.heartPeriod = 0.86
		else:
			self.heartPeriod = 0.4

		self.clearBroadcastData()
		self.addBroadcastVariable('heartPeriod', self.heartPeriod)

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