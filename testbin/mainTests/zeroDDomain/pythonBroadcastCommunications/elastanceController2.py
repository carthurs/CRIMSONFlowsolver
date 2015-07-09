from CRIMSONPython import *
from math import pi, cos

# The parameter controller must have exactly this name
class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile):
		# import io
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile)
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_timeToMaximumElastance = 0.2782;
		self.m_timeToRelax = 0.1391;
		self.m_minimumElastance = 4.10246e-3;
		self.m_maximumElastance = 3.0827e-1;
		self.m_heartPeriod = 0.86;
		self.finishSetup()

	# This method must have exactly this name
	def setFirstTimestepBroadcastValues(self):
		self.clearBroadcastData()
		self.addBroadcastVariable('five', 5)
		self.addBroadcastVariable('six', 6)
		self.addBroadcastVariable('LeftVentricularVolume', 100000.0)
		self.addBroadcastVariable('AorticValveFlow', 0.0)
		
	# This method must have exactly this name
	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByComponentIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.clearBroadcastData()
		self.addBroadcastVariable('five', 5)
		self.addBroadcastVariable('six', 6)
		LeftVentricularVolume = dictionaryOfVolumesByComponentIndex[5]
		self.addBroadcastVariable('LeftVentricularVolume', LeftVentricularVolume)
		AorticValveFlow = dictionaryOfFlowsByComponentIndex[1]
		self.addBroadcastVariable('AorticValveFlow', AorticValveFlow)
		# print "Elastance Controller Reporting!"
		# self.printAllRecievedData()

		controlSignal = self.getRecievedBroadcastValue('masterController','masterControlSignal')
		# print "in elastance:", controlSignal

		# Only update the time if this controller is receiving the (otherwise-unused)
		# broadcasts from other controllers. The only purpose of this is to
		# make the test fail (due to the time not being updated properly) if
		# there is a problem with the broadcast reception.
		if self.getRecievedBroadcastValue('nodeController2','four') == 4 and self.getRecievedBroadcastValue('nodeController_downstream','seven') == 7:
			self.updatePeriodicTime(delt)

		elastance = self.getElastance(currentParameterValue) * (abs(controlSignal) + 0.5)

		# for key in dictionaryOfPressuresByComponentIndex:
		# 	print "Pressure ", key, " was ", dictionaryOfPressuresByComponentIndex[key]
		# for key in dictionaryOfFlowsByComponentIndex:
		# 	print "Flow ", key, " was ", dictionaryOfFlowsByComponentIndex[key]
		# for key in dictionaryOfVolumesByComponentIndex:
		# 	print "Volume", key, "was", dictionaryOfVolumesByComponentIndex[key]

		return elastance


	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod

	def getElastance(self, currentParameterValue):
		# *** analytical elastance function from:
		#     pope, s. r.; ellwein, l. m.; zapata, c. l.; novak, v.; kelley, c. t. & olufsen, m. s.  
		#     estimation and identification of parameters in a lumped cerebrovascular model.
		#     math biosci eng, 2009, 6, 93-115

		# This is the elastance function. It's defined piecewise:
		if self.m_periodicTime <= self.m_timeToMaximumElastance:
			elastance = self.m_minimumElastance \
	        + 0.5*(self.m_maximumElastance - self.m_minimumElastance) \
	        * (1.0 - cos((self.m_periodicTime*pi)/self.m_timeToMaximumElastance))

		elif self.m_periodicTime <= (self.m_timeToMaximumElastance + self.m_timeToRelax):
		 	elastance = self.m_minimumElastance \
		    + 0.5*(self.m_maximumElastance-self.m_minimumElastance) \
		    * (1.0 + cos((self.m_periodicTime-self.m_timeToMaximumElastance)*(pi/self.m_timeToRelax)))

		elif self.m_periodicTime > (self.m_timeToMaximumElastance + self.m_timeToRelax):
			elastance = self.m_minimumElastance

		return elastance;
