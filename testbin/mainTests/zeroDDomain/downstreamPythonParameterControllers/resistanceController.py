from CRIMSONPython import *
from math import pi, cos

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile)
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_timeToMaximumElastance = 0.2782;
		self.m_timeToRelax = 0.1391;
		self.m_minimumElastance = 4.10246e-3;
		self.m_maximumElastance = 3.0827e-1;
		self.m_heartPeriod = 0.86;


	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByComponentIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.updatePeriodicTime(delt)
		pressure = abs(0.001 * cos(pi*100*self.m_periodicTime))

		# for key in dictionaryOfPressuresByComponentIndex:
		# 	print "Pressure ", key, " was ", dictionaryOfPressuresByComponentIndex[key]
		# for key in dictionaryOfFlowsByComponentIndex:
		# 	print "Flow ", key, " was ", dictionaryOfFlowsByComponentIndex[key]

		return pressure

	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod