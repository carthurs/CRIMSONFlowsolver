from CRIMSONPython import *
from math import pi, cos

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.m_baseNameOfThisScript = baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_timeToMaximumElastance = 0.2782;
		self.m_timeToRelax = 0.1391;
		self.m_minimumElastance = 4.10246e-3;
		self.m_maximumElastance = 3.0827e-1;
		self.m_heartPeriod = 0.86;
		self.finishSetup()

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.updatePeriodicTime(delt)
		elastance = self.getElastance(currentParameterValue)

		self.clearBroadcastData()
		self.addBroadcastVariable('foo', 1234.5)
		self.addBroadcastVariable('bar',55646)
		self.addBroadcastVariable('beans','heinz')
		self.addBroadcastVariable('cutlery','useful')

		# for key in dictionaryOfPressuresByNodeIndex:
		# 	print "Pressure ", key, " was ", dictionaryOfPressuresByNodeIndex[key]
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