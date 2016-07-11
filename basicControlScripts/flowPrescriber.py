from CRIMSONPython import *
from math import pi, cos
import numpy
import scipy.interpolate

# disable serialisation
# def loadClassOnRestart(fileName, MPIRank):
# 	return None

# def saveClassForRestart(objectToIgnore):
# 	pass

class parameterController(abstractParameterController):

	def __getstate__(self):
		odict = self.__dict__.copy()
		del odict['flowFunction']
		return odict

	def __setstate__(self, odict):
		self.__dict__.update(odict)
		self.getPeriodicFlowPrescriberData()


	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_nameOfThisScript = baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile
		self.getPeriodicFlowPrescriberData()
		self.finishSetup()

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.updatePeriodicTime(delt)	
		prescribedFlow = self.flowFunction(self.m_periodicTime)

		return prescribedFlow.astype(float)

	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,self.endTime):
		if self.m_periodicTime >= self.endTime:
			self.m_periodicTime = self.m_periodicTime - self.endTime

	def getPeriodicFlowPrescriberData(self):
		# Load a flow file which has the same name as this controller script, but with extension ".dat".
		# It should contain two columns in plain text: the first gives the time,
		# and the second gives the flow to prescribe at that time.
		#
		# The time should start at zero.
		#
		# This script will loop the flow data when it reaches the end of the time.
		self.flowFileData = numpy.loadtxt(self.m_nameOfThisScript+'.dat')
		self.flowFunction = scipy.interpolate.interp1d(self.flowFileData[:,0],self.flowFileData[:,1])
		self.endTime = self.flowFileData[-1,0].astype(float)