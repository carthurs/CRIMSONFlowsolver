from math import pi, cos
import numpy
import scipy.interpolate

class parameterController:

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile):

		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_nameOfThisScript = baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile
		self.getPeriodicPressurePrescriberData()

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByComponentIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.updatePeriodicTime(delt)	
		prescribedFlow = self.pressureFunction(self.m_periodicTime)

		return prescribedFlow.astype(float)

	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,self.endTime):
		if self.m_periodicTime >= self.endTime:
			self.m_periodicTime = self.m_periodicTime - self.endTime

	def getPeriodicPressurePrescriberData(self):
		# Load a flow file which has the same name as this controller script, but with extension ".dat".
		# It should contain two columns in plain text: the first gives the time,
		# and the second gives the flow to prescribe at that time.
		#
		# The time should start at zero.
		#
		# This script will loop the flow data when it reaches the end of the time.
		self.pressureFileData = numpy.loadtxt(self.m_nameOfThisScript+'.dat')
		self.pressureFunction = scipy.interpolate.interp1d(self.pressureFileData[:,0],self.pressureFileData[:,1])
		self.endTime = self.pressureFileData[-1,0].astype(float)