from CRIMSONPython import *

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = -1
		self.finishSetup()

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		# format is getReceivedBroadcastValue(nameOfBroadcastingScript, nametagOfBroadcastValue)
		leftVentricularPressureToSet = self.getRecievedBroadcastValue('ventricularPressureObserver','leftVentricularPressure')
	
		return leftVentricularPressureToSet