from CRIMSONPython import *

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.finishSetup()

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		# format is getReceivedBroadcastValue(nameOfBroadcastingScript, nametagOfBroadcastValue)
		R_pResistanceToSet = self.getRecievedBroadcastValue('coronaryRdController','Rp')

		#print "R_p was:", R_pResistanceToSet
		
		return R_pResistanceToSet