from CRIMSONPython import *

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = -2
		self.finishSetup()

	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		leftVentricularPressureNodeIndex = 3
		observedLeftVentricularPressure = dictionaryOfPressuresByNodeIndex[leftVentricularPressureNodeIndex]

		# Having observed the pressure in the LV, we broadcast it so that the coronary can catch it.
		# The coronary model control script will look for the name 'leftVentricularPressure' in the broadcast data, coming from
		# this script ('ventricularPressureObserver')
		self.clearBroadcastData()
		self.addBroadcastVariable('leftVentricularPressure', observedLeftVentricularPressure)

		observedLeftVentricularVolume = dictionaryOfVolumesByComponentIndex[5]
		self.addBroadcastVariable('leftVentricularVolume', observedLeftVentricularVolume)

		# This is a little bit of a hack: we just use this script to observe the LV pressure, but it
		# must return a new pressure for the node it's attached to in the netlist_surfaces.dat.
		#
		# In this case, this script is attached to a ground node with pressure zero (node 6).
		# We just return it's original value without chaing it.
		unchangedPressureAtNode6 = dictionaryOfPressuresByNodeIndex[6]
		return unchangedPressureAtNode6