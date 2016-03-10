from CRIMSONPython import *
import numpy

# This is an example general control script. Edit it, and save it with a meaningful name
# (e.g. downstreamResistanceController.py) in your simulation directory. You'll need to refer to it in netlist_surfaces.dat
# or netlist_closed_loop_downstream.dat.
#
# Its purpose is to adjust a parameter (flow, pressure or component parameter (resistance,compliance, etc.))

class parameterController(abstractParameterController): #NECESSARY

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank): #NECESSARY
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank) #NECESSARY
		self.controllerPriority = 1
		self.finishSetup() #NECESSARY
		# The baseline resistance will be set the first time the controller steps, 
		# using the value set in netlist_surfaces.xml for the resistor its attached to
		self.initialUnstressedVolumeHasBeenSet = False


	# def setFirstTimestepBroadcastValues(self):
	# 	self.clearBroadcastData()
	# 	self.addBroadcastVariable('three', 3) # just a non-functional example broadcast
	# 	self.addBroadcastVariable('four', 4) # just a non-functional example broadcast

	# This method returns a new value which is set as the parameter value of the associated component (e.g. resistance, compliance, etc.)
	# OR the prescribed pressure at a node
	# OR the prescribed flow through the component
	#
	# This depends on how this controller was described in netlist_surfaces.dat or netlist_closed_loop_downstream.dat
	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex): #NECESSARY
		# Uncomment to see what the recieved pressures, flows and volumes in the associated netlist
		# circuit were:
		#
		# for key in dictionaryOfPressuresByNodeIndex:
		# 	print "Pressure ", key, " was ", dictionaryOfPressuresByNodeIndex[key]
		# for key in dictionaryOfFlowsByComponentIndex:
		# 	print "Flow ", key, " was ", dictionaryOfFlowsByComponentIndex[key]

		# Get a baseline resistance on the first time-step taken
		if (not self.initialUnstressedVolumeHasBeenSet):
			self.baselineUnstressedVolume = currentParameterValue
			self.initialUnstressedVolumeHasBeenSet = True

		controlSignal = self.getRecievedBroadcastValue('baroreceptor','venousUnstressedVolumeControlSignal')
		unstressedVolumeToSet = self.baselineUnstressedVolume * controlSignal
		print "python setting unstresed volume", unstressedVolumeToSet
		
		return unstressedVolumeToSet