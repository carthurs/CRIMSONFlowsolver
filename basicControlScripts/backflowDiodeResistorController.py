# This script can be used to control a resistor component located where we have a diode/valve
# which we wish to allow backflow through, as opposed to zero flow, when the pressure
# gradient across it is negative.

from CRIMSONPython import *
import sys

class parameterController(abstractParameterController): #NECESSARY

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank): #NECESSARY
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank) #NECESSARY
		self.finishSetup() #NECESSARY

		# SET THSE YOURSELF!:
		self.indexOfStartNodeOfControlledResistor = 
		self.indexOfEndNodeOfControlledResistor = 
		self.resistanceToForwardFlow = 0.0001
		self.maxResistanceWithClosedValve = 0.1#sys.float_info.max

	# def setFirstTimestepBroadcastValues(self): #NECESSARY - note that the masterController.py (i.e. this example controller)
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

		startNodePressure = dictionaryOfPressuresByNodeIndex[self.indexOfStartNodeOfControlledResistor]
		endNodePressure = dictionaryOfPressuresByNodeIndex[self.indexOfEndNodeOfControlledResistor]

		pressureGradientAcrossResistor = startNodePressure - endNodePressure

		if (pressureGradientAcrossResistor > 0.0):
			resistanceToSet = self.resistanceToForwardFlow
		elif (pressureGradientAcrossResistor <= 0.0):
			resistanceToSet = self.maxResistanceWithClosedValve


		return resistanceToSet

