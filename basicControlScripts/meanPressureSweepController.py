# This script can be used to control a resistor component located where we have a diode/valve
# which we wish to allow backflow through, as opposed to zero flow, when the pressure
# gradient across it is negative.

# !!!!!!!!!THIS IS UNFINISHED!!!!!!!!!!!!!!
# I realise that there's no obvious way of reading the pressure we want to compute the mean of - all we have is outlet surface pressures!

from CRIMSONPython import *
import numpy

class meanPressureSweepController(abstractParameterController): #NECESSARY

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank): #NECESSARY
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank) #NECESSARY
		self.finishSetup() #NECESSARY

		# SET THSE YOURSELF!:
		self.targetMeanPressure_mmHg = 
		self.resistanceScalingToSetAfterEachStep = 0.9
		self.cyclePeriodInSteps = 1214
		self.timestepCounter = 0
		self.cardiacCycleCounter = 0
		self.indexOfNodeAt3dInterface = 1

		# these are computed for you:
		self.targetMeanPressure = self.targetMeanPressure_mmHg / 133.32
		self.pressureHistory = numpy.zeros(self.cyclePeriodInSteps)
		self.meanPressureHistory = 

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

		self.pressureHistory[self.timestepCounter] = dictionaryOfPressuresByNodeIndex[self.indexOfNodeAt3dInterface]

		self.timestepCounter += 1

		if (self.timestepCounter % self.cyclePeriodInSteps == 0):
			self.cardiacCycleCounter += 1
			computeAndStoreMeanPressure()

		if (self.timestepCounter % self.cyclePeriodInSteps * self.numberOfCyclesBetweenParameterChanges == 0 ):
			newCompliance = currentParameterValue * self.complianceScalingToSetAfterEachStep
		else:
			newCompliance = currentParameterValue



		return newCompliance

	def computeAndStoreMeanPressure():

