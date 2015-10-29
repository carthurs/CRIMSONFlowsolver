from CRIMSONPython import *
from math import pi, cos

# This is an example general control script. Edit it, and save it with a meaningful name
# (e.g. downstreamResistanceController.py) in your simulation directory. You'll need to refer to it in netlist_surfaces.dat
# or netlist_closed_loop_downstream.dat.
#
# Its purpose is to adjust a parameter (flow, pressure or component parameter (resistance,compliance, etc.))
#
# Lines which must not be changed are marked with "#NECESSARY" at the end.

class parameterController(abstractParameterController): #NECESSARY

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank): #NECESSARY
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank) #NECESSARY
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_heartPeriod = 0.86;
		self.finishSetup() #NECESSARY

	def setFirstTimestepBroadcastValues(self):
		self.clearBroadcastData()
		self.addBroadcastVariable('three', 3) # just a non-functional example broadcast
		self.addBroadcastVariable('four', 4) # just a non-functional example broadcast

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

		self.clearBroadcastData()
		self.addBroadcastVariable('three', 3) # just a non-functional example broadcast
		self.addBroadcastVariable('four', 4) # just a non-functional example broadcast
		# A more meaningful broadcast; if the present script is called myScript.py, 
		# then this will be accessible in other scripts via a call 
		# self.getRecievedBroadcastValue('myScript','firstComponentsFlow') in all other
		# control scripts.
		flowFromAssociatedNetlistToBroadcast = dictionaryOfFlowsByComponentIndex[1]
		self.addBroadcastVariable('firstComponentsFlow', flowFromAssociatedNetlistToBroadcast)

		# Can use this to see the data recieved by broadcast:
		#
		# print "Custom Controller Reporting!"
		# self.printAllRecievedData()

		# Do something with some data:
		# use a variable tagged "six", broadcast by a controller called elastanceController.py (which we expect to have the value 6 here), and
		# use a variable tagged "eight", broadcast by a controller called nodeController_downstream.py (which we expect to have the value 8 here)
		if self.getRecievedBroadcastValue('elastanceController','six') == 6 and self.getRecievedBroadcastValue('nodeController_downstream','eight') == 8:
			self.updatePeriodicTime(delt)

		# Do something with some data:
		previousResistance = currentParameterValue # we rename here just to make it clear that we're imagining this script to be controlling a resistance
		resistanceToSet = previousResistance * self.getRecievedBroadcastValue('masterController','masterControlSignal')

		return resistanceToSet

	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod
