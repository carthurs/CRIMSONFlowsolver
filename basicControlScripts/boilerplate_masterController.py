from CRIMSONPython import *
from math import pi, cos

# This is an example masterController.py script. Edit it, and save it with exactly that
# name in your simulation directory.
#
# Its purpose is to broadcast master control signals to the other parameter controllers
# (e.g. an exercise state or a nerual command signal).
#
# Lines which must not be changed are marked with "#NECESSARY" at the end.

class parameterController(abstractParameterController): #NECESSARY

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile): #NECESSARY
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile) #NECESSARY
		self.m_periodicTime = 0.0; #\todo think about this for restarts!
		self.m_heartPeriod = 0.86;
		self.controlSignal = 0.0 #an intiial value
		self.finishSetup() #NECESSARY

	def setFirstTimestepBroadcastValues(self): #NECESSARY - IF AND ONLY IF you have control broadcasts, you must initialise their values here
		self.clearBroadcastData()
		self.addBroadcastVariable('masterControlSignal', self.controlSignal)
		self.someValue = 3.21 # we ensure that someValue is initialised before we try to broadcast it on the first tep
		self.addBroadcastVariable('anotherBroadcastValue', self.someValue)
		
	def updateControl(self, delt): #NECESSARY - note that the masterController.py (i.e. this example controller)
		
		# Can use this to see the data recieved by broadcast:
		#
		# print "Master Controller Reporting!"
		# self.printAllRecievedData()

		# Do something with some data:
		# use a variable tagged "six", broadcast by a controller called elastanceController.py (which we expect to have the value 6 here), and
		# use a variable tagged "eight", broadcast by a controller called nodeController_downstream.py (which we expect to have the value 8 here)
		if self.getRecievedBroadcastValue('elastanceController','six') == 6 and self.getRecievedBroadcastValue('nodeController_downstream','eight') == 8:
			self.updatePeriodicTime(delt)

		self.controlSignal = cos(self.m_periodicTime)

		# Update any which change each iteration:
		self.clearBroadcastData()
		self.addBroadcastVariable('masterControlSignal', self.controlSignal)

		# Note that masterController.py returns nothing (because it doesn't have a particular node/component whose)
		# parameter it is adjusting
		

	# This is just a convenience function - it's not fundamentally necessary!
	def updatePeriodicTime(self, delt):
	
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod