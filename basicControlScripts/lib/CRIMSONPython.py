class stateDataContainer:
	
	def __init__(self, containerNameTag):
		self.containerNameTag = containerNameTag
		self.stateDataInternal = dict() # empty dictionary to store state data itself

	def clear(self):
		self.stateDataInternal.clear()

	def addValue(self,humanReadableNameString, value):
		self.stateDataInternal[humanReadableNameString] = value

	def broadcast(self):
		# add another level to the state data dictionary before returning it, which
		# identifies this script as its source:
		returnValue = dict()
		returnValue[self.containerNameTag] = self.stateDataInternal
		return returnValue

class abstractParameterController:
	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile):
		self.m_baseNameOfThisScript = baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile
		self.m_stateDataContainer = stateDataContainer(self.m_baseNameOfThisScript)

	def finishSetup(self):
		# Handle the case where the user doesn't need broadcasts
		# (and so hasn't implemented setFirstTimestepBroadcastValues)
		if hasattr(self.__class__, 'setFirstTimestepBroadcastValues') and callable(getattr(self.__class__, 'setFirstTimestepBroadcastValues')):
			self.setFirstTimestepBroadcastValues()

	def clearBroadcastData(self):
		self.m_stateDataContainer.clear()

	def addBroadcastVariable(self, humanReadableNameString, value):
		self.m_stateDataContainer.addValue(humanReadableNameString, value)

	def getRecievedBroadcastValue(self, broadcastingScriptBaseName, humanReadableNameString):
		return self.m_stateDataFromAllPythonControllers[broadcastingScriptBaseName][humanReadableNameString]

	def printAllRecievedData(self):
		print "Recieved broadcast from other Python controllers:", self.m_stateDataFromAllPythonControllers

	# This gets called by C++ to gather the boradcast data from this parameter controller.
	def broadcastStateDataToOtherParameterControllers(self):
		return self.m_stateDataContainer.broadcast()

	# This gets called by C++ once it has gathered and bundled the broadcast data from all 
	# the parameter controllers.
	def receiveStateDataFromAllOtherParameterControllers(self, stateDataFromAllPythonControllers):
		self.m_stateDataFromAllPythonControllers = stateDataFromAllPythonControllers

	# If you want to broadcast anything, add a class like this in your controller.
	# Some example entries are commented out here.
	#
	# Note that this gets called before the concrete controller's constructor (__init__() method),
	# so it cannot use any variables defined therein. Define any needed member variables
	# (at least by giving them initial values) in setFirstTimestepBroadcastValues
	#
	# def setFirstTimestepBroadcastValues(self):
	# 	self.clearBroadcastData()
	# 	self.addbroadcastVariable('variableName1', 0.0)
	# 	self.addbroadcastVariable('variableName2', 0.457)
	# 	self.initialValueForSomeVariable = 1.34
	# 	self.addBroadcastVariable('myVariable3', self.initialValueForSomeVariable)