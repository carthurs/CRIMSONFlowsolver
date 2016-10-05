import copy
import cPickle

# This serializes the class so it can be re-loaded on restart.
# OVERRIDE IN YOUR SUBCLASS IF YOU DONT WANT SERIALISATION
# (I.E. IF YOU'RE USING UNSERIALISEABLE MEMBERS)
def saveClassForRestart(objectToSave):
	if objectToSave.MPIRank == 0:
		fullFileName = objectToSave.m_baseNameOfThisScript + objectToSave.controllerNameQualification + ".pickle"
		outputFile = open(fullFileName, "w")

		cPickle.dump(objectToSave, outputFile, cPickle.HIGHEST_PROTOCOL)
		outputFile.close()

# OVERRIDE (-READ THIS WHOLE COMMENT BEFORE DOING THAT-) IN YOUR SUBCLASS IF YOU DONT WANT SERIALISATION
# (I.E. IF YOU'RE USING UNSERIALISEABLE MEMBERS)
# RETURN None IN THAT CASE, SO C++ WILL CONSTRUCT THE OBJECT
# AFRESH ON RESTART.
#
# Alternatively, remove the things that you can't serialise, then load them again manually on unpickle using something like
# this in your subclass:
#
# def __getstate__(self):
# 	odict = self.__dict__.copy()
# 	del odict['odeSolver'] # remove things that can't be pickled
# 	return odict
#
# def __setstate__(self, odict):
# 	self.__dict__.update(odict)
# 	self.setupOdeSolver(self.states) # On unpickle, manually recreate the things that you couldn't pickle
#
def loadClassOnRestart(fileName, MPIRank):
	inputFile = open(fileName + ".pickle")
	loadedObject = cPickle.load(inputFile)
	inputFile.close()
	# correct the MPI rank (the loaded class was pickled by the rank-zero thread, so need to reset it appropriately for each rank now).
	loadedObject.MPIRank = MPIRank
	return loadedObject


class stateDataContainer(object):
	
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
		# deepcopy it so the data can't be changed after broadcast.
		return copy.deepcopy(returnValue)


class abstractParameterController(object):
	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		self.m_baseNameOfThisScript = baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile
		self.MPIRank = MPIRank
		self.m_stateDataContainer = stateDataContainer(self.m_baseNameOfThisScript)
		self.controllerPriority = 0 # See the method getControllerPriority, below. Override this in your controller if you want to change priorities.

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

	# Priorities are used by c++ to determine the order in which controllers get updated
	# on each time step. Controllers with values closer to -infinity have higher priortiy,
	# numbers closer to +infinity are lower priority.
	#
	# This is useful e.g. when using the control broadcast system, to ensure the value
	# you want has been broadcast before the time when you want to use it in another
	# controller.
	def getControllerPriority(self):
		return self.controllerPriority

	# Package any other new data that you need to give to the controllers here.
	# (and make corresponding changes in the C++)
	#
	# This saves us from putting any new variables in the constructor (in Python),
	# which would break backward-compatibility.
	def recieveExtraData(self, controllerNameQualification):
		self.controllerNameQualification = controllerNameQualification

	# If you want to broadcast anything fixed on the first timestep, add a class like this in your controller.
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