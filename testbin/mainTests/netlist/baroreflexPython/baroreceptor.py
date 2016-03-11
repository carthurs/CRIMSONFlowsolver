from CRIMSONPython import *
from scipy.integrate import ode
import numpy

def loadClassOnRestart(fileName, MPIRank):
	print "Baroreflex controller is not pickleable. Restart will not work properly."
	pass

def saveClassForRestart(objectToSave):
	pass

class parameterController(abstractParameterController):

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank):
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank)
		self.controllerPriority = 0
		self.finishSetup()
		
		# Time parameters
		self.timestepsPerCardiacCycle = 1000
		self.currentTimestepIndex = 0
		startTime = 0.0
		self.currentTime = startTime
		self.periodicTime = 0.0

		# Initial state
		self.initialHeartRate = 60.0 #85.7
		self.heartRate = self.initialHeartRate # we manage the actual heartRate (rather than just a scaling of it) here
		self.heartPeriod = self.heartRate / 60

		initialHeartRateControlSignal = 1.0
		initialMaxElastanceControlSignal = 1.0
		initialPeripheralResistanceControlSignal = 1.0
		initialVenousComplianceControlSignal = 1.0
		initialVenousUnstressedVolumeControlSignal = 1.0
		self.initialStates = numpy.array([initialHeartRateControlSignal, initialMaxElastanceControlSignal, initialPeripheralResistanceControlSignal, initialVenousComplianceControlSignal, initialVenousUnstressedVolumeControlSignal])

		# Baroreflex model parameters
		self.baroreflexOn = True
		self.targetMeanPressure = 100.0 * 133.32
		self.currentMeanPressure = self.targetMeanPressure
		self.autonomicSigmoidSteepness = 5.0

		self.pressureHistoryArray = numpy.array([])

		# ODE solver
		self.odeSolver = ode(self.computeOdeDerivatives)
		self.odeSolver.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
		self.odeSolver.set_initial_value(self.initialStates, self.currentTime)
		self.odeSolver.set_f_params(self.currentMeanPressure)
		# self.odeSolver.set_f_params(constants)


	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.currentTime += delt
		print "mean pressure: ", self.currentMeanPressure/133.32

		indexOfLpnNodeAtAorticSurface = 1
		baroreceptorObservedPressure = dictionaryOfPressuresByNodeIndex[indexOfLpnNodeAtAorticSurface]

		# if the script is slow, consider doing some preallocation and recycling on this array:
		self.pressureHistoryArray = numpy.append(self.pressureHistoryArray, baroreceptorObservedPressure)

		# self.odeSolver.set_f_params()

		# integrate the ODEs up to the current time
		self.odeSolver.integrate(self.currentTime)
		if (not self.odeSolver.successful()):
			print "============================FAILURE IN BAROREFLEX ODE INTEGRATION============================"
		# get the state after the integration
		states = self.odeSolver.y
		
		# Extract the state variables using friendly names, for clarity
		heartRateControlSignal = states[0]
		maximumElastanceControlSignal = states[1]
		peripheralResistanceControlSignal = states[2]
		venousComplianceControlSignal = states[3]
		venousUnstressedVolumeControlSignal = states[4]

		if (self.heartPeriodCompleted(delt)):
			# this is really the mean pressure over the /previous/ cardiac cycle
			self.currentMeanPressure = numpy.mean(self.pressureHistoryArray)
			self.odeSolver.set_f_params(self.currentMeanPressure)
			# reset the array to empty, ready to begin gathering pressures again over the new cycle
			self.pressureHistoryArray = numpy.array([])

			self.heartRate = self.initialHeartRate * heartRateControlSignal
			self.heartPeriod = self.heartRate / 60

		self.clearBroadcastData()
		if (self.baroreflexOn):
			self.addBroadcastVariable('heartRate', self.heartRate)
			self.addBroadcastVariable('maximumElastanceControlSignal', states[1])
			self.addBroadcastVariable('peripheralResistanceControlSignal', states[2])
			self.addBroadcastVariable('venousComplianceControlSignal', states[3])
			self.addBroadcastVariable('venousUnstressedVolumeControlSignal', states[4])
			print "broadcasting ", self.heartRate, states[1:4]
		else:
			self.addBroadcastVariable('heartRate', self.initialHeartRate)
			self.addBroadcastVariable('maximumElastanceControlSignal', 1.0)
			self.addBroadcastVariable('peripheralResistanceControlSignal', 1.0)
			self.addBroadcastVariable('venousComplianceControlSignal', 1.0)
			self.addBroadcastVariable('venousUnstressedVolumeControlSignal', 1.0)
			print "broadcasting ", [self.initialHeartRate, 1.0, 1.0, 1.0, 1.0]

		# This is a little bit of a hack: we just use this script to observe the LV pressure, but it
		# must return a new pressure for the node it's attached to in the netlist_surfaces.dat.
		#
		# In this case, this script is attached to a ground node with pressure zero (node 6).
		# We just return it's original value without chaing it.
		unchangedPressureAtNode6 = dictionaryOfPressuresByNodeIndex[6]
		return unchangedPressureAtNode6

	def computeOdeDerivatives(self, t, y, currentMeanPressure):
		proportionalDistanceFromTargetMeanPressure = currentMeanPressure / self.targetMeanPressure
		sympatheticSignal = self.sympatheticActivation(proportionalDistanceFromTargetMeanPressure)
		parasympatheticSignal = self.parasympatheticActivation(proportionalDistanceFromTargetMeanPressure)
		
		heartRateDeriv = (1.75*sympatheticSignal - 0.25*parasympatheticSignal - y[0])/3.0
		elastanceDeriv = (0.4*sympatheticSignal + 0.8 - y[1])/3.0
		peripheralResistanceDeriv = (0.8*sympatheticSignal + 0.6 - y[2])/3.0
		venousComplianceDeriv = (-0.2*sympatheticSignal + 1.10 - y[3])/30.0
		venousUnstressedVolumeDeriv = (-0.42*sympatheticSignal + 1.21 - y[4])/30.0
		return numpy.array([heartRateDeriv, elastanceDeriv, peripheralResistanceDeriv, venousComplianceDeriv, venousUnstressedVolumeDeriv])

	def sympatheticActivation(self, delta):
		return 1.0/(1 + pow(delta, self.autonomicSigmoidSteepness))

	def parasympatheticActivation(self, delta):
		return 1.0/(1 + pow(delta, - self.autonomicSigmoidSteepness))

	def heartPeriodCompleted(self, delt):
		self.periodicTime += delt/self.heartPeriod

		periodJustReset = False
		if self.periodicTime > 1.0:
			self.periodicTime = 0.0
			periodJustReset = True
		return periodJustReset
