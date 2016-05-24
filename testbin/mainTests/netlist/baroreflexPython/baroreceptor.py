from CRIMSONPython import *
from scipy.integrate import ode
import numpy
from math import cosh, tanh

# def loadClassOnRestart(fileName, MPIRank):
# 	print "Baroreflex controller is not pickleable. Restart will not work properly."
# 	pass

# def saveClassForRestart(objectToSave):
# 	pass

class parameterController(abstractParameterController):

	def __getstate__(self):
		odict = self.__dict__.copy()
		del odict['odeSolver'] # remove things that can't be pickled
		return odict

	def __setstate__(self, odict):
		self.__dict__.update(odict)
		self.setupOdeSolver(self.states)


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
		self.periodIndex = 0

		# Initial state
		self.initialHeartRate = 60.0 #85.7
		self.heartRate = self.initialHeartRate # we manage the actual heartRate (rather than just a scaling of it) here
		self.heartPeriod = self.heartRate / 60

		initialHeartRateControlSignal = 1.0
		initialMaxElastanceControlSignal = 1.0
		initialPeripheralResistanceControlSignal = 1.0
		initialVenousComplianceControlSignal = 1.0
		initialVenousUnstressedVolumeControlSignal = 1.0
		initialHeartRateCentringValue = 1.0
		initialPeripheralResistanceCentringValue = 1.0
		
		initialStates = numpy.array([initialHeartRateControlSignal, initialMaxElastanceControlSignal, initialPeripheralResistanceControlSignal, initialVenousComplianceControlSignal, initialVenousUnstressedVolumeControlSignal, initialHeartRateCentringValue, initialPeripheralResistanceCentringValue])

		# Baroreflex model parameters
		self.baroreflexOn = True
		self.targetMeanPressure = 100.0 * 133.32
		self.currentMeanPressure = self.targetMeanPressure
		self.autonomicSigmoidSteepness = 5.0

		self.pressureHistoryArray = numpy.array([])

		# ODE solver
		self.setupOdeSolver(initialStates)
		self.lastHeartRateDeriv = 0.0
		self.lastPeripheralResistanceDeriv = 0.0

		if self.MPIRank == 0:
			self.fileWriteBufferLength = 1000
			self.fileWriteBuffer = numpy.zeros((self.fileWriteBufferLength,8))
			self.fileWriteBufferNextWriteIndex = 0

	def setupOdeSolver(self, startingState):
		self.odeSolver = ode(self.computeOdeDerivatives)
		self.odeSolver.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=1)
		self.odeSolver.set_initial_value(startingState, self.currentTime)
		self.odeSolver.set_f_params(self.currentMeanPressure)


	def updateControl(self, currentParameterValue, delt, dictionaryOfPressuresByNodeIndex, dictionaryOfFlowsByComponentIndex, dictionaryOfVolumesByComponentIndex):

		self.currentTime += delt
		self.clearBroadcastData()
		print "mean pressure: ", self.currentMeanPressure/133.32

		if self.currentTime < 15.0:
			self.exerciseIntensityInWatts = 0.0
		elif self.currentTime < 60.0:
			self.exerciseIntensityInWatts = 300.0
		else:
			self.exerciseIntensityInWatts = 0.0

		self.addBroadcastVariable('exerciseIntensityInWatts', self.exerciseIntensityInWatts)

		self.targetMeanPressure = 133.32 * (80.0 + 40.0/300.0 * self.exerciseIntensityInWatts) # see "physical exercise - hemodynamics" in Evernote
		self.targetHeartRate = 60 + 130.0/300.0 * self.exerciseIntensityInWatts # see "physical exercise - hemodynamics" in Evernote
		self.targetHeartRateNormalised = self.targetHeartRate / self.initialHeartRate

		initialPeripheralResistance = 0.1 # should match whatever you set in the boundary conditions
		initialPeripheralConductance = 1.0/initialPeripheralResistance
		targetPeripheralConductanceDelta = max(-5.0 * self.exerciseIntensityInWatts / 150.0, -5.0)
		targetPeripheralConductance = initialPeripheralConductance + targetPeripheralConductanceDelta
		self.targetPeripheralResistanceNormalised = (1.0 / targetPeripheralConductance) / initialPeripheralResistance

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
		self.states = self.odeSolver.y # we store self.states as a member variable so that it gets pickled, and can be used as an initial condition upon unpickling
		
		# Extract the state variables using friendly names, for clarity
		heartRateControlSignal = self.states[0]
		maximumElastanceControlSignal = self.states[1]
		peripheralResistanceControlSignal = self.states[2]
		venousComplianceControlSignal = self.states[3]
		venousUnstressedVolumeControlSignal = self.states[4]

		print {"targetHeartRate": self.targetHeartRate, "targetHeartRateNormalised": self.targetHeartRateNormalised, "heartRateControlSignal": heartRateControlSignal}

		if (self.heartPeriodCompleted(delt)):
			# this is really the mean pressure over the /previous/ cardiac cycle
			self.currentMeanPressure = numpy.mean(self.pressureHistoryArray)
			self.odeSolver.set_f_params(self.currentMeanPressure)
			# reset the array to empty, ready to begin gathering pressures again over the new cycle
			self.pressureHistoryArray = numpy.array([])

			self.heartRate = self.initialHeartRate * heartRateControlSignal
			self.heartPeriod = self.heartRate / 60


		if (self.baroreflexOn):
			self.addBroadcastVariable('heartRate', self.heartRate)
			self.addBroadcastVariable('maximumElastanceControlSignal', self.states[1])
			self.addBroadcastVariable('peripheralResistanceControlSignal', self.states[2])
			self.addBroadcastVariable('venousComplianceControlSignal', self.states[3])
			self.addBroadcastVariable('venousUnstressedVolumeControlSignal', self.states[4])
			print "broadcasting ", [self.heartRate, self.states[1], self.states[2], self.states[3], self.states[4]], "not: ", self.states[5], self.states[6]
		else:
			self.addBroadcastVariable('heartRate', self.initialHeartRate)
			self.addBroadcastVariable('maximumElastanceControlSignal', 1.0)
			self.addBroadcastVariable('peripheralResistanceControlSignal', 1.0)
			self.addBroadcastVariable('venousComplianceControlSignal', 1.0)
			self.addBroadcastVariable('venousUnstressedVolumeControlSignal', 1.0)
			print "broadcasting ", [self.initialHeartRate, 1.0, 1.0, 1.0, 1.0], "not: ", self.states[5], self.states[6]

		if self.MPIRank == 0:
			self.fileWriteBuffer[self.fileWriteBufferNextWriteIndex, 0] = self.currentTimestepIndex
			self.fileWriteBuffer[self.fileWriteBufferNextWriteIndex, 1:-1] = self.states[0:-1]
			# irrifating that we have to do this, but can't slice to last element:
			# see e.g. http://stackoverflow.com/questions/15627312/what-value-do-i-use-in-a-slicing-range-to-include-the-last-value-in-a-numpy-arra
			self.fileWriteBuffer[self.fileWriteBufferNextWriteIndex, -1] = self.states[-1]
			self.fileWriteBufferNextWriteIndex += 1
			if self.fileWriteBufferNextWriteIndex % self.fileWriteBufferLength == 0:
				self.fileWriteBufferNextWriteIndex = 0
				outputFile = open('states.dat','a')
				numpy.savetxt(outputFile, self.fileWriteBuffer, delimiter=' ')
				outputFile.close()

		self.currentTimestepIndex += 1

		# This is a little bit of a hack: we just use this script to observe the LV pressure, but it
		# must return a new pressure for the node it's attached to in the netlist_surfaces.dat.
		#
		# In this case, this script is attached to a ground node with pressure zero (node 6).
		# We just return it's original value without chaing it.
		unchangedPressureAtNode6 = dictionaryOfPressuresByNodeIndex[6]
		return unchangedPressureAtNode6

	def computeOdeDerivatives(self, t, y, currentMeanPressure):
		# unpack extraPArameters:
		# currentMeanPressure = extraParameters[0]


		proportionalDistanceFromTargetMeanPressure = currentMeanPressure / self.targetMeanPressure
		sympatheticSignal = self.sympatheticActivation(proportionalDistanceFromTargetMeanPressure)
		parasympatheticSignal = self.parasympatheticActivation(proportionalDistanceFromTargetMeanPressure)
		
		# Kev's original formulation
		# heartRateDeriv = (1.75*sympatheticSignal + 0.25*parasympatheticSignal - y[0])/3.0
		# elastanceDeriv = (0.4*sympatheticSignal + 0.8 - y[1])/3.0
		# peripheralResistanceDeriv = (0.8*sympatheticSignal + 0.6 - y[2])/3.0
		# venousComplianceDeriv = (-0.2*sympatheticSignal + 1.10 - y[3])/30.0
		# venousUnstressedVolumeDeriv = (-0.42*sympatheticSignal + 1.21 - y[4])/30.0

		# My formulation:
		heartRateCentringDeriv = (self.targetHeartRateNormalised - y[5]) / 5.0#numpy.absolute(heartRateDeriv)
		peripheralResistanceCentringDeriv = (self.targetPeripheralResistanceNormalised - y[6]) / 5.0#numpy.absolute(peripheralResistanceDeriv)

		gradientDropoffRate = 2.0 # configures the value of 1/cosh^2 (=sech^2) where we consider further movement to be "negligible" (a value of 2.0 means the derivatives will be scaled by 0.07 when their parameters reach 100% of their "extremal" range - note that this makes the extrema soft.)

		maxProportionalHeartRate = 1.75
		heartRateDeriv = (((sympatheticSignal-0.5) - (parasympatheticSignal-0.5)) / pow(cosh((y[0] - y[5])*gradientDropoffRate/(maxProportionalHeartRate - 1.0)),2) + (y[5]-y[0]) * (1.0 + tanh(-(maxProportionalHeartRate - 1.0) + 2000.0 * (self.lastHeartRateDeriv*pow(heartRateCentringDeriv,2) - pow(heartRateCentringDeriv,3) ))))/3.0
		self.lastHeartRateDeriv = heartRateDeriv

		elastanceDeriv = ((sympatheticSignal - 0.5)) / pow(cosh((y[1] - 1.0)*gradientDropoffRate/0.2),2)/3.0

		maxProportionalResistance = 1.8
		peripheralResistanceDeriv = (((sympatheticSignal - 0.5)) / pow(cosh((y[2] - y[6])*gradientDropoffRate/(maxProportionalResistance - 1.0)),2) + (y[6] - y[2]) * (1.0 + tanh(-(maxProportionalResistance-1.0) + 2000.0 * (self.lastPeripheralResistanceDeriv*pow(peripheralResistanceCentringDeriv,2) - pow(peripheralResistanceCentringDeriv,3) )))  )/3.0 # 0.8 from 0.4
		self.lastPeripheralResistanceDeriv = peripheralResistanceDeriv
		venousComplianceDeriv = (-(sympatheticSignal - 0.5)) / pow(cosh((y[3] - 1.0)*gradientDropoffRate/0.5),2)/30.0 #0.5 from 0.1
		venousUnstressedVolumeDeriv = (-(sympatheticSignal - 0.5)) / pow(cosh((y[4] - 1.0)*gradientDropoffRate/0.5),2)/30.0 # 0.5 from 0.21

		return numpy.array([heartRateDeriv, elastanceDeriv, peripheralResistanceDeriv, venousComplianceDeriv, venousUnstressedVolumeDeriv, heartRateCentringDeriv, peripheralResistanceCentringDeriv])

	def sympatheticActivation(self, delta):
		return 1.0/(1.0 + pow(delta, self.autonomicSigmoidSteepness))

	def parasympatheticActivation(self, delta):
		return 1.0/(1.0 + pow(delta, - self.autonomicSigmoidSteepness))

	def heartPeriodCompleted(self, delt):
		self.periodicTime += delt/self.heartPeriod
		print "periodIndex", self.periodIndex

		periodJustReset = False
		if self.periodicTime > 1.0:
			self.periodicTime = 0.0
			self.periodIndex += 1
			periodJustReset = True
		return periodJustReset
