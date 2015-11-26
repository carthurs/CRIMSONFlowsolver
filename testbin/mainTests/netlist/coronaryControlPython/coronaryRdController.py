from CRIMSONPython import *
import numpy

# This is an example general control script. Edit it, and save it with a meaningful name
# (e.g. downstreamResistanceController.py) in your simulation directory. You'll need to refer to it in netlist_surfaces.dat
# or netlist_closed_loop_downstream.dat.
#
# Its purpose is to adjust a parameter (flow, pressure or component parameter (resistance,compliance, etc.))

# Coronary Structure with parameter names:
#
# (3D)--[==R_a==]------[==R_p==]------[==R_d==]-----.
#                   |              |                |
#                 -----          -----              |
#             C_a -----     C_im -----              |
#                   |              |                |
#                  ---           (PLV)			   ---
#                   -                               -



class parameterController(abstractParameterController): #NECESSARY

	def __init__(self, baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank): #NECESSARY
		abstractParameterController.__init__(self,baseNameOfThisScriptAndOfRelatedFlowOrPressureDatFile, MPIRank) #NECESSARY
		self.controllerPriority = -1
		self.m_periodicTime = 0.0
		# self.m_heartPeriod = 0.86
		self.finishSetup() #NECESSARY
		self.setupControlSystemVariables()

	def setupControlSystemVariables(self):
		self.proportionOfLeftMyocardiumThisCoronaryPerfuses = 0.1
		self.initialMVO2 = 2720.0 * self.proportionOfLeftMyocardiumThisCoronaryPerfuses # mm^3. 2720.0 should be set to some resting-state value for your particular application. Simulate to discover it.
		self.previousMVO2 = self.initialMVO2

		self.MVO2History = numpy.array(self.initialMVO2) #empty array
		self.leftVentricularVolumeHistory = numpy.array([]) #empty array
		self.leftVentricularPressureHistory = numpy.array([])

		self.R_d_history = numpy.array([])
		self.R_p_history = numpy.array([])
		self.hungerHistory = numpy.array([])

		self.currentTimestepIndex = 0
		self.beatStartTimestepIndices = numpy.array([self.currentTimestepIndex])

		self.currentBeatIndex = 0

		self.R_a = 0.32155 / self.proportionOfLeftMyocardiumThisCoronaryPerfuses
		self.R_p = 0.64310 / self.proportionOfLeftMyocardiumThisCoronaryPerfuses
		self.R_d = 7.50000 / self.proportionOfLeftMyocardiumThisCoronaryPerfuses

		self.oneOverTotalResistance = 1.0/(self.R_a + self.R_p + self.R_d)

		self.deltaMVO2 = 0.0
		self.feedback_gain = 0.4
		self.alphaFeedforwarGain = -0.1
		self.arterial_O2_volume_proportion = 0.125 # volume percent
		self.dampingCoefficient = 0.2625
		self.minimumAllowedTotalResistance = 0.00010 / self.proportionOfLeftMyocardiumThisCoronaryPerfuses
		self.maximumAllowedTotalResistance = 10000.00000 / self.proportionOfLeftMyocardiumThisCoronaryPerfuses
		self.meanPerfusionPressure = 100.0 * 133.3 # 100 mmHg
		self.hungerSignal = 0.0
		self.debtDelta = 0.0
		self.maximumAllowedNegativeHunger = -250.0 * self.proportionOfLeftMyocardiumThisCoronaryPerfuses # mm^3


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

		leftVentricularPressure = self.getRecievedBroadcastValue('ventricularPressureObserver','leftVentricularPressure')
		leftVentricularVolume = self.getRecievedBroadcastValue('ventricularPressureObserver','leftVentricularVolume')
		self.recordPressureAndVolume(leftVentricularPressure, leftVentricularVolume)

		self.m_heartPeriod = self.getRecievedBroadcastValue('lvElastanceController','heartPeriod')

		if self.currentBeatIndex > 0:
			# only if there's enough data to interpolate:
			lastCompleteBeatStartTime = self.beatStartTimestepIndices[-2] * delt
			lastCompleteBeatEndTime = self.beatStartTimestepIndices[-1] * delt
			xValuesForInterpolator = [lastCompleteBeatStartTime, lastCompleteBeatEndTime]
			
			mvo2AtStartOfLastCompleteBeat = self.MVO2History[-2]
			mvo2AtEndOfLastCompleteBeat = self.MVO2History[-1]
			yValuesForInterpolator = [mvo2AtStartOfLastCompleteBeat, mvo2AtEndOfLastCompleteBeat]

			# We interpolate MVO2 at n% through the previous heart period, and use it when we are n%
			# of the way through the current heart period. This allows us to have a realistic MVO2 gradient,
			# and it introduces a delat of one beat to the MVO2 processing.
			interpolationTime = self.currentTimestepIndex * delt - self.m_heartPeriod
			self.currentMVO2 = numpy.interp(interpolationTime, xValuesForInterpolator, yValuesForInterpolator)
		else:
			#else there's not enough data to interpolate yet
			self.currentMVO2 = self.initialMVO2


		currentBloodFlow = dictionaryOfFlowsByComponentIndex[1] # flow through R_p
		self.debtDelta = self.currentMVO2 - currentBloodFlow * self.arterial_O2_volume_proportion
		self.hungerSignal = max(self.maximumAllowedNegativeHunger, self.hungerSignal + self.debtDelta * delt)
		self.hungerHistory = numpy.append(self.hungerHistory, self.hungerSignal)

		self.stepControlODE(delt)
		self.restrictResistanceToWithinPhysiologicalRange()

		# Set the new resistances:
		self.R_d = 1.0/self.oneOverTotalResistance - self.R_p - self.R_a
		self.R_d_history = numpy.append(self.R_d_history, self.R_d)
		self.R_p = 1.0/(1.0/self.R_p + delt * self.alphaFeedforwarGain * self.deltaMVO2 / self.arterial_O2_volume_proportion / self.meanPerfusionPressure)
		self.R_p_history = numpy.append(self.R_p_history, self.R_p)

		# broadcast the new R_p, to be picked up by the controller for R_p:
		self.clearBroadcastData()
		self.addBroadcastVariable('Rp', self.R_p)

		self.updatePeriodicTime(delt)
		self.finalise()

		# Return the resistance to CRIMSON flowsolver, to be applied to the component this control script is attached to in netlist_surfaces.dat
		# (should be R_d).
		resistanceToSetAtResistorRd = self.R_d
		return resistanceToSetAtResistorRd

	def finalise(self):
		self.previousMVO2 = self.currentMVO2

	def updatePeriodicTime(self, delt):
		
		self.m_periodicTime = self.m_periodicTime + delt
		# Keep m_periodicTime in the range [0,m_heartPeriod):
		if self.m_periodicTime >= self.m_heartPeriod:
			self.m_periodicTime = self.m_periodicTime - self.m_heartPeriod
			self.currentBeatIndex += 1
			self.computeMVO2OverPreviousBeat(self.beatStartTimestepIndices[-1])
			self.beatStartTimestepIndices = numpy.append(self.beatStartTimestepIndices,self.currentTimestepIndex)
			if self.MPIRank == 0:
				numpy.savetxt('R_d_history.dat', self.R_d_history)
				numpy.savetxt('R_p_history.dat', self.R_p_history)
				numpy.savetxt('MVO2History.dat', self.MVO2History)
				numpy.savetxt('myocardialHungerHistory.dat', self.hungerHistory)
		self.currentTimestepIndex += 1

	def stepControlODE(self, delt):
		self.deltaMVO2 = ( self.currentMVO2 - self.previousMVO2 ) / delt
		print "deltaMVO2", self.deltaMVO2

		self.oneOverTotalResistance = self.oneOverTotalResistance + \
						 delt*(
						 (
						 self.feedback_gain/self.arterial_O2_volume_proportion * self.hungerSignal + \
						 self.dampingCoefficient / self.arterial_O2_volume_proportion * self.debtDelta + \
						 self.deltaMVO2 / self.arterial_O2_volume_proportion \
						                  )/ \
						       self.meanPerfusionPressure)
		
		# self.oneOverTotalResistance = self.oneOverTotalResistance + \
		# 				 delt*(
		# 				 (
		# 				 self.feedback_gain/self.arterial_O2_volume_proportion * self.hungerSignal + \
		# 				 self.dampingCoefficient / self.arterial_O2_volume_proportion * self.debtDelta + \
		# 				 self.deltaMVO2 / self.arterial_O2_volume_proportion - \
		# 				 0.0*(part-part_previous) / delt*self.oneOverTotalResistance + \
		# 				 self.feedforward_gain * self.deltaMVO2 / self.arterial_O2_volume_proportion + \
		# 				                  )/ \
		# 				       self.meanPerfusionPressure);

	def restrictResistanceToWithinPhysiologicalRange(self):
		if self.oneOverTotalResistance > 1/(self.minimumAllowedTotalResistance + self.R_p + self.R_a):
			self.oneOverTotalResistance = 1.0/(self.minimumAllowedTotalResistance + self.R_p + self.R_a)
		elif self.oneOverTotalResistance < 1.0/(self.maximumAllowedTotalResistance + self.R_p + self.R_a):
			self.oneOverTotalResistance = 1.0/(self.maximumAllowedTotalResistance + self.R_p + self.R_a)

	def computeMVO2OverPreviousBeat(self, beatStartTimeIndex):
		pressureVolumeAreaOverPreviousBeat = self.computePV_Area(beatStartTimeIndex, self.currentTimestepIndex)
		# units - Joules (converted volume to m^3, and pressure is already in Pa)
		pressureVolumeAreaOverPreviousBeat = pressureVolumeAreaOverPreviousBeat / 1e9
		joulesFrom1cubicMillimetreOfOxygen = 0.020
		pressureVolumeAreaToOxygenRequirementsScaleFactor = 3.0
		MVO2OverPreviousBeat = pressureVolumeAreaOverPreviousBeat * self.proportionOfLeftMyocardiumThisCoronaryPerfuses * pressureVolumeAreaToOxygenRequirementsScaleFactor / \
																															joulesFrom1cubicMillimetreOfOxygen / \
																															 self.m_heartPeriod
		self.MVO2History = numpy.append(self.MVO2History, MVO2OverPreviousBeat)

	def recordPressureAndVolume(self, currentLVPressure, currentLVVolume):
		self.leftVentricularVolumeHistory = numpy.append(self.leftVentricularVolumeHistory, currentLVVolume)
		self.leftVentricularPressureHistory = numpy.append(self.leftVentricularPressureHistory, currentLVPressure)


	#  This function works with a polygon as a set of points (x-coords coords in leftVentricularVolumeHistory, y-coords in leftVentricularPressureHistory)
	#  and returns the area enclosed by the convex hull consisting of the polygoin points PLUS THE ORIGIN.
	#  The sides of the polygon should not intersect. (0,0) should not be included.
	#  This is a MODIFICATION of a standard formula; see eg. http://www.mathopenref.com/coordpolygonarea.html
	#  The modification is to include the origin in the convex hull.
	#  Note that the points MUST be numbered anticlockwise!
	def computePV_Area(self, arrayDataStartIndex, arrayDataEndIndex):
		# Form the sum at the core of the formula:
		accumulator = 0.0
		ii=arrayDataStartIndex
		while (ii<arrayDataEndIndex):
			potentialNextAddition = self.leftVentricularVolumeHistory[ii]*self.leftVentricularPressureHistory[ii+1] - self.leftVentricularVolumeHistory[ii+1]*self.leftVentricularPressureHistory[ii]
			# It is this if-guard that allows the inclusion of the area to the origin.
			if (potentialNextAddition > 0):
			   accumulator = accumulator + potentialNextAddition
			ii += 1
		# Add the last sum term, which doesn't work nicely with the above loop:
		potentialNextAddition = self.leftVentricularVolumeHistory[arrayDataEndIndex]*self.leftVentricularPressureHistory[1] - self.leftVentricularVolumeHistory[1]*self.leftVentricularPressureHistory[arrayDataEndIndex]
		if (potentialNextAddition > 0):
			accumulator += potentialNextAddition

		# Just need half now, to complete the formula:
		# (the absolute value here could be taken if we didn't want the origin included, which would relax the constraint that the points must be anticlockwise.)
		computedPV_Area = accumulator/2.0
		return computedPV_Area