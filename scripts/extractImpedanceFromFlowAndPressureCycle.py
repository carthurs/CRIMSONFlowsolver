#!/usr/bin/python

# This file makes an impedance input function and the corresponding initial pressure history
# for the C++ impedance boundary condition.
#
# Inputs:
# Pressure and flow waveforms, ideally for several full cycles.
# (This script will extract these from the standard PressHist.dat and
# and FlowHist.dat, although you should change to use the .gplot files if you have more than
# 3 flow surfaces, due to Fortran's output stacking in 3 columns)
#
#
# Outputs:
#
# impedance.dat
# -> Contains the time-domain impedance computed from the input pressure and flow waveforms
# -> Should be renamed time_varying_impedance_surface_X.dat, where X is the index of the
#	 associated flow surface in your geometry.
# 
# impedanceBCFlowHist.dat
# -> Contains the initial flow history, which will be used by the impedance boundary condition.
#	 This will be updated by the simulation as the recent flow history changes.
# -> Should be renamed Qhistor_impedance_surface_X.dat, where X is the index of the
#	 associated flow surface in your geometry.
#
# impedancePressure.dat
# -> Contains the pressure history corresponding to the flow history used in the impedance
#	 creation.
# -> Not required for simulation.
#
# impedanceAnalyticInFrequencyDomain.dat
# -> If you set the Windkessel parameters Rp, Rd and C in the below code, this will contain the
# 	 analytic frequency domain impedance of that idealised Windkessel
# -> Not required for simulation.
#
# impedanceAnalyticInTimeDomain.dat
# -> If you set the Windkessel parameters Rp, Rd and C in the below code, this will contain the
# 	 analytic time domain impedance of that idealised Windkessel
# -> This can be used in place of impedance.dat if you wish to use an analytic windkessel
#	 implemented in terms of its impedance
#
# analyticalConvolution.dat
# -> Contains the convolution of the analytical Windkessel impedance (as in
#    impedanceAnalyticInTimeDomain.dat) with the input flow waveform.
#    Its values indicate roughly what the pressure waveform will look like if 
#    impedanceAnalyticInTimeDomain.dat is used for simulation.
# -> Useful for understanding what the convolution is doing, and for troubleshooting.
# -> Not required for simulation.
#
# convolution.dat
# -> Contains the convolution of the numerically-computed impedance (as in impedance.dat) with the
#	 input flow waveform, Its values indicate roughly what the pressure waveform will look like if 
#    impedance.dat is used for simulation.
# -> Useful for understanding what the convolution is doing, and for troubleshooting.
# -> Not required for simulation.
import numpy

# Define the cardiac cycles to extract:
cycleStartTimestepIndex = 15000#917
cycleEndTimestepIndex = 16720#15860#16720#1781

computeAnalyticWindkesselImpedance = True

timestep = 0.001 #seconds - as in solver.inp for the original simulation
print "\nUsing a timestep of ", timestep, "This should match the solver.inp value of the simulation which generated the input data.\n"

# Define which surface's flow and pressure data to use: zero-indexed column to use in the data
columnToUse = 1 # ZERO-INDEXED! 0 = first column, 1 = second column, ....

pressureData = numpy.loadtxt('PressHist.dat', skiprows=1) # skiprows to ignore the first line of the file (which contains metadata)
flowData = numpy.loadtxt('FlowHist.dat', skiprows=1)

extractedPressureCycles = pressureData[cycleStartTimestepIndex-1:cycleEndTimestepIndex-1,columnToUse]
extractedFlowCycles = flowData[cycleStartTimestepIndex-1:cycleEndTimestepIndex-1,columnToUse]

pressureFft = numpy.fft.fft(extractedPressureCycles)
flowFft = numpy.fft.fft(extractedFlowCycles)

impedanceInFrequencyDomain = pressureFft / flowFft

numpy.savetxt('impedanceInFreqencyDomain.dat', numpy.abs(impedanceInFrequencyDomain))

print "Computed impedanceInFrequencyDomain shape: ", numpy.shape(impedanceInFrequencyDomain)
print "zeroth mode numerical impedance: ", impedanceInFrequencyDomain[0]

impedanceInTimeDomain = numpy.fft.ifft(impedanceInFrequencyDomain)

numberOfSteps = cycleEndTimestepIndex - cycleStartTimestepIndex
endTime = numberOfSteps * timestep
timestepArray = numpy.linspace(0.0, endTime, numberOfSteps)

# print numpy.shape(timestepArray), numpy.shape(impedanceInTimeDomain)

stack = numpy.column_stack((timestepArray, numpy.abs(impedanceInTimeDomain)))

numpy.savetxt('impedance.dat', stack)
numpy.savetxt('impedanceBCFlowHist.dat', extractedFlowCycles)
numpy.savetxt('impedancePressure.dat', extractedPressureCycles)

# We must use "valid" here, because otherwise the discrete convolution will contain values where the two arrays do not entirely
# overlap (including values where they only overlap in a single point, 2 pts, 3 pts,... etc) which gives pressures v. close to zero
# initially.
#
# We also need the flow cycles trace to be longer than the impedance for this to work, so we just append to the input flow data
# a copy of itself. This is fine as we're assuming periodicity anyway.
doubleRepeatOfInputFlowCycles = numpy.append(extractedFlowCycles, extractedFlowCycles)
numpy.savetxt('convolution.dat', numpy.abs(numpy.convolve(impedanceInTimeDomain, doubleRepeatOfInputFlowCycles, mode='valid')))

if computeAnalyticWindkesselImpedance:
	# Analytic version for a Windkessel:
	Rp = 0.3815889462107935
	Rd = 1.740713618386789
	C = 0.1885526988247284

	print "\nComputing analytic Windkessel impedance too..."
	print "Using Rp, Rd and C: ", Rp, Rd, C, "\n"

	# Deteermine the frequencies at which we need to compute the value of the frequency-domain impedance.
	# We match the length with the numerical impedance we computed above, just for convenience.
	frequencyDomainFrequencies = numpy.fft.fftfreq(numpy.size(impedanceInTimeDomain), d=timestep/2/numpy.pi)
	print "First few analytic Windkessel frequencies: ", frequencyDomainFrequencies[0:10], "total: ", numpy.shape(frequencyDomainFrequencies)

	# Compute the impedance at each frequency of interest:
	analyticImpedance = numpy.array([])
	for frequencyIndex in range(0,numpy.size(frequencyDomainFrequencies)):
		impedanceAtThisFrequency = (Rp + Rd + frequencyDomainFrequencies[frequencyIndex] * 1j * Rp * Rd * C) / (1 + frequencyDomainFrequencies[frequencyIndex] * 1j * Rd * C)
		analyticImpedance = numpy.append(analyticImpedance, impedanceAtThisFrequency)

	numpy.savetxt('impedanceAnalyticInFrequencyDomain.dat', numpy.abs(analyticImpedance))

	analyticImpedanceInTimeDomain = numpy.fft.ifft(analyticImpedance)

	numberOfAnalyticSteps = numpy.size(frequencyDomainFrequencies)
	timestepArrayAnalytic = numpy.linspace(0.0, numberOfAnalyticSteps * timestep, numberOfAnalyticSteps)
	# print numpy.shape(timestepArrayAnalytic), numpy.shape(analyticImpedanceInTimeDomain)

	analyticStack = numpy.column_stack((timestepArrayAnalytic, numpy.abs(analyticImpedanceInTimeDomain)))
	numpy.savetxt('impedanceAnalyticInTimeDomain.dat', analyticStack)

	# We must use "valid" here, because otherwise the discrete convolution will contain values where the two arrays do not entirely
	# overlap (including values where they only overlap in a single point, 2 pts, 3 pts,... etc) which gives pressures v. close to zero
	# initially.
	analyticalConvolution = numpy.convolve(doubleRepeatOfInputFlowCycles, analyticImpedanceInTimeDomain, mode='valid')
	numpy.savetxt('analyticalConvolution.dat', numpy.abs(analyticalConvolution))