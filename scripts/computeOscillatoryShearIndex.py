#!/usr/bin/python
from vtk import *

# To use this script, first use the Merge Blocks filter in Paraview on a loaded time-data series, containing
# Wall Shear Stress. Then use Extract Surface filter, then save the data as a .vtk file, using write all timesteps as file series in the dialogue.
# Set that file's base name as fileBaseName in this script, set the startStep and endStep indices to integrate in time over
# (as indexed in the Paraview scene stepping tool), then ensure that indexOfWSSArrayInInput is set to the
# array index containing the wall shear stress in the input data (the script will give you a hint for this if
# you try running it - it loads the data and prints some metadata regarding the input file to the console).
#
# Then run this script from the folder containing that .vtk file.
#
# Then load the file osi.vtk in Paraview.

def readVtkFile(fileName):
	reader = vtkUnstructuredGridReader()
	reader.SetFileName(fileName)

	reader.Update()
	fullReaderResultOutput = reader.GetOutput()
	# scalar_range = fullReaderResultOutput.GetScalarRange()

	return fullReaderResultOutput

def vectorMagnitude(inputThreeTuple):
	sumOfSquares = pow(inputThreeTuple[0], 2) + pow(inputThreeTuple[1], 2) + pow(inputThreeTuple[2], 2)
	return pow(sumOfSquares, 0.5)

def addTimestepDataToIntegral(timestepIndex, integrandScaling):
	currentFileName = fileNameBase + str(timestepIndex) + '.vtk'

	readDataStructure = readVtkFile(currentFileName)
	wssdata = readDataStructure.GetPointData().GetArray(indexOfWSSArrayInInput)

	for nodalDataVectorIndex in range(0, dataSize):
		# Add the 3-element vector tuples componentwise:
		sumData = [integrandScaling * (wssdata.GetTuple3(nodalDataVectorIndex)[i] +  timeIntegratedWSSArray.GetTuple3(nodalDataVectorIndex)[i]) for i in range(0,3)]
		timeIntegratedWSSArray.SetTuple3(nodalDataVectorIndex, sumData[0], sumData[1], sumData[2])

	# use the pressure scalar array to store the magnitude integral:
	for nodalDataVectorIndex in range(0, dataSize):
		WSSvector = wssdata.GetTuple3(nodalDataVectorIndex)
		summandSoFar = timeIntegratedMagnitudeOfWSSArray.GetTuple1(nodalDataVectorIndex)
		timeIntegratedMagnitudeOfWSSArray.SetTuple1(nodalDataVectorIndex, summandSoFar + integrandScaling * vectorMagnitude(WSSvector))


def createVtkFloatArrayAndAddToUnstructuredGridData(dataSize, arrayName):
	newArray = vtkDoubleArray()
	newArray.SetName(arrayName)
	newArray.SetNumberOfValues(dataSize)
	firstStepVtkFileDataStructure.GetPointData().AddArray(newArray)
	return newArray

def createVtk3VectorArrayAndAddToUnstructuredGridData(dataSize, arrayName):
	newArray = vtkDoubleArray()
	newArray.SetNumberOfComponents(3)
	newArray.SetNumberOfTuples(dataSize)
	newArray.SetName(arrayName)
	firstStepVtkFileDataStructure.GetPointData().AddArray(newArray)
	return newArray

if __name__ == '__main__':
	# fileNameBase = 'vtkFileTestPulmonaries_'
	fileNameBase = 'longWideLandingZone_'

	startStep = 0
	endStep = 67

	stepsBetweenDataFiles = 100
	simulationDeltaT = 0.0001
	deltaT = simulationDeltaT * stepsBetweenDataFiles

	firstStepFileName = fileNameBase + str(startStep) + '.vtk'
	firstTimestepOutput = None
	firstStepVtkFileDataStructure = readVtkFile(firstStepFileName)

	# Can work this out using:
	print "\nInput file metadata (can use to find the Array index of the Wall Shear Stress data in the file:\n"
	print firstStepVtkFileDataStructure.GetAttributes(0)
	indexOfWSSArrayInInput = 3

	firstTimestepWSS = firstStepVtkFileDataStructure.GetPointData().GetArray(indexOfWSSArrayInInput)
	dataSize = firstTimestepWSS.GetNumberOfTuples()
	print "dataSize", dataSize

	# Create the array to use to store the time integrated WSS, and copy the first timestep data into it:
	timeIntegratedWSSArray = createVtk3VectorArrayAndAddToUnstructuredGridData(dataSize, 'Time-Integrated Wall Shear Stress')
	for nodalDataVectorIndex in range(0, dataSize):
		firstStepWSSPointVector = firstTimestepWSS.GetTuple3(nodalDataVectorIndex)
		# Here we copy, but also we prepare the WSS vector data itself for integration, by scaling the first step by 0.5 (for Trapesium rule integration)
		timeIntegratedWSSArray.SetTuple3(nodalDataVectorIndex, firstStepWSSPointVector[0] * 0.5, firstStepWSSPointVector[1] * 0.5, firstStepWSSPointVector[2] * 0.5)

	timeIntegratedMagnitudeOfWSSArray = createVtkFloatArrayAndAddToUnstructuredGridData(dataSize, 'Time-Integrated Magnitude of WSS')

	# use the pressure scalar array to store the magnitude integral:
	for nodalDataVectorIndex in range(0, dataSize):
		firstStepWSSvector = timeIntegratedWSSArray.GetTuple3(nodalDataVectorIndex)
		timeIntegratedMagnitudeOfWSSArray.SetTuple1(nodalDataVectorIndex, 0.5 * vectorMagnitude(firstStepWSSvector)) # 0.5 because it's the first step in a Trapesium rule integration

	# Get the maximum WSS magnitude in and time:
	spaceTimeWSSMagnitudeMaximum = 0.0
	for timestepIndex in range(startStep, endStep+1):
		currentFileName = fileNameBase + str(timestepIndex) + '.vtk'
		readDataStructure = readVtkFile(currentFileName)
		wssdata = readDataStructure.GetPointData().GetArray(indexOfWSSArrayInInput)

		for nodalDataVectorIndex in range(0, dataSize):
			wssPointMagnitude = vectorMagnitude(wssdata.GetTuple3(nodalDataVectorIndex))
			if wssPointMagnitude > spaceTimeWSSMagnitudeMaximum:
				spaceTimeWSSMagnitudeMaximum = wssPointMagnitude
	print fileNameBase, "spaceTimeWSSMagnitudeMaximum:", spaceTimeWSSMagnitudeMaximum



	# Actually loop the timesteps and do the core of the integration:
	for timestepIndex in range(startStep + 1, endStep):
		integrandScaling = 1.0
		addTimestepDataToIntegral(timestepIndex, integrandScaling)


		# Can use wssdata.SetTuple3(0,1.0,2.0,3.0) to reset the 0th tuple to (1.0,2.0,3.0)

	# do the final timestep (which differs in that it has 0.5 scaling on the integrand values, for the Trapesium Rule)
	integrandScaling = 0.5
	addTimestepDataToIntegral(timestepIndex, integrandScaling)


	# Include the (uniform) time scaling in the integrals:
	for nodalDataVectorIndex in range(0, dataSize):
		timeScaledData = [deltaT * timeIntegratedWSSArray.GetTuple3(nodalDataVectorIndex)[i] for i in range(0,3)]
		timeIntegratedWSSArray.SetTuple3(nodalDataVectorIndex, timeScaledData[0], timeScaledData[1], timeScaledData[2])

		timeScaledWssMagnitudeIntegral = deltaT * timeIntegratedMagnitudeOfWSSArray.GetTuple1(nodalDataVectorIndex)
		timeIntegratedMagnitudeOfWSSArray.SetTuple1(nodalDataVectorIndex, timeScaledWssMagnitudeIntegral)

	osi_array = createVtkFloatArrayAndAddToUnstructuredGridData(dataSize, 'Oscillatory Shear Index')
	# Compute the OSI:
	for dataIndex in range(0, dataSize):
		wssMagnitudeIntegral = timeIntegratedMagnitudeOfWSSArray.GetTuple1(dataIndex)
		magnitudeOfWssIntegral = vectorMagnitude(timeIntegratedWSSArray.GetTuple3(dataIndex))
		if wssMagnitudeIntegral == 0:
			OSI_point_value = 0.0 # float('NaN')
		else:
			OSI_point_value = 0.5 * (1.0 -  magnitudeOfWssIntegral / wssMagnitudeIntegral)
		osi_array.SetTuple1(dataIndex, OSI_point_value)

	writer = vtkUnstructuredGridWriter()
	writer.SetFileName('osi.vtk')
	writer.SetInputData(firstStepVtkFileDataStructure)
	writer.Write()