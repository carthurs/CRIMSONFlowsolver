#!/sw/lsa/centos7/python-anaconda2/201607/bin/python

import numpy
import sys
import os

try:
    import pandas
    prettyPrintingAvailable = True
except ImportError:
    prettyPrintingAvailable = False
    print "\n\n ===> (II) Python module pandas not found. Pretty output formatting will not be possible, but the script will still work.\n\n"

if len(sys.argv) < 3:
    print "Computes pressure and flow summary at the outlets from PressHist.gplot and FlowHist.gplot.\n"
    print "\nIf faceInfo.dat, solver.inp and Python module pandas are available, the output will be easier to read.\n"
    print "\nResults will be best when a single pressure cycle is provided.\n"
    print "\nUsage: computeHaemodynamicSummary.py <startTimestepIndex> <endTimestepIndex>\n"
    print "\nThis will compute the pressure mean for each outlet over the timestep index range [startTimestepIndex, endTimestepIndex].\n\n"
    sys.exit(0)

try:
    with open('solver.inp','r') as solverInpFile:
        for lineFromFile in solverInpFile:
            if 'list of output surfaces' in lineFromFile.lower():
                # Split on colon, then split the part after the colon on whitespace (to retrieve the list of output surface indices)
                outputSurfacesList = lineFromFile.split(':')[1].split()
            elif 'time step size' in lineFromFile.lower():
                timestep = float(lineFromFile.split(':')[1])
    solverInpAvailable = True
except IOError:
    prettyPrintingAvailable = False
    solverInpAvailable = False
    print "\n\n ===> (II) No solver.inp found. Pretty output formatting will not be possible, but the script will still work.\n\n"

try:
    with open('faceInfo.dat','r') as faceInfoFile:
        faceIdToFaceNameDictionary = {}
        for lineFromFile in faceInfoFile:
            spaceDelimitedFileLineSplit = lineFromFile.split()
            faceIdToFaceNameDictionary[spaceDelimitedFileLineSplit[-2]] = spaceDelimitedFileLineSplit[-1]
except IOError:
    prettyPrintingAvailable = False
    print "\n\n ===> (II) No faceInfo.dat found. Pretty output formatting will not be possible, but the script will still work.\n\n"

if prettyPrintingAvailable:
    orderedDataLabels = [faceIdToFaceNameDictionary[surfaceIndex] for surfaceIndex in outputSurfacesList]
   
startTimestepIndex = int(sys.argv[1])
endTimestepIndex = int(sys.argv[2])

print "\nMake sure you've generated the PressHist.gplot and FlowHist.gplot files using write_rcr_data.gpi before running this script.\n"
print "\nI'm assuming the model is in mm, so scaling the pressures to mmHg by dividing by 133.32.\n"

pressureFile = numpy.loadtxt("PressHist.gplot")
pressureScaling = 133.32
userRequestedPressureDataRange = pressureFile[startTimestepIndex:endTimestepIndex,:]/pressureScaling

pressureMeansArray = numpy.mean(userRequestedPressureDataRange,axis=0)
pressureMeansArray = numpy.delete(pressureMeansArray, 0) # remove the timestep index mean from the start (it's useless data)
pressureMeansArray = numpy.transpose(pressureMeansArray)

peakPressureArray = numpy.max(userRequestedPressureDataRange, axis=0)
peakPressureArray = numpy.delete(peakPressureArray, 0)
peakPressureArray = numpy.transpose(peakPressureArray)

minimumPressureArray = numpy.min(userRequestedPressureDataRange, axis=0)
minimumPressureArray = numpy.delete(minimumPressureArray, 0)
minimumPressureArray = numpy.transpose(minimumPressureArray)

pulsePressureArray = peakPressureArray - minimumPressureArray

if solverInpAvailable: # (incidentally, if this is true, so is prettyPrentingAvailable)
    print "Detected timestep of ===>", timestep, "<=== from solver.inp. CHECK! If this is incorrect, the reported data duration will be incorrect."

    dataDurationInSeconds = timestep * (endTimestepIndex - startTimestepIndex)

    print "The duration of the requested data range is therefore", dataDurationInSeconds, "seconds."
    
print "\nI'm assuming the model is in mm, so scaling the flows to litres per minute.\n"

flowFile = numpy.loadtxt("FlowHist.gplot")
userRequestedFlowDataRange = flowFile[startTimestepIndex:endTimestepIndex,:] / 1.0e6 * 60

flowMeansArray = numpy.mean(userRequestedFlowDataRange, axis=0)
flowMeansArray = numpy.delete(flowMeansArray, 0)
flowMeansArray = numpy.transpose(flowMeansArray)
    
indexOfInflowInData = numpy.argmin(flowMeansArray) # assume that the inflow is the most negative. Inform the user and warn them that this must be correct.
surfaceIndexOfInflow = outputSurfacesList[indexOfInflowInData]
print "\nI think the inflow is surface ===>", surfaceIndexOfInflow, "<=== with name", faceIdToFaceNameDictionary[surfaceIndexOfInflow], "CHECK! If This is incorrect, the flow splits will be incorrect.\n\n"

meanInflowRateLitresPerMinute = flowMeansArray[indexOfInflowInData]

flowSplits = flowMeansArray / -meanInflowRateLitresPerMinute



print "\nPressures are in mmHg, flows are in litres per minute.\n"
if prettyPrintingAvailable:
    stackedData = numpy.transpose(numpy.row_stack((pressureMeansArray, pulsePressureArray,  peakPressureArray, minimumPressureArray, flowMeansArray, flowSplits)))
    dataFrame = pandas.DataFrame(stackedData, index=orderedDataLabels, columns=['Pressure Means', 'Pulse Pressure', 'Systolic Pressure', 'Diastolic Pressure', 'Mean Flow', 'Proportion Of Inflow'])
    print dataFrame
else:
    print "Means of each column of PressHist:"
    print "----"
    print pressureMeansArray
    print "\n\nMeans of each column of FlowHist:"
    print "----"
    print flowMeansArray

print "\n"


# Gather any interesting files that we also want to summarise:
if prettyPrintingAvailable:
    downstreamPressureFiles = []
    for file in os.listdir('.'):
        if 'netlistPressures_downstream' in file:
           downstreamPressureFiles += [file]

    for downstreamPressureFile in downstreamPressureFiles:
        pressureFile = numpy.loadtxt(downstreamPressureFile)
                
        userRequestedPressureDataRange = pressureFile[startTimestepIndex:endTimestepIndex,:]/pressureScaling

        pressureMeansArray = numpy.mean(userRequestedPressureDataRange,axis=0)
        pressureMeansArray = numpy.delete(pressureMeansArray, 0) # remove the timestep index mean from the start (it's useless data)
        pressureMeansArray = numpy.transpose(pressureMeansArray)

        peakPressureArray = numpy.max(userRequestedPressureDataRange, axis=0)
        peakPressureArray = numpy.delete(peakPressureArray, 0)
        peakPressureArray = numpy.transpose(peakPressureArray)

        minimumPressureArray = numpy.min(userRequestedPressureDataRange, axis=0)
        minimumPressureArray = numpy.delete(minimumPressureArray, 0)
        minimumPressureArray = numpy.transpose(minimumPressureArray)

        pulsePressureArray = peakPressureArray - minimumPressureArray

        stackedData = numpy.transpose(numpy.row_stack((pressureMeansArray, pulsePressureArray,  peakPressureArray, minimumPressureArray)))
        dataFrame = pandas.DataFrame(stackedData, index=range(1, len(pressureMeansArray)+1), columns=['Pressure Means', 'Pulse Pressure', 'Systolic Pressure', 'Diastolic Pressure'])
        print "Nodal pressure data from", downstreamPressureFile + ":\n"
        print dataFrame

if prettyPrintingAvailable:
    downstreamFlowFiles = []
    for file in os.listdir('.'):
        if 'netlistFlows_downstream' in file and not 'componentLabels_' in file:
           downstreamFlowFiles += [file]

    for downstreamFlowFile in downstreamFlowFiles:
        flowFile = numpy.loadtxt(downstreamFlowFile)
                
        userRequestedFlowDataRange = flowFile[startTimestepIndex:endTimestepIndex,:]/ 1.0e6 * 60.0

        flowMeansArray = numpy.mean(userRequestedFlowDataRange,axis=0)
        flowMeansArray = numpy.delete(flowMeansArray, 0) # remove the timestep index mean from the start (it's useless data)
        flowMeansArray = numpy.transpose(flowMeansArray)

        maxFlowArray = numpy.max(userRequestedFlowDataRange, axis=0)
        maxFlowArray = numpy.delete(maxFlowArray, 0)
        maxFlowArray = numpy.transpose(maxFlowArray)

        minimumFlowArray = numpy.min(userRequestedFlowDataRange, axis=0)
        minimumFlowArray = numpy.delete(minimumFlowArray, 0)
        minimumFlowArray = numpy.transpose(minimumFlowArray)

        try:
            with open('componentLabels_'+downstreamFlowFile, 'r') as componentLabelsFile:
                componentLabelsDictionary ={}
                for line in componentLabelsFile:
                    spaceDelimitedLine = line.split()
                    componentLabelsDictionary[spaceDelimitedLine[0]] = spaceDelimitedLine[1]
            
            flowComponentRowLabels = [componentLabelsDictionary[str(index)] for index in range(1,len(flowMeansArray)+1)]
        except IOError:
            flowComponentRowLabels = range(1,len(flowMeansArray)+1)
            print "(II) No custom labels for the components found for file", downstreamFlowFiles, "create one with two columns: component index and friendly nametag, with name", "componentLabels_"+downstreamFlowFile



        stackedData = numpy.transpose(numpy.row_stack((flowMeansArray, maxFlowArray, minimumFlowArray)))
        dataFrame = pandas.DataFrame(stackedData, index=flowComponentRowLabels, columns=['Component Flow Means', 'Peak Flow', 'Minimum Flow'])
        print "\n\nComponent flow data from", downstreamFlowFile + ":\n"
        print "Note that flow sign is dependent upon the component orientation in your netlist circuit specification."
        print dataFrame

