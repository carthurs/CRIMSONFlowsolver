#!/sw/lsa/centos7/python-anaconda2/201607/bin/python

import numpy
import sys

try:
    import pandas
    prettyPrintingAvailable = True
except ImportError:
    prettyPrintingAvailable = False
    print "\n\n ===> (II) Python module pandas not found. Pretty output formatting will not be possible, but the script will still work.\n\n"

if len(sys.argv) < 3:
    print "Computes pressure mean from PressHist.gplot.\n"
    print "\nIf faceInfo.dat, solver.inp and Python module pandas are available, the output will be easier to read.\n"
    print "\nResults will be best when a single pressure cycle is provided.\n"
    print "\nUsage: computePressureMeans.py <startTimestepIndex> <endTimestepIndex>\n"
    print "\nThis will compute the pressure mean for each outlet over the timestep index range [startTimestepIndex, endTimestepIndex].\n\n"
    print "\n~~~~OR~~~~\n\n"
    print "\nUsage: computePressureMeans.py <startTimestepIndex> <endTimestepIndex> <pressureFileName>\n\n\n"
    sys.exit(0)

try:
    with open('solver.inp','r') as solverInpFile:
        for lineFromFile in solverInpFile:
            if 'list of output surfaces' in lineFromFile.lower():
                # Split on colon, then split the part after the colon on whitespace (to retrieve the list of output surface indices)
                outputSurfacesList = lineFromFile.split(':')[1].split()
except IOError:
    prettyPrintingAvailable = False
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

print "\nMake sure you've generated the PressHist.gplot file using write_rcr_data.gpi before running this script.\n"
print "\nI'm assuming the model is in mm, so scaling the pressures to mmHg by dividing by 133.32.\n"

if (len(sys.argv) == 3):
    file = numpy.loadtxt("PressHist.gplot")
else:
    file =  numpy.loadtxt(sys.argv[3])

pressureScaling = 133.32
userRequestedDataRange = file[startTimestepIndex:endTimestepIndex,:]/pressureScaling

meansArray = numpy.mean(userRequestedDataRange,axis=0)
meansArray = numpy.delete(meansArray, 0) # remove the timestep index mean from the start (it's useless data)
meansArray = numpy.transpose(meansArray)

peakPressureArray = numpy.max(userRequestedDataRange, axis=0)
peakPressureArray = numpy.delete(peakPressureArray, 0)
peakPressureArray = numpy.transpose(peakPressureArray)

minimumPressureArray = numpy.min(userRequestedDataRange, axis=0)
minimumPressureArray = numpy.delete(minimumPressureArray, 0)
minimumPressureArray = numpy.transpose(minimumPressureArray)

pulsePressureArray = peakPressureArray - minimumPressureArray


if prettyPrintingAvailable:
    stackedData = numpy.transpose(numpy.row_stack((meansArray, pulsePressureArray,  peakPressureArray, minimumPressureArray)))
    dataFrame = pandas.DataFrame(stackedData, index=orderedDataLabels, columns=['Pressure Means', 'Pulse Pressure', 'Systolic Pressure', 'Diastolic Pressure'])
    print dataFrame
else:
    print "Means of each column of PressHist (including the time index column!) :"
    print "----"
    print meansArray

print "\n"
