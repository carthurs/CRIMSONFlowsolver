#!/usr/bin/python
import numpy
import sys

if len(sys.argv) < 5:
    print "Computes flow means in a flow file.\n"
    print "\nUsage: computeFlowMeans.py <startTimestepIndex> <endTimestepIndex> <timestepInSeconds> <fileName>\n"
    print "\nThis will compute the flow mean for each outlet over the timestep index range [startTimestepIndex, endTimestepIndex].\n\n"
    sys.exit(0)

fileName = sys.argv[4]

print "\n\n ===> Means of each column of " + fileName + " FlowHist (including the time index column!):\n"

startTimestepIndex = int(sys.argv[1])
endTimestepIndex = int(sys.argv[2])
timestep = float(sys.argv[3])

file = numpy.loadtxt(fileName)
columnMeans=numpy.mean(file[startTimestepIndex:endTimestepIndex,:],axis=0)
print columnMeans
print "\n  ===> Sum of means of of all except the first column:"
#print columnMeans[2:len(columnMeans)]
print numpy.sum(columnMeans[1:len(columnMeans)])

columnIntegrals=numpy.sum(file[startTimestepIndex:endTimestepIndex,:],axis=0)*timestep
print "\n\n ===> Column integrals:\n\n", columnIntegrals, "\n"
