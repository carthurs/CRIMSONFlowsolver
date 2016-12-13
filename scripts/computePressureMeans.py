#!/apps/anaconda/bin/python
import numpy
import sys

if len(sys.argv) < 3:
    print "Computes pressure mean from PressHist.gplot.\n"
    print "\nUsage: computePressureMeans.py <startTimestepIndex> <endTimestepIndex>\n"
    print "\nThis will compute the pressure mean for each outlet over the timestep index range [startTimestepIndex, endTimestepIndex].\n\n"
    print "\n~~~~OR~~~~\n\n"
    print "\nUsage: computePressureMeans.py <startTimestepIndex> <endTimestepIndex> <pressureFileName>\n\n\n"
    sys.exit(0)

startTimestepIndex = int(sys.argv[1])
endTimestepIndex = int(sys.argv[2])

print "Make sure you've generated the PressHist.gplot file using write_rcr_data.gpi before running this script."
print "Means of each column of PressHist (including the time index column!) :"
print "----"
print "I'm assuming the model is in mm, so scaling the pressures to mmHg by dividing by 133.32."

if (len(sys.argv) == 3):
    file = numpy.loadtxt("PressHist.gplot")
else:
    file =  numpy.loadtxt(sys.argv[3])

pressureScaling = 133.32
print numpy.mean(file[startTimestepIndex:endTimestepIndex,:]/pressureScaling,axis=0)
