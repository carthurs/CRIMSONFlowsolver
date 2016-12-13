#!/apps/anaconda/bin/python

# This script adds up rows of the FlowHist.gplot (which you must have created first).
# This allows you to check for mass conservation (the sums should be almost zero).

import numpy

fileData = numpy.loadtxt('FlowHist.gplot')

print "Sums: "
print fileData.sum(1)


