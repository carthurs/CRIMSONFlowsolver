#!/apps/anaconda/bin/python
import numpy
import matplotlib.pyplot

dataFile=numpy.loadtxt('FlowHist.gplot')

numberOfTimesteps=dataFile.shape[0]
print "Found %d timesteps in file." % numberOfTimesteps

netFlowsAcrossBoundaries = numpy.ndarray(shape=(numberOfTimesteps,1), dtype=float, order='C')

for ii in range(0,numberOfTimesteps):
	netFlowsAcrossBoundaries[ii] = numpy.sum(dataFile[ii,1:-1])

matplotlib.pyplot.plot(numpy.array(range(0,numberOfTimesteps)), netFlowsAcrossBoundaries)
matplotlib.pyplot.draw()

print "sum of flows on final line of gplot file:"
print numpy.sum(dataFile[-1,1:-1])

matplotlib.pyplot.show()
