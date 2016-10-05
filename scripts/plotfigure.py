#!/usr/bin/python
import sys
import matplotlib.pyplot
matplotlib.pyplot.ion() # enable interactive plotting (not really necessary if you're running this as a script.. neither are all the draw() calls!)
import numpy
import copy
import time

print ""
print "=============="
print "Full usage: plotfigure.py <datafiley> <ycolumn in file (zero-indexed)> <ylabel> <yscalefactor> <title> <datafilex> <xcolumn in file (zero-indexed) OR 'auto' for consecutive integers on the x axis> <xlabel> <xscalefactor>"
print "Use quotes if you need spaces in any of these arguments."
print ""
print "You can save the plot as a .png from the GUI, if you wish."
print "If you want text annotations, see the commented-out lines at the bottom of this script."
print "=============="
print ""
print ""
commandLineArgs = copy.deepcopy(sys.argv)
if len(commandLineArgs) < 3:
	print "Not enough command line arguments. Please provide at least the first two."
	print ""
	sys.exit(0)
if len(commandLineArgs) > 10:
	print "Too many command line arguments."
	print ""
	sys.exit(0)
if len(commandLineArgs) < 10:
	print ""
	print "=============="
	print "Proceeding with missing arguments. This is okay, you just won't get some labels."
	print "=============="
	print ""
	print ""
	for nullEntryIndex in range(len(commandLineArgs),10):
		commandLineArgs.append('command line option missing')
	# Make xscalefactor and yscalefactor 1.0 instead of "null"
	for scalingEntryIndex in [4,9]:
		if len(sys.argv) < scalingEntryIndex+1:
			commandLineArgs[scalingEntryIndex] = 1.0

# Put the input data into a string-keyed dictionary:
commandLineArgsNamestrings = ['script','datafiley','ycolumn','ylabel','yscalefactor','title','datafilex','xcolumn','xlabel','xscalefactor']
commandLineDataZipped = zip(commandLineArgsNamestrings, commandLineArgs)
commandLineData = dict(commandLineDataZipped)

# Do some type conversions:
commandLineData['yscalefactor'] = float(commandLineData['yscalefactor'])
commandLineData['xscalefactor'] = float(commandLineData['xscalefactor'])

datax=numpy.loadtxt(commandLineData['datafilex'])
if commandLineData['datafilex'] != 'command line option missing':
	datay = numpy.loadtxt(commandLineData['datafiley'])
else:
	datay = datax

if commandLineData['xcolumn'] == 'auto':
	myplot=matplotlib.pyplot.plot(numpy.array(range(1,len(datax)+1))*commandLineData['xscalefactor'],datay[:,commandLineData['ycolumn']]*commandLineData['yscalefactor'])
elif commandLineData['xcolumn'] != 'command line option missing':
	myplot=matplotlib.pyplot.plot(datax[:,commandLineData['xcolumn']]*commandLineData['xscalefactor'],datay[:,commandLineData['ycolumn']]*commandLineData['yscalefactor'])
else:
	myplot=matplotlib.pyplot.plot(datay[:,commandLineData['ycolumn']]*commandLineData['yscalefactor'])

xlabel=matplotlib.pyplot.xlabel(commandLineData['xlabel'])
matplotlib.pyplot.draw()
ylabel=matplotlib.pyplot.ylabel(commandLineData['ylabel'])
xlabel.set_fontsize(32.0)
ylabel.set_fontsize(32.0)
matplotlib.pyplot.draw()
title=matplotlib.pyplot.title(commandLineData['title'])
title.set_fontsize(34.0)
matplotlib.pyplot.draw()
line,=myplot
line.set_linewidth(2.0)
matplotlib.pyplot.draw()
matplotlib.rc('xtick', labelsize=28.0)
matplotlib.pyplot.draw()
myseries,=myplot
myseries.set_linewidth(2.0)
matplotlib.pyplot.draw()
ax=myseries.get_axes()
ax.tick_params(axis='both',labelsize=28.0)
ax.get_xticks()
ax.set_xticks( numpy.array([ 70., 80.,   90.,  100.,  110.,  120.,  130.]))
myseries.get_figure().set_facecolor('white')
ax.spines['top'].set_visible(False)
matplotlib.pyplot.draw()
ax.spines['right'].set_visible(False)
matplotlib.pyplot.draw()
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')


#############################################
#############################################
#				ANNOTATIONS					#
#############################################
#############################################
# line.axes.text(8, 50, 'boxed italics text in data coords', style='italic',bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
# ax.annotate('Electrical stimulus period changed from 800 to 600 ms', xy=(4.15, 70), xytext=(4.15, 50), fontsize=32.0, arrowprops=dict(facecolor='black', shrink=0.1))
#############################################
#############################################
#				ANNOTATIONS					#
#############################################
#############################################

matplotlib.pyplot.draw()

matplotlib.pyplot.ioff()
matplotlib.pyplot.show()
