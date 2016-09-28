#!/home/carthurs/anaconda2/bin/python

import sys

if len(sys.argv) != 3:
	print "\nScript to extract nodal pressure and velocity from a nodalData.dat.X file."
	print "Usage: extractNodalData.py <inputFileName> <globalNodeIndex>"
	print "Outputs to <inputFileName>.<globalNodeIndex>.out\n"
	sys.exit(1)

fileName = sys.argv[1]
globalNodeIndexToExtract = int(sys.argv[2])

nodalDataInputFile = open(fileName, 'r')

# Skip over the header text:
for lineIndex in range(0,3):
	nodalDataInputFile.readline()

numberOfNodesInThisFile = int(nodalDataInputFile.readline())

print "Nodes in file: ", numberOfNodesInThisFile

# Skip over the comment line:
nodalDataInputFile.readline()

# Find out which node within the file is wanted for extraction:
nodeIndexWithinFileToExtract = -1 # unfound value
for nodeIndexWithinFile in range(0, numberOfNodesInThisFile):
	globalNodeIndex = int(nodalDataInputFile.readline())
	if globalNodeIndex == globalNodeIndexToExtract:
		nodeIndexWithinFileToExtract = nodeIndexWithinFile

if nodeIndexWithinFileToExtract == -1:
	print "Node", nodeIndexWithinFileToExtract, "was not found in this file."
	sys.exit(1)
else:
	print "The output node has internal index", nodeIndexWithinFileToExtract

# Skip over the comment line:
nodalDataInputFile.readline()

allFileDataLines = nodalDataInputFile.read().split('\n')

outputFileName = fileName + "." + str(globalNodeIndexToExtract) + ".out"
with open(outputFileName, 'w') as outputFile:
	for lineIndex in range(0, len(allFileDataLines)-1):
		if (lineIndex - 2*nodeIndexWithinFileToExtract) % (2*numberOfNodesInThisFile) == 0:
			print allFileDataLines[lineIndex] + allFileDataLines[lineIndex+1]
			outputFile.write(str(allFileDataLines[lineIndex] + allFileDataLines[lineIndex+1])+"\n")

nodalDataInputFile.close()