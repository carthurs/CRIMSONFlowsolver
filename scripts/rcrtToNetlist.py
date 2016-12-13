#!/usr/bin/python

rcrNetlistData = """### Begin first netlist boundary condition model
# Number of Components
3
# Component 1 type (resistor here):
r
# Component 1 details (start-node index, end-node index, associated parameter (resistance for resistors, capacitance for capacitors):
1
2
%s
# Component 2 type (capacitor here):
c
# Component 2 details:
2
3 
%s
# Component 3 type:
r
# Component 3 details:
2
4
%s
# Number of prescribed pressure nodes:
2
# Indices of nodes with prescribed pressures:
3
4
# Prescribed pressure values / scalings (dependent on types, given by next component):
0.0
%s
# Prescribed pressure types (f=fixed to value given in previous line, l=left ventricular pressure, scaled by value given in previous line):
f
f
# Number of prescribed flows:
1
# Indices of components with prescribed flows
1
# Values of prescribed flows (3D interface set to -1; this value is irrelevant and unused by the code):
-1.0e0
# Types of prescribed flows (t=threeD domain interface)
t
# Number of pressure nodes (including everything- 3D interface, zero-pressure points, internal nodes, etc.):
4
# Initial pressures at the pressure nodes:
1 10585.0
2 10585.0
3 0.0e0
4 %s
# Index of node at 3D Interface
1
# Number of components with control
0
# List of components with control, with the type of control
# number of nodes with control
0
# List of nodes with control, with the type of control (l = LV elastance)
"""

def parseSingleRCR(rcrtFile_local):
	numDistalPressureDatapoints = rcrtFile_local.readline()
	proximalResistance = rcrtFile_local.readline()
	compliance = rcrtFile_local.readline()
	distalResistance = rcrtFile_local.readline()
	distalPressure_fullLine = rcrtFile_local.readline()
	distalPressure = distalPressure_fullLine.split()[1]

	# skip the remaining distal pressure prescriptions (yes, this file only works with the first of them if there are multiple ones)
	for i in range(1,int(numDistalPressureDatapoints)):
		rcrtFile_local.readline()

	rcrParameterPackage = (proximalResistance, compliance, distalResistance, distalPressure, distalPressure) # last one twice as it appears in two places in the output file
	# remove any carriage returns or line feeds:
	rcrParameterPackage = tuple([i.strip() for i in rcrParameterPackage])

	return rcrParameterPackage

def countRCRS(rcrtFile_local):
	# save the file pointer, so we can reset it when we're done
	initialFilePointerLocation = rcrtFile_local.tell()
	#rewind the file
	rcrtFile_local.seek(0)

	rcrtFile_local.readline() # ignore the first line of the file	

	numberOfPressurePrescriptionsForThisRCR = rcrtFile_local.readline() # initialise by getting the first RCR's number of distal pressure prescription lines
	numberOfRCRs = 0
	while numberOfPressurePrescriptionsForThisRCR != "":
		numberOfRCRs = numberOfRCRs + 1
		numberOfDistalPressurePrescriptionLinesForThisRCR = numberOfPressurePrescriptionsForThisRCR
		rcrtFile_local.readline() # skip proximal resistance
		rcrtFile_local.readline() # skip compliance
		rcrtFile_local.readline() # skip distal resistance
		for i in range(0,int(numberOfDistalPressurePrescriptionLinesForThisRCR)):
			rcrtFile_local.readline() # skip the pressure prescriptions
		numberOfPressurePrescriptionsForThisRCR = rcrtFile_local.readline()

	# Reset the file pointer:
	rcrtFile_local.seek(initialFilePointerLocation)

	return numberOfRCRs


rcrtFile = open("rcrt.dat", "r")
rcrtFile.readline() # ditch the first line of the file as we don't need it

numberOfRCRs = countRCRS(rcrtFile)

netlistFile = open("netlist_surfaces.dat", "w")

for i in range(0,numberOfRCRs):
	rcrtDataForThisRCR = parseSingleRCR(rcrtFile)#('a','a','a','a')
	print "Detected proximal resistance, compliance, distal resistance and distal pressure for this surface:"
	print rcrtDataForThisRCR
	netlistFile.write(rcrNetlistData % rcrtDataForThisRCR)

netlistFile.close()
rcrtFile.close()

multidomainFile = open("multidomain.dat","w")
multidomainFile.write("#\n0\n#\n0\n")
multidomainFile.close()
