#include "fileReaders.hxx"
#include "debuggingToolsForCpp.hxx"
#include "common_c.h"
#include "datatypesInCpp.hxx"
#include <iterator>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include "indexShifters.hxx"
#include "NetlistXmlReader.hxx"

rcrtReader* rcrtReader::instance = 0;
controlledCoronaryReader* controlledCoronaryReader::instance = 0;
NetlistReader* NetlistReader::instance = 0;
NetlistDownstreamCircuitReader* NetlistDownstreamCircuitReader::downstreamReaderInstance = 0;


// Reads a file line, returns a successful-read bool.
// Stores the read data in a member std::vector<std::string>, which is
// the read line from the file, split by spaces.
bool abstractFileReader::readNextLine()
{
	if (mp_file->eof())
	{
		std::stringstream error;
		error << "Reached the end of " << m_fileName << " before all the required data was read.";
		throw std::runtime_error(error.str());
	}

	if (mp_file->fail())
	{
		std::stringstream error;
		error << "File " << m_fileName << " appears to be malformed.";
		throw std::runtime_error(error.str());
	}

	if(!mp_file->is_open())
	{
		std::stringstream error;
		error << "File " << m_fileName << " was not opened properly.";
		throw std::runtime_error(error.str());	
	}

	// Read the next line from the file
	m_currentLine.clear();
	bool fileNotEnded;

	fileNotEnded = !(std::getline(*mp_file,m_currentLine).eof());

	// If the end of the file had not been reached before the above read:
	if (fileNotEnded)
	{
		// See if we have a hash-commented line.
		// If we do, try reading the next line
		while(m_currentLine.compare(0,1,"#") == int(0))
		{
			m_currentLine.clear();
			fileNotEnded = !(std::getline(*mp_file,m_currentLine).eof());
		}

		if (fileNotEnded)
		{
			// Now we've found a non-commented line, actually read the data from it.
			std::stringstream lineSplitBuffer;
			lineSplitBuffer << m_currentLine;

			std::string substring;

			mp_currentLineSplitBySpaces->clear();

			while(std::getline(lineSplitBuffer,substring,' '))
			{
				// ignore whitespace:
				if (substring.empty())
				{
					continue;
				}
				mp_currentLineSplitBySpaces->push_back(substring);
			}
		}
	}

	return fileNotEnded;
}

bool abstractFileReader::readNextLineWithKnownNumberOfColumns()
{
	assert(m_hasNumberOfColumns);
	int index;
	double value;

	m_dataReadFromFile_line.clear();
	// Get all the data entries from a single line of the file.
	// Return false if the end of file is reached, otherwise, return true.
	// Read data is sequentially placed in the vector m_dataReadFromFile_line,
	// which is newly cleared on each call to this function.
 	for (int currentColumn=0; currentColumn<m_numColumns; currentColumn++)
 	{
 		value = 0.0;
 		if (!mp_file->fail())
 		{
 			*mp_file >> value;
 		}
 		else
		{
			std::stringstream error;
			error << "File " << m_fileName << " appears to be malformed.";
			throw std::runtime_error(error.str());
 		}
 		if (mp_file->eof())
	 	{
	 		if (currentColumn>0)
	 		{
	 			std::stringstream error;
	 			error << "File " << m_fileName << " terminated early.";
				throw std::runtime_error(error.str());
	 		}
	 		return false;
	 	}
 		m_dataReadFromFile_line.push_back(value);
 	}

 	// return false case is guarded by an if above
 	if (m_numColumns > 0)
 	{
 		return true;
 	} else {
 		return false; // for the case where the file is empty anyway
 	}

}

double abstractFileReader::getNextDatum()
{
	assert(m_hasNumberOfColumns);

	if (!m_fileHasBeenRead)
	{
		std::stringstream error;
		error << "Attempted to access data in file " << m_fileName << " before it has been read. Terminating.";
		throw std::runtime_error(error.str());
	}

	double returnValue = getReadFileData(m_nextColumnReadLocation, m_nextRowReadLoacation);
	m_nextColumnReadLocation = (m_nextColumnReadLocation + 1) % m_numColumns;
	
	if (m_nextColumnReadLocation == 0) // if the row counter was just reset
	{
		m_nextRowReadLoacation++;
	}

	return returnValue;
}

// The columnIndex refers to the file columns. It's zero-indexed.
double abstractFileReader::getReadFileData(int columnIndex, int timestepNumber)
{
	if (!m_fileHasBeenRead)
	{
		std::stringstream error;
		error << "Attempted to access data in file " << m_fileName << " before it has been read. Terminating.";
		throw std::runtime_error(error.str());
	}

	double returnValue;
	try {
		returnValue = ((m_dataReadFromFile.find(timestepNumber))->second).at(columnIndex);
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw e;
	}

	return returnValue;
}

void abstractFileReader::readFileInternalMetadata()
{
	if (!m_metadataOnNumberOfLinesInFileAvailable)
	{
		// e.g. hist files start with a single integer on the first line, giving the length of the file.
		// We read this first, as a special case
		*mp_file >> m_expectedNumberOfLinesInFile;
		m_metadataOnNumberOfLinesInFileAvailable = true;
	}
	else
	{
		std::stringstream error;
		error << "Attempted to read metadata from file " << m_fileName << " twice. Don't do this. Exiting";
		throw std::runtime_error(error.str());
	}
}

void histFileReader::readAndSplitMultiSurfaceRestartFile()
{
	int lineIndex = 0;
	while(readNextLineWithKnownNumberOfColumns())
	{
		lineIndex++;
		m_dataReadFromFile.insert(std::pair<int,std::vector<double>> (lineIndex, m_dataReadFromFile_line));
	}
	if (m_metadataOnNumberOfLinesInFileAvailable && (lineIndex != m_expectedNumberOfLinesInFile))
	{
		std::cerr << "WARNING: Failed to read as many lines from  " << m_fileName << " as the integer on its first line suggests it should have!" << std::endl;
		std::cerr << "Expected " << m_expectedNumberOfLinesInFile << " but read " << lineIndex << " lines!" << std::endl;
	}

	m_fileHasBeenRead = true;
}

// The output objects are std::vectors at the top level, with each vector
// entry corresponding to one surface
void rcrtReader::readAndSplitMultiSurfaceInputFile()
{
	std::pair<double,double> tempTimeAndPdistval;
	std::vector<std::pair<double,double>> tempTimeDataPdist;

	// Get the pdmax (first line of the file)
	readNextLine();
	pdmax = std::atoi((*mp_currentLineSplitBySpaces)[0].c_str());

	// Loop over the rest of the file to get the relevant RCR data for this boundary:
	while(readNextLine())
	{
		tempTimeDataPdist.clear();

		numDataRCR.push_back(atoi((*mp_currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		r1.push_back(atof((*mp_currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		c.push_back(atof((*mp_currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		r2.push_back(atof((*mp_currentLineSplitBySpaces)[0].c_str()));

		for(int ii=0; ii<numDataRCR.back(); ii++)
		{
			readNextLine();
			tempTimeAndPdistval.first = atof((*mp_currentLineSplitBySpaces)[0].c_str());
			tempTimeAndPdistval.second = atof((*mp_currentLineSplitBySpaces)[1].c_str());
			tempTimeDataPdist.push_back(tempTimeAndPdistval);
		}
		timeDataPdist.push_back(tempTimeDataPdist);
	}
	m_fileHasBeenRead = true;
}

int rcrtReader::getPdmax()
{
	return pdmax;
}

std::vector<double> rcrtReader::getR1()
{
	return r1;
}

std::vector<double> rcrtReader::getC()
{
	return c;
}

std::vector<double> rcrtReader::getR2()
{
	return r2;
}

std::vector<std::vector<std::pair<double,double>>> rcrtReader::getTimeDataPdist()
{
	return timeDataPdist;
}

std::vector<int> rcrtReader::getNumDataRCR()
{
    return numDataRCR;
}

void controlledCoronaryReader::readAndSplitMultiSurfaceInputFile()
{

	// Loop over the rest of the file to get the relevant RCR data for this boundary:
	while(readNextLine())
	{
		resistanceNearAorta.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));
		
		readNextLine();
		midResistance.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));
		
		readNextLine();
		distalResistance.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		complianceNearAorta.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		intramyocardialCompliance.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		minimumAllowedResistance.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		maximumAllowedResistance.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		perfusionBedMVO2_previous.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		perfusionBedMVO2_current.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		proportionOfMyocardiumPerfusedByThisSurface.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		metabolicFeedbackGain.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		alphaAdrenergicFeedforwardGain.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		betaAdrenergicFeedforwardGain.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		feedbackDamping.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		O2DemandIntegrationWindow.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		capacitorNearAortaTopPressure.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		intramyocardialCapacitorTopPressure.push_back(atof((*mp_currentLineSplitBySpaces).at(0).c_str()));
	}

	m_fileHasBeenRead = true;
}

std::vector<double> controlledCoronaryReader::getResistanceNearAorta()
{
	return resistanceNearAorta;
}

std::vector<double> controlledCoronaryReader::getComplianceNearAorta()
{
	return complianceNearAorta;
}

std::vector<double> controlledCoronaryReader::getMidResistance()
{
	return midResistance;
}

std::vector<double> controlledCoronaryReader::getIntramyocardialCompliance()
{
	return intramyocardialCompliance;
}

std::vector<double> controlledCoronaryReader::getDistalResistance()
{
	return distalResistance;
}

std::vector<double> controlledCoronaryReader::getMinimumAllowedResistance()
{
	return minimumAllowedResistance;
}

std::vector<double> controlledCoronaryReader::getMaximumAllowedResistance()
{
	return maximumAllowedResistance;
}

std::vector<double> controlledCoronaryReader::getPerfusionBedMVO2_previous()
{
	return perfusionBedMVO2_previous;
}

std::vector<double> controlledCoronaryReader::getPerfusionBedMVO2_current()
{
	return perfusionBedMVO2_current;
}

std::vector<double> controlledCoronaryReader::getProportionOfMyocardiumPerfusedByThisSurface()
{
	return proportionOfMyocardiumPerfusedByThisSurface;
}

std::vector<double> controlledCoronaryReader::getMetabolicFeedbackGain()
{
	return metabolicFeedbackGain;
}

std::vector<double> controlledCoronaryReader::getAlphaAdrenergicFeedforwardGain()
{
	return alphaAdrenergicFeedforwardGain;
}

std::vector<double> controlledCoronaryReader::getBetaAdrenergicFeedforwardGain()
{
	return betaAdrenergicFeedforwardGain;
}

std::vector<double> controlledCoronaryReader::getFeedbackDamping()
{
	return feedbackDamping;
}

std::vector<double> controlledCoronaryReader::getO2DemandIntegrationWindow()
{
	return O2DemandIntegrationWindow;
}

std::vector<double> controlledCoronaryReader::getCapacitorNearAortaTopPressure()
{
	return capacitorNearAortaTopPressure;
}

std::vector<double> controlledCoronaryReader::getIntramyocardialCapacitorTopPressure()
{
	return intramyocardialCapacitorTopPressure;
}

void NetlistReader::readAndSplitMultiSurfaceInputFile()
{
	m_indexOfNetlistCurrentlyBeingReadInFile = 0;
	while(readNextLine())
	{
		// This is for error reporting purposes, to inform the user of errors in particular netlists during the read.
		m_indexOfNetlistCurrentlyBeingReadInFile++;

		// These member functions are used just to break
		// the file read up into more human-understandable chunks.
		// Don't change their order.
		readCircuitStructure();
		readPrescribedPressureNodes();
		readPrescribedPressureValues();
		readPrescribedPressureTypes();
		readPrescribedFlowComponents();
		readPrescribedFlowValues();
		readPrescribedFlowTypes();
		readInitialPressures();

		// Get the tag for the node at the 3D interface:
		readNextLine();
		m_indicesOfNodesAt3DInterface.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));

		readControlSystemPrescriptions();

	}

	m_numberOfNetlistSurfacesIn_netlist_surfacesdat = m_indexOfNetlistCurrentlyBeingReadInFile;

	m_fileHasBeenRead = true;
}

void NetlistReader::readCircuitStructure()
{
	std::vector<circuit_component_t> tempComponentTypes;
	std::vector<int> tempComponentStartNodes;
	std::vector<int> tempComponentEndNodes;
	
	// Two nested vectors to allow for when componets have multiple parameters, given on the same netlist_surfaces.dat (or related file) and separated by spaces.
	// Multiple parameter info should be documented just below.
	std::vector<ComponentParameterContainer> tempComponentParameterValues;

	m_numberOfComponents.push_back(atoi((*mp_currentLineSplitBySpaces).at(0).c_str()));

	// Get the netlist-format-style data for each component in the circuit:
	for (int componentIndex=0; componentIndex < m_numberOfComponents.back(); componentIndex++)
	{
		readNextLine();
		if (mp_currentLineSplitBySpaces->at(0).compare("r") == 0)
		{
			tempComponentTypes.push_back(Component_Resistor);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("c") == 0)
		{
			tempComponentTypes.push_back(Component_Capacitor);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("i") == 0)
		{
			tempComponentTypes.push_back(Component_Inductor);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("d") == 0)
		{
			tempComponentTypes.push_back(Component_Diode);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("t") == 0)
		{
			tempComponentTypes.push_back(Component_VolumeTracking);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("v") == 0)
		{
			tempComponentTypes.push_back(Component_VolumeTrackingPressureChamber);
		}
		else
		{
			throw std::runtime_error("ERROR: Unknown netlist component type. This often indicates a malformed netlist_surfaces.dat.\n");
		}

		readNextLine();
		tempComponentStartNodes.push_back(atoi((*mp_currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		tempComponentEndNodes.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));

		readNextLine();
		ComponentParameterContainer tempContainer;
		tempContainer.setParameter(atof(mp_currentLineSplitBySpaces->at(0).c_str()));
		// If this component has also an initial volume, put it in the ComponentParameterContainer:
		if (tempComponentTypes.back() == Component_VolumeTracking || tempComponentTypes.back() == Component_VolumeTrackingPressureChamber)
		{
			if (mp_currentLineSplitBySpaces->size() < 2)
			{ 
				throw std::runtime_error("EE: Insufficient parameters given for one of the netlist volume-tracking components.");
			}

			tempContainer.setInitialVolume(atof(mp_currentLineSplitBySpaces->at(1).c_str()));
		}
		tempComponentParameterValues.push_back(tempContainer);
	}
	m_componentTypes.push_back(tempComponentTypes);
	m_componentStartNodes.push_back(tempComponentStartNodes);
	m_componentEndNodes.push_back(tempComponentEndNodes);
	m_componentParameterValues.push_back(tempComponentParameterValues);
}

void NetlistReader::readPrescribedPressureNodes()
{
	std::vector<int> tempListOfPrescribedPressures;
	
	readNextLine();
	m_numberOfPrescribedPressures.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));

	for (int prescribedPressureNodeIndex=0; prescribedPressureNodeIndex < m_numberOfPrescribedPressures.back(); prescribedPressureNodeIndex++)
	{
		readNextLine();
		tempListOfPrescribedPressures.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));
	}
	m_listOfPrescribedPressures.push_back(tempListOfPrescribedPressures);
}

void NetlistReader::readPrescribedPressureValues()
{
	std::vector<double> tempValueOfPrescribedPressures;
	for (int prescribedPressureNodeIndex=0; prescribedPressureNodeIndex < m_numberOfPrescribedPressures.back(); prescribedPressureNodeIndex++)
	{
		readNextLine();
		tempValueOfPrescribedPressures.push_back(atof(mp_currentLineSplitBySpaces->at(0).c_str()));
	}
	m_valueOfPrescribedPressures.push_back(tempValueOfPrescribedPressures);
}

void NetlistReader::readPrescribedPressureTypes()
{
	std::vector<circuit_nodal_pressure_prescription_t> tempTypeOfPrescribedPressures;
	
	for (int prescribedPressureNodeIndex=0; prescribedPressureNodeIndex < m_numberOfPrescribedPressures.back(); prescribedPressureNodeIndex++)
	{
		readNextLine();
		if (mp_currentLineSplitBySpaces->at(0).compare("f") == 0)
		{
			tempTypeOfPrescribedPressures.push_back(Pressure_Fixed);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("l") == 0)
		{
			tempTypeOfPrescribedPressures.push_back(Pressure_LeftVentricular);
		}
		else
		{
			throw std::runtime_error("EE: Unknown netlist nodal pressure prescription. This often indicates a malformed netlist_surfaces.dat.");
		}
	}
	m_typeOfPrescribedPressures.push_back(tempTypeOfPrescribedPressures);
}

void NetlistReader::readPrescribedFlowComponents()
{
	std::vector<int> tempListOfPrescribedFlows;
	readNextLine();
	m_numberOfPrescribedFlows.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));

	for (int prescribedFlowComponentIndex=0; prescribedFlowComponentIndex < m_numberOfPrescribedFlows.back(); prescribedFlowComponentIndex++)
	{
		readNextLine();
		tempListOfPrescribedFlows.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));
	}
	m_listOfPrescribedFlows.push_back(tempListOfPrescribedFlows);
}

void NetlistReader::readPrescribedFlowValues()
{
	std::vector<double> tempValueOfPrescribedFlows;
	for (int prescribedFlowComponentIndex=0; prescribedFlowComponentIndex < m_numberOfPrescribedFlows.back(); prescribedFlowComponentIndex++)
	{
		readNextLine();
		tempValueOfPrescribedFlows.push_back(atof(mp_currentLineSplitBySpaces->at(0).c_str()));
	}
	m_valueOfPrescribedFlows.push_back(tempValueOfPrescribedFlows);
}

void NetlistReader::readPrescribedFlowTypes()
{
	std::vector<circuit_component_flow_prescription_t> tempTypeOfPrescribedFlows;
	for (int prescribedFlowComponentIndex=0; prescribedFlowComponentIndex < m_numberOfPrescribedFlows.back(); prescribedFlowComponentIndex++)
	{
		readNextLine();
		if (mp_currentLineSplitBySpaces->at(0).compare("f") == 0)
		{
			tempTypeOfPrescribedFlows.push_back(Flow_Fixed);
		}
		else if (mp_currentLineSplitBySpaces->at(0).compare("t") == 0)
		{
			tempTypeOfPrescribedFlows.push_back(Flow_3DInterface);
		}
		else
		{
			throw std::runtime_error("ERROR: Unknown netlist component flow prescription. This often indicates a malformed netlist_surfaces.dat.");
		}
	}
	m_typeOfPrescribedFlows.push_back(tempTypeOfPrescribedFlows);
}

void NetlistReader::readInitialPressures()
{
	std::map<int,double> tempInitialPressures;
	
	readNextLine();
	m_numberOfPressureNodes.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));
	for (int pressureNode=0; pressureNode < m_numberOfPressureNodes.back(); pressureNode++)
	{
		readNextLine();
		tempInitialPressures.insert( std::pair<int,double> (atoi(mp_currentLineSplitBySpaces->at(0).c_str()), atof(mp_currentLineSplitBySpaces->at(1).c_str()))  );
	}
	m_initialPressures.push_back(tempInitialPressures);
}

void NetlistReader::readControlSystemPrescriptions()
{
	// Get the number, list and control types of components which have attached control systems:
	readNextLine();
	m_numberOfComponentsWithControl.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));
	std::map<int,parameter_controller_t> componentControlTypesForThisSurface;
	std::map<int,std::string> userDefinedComponentControllersAndPythonNamesForThisSurface;
	for (int controlledComponent=0; controlledComponent<m_numberOfComponentsWithControl.back(); controlledComponent++)
	{
		parameter_controller_t controlType;
		readNextLine();
		int componentIndex = atoi(mp_currentLineSplitBySpaces->at(0).c_str());
		// Work out what sort of control is required for this component:
		if (mp_currentLineSplitBySpaces->at(1).compare("pc") == 0) // proximal coronary resistor control
		{
			throw std::logic_error("EE: This control hasn't been implemented yet!"); // need to do the push_backs below, too.
			// controlType = 
		}
		else if (mp_currentLineSplitBySpaces->at(1).compare("l") == 0)
		{
			controlType = Controller_LeftVentricularElastance;
		}
		else if (mp_currentLineSplitBySpaces->at(1).compare("bleed") == 0)
		{
			controlType = Controller_BleedResistance;
		}
		else if (mp_currentLineSplitBySpaces->at(1).compare("customPython") == 0)
		{
			controlType = Controller_CustomPythonComponentParameter;
			// Store the name of the Python script that controls this surface:
			userDefinedComponentControllersAndPythonNamesForThisSurface.insert(std::make_pair(componentIndex, mp_currentLineSplitBySpaces->at(2)));
		}
		else if (mp_currentLineSplitBySpaces->at(1).compare("prescribedPeriodicFlow") == 0)
		{
			controlType = Controller_CustomPythonComponentFlowFile;
			// Store the name of the dat file that controls this surface:
			userDefinedComponentControllersAndPythonNamesForThisSurface.insert(std::make_pair(componentIndex, mp_currentLineSplitBySpaces->at(2)));	
		}
		else
		{
			std::stringstream error;
			error << "EE: Unknown component control type found during read of netlist surface " << m_indexOfNetlistCurrentlyBeingReadInFile << ", as indexed by order of appearance in netlist_surfaces.dat." << std::endl;
			error << "This may just indicate a malformed netlist_surfaces.dat" << std::endl;
			throw std::runtime_error(error.str());
		}
		componentControlTypesForThisSurface.insert( std::make_pair(componentIndex,controlType) );
	}
	m_mapsOfComponentControlTypesForEachSurface.push_back(componentControlTypesForThisSurface);
	m_userDefinedComponentControllersAndPythonNames.push_back(userDefinedComponentControllersAndPythonNamesForThisSurface);

	// Get the number, list and control types of nodes which have attached control systems:
	readNextLine();
	m_numberOfNodesWithControl.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));
	std::map<int,parameter_controller_t> nodalControlTypesForThisSurface;
	std::map<int,std::string> userDefinedNodeControllersAndPythonNamesForThisSurface;
	for (int controlledNode=0; controlledNode < m_numberOfNodesWithControl.back(); controlledNode++)
	{
		parameter_controller_t controlType;
		readNextLine();
		// Work out what sort of control is required for this component:
		if (mp_currentLineSplitBySpaces->at(1).compare("customPython") == 0)
		{
			controlType = Controller_CustomPythonNode;
			// Store the name of the Python script that controls this surface:
			int nodeIndex = atoi(mp_currentLineSplitBySpaces->at(0).c_str());
			userDefinedNodeControllersAndPythonNamesForThisSurface.insert(std::make_pair(nodeIndex, mp_currentLineSplitBySpaces->at(2)));
		}
		else if (mp_currentLineSplitBySpaces->at(1).compare("prescribedPeriodicPressure") == 0)
		{
			controlType = Controller_CustomPythonNodePressureFile;
			// Store the name of the dat file that controls this surface:
			int nodeIndex = atoi(mp_currentLineSplitBySpaces->at(0).c_str());
			userDefinedNodeControllersAndPythonNamesForThisSurface.insert(std::make_pair(nodeIndex, mp_currentLineSplitBySpaces->at(2)));	
		}
		else
		{
			std::stringstream error;
			error << "EE: Unknown node control type found during read of netlist surface " << m_indexOfNetlistCurrentlyBeingReadInFile << ", as indexed by order of appearance in netlist_surfaces.dat." << std::endl;
			error << "This may just indicate a malformed netlist_surfaces.dat" << std::endl;
			throw std::runtime_error(error.str());
		}
		int nodeIndex = atoi(mp_currentLineSplitBySpaces->at(0).c_str());
		nodalControlTypesForThisSurface.insert( std::make_pair(nodeIndex,controlType) );
	}
	m_mapsOfNodalControlTypesForEachSurface.push_back(nodalControlTypesForThisSurface);
	m_userDefinedNodeControllersAndPythonNames.push_back(userDefinedNodeControllersAndPythonNamesForThisSurface);
}

// This is designed only for use in converting old-format netlist_surfaces.dat into 
// the new netlist_surfaces.xml format
void NetlistReader::writeCircuitSpecificationInXmlFormat() const
{
	using boost::property_tree::ptree;
	ptree pt;


	for (int circuitIndex = 0; circuitIndex < m_numberOfNetlistSurfacesIn_netlist_surfacesdat; circuitIndex++)
	{
		ptree currentCircuit;
		currentCircuit.put("circuitIndex", toOneIndexing(circuitIndex));

		// Do the components:
		for (int componentIndex = 0; componentIndex < m_numberOfComponents.at(circuitIndex); componentIndex++)
		{
			ptree currentComponent;
			currentComponent.put("index", toOneIndexing(componentIndex));

			circuit_component_t componentType = m_componentTypes.at(circuitIndex).at(componentIndex);
			currentComponent.put("type", NetlistXmlReader::getXmlComponentNameFromComponentType(componentType));

			currentComponent.put("startNodeIndex", m_componentStartNodes.at(circuitIndex).at(componentIndex));
			currentComponent.put("endNodeIndex", m_componentEndNodes.at(circuitIndex).at(componentIndex));
			currentComponent.put("parameterValue", m_componentParameterValues.at(circuitIndex).at(componentIndex).getParameter());

			if (componentType == Component_VolumeTracking || componentType == Component_VolumeTrackingPressureChamber)
			{
				double initialVolume = getComponentInitialVolume(circuitIndex, componentIndex);
				currentComponent.put("initialVolume", initialVolume);
			}


			for (int indexAmongstPrescribedFlowComponents = 0; indexAmongstPrescribedFlowComponents < m_numberOfPrescribedFlows.at(circuitIndex); indexAmongstPrescribedFlowComponents++)
			{
				if (m_listOfPrescribedFlows.at(circuitIndex).at(indexAmongstPrescribedFlowComponents) == toOneIndexing(componentIndex))
				{
					currentComponent.put("prescribedFlowType", NetlistXmlReader::getXmlFlowPrescritpionNameFromFlowPrescriptionType(m_typeOfPrescribedFlows.at(circuitIndex).at(indexAmongstPrescribedFlowComponents)));
				}
			}

			// Gather the component control info into the property tree
			for (auto controlledComponentInfo : m_mapsOfComponentControlTypesForEachSurface.at(circuitIndex))
			{
				int controlledComponentIndex = controlledComponentInfo.first;
				if (controlledComponentIndex == toOneIndexing(componentIndex))
				{
					ptree controlInfoForThisComponent;

					parameter_controller_t controlType = controlledComponentInfo.second;
					controlInfoForThisComponent.put("type", NetlistXmlReader::getXmlControlNameFromControlType(controlType));
					if (controlType == Controller_CustomPythonComponentParameter || controlType == Controller_CustomPythonComponentFlowFile)
					{
						std::string controlSourceInfo = getUserDefinedComponentControllersAndPythonNames(circuitIndex).at(toOneIndexing(componentIndex));
						controlInfoForThisComponent.put("source", controlSourceInfo);
					}
					currentComponent.add_child("control", controlInfoForThisComponent);
					break; // the can't be more than one control specification for this component, and we've found it
				}
			}

			currentCircuit.add_child("components.component", currentComponent);
		}

		// Do the nodes:
		for (int nodeIndex = 0; nodeIndex < m_numberOfPressureNodes.at(circuitIndex); nodeIndex++)
		{
			ptree currentNode;
			currentNode.put("index", toOneIndexing(nodeIndex));
			currentNode.put("initialPressure", m_initialPressures.at(circuitIndex).at(toOneIndexing(nodeIndex)));
			if (toOneIndexing(nodeIndex) == m_indicesOfNodesAt3DInterface.at(circuitIndex))
			{
				currentNode.put("isAt3DInterface", "true");
			}

			for (int indexAmongstPrescribedPressureNodes = 0; indexAmongstPrescribedPressureNodes < m_numberOfPrescribedPressures.at(circuitIndex); indexAmongstPrescribedPressureNodes++)
			{
				if (m_listOfPrescribedPressures.at(circuitIndex).at(indexAmongstPrescribedPressureNodes) == toOneIndexing(nodeIndex))
				{
					currentNode.put("prescribedPressureType", NetlistXmlReader::getXmlPressurePrescriptionNameFromPressurePrescriptionType(m_typeOfPrescribedPressures.at(circuitIndex).at(indexAmongstPrescribedPressureNodes)));
					
					// during the switch to xml input files, we've consolidated the fixed pressure node's pressure with the same node's initial pressure
					// since having both didn't make sense. Therefore, overwrite the initialPressure with the prescribedPressure
					// (the user should never have given them different values anyway!)
					currentNode.put("initialPressure", m_valueOfPrescribedPressures.at(circuitIndex).at(indexAmongstPrescribedPressureNodes));
				}
			}
			
			// Gather the nodal control info:
			for (auto controlledNodeInfo : m_mapsOfNodalControlTypesForEachSurface.at(circuitIndex))
			{
				int controlledNodeIndex = controlledNodeInfo.first;
				if (controlledNodeIndex == toOneIndexing(nodeIndex))
				{
					ptree controlInfoForThisNode;

					parameter_controller_t controlType = controlledNodeInfo.second;
					controlInfoForThisNode.put("type", NetlistXmlReader::getXmlControlNameFromControlType(controlType));
					if (controlType == Controller_CustomPythonNode || controlType == Controller_CustomPythonNodePressureFile)
					{
						std::string controlSourceInfo = getUserDefinedNodeControllersAndPythonNames(circuitIndex).at(toOneIndexing(nodeIndex));
						controlInfoForThisNode.put("source", controlSourceInfo);
					}
					currentNode.add_child("control", controlInfoForThisNode);
					break; // the can't be more than one control specification for this node, and we've found it
				}
			}

			currentCircuit.add_child("nodes.node", currentNode);
		}


		pt.add_child("netlistCircuits.circuit", currentCircuit);
	}

	write_xml("testlist_surfaces.xml", pt);
}

std::map<int,std::string> NetlistReader::getUserDefinedComponentControllersAndPythonNames(const int surfaceIndex) const
{
	return m_userDefinedComponentControllersAndPythonNames.at(surfaceIndex);
}

std::map<int,std::string> NetlistReader::getUserDefinedNodeControllersAndPythonNames(const int surfaceIndex) const
{
	return m_userDefinedNodeControllersAndPythonNames.at(surfaceIndex);
}

std::vector<std::vector<circuit_component_t>> NetlistReader::getComponentTypes()
{
	assert(m_fileHasBeenRead);
	return m_componentTypes;
}
std::vector<std::vector<int>> NetlistReader::getComponentStartNodes()
{
	assert(m_fileHasBeenRead);
	return m_componentStartNodes;
}
std::vector<std::vector<int>> NetlistReader::getComponentEndNodes()
{
	assert(m_fileHasBeenRead);
	return m_componentEndNodes;
}

std::vector<double> NetlistReader::getComponentParameterValues(const int indexOfRequestedNetlistLPNDataInInputFile) const
{
	assert(m_fileHasBeenRead);

	// We extract the component parameter values from their containers and place them
	// in a vector in the order in which they appear in the netlist_surfaces.dat (or whichever input file they come from)
	std::vector<double> allParametersForThisSurface;

	for (auto parameterContainer = m_componentParameterValues.at(indexOfRequestedNetlistLPNDataInInputFile).begin(); parameterContainer != m_componentParameterValues.at(indexOfRequestedNetlistLPNDataInInputFile).end(); parameterContainer++)
	{
		allParametersForThisSurface.push_back(parameterContainer->getParameter());
	}
	
	return allParametersForThisSurface;
}

double NetlistReader::getComponentInitialVolume(const int indexOfRequestedNetlistLPNDataInInputFile, const int componentIndexWithinNetlist) const
{
	return m_componentParameterValues.at(indexOfRequestedNetlistLPNDataInInputFile).at(componentIndexWithinNetlist).getInitialVolume();
}

std::vector<int> NetlistReader::getNumberOfComponents()
{
	assert(m_fileHasBeenRead);
	return m_numberOfComponents;
}
std::vector<int> NetlistReader::getNumberOfPrescribedPressures()
{
	assert(m_fileHasBeenRead);
	return m_numberOfPrescribedPressures;
}
std::vector<int> NetlistReader::getNumberOfPrescribedFlows()
{
	assert(m_fileHasBeenRead);
	return m_numberOfPrescribedFlows;
}
std::vector<std::vector<int>> NetlistReader::getListOfPrescribedPressures()
{
	assert(m_fileHasBeenRead);
	return m_listOfPrescribedPressures;
}
std::vector<std::vector<int>> NetlistReader::getListOfPrescribedFlows()
{
	assert(m_fileHasBeenRead);
	return m_listOfPrescribedFlows;
}
std::vector<std::vector<double>> NetlistReader::getValueOfPrescribedPressures()
{
	assert(m_fileHasBeenRead);
	return m_valueOfPrescribedPressures;
}
std::vector<std::vector<double>> NetlistReader::getValueOfPrescribedFlows()
{
	assert(m_fileHasBeenRead);
	return m_valueOfPrescribedFlows;
}
std::vector<std::vector<circuit_nodal_pressure_prescription_t>> NetlistReader::getTypeOfPrescribedPressures()
{
	assert(m_fileHasBeenRead);
	return m_typeOfPrescribedPressures;
}
std::vector<std::vector<circuit_component_flow_prescription_t>> NetlistReader::getTypeOfPrescribedFlows()
{
	assert(m_fileHasBeenRead);
	return m_typeOfPrescribedFlows;
}
std::vector<int> NetlistReader::getNumberOfPressureNodes()
{
	assert(m_fileHasBeenRead);
	return m_numberOfPressureNodes;
}
std::vector<std::map<int,double>> NetlistReader::getInitialPressures()
{
	assert(m_fileHasBeenRead);
	return m_initialPressures;
}
std::vector<int> NetlistReader::getIndicesOfNodesAt3DInterface()
{
	assert(m_fileHasBeenRead);
	return m_indicesOfNodesAt3DInterface;
}
std::vector<int>& NetlistReader::getNumberOfComponentsWithControl()
{
	assert(m_fileHasBeenRead);
	return m_numberOfComponentsWithControl;
}
std::vector<std::map<int,parameter_controller_t>>& NetlistReader::getMapsOfComponentControlTypesForEachSurface()
{
	assert(m_fileHasBeenRead);
	return m_mapsOfComponentControlTypesForEachSurface;
}
std::vector<int>& NetlistReader::getNumberOfNodesWithControl()
{
	assert(m_fileHasBeenRead);
	return m_numberOfNodesWithControl;
}
std::vector<std::map<int,parameter_controller_t>>& NetlistReader::getMapsOfNodalControlTypesForEachSurface()
{
	assert(m_fileHasBeenRead);
	return m_mapsOfNodalControlTypesForEachSurface;
}

int NetlistReader::getNumberOfNetlistSurfaces()
{
	assert(m_fileHasBeenRead);
	return m_numberOfNetlistSurfacesIn_netlist_surfacesdat;
}

void NetlistDownstreamCircuitReader::readAndSplitMultiSurfaceInputFile()
{
	m_indexOfNetlistCurrentlyBeingReadInFile = 0;
	while(readNextLine())
	{
		// This is for error reporting purposes, to inform the user of errors in particular netlists during the read.
		m_indexOfNetlistCurrentlyBeingReadInFile++;

		// These member functions are used just to break
		// the file read up into more human-understandable chunks.
		// Don't change their order.
		readCircuitStructure();
		readPrescribedPressureNodes();
		readPrescribedPressureValues();
		readPrescribedPressureTypes();
		readPrescribedFlowComponents();
		readPrescribedFlowValues();
		readPrescribedFlowTypes();
		readInitialPressures();
		readBoundaryConditionConnectivity();
		readControlSystemPrescriptions();

	}

	checkForBadCircuitDesign();

	m_fileHasBeenRead = true;
}

void NetlistDownstreamCircuitReader::readBoundaryConditionConnectivity()
{
	readNextLine();
	m_numberOfBoundaryConditionsConnectedTo.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));

	{
		std::vector<int> tempConnectedCircuitSurfaceIndces;
		std::vector<int> tempLocalBoundaryConditionInterfaceNodes;
		std::vector<int> tempRemoteBoundaryConditionInterfaceNodes;
		for (int connectedBoundaryCondition = 0; connectedBoundaryCondition < m_numberOfBoundaryConditionsConnectedTo.back(); connectedBoundaryCondition++)
		{
			readNextLine();
			tempConnectedCircuitSurfaceIndces.push_back(atoi(mp_currentLineSplitBySpaces->at(0).c_str()));
			tempLocalBoundaryConditionInterfaceNodes.push_back(atoi(mp_currentLineSplitBySpaces->at(1).c_str()));
			tempRemoteBoundaryConditionInterfaceNodes.push_back(atoi(mp_currentLineSplitBySpaces->at(2).c_str()));
		}
		m_connectedCircuitSurfaceIndices.push_back(tempConnectedCircuitSurfaceIndces);
		m_localBoundaryConditionInterfaceNodes.push_back(tempLocalBoundaryConditionInterfaceNodes);
		m_remoteBoundaryConditionInterfaceNodes.push_back(tempRemoteBoundaryConditionInterfaceNodes);
	}


}

int NetlistDownstreamCircuitReader::getNumberOfBoundaryConditionsConnectedTo(const int downstreamCircuitIndex) const
{
	assert(m_fileHasBeenRead);
	return m_numberOfBoundaryConditionsConnectedTo.at(downstreamCircuitIndex);
}
std::vector<int> NetlistDownstreamCircuitReader::getConnectedCircuitSurfaceIndices(const int downstreamCircuitIndex) const
{
	assert(m_fileHasBeenRead);
	return m_connectedCircuitSurfaceIndices.at(downstreamCircuitIndex);
}
std::vector<int> NetlistDownstreamCircuitReader::getLocalBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const
{
	assert(m_fileHasBeenRead);
	return m_localBoundaryConditionInterfaceNodes.at(downstreamCircuitIndex);
}
std::vector<int> NetlistDownstreamCircuitReader::getRemoteBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const
{
	assert(m_fileHasBeenRead);
	return m_remoteBoundaryConditionInterfaceNodes.at(downstreamCircuitIndex);
}

// boundaryConditionIndex here should be as in the solver.inp.
std::set<int> NetlistDownstreamCircuitReader::getSetOfNodesInBoundaryConditionWhichConnectToDownstreamCircuit(const int boundaryConditionIndex) const
{
	assert(m_fileHasBeenRead);
	std::set<int> nodesInBoundaryConditionWhichConnectToSomeDownstreamCircuit;

	// Find all the downstream circuits which connect to the requested boundary condition (boundaryConditionIndex)
	for (int downstreamCircuitIndex = 0; downstreamCircuitIndex < m_connectedCircuitSurfaceIndices.size(); downstreamCircuitIndex++)
	{
		for (int attachedBoundaryConditionIndex = 0; attachedBoundaryConditionIndex < m_connectedCircuitSurfaceIndices.at(downstreamCircuitIndex).size(); attachedBoundaryConditionIndex++)
		{
			if (m_connectedCircuitSurfaceIndices.at(downstreamCircuitIndex).at(attachedBoundaryConditionIndex) == boundaryConditionIndex)
			{
				// Gather the nodes in this boundary condition which connect to downstream circuits:
				nodesInBoundaryConditionWhichConnectToSomeDownstreamCircuit.insert(m_remoteBoundaryConditionInterfaceNodes.at(downstreamCircuitIndex).at(attachedBoundaryConditionIndex));
			}
		}
	}

	return nodesInBoundaryConditionWhichConnectToSomeDownstreamCircuit;
}

void NetlistDownstreamCircuitReader::checkForBadCircuitDesign()
{
	// Ensure that each boundary condition connects to at most one
	// downstream circuit (i.e. that it is part of a unique closed
	// loop connected component (in the topological sense)).
	//
	// Note that this requirement is not fundamental; it would just
	// require lots of work in the closed loop matrix construction
	// code to make it viable.

	// To test for uniqueness, we assert that the intersection of
	// the lists of boundary conditions connected to the downstream
	// circuits is empty OR that there is a unique closed loop 
	// circuit.
	if (m_connectedCircuitSurfaceIndices.size() != 1)
	{
		// Work out the size of the largest 
		std::size_t maxNumberOfSurfacesConnectedToAnyClosedLoop = 0; // size_t because std::max expects it as the size() operator in the second argument also returns a size_t
		for (auto connectedSurfacesForThisClosedLoop = m_connectedCircuitSurfaceIndices.begin(); connectedSurfacesForThisClosedLoop != m_connectedCircuitSurfaceIndices.end(); connectedSurfacesForThisClosedLoop++)
		{
			maxNumberOfSurfacesConnectedToAnyClosedLoop = std::max(maxNumberOfSurfacesConnectedToAnyClosedLoop, connectedSurfacesForThisClosedLoop->size());
		}

		// Intersect the lists of connected boundary conditions for each closed loop:
		std::vector<int> intersection(maxNumberOfSurfacesConnectedToAnyClosedLoop);
		// Initialise by intersecting for the first two closed loops:
		intersection = intersectVectors(m_connectedCircuitSurfaceIndices.at(0), m_connectedCircuitSurfaceIndices.at(1));
		// Do all the remaining intersections
		// Begin by moving a start-iterator to m_connectedCircuitSurfaceIndices.at(2):
		auto startingIterator = m_connectedCircuitSurfaceIndices.begin();
		std::advance(startingIterator,2);
		for (auto connectedSurfacesForThisClosedLoop = startingIterator; connectedSurfacesForThisClosedLoop != m_connectedCircuitSurfaceIndices.end(); connectedSurfacesForThisClosedLoop++)
		{
			intersection = intersectVectors(intersection,*connectedSurfacesForThisClosedLoop);
		}

		if (intersection.size() > 0)
		{
			throw std::runtime_error("EE: Each boundary condition should connect to at most one closed-loop downstream circuit.");
		}
	}
	else if (m_connectedCircuitSurfaceIndices.size() == 0)
	{
		throw std::runtime_error("EE: Closed loop code has run, but no closed loop circuits were found. This is probably due some disagreement in the input data.");
	}

}

// Disable unwanted methods:
std::vector<int> NetlistDownstreamCircuitReader::getIndicesOfNodesAt3DInterface()
{
	throw std::logic_error("Method getIndicesOfNodesAt3DInterface() should not be called on the NetlistDownstreamCircuitReader, as it has no 3D interface.");
}
int NetlistDownstreamCircuitReader::getNumberOfNetlistSurfaces()
{
	throw std::logic_error("Method getNumberOfNetlistSurfaces() should not be called on the NetlistDownstreamCircuitReader, as it has no surfaces.");
}