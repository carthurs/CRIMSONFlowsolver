#include "fileReaders.hxx"
#include "debuggingToolsForCpp.hxx"
#include "common_c.h"

rcrtReader* rcrtReader::instance = 0;
controlledCoronaryReader* controlledCoronaryReader::instance = 0;
netlistReader* netlistReader::instance = 0;


// Reads a file line, returns a successful-read bool.
// Stores the read data in a member std::vector<std::string>, which is
// the read line from the file, split by spaces.
bool abstractFileReader::readNextLine()
{
	if (fileHandle->fail())
	{
		std::stringstream error;
		error << "File " << fileName << " appears to be malformed.";
		throw std::runtime_error(error.str());
	}

	// Read the next line from the file
	currentLine.clear();
	bool fileNotEnded;
	fileNotEnded = !(std::getline(*fileHandle,currentLine).eof());

	// If the end of the file had not been reached before the above read:
	if (fileNotEnded)
	{
		// See if we have a hash-commented line.
		// If we do, try reading the next line
		while(currentLine.compare(0,1,"#") == int(0))
		{
			currentLine.clear();
			fileNotEnded = !(std::getline(*fileHandle,currentLine).eof());
		}

		if (fileNotEnded)
		{
			// Now we've found a non-commented line, actually read the data from it.
			std::stringstream lineSplitBuffer;
			lineSplitBuffer << currentLine;

			std::string substring;

			currentLineSplitBySpaces->clear();

			while(std::getline(lineSplitBuffer,substring,' '))
			{
				currentLineSplitBySpaces->push_back(substring);
			}
		}
	}

	return fileNotEnded;
}

bool abstractFileReader::readNextLineWithKnownNumberOfColumns()
{

	int index;
	double value;

	dataReadFromFile_line.clear();
	// Get all the data entries from a single line of the file.
	// Return false if the end of file is reached, otherwise, return true.
	// Read data is sequentially placed in the vector dataReadFromFile_line,
	// which is newly cleared on each call to this function.
 	for (int currentColumn=0; currentColumn<numColumns; currentColumn++)
 	{
 		value = 0.0;
 		if (!fileHandle->fail())
 		{
 			*fileHandle >> value;
 		}
 		else
		{
			std::stringstream error;
			error << "File " << fileName << " appears to be malformed.";
			throw std::runtime_error(error.str());
 		}
 		if (fileHandle->eof())
	 	{
	 		if (currentColumn>0)
	 		{
	 			std::stringstream error;
	 			error << "File " << fileName << " terminated early.";
				throw std::runtime_error(error.str());
	 		}
	 		return false;
	 	}
 		dataReadFromFile_line.push_back(value);
 	}

 	// return false case is guarded by an if above
 	return true;

}

// The columnIndex refers to the file columns. It's zero-indexed.
double abstractFileReader::getReadFileData(int columnIndex, int timestepNumber)
{
	if (!fileHasBeenRead)
	{
		std::stringstream error;
		error << "Attempted to access data in file " << fileName << " before it has been read. Terminating.";
		throw std::runtime_error(error.str());
	}

	return ((dataReadFromFile.find(timestepNumber))->second).at(columnIndex);
}

void abstractFileReader::readFileInternalMetadata()
{
	if (!metadataOnNumberOfLinesInFileAvailable)
	{
		// e.g. hist files start with a single integer on the first line, giving the length of the file.
		// We read this first, as a special case
		*fileHandle >> expectedNumberOfLinesInFile;
		metadataOnNumberOfLinesInFileAvailable = true;
	}
	else
	{
		std::stringstream error;
		error << "Attempted to read metadata from file " << fileName << " twice. Don't do this. Exiting";
		throw std::runtime_error(error.str());
	}
}

void histFileReader::readAndSplitMultiSurfaceRestartFile()
{
	int lineIndex = 0;
	while(readNextLineWithKnownNumberOfColumns())
	{
		lineIndex++;
		dataReadFromFile.insert(std::pair<int,std::vector<double>> (lineIndex, dataReadFromFile_line));
	}
	if (metadataOnNumberOfLinesInFileAvailable && (lineIndex != expectedNumberOfLinesInFile))
	{
		std::cout << "WARNING: Failed to read as many lines from  " << fileName << " as the integer on its first line suggests it should have!" << std::endl;
	}

	fileHasBeenRead = int(1);
}

// The output objects are std::vectors at the top level, with each vector
// entry corresponding to one surface
void rcrtReader::readAndSplitMultiSurfaceInputFile()
{
	std::pair<double,double> tempTimeAndPdistval;
	std::vector<std::pair<double,double>> tempTimeDataPdist;

	// Get the pdmax (first line of the file)
	readNextLine();
	pdmax = std::atoi((*currentLineSplitBySpaces)[0].c_str());

	// Loop over the rest of the file to get the relevant RCR data for this boundary:
	while(readNextLine())
	{
		tempTimeDataPdist.clear();

		numDataRCR.push_back(atoi((*currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		r1.push_back(atof((*currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		c.push_back(atof((*currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		r2.push_back(atof((*currentLineSplitBySpaces)[0].c_str()));

		for(int ii=0; ii<numDataRCR.back(); ii++)
		{
			readNextLine();
			tempTimeAndPdistval.first = atof((*currentLineSplitBySpaces)[0].c_str());
			tempTimeAndPdistval.second = atof((*currentLineSplitBySpaces)[1].c_str());
			tempTimeDataPdist.push_back(tempTimeAndPdistval);
		}
		timeDataPdist.push_back(tempTimeDataPdist);
	}
	fileHasBeenRead = int(1);
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
		resistanceNearAorta.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));
		
		readNextLine();
		midResistance.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));
		
		readNextLine();
		distalResistance.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		complianceNearAorta.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		intramyocardialCompliance.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		minimumAllowedResistance.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		maximumAllowedResistance.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		perfusionBedMVO2_previous.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		perfusionBedMVO2_current.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		proportionOfMyocardiumPerfusedByThisSurface.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		metabolicFeedbackGain.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		alphaAdrenergicFeedforwardGain.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		betaAdrenergicFeedforwardGain.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		feedbackDamping.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		O2DemandIntegrationWindow.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		capacitorNearAortaTopPressure.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));

		readNextLine();
		intramyocardialCapacitorTopPressure.push_back(atof((*currentLineSplitBySpaces).at(0).c_str()));
	}

	fileHasBeenRead = int(1);
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

void netlistReader::readAndSplitMultiSurfaceInputFile()
{

	std::vector<char> tempComponentTypes;
	std::vector<int> tempComponentStartNodes;
	std::vector<int> tempComponentEndNodes;
	std::vector<double> tempComponentParameterValues;
	std::vector<int> tempListOfPrescribedPressures;
	std::vector<double> tempValueOfPrescribedPressures;
	std::vector<double> tempValueOfPrescribedFlows;
	std::vector<char> tempTypeOfPrescribedPressures;
	std::vector<int> tempListOfPrescribedFlows;
	std::vector<char> tempTypeOfPrescribedFlows;
	std::vector<double> tempInitialPressures;

	while(readNextLine())
	{
		tempComponentTypes.clear();
		tempComponentStartNodes.clear();
		tempComponentEndNodes.clear();
		tempComponentParameterValues.clear();
		tempListOfPrescribedPressures.clear();
		tempValueOfPrescribedPressures.clear();
		tempValueOfPrescribedFlows.clear();
		tempTypeOfPrescribedPressures.clear();
		tempListOfPrescribedFlows.clear();
		tempTypeOfPrescribedFlows.clear();
		tempInitialPressures.clear();

		numberOfComponents.push_back(atoi((*currentLineSplitBySpaces).at(0).c_str()));

		// Get the netlist-format-style data for each component in the circuit:
		for (int componentIndex=0; componentIndex < numberOfComponents.back(); componentIndex++)
		{
			readNextLine();
			tempComponentTypes.push_back((currentLineSplitBySpaces->at(0)).at(0));

			readNextLine();
			tempComponentStartNodes.push_back(atoi((*currentLineSplitBySpaces).at(0).c_str()));

			readNextLine();
			tempComponentEndNodes.push_back(atoi(currentLineSplitBySpaces->at(0).c_str()));

			readNextLine();
			tempComponentParameterValues.push_back(atof(currentLineSplitBySpaces->at(0).c_str()));
		}
		componentTypes.push_back(tempComponentTypes);
		componentStartNodes.push_back(tempComponentStartNodes);
		componentEndNodes.push_back(tempComponentEndNodes);
		componentParameterValues.push_back(tempComponentParameterValues);

		readNextLine();
		numberOfPrescribedPressures.push_back(atoi(currentLineSplitBySpaces->at(0).c_str()));

		for (int prescribedPressureNodeIndex=0; prescribedPressureNodeIndex < numberOfPrescribedPressures.back(); prescribedPressureNodeIndex++)
		{
			readNextLine();
			tempListOfPrescribedPressures.push_back(atoi(currentLineSplitBySpaces->at(0).c_str()));
		}
		listOfPrescribedPressures.push_back(tempListOfPrescribedPressures);

		for (int prescribedPressureNodeIndex=0; prescribedPressureNodeIndex < numberOfPrescribedPressures.back(); prescribedPressureNodeIndex++)
		{
			readNextLine();
			tempValueOfPrescribedPressures.push_back(atof(currentLineSplitBySpaces->at(0).c_str()));
		}
		valueOfPrescribedPressures.push_back(tempValueOfPrescribedPressures);


		for (int prescribedPressureNodeIndex=0; prescribedPressureNodeIndex < numberOfPrescribedPressures.back(); prescribedPressureNodeIndex++)
		{
			readNextLine();
			tempTypeOfPrescribedPressures.push_back(currentLineSplitBySpaces->at(0).at(0));
		}
		typeOfPrescribedPressures.push_back(tempTypeOfPrescribedPressures);
		

		readNextLine();
		numberOfPrescribedFlows.push_back(atoi(currentLineSplitBySpaces->at(0).c_str()));

		for (int prescribedFlowComponentIndex=0; prescribedFlowComponentIndex < numberOfPrescribedFlows.back(); prescribedFlowComponentIndex++)
		{
			readNextLine();
			tempListOfPrescribedFlows.push_back(atoi(currentLineSplitBySpaces->at(0).c_str()));
		}
		listOfPrescribedFlows.push_back(tempListOfPrescribedFlows);

		for (int prescribedFlowComponentIndex=0; prescribedFlowComponentIndex < numberOfPrescribedFlows.back(); prescribedFlowComponentIndex++)
		{
			readNextLine();
			tempValueOfPrescribedFlows.push_back(atof(currentLineSplitBySpaces->at(0).c_str()));
		}
		valueOfPrescribedFlows.push_back(tempValueOfPrescribedFlows);

		for (int prescribedFlowComponentIndex=0; prescribedFlowComponentIndex < numberOfPrescribedFlows.back(); prescribedFlowComponentIndex++)
		{
			readNextLine();
			tempTypeOfPrescribedFlows.push_back(currentLineSplitBySpaces->at(0).at(0));
		}
		typeOfPrescribedFlows.push_back(tempTypeOfPrescribedFlows);

		readNextLine();
		numberOfPressureNodes.push_back(atoi(currentLineSplitBySpaces->at(0).c_str()));
		for (int pressureNode=0; pressureNode < numberOfPressureNodes.back(); pressureNode++)
		{
			readNextLine();
			tempInitialPressures.push_back(atof(currentLineSplitBySpaces->at(0).c_str()));
		}
		initialPressures.push_back(tempInitialPressures);


	}

	fileHasBeenRead = int(1);
}

std::vector<std::vector<char>> netlistReader::getComponentTypes()
{
	return componentTypes;
}
std::vector<std::vector<int>> netlistReader::getComponentStartNodes()
{
	return componentStartNodes;
}
std::vector<std::vector<int>> netlistReader::getComponentEndNodes()
{
	return componentEndNodes;
}
std::vector<std::vector<double>> netlistReader::getComponentParameterValues()
{
	return componentParameterValues;
}
std::vector<int> netlistReader::getNumberOfComponents()
{
	return numberOfComponents;
}
std::vector<int> netlistReader::getNumberOfPrescribedPressures()
{
	return numberOfPrescribedPressures;
}
std::vector<int> netlistReader::getNumberOfPrescribedFlows()
{
	return numberOfPrescribedFlows;
}
std::vector<std::vector<int>> netlistReader::getListOfPrescribedPressures()
{
	return listOfPrescribedPressures;
}
std::vector<std::vector<int>> netlistReader::getListOfPrescribedFlows()
{
	return listOfPrescribedFlows;
}
std::vector<std::vector<double>> netlistReader::getValueOfPrescribedPressures()
{
	return valueOfPrescribedPressures;
}
std::vector<std::vector<double>> netlistReader::getValueOfPrescribedFlows()
{
	return valueOfPrescribedFlows;
}
std::vector<std::vector<char>> netlistReader::getTypeOfPrescribedPressures()
{
	return typeOfPrescribedPressures;
}
std::vector<std::vector<char>> netlistReader::getTypeOfPrescribedFlows()
{
	return typeOfPrescribedFlows;
}
std::vector<int> netlistReader::getNumberOfPressureNodes()
{
	return numberOfPressureNodes;
}
std::vector<std::vector<double>> netlistReader::getInitialPressures()
{
	return initialPressures;
}