#include "fileReaders.hxx"

// Reads a file line, returns a successful-read bool.
// Stores the read data in a member std::vector<std::string>, which is
// the read line from the file, split by spaces.
bool abstractFileReader::readNextLine()
{
	if (fileHandle->fail())
	{
		std::cout << "Failed to open " << fileName << "!" << std::endl;
		std::exit(1);
	}

	bool returnValue = !(std::getline(*fileHandle,currentLine).eof());

	std::stringstream lineSplitter;
	lineSplitter << currentLine;

	std::string substring;

	currentLineSplitBySpaces->clear();

	while(std::getline(lineSplitter,substring,' '))
	{
		currentLineSplitBySpaces->push_back(substring);
	}

	return returnValue;
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
		r1.push_back(std::stod((*currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		c.push_back(std::stod((*currentLineSplitBySpaces)[0].c_str()));
		readNextLine();
		r2.push_back(std::stod((*currentLineSplitBySpaces)[0].c_str()));

		for(int ii=0; ii<numDataRCR.back(); ii++)
		{
			readNextLine();
			tempTimeAndPdistval.first = std::stod((*currentLineSplitBySpaces)[0].c_str());
			tempTimeAndPdistval.second = std::stod((*currentLineSplitBySpaces)[1].c_str());
			tempTimeDataPdist.push_back(tempTimeAndPdistval);
		}
		timeDataPdist.push_back(tempTimeDataPdist);
	}
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