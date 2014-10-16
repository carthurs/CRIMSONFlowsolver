#ifndef FILEREADERS_H_ 
#define FILEREADERS_H_

#include <fstream>
#include <vector>
#include <utility>
#include <iostream>
#include <memory>
#include <sstream>

class abstractFileReader
{
public:
	abstractFileReader(std::string fileNameIn)
	{
		fileName = fileNameIn;
		std::cout << "filename is: " <<fileName << std::endl;
		fileHandle = new std::ifstream(fileName);
		currentLineSplitBySpaces = new std::vector<std::string>;
		if (fileHandle->fail())
		{
			std::cout << "Failed to open " << fileName << "!" << std::endl;
			std::exit(1);
		}
	}
	~abstractFileReader()
	{
		delete fileHandle;
		delete currentLineSplitBySpaces;
		// delete fileName;
	}
protected:
	abstractFileReader()
	{
	}
	std::ifstream* fileHandle;
	std::vector<std::string>* currentLineSplitBySpaces;
	std::string currentLine;
	std::string fileName;
	bool readNextLine();
private:
};

// This abstract class is desigend for extenntion to classes which read a single file which contains data for multiple boundaries
// It also acts as a container for the data, from which it can be easily copied to the relevant boundary condition classes
class abstractMultipleSurfaceFileReader : public abstractFileReader
{
public:
	virtual void readAndSplitMultiSurfaceInputFile() = 0;
	abstractMultipleSurfaceFileReader(std::string filename, int surfaces)
	: abstractFileReader(filename)
	{
		numberOfSurfacesOfThisType = surfaces;
	}

	abstractMultipleSurfaceFileReader()
	{
	}
protected:
	int numberOfSurfacesOfThisType;
};


class rcrtReader : public abstractMultipleSurfaceFileReader
{
public:
	static rcrtReader& Instance(std::string filename, int surfaces)
	{
		static rcrtReader *instance = new rcrtReader(filename, surfaces);
		return *instance;
	}

	static rcrtReader& Instance()
	{
		static rcrtReader *instance = new rcrtReader();
		return *instance;
	}

	void readAndSplitMultiSurfaceInputFile();
	int getPdmax();
	std::vector<double> getR1();
	std::vector<double> getC();
	std::vector<double> getR2();
	std::vector<std::vector<std::pair<double,double>>> getTimeDataPdist();
private:
	// Make the constructor private; it's only ever called as a static method
	// via the public Instance().
	rcrtReader(std::string filename, int surfaces)
	: abstractMultipleSurfaceFileReader(filename,surfaces)
	{
	}

	rcrtReader()
	{
	}

	int pdmax;
	std::vector<int> numDataRCR;
	std::vector<double> r1;
	std::vector<double> c;
	std::vector<double> r2;
	std::vector<std::vector<std::pair<double,double>>> timeDataPdist;
};

#endif