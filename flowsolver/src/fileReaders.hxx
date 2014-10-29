#ifndef FILEREADERS_H_ 
#define FILEREADERS_H_

#include <fstream>
#include <vector>
#include <utility>
#include <iostream>
#include <memory>
#include <sstream>
#include <map>

class abstractFileReader
{
public:
	abstractFileReader()
	{
		currentLineSplitBySpaces = new std::vector<std::string>;
	}
	
	void setFileName(std::string fileNameIn)
	{
	    fileName = fileNameIn;
		fileHandle = new std::ifstream(fileName);
		
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
	std::ifstream* fileHandle;
	std::vector<std::string>* currentLineSplitBySpaces;
	std::string currentLine;
	std::string fileName;
	bool readNextLine();
private:
};

class histFileReader : public abstractFileReader
{
public:
	histFileReader()
	{
	}
	void readAndSplitMultiSurfaceRestartFile();
	// The data as a map of timestep index to vector of all data for that timestep,
	// for all relevant surfaces, in the order in which they appear in the file.
	//\todo: this is buggy, because FORTRAN writes three-in-a-row!
	std::map<int,std::vector<double>> dataReadFromFile;
};

// This abstract class is desigend for extenntion to classes which read a single file which contains data for multiple boundaries
// It also acts as a container for the data, from which it can be easily copied to the relevant boundary condition classes
class abstractMultipleSurfaceFileReader : public abstractFileReader
{
public:
	virtual void readAndSplitMultiSurfaceInputFile() = 0;
	abstractMultipleSurfaceFileReader()
	{
	}
};


class rcrtReader : public abstractMultipleSurfaceFileReader
{
public:
	static rcrtReader* Instance()
	{
		static rcrtReader* instance = new rcrtReader();
		return instance;
	}


	void readAndSplitMultiSurfaceInputFile();
	int getPdmax();
	std::vector<double> getR1();
	std::vector<double> getC();
	std::vector<double> getR2();
	std::vector<std::vector<std::pair<double,double>>> getTimeDataPdist();
	std::vector<int> getNumDataRCR();
private:
	// Make the constructor private; it's only ever called as a static method
	// via the public Instance().
	rcrtReader()
	{
	}

	int pdmax;
	std::vector<int> numDataRCR;
	std::vector<double> r1;
	std::vector<double> c;
	std::vector<double> r2;
	std::vector<std::vector<std::pair<double,double>>> timeDataPdist;
	int lengthOfTimeDataPdist;
};

#endif