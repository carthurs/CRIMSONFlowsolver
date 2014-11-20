#ifndef FILEREADERS_H_ 
#define FILEREADERS_H_

#include <fstream>
#include <vector>
#include <utility>
#include <iostream>
#include <memory>
#include <sstream>
#include <map>
#include <cstdlib>
#include "gtest/gtest_prod.h"
#include "debuggingToolsForCpp.hxx"

class abstractFileReader
{
	friend class testMain;
	friend class testFileReaders;
 	FRIEND_TEST(testMain, checkSimpleShortSimulationWithRCRs);
 	FRIEND_TEST(testMain, checkRestartWorks_simpleShortSimulationWithRCRs);
public:
	abstractFileReader()
	{
		currentLineSplitBySpaces = new std::vector<std::string>;
		std::cout<< "call to abstractFileReader constructor." << std::endl;
	}
	
	void setFileName(std::string fileNameIn)
	{
	    fileName = fileNameIn;
		fileHandle = new std::ifstream(fileName.c_str());
		
		if (fileHandle->fail())
		{
			std::cout << "Failed to open " << fileName << "!" << std::endl;
			std::exit(1);
		}
	}

	void setNumColumns(int numberOfColumns)
	{
		numColumns = numberOfColumns;
	}
	
	virtual ~abstractFileReader()
	{
		fileHandle->close();
		delete fileHandle;
		delete currentLineSplitBySpaces;
		// delete fileName;
	}
protected:
	std::ifstream* fileHandle;
	std::vector<std::string>* currentLineSplitBySpaces;
	std::string currentLine;
	std::string fileName;
	// The data as a map of timestep index to vector of all data for that timestep,
	// for all relevant surfaces, in the order in which they appear in the file.
	//\todo: this is buggy, because FORTRAN writes three-in-a-row!
	std::map<int,std::vector<double>> dataReadFromFile;
	bool readNextLine();

	int numColumns;
	std::vector<double> dataReadFromFile_line;
	bool readNextLineWithKnownNumberOfColumns();
private:
};

class histFileReader : public abstractFileReader
{
public:
	histFileReader()
	{
	}
	void readAndSplitMultiSurfaceRestartFile();
	
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

	virtual ~abstractMultipleSurfaceFileReader()
	{
	}
};


class rcrtReader : public abstractMultipleSurfaceFileReader
{
	friend class testFileReaders;
	friend class testMultidom;
public:
	static rcrtReader* Instance()
	{
		if (!instance)
		{
			instance = new rcrtReader();
		}
		return instance;
	}

	static void Term()
	{
		delete instance;
		instance = 0;
	}

	void readAndSplitMultiSurfaceInputFile();
	int getPdmax();
	std::vector<double> getR1();
	std::vector<double> getC();
	std::vector<double> getR2();
	std::vector<std::vector<std::pair<double,double>>> getTimeDataPdist();
	std::vector<int> getNumDataRCR();
private:
	static rcrtReader* instance;
	// Make the constructor private; it's only ever called as a static method
	// via the public Instance().
	rcrtReader()
	{
	}

	// For testing purposes, to clear the static class out before the next test begins
    // Note that you'll have to make the test class a friend in order to use this..
    // I've made it private on purpose!
    void tearDown()
    {
    	Term();
    	// numDataRCR.clear();
    	// r1.clear();
    	// c.clear();
    	// r2.clear();
    	// timeDataPdist.clear();
    	// currentLineSplitBySpaces->clear();
    	// dataReadFromFile.clear();
    	// dataReadFromFile_line.clear();
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
