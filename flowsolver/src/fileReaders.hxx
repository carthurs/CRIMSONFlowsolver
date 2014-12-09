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
 	FRIEND_TEST(testMain, checkRCRSimpleShortSimulation);
 	FRIEND_TEST(testMain, checkRestartWorks_RCRSimpleShortSimulation);
public:
	abstractFileReader()
	{
		fileHasBeenRead = 0;
		currentLineSplitBySpaces = new std::vector<std::string>;
		fileHandle = new std::ifstream();
	}
	
	void setFileName(std::string fileNameIn)
	{
	    fileName = fileNameIn;
		fileHandle->open(fileName.c_str());
		
		if (fileHandle->fail())
		{
			std::cout << "Failed to open " << fileName << "!" << std::endl;
			throw std::runtime_error("Failed to open file.");
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

	double getReadFileData(int columnIndex, int timestepNumber);

	void readFileInternalMetadata();

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
	int fileHasBeenRead;
	bool readNextLineWithKnownNumberOfColumns();

	int expectedNumberOfLinesInFile;
	bool metadataOnNumberOfLinesInFileAvailable = false;
private:
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


class histFileReader : public abstractFileReader
{
public:
	histFileReader()
	{
	}
	void readAndSplitMultiSurfaceRestartFile();
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
		if (instance!=0)
		{
			delete instance;
			instance = 0;
		}
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

class controlledCoronaryReader : public abstractMultipleSurfaceFileReader
{
public:
	static controlledCoronaryReader* Instance()
	{
		if (!instance)
		{
			instance = new controlledCoronaryReader();
		}
		return instance;
	}

	static void Term()
	{
		if (instance!=0)
		{
			delete instance;
			instance = 0;
		}
	}

	void readAndSplitMultiSurfaceInputFile();

	std::vector<double> getResistanceNearAorta();
	std::vector<double> getComplianceNearAorta();
	std::vector<double> getMidResistance();
	std::vector<double> getIntramyocardialCompliance();
	std::vector<double> getDistalResistance();

	std::vector<double> getMinimumAllowedResistance();
	std::vector<double> getMaximumAllowedResistance();
	std::vector<double> getPerfusionBedMVO2_previous();
	std::vector<double> getPerfusionBedMVO2_current();
	std::vector<double> getProportionOfMyocardiumPerfusedByThisSurface(); // Just this ventricle
	std::vector<double> getMetabolicFeedbackGain();
	std::vector<double> getAlphaAdrenergicFeedforwardGain();
	std::vector<double> getBetaAdrenergicFeedforwardGain();
	std::vector<double> getFeedbackDamping();
	std::vector<double> getO2DemandIntegrationWindow();

	std::vector<double> getCapacitorNearAortaTopPressure();
	std::vector<double> getIntramyocardialCapacitorTopPressure();

private:
	controlledCoronaryReader()
	{
	}

	static controlledCoronaryReader* instance;

	std::vector<double> resistanceNearAorta;
	std::vector<double> complianceNearAorta;
	std::vector<double> midResistance;
	std::vector<double> intramyocardialCompliance;
	std::vector<double> distalResistance;

	std::vector<double> minimumAllowedResistance;
	std::vector<double> maximumAllowedResistance;
	std::vector<double> perfusionBedMVO2_previous;
	std::vector<double> perfusionBedMVO2_current;
	std::vector<double> proportionOfMyocardiumPerfusedByThisSurface; // Just this ventricle
	std::vector<double> metabolicFeedbackGain;
	std::vector<double> alphaAdrenergicFeedforwardGain;
	std::vector<double> betaAdrenergicFeedforwardGain;
	std::vector<double> feedbackDamping;
	std::vector<double> O2DemandIntegrationWindow;

	std::vector<double> capacitorNearAortaTopPressure;
	std::vector<double> intramyocardialCapacitorTopPressure;
};

class netlistReader : public abstractMultipleSurfaceFileReader
{
public:
	static netlistReader* Instance()
	{
		if (!instance)
		{
			instance = new netlistReader();
		}
		return instance;
	}

	static void Term()
	{
		if (instance!=0)
		{
			delete instance;
			instance = 0;
		}
	}

	void readAndSplitMultiSurfaceInputFile();

	std::vector<std::vector<char>> getComponentTypes();
	std::vector<std::vector<int>> getComponentStartNodes();
	std::vector<std::vector<int>> getComponentEndNodes();
	std::vector<std::vector<double>> getComponentParameterValues();
	std::vector<int> getNumberOfComponents();
	std::vector<int> getNumberOfPrescribedPressures();
	std::vector<int> getNumberOfPrescribedFlows();
	std::vector<std::vector<int>> getListOfPrescribedPressures();
	std::vector<std::vector<int>> getListOfPrescribedFlows();
	std::vector<std::vector<int>> getListOfHistoryPressures();
	std::vector<std::vector<int>> getListOfHistoryFlows();
	std::vector<std::vector<double>> getValueOfPrescribedPressures();
	std::vector<std::vector<double>> getValueOfPrescribedFlows();
	std::vector<std::vector<char>> getTypeOfPrescribedPressures();
	std::vector<std::vector<char>> getTypeOfPrescribedFlows();
	std::vector<int> getNumberOfPressureNodes();
	std::vector<std::vector<double>> getInitialPressures();
protected:
private:
	netlistReader()
	{
	}

	static netlistReader* instance;

	std::vector<std::vector<char>> componentTypes; // the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	std::vector<std::vector<int>> componentStartNodes;
	std::vector<std::vector<int>> componentEndNodes;
	std::vector<std::vector<double>> componentParameterValues;

	std::vector<int> numberOfComponents;
	std::vector<int> numberOfPrescribedPressures;
	std::vector<int> numberOfPrescribedFlows;

	std::vector<std::vector<int>> listOfPrescribedPressures;
	std::vector<std::vector<int>> listOfPrescribedFlows;

	std::vector<std::vector<double>> valueOfPrescribedPressures;
	std::vector<std::vector<double>> valueOfPrescribedFlows;
	std::vector<std::vector<char>> typeOfPrescribedPressures;
	std::vector<std::vector<char>> typeOfPrescribedFlows;

	std::vector<int> numberOfPressureNodes;
	std::vector<std::vector<double>> initialPressures;

};

#endif
