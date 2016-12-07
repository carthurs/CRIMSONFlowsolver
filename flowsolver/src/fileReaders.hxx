#ifndef FILEREADERS_H_ 
#define FILEREADERS_H_

#include <fstream>
#include <vector>
#include <utility>
#include <iostream>
#include <memory>
#include <sstream>
#include <map>
#include <set>
#include <cstdlib>
#include <algorithm>
#include "gtest/gtest_prod.h"
#include "debuggingToolsForCpp.hxx"
#include "datatypesInCpp.hxx"
#include "customCRIMSONContainers.hxx"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
// #include <boost/property_tree/json_parser.hpp>

class AbstractFileReader
{
	friend class testMain;
	friend class testFileReaders;
 	FRIEND_TEST(testMain, checkRCRSimpleShortSimulation);
 	FRIEND_TEST(testMain, checkImpedanceMimickingWindkessel);
 	FRIEND_TEST(testMain, checkRestartWorks_RCRSimpleShortSimulation);
public:
	AbstractFileReader()
	{
		m_fileHasBeenRead = false;
		m_metadataOnNumberOfLinesInFileAvailable = false;
		m_hasNumberOfColumns = false;
		m_nextColumnReadLocation = 0;
		m_nextRowReadLoacation = 0;
		mp_currentLineSplitBySpaces = new std::vector<std::string>;
		mp_file = new std::ifstream();
		m_skipHeaderLinesAlreadyCalled = false;
	}
	
	void setFileName(std::string fileNameIn)
	{
	    m_fileName = fileNameIn;
		mp_file->open(m_fileName.c_str());
		
		if (mp_file->fail())
		{
			std::stringstream error;
			error << "Failed to open " << m_fileName << " for reading!" << std::endl;
			throw std::runtime_error(error.str());
		}
		mp_file->clear(); // reset error state flags
		mp_file->seekg(0,mp_file->beg); // ensure the file is rewound
	}

	void setNumColumns(int numberOfColumns)
	{
		m_numColumns = numberOfColumns;
		m_hasNumberOfColumns = true;
	}

	void skipHeaderLines(const int numberOfHeaderLinesToSkip);
	
	virtual ~AbstractFileReader()
	{
		mp_file->close();
		delete mp_file;
		delete mp_currentLineSplitBySpaces;
		// delete m_fileName;
	}

	double getReadFileData(int columnIndex, int timestepNumber);
	double getNextDatum();

	void readFileInternalMetadata();

protected:
	std::ifstream* mp_file;
	std::vector<std::string>* mp_currentLineSplitBySpaces;
	std::string m_currentLine;
	std::string m_fileName;
	// The data as a map of timestep index to vector of all data for that timestep,
	// for all relevant surfaces, in the order in which they appear in the file.
	//\todo: this is buggy, because FORTRAN writes three-in-a-row!
	std::map<int,std::vector<double>> m_dataReadFromFile;
	bool readNextLine();

	int m_numColumns;
	bool m_hasNumberOfColumns;
	std::vector<double> m_dataReadFromFile_line;
	bool m_fileHasBeenRead;
	bool readNextLineWithKnownNumberOfColumns();

	int m_expectedNumberOfLinesInFile;
	bool m_metadataOnNumberOfLinesInFileAvailable;
	int m_nextColumnReadLocation;
	int m_nextRowReadLoacation;
private:
	bool m_skipHeaderLinesAlreadyCalled;
};


// This abstract class is desigend for extenntion to classes which read a single file which contains data for multiple boundaries
// It also acts as a container for the data, from which it can be easily copied to the relevant boundary condition classes
class abstractMultipleSurfaceFileReader : public AbstractFileReader
{
public:
	virtual void readAndSplitMultiSurfaceInputFile() = 0;
	abstractMultipleSurfaceFileReader()
	{
	}

	virtual ~abstractMultipleSurfaceFileReader() override
	{
	}
};

class SimpleFileReader : public AbstractFileReader
{
public:
	SimpleFileReader(std::string fileNameIn)
	{
		setFileName(fileNameIn);
	}

	std::string getNextDataSplitBySpacesOrEndOfLine(bool& success)
	{
		success = false;
		// If there's no data left on the current line from the file, get a new line:
		if (mp_currentLineSplitBySpaces->size()==0)
		{
			success = readNextLine();
			if (!success)
			{
				return std::string("fail");
			}
			// Reverse so we can pop off from the end of the vector in the order
			// that the data appears on the line in the file
			std::reverse(mp_currentLineSplitBySpaces->begin(), mp_currentLineSplitBySpaces->end());
		}

		std::string nextSpaceDelimitedString = mp_currentLineSplitBySpaces->back();
		mp_currentLineSplitBySpaces->pop_back();
		
		// strip any stray carriage return, which will exist at the end of the string if solver.inp was written
		// in Windows and is being read on Linux.
		nextSpaceDelimitedString.erase(std::remove(nextSpaceDelimitedString.begin(), nextSpaceDelimitedString.end(), '\r'), nextSpaceDelimitedString.end());

		success = true;
		return nextSpaceDelimitedString;
	}
private:

};

class HistFileReader : public AbstractFileReader
{
public:
	HistFileReader()
	{
		m_nextRowReadLoacation = 1; // indexing is 1-indexed for hist files, so we override this value from the AbstractFileReader constructor.
	}
	void readAndSplitMultiSurfaceRestartFile();
};


class RcrtReader : public abstractMultipleSurfaceFileReader
{
	friend class testFileReaders;
	friend class testMultidom;
public:
	static RcrtReader* Instance()
	{
		if (!instance)
		{
			instance = new RcrtReader();
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

	void readAndSplitMultiSurfaceInputFile() override;
	int getPdmax();
	std::vector<double> getR1();
	std::vector<double> getC();
	std::vector<double> getR2();
	std::vector<std::vector<std::pair<double,double>>> getTimeDataPdist();
	std::vector<int> getNumDataRCR();
private:
	static RcrtReader* instance;
	// Make the constructor private; it's only ever called as a static method
	// via the public Instance().
	RcrtReader()
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
    	// mp_currentLineSplitBySpaces->clear();
    	// m_dataReadFromFile.clear();
    	// m_dataReadFromFile_line.clear();
    }

	int pdmax;
	std::vector<int> numDataRCR;
	std::vector<double> r1;
	std::vector<double> c;
	std::vector<double> r2;
	std::vector<std::vector<std::pair<double,double>>> timeDataPdist;
	int lengthOfTimeDataPdist;
};

class ControlledCoronaryReader : public abstractMultipleSurfaceFileReader
{
public:
	static ControlledCoronaryReader* Instance()
	{
		if (!instance)
		{
			instance = new ControlledCoronaryReader();
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

	void readAndSplitMultiSurfaceInputFile() override;

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
	ControlledCoronaryReader()
	{
	}

	static ControlledCoronaryReader* instance;

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

using boost::property_tree::ptree;

class NetlistReader : public abstractMultipleSurfaceFileReader
{
public:
	static NetlistReader* Instance()
	{
		if (!instance)
		{
			instance = new NetlistReader();
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

	virtual void readAndSplitMultiSurfaceInputFile() override;

	std::vector<std::vector<circuit_component_t>> getComponentTypes();
	std::vector<std::vector<int>> getComponentStartNodes();
	std::vector<std::vector<int>> getComponentEndNodes();
	std::vector<double> getComponentParameterValues(const int indexOfRequestedNetlistLPNDataInInputFile) const;
	double getComponentInitialVolume(const int indexOfRequestedNetlistLPNDataInInputFile, const int componentIndexWithinNetlist) const;
	std::vector<int> getNumberOfComponents();
	std::vector<int> getNumberOfPrescribedPressures();
	std::vector<int> getNumberOfPrescribedFlows();
	std::vector<std::vector<int>> getListOfPrescribedPressures();
	std::vector<std::vector<int>> getListOfPrescribedFlows();
	std::vector<std::vector<double>> getValueOfPrescribedPressures();
	std::vector<std::vector<double>> getValueOfPrescribedFlows();
	std::vector<std::vector<circuit_nodal_pressure_prescription_t>> getTypeOfPrescribedPressures();
	std::vector<std::vector<circuit_component_flow_prescription_t>> getTypeOfPrescribedFlows();
	std::vector<int> getNumberOfPressureNodes();
	std::vector<std::map<int,double>> getInitialPressures();
	virtual std::vector<int> getIndicesOfNodesAt3DInterface();
	
	std::vector<int>& getNumberOfComponentsWithControl();
	std::vector<std::map<int,parameter_controller_t>>& getMapsOfComponentControlTypesForEachSurface();
	std::vector<int>& getNumberOfNodesWithControl();
	std::vector<std::map<int,parameter_controller_t>>& getMapsOfNodalControlTypesForEachSurface();
	virtual int getNumberOfNetlistSurfaces();
	std::map<int,std::string> getUserDefinedComponentControllersAndPythonNames(const int surfaceIndex) const;
	std::map<int,std::string> getUserDefinedNodeControllersAndPythonNames(const int surfaceIndex) const;

	void writeCircuitSpecificationInXmlFormat();

protected:
	int m_indexOfNetlistCurrentlyBeingReadInFile;
	virtual ~NetlistReader() override
	{
	}

	void readCircuitStructure();
	void readPrescribedPressureNodes();
	void readPrescribedPressureValues();
	void readPrescribedPressureTypes();
	void readPrescribedFlowComponents();
	void readPrescribedFlowValues();
	void readPrescribedFlowTypes();
	void readInitialPressures();
	void readControlSystemPrescriptions();

	void generateBasicPropertyTreeForBoundaryConditionCircuits();
	void writePropertyTreeToDisk(const std::string fileName) const;

	ptree m_propertyTreeRepresentationOfCircuitData;

	std::vector<std::vector<circuit_component_t>> m_componentTypes; // the data in here will be the stripped first column of the netlist, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	std::vector<std::vector<int>> m_componentStartNodes;
	std::vector<std::vector<int>> m_componentEndNodes;
	std::vector<std::vector<ComponentParameterContainer>> m_componentParameterValues;

	std::vector<int> m_numberOfComponents;
	std::vector<int> m_numberOfPrescribedPressures;
	std::vector<int> m_numberOfPrescribedFlows;

	std::vector<std::vector<int>> m_listOfPrescribedPressures;
	std::vector<std::vector<int>> m_listOfPrescribedFlows;

	std::vector<std::vector<double>> m_valueOfPrescribedPressures;
	std::vector<std::vector<double>> m_valueOfPrescribedFlows;
	std::vector<std::vector<circuit_nodal_pressure_prescription_t>> m_typeOfPrescribedPressures;
	std::vector<std::vector<circuit_component_flow_prescription_t>> m_typeOfPrescribedFlows;

	std::vector<int> m_numberOfPressureNodes;
	std::vector<std::map<int,double>> m_initialPressures;

	std::vector<int> m_numberOfComponentsWithControl;
	std::vector<std::map<int,parameter_controller_t>> m_mapsOfComponentControlTypesForEachSurface;

	std::vector<int> m_numberOfNodesWithControl;
	std::vector<std::map<int,parameter_controller_t>> m_mapsOfNodalControlTypesForEachSurface;

	std::vector<std::map<int,std::string>> m_userDefinedComponentControllersAndPythonNames;
	std::vector<std::map<int,std::string>> m_userDefinedNodeControllersAndPythonNames;

	int m_numberOfNetlistSurfacesInDatFile;

	NetlistReader()
	{
	}
private:

	static NetlistReader* instance;

	std::vector<int> m_indicesOfNodesAt3DInterface;

};

class NetlistDownstreamCircuitReader : public NetlistReader
{
public:
	static NetlistDownstreamCircuitReader* Instance()
	{
		if (!downstreamReaderInstance)
		{
			downstreamReaderInstance = new NetlistDownstreamCircuitReader();
		}
		return downstreamReaderInstance;
	}

	static void Term()
	{
		if (downstreamReaderInstance!=0)
		{
			delete downstreamReaderInstance;
			downstreamReaderInstance = 0;
		}
	}

	void readAndSplitMultiSurfaceInputFile() override;

	void writeDownstreamCircuitSpecificationInXmlFormat();

	int getNumberOfBoundaryConditionsConnectedTo(const int downstreamCircuitIndex) const;
	std::vector<int> getConnectedCircuitSurfaceIndices(const int downstreamCircuitIndex) const;
	std::vector<int> getLocalBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const;
	std::vector<int> getRemoteBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const;
	std::set<int> getSetOfNodesInBoundaryConditionWhichConnectToDownstreamCircuit(const int boundaryConditionIndex) const; // boundaryConditionIndex here should be as in the solver.inp.

private:
	static NetlistDownstreamCircuitReader* downstreamReaderInstance;
	std::vector<int> getIndicesOfNodesAt3DInterface() override;
	int getNumberOfNetlistSurfaces() override;
	void readBoundaryConditionConnectivity();
	void checkForBadCircuitDesign();
	void addUpstreamConnectivityInfoToPropertyTree();

	std::vector<int> m_numberOfBoundaryConditionsConnectedTo;
	std::vector<std::vector<int>> m_connectedCircuitSurfaceIndices;
	std::vector<std::vector<int>> m_localBoundaryConditionInterfaceNodes;
	std::vector<std::vector<int>> m_remoteBoundaryConditionInterfaceNodes;

	template <typename Type>
	std::vector<Type> intersectVectors(std::vector<Type> vector1, std::vector<Type> vector2)
	{
		std::vector<Type> intersection(std::max(vector1.size(), vector2.size()));
		auto iteratorToEndOfIntersection = std::set_intersection(vector1.begin(), vector1.end(), vector2.begin(), vector2.end(), intersection.begin());

		// shrink intersection to the size of just the intersecting elements (it's usually going to be longer than that before doing this.)
		intersection.resize(iteratorToEndOfIntersection - intersection.begin());

		return intersection;
	}
};

#endif
