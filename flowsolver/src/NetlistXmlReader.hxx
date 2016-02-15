#ifndef NETLISTXMLREADER_H_ 
#define NETLISTXMLREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include "datatypesInCpp.hxx"
#include "customCRIMSONContainers.hxx"

class NetlistXmlReader
{
public:
	static NetlistXmlReader* Instance() 
	{
		if (!msp_instance)
		{
			msp_instance = new NetlistXmlReader("netlist_surfaces.xml");
		}
		return msp_instance;
	}

	static void Term()
	{
		if (msp_instance)
		{
			delete msp_instance;
			msp_instance = 0;
		}
	}

	const std::map<int, int>& getIndicesOfNodesAt3DInterface() const;
	const std::map<int, std::string>& getUserDefinedComponentControllersAndPythonNames(const int surfaceIndex) const;
	const std::map<int, std::map<int,parameter_controller_t>>& getMapsOfComponentControlTypesForEachSurface() const;
	const std::map<int, std::string>& getUserDefinedNodeControllersAndPythonNames(const int surfaceIndex) const;
	const std::map<int, std::map<int,parameter_controller_t>>& getMapsOfNodalControlTypesForEachSurface() const;
	const std::map<int, int>& getNumberOfComponents() const;
	const std::map<int, int>& getNumberOfPressureNodes() const;
	const std::map<int, int>& getNumberOfPrescribedPressures() const;
	const std::map<int, int>& getNumberOfPrescribedFlows() const;
	const std::map<int, std::vector<circuit_component_t>>& getComponentTypes() const;
	const std::map<int, std::vector<int>>& getComponentStartNodes() const;
	const std::map<int, std::vector<int>>& getComponentEndNodes() const;
	const std::vector<double> getComponentParameterValues(const int indexOfRequestedNetlistLPNDataInInputFile) const;
	double getComponentInitialVolume(const int indexOfRequestedNetlistLPNDataInInputFile, const int componentIndexWithinNetlist) const;
	const std::set<int>& getKalmanFilteredComponentIndicesByCircuitIndex(const int circuitIndex) const;

	const std::map<int, std::vector<int>>& getListOfPrescribedFlows() const;
	const std::map<int, std::vector<double>>& getValueOfPrescribedFlows() const;
	const std::map<int, std::vector<circuit_component_flow_prescription_t>>& getTypeOfPrescribedFlows() const;
	const std::map<int, std::map<int,double>>& getInitialPressures() const;
	const std::map<int, std::vector<int>> getListOfPrescribedPressures() const;
	const std::map<int, std::vector<double>> getValueOfPrescribedPressures() const;
	const std::map<int, std::vector<circuit_nodal_pressure_prescription_t>> getTypeOfPrescribedPressures() const;
	static const std::string getXmlComponentNameFromComponentType(const circuit_component_t componentType);
	static const std::string getXmlFlowPrescritpionNameFromFlowPrescriptionType(const circuit_component_flow_prescription_t flowPrescriptionType);
	static const std::string getXmlControlNameFromControlType(const parameter_controller_t controlType);
	static const std::string getXmlPressurePrescriptionNameFromPressurePrescriptionType(const circuit_nodal_pressure_prescription_t pressurePrescriptionType);

	virtual ~NetlistXmlReader() {};

protected:
	NetlistXmlReader(const std::string fileName)
	{
		initialise(fileName);
	}

	void initialise(const std::string fileName)
	{
		try {
			read_xml(fileName, m_netlistDataFromFile);
			parseReadData();
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw e;
		}
	}

	void parseReadData();
	void readKalmanFilteringTags();
	void gatherNodeIndicesAt3DInterface();
	void gatherComponentControllerNames();
	void gatherNodalControllerNames();
	void countComponentsForEachCircuit();
	void countPressureNodesForEachCircuit();
	void countPrescribedPressureNodesForEachCircuit();
	void countPrescribedFlowsForEachCircuit();
	void readCircuitStructure();
	void readPrescribedFlowComponents();
	circuit_component_flow_prescription_t convertToFlowPrescriptionType(const std::string inputFileFlowPrescriptionKeyword) const;
	void readInitialPressures();
	void readPrescribedPressureNodes();
	circuit_nodal_pressure_prescription_t convertToPressurePrescriptionType(const std::string inputFilePressurePrescriptionKeyword) const;

	static NetlistXmlReader* msp_instance;

	std::map<int, int> m_nodeAt3DInterfaceForEachNetlist;

	boost::property_tree::ptree m_netlistDataFromFile;
	std::map<int, std::map<int, std::string>> m_userDefinedComponentControllersAndPythonNames;
	std::map<int, std::map<int,parameter_controller_t>> m_mapsOfComponentControlTypesForEachSurface;
	std::map<int, std::map<int,std::string>> m_userDefinedNodeControllersAndPythonNames;
	std::map<int, std::map<int,parameter_controller_t>> m_mapsOfNodalControlTypesForEachSurface;
	std::map<int, int> m_numberOfComponents;
	std::map<int, int> m_numberOfNodes;
	std::map<int, int> m_numberOfPrescribedPressures;
	std::map<int, int> m_numberOfPrescribedFlows;
	std::map<int, std::vector<circuit_component_t>> m_componentTypes; // the data in here will be the stripped first column of the netlist, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	std::map<int, std::vector<int>> m_componentStartNodes;
	std::map<int, std::vector<int>> m_componentEndNodes;
	std::map<int, std::vector<ComponentParameterContainer>> m_componentParameterValues;
	std::map<int, std::map<int,double>> m_initialPressures;
	std::map<int, int> m_numberOfPressureNodes;
	std::map<int, std::vector<int>> m_listOfPrescribedFlows;
	std::map<int, std::vector<double>> m_valueOfPrescribedFlows;
	std::map<int, std::vector<circuit_component_flow_prescription_t>> m_typeOfPrescribedFlows;
	std::map<int, std::vector<int>> m_listOfPrescribedPressures;
	std::map<int, std::vector<double>> m_valueOfPrescribedPressures;
	std::map<int, std::vector<circuit_nodal_pressure_prescription_t>> m_typeOfPrescribedPressures;
	std::map<int, std::set<int>> m_listOfKalmanFilteredComponents;
};

class NetlistDownstreamXmlReader : public NetlistXmlReader
{
public:

	static NetlistDownstreamXmlReader* Instance() 
	{
		if (!msp_downstreamReaderInstance)
		{
			msp_downstreamReaderInstance = new NetlistDownstreamXmlReader("netlist_closed_loop_downstream.xml");
		}
		return msp_downstreamReaderInstance;
	}

	static void Term()
	{
		if (msp_downstreamReaderInstance)
		{
			delete msp_downstreamReaderInstance;
			msp_downstreamReaderInstance = 0;
		}
	}

	int getNumberOfBoundaryConditionsConnectedTo(const int downstreamCircuitIndex) const;
	const std::vector<int>& getConnectedCircuitSurfaceIndices(const int downstreamCircuitIndex) const;
	const std::vector<int>& getLocalBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const;
	const std::vector<int>& getRemoteBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const;
	const std::set<int> getSetOfNodesInBoundaryConditionWhichConnectToDownstreamCircuit(const int boundaryConditionIndex) const; // boundaryConditionIndex here should be as in the solver.inp.
	void writeDownstreamCircuitSpecificationInXmlFormat() const;
private:
	
	static NetlistDownstreamXmlReader* msp_downstreamReaderInstance;

	NetlistDownstreamXmlReader(const std::string fileName)
	: NetlistXmlReader(fileName)
	{
		try {
			initialise_downstream();
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw e;
		}
	}

	void initialise_downstream()
	{
		try {
			parseBoundaryConditionConnectivity();
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw e;
		}
	}
	void parseBoundaryConditionConnectivity();

	std::map<int, int> m_numberOfBoundaryConditionsConnectedTo;
	std::map<int, std::vector<int>> m_connectedCircuitSurfaceIndices;
	std::map<int, std::vector<int>> m_localBoundaryConditionInterfaceNodes;
	std::map<int, std::vector<int>> m_remoteBoundaryConditionInterfaceNodes;
};

#endif