#ifndef NETLISTXMLREADER_H_ 
#define NETLISTXMLREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <vector>
#include <map>
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
			msp_instance = new NetlistXmlReader();
		}
		return msp_instance;
	}

	static void Terminante()
	{
		if (msp_instance)
		{
			delete msp_instance;
			msp_instance = 0;
		}
	}

	const std::map<int, int>& getIndicesOfNodesAt3DInterface() const;
	const std::vector<std::pair<int,std::string>>& getUserDefinedComponentControllersAndPythonNames(const int surfaceIndex) const;
	const std::map<int, std::map<int,parameter_controller_t>>& getMapsOfComponentControlTypesForEachSurface() const;
	const std::vector<std::pair<int,std::string>>& getUserDefinedNodeControllersAndPythonNames(const int surfaceIndex) const;
	const std::map<int, std::map<int,parameter_controller_t>>& getMapsOfNodalControlTypesForEachSurface() const;
	const std::map<int, int>& getNumberOfComponents() const;
	const std::map<int, int>& getNumberOfPressureNodes() const;
	const std::map<int, int>& getNumberOfPrescribedPressures() const;
	const std::map<int, int>& getNumberOfPrescribedFlows() const;
	const std::map<int, std::vector<circuit_component_t>>& getComponentTypes() const;
	const std::map<int, std::vector<int>>& getComponentStartNodes() const;
	const std::map<int, std::vector<int>>& getComponentEndNodes() const;
	const std::vector<double> getComponentParameterValues(const int indexOfRequestedNetlistLPNDataInInputFile) const;
	const double getComponentInitialVolume(const int indexOfRequestedNetlistLPNDataInInputFile, const int componentIndexWithinNetlist) const;
private:
	NetlistXmlReader()
	{
		try {
			read_xml("netlist_surfaces.xml", m_netlistDataFromFile);
			parseReadData();
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw e;
		}
	}

	void parseReadData();
	void gatherNodeIndicesAt3DInterface();
	void gatherComponentControllerNames();
	void gatherNodalControllerNames();
	void countComponentsForEachCircuit();
	void countPressureNodesForEachCircuit();
	void countPrescribedPressureNodesForEachCircuit();
	void countPrescribedFlowsForEachCircuit();
	void readCircuitStructure();

	static NetlistXmlReader* msp_instance;

	std::map<int, int> m_nodeAt3DInterfaceForEachNetlist;

	boost::property_tree::ptree m_netlistDataFromFile;
	std::map<int, std::vector<std::pair<int, std::string>>> m_userDefinedComponentControllersAndPythonNames;
	std::map<int, std::map<int,parameter_controller_t>> m_mapsOfComponentControlTypesForEachSurface;
	std::map<int, std::vector<std::pair<int,std::string>>> m_userDefinedNodeControllersAndPythonNames;
	std::map<int, std::map<int,parameter_controller_t>> m_mapsOfNodalControlTypesForEachSurface;
	std::map<int, int> m_numberOfComponents;
	std::map<int, int> m_numberOfNodes;
	std::map<int, int> m_numberOfPrescribedPressures;
	std::map<int, int> m_numberOfPrescribedFlows;
	std::map<int, std::vector<circuit_component_t>> m_componentTypes; // the data in here will be the stripped first column of the netlist, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	std::map<int, std::vector<int>> m_componentStartNodes;
	std::map<int, std::vector<int>> m_componentEndNodes;
	std::map<int, std::vector<ComponentParameterContainer>> m_componentParameterValues;
};

#endif