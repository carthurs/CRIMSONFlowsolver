#ifndef NETLISTXMLREADER_H_ 
#define NETLISTXMLREADER_H_

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <vector>
#include <map>
#include <utility>
#include <iostream>
#include "datatypesInCpp.hxx"

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
	const std::map<int, std::map<int,parameter_controller_t>>& getMapsOfComponentControlTypesForEachSurface();	
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

	static NetlistXmlReader* msp_instance;

	std::map<int, int> m_nodeAt3DInterfaceForEachNetlist;

	boost::property_tree::ptree m_netlistDataFromFile;
	std::map<int, std::vector<std::pair<int, std::string>>> m_userDefinedComponentControllersAndPythonNames;
	std::map<int, std::map<int,parameter_controller_t>> m_mapsOfComponentControlTypesForEachSurface;
};

#endif