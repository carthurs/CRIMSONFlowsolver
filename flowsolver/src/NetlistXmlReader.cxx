#include "NetlistXmlReader.hxx"
#include <boost/foreach.hpp>
#include "indexShifters.hxx"

NetlistXmlReader* NetlistXmlReader::msp_instance = 0;

void NetlistXmlReader::parseReadData()
{
	gatherNodeIndicesAt3DInterface();
	gatherComponentControllerNames();
}

void NetlistXmlReader::gatherNodeIndicesAt3DInterface()
{
	// Search each circuit in the netlist_surfaces.xml file to find which node is at the 3D interface in each circuit.
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		// loop the nodes of this circuit
		for (auto node : circuit.second.get_child("nodes"))
		{
			boost::optional<std::string> isAt3DInterface = node.second.get_optional<std::string>("isAt3DInterface");
			if (isAt3DInterface && (*isAt3DInterface).compare("true") == 0)
			{
				int componentIndex = node.second.get<int>("index");
				m_nodeAt3DInterfaceForEachNetlist.insert(std::make_pair(circuitIndex, componentIndex));
			}
		}
	}
}

void NetlistXmlReader::gatherComponentControllerNames()
{
	// Loop the circuits
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));

		// Some containers to fill:
		std::map<int,parameter_controller_t> componentControlTypesForThisSurface;
		std::vector<std::pair<int,std::string>> userDefinedComponentControllersAndPythonNamesForThisSurface;

		// loop the components of this circuit
		for (auto component : circuit.second.get_child("components"))
		{
			int componentIndex = component.second.get<int>("index");
			// If there's a controller defined for this component, get its details:
			boost::optional<boost::property_tree::ptree&> componentControlSpecification = component.second.get_child_optional("control");
			if (componentControlSpecification)
			{
				std::string controlType_raw = componentControlSpecification->get<std::string>("type");
				// Work out what sort of control is required for this component:
				parameter_controller_t controlType;
				if (controlType_raw.compare("pc") == 0) // proximal coronary resistor control
				{
					throw std::logic_error("EE: This control hasn't been implemented yet!"); // need to do the push_backs below, too.
				}
				else if (controlType_raw.compare("l") == 0)
				{
					controlType = Controller_LeftVentricularElastance;
				}
				else if (controlType_raw.compare("bleed") == 0)
				{
					controlType = Controller_BleedResistance;
				}
				else if (controlType_raw.compare("customPython") == 0)
				{
					controlType = Controller_CustomPythonComponentParameter;
					// Store the name of the Python script that controls this surface:
					std::string controlScriptName = componentControlSpecification->get<std::string>("source");

					userDefinedComponentControllersAndPythonNamesForThisSurface.push_back(std::make_pair(componentIndex, controlScriptName));
				}
				else if (controlType_raw.compare("prescribedPeriodicFlow") == 0)
				{
					controlType = Controller_CustomPythonComponentFlowFile;
					// Store the name of the dat file that controls this surface:
					std::string controlScriptName = componentControlSpecification->get<std::string>("source");

					userDefinedComponentControllersAndPythonNamesForThisSurface.push_back(std::make_pair(componentIndex, controlScriptName));
				}
				else
				{
					std::stringstream error;
					error << "EE: Unknown component control type found during read of netlist surface " << toOneIndexing(circuitIndex) << ", as indexed by order of appearance in netlist_surfaces.dat." << std::endl;
					error << "This may just indicate a malformed netlist_surfaces.dat" << std::endl;
					throw std::runtime_error(error.str());
				}
				componentControlTypesForThisSurface.insert( std::make_pair(componentIndex,controlType) );
			}
		}
		// gather the info for this circuit into the all-circuits data structures:
		m_mapsOfComponentControlTypesForEachSurface.insert(std::make_pair(circuitIndex, componentControlTypesForThisSurface));
		m_userDefinedComponentControllersAndPythonNames.insert(std::make_pair(circuitIndex, userDefinedComponentControllersAndPythonNamesForThisSurface));
	}
}

const std::map<int, int>& NetlistXmlReader::getIndicesOfNodesAt3DInterface() const
{
	return m_nodeAt3DInterfaceForEachNetlist;
}

const std::vector<std::pair<int,std::string>>& NetlistXmlReader::getUserDefinedComponentControllersAndPythonNames(const int surfaceIndex) const
{
	return m_userDefinedComponentControllersAndPythonNames.at(surfaceIndex);
}

const std::map<int, std::map<int,parameter_controller_t>>& NetlistXmlReader::getMapsOfComponentControlTypesForEachSurface()
{
	return m_mapsOfComponentControlTypesForEachSurface;
}

