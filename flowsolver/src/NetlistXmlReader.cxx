#include "NetlistXmlReader.hxx"
#include <boost/foreach.hpp>
#include "indexShifters.hxx"
#include <boost/algorithm/string.hpp>

NetlistXmlReader* NetlistXmlReader::msp_instance = 0;

// Netlist data is read from netlist_surfaces.xml and stored in a boost property tree.
// The main task of this class is to extract this data into useable internal data structures,
// store it, and provide it to clients on demand.
void NetlistXmlReader::parseReadData()
{
	gatherNodeIndicesAt3DInterface();
	gatherComponentControllerNames();
	gatherNodalControllerNames();
	countComponentsForEachCircuit();
	countPressureNodesForEachCircuit();
	countPrescribedPressureNodesForEachCircuit();
	countPrescribedFlowsForEachCircuit();
	readCircuitStructure();
	readPrescribedFlowComponents();
	readInitialPressures();
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

void NetlistXmlReader::gatherNodalControllerNames()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));

		// Some containers to fill
		std::map<int,parameter_controller_t> nodalControlTypesForThisSurface;
		std::vector<std::pair<int,std::string>> userDefinedNodeControllersAndPythonNamesForThisSurface;

		// Loop the nodes of this circuit:
		for (auto node : circuit.second.get_child("nodes"))
		{
			int nodeIndex = node.second.get<int>("index");
			boost::optional<boost::property_tree::ptree&> nodalControlSpecification = node.second.get_child_optional("control");
			if (nodalControlSpecification)
			{
				std::string controlType_raw = nodalControlSpecification->get<std::string>("type");
				parameter_controller_t controlType;
				// Work out what sort of control is required for this node:
				if (controlType_raw.compare("customPython") == 0)
				{
					controlType = Controller_CustomPythonNode;
					std::string controlScriptName = nodalControlSpecification->get<std::string>("source");
					// Store the name of the Python script that controls this surface:
					userDefinedNodeControllersAndPythonNamesForThisSurface.push_back(std::make_pair(nodeIndex, controlScriptName));
				}
				else if (controlType_raw.compare("prescribedPeriodicPressure") == 0)
				{
					controlType = Controller_CustomPythonNodePressureFile;
					std::string controlScriptName = nodalControlSpecification->get<std::string>("source");
					// Store the name of the dat file that controls this surface:
					userDefinedNodeControllersAndPythonNamesForThisSurface.push_back(std::make_pair(nodeIndex, controlScriptName));	
				}
				else
				{
					std::stringstream error;
					error << "EE: Unknown node control type found during read of netlist surface " << toOneIndexing(circuitIndex) << ", as indexed by order of appearance in netlist_surfaces.dat." << std::endl;
					error << "This may just indicate a malformed netlist_surfaces.dat" << std::endl;
					throw std::runtime_error(error.str());
				}
				nodalControlTypesForThisSurface.insert( std::make_pair(nodeIndex,controlType) );
			}
		}
		// gather the info for this circuit into the all-circuits data structures:
		m_mapsOfNodalControlTypesForEachSurface.insert(std::make_pair(circuitIndex, nodalControlTypesForThisSurface));
		m_userDefinedNodeControllersAndPythonNames.insert(std::make_pair(circuitIndex, userDefinedNodeControllersAndPythonNamesForThisSurface));
	}
}

void NetlistXmlReader::countComponentsForEachCircuit()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		int numberOfComponentsForThisCircuit = 0;
		for (auto component : circuit.second.get_child("components"))
		{
			numberOfComponentsForThisCircuit++;
		}
		m_numberOfComponents.insert(std::make_pair(circuitIndex, numberOfComponentsForThisCircuit));
	}
}

void NetlistXmlReader::countPressureNodesForEachCircuit()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		int numberOfNodesForThisCircuit = 0;
		for (auto node : circuit.second.get_child("nodes"))
		{
			numberOfNodesForThisCircuit++;
		}
		m_numberOfNodes.insert(std::make_pair(circuitIndex, numberOfNodesForThisCircuit));
	}
}

void NetlistXmlReader::countPrescribedPressureNodesForEachCircuit()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		int numberOfPrescribedPressuresForThisSurface = 0;
		for (auto node : circuit.second.get_child("nodes"))
		{
			boost::optional<std::string> prescribedPressure = node.second.get_optional<std::string>("prescribedPressureType");
			if (prescribedPressure)
			{
				numberOfPrescribedPressuresForThisSurface++;
			}
		}
		m_numberOfPrescribedPressures.insert(std::make_pair(circuitIndex, numberOfPrescribedPressuresForThisSurface));
	}
}
void NetlistXmlReader::countPrescribedFlowsForEachCircuit()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		int numberOfPrescribedFlowsForThisSurface = 0;
		for (auto component : circuit.second.get_child("components"))
		{
			boost::optional<std::string> prescribedFlow = component.second.get_optional<std::string>("prescribedFlowType");
			if (prescribedFlow)
			{
				numberOfPrescribedFlowsForThisSurface++;
			}
		}
		m_numberOfPrescribedFlows.insert(std::make_pair(circuitIndex, numberOfPrescribedFlowsForThisSurface));
	}
}

void NetlistXmlReader::readCircuitStructure()
{
	// Get the netlist-format-style data for each component in the circuit:
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		std::vector<ComponentParameterContainer> tempComponentParameterValues;
		std::vector<circuit_component_t> tempComponentTypes;
		std::vector<int> tempComponentStartNodes;
		std::vector<int> tempComponentEndNodes;

		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		for (auto component : circuit.second.get_child("components"))
		{
			int componentIndex = component.second.get<int>("index");
			std::string componentTypeString = component.second.get<std::string>("type");
			if (boost::iequals(componentTypeString, "resistor"))
			{
				tempComponentTypes.push_back(Component_Resistor);
			}
			else if (boost::iequals(componentTypeString, "capacitor"))
			{
				tempComponentTypes.push_back(Component_Capacitor);
			}
			else if (boost::iequals(componentTypeString, "inductor"))
			{
				tempComponentTypes.push_back(Component_Inductor);
			}
			else if (boost::iequals(componentTypeString, "diode"))
			{
				tempComponentTypes.push_back(Component_Diode);
			}
			else if (boost::iequals(componentTypeString, "volumeTracking"))
			{
				tempComponentTypes.push_back(Component_VolumeTracking);
			}
			else if (boost::iequals(componentTypeString, "volumeTrackingPressureChamber"))
			{
				tempComponentTypes.push_back(Component_VolumeTrackingPressureChamber);
			}
			else
			{
				throw std::runtime_error("ERROR: Unknown netlist component type. This often indicates a malformed netlist_surfaces.dat.\n");
			}

			int startNodeIndex = component.second.get<int>("startNodeIndex");
			tempComponentStartNodes.push_back(startNodeIndex);

			int endNodeIndex = component.second.get<int>("endNodeIndex");
			tempComponentEndNodes.push_back(endNodeIndex);

			ComponentParameterContainer tempContainer;
			tempContainer.setParameter(component.second.get<double>("parameterValue"));
			// If this component has also an initial volume, put it in the ComponentParameterContainer:
			if (tempComponentTypes.back() == Component_VolumeTracking || tempComponentTypes.back() == Component_VolumeTrackingPressureChamber)
			{
				boost::optional<double> initialVolume = component.second.get_optional<double>("initialVolume");
				if (!initialVolume)
				{ 
					throw std::runtime_error("EE: Insufficient parameters given for one of the netlist volume-tracking components.");
				}

				tempContainer.setInitialVolume(*initialVolume);
			}
			tempComponentParameterValues.push_back(tempContainer);
		}
		m_componentTypes.insert(std::make_pair(circuitIndex, tempComponentTypes));
		m_componentStartNodes.insert(std::make_pair(circuitIndex, tempComponentStartNodes));
		m_componentEndNodes.insert(std::make_pair(circuitIndex, tempComponentEndNodes));
		m_componentParameterValues.insert(std::make_pair(circuitIndex, tempComponentParameterValues));
	}
}

void NetlistXmlReader::readPrescribedFlowComponents()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		std::vector<int> tempListOfPrescribedFlows;
		std::vector<double> tempValueOfPrescribedFlows;
		std::vector<circuit_component_flow_prescription_t> tempTypeOfPrescribedFlows;

		for (auto component : circuit.second.get_child("components"))
		{
			boost::optional<std::string> prescribedFlowType = component.second.get_optional<std::string>("prescribedFlowType");
			if (prescribedFlowType)
			{
				tempListOfPrescribedFlows.push_back(component.second.get<int>("index"));
				circuit_component_flow_prescription_t typeOfPrescribedFlow = convertToFlowPrescriptionType(*prescribedFlowType);
				tempTypeOfPrescribedFlows.push_back(typeOfPrescribedFlow);
				if (typeOfPrescribedFlow == Flow_Fixed)
				{
					tempValueOfPrescribedFlows.push_back(component.second.get<int>("prescribedFlowValue"));
				}
				else if (typeOfPrescribedFlow == Flow_3DInterface)
				{
					// just a dummy value which should never be read (the flow will be passed by the 3D interface here):
					tempValueOfPrescribedFlows.push_back(-1.0);	
				}
			}
		}
		m_listOfPrescribedFlows.insert(std::make_pair(circuitIndex, tempListOfPrescribedFlows));
		m_valueOfPrescribedFlows.insert(std::make_pair(circuitIndex, tempValueOfPrescribedFlows));
		m_typeOfPrescribedFlows.insert(std::make_pair(circuitIndex, tempTypeOfPrescribedFlows));
	}
}

circuit_component_flow_prescription_t NetlistXmlReader::convertToFlowPrescriptionType(const std::string inputFileFlowPrescriptionKeyword) const
{
	circuit_component_flow_prescription_t returnValue;
	if (boost::iequals(inputFileFlowPrescriptionKeyword, "fixed"))
	{
		returnValue = Flow_Fixed;
	}
	else if (boost::iequals(inputFileFlowPrescriptionKeyword, "threeDInterface"))
	{
		returnValue = Flow_3DInterface;
	}
	else
	{
		throw std::runtime_error("ERROR: Unknown netlist component flow prescription. This often indicates a malformed netlist_surfaces.dat.");
	}
	return returnValue;
}

void NetlistXmlReader::readInitialPressures()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		std::map<int,double> tempInitialPressures;

		int numberOfPressureNodesCounter = 0;
		for (auto node : circuit.second.get_child("nodes"))
		{
			numberOfPressureNodesCounter++;
			tempInitialPressures.insert( std::make_pair(node.second.get<int>("index"), node.second.get<double>("initialPressure")) );
		}
		m_numberOfPressureNodes.insert(std::make_pair(circuitIndex, numberOfPressureNodesCounter));
		m_initialPressures.insert(std::make_pair(circuitIndex, tempInitialPressures));
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

const std::map<int, std::map<int,parameter_controller_t>>& NetlistXmlReader::getMapsOfComponentControlTypesForEachSurface() const
{
	return m_mapsOfComponentControlTypesForEachSurface;
}

const std::vector<std::pair<int,std::string>>& NetlistXmlReader::getUserDefinedNodeControllersAndPythonNames(const int surfaceIndex) const
{
	return m_userDefinedNodeControllersAndPythonNames.at(surfaceIndex);
}

const std::map<int, std::map<int,parameter_controller_t>>& NetlistXmlReader::getMapsOfNodalControlTypesForEachSurface() const
{
	return m_mapsOfNodalControlTypesForEachSurface;
}

const std::map<int,int>& NetlistXmlReader::getNumberOfComponents() const
{
	return m_numberOfComponents;
}

const std::map<int, int>& NetlistXmlReader::getNumberOfPressureNodes() const
{
	return m_numberOfNodes;
}

const std::map<int, int>& NetlistXmlReader::getNumberOfPrescribedPressures() const
{
	return m_numberOfPrescribedPressures;
}

const std::map<int, int>& NetlistXmlReader::getNumberOfPrescribedFlows() const
{
	return m_numberOfPrescribedFlows;
}

const std::map<int, std::vector<circuit_component_t>>& NetlistXmlReader::getComponentTypes() const
{
	return m_componentTypes;
}

const std::map<int, std::vector<int>>& NetlistXmlReader::getComponentStartNodes() const
{
	return m_componentStartNodes;
}

const std::map<int, std::vector<int>>& NetlistXmlReader::getComponentEndNodes() const
{
	return m_componentEndNodes;
}

const std::vector<double> NetlistXmlReader::getComponentParameterValues(const int indexOfRequestedNetlistLPNDataInInputFile) const
{
	// We extract the component parameter values from their containers and place them
	// in a vector in the order in which they appear in the netlist_surfaces.dat (or whichever input file they come from)
	std::vector<double> allParametersForThisSurface;

	for (auto parameterContainer = m_componentParameterValues.at(indexOfRequestedNetlistLPNDataInInputFile).begin(); parameterContainer != m_componentParameterValues.at(indexOfRequestedNetlistLPNDataInInputFile).end(); parameterContainer++)
	{
		allParametersForThisSurface.push_back(parameterContainer->getParameter());
	}
	
	return allParametersForThisSurface;
}

double NetlistXmlReader::getComponentInitialVolume(const int indexOfRequestedNetlistLPNDataInInputFile, const int componentIndexWithinNetlist) const
{
	return m_componentParameterValues.at(indexOfRequestedNetlistLPNDataInInputFile).at(componentIndexWithinNetlist).getInitialVolume();
}

const std::map<int, std::vector<int>>& NetlistXmlReader::getListOfPrescribedFlows() const
{
	return m_listOfPrescribedFlows;
}
const std::map<int, std::vector<double>>& NetlistXmlReader::getValueOfPrescribedFlows() const
{
	return m_valueOfPrescribedFlows;
}
const std::map<int, std::vector<circuit_component_flow_prescription_t>>& NetlistXmlReader::getTypeOfPrescribedFlows() const
{
	return m_typeOfPrescribedFlows;
}
const std::map<int, std::map<int,double>>& NetlistXmlReader::getInitialPressures() const
{
	return m_initialPressures;
}