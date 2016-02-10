#include "NetlistXmlReader.hxx"
#include <boost/foreach.hpp>
#include "indexShifters.hxx"
#include <boost/algorithm/string.hpp>

NetlistXmlReader* NetlistXmlReader::msp_instance = 0;
NetlistDownstreamXmlReader* NetlistDownstreamXmlReader::msp_downstreamReaderInstance = 0;

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
	readPrescribedPressureNodes();
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
		std::map<int,std::string> userDefinedComponentControllersAndPythonNamesForThisSurface;

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
				else if (boost::iequals(controlType_raw, "leftVentricularElastance"))
				{
					controlType = Controller_LeftVentricularElastance;
				}
				else if (controlType_raw.compare("bleedResistance") == 0)
				{
					controlType = Controller_BleedResistance;
				}
				else if (controlType_raw.compare("customPython") == 0)
				{
					controlType = Controller_CustomPythonComponentParameter;
					// Store the name of the Python script that controls this surface:
					std::string controlScriptName = componentControlSpecification->get<std::string>("source");

					userDefinedComponentControllersAndPythonNamesForThisSurface.insert(std::make_pair(componentIndex, controlScriptName));
				}
				else if (controlType_raw.compare("prescribedPeriodicFlow") == 0)
				{
					controlType = Controller_CustomPythonComponentFlowFile;
					// Store the name of the dat file that controls this surface:
					std::string controlScriptName = componentControlSpecification->get<std::string>("source");

					userDefinedComponentControllersAndPythonNamesForThisSurface.insert(std::make_pair(componentIndex, controlScriptName));
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
		std::map<int,std::string> userDefinedNodeControllersAndPythonNamesForThisSurface;

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
					userDefinedNodeControllersAndPythonNamesForThisSurface.insert(std::make_pair(nodeIndex, controlScriptName));
				}
				else if (controlType_raw.compare("prescribedPeriodicPressure") == 0)
				{
					controlType = Controller_CustomPythonNodePressureFile;
					std::string controlScriptName = nodalControlSpecification->get<std::string>("source");
					// Store the name of the dat file that controls this surface:
					userDefinedNodeControllersAndPythonNamesForThisSurface.insert(std::make_pair(nodeIndex, controlScriptName));	
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
					tempValueOfPrescribedFlows.push_back(component.second.get<double>("prescribedFlowValue"));
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

void NetlistXmlReader::readPrescribedPressureNodes()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int circuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));
		std::vector<int> tempListOfPrescribedPressures;
		std::vector<double> tempValueOfPrescribedPressures;
		std::vector<circuit_nodal_pressure_prescription_t> tempTypeOfPrescribedPressures;

		for (auto node : circuit.second.get_child("nodes"))
		{
			boost::optional<std::string> prescribedPressureType = node.second.get_optional<std::string>("prescribedPressureType");
			if (prescribedPressureType)
			{
				tempListOfPrescribedPressures.push_back(node.second.get<int>("index"));
				circuit_nodal_pressure_prescription_t typeOfPrescribedPressure = convertToPressurePrescriptionType(*prescribedPressureType);
				tempTypeOfPrescribedPressures.push_back(typeOfPrescribedPressure);
				if (typeOfPrescribedPressure == Pressure_Fixed)
				{
					tempValueOfPrescribedPressures.push_back(node.second.get<double>("initialPressure"));
				}
				else if (typeOfPrescribedPressure == Pressure_LeftVentricular)
				{
					// this 1.0 really gets used as a LV pressure scaling value before applying it to the node (in the case of Pressure_LeftVentricular only)
					tempValueOfPrescribedPressures.push_back(1.0);	
				}
			}
		}
		m_listOfPrescribedPressures.insert(std::make_pair(circuitIndex, tempListOfPrescribedPressures));
		m_valueOfPrescribedPressures.insert(std::make_pair(circuitIndex, tempValueOfPrescribedPressures));
		m_typeOfPrescribedPressures.insert(std::make_pair(circuitIndex, tempTypeOfPrescribedPressures));
	}
}

circuit_nodal_pressure_prescription_t NetlistXmlReader::convertToPressurePrescriptionType(const std::string inputFilePressurePrescriptionKeyword) const
{
	circuit_nodal_pressure_prescription_t returnValue;
	if (boost::iequals(inputFilePressurePrescriptionKeyword, "fixed"))
	{
		returnValue = Pressure_Fixed;
	}
	else if (boost::iequals(inputFilePressurePrescriptionKeyword, "leftVentricular"))
	{
		returnValue = Pressure_LeftVentricular;
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

const std::map<int,std::string>& NetlistXmlReader::getUserDefinedComponentControllersAndPythonNames(const int surfaceIndex) const
{
	return m_userDefinedComponentControllersAndPythonNames.at(surfaceIndex);
}

const std::map<int, std::map<int,parameter_controller_t>>& NetlistXmlReader::getMapsOfComponentControlTypesForEachSurface() const
{
	return m_mapsOfComponentControlTypesForEachSurface;
}

const std::map<int,std::string>& NetlistXmlReader::getUserDefinedNodeControllersAndPythonNames(const int surfaceIndex) const
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

const std::map<int, std::vector<int>> NetlistXmlReader::getListOfPrescribedPressures() const
{
	return m_listOfPrescribedPressures;
}
const std::map<int, std::vector<double>> NetlistXmlReader::getValueOfPrescribedPressures() const
{
	return m_valueOfPrescribedPressures;
}
const std::map<int, std::vector<circuit_nodal_pressure_prescription_t>> NetlistXmlReader::getTypeOfPrescribedPressures() const
{
	return m_typeOfPrescribedPressures;
}

const std::string NetlistXmlReader::getXmlComponentNameFromComponentType(const circuit_component_t componentType)
{
	std::string xmlNameString;
	switch (componentType)
	{
		case Component_Resistor:
			xmlNameString = "resistor";
			break;
		case Component_Capacitor:
			xmlNameString = "capacitor";
			break;
		case Component_Inductor:
			xmlNameString = "inductor";
			break;
		case Component_Diode:
			xmlNameString = "diode";
			break;
		case Component_MonopolePressureNode:
			xmlNameString = "monopole";
			break;
		case Component_VolumeTracking:
			xmlNameString = "volumeTracking";
			break;
		case Component_VolumeTrackingPressureChamber:
			xmlNameString = "volumeTrackingPressureChamber";
			break;
		default :
			throw std::runtime_error("Unknown circuit component type in getXmlComponentNameFromComponentType");
	}
	return xmlNameString;
}

const std::string NetlistXmlReader::getXmlFlowPrescritpionNameFromFlowPrescriptionType(const circuit_component_flow_prescription_t flowPrescriptionType)
{
	std::string xmlNameString;
	switch(flowPrescriptionType)
	{
		case Flow_NotPrescribed:
			xmlNameString = "notPrescribed";
			break;
		case Flow_MonopoleSoUndefined:
			xmlNameString = "monopole";
			break;
		case Flow_Fixed:
			xmlNameString = "fixed";
			break;
		case Flow_3DInterface:
			xmlNameString = "threeDInterface";
			break;
		default:
			throw std::runtime_error("Unknown flow prescription type in getXmlFlowPrescritpionNameFromFlowPrescriptionType");
	}

	return xmlNameString;
}

const std::string NetlistXmlReader::getXmlControlNameFromControlType(const parameter_controller_t controlType)
{
	std::string xmlNameString;
	switch (controlType)
	{
		case Controller_LeftVentricularElastance:
			xmlNameString = "leftVentricularElastance";
			break;
		case Controller_BleedResistance:
			xmlNameString = "bleedResistance";
			break;
		case Controller_BleedCompliance:
			xmlNameString = "bleedCompliance";
			break;
		case Controller_CustomPythonComponentParameter:
			xmlNameString = "customPython";
			break;
		case Controller_CustomPythonComponentFlowFile:
			xmlNameString = "prescribedPeriodicFlow";
			break;
		case Controller_CustomPythonNode:
			xmlNameString = "customPython";
			break;
		case Controller_CustomPythonNodePressureFile:
			xmlNameString = "prescribedPeriodicPressure";
			break;
		default:
			throw std::runtime_error("Unknown control prescription type in getXmlControlNameFromControlType");
	}
	return xmlNameString;
}

const std::string NetlistXmlReader::getXmlPressurePrescriptionNameFromPressurePrescriptionType(const circuit_nodal_pressure_prescription_t pressurePrescriptionType)
{
	std::string xmlNameString;
	switch (pressurePrescriptionType)
	{
		case Pressure_NotPrescribed:
			xmlNameString = "pressureNotPrescribed";
			break;
		case Pressure_Fixed:
			xmlNameString = "fixed";
			break;
		case Pressure_LeftVentricular:
			xmlNameString = "leftVentricularPressure";
			break;
		case Pressure_3DInterface:
			xmlNameString = "threeDInterface";
			break;
	}
	return xmlNameString;
}

void NetlistDownstreamXmlReader::parseBoundaryConditionConnectivity()
{
	for (auto circuit : m_netlistDataFromFile.get_child("netlistCircuits"))
	{
		int downstreamCircuitIndex = toZeroIndexing(circuit.second.get<int>("circuitIndex"));

		std::vector<int> tempConnectedCircuitSurfaceIndces;
		std::vector<int> tempLocalBoundaryConditionInterfaceNodes;
		std::vector<int> tempRemoteBoundaryConditionInterfaceNodes;

		int numberOfBoundaryConditionsConnectedTo = 0;
		std::set<int> allUpstreamCircuitSurfaceIndices;
		for (auto node : circuit.second.get_child("nodes"))
		{
			boost::optional<boost::property_tree::ptree&> nodalConenctivitySpecification = node.second.get_child_optional("upstreamConnectivity");
			if (nodalConenctivitySpecification)
			{
				for (auto upstreamConnection : node.second.get_child("upstreamConnectivity"))
				{
					// upstream circuit surface index info:
					int upstreamCircuitSurfaceIndex = upstreamConnection.second.get<int>("upstreamCircuitSurfaceIndex");
					allUpstreamCircuitSurfaceIndices.insert(upstreamCircuitSurfaceIndex);
					tempConnectedCircuitSurfaceIndces.push_back(upstreamCircuitSurfaceIndex);

					// local node of the connection
					tempLocalBoundaryConditionInterfaceNodes.push_back(node.second.get<int>("index"));

					// upstream node of the connection
					tempRemoteBoundaryConditionInterfaceNodes.push_back(upstreamConnection.second.get<int>("connectsToUpstreamNode"));
				}
			}
		}
		// m_numberOfPressureNodes.insert(std::make_pair(circuitIndex, numberOfBoundaryConditionsConnectedTo));
		// m_initialPressures.insert(std::make_pair(circuitIndex, tempInitialPressures));
		m_numberOfBoundaryConditionsConnectedTo.insert(std::make_pair(downstreamCircuitIndex, allUpstreamCircuitSurfaceIndices.size()));
		m_connectedCircuitSurfaceIndices.insert(std::make_pair(downstreamCircuitIndex, tempConnectedCircuitSurfaceIndces));
		m_localBoundaryConditionInterfaceNodes.insert(std::make_pair(downstreamCircuitIndex, tempLocalBoundaryConditionInterfaceNodes));
		m_remoteBoundaryConditionInterfaceNodes.insert(std::make_pair(downstreamCircuitIndex, tempRemoteBoundaryConditionInterfaceNodes));
	}
}

int NetlistDownstreamXmlReader::getNumberOfBoundaryConditionsConnectedTo(const int downstreamCircuitIndex) const
{
	return m_numberOfBoundaryConditionsConnectedTo.at(downstreamCircuitIndex);
}

const std::vector<int>& NetlistDownstreamXmlReader::getConnectedCircuitSurfaceIndices(const int downstreamCircuitIndex) const
{
	return m_connectedCircuitSurfaceIndices.at(downstreamCircuitIndex);
}

const std::vector<int>& NetlistDownstreamXmlReader::getLocalBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const
{
	return m_localBoundaryConditionInterfaceNodes.at(downstreamCircuitIndex);
}

const std::vector<int>& NetlistDownstreamXmlReader::getRemoteBoundaryConditionInterfaceNodes(const int downstreamCircuitIndex) const
{
	return m_remoteBoundaryConditionInterfaceNodes.at(downstreamCircuitIndex);
}

const std::set<int> NetlistDownstreamXmlReader::getSetOfNodesInBoundaryConditionWhichConnectToDownstreamCircuit(const int boundaryConditionIndex) const // boundaryConditionIndex here should be as in the solver.inp
{
	std::set<int> nodesInBoundaryConditionWhichConnectToSomeDownstreamCircuit;

	// Find all the downstream circuits which connect to the requested boundary condition (boundaryConditionIndex)
	for (int downstreamCircuitIndex = 0; downstreamCircuitIndex < m_connectedCircuitSurfaceIndices.size(); downstreamCircuitIndex++)
	{
		for (int attachedBoundaryConditionIndex = 0; attachedBoundaryConditionIndex < m_connectedCircuitSurfaceIndices.at(downstreamCircuitIndex).size(); attachedBoundaryConditionIndex++)
		{
			if (m_connectedCircuitSurfaceIndices.at(downstreamCircuitIndex).at(attachedBoundaryConditionIndex) == boundaryConditionIndex)
			{
				// Gather the nodes in this boundary condition which connect to downstream circuits:
				nodesInBoundaryConditionWhichConnectToSomeDownstreamCircuit.insert(m_remoteBoundaryConditionInterfaceNodes.at(downstreamCircuitIndex).at(attachedBoundaryConditionIndex));
			}
		}
	}

	return nodesInBoundaryConditionWhichConnectToSomeDownstreamCircuit;
}
