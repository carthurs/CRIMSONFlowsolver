#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include <set>
#include <vector>
#include "datatypesInCpp.hxx"
#include "netlistSubcircuit.hxx"

class CircuitComponent
{
public:
	circuit_component_t type;
	CircuitPressureNode startNode;
	CircuitPressureNode endNode;
	double parameterValue; // resistance or compliance or inductance etc.
	int indexInCircuit;
	circuit_component_flow_prescription_t prescribedFlowType
	double valueOfPrescribedFlow;
	CircuitComponent()
	{
		type = Component_Null;
		prescribedFlowType = Flow_Null;
	}
};

class CircuitPressureNode
{
public:
	double valueOfPrescribedPressure;
	double pressure;
	circuit_nodal_pressure_prescription_t prescribedPressureType;
	int indexInCircuit;
	bool nodeIsAt3DInterface;
	CircuitPressureNode()
	{
		prescribedPressureType = Pressure_Null;
		nodeIsAt3DInterface = false;
	}
};

class CircuitData
{
public:
	std::vector<boost::shared_ptr<CircuitComponent>> components;
	// std::vector<circuit_component_t> componentTypes; // the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	// std::vector<int> componentStartNodes;
	// std::vector<int> componentEndNodes;
	// std::vector<double> componentParameterValues;
	// std::vector<int> componentIndices; // This is just a list 1,2,3,..,numberOfComponents
	// std::vector<circuit_nodal_pressure_prescription_t> typeOfPrescribedPressures;
	// std::vector<circuit_component_flow_prescription_t> typeOfPrescribedFlows;
	int numberOfPrescribedPressures;
	int numberOfPrescribedFlows;
	int numberOfPressureNodes;
	int numberOfComponents;
	// std::vector<int> listOfPrescribedPressures;         // input data, listing node indices of prescribed pressures (e.g. LV pressure on a capacitor, venous pressure in open loop, etc.).
	// std::vector<int> listOfPrescribedFlows;             // input data, listing node indices of prescribed flows, listing the component indices with prescribed flows (e.g. 3D outlet flow)
	// std::vector<double> pressuresInSubcircuit;        // Pressure at each LPN node, using the same node indexing as in the netlist
	// std::vector<double> valueOfPrescribedPressures;
	// std::vector<double> valueOfPrescribedFlows;
	std::set<int> setOfPressureNodes; // Utility data structure, containing all the pressure nodes of the circuit, exactly once each. If needed, initialise this manually.
	std::set<int> setOfPrescribedPressureNodes;
	std::set<int> setOfComponentIndices; // Utility data structure containing all the component indices of the circuit, exactly once each. If needed, initialise this manually.
	
	void rebuildCircuitMetadata()
	{
		// This subroutine builds everything in the CircuitData class that isn't the CircuitComponents themselves.
		// It parses the CircuitComponents to do this, so it assumes the CircuitComponents are already fully initialised.

		numberOfPrescribedFlows = 0;

		setOfPressureNodes.clear();
		setOfComponentIndices.clear();
		setOfPrescribedPressureNodes.clear();
		for (auto component = components.begin(); component!=components.end(); component++)
		{
			assert(component->startNode.prescribedPressureType!=Pressure_Null)
			if(component->startNode.prescribedPressureType != Pressure_NotPrescribed)
			{
				setOfPrescribedPressureNodes.insert(component->startNode.indexInCircuit);
			}
			setOfPressureNodes.insert(component->startNode.indexInCircuit);

			assert(component->endNode.prescribedPressureType!=Pressure_Null)
			if(component->endNode.prescribedPressureType != Pressure_NotPrescribed)
			{
				setOfPrescribedPressureNodes.insert(component->endNode.indexInCircuit);
			}
			setOfPressureNodes.insert(component->endNode.indexInCircuit);

			setOfComponentIndices.insert(component->indexInCircuit);

			assert(component->prescribedFlowType != Flow_Null);
			if(component->prescribedFlowType != Flow_NotPrescribed)
			{
				numberOfPrescribedFlows++;
			}
		}

		numberOfPressureNodes = setOfPressureNodes.size();
		numberOfPrescribedPressures = setOfPrescribedPressureNodes.size();
		numberOfComponents = setOfComponentIndices.size();

	}
};

// The particular pattern of open and closed valves in the circuit breaks it into
// varying isolated subcircuits. All possible circuits are enumerated by the
// netlistBoundaryCondition; the struct activeSubcircuits tracks those which
// are currently active due to the current state of all valves.
struct ActiveSubcircuits
{
	std::vector<int> list; // The list of indices of currently active subcircuits
	int indexOfSubcircuitAt3DInterface;
};

class AtomicSubcircuitConnectionManager
{
public:
	AtomicSubcircuitConnectionManager::AtomicSubcircuitConnectionManager(int numberOfDiodes)
	{
		startNodeType.insert(startNodeType.begin(),numberOfDiodes,Node_Null);
		endNodeType.insert(endNodeType.begin(),numberOfDiodes,Node_Null);
		circuitConnectedToStartNode.insert(circuitConnectedToStartNode.begin(),numberOfDiodes,NULL);
		circuitConnectedToEndNode.insert(circuitConnectedToEndNode.begin(),numberOfDiodes,NULL);
	}

	void setDiodeInfo(int diodeIndex, circuit_diode_node_t startNodeType_in, circuit_diode_node_t endNodeType_in, boost::shared_ptr<CircuitData> circuitConnectedToStartNode_in, boost::shared_ptr<CircuitData> circuitConnectedToEndNode_in)
	{
		assert(startNodeType.at(diodeIndex) == Node_Null);
		startNodeType.at(diodeIndex) = startNodeType_in;
		endNodeType.at(diodeIndex) = endNodeType_in;
		circuitConnectedToStartNode.at(diodeIndex) = circuitConnectedToStartNode_in;
		circuitConnectedToEndNode.at(diodeIndex) = circuitConnectedToEndNode_in;
	}

	circuit_diode_node_t getStartNodeType(int diodeIndex)
	{
		assert(startNodeType.at(diodeIndex) != Node_Null);
		return startNodeType.at(diodeIndex);
	}
	circuit_diode_node_t getEndNodeType(int diodeIndex)
	{
		assert(endNodeType.at(diodeIndex) != Node_Null);
		return endNodeType.at(diodeIndex);
	}
	boost::shared_ptr<CircuitData> getCircuitConnectedToStartNode(int diodeIndex)
	{
		assert(circuitConnectedToStartNode.at(diodeIndex) != NULL);
		return circuitConnectedToStartNode.at(diodeIndex);
	}
	boost::shared_ptr<CircuitData> getCircuitConnectedToEndNode(int diodeIndex)
	{
		assert(circuitConnectedToEndNode.at(diodeIndex) != NULL);
		return circuitConnectedToEndNode.at(diodeIndex);
	}
private:
	std::vector<circuit_diode_node_t> startNodeType;
	std::vector<circuit_diode_node_t> endNodeType;
	std::vector<boost::shared_ptr<CircuitData>> circuitConnectedToStartNode;
	std::vector<boost::shared_ptr<CircuitData>> circuitConnectedToEndNode;
};

// The netlistBoundaryCondition is really a manager class for a collection
// of subcircuits, divided by the diodes/valves in the input data.
// It determines which subcircuits should be used on any given time-step,
// and presents to the boundaryConditionManager an interface which is consistent
// with the other boundary condition types (such as the RCR), so boundaryConditionManager
// needn't concern itself with the fine details.
class NetlistBoundaryCondition : public abstractBoundaryCondition
{
public:
	NetlistBoundaryCondition(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		m_IndexOfThisNetlistLPN = numberOfInitialisedNetlistLPNs;
		numberOfInitialisedNetlistLPNs++;
		initialiseModel();
	}

 	void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);

	void updateLPN();

	~NetlistBoundaryCondition()
	{
		numberOfInitialisedNetlistLPNs--;
	}

private:
	static int numberOfInitialisedNetlistLPNs;
	int m_IndexOfThisNetlistLPN;

	// This vector contains all possible legal subcircuits. The exact ones
	// in use on a given time-step will be determined by diode (valve) state checks..
	// .. a closed valve effectively splits the circuit into two subcircuits
	std::vector<boost::shared_ptr<NetlistSubcircuit>> m_SubcircuitsAsDelineatedByDiodes;
	int m_NumberOfAtomicSubcircuits; // These are what you get with all valves closed: no subcircuit divides an atomic subcircuit.
	int m_NumberOfSubcircuits;
	int m_numberOfDiodes;
	std::map<int,boost:shared_ptr<CircuitComponent>> m_diodeIndexingMap;

	boost::shared_ptr<AtomicSubcircuitConnectionManager> m_atomicSubcircuitConnectionManager;

	void initialiseModel();
	void createCircuitDescription();
	void identifySubciruits();
	void passPressuresAndFlowsToAllSubcircuits();
	void selectActiveSubcircuits();
	void createInitialCircuitDescriptionWithoutDiodes();
	void assignComponentsToAtomicSubcircuits();
	void createSubcircuitDescriptions();
	void buildDiodeMetadata();

	CircuitData m_CircuitDescription;
	CircuitData m_CircuitDescriptionWithoutDiodes;
	std::vector<boost::shared_ptr<CircuitData>> m_CircuitDataForSubcircuits;
	std::vector<boost::shared_ptr<CircuitData>> m_CircuitDataForAtomicSubcircuits;

	std::vector<int> m_AtomicSubcircuitsComponentsBelongsTo; // This is indexed by component, as they appear in m_CircuitDescriptionWithoutDiodes

	std::vector<double> m_PressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	std::vector<double> m_HistoryPressuresInLPN;                // As m_PressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	std::vector<double> m_FlowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	std::vector<double> m_HistoryFlowsInLPN;					  // As m_FlowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	ActiveSubcircuits m_ActiveSubcircuits;
	
	// double P_a;

};

#endif
