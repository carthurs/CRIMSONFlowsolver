#ifndef CIRCUITDATA_HXX_
#define CIRCUITDATA_HXX_

#include <map>
#include <vector>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include "datatypesInCpp.hxx"

class CircuitPressureNode
{
public:
	double valueOfPrescribedPressure;
	double pressure;
	double historyPressure;
	bool hasHistoryPressure;
	circuit_nodal_pressure_prescription_t prescribedPressureType;
	int indexInInputData;
	int indexLocalToSubcircuit;
	bool m_connectsTo3DDomain;
	CircuitPressureNode()
	{
		prescribedPressureType = Pressure_Null;
		hasHistoryPressure = false;
	    m_connectsTo3DDomain = false;
	}
};


class CircuitComponent
{
public:
	circuit_component_t type;
	boost::shared_ptr<CircuitPressureNode> startNode;
	boost::shared_ptr<CircuitPressureNode> endNode;
	double parameterValue; // resistance or compliance or inductance etc.
	int indexInInputData;
	int indexLocalToSubcircuit;
	circuit_component_flow_prescription_t prescribedFlowType;
	double valueOfPrescribedFlow;
	double flow;
	double historyFlow;
	bool hasHistoryFlow;
	CircuitComponent() //:
		// startNode(new CircuitPressureNode),
		// endNode(new CircuitPressureNode)
	{
		type = Component_Null;
		prescribedFlowType = Flow_Null;
		hasHistoryFlow = false;
	}
};

class CircuitData
{
public:
	std::vector<boost::shared_ptr<CircuitComponent>> components;
	int index;
	
	// Begin metadata, updated with rebuildCircuitMetadata.
	int numberOfPrescribedPressures;
	int numberOfPrescribedFlows;
	int numberOfPressureNodes;
	int numberOfComponents;

	// These maps have indexInInputData as the key, and a shared_ptr to the relevant node/circuit as the mapped value.
	// Although this provides useful random access by indexInInputData, it is often useful to just use the std::map iterator to process them all.
	std::map<int,boost::shared_ptr<CircuitPressureNode>> mapOfPressureNodes; // Utility data structure, containing all the pressure nodes of the circuit, exactly once each.
	std::map<int,boost::shared_ptr<CircuitPressureNode>> mapOfPrescribedPressureNodes;
	std::map<int,boost::shared_ptr<CircuitComponent>> mapOfComponents; // Utility data structure containing all the component indices of the circuit, exactly once each.
	std::map<int,boost::shared_ptr<CircuitComponent>> mapOfPrescribedFlowComponents;
	// End of medatata
	
	void rebuildCircuitMetadata();
	bool connectsTo3DDomain() const;
	void generateNodeAndComponentIndicesLocalToSubcircuit();
	void tagNodeAt3DInterface();

	boost::shared_ptr<CircuitPressureNode> ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in);
private:
	int toOneIndexing(const int oneIndexedValue);
	void rebuildCircuitPressureNodeMap();
};

#endif
