#ifndef ATOMICSUBCIRCUITCONNECTIONMANAGER_HXX_
#define ATOMICSUBCIRCUITCONNECTIONMANAGER_HXX_

#include <boost/shared_ptr.hpp>
#include <vector>
#include <map>
#include "CircuitData.hxx"

class AtomicSubcircuitConnectionManager
{
public:

	std::map<int,boost::shared_ptr<CircuitComponent>> m_diodeIndexingMap;

	AtomicSubcircuitConnectionManager(const CircuitData& circuitDescription, const std::vector<boost::shared_ptr<CircuitData>>& circuitDataForAtomicSubcircuits)
	{

		// Count the diodes, get a map to them:
	    m_numberOfDiodes = 0;
	    for (auto component = circuitDescription.components.begin(); component != circuitDescription.components.end(); component++)
	    {
	        if ((*component)->type == Component_Diode)
	        {
	            m_diodeIndexingMap.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> (m_numberOfDiodes, *component));
	            m_numberOfDiodes++;
	        }
	    }

	    // Get a map from diode index to the two atomic subcircuits that it joins, and store that info
	    discoverAtomicSubcircuitsJoinedByEachDiode(circuitDataForAtomicSubcircuits);

		m_startNodeType.insert(m_startNodeType.begin(),m_numberOfDiodes,Node_Null);
		m_endNodeType.insert(m_endNodeType.begin(),m_numberOfDiodes,Node_Null);

		boost::shared_ptr<CircuitData> nullPtr; // NULL doesn't work with boost::shared_ptr, so we make our own NULL.
		m_circuitConnectedToStartNode.insert(m_circuitConnectedToStartNode.begin(),m_numberOfDiodes,nullPtr);
		m_circuitConnectedToEndNode.insert(m_circuitConnectedToEndNode.begin(),m_numberOfDiodes,nullPtr);
		m_diodesOpen.insert(m_diodesOpen.begin(),m_numberOfDiodes,false);
	}

	circuit_diode_node_t getStartNodeType(int diodeIndex);
	circuit_diode_node_t getEndNodeType(int diodeIndex);
	boost::shared_ptr<CircuitData> getCircuitConnectedToStartNode(int diodeIndex);
	boost::shared_ptr<CircuitData> getCircuitConnectedToEndNode(int diodeIndex);
	int getNumberOfDiodes();
	bool diodeIsOpen(int diodeIdx);
	void setDiodeOpen(int diodeIdx, bool isOpen);
	bool diodeStartNodeConnectsCircuit(int diodeIdx);
	bool diodeEndNodeConnectsCircuit(int diodeIdx);
private:
	int m_numberOfDiodes;
	std::vector<circuit_diode_node_t> m_startNodeType;
	std::vector<circuit_diode_node_t> m_endNodeType;
	std::vector<boost::shared_ptr<CircuitData>> m_circuitConnectedToStartNode;
	std::vector<boost::shared_ptr<CircuitData>> m_circuitConnectedToEndNode;
	std::vector<bool> m_diodesOpen;

	void setDiodeInfo(int diodeIndex, circuit_diode_node_t startNodeType_in, circuit_diode_node_t endNodeType_in, boost::shared_ptr<CircuitData> circuitConnectedToStartNode_in, boost::shared_ptr<CircuitData> circuitConnectedToEndNode_in);
	void discoverAtomicSubcircuitsJoinedByEachDiode(const std::vector<boost::shared_ptr<CircuitData>>& circuitDataForAtomicSubcircuits);
};

#endif