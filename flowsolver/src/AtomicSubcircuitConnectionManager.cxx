#include "AtomicSubcircuitConnectionManager.hxx"
#include "datatypesInCpp.hxx"
#include <stdexcept>

circuit_diode_node_t AtomicSubcircuitConnectionManager::getStartNodeType(int diodeIndex)
{
	assert(m_startNodeType.at(diodeIndex) != Node_Null);
	return m_startNodeType.at(diodeIndex);
}
circuit_diode_node_t AtomicSubcircuitConnectionManager::getEndNodeType(int diodeIndex)
{
	assert(m_endNodeType.at(diodeIndex) != Node_Null);
	return m_endNodeType.at(diodeIndex);
}
boost::shared_ptr<CircuitData> AtomicSubcircuitConnectionManager::getCircuitConnectedToStartNode(int diodeIndex)
{
	assert(m_circuitConnectedToStartNode.at(diodeIndex) != NULL);
	return m_circuitConnectedToStartNode.at(diodeIndex);
}
boost::shared_ptr<CircuitData> AtomicSubcircuitConnectionManager::getCircuitConnectedToEndNode(int diodeIndex)
{
	assert(m_circuitConnectedToEndNode.at(diodeIndex) != NULL);
	return m_circuitConnectedToEndNode.at(diodeIndex);
}
int AtomicSubcircuitConnectionManager::getNumberOfDiodes()
{
	return m_numberOfDiodes;
}
bool AtomicSubcircuitConnectionManager::diodeIsOpen(int diodeIdx)
{
	return m_diodesOpen.at(diodeIdx);
}
void AtomicSubcircuitConnectionManager::setDiodeOpen(int diodeIdx, bool isOpen)
{
	m_diodesOpen.at(diodeIdx) = isOpen;
}

void AtomicSubcircuitConnectionManager::setDiodeInfo(int diodeIndex, circuit_diode_node_t startNodeType_in, circuit_diode_node_t endNodeType_in, boost::shared_ptr<CircuitData> circuitConnectedToStartNode_in, boost::shared_ptr<CircuitData> circuitConnectedToEndNode_in)
{
	assert(diodeIndex < m_numberOfDiodes);
	assert(startNodeType_in == Node_Null);
	m_startNodeType.at(diodeIndex) = startNodeType_in;
	m_endNodeType.at(diodeIndex) = endNodeType_in;
	m_circuitConnectedToStartNode.at(diodeIndex) = circuitConnectedToStartNode_in;
	m_circuitConnectedToEndNode.at(diodeIndex) = circuitConnectedToEndNode_in;
}

void AtomicSubcircuitConnectionManager::discoverAtomicSubcircuitsJoinedByEachDiode(const std::vector<boost::shared_ptr<CircuitData>>& circuitDataForAtomicSubcircuits)
{
	for (int diodeIdx = 0; diodeIdx < m_numberOfDiodes; diodeIdx++)
    {
        boost::shared_ptr<CircuitComponent> currentDiode = m_diodeIndexingMap.at(diodeIdx);

        // Prepare a place to temporarily package the two adjoining circuits for this diode:
        std::pair<boost::shared_ptr<CircuitData>, boost::shared_ptr<CircuitData>> adjoiningAtomicSubcircuitsPair;

        // Loop the atomic subcircuits, looking for the subcircuit that the diode connects to:
        for (auto atomicSubcircuit = circuitDataForAtomicSubcircuits.begin(); atomicSubcircuit != circuitDataForAtomicSubcircuits.end(); atomicSubcircuit++)
        {
            if ((*atomicSubcircuit)->mapOfPressureNodes.count(currentDiode->startNode->indexInInputData) == 1)
            {
                adjoiningAtomicSubcircuitsPair.first = *atomicSubcircuit;
            }

            if ((*atomicSubcircuit)->mapOfPressureNodes.count(currentDiode->endNode->indexInInputData) == 1)
            {
                adjoiningAtomicSubcircuitsPair.second = *atomicSubcircuit;
            }
        }

        // Note that in the case where a diode is "hanging" (i.e. the diode has no further components on one side of it)
        // the relevant entry of adjoiningAtomicSubcircuitsPair will be NULL. Careful with this! The AtomicSubcircuitConnectionManager
        // is designed to help avoid making mistakes here.
        circuit_diode_node_t startNodeSubcircuitType;
        if (adjoiningAtomicSubcircuitsPair.first == NULL)
        {
            startNodeSubcircuitType = Node_IsMonopolar;
        }
        else
        {
            startNodeSubcircuitType = Node_ConnectsCircuit;
        }

        circuit_diode_node_t endNodeSubcircuitType;
        if (adjoiningAtomicSubcircuitsPair.second == NULL)
        {
            endNodeSubcircuitType = Node_IsMonopolar;
        }
        else
        {
            endNodeSubcircuitType = Node_ConnectsCircuit;
        }

        setDiodeInfo(diodeIdx,startNodeSubcircuitType,endNodeSubcircuitType,adjoiningAtomicSubcircuitsPair.first,adjoiningAtomicSubcircuitsPair.second);
    }
}

bool AtomicSubcircuitConnectionManager::diodeStartNodeConnectsCircuit(int diodeIdx)
{
    if (m_startNodeType.at(diodeIdx) == Node_ConnectsCircuit)   
    {
        assert(m_circuitConnectedToStartNode.at(diodeIdx) != NULL);
        return true;
    }
    else if(m_startNodeType.at(diodeIdx) == Node_IsMonopolar)
    {
        return false;
    }
    else
    {
        throw std::logic_error("EE: Unknown diode node type.");
    }
}

bool AtomicSubcircuitConnectionManager::diodeEndNodeConnectsCircuit(int diodeIdx)
{
    if (m_endNodeType.at(diodeIdx) == Node_ConnectsCircuit)   
    {
        assert(m_circuitConnectedToEndNode.at(diodeIdx) != NULL);
        return true;
    }
    else if(m_endNodeType.at(diodeIdx) == Node_IsMonopolar)
    {
        return false;
    }
    else
    {
        throw std::logic_error("EE: Unknown diode node type.");
    }
}