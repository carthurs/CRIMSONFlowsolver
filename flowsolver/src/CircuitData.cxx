#include "CircuitData.hxx"
#include <stdexcept>

void CircuitData::rebuildCircuitMetadata()
{
	// This subroutine builds everything in the CircuitData class that isn't the CircuitComponents themselves.
	// It parses the CircuitComponents to do this, so it assumes the CircuitComponents are already fully initialised.

	numberOfPrescribedFlows = 0;

	// mapOfPressureNodes.clear();
	mapOfComponents.clear();
	mapOfPrescribedPressureNodes.clear();
	for (auto component = components.begin(); component!=components.end(); component++)
	{
		assert((*component)->startNode->prescribedPressureType!=Pressure_Null);
		if((*component)->startNode->prescribedPressureType != Pressure_NotPrescribed)
		{
			mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->indexInInputData, (*component)->startNode));
		}
		// mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->indexInInputData, (*component)->startNode));

		assert((*component)->endNode->prescribedPressureType!=Pressure_Null);
		if((*component)->endNode->prescribedPressureType != Pressure_NotPrescribed)
		{
			mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->indexInInputData, (*component)->endNode));
		}
		// mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->indexInInputData, (*component)->endNode));

		mapOfComponents.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> ((*component)->indexInInputData,(*component)));

		assert((*component)->prescribedFlowType != Flow_Null);
		if((*component)->prescribedFlowType != Flow_NotPrescribed)
		{
			numberOfPrescribedFlows++;
			mapOfPrescribedFlowComponents.insert(std::pair<int, boost::shared_ptr<CircuitComponent>> ((*component)->indexInInputData, *component));
		}
	}

	rebuildCircuitPressureNodeMap();

	// numberOfPressureNodes = mapOfPressureNodes.size();
	numberOfPrescribedPressures = mapOfPrescribedPressureNodes.size();
	numberOfComponents = mapOfComponents.size();

}

void CircuitData::tagNodeAt3DInterface()
{
	int numberOfComponentsTaggedFor3DFlow = 0; //a counter to verify there exists a unique 3D flow-tagged component
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		// Find the component tagged to recieve flow from the 3D interface:
		if ((*component)->prescribedFlowType == Flow_3DInterface)
		{
			numberOfComponentsTaggedFor3DFlow++;
			// Find which end of the component doesn't link to any other components (and so therefore is at the 3D interface)
			int baseComponentStartNodeIdx=(*component)->startNode->indexInInputData;
			int baseComponentEndNodeIdx=(*component)->endNode->indexInInputData;
			int numberOfTimesComponentStartNodeAppearsInCircuit = 0;
			int numberOfTimesComponentEndNodeAppearsInCircuit = 0;
			for (auto otherComponent=components.begin(); otherComponent!=components.end(); otherComponent++)
			{
				if ((*otherComponent)->startNode->indexInInputData == baseComponentStartNodeIdx)
				{
					numberOfTimesComponentStartNodeAppearsInCircuit++;
				}
				if ((*otherComponent)->endNode->indexInInputData == baseComponentStartNodeIdx)
				{
					numberOfTimesComponentStartNodeAppearsInCircuit++;
				}
				if ((*otherComponent)->startNode->indexInInputData == baseComponentEndNodeIdx)
				{
					numberOfTimesComponentEndNodeAppearsInCircuit++;
				}
				if ((*otherComponent)->endNode->indexInInputData == baseComponentEndNodeIdx)
				{
					numberOfTimesComponentEndNodeAppearsInCircuit++;
				}
			}

			// Make sure both the start and end nodes were actually found in the circuit (or else something has gone badly wrong!)
			assert(numberOfTimesComponentStartNodeAppearsInCircuit>=1);
			assert(numberOfTimesComponentEndNodeAppearsInCircuit>=1);

			if (numberOfTimesComponentStartNodeAppearsInCircuit>1 && numberOfTimesComponentEndNodeAppearsInCircuit>1)
			{
				throw std::runtime_error("EE: Only one component may be directly connected to the 3D interface in the Netlist.");
			}
			else if (numberOfTimesComponentStartNodeAppearsInCircuit>1)
			{
				throw std::runtime_error("EE: The netlist component at the 3D interface must be connected to it via its start node, not its end node.");
			}
			else // actually tag the 3D interface node, now we've checked for errors...
			{
				assert(numberOfTimesComponentStartNodeAppearsInCircuit==1);
				(*component)->startNode->m_connectsTo3DDomain = true; // Setting this bool tag here is the point of the whole subroutine!
			}
		}
	}
	assert(numberOfComponentsTaggedFor3DFlow==1); // May fail if this subroutine gets called on a subcircuit, instead of the whole input data circuit, if the subcircuit doesn't have the 3D flow interface component. Don't call it in this case!
}

void CircuitData::rebuildCircuitPressureNodeMap()
{
	mapOfPressureNodes.clear();
	for (auto component = components.begin(); component!=components.end(); component++)
	{
		if ((*component)->startNode) // ensure the pointer startNode is not null 
		{
			mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->indexInInputData, (*component)->startNode));
		}
		if ((*component)->endNode) // ensure the pointer endNode is not null
		{
			mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->indexInInputData, (*component)->endNode));
		}
	}
	numberOfPressureNodes = mapOfPressureNodes.size();
}

boost::shared_ptr<CircuitPressureNode> CircuitData::ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in)
{
	rebuildCircuitPressureNodeMap();

	bool nodeAlreadyConstructed = (mapOfPressureNodes.count(indexInInputData_in) == 1);
	if (nodeAlreadyConstructed)
	{
		return mapOfPressureNodes.at(indexInInputData_in);
	}
	else // node not already constructed
	{
		CircuitPressureNode* ptrToMakeShared = new CircuitPressureNode;
		return boost::shared_ptr<CircuitPressureNode> (ptrToMakeShared);
	}

}

void CircuitData::generateNodeAndComponentIndicesLocalToSubcircuit()
{
	// Ensure the metadata is up-to-date before we start
	rebuildCircuitMetadata();

	// we take advantage of the fact that std::map<int,*> keeps elements sorted 
	// in increasing value order of int keys
	{
		// Get the circuit-local pressure node indexing
		int indexLocalToSubcircuit_loopCounter = 1;
		for (auto pressureNode = mapOfPressureNodes.begin(); pressureNode != mapOfPressureNodes.end(); pressureNode++)
		{
			// (*pressureNode)->indexLocalToSubcircuit = indexLocalToSubcircuit_loopCounter;
			pressureNode->second->indexLocalToSubcircuit = indexLocalToSubcircuit_loopCounter;
			indexLocalToSubcircuit_loopCounter++;
		}
	}
	// for (int ii=0; ii<numberOfPressureNodes; ii++)
	// {
	// 	for (auto component=components.begin(); component!=components.end(); component++)
	// 	{
	// 		// If the component's index from the input data matches with the
	// 		// ii-th entry in the mapOfPressureNodes for this subcircuit,
	// 		// then we know that the component is the ii-th of the circuit,
	// 		// - this is just because the std::set<int> keeps its entries
	// 		// sorted.
	// 		// We also converto from zero- to one-indexing, for consistency
	// 		// with the usual semantics of the input data.
	// 		if (*component->startNode->indexInInputData == mapOfPressureNodes.at(ii))
	// 		{
	// 			*component->startNode->indexLocalToSubcircuit = toOneIndexing(ii);
	// 		}
	// 		if (*component->endNode->indexInInputData == mapOfPressureNodes.at(ii))
	// 		{
	// 			*component->endNode->indexLocalToSubcircuit = toOneIndexing(ii);
	// 		}
	// 	}
	// }

	{
		// Get the circuit-local component indexing
		int indexLocalToSubcircuit_loopCounter = 1;
		for (auto component = mapOfComponents.begin(); component!=mapOfComponents.end(); component++)
		{
			component->second->indexLocalToSubcircuit = indexLocalToSubcircuit_loopCounter;
		}
	}
	// for (int ii=0; ii<numberOfComponents; ii++)
	// {
	// 	for (auto component=components.begin(); component!=components.end(); component++)
	// 	{
	// 		if(*component->indexInInputData == mapOfComponents.at(ii))
	// 		{
	// 			*component->indexLocalToSubcircuit = toOneIndexing(ii);
	// 		}
	// 	}
	// }

}

bool CircuitData::connectsTo3DDomain() const
{
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		if ((*component)->startNode->m_connectsTo3DDomain)
		{
			return true;
		}
		if ((*component)->endNode->m_connectsTo3DDomain)
		{
			return true;
		}
	}
	return false;
}

inline int CircuitData::toOneIndexing(const int zeroIndexedValue)
{
	int oneIndexedValue = zeroIndexedValue + 1;
	return oneIndexedValue;
}