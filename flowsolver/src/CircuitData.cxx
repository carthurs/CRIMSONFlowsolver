#include "CircuitData.hxx"

void CircuitData::rebuildCircuitMetadata()
{
	// This subroutine builds everything in the CircuitData class that isn't the CircuitComponents themselves.
	// It parses the CircuitComponents to do this, so it assumes the CircuitComponents are already fully initialised.

	numberOfPrescribedFlows = 0;

	mapOfPressureNodes.clear();
	mapOfComponents.clear();
	mapOfPrescribedPressureNodes.clear();
	for (auto component = components.begin(); component!=components.end(); component++)
	{
		assert((*component)->startNode->prescribedPressureType!=Pressure_Null);
		if((*component)->startNode->prescribedPressureType != Pressure_NotPrescribed)
		{
			mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->indexInInputData, (*component)->startNode));
		}
		mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->indexInInputData, (*component)->startNode));

		assert((*component)->endNode->prescribedPressureType!=Pressure_Null);
		if((*component)->endNode->prescribedPressureType != Pressure_NotPrescribed)
		{
			mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->indexInInputData, (*component)->endNode));
		}
		mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->indexInInputData, (*component)->endNode));

		mapOfComponents.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> ((*component)->indexInInputData,(*component)));

		assert((*component)->prescribedFlowType != Flow_Null);
		if((*component)->prescribedFlowType != Flow_NotPrescribed)
		{
			numberOfPrescribedFlows++;
			mapOfPrescribedFlowComponents.insert(std::pair<int, boost::shared_ptr<CircuitComponent>> ((*component)->indexInInputData, *component));
		}
	}

	numberOfPressureNodes = mapOfPressureNodes.size();
	numberOfPrescribedPressures = mapOfPrescribedPressureNodes.size();
	numberOfComponents = mapOfComponents.size();

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