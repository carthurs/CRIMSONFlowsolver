#include "CircuitData.hxx"
#include "indexShifters.hxx"
#include <stdexcept>
#include <float.h>
#include <algorithm>
#include <stack>
#include <sstream>

void VolumeTrackingPressureChamber::passPressureToStartNode()
{
	m_pressure = (m_storedVolume - m_unstressedVolume)/currentParameterValue;
	std::cout << "compliance set to: " << currentParameterValue << std::endl;
	std::cout << "pressure: " << m_pressure << std::endl;
	startNode->setPressure(m_pressure);
}

bool CircuitComponent::hasNonnegativePressureGradientOrForwardFlow() // whether the diode should be open
{
	bool hasNonnegativePressureGradient = (startNode->getPressure() >= endNode->getPressure());
	bool hasForwardFlow = (signForPrescribed3DInterfaceFlow*flow >= 1e-16); // We use 1e-16 because it's essentially zero. Diode closure is enforced by setting diode resistance to DBL_MAX, so there remains a small flow on the order 1e-308 across a closed diode.
	return (hasNonnegativePressureGradient || hasForwardFlow);
}

void CircuitData::rebuildCircuitMetadata()
{
	// This subroutine builds everything in the CircuitData class that isn't the CircuitComponents themselves.
	// It parses the CircuitComponents to do this, so it assumes the CircuitComponents are already fully initialised.

	numberOfPrescribedFlows = 0;

	// mapOfPressureNodes.clear();
	mapOfComponents.clear();
	mapOfPrescribedPressureNodes.clear();
	mapOfPrescribedFlowComponents.clear();
	mapOfVolumeTrackingComponents.clear();
	mapOfPrescribedVolumeTrackingComponents.clear();
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
		// bool isDiode = ((*component)->type == Component_Diode);
		// bool diodeIsClosedSoPrescribeZeroFlow = !((*component)->hasNonnegativePressureGradientAndNoBackflow());
		if((*component)->prescribedFlowType != Flow_NotPrescribed) // || (diodeIsClosedSoPrescribeZeroFlow && isDiode))
		{
			numberOfPrescribedFlows++;
			mapOfPrescribedFlowComponents.insert(std::pair<int, boost::shared_ptr<CircuitComponent>> ((*component)->indexInInputData, *component));
		}
	}

	rebuildCircuitPressureNodeMap();

	setupComponentNeighbourPointers();

	numberOfPrescribedPressures = mapOfPrescribedPressureNodes.size();
	numberOfComponents = mapOfComponents.size();

	// Get the number of VolumeTrackingPressureChambers,
	// and generate mapOfVolumeTrackingComponents:
	m_numberOfVolumeTrackingPressureChambers = 0;
	for (auto component=mapOfComponents.begin(); component!=mapOfComponents.end(); component++)
	{
		if (component->second->type == Component_VolumeTrackingPressureChamber)
		{
			mapOfVolumeTrackingComponents.insert(std::make_pair(component->first,component->second));
			m_numberOfVolumeTrackingPressureChambers++;
		}
	}

	// mapOfPrescribedVolumeTrackingComponents
	// for (auto component=components.begin(); component!=components.end(); component++)
	// {
	// 	if ((*component)->type == Component_VolumeTrackingPressureChamber)
	// 	{
	// 		VolumeTrackingPressureChamber* volumeTrackingPressureChamber = dynamic_cast<VolumeTrackingPressureChamber> ((*component).get());
	// 		if (volumeTrackingPressureChamber->
	// 	}
	// }

}

void CircuitData::setupComponentNeighbourPointers()
{
	for(auto component = components.begin(); component!=components.end(); component++)
	{
		(*component)->neighbouringComponentsAtStartNode.clear();
		(*component)->neighbouringComponentsAtEndNode.clear();

		int startNodeIdx = (*component)->startNode->indexInInputData;
		int endNodeIdx = (*component)->endNode->indexInInputData;
		for (auto possibleNeighbouringComponent=components.begin(); possibleNeighbouringComponent!=components.end(); possibleNeighbouringComponent++)
		{
			// Avoid giving a component itself as a neighbour:
			if (component!=possibleNeighbouringComponent)
			{
				int neighbourStartNodeIdx = (*possibleNeighbouringComponent)->startNode->indexInInputData;
				int neighbourEndNodeIdx = (*possibleNeighbouringComponent)->endNode->indexInInputData;
				// check whether these nodes match, if so, give the component a pointer to its newly-discovered neighbour
				
				bool isANeighbourViaComponentStartNode = (neighbourStartNodeIdx == startNodeIdx ||
														    neighbourEndNodeIdx == startNodeIdx );

				bool isANeighbourViaComponentEndNode = (neighbourStartNodeIdx == endNodeIdx ||
														  neighbourEndNodeIdx == endNodeIdx );

				if(isANeighbourViaComponentStartNode)
				{
					boost::weak_ptr<CircuitComponent> toPushBack(*possibleNeighbouringComponent);
					(*component)->neighbouringComponentsAtStartNode.push_back(toPushBack);
				}

				if(isANeighbourViaComponentEndNode)
				{
					boost::weak_ptr<CircuitComponent> toPushBack(*possibleNeighbouringComponent);
					(*component)->neighbouringComponentsAtEndNode.push_back(toPushBack);
				}
			}

		}

		// Reverse the neighbour vectors so that components are listed in the same order as they appear in the input data netlist_surfaces.dat
		// This keeps things consistent.
		std::reverse((*component)->neighbouringComponentsAtStartNode.begin(),(*component)->neighbouringComponentsAtStartNode.end());
		std::reverse((*component)->neighbouringComponentsAtEndNode.begin(),(*component)->neighbouringComponentsAtEndNode.end());
	}
}

void CircuitData::initialiseNodeAndComponentAt3DInterface(int threeDInterfaceNodeIndex)
{
	setIndexOfNodeAt3DInterface(threeDInterfaceNodeIndex);
	// tag the node at the 3D interface:
    mapOfPressureNodes.at(threeDInterfaceNodeIndex)->m_connectsTo3DDomain = true;

	int numberOfComponentsTaggedFor3DFlow = 0; //a counter to verify there exists a unique 3D flow-tagged component
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		// Find the component tagged to recieve flow from the 3D interface:
		if ((*component)->prescribedFlowType == Flow_3DInterface)
		{
			numberOfComponentsTaggedFor3DFlow++;
			// Ensure one of this component's nodes is the one at the 3D interface
			bool startNodeIsAt3Dinterface = ((*component)->startNode->indexInInputData  ==  m_indexOfNodeAt3DInterface);
			bool endNodeIsAt3Dinterface = ((*component)->endNode->indexInInputData  ==  m_indexOfNodeAt3DInterface);
			if (startNodeIsAt3Dinterface || endNodeIsAt3Dinterface)
			{
				// tag the component as being at the 3D interface
				(*component)->setConnectsToNodeAt3DInterface();
			}
			else
			{
				throw std::runtime_error("EE: Component with prescribed flow from the 3D domain does not have a node tagged as being at the 3D interface.");
			}

			// Set up the sign for the flow prescription, depending on whether it's the start or end node of the
			// component which is at the 3D interface (if the user has listed the end node before the start node for
			// the component at the 3D interface in the netlist_surfaces.dat, we need to flip the sign of the flow
			// coming from the 3D interface. This affects all component types, but is needed for diodes in particular,
			// where the user has a legitimate reason to want to reverse those nodes, depending on the orientation they
			// want to give to the diode.)
			if (startNodeIsAt3Dinterface)
			{
				(*component)->signForPrescribed3DInterfaceFlow = 1.0;
			}
			else if (endNodeIsAt3Dinterface)
			{
				(*component)->signForPrescribed3DInterfaceFlow = -1.0;
			}
			else
			{
				throw std::runtime_error("EE: Component with prescribed flow from the 3D domain does not have a node tagged as being at the 3D interface.");
			}

			// // Find which end of the component doesn't link to any other components (and so therefore is at the 3D interface)
			// int baseComponentStartNodeIdx=(*component)->startNode->indexInInputData;
			// int baseComponentEndNodeIdx=(*component)->endNode->indexInInputData;
			// int numberOfTimesComponentStartNodeAppearsInCircuit = 0;
			// int numberOfTimesComponentEndNodeAppearsInCircuit = 0;
			// for (auto otherComponent=components.begin(); otherComponent!=components.end(); otherComponent++)
			// {
			// 	if ((*otherComponent)->startNode->indexInInputData == baseComponentStartNodeIdx)
			// 	{
			// 		numberOfTimesComponentStartNodeAppearsInCircuit++;
			// 	}
			// 	if ((*otherComponent)->endNode->indexInInputData == baseComponentStartNodeIdx)
			// 	{
			// 		numberOfTimesComponentStartNodeAppearsInCircuit++;
			// 	}
			// 	if ((*otherComponent)->startNode->indexInInputData == baseComponentEndNodeIdx)
			// 	{
			// 		numberOfTimesComponentEndNodeAppearsInCircuit++;
			// 	}
			// 	if ((*otherComponent)->endNode->indexInInputData == baseComponentEndNodeIdx)
			// 	{
			// 		numberOfTimesComponentEndNodeAppearsInCircuit++;
			// 	}
			// }

			// // Make sure both the start and end nodes were actually found in the circuit (or else something has gone badly wrong!)
			// assert(numberOfTimesComponentStartNodeAppearsInCircuit>=1);
			// assert(numberOfTimesComponentEndNodeAppearsInCircuit>=1);

			// if (numberOfTimesComponentStartNodeAppearsInCircuit>1 && numberOfTimesComponentEndNodeAppearsInCircuit>1)
			// {
			// 	throw std::runtime_error("EE: Only one component may be directly connected to the 3D interface in the Netlist.");
			// }
			// else if (numberOfTimesComponentStartNodeAppearsInCircuit>1)
			// {
			// 	throw std::runtime_error("EE: The netlist component at the 3D interface must be connected to it via its start node, not its end node.");
			// }
			// else // actually tag the 3D interface node, now we've checked for errors...
			// {
			// 	assert(numberOfTimesComponentStartNodeAppearsInCircuit==1);
			// 	// Setting these two bool tags here is the purpose of the whole subroutine!
			// 	// (*component)->startNode->m_connectsTo3DDomain = true;
			// 	(*component)->setConnectsToNodeAt3DInterface();
			// }
		}
	}
	if (numberOfComponentsTaggedFor3DFlow!=1)
	{
		// May fail if this subroutine gets called on a subcircuit, instead of the whole input data circuit, if the subcircuit doesn't have the 3D flow interface component. Don't call it in this case!
		std::stringstream errorMessage;
		errorMessage << "EE: Expected exactly one component at the 3D interface, but found " << numberOfComponentsTaggedFor3DFlow << "." << std::endl;
		throw std::runtime_error(errorMessage.str());
	}

	// Ensure that the 3D interface node belongs to a unique component:
	int numberOfComponentsConnectingTo3DInterfaceNode = 0; //a counter to verify there exists a unique component connecting to the 3D interface node
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		bool startNodeIsAt3Dinterface = ((*component)->startNode->indexInInputData  ==  m_indexOfNodeAt3DInterface);
		bool endNodeIsAt3Dinterface = ((*component)->endNode->indexInInputData  ==  m_indexOfNodeAt3DInterface);
		if (startNodeIsAt3Dinterface)
		{
			numberOfComponentsConnectingTo3DInterfaceNode++;
		}
		if (endNodeIsAt3Dinterface)
		{
			numberOfComponentsConnectingTo3DInterfaceNode++;
		}
	}
	if (numberOfComponentsConnectingTo3DInterfaceNode!=1)
	{
		// May fail if this subroutine gets called on a subcircuit, instead of the whole input data circuit, if the subcircuit doesn't have the 3D flow interface component. Don't call it in this case!
		std::stringstream errorMessage;
		errorMessage << "EE: Expected exactly one appearance of the 3D interface node in the circuit data, but found " << numberOfComponentsConnectingTo3DInterfaceNode << "." << std::endl;
		throw std::runtime_error(errorMessage.str());
	}

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

boost::shared_ptr<CircuitPressureNode> CircuitData::ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode)
{
	rebuildCircuitPressureNodeMap();

	bool nodeAlreadyConstructed = (mapOfPressureNodes.count(indexInInputData_in) == 1);
	if (nodeAlreadyConstructed)
	{
		// add the pressure node to the list of neighbours, so the node knows which components are attached to it:
		boost::weak_ptr<CircuitComponent> componentToPushBack(componentNeighbouringThisNode);
		mapOfPressureNodes.at(indexInInputData_in)->listOfComponentstAttachedToThisNode.push_back(componentToPushBack);
		// Return this existing node:
		return mapOfPressureNodes.at(indexInInputData_in);
	}
	else // node not already constructed
	{
		// Check for a non-set pressure type first:
		if (typeOfPrescribedPressure == Pressure_Null)
		{
			std::stringstream message;
			message << "Invalid pressure prescription for netlist circuit node number " <<  indexInInputData_in << std::endl;
			throw std::runtime_error(message.str());
		}

		// Build the correct sort of node:
		CircuitPressureNode* ptrToNewNode;
		ptrToNewNode = new CircuitPressureNode(indexInInputData_in, typeOfPrescribedPressure, m_hstep);

		// add the pressure node to the list of neighbours, so the node knows which components are attached to it:
		boost::weak_ptr<CircuitComponent> componentToPushBack(componentNeighbouringThisNode);
		ptrToNewNode->listOfComponentstAttachedToThisNode.push_back(componentToPushBack);

		// Finally, make a shared pointer to the new node and return it.
		return boost::shared_ptr<CircuitPressureNode> (ptrToNewNode);
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

void CircuitData::switchDiodeStatesIfNecessary()
{
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		if ((*component)->type == Component_Diode)
		{
			// (*component)->prescribedFlowType = Flow_Diode_FixedWhenClosed;
			// (*component)->valueOfPrescribedFlow = 0.0; // For enforcing zero flow when the diode is closed
			bool diodeIsOpen = (*component)->hasNonnegativePressureGradientOrForwardFlow();
			if (diodeIsOpen)
			{
				(*component)->currentParameterValue = (*component)->parameterValueFromInputData; // For enforcing zero resistance when the diode is open
				(*component)->permitsFlow = true;
			}
			else
			{
				(*component)->currentParameterValue = DBL_MAX; // For enforcing "infinite" resistance when the diode is open
				(*component)->permitsFlow = false;
			}
		}
	}
}

void CircuitData::detectWhetherClosedDiodesStopAllFlowAt3DInterface()
{
	bool previousStateOf_m_flowPermittedAcross3DInterface = m_flowPermittedAcross3DInterface;
	m_flowPermittedAcross3DInterface = false;
	// Find the index of the component connected to the 3D domain, so we can use it as a starting point for walking the component tree::
	int indexOfComponentAt3DInterface = -1;
	for (auto component = components.begin(); component!= components.end(); component++)
	{
		if ((*component)->connectsToNodeAt3DInterface())
		{
			indexOfComponentAt3DInterface = (*component)->indexInInputData;
		}
	}

	assert(indexOfComponentAt3DInterface != -1);

	std::stack<boost::weak_ptr<CircuitComponent>> componentsNeedingChecking;
	boost::weak_ptr<CircuitComponent> toPushOntoStack(mapOfComponents.at(indexOfComponentAt3DInterface));
	componentsNeedingChecking.push(toPushOntoStack);
	// To keep track of which components have been already checked:
	std::vector<bool> componentsWhichHaveBeenChecked(components.size(),false);
	while(!componentsNeedingChecking.empty())
	{
		// pop the stack to set the current component:
		boost::weak_ptr<CircuitComponent> currentComponent(componentsNeedingChecking.top());
		componentsNeedingChecking.pop();
		// For the current component:
		// Ensure we've not done this component yet (avoids circular problems)
		if (componentsWhichHaveBeenChecked.at(toZeroIndexing(currentComponent.lock()->indexInInputData)) == false)
		{
			// note that we're checking this component:
			componentsWhichHaveBeenChecked.at(toZeroIndexing(currentComponent.lock()->indexInInputData)) = true;

			// Don't parse the neighbours if flow is banned (i.e. if there's a closed diode)
			if (currentComponent.lock()->permitsFlow)
			{
				// Discover whether currentComponent has no neighbours at one end (i.e. it's a flow sink), so there is
				// somewhere for flow coming in at the 3D domain to go (i.e. it's OK for the surface to have a Dirichlet boundary condition).
				// We also ensure that we haven't accidentally detected the 3D interface node itself, using the bools.
				int numberOfStartNodeNeighbours = currentComponent.lock()->neighbouringComponentsAtStartNode.size();
				bool startNodeNotAt3DInterface = !(currentComponent.lock()->startNode->m_connectsTo3DDomain);
				int numberOfEndNodeNeighbours = currentComponent.lock()->neighbouringComponentsAtEndNode.size();
				bool endNodeNotAt3DInterface = !(currentComponent.lock()->endNode->m_connectsTo3DDomain);
				if ((numberOfStartNodeNeighbours == 0 && startNodeNotAt3DInterface) || (numberOfEndNodeNeighbours == 0 && endNodeNotAt3DInterface))
				{
					m_flowPermittedAcross3DInterface = true;
					// we've found  what we were looking for, so break out
					break;
				}

				// Put all the neighbours on a stack
				for (auto neighbouringComponent=currentComponent.lock()->neighbouringComponentsAtEndNode.begin(); neighbouringComponent!=currentComponent.lock()->neighbouringComponentsAtEndNode.end(); neighbouringComponent++)
				{
					componentsNeedingChecking.push(*neighbouringComponent);
				}
				for (auto neighbouringComponent=currentComponent.lock()->neighbouringComponentsAtStartNode.begin(); neighbouringComponent!=currentComponent.lock()->neighbouringComponentsAtStartNode.end(); neighbouringComponent++)
				{
					componentsNeedingChecking.push(*neighbouringComponent);
				}
			}
		}
	}
	// Check whether flow has just become possible across this boundary (or vice-versa), when compared to the previous time-step.
	// If this is the case, the boolean will reach the linear system to tell it that it needs rebuilding
	// for the new boundary condition type (Dirichlet/Neumann).
	if (previousStateOf_m_flowPermittedAcross3DInterface != m_flowPermittedAcross3DInterface)
	{
		m_boundaryConditionTypeHasJustChanged = true;
	}

	// If necessary, adjust the circuit data as appropriate for a change in boundary condition type
	// between Neumann and Dirichlet.
	if (m_boundaryConditionTypeHasJustChanged)
	{
		switchBetweenDirichletAndNeumannCircuitDesign();
		m_boundaryConditionTypeHasJustChanged = false;
	}
}

// Adjusts the circuit data as appropriate for a change in boundary condition type
// between Neumann and Dirichlet.
void CircuitData::switchBetweenDirichletAndNeumannCircuitDesign()
{
	if (m_flowPermittedAcross3DInterface) // The condition we have just changed to is Neumann from Dirichlet:
	{
		std::cout << "Switching to Neumann in CircuitData.cxx!" << std::endl;
		// Remove the node at the 3D interface from the list of those with prescribed pressure:
		for (auto node=mapOfPrescribedPressureNodes.begin(); node!=mapOfPrescribedPressureNodes.end(); node++)
		{
			if (node->second->m_connectsTo3DDomain)
			{
				node->second->prescribedPressureType=Pressure_NotPrescribed;
				mapOfPrescribedPressureNodes.erase(node);
				break;
			}
		}

		// Add the component at the 3D interface to the list of those with prescribed flow:
		for (auto component=mapOfComponents.begin(); component!=mapOfComponents.end(); component++)
		{
			if (component->second->connectsToNodeAt3DInterface())
			{
				component->second->prescribedFlowType=Flow_3DInterface;
				mapOfPrescribedFlowComponents.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> (component->second->indexInInputData,component->second));
				break;
			}
		}

		// Update the counts of each type of prescription (in switching to a Neumann condition, we're now
		// giving a prescribed pressure to the 3D domain, which means this boundary condition will be
		// receiving a flow back from the 3D domain to prescribe at the LPN interface, instead of the
		// pressure it was previously receiving under Dirichlet conditions. These counters therefore
		// need updating.)
		numberOfPrescribedPressures--;
		numberOfPrescribedFlows++;
	}
	else // Else we have just changed to a Dirichlet condition
	{
		std::cout << "Switching to Dirichlet in CircuitData.cxx!" << std::endl;
		// Add the node at the 3D interface to the list of those with prescribed pressure:
		for (auto node=mapOfPressureNodes.begin(); node!=mapOfPressureNodes.end(); node++)
		{
			if (node->second->m_connectsTo3DDomain)
			{
				node->second->prescribedPressureType=Pressure_3DInterface;
				mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> (node->second->indexInInputData,node->second));
				break;
			}
		}

		// Remove the component at the 3D interface from the list of those with prescribed flow:
		for (auto component=mapOfPrescribedFlowComponents.begin(); component!=mapOfPrescribedFlowComponents.end(); component++)
		{
			if (component->second->connectsToNodeAt3DInterface())
			{
				component->second->prescribedFlowType=Flow_NotPrescribed;
				mapOfPrescribedFlowComponents.erase(component);
				break;
			}
		}

		// Update the counts of each type of prescription (in switching to a Neumann condition, we're now
		// giving a prescribed flow to the 3D domain, which means this boundary condition will be
		// receiving a pressure back from the 3D domain to prescribe at the LPN interface, instead of the
		// flow it was previously receiving under Neumann conditions. These counters therefore
		// need updating.)
		numberOfPrescribedPressures++;
		numberOfPrescribedFlows--;
	}	
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

bool CircuitData::flowPermittedAcross3DInterface() const
{
	return m_flowPermittedAcross3DInterface;
}

bool CircuitData::boundaryConditionTypeHasJustChanged()
{
	bool returnVal = m_boundaryConditionTypeHasJustChanged;
	// Reset m_boundaryConditionTypeHasJustChanged before returning:
	m_boundaryConditionTypeHasJustChanged = false;
	return m_boundaryConditionTypeHasJustChanged;
}

void CircuitData::setIndexOfNodeAt3DInterface(int indexToSet)
{
	m_indexOfNodeAt3DInterface = indexToSet;
}

int CircuitData::getIndexOfNodeAt3DInterface()
{
	return m_indexOfNodeAt3DInterface;
}

bool CircuitComponent::connectsToNodeAt3DInterface()
{
	return m_connectsToNodeAt3DInterface;
}

void CircuitComponent::setConnectsToNodeAt3DInterface()
{
	m_connectsToNodeAt3DInterface = true;
}

// inline int CircuitData::toOneIndexing(const int zeroIndexedValue)
// {
// 	int oneIndexedValue = zeroIndexedValue + 1;
// 	return oneIndexedValue;
// }

