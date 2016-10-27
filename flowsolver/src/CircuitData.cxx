#include "CircuitData.hxx"
#include "indexShifters.hxx"
#include <stdexcept>
#include <algorithm>
#include <stack>
#include <sstream>
#include <boost/make_shared.hpp>

#include "common_c.h"


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
		assert((*component)->startNode->getPressurePrescriptionType() != Pressure_Null);
		if((*component)->startNode->getPressurePrescriptionType() != Pressure_NotPrescribed)
		{
			mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->getIndex(), (*component)->startNode));
		}
		// mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->indexInInputData, (*component)->startNode));

		assert((*component)->endNode->getPressurePrescriptionType() != Pressure_Null);
		if((*component)->endNode->getPressurePrescriptionType() != Pressure_NotPrescribed)
		{
			mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->getIndex(), (*component)->endNode));
		}
		// mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->endNode->indexInInputData, (*component)->endNode));

		mapOfComponents.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> ((*component)->getIndex(),(*component)));

		assert((*component)->prescribedFlowType != Flow_Null);
		// bool isDiode = ((*component)->type == Component_Diode);
		// bool diodeIsClosedSoPrescribeZeroFlow = !((*component)->hasNonnegativePressureGradientAndNoBackflow());
		if((*component)->prescribedFlowType != Flow_NotPrescribed) // || (diodeIsClosedSoPrescribeZeroFlow && isDiode))
		{
			numberOfPrescribedFlows++;
			mapOfPrescribedFlowComponents.insert(std::pair<int, boost::shared_ptr<CircuitComponent>> ((*component)->getIndex(), *component));
		}
	}

	rebuildCircuitPressureNodeMap();

	// setupComponentNeighbourPointers();

	numberOfPrescribedPressures = mapOfPrescribedPressureNodes.size();
	numberOfComponents = mapOfComponents.size();

	// Get the number of VolumeTrackingComponents,
	// and generate mapOfVolumeTrackingComponents:
	m_numberOfVolumeTrackingComponenets = 0;
	for (auto component=mapOfComponents.begin(); component!=mapOfComponents.end(); component++)
	{
		if (boost::dynamic_pointer_cast<VolumeTrackingComponent> (component->second))
		{
			mapOfVolumeTrackingComponents.insert(std::make_pair(component->first,component->second));
			m_numberOfVolumeTrackingComponenets++;
		}
	}

	// mapOfPrescribedVolumeTrackingComponents
	// for (auto component=components.begin(); component!=components.end(); component++)
	// {
	// 	if ((*component)->type == Component_VolumeTrackingPressureChamber)
	// 	{
	// 		VolumeTrackingComponent* volumeTrackingPressureChamber = dynamic_cast<VolumeTrackingComponent> ((*component).get());
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

		int startNodeIdx = (*component)->startNode->getIndex();
		int endNodeIdx = (*component)->endNode->getIndex();
		for (auto possibleNeighbouringComponent=components.begin(); possibleNeighbouringComponent!=components.end(); possibleNeighbouringComponent++)
		{
			// Avoid giving a component itself as a neighbour:
			if (component!=possibleNeighbouringComponent)
			{
				int neighbourStartNodeIdx = (*possibleNeighbouringComponent)->startNode->getIndex();
				int neighbourEndNodeIdx = (*possibleNeighbouringComponent)->endNode->getIndex();
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

void CircuitData::initialiseNodeAndComponentAtInterface(int threeDInterfaceNodeIndex)
{
	{
		std::vector<int> vectorToSet;
		vectorToSet.push_back(threeDInterfaceNodeIndex);
		setIndicesOfNodesAtInterface(vectorToSet);
	}
	// tag the node at the 3D interface:
	try {
    	mapOfPressureNodes.at(threeDInterfaceNodeIndex)->setIsAtBoundary();
    } catch (const std::exception& e) {
    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
    	throw;
    }

	int numberOfComponentsTaggedFor3DFlow = 0; //a counter to verify there exists a unique 3D flow-tagged component
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		// Find the component tagged to recieve flow from the 3D interface:
		if ((*component)->prescribedFlowType == Flow_3DInterface)
		{
			numberOfComponentsTaggedFor3DFlow++;
			// Ensure one of this component's nodes is the one at the 3D interface
			bool startNodeIsAt3Dinterface = ((*component)->startNode->getIndex()  ==  getIndexOfNodeAtInterface());
			bool endNodeIsAt3Dinterface = ((*component)->endNode->getIndex()  ==  getIndexOfNodeAtInterface());
			if (startNodeIsAt3Dinterface || endNodeIsAt3Dinterface)
			{
				// tag the component as being at the 3D interface
				(*component)->setConnectsToNodeAtInterface();
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
				m_signForPrescribed3DInterfaceFlow = 1.0;
				std::cout << "start node is at 3d interface" << std::endl;
			}
			else if (endNodeIsAt3Dinterface)
			{
				m_signForPrescribed3DInterfaceFlow = -1.0;
				std::cout << "end node is at 3d interface" << std::endl;
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
			// 	// (*component)->startNode->m_isAtBoundary = true;
			// 	(*component)->setConnectsToNodeAtInterface();
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
		bool startNodeIsAt3Dinterface = ((*component)->startNode->getIndex()  ==  getIndexOfNodeAtInterface());
		bool endNodeIsAt3Dinterface = ((*component)->endNode->getIndex()  ==  getIndexOfNodeAtInterface());
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

double CircuitData::getSignForPrescribed3DInterfaceFlow() const
{
	return m_signForPrescribed3DInterfaceFlow;
}

void CircuitData::rebuildCircuitPressureNodeMap()
{
	mapOfPressureNodes.clear();
	for (auto component = components.begin(); component!=components.end(); component++)
	{
		if ((*component)->startNode) // ensure the pointer startNode is not null 
		{
			mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> ((*component)->startNode->getIndex(), (*component)->startNode));
		}
		if ((*component)->endNode) // ensure the pointer endNode is not null
		{
			int index = (*component)->endNode->getIndex();
			boost::shared_ptr<CircuitPressureNode> node = (*component)->endNode;
			mapOfPressureNodes[index] = node;
			// mapOfPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> (index, node));
		}
	}
	numberOfPressureNodes = mapOfPressureNodes.size();
}

int CircuitData::getIndexOfComponentConnectingToNodeAtInterface()
{
	int componentIndexToReturn = 0;
	for (auto component=components.begin(); component != components.end(); component++)
	{
		componentIndexToReturn++;
		if ((*component)->connectsToNodeAtInterface())
		{
			return componentIndexToReturn;
		}
	}
	std::stringstream errorMessage;
	errorMessage << "EE: Internal error: unreachable location in getIndexOfComponentConnectingToNodeAtInterface reached." << std::endl;
	throw std::logic_error(errorMessage.str());
}

boost::shared_ptr<CircuitPressureNode> CircuitData::ifExistsGetNodeOtherwiseConstructNode(const int indexInInputData_in, const circuit_nodal_pressure_prescription_t typeOfPrescribedPressure, const boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode)
{
	rebuildCircuitPressureNodeMap();

	bool nodeAlreadyConstructed = (mapOfPressureNodes.count(indexInInputData_in) == 1);
	if (nodeAlreadyConstructed)
	{
		// add the pressure node to the list of neighbours, so the node knows which components are attached to it:
		boost::weak_ptr<CircuitComponent> componentToPushBack(componentNeighbouringThisNode);
		try {
			mapOfPressureNodes.at(indexInInputData_in)->listOfComponentstAttachedToThisNode.push_back(componentToPushBack);
		} catch (const std::exception& e) {
	    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    	throw;
	    }
		
		// Return this existing node:
		boost::shared_ptr<CircuitPressureNode> returnValue;
		try {
			returnValue = mapOfPressureNodes.at(indexInInputData_in);
		} catch (const std::exception& e) {
	    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    	throw;
	    }

		return returnValue;
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

void CircuitData::switchDiodeStatesIfNecessary()
{
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		if ((*component)->getType() == Component_Diode)
		{
			// (*component)->prescribedFlowType = Flow_Diode_FixedWhenClosed;
			// (*component)->valueOfPrescribedFlow = 0.0; // For enforcing zero flow when the diode is closed
			bool diodeIsOpen = (*component)->hasNonnegativePressureGradientOrForwardFlow();
			if (diodeIsOpen)
			{
				(*component)->enableDiodeFlow();
			}
			else
			{
				(*component)->disableDiodeFlow();
			}
		}
	}
}

// In particular, this is for starting a simulation with all diodes shut, for stability.
void CircuitData::closeAllDiodes()
{
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		if ((*component)->getType() == Component_Diode)
		{
			(*component)->disableDiodeFlow();
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
		if ((*component)->connectsToNodeAtInterface())
		{
			indexOfComponentAt3DInterface = (*component)->getIndex();
		}
	}

	assert(indexOfComponentAt3DInterface != -1);

	std::stack<boost::weak_ptr<CircuitComponent>> componentsNeedingChecking;
	boost::weak_ptr<CircuitComponent> toPushOntoStack;
	try {
		toPushOntoStack = boost::weak_ptr<CircuitComponent>(mapOfComponents.at(indexOfComponentAt3DInterface));
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw;
	}
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
		try {
			if (componentsWhichHaveBeenChecked.at(toZeroIndexing(currentComponent.lock()->getIndex())) == false)
				{
					// note that we're checking this component:
					componentsWhichHaveBeenChecked.at(toZeroIndexing(currentComponent.lock()->getIndex())) = true;
		
					// Don't parse the neighbours if flow is banned (i.e. if there's a closed diode)
					if (currentComponent.lock()->permitsFlow())
					{
						// Discover whether currentComponent has no neighbours at one end (i.e. it's a flow sink), so there is
						// somewhere for flow coming in at the 3D domain to go (i.e. it's OK for the surface to have a Dirichlet boundary condition).
						// We also ensure that we haven't accidentally detected the 3D interface node itself, using the bools.
						int numberOfStartNodeNeighbours = currentComponent.lock()->neighbouringComponentsAtStartNode.size();
						bool startNodeNotAt3DInterface = !(currentComponent.lock()->startNode->isAtBoundary());
						int numberOfEndNodeNeighbours = currentComponent.lock()->neighbouringComponentsAtEndNode.size();
						bool endNodeNotAt3DInterface = !(currentComponent.lock()->endNode->isAtBoundary());
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
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw;
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
			if (node->second->isAtBoundary())
			{
				node->second->setPressurePrescriptionType(Pressure_NotPrescribed);
				mapOfPrescribedPressureNodes.erase(node);
				break;
			}
		}

		// Add the component at the 3D interface to the list of those with prescribed flow:
		for (auto component=mapOfComponents.begin(); component!=mapOfComponents.end(); component++)
		{
			if (component->second->connectsToNodeAtInterface())
			{
				component->second->prescribedFlowType=Flow_3DInterface;
				mapOfPrescribedFlowComponents.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> (component->second->getIndex(),component->second));
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
		assert(numberOfPrescribedPressures>=0);
	}
	else // Else we have just changed to a Dirichlet condition
	{
		std::cout << "Switching to Dirichlet in CircuitData.cxx!" << std::endl;
		// Add the node at the 3D interface to the list of those with prescribed pressure:
		for (auto node=mapOfPressureNodes.begin(); node!=mapOfPressureNodes.end(); node++)
		{
			if (node->second->isAtBoundary())
			{
				node->second->setPressurePrescriptionType(Pressure_3DInterface);
				mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> (node->second->getIndex(),node->second));
				break;
			}
		}

		// Remove the component at the 3D interface from the list of those with prescribed flow:
		for (auto component=mapOfPrescribedFlowComponents.begin(); component!=mapOfPrescribedFlowComponents.end(); component++)
		{
			if (component->second->connectsToNodeAtInterface())
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
		assert(numberOfPrescribedFlows>=0);
	}	
}

bool CircuitData::connectsTo3DDomain() const
{
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		if ((*component)->startNode->isAtBoundary())
		{
			return true;
		}
		if ((*component)->endNode->isAtBoundary())
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

void CircuitData::setIndicesOfNodesAtInterface(std::vector<int> indicesToSet)
{
	m_indexOfNodeAt3DInterface = indicesToSet;
}

int CircuitData::getIndexOfNodeAtInterface()
{
	try {
		return m_indexOfNodeAt3DInterface.at(0); // The basic netlist boundary condition currently only uses the 0th entry of m_indexOfNodeAt3DInterface.
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw;
	}
}

inline int CircuitData::toOneIndexing(const int zeroIndexedValue)
{
	int oneIndexedValue = zeroIndexedValue + 1;
	return oneIndexedValue;
}

bool CircuitData::hasPrescribedFlowAcrossInterface() const
{
	return m_flowPermittedAcross3DInterface;
}

bool CircuitData::hasPrescribedPressureAcrossInterface() const
{
	// Negate and return:
	return !m_flowPermittedAcross3DInterface;
}

boost::shared_ptr<CircuitComponent> CircuitData::getComponentByInputDataIndex(const int componentIndex)
{
	try {
		return mapOfComponents.at(componentIndex);
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw;
	}
}

boost::shared_ptr<CircuitPressureNode> CircuitData::getNodeByInputDataIndex(const int componentIndex)
{
	try {
		return mapOfPressureNodes.at(componentIndex);
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw;
	}
}

std::vector<std::pair<int,double*>> CircuitData::getComponentInputDataIndicesAndFlows() const
{
	std::vector<std::pair<int,double*>> returnValue;
	for (auto componentPair = mapOfComponents.begin(); componentPair != mapOfComponents.end(); componentPair++)
	{
		returnValue.push_back(std::make_pair(componentPair->first, componentPair->second->getFlowPointer()));
	}
	return returnValue;
}

std::vector<std::pair<int,double*>> CircuitData::getNodeInputDataIndicesAndPressures() const
{
	std::vector<std::pair<int,double*>> returnValue;
	for (auto nodePair = mapOfPressureNodes.begin(); nodePair != mapOfPressureNodes.end(); nodePair++)
	{
		returnValue.push_back(std::make_pair(nodePair->first, nodePair->second->getPressurePointer()));
	}
	return returnValue;
}

std::vector<std::pair<int,double*>> CircuitData::getVolumeTrackingComponentInputDataIndicesAndVolumes() const
{
	std::vector<std::pair<int,double*>> returnValue;
	for (auto componentPair = mapOfVolumeTrackingComponents.begin(); componentPair != mapOfVolumeTrackingComponents.end(); componentPair++)
	{
		boost::shared_ptr<VolumeTrackingComponent> downcastPressureChamber = boost::dynamic_pointer_cast<VolumeTrackingComponent> (componentPair->second);
		returnValue.push_back(std::make_pair(componentPair->first, downcastPressureChamber->getVolumePointer()));
	}
	return returnValue;
}

std::vector<double*> CircuitData::getCapacitorNodalHistoryPressurePointers() const
{
	std::vector<double*> capacitorNodalHistoryPressurePointers;

	for (auto componentMapEntry : mapOfComponents)
	{
		CircuitComponent component = *(componentMapEntry.second);
		if (component.getType() == Component_Capacitor)
		{
			boost::shared_ptr<CircuitPressureNode> node = component.getStartNode();
			// Fixed pressures should not be allowed to be varied by the Kalman filter
			// (for which this method was written)
			if (node->getPressurePrescriptionType() != Pressure_Fixed)
			{
				capacitorNodalHistoryPressurePointers.push_back(node->getPressurePointer());
			}
			
		}
	}
	return capacitorNodalHistoryPressurePointers;
}

boost::shared_ptr<std::vector<std::pair<parameter_controller_t, int>>> CircuitData::getControlTypesAndComponentIndices() const
{
	boost::shared_ptr<std::vector<std::pair<parameter_controller_t, int>>> controlTypesAndComponentIndices(new std::vector<std::pair<parameter_controller_t, int>>());
    for (auto& component : components)
    {
        if (component->hasUserDefinedExternalPythonScriptParameterController())
        {
        	controlTypesAndComponentIndices->push_back(std::make_pair(component->getControlType(), component->getIndex()));
        }
    }
    return controlTypesAndComponentIndices;
}

bool Netlist3DDomainReplacementCircuitData::hasPrescribedFlowAcrossInterface() const
{
	// Negate and return:
	return !m_flowPermittedAcross3DInterface;
}

bool Netlist3DDomainReplacementCircuitData::hasPrescribedPressureAcrossInterface() const
{
	return m_flowPermittedAcross3DInterface;
}

bool Netlist3DDomainReplacementCircuitData::isNodeAtBoundaryInterface(int nodeIndex)
{
	bool returnValue = false;
	for (auto boundaryNodeIndex = m_indexOfNodeAt3DInterface.begin(); boundaryNodeIndex != m_indexOfNodeAt3DInterface.end(); boundaryNodeIndex++)
	{
		if (nodeIndex == *boundaryNodeIndex)
		{
			returnValue = true;
		}
	}
	return returnValue;
}

void Netlist3DDomainReplacementCircuitData::initialiseNodesAndComponentsAtInterface_vector(std::vector<int> threeDInterfaceNodeIndices)
{
	setIndicesOfNodesAtInterface(threeDInterfaceNodeIndices);
	// tag the nodes at the 3D interface:
	for (auto threeDInterfaceNodeIndex = threeDInterfaceNodeIndices.begin(); threeDInterfaceNodeIndex != threeDInterfaceNodeIndices.end(); threeDInterfaceNodeIndex++)
	{
    	try {
    		mapOfPressureNodes.at(*threeDInterfaceNodeIndex)->setIsAtBoundary();
    	} catch (const std::exception& e) {
    	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
    	    throw;
    	}
    }

	int numberOfComponentsTaggedAsBeingAtInterfaces = 0; //a counter to verify we get m_numberOfNetlistsUsedAsBoundaryConditions of these.
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		// Find the component tagged to recieve flow from the 3D interface:
		// if ((*component)->prescribedFlowType == Flow_3DInterface)
		// {
		// Ensure one of this component's nodes is the one at the 3D interface
		bool startNodeIsAt3Dinterface = isNodeAtBoundaryInterface((*component)->startNode->getIndex());
		bool endNodeIsAt3Dinterface = isNodeAtBoundaryInterface((*component)->endNode->getIndex());
		if (startNodeIsAt3Dinterface || endNodeIsAt3Dinterface)
		{
			// tag the component as being at the 3D interface
			(*component)->setConnectsToNodeAtInterface();
			numberOfComponentsTaggedAsBeingAtInterfaces++;
		}

		// Set up the sign for the flow prescription, depending on whether it's the start or end node of the
		// component which is at the 3D interface (if the user has listed the end node before the start node for
		// the component at the 3D interface in the netlist_surfaces.dat, we need to flip the sign of the flow
		// coming from the 3D interface. This affects all component types, but is needed for diodes in particular,
		// where the user has a legitimate reason to want to reverse those nodes, depending on the orientation they
		// want to give to the diode.)
		if (startNodeIsAt3Dinterface)
		{
			m_signForPrescribed3DInterfaceFlow = -1.0;
		}
		else if (endNodeIsAt3Dinterface)
		{
			m_signForPrescribed3DInterfaceFlow = 1.0;
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
		// 	// (*component)->startNode->m_isAtBoundary = true;
		// 	(*component)->setConnectsToNodeAtInterface();
		// }
	}
	if (numberOfComponentsTaggedAsBeingAtInterfaces!=m_numberOfNetlistsUsedAsBoundaryConditions)
	{
		// May fail if this subroutine gets called on a subcircuit, instead of the whole input data circuit, if the subcircuit doesn't have the 3D flow interface component. Don't call it in this case!
		std::stringstream errorMessage;
		errorMessage << "EE: Expected " << m_numberOfNetlistsUsedAsBoundaryConditions << " components at the boundary interfaces in the 0D replacement" << std::endl;
		errorMessage << "               for the 3D domain, but found " << numberOfComponentsTaggedAsBeingAtInterfaces << ". You probably can't fix this." << std::endl;
		throw std::logic_error(errorMessage.str());
	}

	// Ensure that the 3D interface nodes belongs to unique components:
	int numberOfComponentsConnectingTo3DInterfaceNode = 0; //a counter to verify there exists a unique component connecting to the 3D interface node
	for (auto component=components.begin(); component!=components.end(); component++)
	{
		bool startNodeIsAt3Dinterface = isNodeAtBoundaryInterface((*component)->startNode->getIndex());
		bool endNodeIsAt3Dinterface = isNodeAtBoundaryInterface((*component)->endNode->getIndex());
		if (startNodeIsAt3Dinterface)
		{
			numberOfComponentsConnectingTo3DInterfaceNode++;
		}
		if (endNodeIsAt3Dinterface)
		{
			numberOfComponentsConnectingTo3DInterfaceNode++;
		}
	}
	if (numberOfComponentsConnectingTo3DInterfaceNode!=m_numberOfNetlistsUsedAsBoundaryConditions)
	{
		// May fail if this subroutine gets called on a subcircuit, instead of the whole input data circuit, if the subcircuit doesn't have the 3D flow interface component. Don't call it in this case!
		std::stringstream errorMessage;
		errorMessage << "EE: Expected exactly one appearance of the 3D interface node in the circuit data, but found " << numberOfComponentsConnectingTo3DInterfaceNode << "." << std::endl;
		throw std::runtime_error(errorMessage.str());
	}

}

// Adjusts the circuit data as appropriate for a change in boundary condition type
// between Neumann and Dirichlet.
void Netlist3DDomainReplacementCircuitData::setBoundaryPrescriptionsAndBoundaryConditionTypes(std::vector<std::pair<boundary_data_t,double>>& boundaryFlowsOrPressuresAsAppropriate)
{
	// std::cout << "map of prescribed flow components (pre): " << std::endl;
	// for (auto flowcomp=mapOfPrescribedFlowComponents.begin(); flowcomp!=mapOfPrescribedFlowComponents.end(); flowcomp++)
	// {
	// 	std::cout << flowcomp->first << " " << flowcomp->second->indexInInputData << std::endl;
	// }

	// std::cout << "map of prescribed pressure components (pre): " << std::endl;
	// for (auto flowcomp=mapOfPrescribedPressureNodes.begin(); flowcomp!=mapOfPrescribedPressureNodes.end(); flowcomp++)
	// {
	// 	std::cout << flowcomp->first << " " << flowcomp->second->indexInInputData << std::endl;
	// }

	for (int componentAtBoundaryIndex = 0; componentAtBoundaryIndex<m_numberOfNetlistsUsedAsBoundaryConditions; componentAtBoundaryIndex++)
	{
		// Scoping unit to avoid errors:
		{
			bool flowIsGivenTo0DReplacementDomainAtThisBoundary;
			try {
				flowIsGivenTo0DReplacementDomainAtThisBoundary = (boundaryFlowsOrPressuresAsAppropriate.at(componentAtBoundaryIndex).first == Boundary_Flow);
			} catch (const std::exception& e) {
			    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
			    throw;
			}
			if (flowIsGivenTo0DReplacementDomainAtThisBoundary)
			{
				try {
					givePrescribedFlowToBoundaryComponent(toOneIndexing(componentAtBoundaryIndex),boundaryFlowsOrPressuresAsAppropriate.at(componentAtBoundaryIndex).second);
				} catch (const std::exception& e) {
				    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
				    throw;
				}
			}
			bool flowWasPreviouslyNotGivenTo0DReplacementDomainAtThisBoundary = (mapOfPrescribedFlowComponents.find(toOneIndexing(componentAtBoundaryIndex)) == mapOfPrescribedFlowComponents.end());
			bool flowNewlyPrescribedAtThisBoundary = (flowIsGivenTo0DReplacementDomainAtThisBoundary && flowWasPreviouslyNotGivenTo0DReplacementDomainAtThisBoundary);
			if (flowNewlyPrescribedAtThisBoundary)
			{
				std::cout << "Switching to Dirichlet in zero-D domain replacement!" << std::endl;
				// Remove the node at the interface from the list of those with prescribed pressure:
				for (auto node=mapOfPrescribedPressureNodes.begin(); node!=mapOfPrescribedPressureNodes.end(); node++)
				{
					if (node->first == toOneIndexing(componentAtBoundaryIndex))
					{
						node->second->setPressurePrescriptionType(Pressure_NotPrescribed);
						mapOfPrescribedPressureNodes.erase(node);
						break;
					}
				}

				// Add the component at the interface to the list of those with prescribed flow:
				for (auto component=mapOfComponents.begin(); component!=mapOfComponents.end(); component++)
				{
					if (component->first == toOneIndexing(componentAtBoundaryIndex))
					{
						component->second->prescribedFlowType=Flow_3DInterface;
						mapOfPrescribedFlowComponents.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> (component->second->getIndex(),component->second));
						break;
					}
				}

				// Update the counts of each type of prescription (in switching to a Neumann condition, we're now
				// giving a prescribed pressure to the domain, which means this boundary condition will be
				// receiving a flow back from the domain to prescribe at the LPN interface, instead of the
				// pressure it was previously receiving under Dirichlet conditions. These counters therefore
				// need updating.)
				numberOfPrescribedPressures--;
				numberOfPrescribedFlows++;
				assert(numberOfPrescribedPressures>=0);
			}
		}
		// Scoping unit to avoid errors:
		{
			bool pressureIsGivenTo0DReplacementDomainAtThisBoundary;
			try {
				pressureIsGivenTo0DReplacementDomainAtThisBoundary = (boundaryFlowsOrPressuresAsAppropriate.at(componentAtBoundaryIndex).first == Boundary_Pressure);
			} catch (const std::exception& e) {
			    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
			    throw;
			}
			if(pressureIsGivenTo0DReplacementDomainAtThisBoundary)
			{
				try {
					givePrescribedPressureToBoundaryNode(toOneIndexing(componentAtBoundaryIndex),boundaryFlowsOrPressuresAsAppropriate.at(componentAtBoundaryIndex).second);
				} catch (const std::exception& e) {
			    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
			    	throw;
		    	}
			}
			bool pressureWasPreviouslyNotGivenTo0DReplacementDomainAtThisBoundary = (mapOfPrescribedPressureNodes.find(toOneIndexing(componentAtBoundaryIndex)) == mapOfPrescribedPressureNodes.end());
			bool pressureNewlyPrescribedAtThisBoundary = (pressureIsGivenTo0DReplacementDomainAtThisBoundary && pressureWasPreviouslyNotGivenTo0DReplacementDomainAtThisBoundary);
			if (pressureNewlyPrescribedAtThisBoundary)
			{
				std::cout << "Switching to Neumann in zero-D domain replacement!" << std::endl;
				// Add the node at the interface to the list of those with prescribed pressure:
				for (auto node=mapOfPressureNodes.begin(); node!=mapOfPressureNodes.end(); node++)
				{
					if (node->first == toOneIndexing(componentAtBoundaryIndex))
					{
						node->second->setPressurePrescriptionType(Pressure_3DInterface);
						mapOfPrescribedPressureNodes.insert(std::pair<int,boost::shared_ptr<CircuitPressureNode>> (node->second->getIndex(),node->second));
						break;
					}
				}

				// Remove the component at the interface from the list of those with prescribed flow:
				for (auto component=mapOfPrescribedFlowComponents.begin(); component!=mapOfPrescribedFlowComponents.end(); component++)
				{
					if (component->first == toOneIndexing(componentAtBoundaryIndex))
					{
						component->second->prescribedFlowType=Flow_NotPrescribed;
						mapOfPrescribedFlowComponents.erase(component);
						break;
					}
				}

				// Update the counts of each type of prescription (in switching to a Neumann condition, we're now
				// giving a prescribed flow to the domain, which means this boundary condition will be
				// receiving a pressure back from the domain to prescribe at the LPN interface, instead of the
				// flow it was previously receiving under Neumann conditions. These counters therefore
				// need updating.)
				numberOfPrescribedPressures++;
				numberOfPrescribedFlows--;
				assert(numberOfPrescribedFlows>=0);
			}
		}

	}

	// std::cout << "map of prescribed flow components (post): " << std::endl;
	// for (auto flowcomp=mapOfPrescribedFlowComponents.begin(); flowcomp!=mapOfPrescribedFlowComponents.end(); flowcomp++)
	// {
	// 	std::cout << flowcomp->first << " " << flowcomp->second->indexInInputData << std::endl;
	// }

	// std::cout << "map of prescribed pressure components (post): " << std::endl;
	// for (auto flowcomp=mapOfPrescribedPressureNodes.begin(); flowcomp!=mapOfPrescribedPressureNodes.end(); flowcomp++)
	// {
	// 	std::cout << flowcomp->first << " " << flowcomp->second->indexInInputData << std::endl;
	// }
}

void Netlist3DDomainReplacementCircuitData::givePrescribedPressureToBoundaryNode(int nodeIndex, double prescribedPressure)
{
	try {
		mapOfPressureNodes.at(nodeIndex)->setPressure(prescribedPressure);
	} catch (const std::exception& e) {
    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
    	throw;
    }
}

void Netlist3DDomainReplacementCircuitData::givePrescribedFlowToBoundaryComponent(int componentIndex, double prescribedFlow)
{
	try {
		mapOfComponents.at(componentIndex)->flow = prescribedFlow;
	} catch (const std::exception& e) {
    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
    	throw;
    }
}

boost::shared_ptr<CircuitComponent> Netlist3DDomainReplacementCircuitData::getDpDqResistorByIndex(int index)
{
	return m_mapOfDpDqResistors.at(index);
}

void Netlist3DDomainReplacementCircuitData::addToMapOfDpDqResistors(int indexOfResistor, boost::shared_ptr<CircuitComponent> dpDqResistor)
{
	m_mapOfDpDqResistors.insert(std::make_pair(indexOfResistor,dpDqResistor));
}

bool Netlist3DDomainReplacementCircuitData::isADpDqResistor(const int componentIndex)
{
	return (m_mapOfDpDqResistors.count(componentIndex) == 1);
}