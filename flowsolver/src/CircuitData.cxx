#include "CircuitData.hxx"
#include "indexShifters.hxx"
#include <stdexcept>
#include <float.h>
#include <algorithm>
#include <stack>
#include <sstream>
#include <boost/make_shared.hpp>

#include "common_c.h"


bool CircuitComponent::hasNonnegativePressureGradientOrForwardFlow() // whether the diode should be open
{
	// We use a tolerance on these floating-point comparisons here
	// because the errors in the netlist linear system solves
	// for the boundary conditions can be of order 1e-9.
	const double floatingPointTolerance = 1e-8;
	// std::cout << "start and end node pressures for switching in hasNonnegativePressureGradientOrForwardFlow: " << startNode->getPressure() << " " << endNode->getPressure() << " difference is: " << startNode->getPressure() - endNode->getPressure() << std::endl;
	// std::cout << "m_connectsToNodeAtInterface: " << m_connectsToNodeAtInterface << std::endl;
	bool hasNonnegativePressureGradient = (startNode->getPressure() - endNode->getPressure() >= 0.0 - floatingPointTolerance);
	bool hasForwardFlow;
	// Diode closure is enforced by setting diode resistance to DBL_MAX, so there remains a small flow on the order 1e-308 across a closed diode.
	hasForwardFlow = flow >= floatingPointTolerance;

	return (hasNonnegativePressureGradient || hasForwardFlow);
}

boost::shared_ptr<CircuitPressureNode> CircuitComponent::getStartNode()
{
	return startNode;
}

boost::shared_ptr<CircuitPressureNode> CircuitComponent::getEndNode()
{
	return endNode;
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

	setupComponentNeighbourPointers();

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
    	throw e;
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
	    	throw e;
	    }
		
		// Return this existing node:
		boost::shared_ptr<CircuitPressureNode> returnValue;
		try {
			returnValue = mapOfPressureNodes.at(indexInInputData_in);
		} catch (const std::exception& e) {
	    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    	throw e;
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

bool CircuitComponent::hasUserDefinedExternalPythonScriptParameterController() const
{
	return m_hasPythonParameterController;
}

std::string CircuitComponent::getPythonControllerName(const parameter_controller_t controllerType) const
{
	assert(m_hasPythonParameterController);
	return m_pythonParameterControllerNames.at(controllerType);
}

// Sets the name "whateverName" of the parameter controller to look for in the
// working directory: whateverName.py, containing class whateverName,
// with class method:
// newParamterValue = updateControl(self, oldParameterValue, delt).
//
// This should be a Python script.
void CircuitComponent::addPythonControllerName(const parameter_controller_t controllerType, const std::string pythonParameterControllerName)
{
	m_hasPythonParameterController = true;
	m_pythonParameterControllerNames.insert(std::make_pair(controllerType, pythonParameterControllerName));
}


bool CircuitComponent::permitsFlow() const
{
	return m_permitsFlow;
}

void CircuitComponent::enableDiodeFlow()
{
	assert(m_type == Component_Diode);
	m_currentParameterValue = parameterValueFromInputData; // For enforcing zero resistance when the diode is open
	m_permitsFlow = true;
}

void CircuitComponent::disableDiodeFlow()
{
	assert(m_type == Component_Diode);
	m_currentParameterValue = DBL_MAX; // For enforcing "infinite" resistance when the diode is open//NOT USED
	m_permitsFlow = false;
}

double* CircuitComponent::getFlowPointer()
{
	double* flowPointer = &flow;
	return flowPointer;
}

void CircuitComponent::setRestartFlowFromHistory()
{
	flow = m_entireFlowHistory.back();
}

double* CircuitPressureNode::getPressurePointer()
{
	double* pressurePointer = &pressure;
	return pressurePointer;
}

double* CircuitPressureNode::getPointerToFixedPressurePrescription()
{
	// we should only be accessing this pointer for modification if it is a fixed-type pressure prescription
	assert(m_prescribedPressureType == Pressure_Fixed);
	double* pressurePointer = &m_fixedPressure;
	return pressurePointer;
}

bool CircuitPressureNode::hasUserDefinedExternalPythonScriptParameterController() const
{
	return m_hasPythonParameterController;
}

bool CircuitComponent::hasPrescribedFlow() const
{
	return m_hasPrescribedFlow;
}

std::string CircuitPressureNode::getPythonControllerName() const
{
	assert(m_hasPythonParameterController);
	return m_pythonParameterControllerName;
}

// Sets the name "whateverName" of the nodal pressure controller to look for in the
// working directory: whateverName.py, containing class whateverName,
// with class method:
// newParamterValue = updateControl(self, oldParameterValue, delt).
//
// This should be a Python script.
void CircuitPressureNode::setPythonControllerName(const std::string pythonParameterControllerName)
{
	m_hasPythonParameterController = true;
	m_pythonParameterControllerName = pythonParameterControllerName;
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

int CircuitComponent::getIndex() const
{
	return m_indexInInputData;
}

void CircuitComponent::setIndex(const int index)
{
	m_indexInInputData = index;
}

double CircuitComponent::getPrescribedFlow() const
{
	assert(m_hasPrescribedFlow);
	return m_valueOfPrescribedFlow;
}

void CircuitComponent::setPrescribedFlow(const double prescribedFlow)
{
	m_valueOfPrescribedFlow = prescribedFlow;
	m_hasPrescribedFlow = true;
}

double* CircuitComponent::getPointerToFixedFlowPrescription()
{
	assert(m_hasPrescribedFlow);
	return &m_valueOfPrescribedFlow;
}

void CircuitComponent::setHasHistoryVolume(const bool hasHistoryVolume)
{
	m_hasHistoryVolume = hasHistoryVolume;
}

bool CircuitComponent::getHasHistoryVolume()
{
	return m_hasHistoryVolume;
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
	    throw e;
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
		    throw e;
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
	    throw e;
	}
}

circuit_component_t& CircuitComponent::getType()
{
	return m_type;
}

double* CircuitComponent::getParameterPointer()
{
	return &m_currentParameterValue;
}

void CircuitComponent::setParameterValue(double const parameterValue)
{
	m_currentParameterValue = parameterValue;
}

bool CircuitComponent::connectsToNodeAtInterface()
{
	return m_connectsToNodeAtInterface;
}

void CircuitComponent::setConnectsToNodeAtInterface()
{
	m_connectsToNodeAtInterface = true;
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
	    throw e;
	}
}

boost::shared_ptr<CircuitPressureNode> CircuitData::getNodeByInputDataIndex(const int componentIndex)
{
	try {
		return mapOfPressureNodes.at(componentIndex);
	} catch (const std::exception& e) {
	    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	    throw e;
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
    	    throw e;
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
			    throw e;
			}
			if (flowIsGivenTo0DReplacementDomainAtThisBoundary)
			{
				try {
					givePrescribedFlowToBoundaryComponent(toOneIndexing(componentAtBoundaryIndex),boundaryFlowsOrPressuresAsAppropriate.at(componentAtBoundaryIndex).second);
				} catch (const std::exception& e) {
				    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
				    throw e;
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
			    throw e;
			}
			if(pressureIsGivenTo0DReplacementDomainAtThisBoundary)
			{
				try {
					givePrescribedPressureToBoundaryNode(toOneIndexing(componentAtBoundaryIndex),boundaryFlowsOrPressuresAsAppropriate.at(componentAtBoundaryIndex).second);
				} catch (const std::exception& e) {
			    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
			    	throw e;
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
    	throw e;
    }
}

void Netlist3DDomainReplacementCircuitData::givePrescribedFlowToBoundaryComponent(int componentIndex, double prescribedFlow)
{
	try {
		mapOfComponents.at(componentIndex)->flow = prescribedFlow;
	} catch (const std::exception& e) {
    	std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
    	throw e;
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

int CircuitPressureNode::getIndex() const
{
	return m_indexInInputData;
}

void CircuitPressureNode::setIsAtBoundary()
{
	m_isAtBoundary = true;
}

bool CircuitPressureNode::isAtBoundary() const
{
	return m_isAtBoundary;
}

double CircuitPressureNode::getPressure()
{
	// If this is a prescribed fixed pressure, ensure we reset it to the original input value.
	// This has the additional benefit of stopping any drift in a supposedly-prescribed value.
	if (m_prescribedPressureType == Pressure_Fixed)
	{
		pressure = m_fixedPressure;
	}
	return pressure;
}

void CircuitPressureNode::setPressure(const double pressure_in)
{
	pressure = pressure_in;
}

void CircuitPressureNode::setPrescribedPressure(const double prescribedPressure)
{
	// We only do anything special with fixed-pressure values. There's nothing
	// to do here with other types of prescribed pressure, as they don't
	// remain fixed at a single value (so we needn't remember it in m_fixedPressure).
	if (m_prescribedPressureType == Pressure_Fixed)
	{
		m_fixedPressure = prescribedPressure;
	}
	else
	{
		pressure = prescribedPressure;
	}
}

void CircuitPressureNode::setRestartPressureFromHistory()
{
	pressure = m_entirePressureHistory.back();
}

void CircuitPressureNode::setHasHistoryPressure(const bool hasHistoryPressure)
{
	m_hasHistoryPressure = hasHistoryPressure;
}

void CircuitPressureNode::copyPressureToHistoryPressure()
{
	std::cout << "setting history pressure to " << getPressure() << " from " << m_historyPressure << std::endl;
	m_historyPressure = getPressure();
}

void CircuitPressureNode::copyHistoryPressureToHistoryHistoryPressure()
{
	// used with the Kalman filter to store a history pressure from the step before the particle
	m_historyHistoryPressure = m_historyPressure;
}

double CircuitPressureNode::getHistoryPressure() const
{
	return m_historyPressure;
}

double CircuitPressureNode::getHistoryHistoryPressure() const
{
	// used with the Kalman filter to store a history pressure from the step before the particle
	return m_historyHistoryPressure;
}

bool CircuitPressureNode::hasHistoryPressure() const
{
	return m_hasHistoryPressure;
}

circuit_nodal_pressure_prescription_t CircuitPressureNode::getPressurePrescriptionType() const
{
	return m_prescribedPressureType;
}

void CircuitPressureNode::setPressurePrescriptionType(const circuit_nodal_pressure_prescription_t prescribedPressureType)
{
	m_prescribedPressureType = prescribedPressureType;
}

int CircuitPressureNode::getPrescribedPressurePointerIndex() const
{
	return m_prescribedPressurePointerIndex;
}

void CircuitPressureNode::setPrescribedPressurePointerIndex(const int prescribedPressurePointerIndex)
{
	m_prescribedPressurePointerIndex = prescribedPressurePointerIndex;
}

double CircuitPressureNode::getFromPressureHistoryByTimestepIndex(const int timestepIndex) const
{
	return m_entirePressureHistory.at(timestepIndex);
}

void CircuitPressureNode::appendToPressureHistory(const double pressure)
{
	m_entirePressureHistory.push_back(pressure);
}

void VolumeTrackingComponent::recordVolumeInHistory()
{
	m_entireVolumeHistory.push_back(m_storedVolume);
}

double VolumeTrackingComponent::getVolumeHistoryAtTimestep(int timestep)
{
	return m_entireVolumeHistory.at(timestep);
}

void VolumeTrackingComponent::setVolumeHistoryAtTimestep(double historyVolume)
{
	m_entireVolumeHistory.push_back(historyVolume);
}

void VolumeTrackingComponent::setStoredVolume(const double newVolume)
{
	m_storedVolume = newVolume;
}
// The /proposed/ volume is the one which gets checked for negative (invalid) values
// so that we can detect such invalid cases, and take steps to remedy.
void VolumeTrackingComponent::setProposedVolume(const double proposedVolume)
{
	m_proposedVolume = proposedVolume;
}
double VolumeTrackingComponent::getVolume()
{
	return m_storedVolume;
}
double* VolumeTrackingComponent::getVolumePointer()
{
	return &m_storedVolume;
}
double VolumeTrackingComponent::getProposedVolume()
{
	std::cout<<"proposed volume was: " << m_proposedVolume << std::endl;
	return m_proposedVolume;	
}
double VolumeTrackingComponent::getHistoryVolume()
{
	return m_historyVolume;
}
void VolumeTrackingComponent::cycleHistoryVolume()
{
	m_historyVolume = m_storedVolume;
}

double VolumeTrackingComponent::getElastance()
{
	return m_currentParameterValue;
}

bool VolumeTrackingComponent::zeroVolumeShouldBePrescribed()
{
	return m_enforceZeroVolumePrescription;
}

void VolumeTrackingComponent::enforceZeroVolumePrescription()
{
	m_enforceZeroVolumePrescription = true;
}

void VolumeTrackingComponent::resetZeroVolumePrescription()
{
	m_enforceZeroVolumePrescription = false;
}

void VolumeTrackingComponent::setRestartVolumeFromHistory()
{
	m_storedVolume = m_entireVolumeHistory.back();
}

void VolumeTrackingPressureChamber::passPressureToStartNode()
{
	m_pressure = (m_storedVolume - m_unstressedVolume)*m_currentParameterValue;
	std::cout << "m_unstressedVolume: " << m_unstressedVolume << std::endl;
	std::cout << "compliance set to: " << m_currentParameterValue << std::endl;
	std::cout << "pressure set to: " << m_pressure << std::endl;
	startNode->setPressure(m_pressure);
}

double VolumeTrackingPressureChamber::getUnstressedVolume() 
{
	return m_unstressedVolume;
}

void VolumeTrackingPressureChamber::setStoredVolume(const double newVolume)
{
	std::cout << "pre m_storedVolume: " << m_storedVolume << std::endl;
	m_storedVolume = newVolume;
	std::cout << "post m_storedVolume: " << m_storedVolume << std::endl;
	passPressureToStartNode();
}

double* VolumeTrackingPressureChamber::getUnstressedVolumePointer()
{
	return &m_unstressedVolume;
}