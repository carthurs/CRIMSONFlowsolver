#include "CircuitComponent.hxx"
#include "CircuitData.hxx"

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

bool CircuitComponent::hasUserDefinedExternalPythonScriptParameterController() const
{
	return m_hasPythonParameterController;
}

std::string CircuitComponent::getPythonControllerName(const parameter_controller_t controllerType) const
{
	assert(m_hasPythonParameterController);
	return m_pythonParameterControllerNames.at(controllerType);
}

parameter_controller_t CircuitComponent::getControlType() const
{
	// There appears to be some old functionality for multiple controllers. I doubt this is actually used, but assert it anyway.
	assert(m_pythonParameterControllerNames.size() == 1);
	return m_pythonParameterControllerNames.begin()->first;
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

void CircuitComponent::appendToFlowHistory(const double flow)
{
	// m_flowHistoryBuffer.at(m_flowHistoryBufferNextWriteIndex) = flow;
	// m_flowHistoryBufferNextWriteIndex++;
	// if (m_flowHistoryBufferNextWriteIndex > m_flowHistoryBufferSize - 1)
	// {
	// 	m_flowHistoryBufferNextWriteIndex = 0;
	// 	m_entireFlowHistory.insert(m_entireFlowHistory.end(), m_flowHistoryBuffer.begin(), m_flowHistoryBuffer.end());
	// }

	m_entireFlowHistory.push_back(flow);
}

double CircuitComponent::getFromFlowHistoryByTimestepIndex(const int timestepIndex) const
{
	return m_entireFlowHistory.at(timestepIndex);
}

int CircuitComponent::getFlowHistoryLength() const
{
	return m_entireFlowHistory.size();
}

bool CircuitComponent::hasPrescribedFlow() const
{
	return m_hasPrescribedFlow;
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
	// std::cout << "just set flow " << prescribedFlow << " for component " << m_indexInInputData << std::endl;
}

double* CircuitComponent::getPointerToFixedFlowPrescription()
{
	assert(m_hasPrescribedFlow);
	return &m_valueOfPrescribedFlow;
}

void CircuitComponent::setHasNoPrescribedFlow()
{
	m_hasPrescribedFlow = false;
}

void CircuitComponent::setHasHistoryVolume(const bool hasHistoryVolume)
{
	m_hasHistoryVolume = hasHistoryVolume;
}

bool CircuitComponent::getHasHistoryVolume()
{
	return m_hasHistoryVolume;
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

void VolumeTrackingPressureChamber::passPressureToStartNode()
{
	m_pressure = (m_storedVolume - m_unstressedVolume)*m_currentParameterValue;
	// std::cout << "m_unstressedVolume: " << m_unstressedVolume << std::endl;
	// std::cout << "compliance set to: " << m_currentParameterValue << std::endl;
	// std::cout << "pressure set to: " << m_pressure << std::endl;
	startNode->setPressure(m_pressure);
}

double VolumeTrackingPressureChamber::getUnstressedVolume() 
{
	return m_unstressedVolume;
}

void VolumeTrackingPressureChamber::setStoredVolume(const double newVolume)
{
	// std::cout << "pre m_storedVolume: " << m_storedVolume << std::endl;
	m_storedVolume = newVolume;
	// std::cout << "post m_storedVolume: " << m_storedVolume << std::endl;
	passPressureToStartNode();
}

double* VolumeTrackingPressureChamber::getUnstressedVolumePointer()
{
	return &m_unstressedVolume;
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