#include "ControlSystemsManager.hxx"
#include <sstream>

void ControlSystemsManager::updateAllControlSystems()
{
	for (auto controlSystem = m_controlSystems.begin(); controlSystem != m_controlSystems.end(); controlSystem++)
	{
		(*controlSystem)->updateControl();
	}
}

void ControlSystemsManager::createParameterController(const parameter_controller_t controllerType, const boost::shared_ptr<NetlistBoundaryCondition> boundaryCondition, const int nodeOrComponentIndex)
{
	if (controllerType == Controller_LeftVentricularElastance)
	{
		// We know that we must be working with a pressure chammber, because Controller_LeftVentricularElastance
		// refers to a controller for a VolumeTrackingPressureChamber, which derives from CircuitComponent
		boost::shared_ptr<CircuitComponent> component = boundaryCondition->getCircuitDescription().mapOfComponents.at(nodeOrComponentIndex);
		// Get the pointer to the compliance which needs to be controlled (in this case, the compliance of the pressure chamber):
		boost::shared_ptr<VolumeTrackingPressureChamber> downcastVolumeComponent = boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (component);

		// Make sure the dynamic cast was successful:
		assert(downcastVolumeComponent!=NULL);

		double* parameterToControl = downcastVolumeComponent->getPointerToElastance();

		boost::shared_ptr<AbstractParameterController> controllerToPushBack(new LeftVentricularElastanceController(parameterToControl,m_delt));
		m_controlSystems.push_back(controllerToPushBack);
	}
	else
	{
		std::stringstream errorMessage;
		errorMessage << "EE: Unknown control parameter type (controllerType) requested for Netlist boundary condition for surface " << boundaryCondition->getSurfaceIndex() << "." << std::endl;
		throw std::runtime_error(errorMessage.str());
	}
}