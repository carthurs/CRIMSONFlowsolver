#include "ControlSystemsManager.hxx"
#include <sstream>

void ControlSystemsManager::updateAllControlSystems()
{
	for (auto controlSystem = m_controlSystems.begin(); controlSystem != m_controlSystems.end(); controlSystem++)
	{
		(*controlSystem)->updateControl();
	}
}

void ControlSystemsManager::createParameterController(const parameter_controller_t controllerType, const boost::shared_ptr<NetlistCircuit> boundaryCondition, const int nodeOrComponentIndex)
{
	switch (controllerType)
	{
		case  Controller_LeftVentricularElastance:
			{
				// We know that we must be working with a pressure chammber, because Controller_LeftVentricularElastance
				// refers to a controller for a VolumeTrackingPressureChamber, which derives from CircuitComponent
				boost::shared_ptr<CircuitComponent> component = boundaryCondition->getComponentByInputDataIndex(nodeOrComponentIndex);
				// Get the pointer to the compliance which needs to be controlled (in this case, the compliance of the pressure chamber):
				boost::shared_ptr<VolumeTrackingPressureChamber> downcastVolumeComponent = boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (component);

				// Make sure the dynamic cast was successful:
				assert(downcastVolumeComponent!=NULL);

				double* parameterToControl = downcastVolumeComponent->getPointerToElastance();

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new LeftVentricularElastanceController(parameterToControl,m_delt));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_BleedResistance:
			{
				// get the resistor:
				boost::shared_ptr<CircuitComponent> resistor = boundaryCondition->getComponentByInputDataIndex(nodeOrComponentIndex);
				assert(resistor->getType() == Component_Resistor);

				double* resistanceToControl = resistor->getParameterPointer();

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new BleedController(resistanceToControl));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_BleedCompliance:
			{
				// get the capacitor:
				boost::shared_ptr<CircuitComponent> capacitor = boundaryCondition->getComponentByInputDataIndex(nodeOrComponentIndex);
				assert(capacitor->getType() == Component_Capacitor);

				double* complianceToControl = capacitor->getParameterPointer();

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new BleedController(complianceToControl));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_CustomPython:
			{
				// get the component or node:
				boost::shared_ptr<CircuitComponent> controlledComponent = boundaryCondition->getComponentByInputDataIndex(nodeOrComponentIndex);
				double* parameterToControl = controlledComponent->getParameterPointer();
				std::string externalPythonControllerName;

				if (controlledComponent->hasUserDefinedExternalPythonScriptParameterController())
				{
					externalPythonControllerName = controlledComponent->getPythonControllerName();
				}
				else
				{
					std::stringstream errorMessage;
					errorMessage << "EE: A component or node of the Netlist circuit at surface " << boundaryCondition->getSurfaceIndex() << 
					errorMessage << " was tagged as having an external Python parameter controller, but none was found." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new UserDefinedCustomPythonParameterController(parameterToControl, m_delt, externalPythonControllerName));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		default:
			std::stringstream errorMessage;
			errorMessage << "EE: Unknown control parameter type (controllerType) requested for Netlist boundary condition for surface " << boundaryCondition->getSurfaceIndex() << "." << std::endl;
			throw std::runtime_error(errorMessage.str());
	}
}