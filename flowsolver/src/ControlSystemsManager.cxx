#include "ControlSystemsManager.hxx"
#include <sstream>
#include <cstdlib>
#include <boost/filesystem.hpp>

void ControlSystemsManager::updateAllControlSystems()
{
	for (auto controlSystem = m_controlSystems.begin(); controlSystem != m_controlSystems.end(); controlSystem++)
	{
		(*controlSystem)->updateControl();
	}
}

void ControlSystemsManager::createParameterController(const parameter_controller_t controllerType, const boost::shared_ptr<NetlistCircuit> netlistCircuit, const int nodeOrComponentIndex)
{
	switch (controllerType)
	{
		case  Controller_LeftVentricularElastance:
			{
				// We know that we must be working with a pressure chammber, because Controller_LeftVentricularElastance
				// refers to a controller for a VolumeTrackingComponent, which derives from CircuitComponent
				boost::shared_ptr<CircuitComponent> component = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				
				// Make sure that the user actually requested this controller on a VolumeTrackingComponent:
				if (boost::dynamic_pointer_cast<VolumeTrackingComponent> (component) == NULL)
				{
					std::stringstream errorMessage;
					errorMessage << "A LV Elastance Controller was requested for Component " << nodeOrComponentIndex;
					errorMessage << " of Netlist circuit number " << netlistCircuit->getIndexAmongstNetlists();
					errorMessage << " but this component is not a VolumeTrackingComponent." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				// Get the pointer to the compliance which needs to be controlled (in this case, the compliance of the pressure chamber):
				double* parameterToControl = component->getParameterPointer();

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new LeftVentricularElastanceController(parameterToControl,m_delt));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_BleedResistance:
			{
				// get the resistor:
				boost::shared_ptr<CircuitComponent> resistor = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				assert(resistor->getType() == Component_Resistor);

				double* resistanceToControl = resistor->getParameterPointer();

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new BleedController(resistanceToControl));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_BleedCompliance:
			{
				// get the capacitor:
				boost::shared_ptr<CircuitComponent> capacitor = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				assert(capacitor->getType() == Component_Capacitor);

				double* complianceToControl = capacitor->getParameterPointer();

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new BleedController(complianceToControl));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_CustomPythonComponentParameter:
			{
				// get the component:
				boost::shared_ptr<CircuitComponent> controlledComponent = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				double* parameterToControl = controlledComponent->getParameterPointer();
				std::string externalPythonControllerName;
				std::vector<std::pair<int,double*>> flowPointerPairs;
				std::vector<std::pair<int,double*>> pressurePointerPairs;
				std::vector<std::pair<int,double*>> volumePointerPairs;

				if (controlledComponent->hasUserDefinedExternalPythonScriptParameterController())
				{
					externalPythonControllerName = controlledComponent->getPythonControllerName();

					// Gather the pressures and flows as pointers, so the CustomPython parameter
					// controller can retrieve the pressure and flow values for each component
					// of its netlist, and pass them to the Python controller for use by the user.
					// They're indexed by the input data indices for the nodes / componnents:
					flowPointerPairs = netlistCircuit->getComponentInputDataIndicesAndFlows();
					pressurePointerPairs = netlistCircuit->getNodeInputDataIndicesAndPressures();
					volumePointerPairs = netlistCircuit->getVolumeTrackingComponentInputDataIndicesAndVolumes();

				}
				else
				{
					std::stringstream errorMessage;
					errorMessage << "EE: A component of the Netlist circuit at surface " << netlistCircuit->getSurfaceIndex() << 
					errorMessage << " was tagged as having an external Python parameter controller, but none was found." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new UserDefinedCustomPythonParameterController(parameterToControl, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_CustomPythonComponentFlow:
			{
				// Begin by getting the Python flow control script and copying it into the working directory, if
				// it doesn't exist there yet:
				char* crimsonFlowsolverHome;
				crimsonFlowsolverHome = getenv("CRIMSON_FLOWSOLVER_HOME");
				if (crimsonFlowsolverHome == NULL)
				{
					throw std::runtime_error("EE: Please set environmental variable CRIMSON_FLOWSOLVER_HOME to the root of the CRIMSON flowsolver source tree\n");
				}

				boost::filesystem::path crimsonFlowsolverHomePath(crimsonFlowsolverHome);
				if (!boost::filesystem::exists(crimsonFlowsolverHomePath))
				{
					throw std::runtime_error("EE: Error relating to environmental variable CRIMSON_FLOWSOLVER_HOME. Please check it is correctly set.\n");
				}

				// Construct a relative path with the location of the python flow control script we need:
				boost::filesystem::path pathOfPythonScriptRelativeToCrimsonFlowsolverHome("basicControlScripts/flowPrescriber.py");
				// Append to crimsonFlowsolverHomePath to get to the location of the python script we need:
				boost::filesystem::path pathToPythonScript = crimsonFlowsolverHomePath /= pathOfPythonScriptRelativeToCrimsonFlowsolverHome;

				boost::filesystem::path workingDirectory( boost::filesystem::current_path() );
				// we'll do the actual copy in a moment once we have aname for our destination file...


				// get the component:
				boost::shared_ptr<CircuitComponent> controlledComponent = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				std::string externalPythonControllerName;
				std::vector<std::pair<int,double*>> flowPointerPairs;
				std::vector<std::pair<int,double*>> pressurePointerPairs;
				std::vector<std::pair<int,double*>> volumePointerPairs;

				if (!controlledComponent->hasPrescribedFlow())
				{
					std::stringstream errorMessage;
					errorMessage << "EE: Component " << nodeOrComponentIndex << " of the netlist boundary condition at surface " << netlistCircuit->getSurfaceIndex() << std::endl;
					errorMessage << "was tagged has having a flow controller, but not as having prescribed flow. Please fix this." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				double* flowToControl = controlledComponent->getPointerToFixedFlowPrescription();

				if (controlledComponent->hasUserDefinedExternalPythonScriptParameterController())
				{
					externalPythonControllerName = controlledComponent->getPythonControllerName();

					// ... finish the copy that we began above, but didn't yet have the correct file name:
					std::string targetFileName_string = externalPythonControllerName;
					targetFileName_string.append(".py");

					boost::filesystem::path targetFileName_path = boost::filesystem::path(targetFileName_string.c_str());

					// Ensure only rank 0 does the copy:
					if (m_rank == 0)
					{
						boost::filesystem::copy_file(pathToPythonScript, workingDirectory /= targetFileName_path);
					}
					MPI_Barrier(MPI_COMM_WORLD);

					// Gather the pressures and flows as pointers, so the CustomPython parameter
					// controller can retrieve the pressure and flow values for each component
					// of its netlist, and pass them to the Python controller for use by the user.
					// They're indexed by the input data indices for the nodes / componnents:
					flowPointerPairs = netlistCircuit->getComponentInputDataIndicesAndFlows();
					pressurePointerPairs = netlistCircuit->getNodeInputDataIndicesAndPressures();
					volumePointerPairs = netlistCircuit->getVolumeTrackingComponentInputDataIndicesAndVolumes();

				}
				else
				{
					std::stringstream errorMessage;
					errorMessage << "EE: A component of the Netlist circuit at surface " << netlistCircuit->getSurfaceIndex() << 
					errorMessage << " was tagged as having an external Python parameter controller, but none was found." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new UserDefinedCustomPythonParameterController(flowToControl, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_CustomPythonNode:
			{
				// get the component:
				boost::shared_ptr<CircuitPressureNode> controlledNode = netlistCircuit->getNodeByInputDataIndex(nodeOrComponentIndex);
				double* pressureToControl = controlledNode->getPointerToFixedPressurePrescription();
				std::string externalPythonControllerName;
				std::vector<std::pair<int,double*>> flowPointerPairs;
				std::vector<std::pair<int,double*>> pressurePointerPairs;
				std::vector<std::pair<int,double*>> volumePointerPairs;

				if (controlledNode->hasUserDefinedExternalPythonScriptParameterController())
				{
					externalPythonControllerName = controlledNode->getPythonControllerName();

					// Gather the pressures and flows as pointers, so the CustomPython parameter
					// controller can retrieve the pressure and flow values for each component
					// of its netlist, and pass them to the Python controller for use by the user.
					// They're indexed by the input data indices for the nodes / componnents:
					flowPointerPairs = netlistCircuit->getComponentInputDataIndicesAndFlows();
					pressurePointerPairs = netlistCircuit->getNodeInputDataIndicesAndPressures();
					volumePointerPairs = netlistCircuit->getVolumeTrackingComponentInputDataIndicesAndVolumes();

				}
				else
				{
					std::stringstream errorMessage;
					errorMessage << "EE: A node of the Netlist circuit at surface " << netlistCircuit->getSurfaceIndex() << 
					errorMessage << " was tagged as having an external Python parameter controller, but none was found." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new UserDefinedCustomPythonParameterController(pressureToControl, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs));
				m_controlSystems.push_back(controllerToPushBack);
			}

			break;
		default:
			std::stringstream errorMessage;
			errorMessage << "EE: Unknown control parameter type (controllerType) requested for Netlist boundary condition for surface " << netlistCircuit->getSurfaceIndex() << "." << std::endl;
			throw std::runtime_error(errorMessage.str());
	}
}