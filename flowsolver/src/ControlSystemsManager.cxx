#include "ControlSystemsManager.hxx"
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <boost/filesystem.hpp>

void ControlSystemsManager::updateBoundaryConditionControlSystems()
{
	updateAndPassStateInformationBetweenPythonParameterControllers();

	// Update the non-python control systems
	for (auto controlSystem = m_nonPythonControlSystems.begin(); controlSystem != m_nonPythonControlSystems.end(); controlSystem++)
	{
		(*controlSystem)->updateControl();
	}
}

void ControlSystemsManager::updateAndPassStateInformationBetweenPythonParameterControllers()
{
	// Create a place to store the broadcast data as we gather it, and to use to send it to all controllers
	PyObject* gatheredBroadcastData = PyDict_New();
	assert(gatheredBroadcastData != NULL);

	// As we loop, the gatheredBroadcastData gets progressively larger, always containing all of the currently-avaialble (i.e. has-been-computed)
	// broadcast data. Control the order of this loop by using controllerPriority in the python control scripts.
	for (auto pythonControlSystem = m_pythonControlSystems.begin(); pythonControlSystem != m_pythonControlSystems.end(); pythonControlSystem++)
	{
		boost::shared_ptr<UserDefinedCustomPythonParameterController> downcastPythonControlSystem = boost::static_pointer_cast<UserDefinedCustomPythonParameterController> (*pythonControlSystem);
		// First we give the current controller all of the broadcast data that is available so far:
		downcastPythonControlSystem->giveStateDataFromOtherPythonControllers(gatheredBroadcastData);

		(*pythonControlSystem)->updateControl();

		PyObject* thisScriptsBroadcastData = NULL;
		downcastPythonControlSystem->getBroadcastStateData(thisScriptsBroadcastData);
		assert(thisScriptsBroadcastData != NULL);
		m_pythonBroadcastDataFromEachController.push_back(thisScriptsBroadcastData);

		int errCode = PyDict_Merge(gatheredBroadcastData, m_pythonBroadcastDataFromEachController.back(), 0); // Final "0" sets the merge not to over-write existing values with same key during dictionary merge
		assert(errCode == 0);
	}

	// Clean up:
	safe_Py_DECREF(gatheredBroadcastData);
	for (auto broadcastData = m_pythonBroadcastDataFromEachController.begin(); broadcastData != m_pythonBroadcastDataFromEachController.end(); broadcastData++)
	{
		safe_Py_DECREF(*broadcastData);
	}
	m_pythonBroadcastDataFromEachController.clear();

	writePythonControlSystemsRestarts();
}

void ControlSystemsManager::writePythonControlSystemsRestarts()
{
	bool theSimulationWasJustRestarted = (m_currentTimestepIndex == m_startingTimestepIndex);
	if (m_currentTimestepIndex % m_timestepsBetweenRestarts == 0 && !theSimulationWasJustRestarted)
	{
		for (auto pythonControlSystem = m_pythonControlSystems.begin(); pythonControlSystem != m_pythonControlSystems.end(); pythonControlSystem++)
		{
			(*pythonControlSystem)->picklePythonController();
		}
	}
	m_currentTimestepIndex++;
}

// void ControlSystemsManager::updateAndPassStateInformationBetweenPythonParameterControllers()
// {
// 	// Gather all the state data broadcasts from the boundary condition controllers
// 	for (auto pythonControlSystem = m_pythonControlSystems.begin(); pythonControlSystem != m_pythonControlSystems.end(); pythonControlSystem++)
// 	{
// 		boost::shared_ptr<UserDefinedCustomPythonParameterController> downcastPythonControlSystem = boost::static_pointer_cast<UserDefinedCustomPythonParameterController> (*pythonControlSystem);
// 		PyObject* thisScriptsBroadcastData = NULL;
// 		downcastPythonControlSystem->getBroadcastStateData(thisScriptsBroadcastData);
// 		assert(thisScriptsBroadcastData != NULL);
// 		m_pythonBroadcastDataFromEachController.push_back(thisScriptsBroadcastData);

// 		// safe_Py_DECREF(thisScriptsBroadcastData);
// 	}

// 	PyObject* gatheredBroadcastData = PyDict_New();
// 	assert(gatheredBroadcastData != NULL);
// 	// Package the broadcast data:
// 	for (auto broadcastData = m_pythonBroadcastDataFromEachController.begin(); broadcastData != m_pythonBroadcastDataFromEachController.end(); broadcastData++)
// 	{
// 		int errCode = PyDict_Merge(gatheredBroadcastData, *broadcastData, 0); // Final "0" sets the merge to not over-write existing values with same key during dictionary merge
// 		assert(errCode == 0);
// 	}

// 	// Update from the master control system script (not associated with any boundary)
// 	if(m_hasMasterPythonController)
// 	{
// 		mp_masterPythonController->giveStateDataFromOtherPythonControllers(gatheredBroadcastData);
// 		mp_masterPythonController->updateControl();

// 		PyObject* thisScriptsBroadcastData = NULL;
// 		mp_masterPythonController->getBroadcastStateData(thisScriptsBroadcastData);
// 		m_pythonBroadcastDataFromEachController.push_back(thisScriptsBroadcastData);
// 		// add the info from the master control script to the broadcast data:
// 		int errCode = PyDict_Merge(gatheredBroadcastData, thisScriptsBroadcastData, 0); // Final "0" sets the merge to not over-write existing values with same key during dictionary merge
// 		assert(errCode == 0);
		
// 		// safe_Py_DECREF(thisScriptsBroadcastData);
// 	}

// 	assert(gatheredBroadcastData != NULL);
// 	// Send the packaged broadcast data to all controllers
// 	for (auto pythonControlSystem = m_pythonControlSystems.begin(); pythonControlSystem != m_pythonControlSystems.end(); pythonControlSystem++)
// 	{
// 		boost::shared_ptr<UserDefinedCustomPythonParameterController> downcastPythonControlSystem = boost::static_pointer_cast<UserDefinedCustomPythonParameterController> (*pythonControlSystem);
// 		downcastPythonControlSystem->giveStateDataFromOtherPythonControllers(gatheredBroadcastData);
// 	}

// 	// Clean up:
// 	safe_Py_DECREF(gatheredBroadcastData);
// 	for (auto broadcastData = m_pythonBroadcastDataFromEachController.begin(); broadcastData != m_pythonBroadcastDataFromEachController.end(); broadcastData++)
// 	{
// 		safe_Py_DECREF(*broadcastData);
// 	}
// 	m_pythonBroadcastDataFromEachController.clear();
// }

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
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new LeftVentricularElastanceController(parameterToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Component, m_delt, m_startingTimestepIndex));
				m_nonPythonControlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_BleedResistance:
			{
				// get the resistor:
				boost::shared_ptr<CircuitComponent> resistor = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				assert(resistor->getType() == Component_Resistor);

				double* resistanceToControl = resistor->getParameterPointer();
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new BleedController(resistanceToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Component, m_startingTimestepIndex));
				m_nonPythonControlSystems.push_back(controllerToPushBack);
			}

			break;

		case Controller_BleedCompliance:
			{
				// get the capacitor:
				boost::shared_ptr<CircuitComponent> capacitor = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				assert(capacitor->getType() == Component_Capacitor);

				double* complianceToControl = capacitor->getParameterPointer();
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<AbstractParameterController> controllerToPushBack(new BleedController(complianceToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Component, m_startingTimestepIndex));
				m_nonPythonControlSystems.push_back(controllerToPushBack);
			}

			break;
		case Controller_CustomPythonComponentUnstressedVolume:
			{
				// get the component:
				boost::shared_ptr<CircuitComponent> controlledComponent = netlistCircuit->getComponentByInputDataIndex(nodeOrComponentIndex);
				boost::shared_ptr<VolumeTrackingPressureChamber> componentWithUnstressedVolumeToControl = boost::static_pointer_cast<VolumeTrackingPressureChamber> (controlledComponent);
				if (!componentWithUnstressedVolumeToControl) 
				{
					std::stringstream errorMessage;
					errorMessage << "EE: A component of the Netlist circuit at surface " << netlistCircuit->getSurfaceIndex() << 
					errorMessage << " was tagged as having unstressed volume control, but the component has no unstressed volume parameter." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}

				//////////////////////////////////////////////////////////////////////
				// this is the difference with normal Python component controllers: //
				//////////////////////////////////////////////////////////////////////
				double* parameterToControl = componentWithUnstressedVolumeToControl->getUnstressedVolumePointer();

				std::string externalPythonControllerName;
				std::vector<std::pair<int,double*>> flowPointerPairs;
				std::vector<std::pair<int,double*>> pressurePointerPairs;
				std::vector<std::pair<int,double*>> volumePointerPairs;

				if (componentWithUnstressedVolumeToControl->hasUserDefinedExternalPythonScriptParameterController())
				{
					externalPythonControllerName = componentWithUnstressedVolumeToControl->getPythonControllerName(controllerType);

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
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<GenericPythonController> controllerToPushBack(new UserDefinedCustomPythonParameterController(parameterToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Component, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs, m_startingTimestepIndex));
				m_pythonControlSystems.push_back(controllerToPushBack);
				controllerToPushBack->initialise();
				// m_controlSystems.push_back(controllerToPushBack);
				sortPythonControlSystemsByPriority();
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
					externalPythonControllerName = controlledComponent->getPythonControllerName(controllerType);

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
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<GenericPythonController> controllerToPushBack(new UserDefinedCustomPythonParameterController(parameterToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Component, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs, m_startingTimestepIndex));
				m_pythonControlSystems.push_back(controllerToPushBack);
				controllerToPushBack->initialise();
				// m_controlSystems.push_back(controllerToPushBack);
				sortPythonControlSystemsByPriority();
			}

			break;

		case Controller_CustomPythonComponentFlowFile:
			{
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
					externalPythonControllerName = controlledComponent->getPythonControllerName(controllerType);

					// Begin by getting the Python flow control script and copying it into the working directory, if
					// it doesn't exist there yet:
					setupPythonBoilerplateScriptPaths();

					std::string targetFileName_string = externalPythonControllerName;
					targetFileName_string.append(".py");
					copyFileToWorkingDirectory(m_pathToBoilerplatePythonFlowPrescriberScript, targetFileName_string);

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
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<GenericPythonController> controllerToPushBack(new UserDefinedCustomPythonParameterController(flowToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Component, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs, m_startingTimestepIndex));
				m_pythonControlSystems.push_back(controllerToPushBack);
				controllerToPushBack->initialise();
				// m_controlSystems.push_back(controllerToPushBack);
				sortPythonControlSystemsByPriority();
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
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<GenericPythonController> controllerToPushBack(new UserDefinedCustomPythonParameterController(pressureToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Node, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs, m_startingTimestepIndex));
				// m_controlSystems.push_back(controllerToPushBack);
				m_pythonControlSystems.push_back(controllerToPushBack);
				controllerToPushBack->initialise();
				sortPythonControlSystemsByPriority();
			}

			break;
		case Controller_CustomPythonNodePressureFile:
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

					// Begin by getting the Python flow control script and copying it into the working directory, if
					// it doesn't exist there yet:
					setupPythonBoilerplateScriptPaths();
					std::string targetFileName_string = externalPythonControllerName;
					targetFileName_string.append(".py");
					copyFileToWorkingDirectory(m_pathToBoilerplatePythonPressurePrescriberScript, targetFileName_string);

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
					errorMessage << "EE: A node of the Netlist circuit at surface " << netlistCircuit->getSurfaceIndex();
					errorMessage << " was tagged as having an external Python parameter controller, but none was found." << std::endl;
					throw std::runtime_error(errorMessage.str());
				}
				int surfaceIndex = netlistCircuit->getSurfaceIndex();
				boost::shared_ptr<GenericPythonController> controllerToPushBack(new UserDefinedCustomPythonParameterController(pressureToControl, surfaceIndex, nodeOrComponentIndex, circuit_item_t::Circuit_Node, m_delt, externalPythonControllerName, flowPointerPairs, pressurePointerPairs, volumePointerPairs, m_startingTimestepIndex));
				m_pythonControlSystems.push_back(controllerToPushBack);
				controllerToPushBack->initialise();
				// m_controlSystems.push_back(controllerToPushBack);
				sortPythonControlSystemsByPriority();
			}
			
			break;

		default:
			std::stringstream errorMessage;
			errorMessage << "EE: Unknown control parameter type (controllerType) requested for Netlist boundary condition for surface " << netlistCircuit->getSurfaceIndex() << "." << std::endl;
			throw std::runtime_error(errorMessage.str());
	}
}

void ControlSystemsManager::sortPythonControlSystemsByPriority()
{
	std::sort(m_pythonControlSystems.begin(), m_pythonControlSystems.end(), 
		              [](boost::shared_ptr<GenericPythonController> controlSystem1ToCompare, boost::shared_ptr<GenericPythonController> controlSystem2ToCompare)
		                 {return (controlSystem1ToCompare->getPriority() < controlSystem2ToCompare->getPriority());});
}

void ControlSystemsManager::createMasterPythonController()
{
	// boost::shared_ptr<CircuitPressureNode> controlledNode = netlistCircuit->getNodeByInputDataIndex(nodeOrComponentIndex);
	// double* pressureToControl = controlledNode->getPointerToFixedPressurePrescription();
	std::string externalPythonControllerName("masterController");
	// std::vector<std::pair<int,double*>> flowPointerPairs;
	// std::vector<std::pair<int,double*>> pressurePointerPairs;
	// std::vector<std::pair<int,double*>> volumePointerPairs;

	// if (controlledNode->hasUserDefinedExternalPythonScriptParameterController())
	// {
		// externalPythonControllerName = controlledNode->getPythonControllerName();

		// Gather the pressures and flows as pointers, so the CustomPython parameter
		// controller can retrieve the pressure and flow values for each component
		// of its netlist, and pass them to the Python controller for use by the user.
		// They're indexed by the input data indices for the nodes / componnents:
		// flowPointerPairs = netlistCircuit->getComponentInputDataIndicesAndFlows();
		// pressurePointerPairs = netlistCircuit->getNodeInputDataIndicesAndPressures();
		// volumePointerPairs = netlistCircuit->getVolumeTrackingComponentInputDataIndicesAndVolumes();

	// }
	// else
	// {
	// 	std::stringstream errorMessage;
	// 	errorMessage << "EE: A node of the Netlist circuit at surface " << netlistCircuit->getSurfaceIndex() << 
	// 	errorMessage << " was tagged as having an external Python parameter controller, but none was found." << std::endl;
	// 	throw std::runtime_error(errorMessage.str());
	// }
	// int surfaceIndex = netlistCircuit->getSurfaceIndex();
	boost::shared_ptr<GenericPythonController> newMasterController(new GenericPythonController(m_delt, externalPythonControllerName, m_startingTimestepIndex));
	mp_masterPythonController = newMasterController;
	mp_masterPythonController->initialise();
	m_pythonControlSystems.push_back(mp_masterPythonController);
	sortPythonControlSystemsByPriority();
}

void ControlSystemsManager::setupPythonBoilerplateScriptPaths()
{
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
	boost::filesystem::path pathOfPythonFlowPrescriberScriptRelativeToCrimsonFlowsolverHome("basicControlScripts/flowPrescriber.py");
	// Append to crimsonFlowsolverHomePath to get to the location of the python script we need:
	m_pathToBoilerplatePythonFlowPrescriberScript = crimsonFlowsolverHomePath;
	m_pathToBoilerplatePythonFlowPrescriberScript /= pathOfPythonFlowPrescriberScriptRelativeToCrimsonFlowsolverHome;

	boost::filesystem::path pathOfPythonPressurePrescriberScriptRelativeToCrimsonFlowsolverHome("basicControlScripts/pressurePrescriber.py");
	m_pathToBoilerplatePythonPressurePrescriberScript = crimsonFlowsolverHomePath;
	m_pathToBoilerplatePythonPressurePrescriberScript /= pathOfPythonPressurePrescriberScriptRelativeToCrimsonFlowsolverHome;
}

void ControlSystemsManager::copyFileToWorkingDirectory(const boost::filesystem::path sourcePath, const std::string targetFileName) const
{
	boost::filesystem::path targetFileName_path = boost::filesystem::path(targetFileName.c_str());

	// Ensure only rank 0 does the copy:
	if (m_rank == 0)
	{
		boost::filesystem::path copyOfWorkingDir = m_workingDirectory;
		try
		{
			boost::filesystem::copy_file(sourcePath, copyOfWorkingDir /= targetFileName_path);
		}
		catch (boost::filesystem::filesystem_error error)
		{
			if (error.code() == boost::system::errc::errc_t::file_exists)
			{
				std::cerr << "WW: Not copying file from " << sourcePath << " as the target " << targetFileName_path << " exists. If this is unwanted, delete the target." << std::endl;
			}
			else
			{
				throw;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

int ControlSystemsManager::getNumberOfControlSystems() const
{
	return m_pythonControlSystems.size() + m_nonPythonControlSystems.size();
}