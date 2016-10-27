#include "NetlistCircuit.hxx"
#include "fileReaders.hxx"
#include "fileWriters.hxx"
#include <boost/make_shared.hpp>
#include "indexShifters.hxx"
#include "SimvascularGlobalArrayTransfer.h"

void NetlistCircuit::initialisePetscArrayNames()
{
    m_RHS = PETSC_NULL;
    m_solutionVector = PETSC_NULL;
    m_systemMatrix = PETSC_NULL;
    m_inverseOfSystemMatrix = PETSC_NULL;
    m_identityMatrixForPetscInversionHack = PETSC_NULL;
}

int NetlistCircuit::getNumberOfDegreesOfFreedom() const
{
    int degreesOfFreedom = mp_circuitData->numberOfPressureNodes;
    degreesOfFreedom += m_numberOfHistoryPressures;
    degreesOfFreedom += mp_circuitData->numberOfComponents;
    degreesOfFreedom += numberOfHistoryFlows;
    degreesOfFreedom += m_numberOfTrackedVolumes;
    degreesOfFreedom += numberOfHistoryVolumes;

    return degreesOfFreedom;
}

bool NetlistCircuit::surfaceIndexMatches(const int surfaceIndexToTest) const
{
    if (m_surfaceIndex == surfaceIndexToTest)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void NetlistCircuit::terminatePetscArrays()
{
    PetscErrorCode errFlag;
    if (m_RHS)
    {
        errFlag = VecDestroy(&m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_solutionVector)
    {
        errFlag = VecDestroy(&m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_systemMatrix)
    {
        errFlag = MatDestroy(&m_systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_inverseOfSystemMatrix)
    {
        errFlag = MatDestroy(&m_inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    if (m_identityMatrixForPetscInversionHack)
    {
        errFlag = MatDestroy(&m_identityMatrixForPetscInversionHack); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
}

void NetlistCircuit::createCircuitDescription()
{
    // This function takes the read-in netlist circuit description and converts it
    // to the internal CircuitData class format.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    mp_netlistFileReader = NetlistReader::Instance();
    mp_netlistXmlReader = NetlistXmlReader::Instance();
    createBasicCircuitDescription();

    // Tell the node at the 3D interface that it connects to the 3D domain:
    {
        int threeDNodeIndex = mp_netlistXmlReader->getIndicesOfNodesAt3DInterface().at(m_IndexOfThisNetlistLPNInInputFile);
        mp_circuitData->initialiseNodeAndComponentAtInterface(threeDNodeIndex);
    }

}

void NetlistCircuit::setupCustomPythonControlSystems()
{
    // Give any components tagged in the input data for custom Python
    // control systems the name of their Python controller script:
    std::map<int, ComponentControlSpecificationContainer> userDefinedComponentControllers = mp_netlistXmlReader->getUserDefinedComponentControllersAndPythonNames(m_IndexOfThisNetlistLPNInInputFile);
    for (auto componentIndexAndPythonControllerSpecification = userDefinedComponentControllers.begin(); componentIndexAndPythonControllerSpecification != userDefinedComponentControllers.end(); componentIndexAndPythonControllerSpecification++)
    {
        // There can be more than one controller attached to this component (e.g. to control both compliance and unstressed volume)
        ComponentControlSpecificationContainer& controlInfo = componentIndexAndPythonControllerSpecification->second;
        for (int controllerIndexForThisComponent = 0; controllerIndexForThisComponent < controlInfo.getNumberOfControlScripts(); controllerIndexForThisComponent++)
        {
            std::string pythonScriptName = controlInfo.getControlScriptNameByIndexLocalToComponent(controllerIndexForThisComponent);
            parameter_controller_t controllerType = controlInfo.getControlTypeByIndexLocalToComponent(controllerIndexForThisComponent);
            mp_circuitData->mapOfComponents.at(componentIndexAndPythonControllerSpecification->first)->addPythonControllerName(controllerType, pythonScriptName);
        }
    }
    // Give any nodes tagged in the input data for custom Python
    // control systems the name of their Python controller script:
    std::map<int,std::string> userDefinedNodeControllers = mp_netlistXmlReader->getUserDefinedNodeControllersAndPythonNames(m_IndexOfThisNetlistLPNInInputFile);
    for (auto nodeIndexAndPythonName = userDefinedNodeControllers.begin(); nodeIndexAndPythonName != userDefinedNodeControllers.end(); nodeIndexAndPythonName++)
    {
        // Nodal parameter (i.e. presure) controllers can only be attached
        // to _prescribed_ pressure nodes - here we try to get the requested
        // node from mapOfPrescribedPressureNodes; if this fails, we throw an error telling
        // the user that they can only use prescribed pressure nodes here.
        try
        {
            mp_circuitData->mapOfPrescribedPressureNodes.at(nodeIndexAndPythonName->first)->setPythonControllerName(nodeIndexAndPythonName->second);
        }
        catch (std::out_of_range& outOfRange)
        {
            bool pressureNodeExistsButIsNotAPrescribedPressureNode = (mp_circuitData->mapOfPressureNodes.count(nodeIndexAndPythonName->first) == 1);
            if (pressureNodeExistsButIsNotAPrescribedPressureNode)
            {
                std::stringstream errorMessage;
                errorMessage << "EE: Nodal parameter controllers can only be attached to prescribed pressure nodes, as specified in the Netlist description." << std::endl;
                throw std::runtime_error(errorMessage.str());
            }
            else
            {
                // If the node doesn't even exist, just rethrow the original out_of_range error
                throw outOfRange;
            }
        }
        
    }
}

void NetlistCircuit::createBasicCircuitDescription()
{
    mp_circuitData->numberOfComponents = mp_netlistXmlReader->getNumberOfComponents().at(m_IndexOfThisNetlistLPNInInputFile);
    mp_circuitData->numberOfPressureNodes = mp_netlistXmlReader->getNumberOfPressureNodes().at(m_IndexOfThisNetlistLPNInInputFile);
    mp_circuitData->numberOfPrescribedPressures = mp_netlistXmlReader->getNumberOfPrescribedPressures().at(m_IndexOfThisNetlistLPNInInputFile);
    mp_circuitData->numberOfPrescribedFlows = mp_netlistXmlReader->getNumberOfPrescribedFlows().at(m_IndexOfThisNetlistLPNInInputFile);

    std::vector<circuit_component_t> retrievedComponentTypes = mp_netlistXmlReader->getComponentTypes().at(m_IndexOfThisNetlistLPNInInputFile);

    // Prepare space for the components in the circuit:
    assert(mp_circuitData->components.empty());
    for (int componentIndex=0; componentIndex<mp_circuitData->numberOfComponents; componentIndex++)
    {
        CircuitComponent* toPushBack;
        if (retrievedComponentTypes.at(componentIndex) == Component_VolumeTrackingPressureChamber)
        {
            double initialVolume = mp_netlistXmlReader->getComponentInitialVolume(m_IndexOfThisNetlistLPNInInputFile, componentIndex);
            double initialUnstressedVolume = mp_netlistXmlReader->getComponentInitialUnstressedVolume(m_IndexOfThisNetlistLPNInInputFile, componentIndex);
            toPushBack = new VolumeTrackingPressureChamber(m_hstep,m_thisIsARestartedSimulation, initialVolume, initialUnstressedVolume);
        }
        else if (retrievedComponentTypes.at(componentIndex) == Component_VolumeTracking)
        {
            double initialVolume = mp_netlistXmlReader->getComponentInitialVolume(m_IndexOfThisNetlistLPNInInputFile, componentIndex);
            toPushBack = new VolumeTrackingComponent(m_hstep,m_thisIsARestartedSimulation, initialVolume);
        }
        else
        {
            toPushBack = new CircuitComponent(m_hstep,m_thisIsARestartedSimulation);
        }

        mp_circuitData->components.push_back(boost::shared_ptr<CircuitComponent> (toPushBack));
        mp_circuitData->components.back()->setIndex(toOneIndexing(componentIndex));
    }


    // Obtain the component- and node-level data for the circuit, for moving into the appropriate data structure CircuitData
    // We want to pop off the component types as we use them, but starting from the beginning of the vector. To do this, we reverse
    // the vector and then pop from the new end.
    std::reverse(retrievedComponentTypes.begin(), retrievedComponentTypes.end());
    std::vector<int> retrievedComponentStartNodes = mp_netlistXmlReader->getComponentStartNodes().at(m_IndexOfThisNetlistLPNInInputFile);
    std::reverse(retrievedComponentStartNodes.begin(), retrievedComponentStartNodes.end());
    std::vector<int> retrievedComponentEndNodes = mp_netlistXmlReader->getComponentEndNodes().at(m_IndexOfThisNetlistLPNInInputFile);
    std::reverse(retrievedComponentEndNodes.begin(), retrievedComponentEndNodes.end());
    std::vector<double> retrievedComponentParameterValues = mp_netlistXmlReader->getComponentParameterValues(m_IndexOfThisNetlistLPNInInputFile);
    std::reverse(retrievedComponentParameterValues.begin(), retrievedComponentParameterValues.end());

    std::vector<int> retrievedListOfPrescribedFlows = mp_netlistXmlReader->getListOfPrescribedFlows().at(m_IndexOfThisNetlistLPNInInputFile);
    std::vector<circuit_component_flow_prescription_t> retrievedTypeOfPrescribedFlows = mp_netlistXmlReader->getTypeOfPrescribedFlows().at(m_IndexOfThisNetlistLPNInInputFile);
    std::vector<double> retrievedValueOfPrescribedFlows = mp_netlistXmlReader->getValueOfPrescribedFlows().at(m_IndexOfThisNetlistLPNInInputFile);
    std::map<int,double> retrievedInitialPressures = mp_netlistXmlReader->getInitialPressures().at(m_IndexOfThisNetlistLPNInInputFile);

    std::set<int> retrievedKalmanFilteredComponentIndices = mp_netlistXmlReader->getKalmanFilteredComponentIndicesByCircuitIndex(m_IndexOfThisNetlistLPNInInputFile);
    
    // Loop over the components, assigning them (and their nodes) the appropriate properties to give the fully-described circuit:
    for (auto component = mp_circuitData->components.begin(); component != mp_circuitData->components.end(); component++)
    {
        (*component)->getType() = retrievedComponentTypes.back();
        retrievedComponentTypes.pop_back();

        int indexOfStartNodeInInputData = retrievedComponentStartNodes.back();
        retrievedComponentStartNodes.pop_back();
        
        // note that we're passing a boost::shared_ptr to the startNode here; if this is currently NULL, it will
        // be constructed during the call to setupPressureNode.
        setupPressureNode(indexOfStartNodeInInputData, (*component)->startNode, *component);

        // (*component)->startNode = mp_circuitData.ifExistsGetNodeOtherwiseConstructNode(indexOfStartNodeInInputData);
        // (*component)->startNode->getIndex() = indexOfStartNodeInInputData;

        // (*component)->startNode->prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        // for (int prescribedPressure=0; prescribedPressure<mp_circuitData.numberOfPrescribedPressures; prescribedPressure++)
        // {
        //     if (retrievedListOfPrescribedPressures.at(prescribedPressure) == (*component)->startNode->getIndex())
        //     {
        //         (*component)->startNode->prescribedPressureType = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->startNode->pressure = retrievedValueOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->startNode->pressure = (*component)->startNode->pressure;
        //     }
        // }

        int indexOfEndNodeInInputData = retrievedComponentEndNodes.back();
        retrievedComponentEndNodes.pop_back();

        // note that we're passing a boost::shared_ptr to the endNode here; if this is currently NULL, it will
        // be constructed during the call to setupPressureNode.
        setupPressureNode(indexOfEndNodeInInputData, (*component)->endNode, *component);

        // int indexOfEndNodeInInputData = retrievedComponentEndNodes.back();
        // retrievedComponentEndNodes.pop_back();
        // (*component)->endNode = mp_circuitData.ifExistsGetNodeOtherwiseConstructNode(indexOfEndNodeInInputData);
        // (*component)->endNode->getIndex() = indexOfEndNodeInInputData;
        // (*component)->endNode->prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        // for (int prescribedPressure=0; prescribedPressure<mp_circuitData.numberOfPrescribedPressures; prescribedPressure++)
        // {
        //     if (retrievedListOfPrescribedPressures.at(prescribedPressure) == (*component)->endNode->getIndex())
        //     {
        //         (*component)->endNode->prescribedPressureType = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->endNode->pressure = retrievedValueOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->endNode->pressure = (*component)->endNode->pressure;
        //     }
        // }

        (*component)->prescribedFlowType = Flow_NotPrescribed;  // initialise as a default, before replacing as necessary
        for (int prescribedFlow=0; prescribedFlow<mp_circuitData->numberOfPrescribedFlows; prescribedFlow++)
        {
            if (retrievedListOfPrescribedFlows.at(prescribedFlow) == (*component)->getIndex())
            {
                (*component)->prescribedFlowType = retrievedTypeOfPrescribedFlows.at(prescribedFlow);
                (*component)->setPrescribedFlow(retrievedValueOfPrescribedFlows.at(prescribedFlow));
            }
        }

        (*component)->setParameterValue(retrievedComponentParameterValues.back());
        // make a copy of this value so we can reset it if necessary. Used for e.g. diode state changes.
        (*component)->parameterValueFromInputData = *((*component)->getParameterPointer());
        retrievedComponentParameterValues.pop_back();

        (*component)->startNode->setPressure(retrievedInitialPressures.at((*component)->startNode->getIndex()));
        (*component)->endNode->setPressure(retrievedInitialPressures.at((*component)->endNode->getIndex()));


        // Kalman filtering setup
        {
            SimvascularGlobalArrayTransfer* gat = SimvascularGlobalArrayTransfer::Get();
            // see if the component's index is on the list of Kalman-filtered components:
            if (retrievedKalmanFilteredComponentIndices.count((*component)->getIndex()) == 1)
            {
                std::stringstream componentNametagBuilder;
                componentNametagBuilder << "circuit_" << m_surfaceIndex << "_component_" << (*component)->getIndex();

                gat->setPointerToFilteredNetlistParameter((*component)->getParameterPointer(), componentNametagBuilder.str());
            }
        }


    }

    // Some metadata is already set-up for this circuit, as a side-effect
    // of the above construction. This call completes the metadata; it's not
    // a problem that it also re-writes some of the existing metadata
    // (rewrites - but does not change - the values are identical!)
    mp_circuitData->rebuildCircuitMetadata();
    mp_circuitData->setupComponentNeighbourPointers();

    // mp_circuitData.switchDiodeStatesIfNecessary();
    // mp_circuitData.detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    // // Component indices are just consecutive integers by default, but sometimes non-consecutive numbering
    // // is needed; componentIndices allows for this.
    // // We initialise it now for the default case.
    // for (int ii=1; ii < mp_circuitData.numberOfComponents + 1; ii++)
    // {
    //     mp_circuitData.componentIndices.push_back(ii);
    // }
    setupCustomPythonControlSystems();
}

// This function gets a pointer to, or creates anew, a pressure node (dependent on whether the node already has been
// constructed, for example because a node is shared by more than one component, and another component has already
// set up the node). It includes the code which decides whether to set up a volume tracking or a simple pressure node.
//
// It also adds the passed-in componentNeighbouringThisNode to the list of components attached to this node, whether or
// not the node is constructed on this call.
//
// Note that the pressure node is not fully set-up until this function has been called once for each component that
// neighbours the node.
void NetlistCircuit::setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode)
{
    // Access the read-in file data:
    std::vector<int> retrievedListOfPrescribedPressures = mp_netlistXmlReader->getListOfPrescribedPressures().at(m_IndexOfThisNetlistLPNInInputFile);
    std::vector<circuit_nodal_pressure_prescription_t> retrievedTypeOfPrescribedPressures = mp_netlistXmlReader->getTypeOfPrescribedPressures().at(m_IndexOfThisNetlistLPNInInputFile);
    std::vector<double> retrievedValueOfPrescribedPressures = mp_netlistXmlReader->getValueOfPrescribedPressures().at(m_IndexOfThisNetlistLPNInInputFile);

    // Discover whether this node has a prescribed pressure, and if so, what type:
    circuit_nodal_pressure_prescription_t typeOfPrescribedPressure = Pressure_NotPrescribed; // initialise, but chnage later if pressure is actually prescribed
    int indexOfPrescribedPressure = -1; // initialise to a nonsense value to detect errors
    for (int prescribedPressure=0; prescribedPressure<mp_circuitData->numberOfPrescribedPressures; prescribedPressure++)
    {
        if (retrievedListOfPrescribedPressures.at(prescribedPressure) == indexOfNodeInInputData)
        {
            typeOfPrescribedPressure = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
            indexOfPrescribedPressure = prescribedPressure;
        }
    }
    // Get the node (or create a new node if this one hasn't been made yet)
    node = mp_circuitData->ifExistsGetNodeOtherwiseConstructNode(indexOfNodeInInputData,typeOfPrescribedPressure,componentNeighbouringThisNode);    
    
    if (node->getPressurePrescriptionType() != Pressure_NotPrescribed && node->getPressurePrescriptionType() != Pressure_Null)
    {
        node->setPrescribedPressure(retrievedValueOfPrescribedPressures.at(indexOfPrescribedPressure));
    }
}

bool NetlistCircuit::flowPermittedAcross3DInterface() const
{
	return mp_circuitData->flowPermittedAcross3DInterface();
}

bool NetlistCircuit::boundaryConditionTypeHasJustChanged()
{
	return mp_circuitData->boundaryConditionTypeHasJustChanged();
}

void NetlistCircuit::closeAllDiodes()
{
	mp_circuitData->closeAllDiodes();
}

void NetlistCircuit::detectWhetherClosedDiodesStopAllFlowAt3DInterface()
{
	mp_circuitData->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
}

void NetlistCircuit::switchDiodeStatesIfNecessary()
{
	mp_circuitData->switchDiodeStatesIfNecessary();
}

void NetlistCircuit::rebuildCircuitMetadata()
{
	mp_circuitData->rebuildCircuitMetadata();
}

void NetlistCircuit::loadPressuresFlowsAndVolumesOnRestart()
{
    if (m_startingTimestepIndex > 0)
    {
        {
          // Write the netlistPressures_zeroDDomainReplacement.dat
          histFileReader boundaryConditionPressureHistoryReader;
          boundaryConditionPressureHistoryReader.setFileName(m_PressureHistoryFileName);
          boundaryConditionPressureHistoryReader.setNumColumns(mp_circuitData->numberOfPressureNodes + 1); // +1 for the timestep indexing column
          boundaryConditionPressureHistoryReader.readAndSplitMultiSurfaceRestartFile();

          for (int stepToRead=0; stepToRead <= m_startingTimestepIndex; stepToRead++)
          {
            boundaryConditionPressureHistoryReader.getNextDatum(); // ditch the timestep index
            for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
            {
              node->second->appendToPressureHistory(boundaryConditionPressureHistoryReader.getNextDatum());
            }
          }

          // Set the restarted nodal pressures:
          for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
            {
              node->second->setRestartPressureFromHistory();
            }
        }

        {
          // Write the netlistFlows_surface_X.dat
          histFileReader boundaryConditionFlowistoryReader;
          boundaryConditionFlowistoryReader.setFileName(m_FlowHistoryFileName);
          boundaryConditionFlowistoryReader.setNumColumns(mp_circuitData->numberOfComponents + 1); // +1 for the timestep indexing column
          boundaryConditionFlowistoryReader.readAndSplitMultiSurfaceRestartFile();

          for (int stepToRead=0; stepToRead <= m_startingTimestepIndex; stepToRead++)
          {
            boundaryConditionFlowistoryReader.getNextDatum(); // ditch the timestep index
            for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
            {
              (*component)->appendToFlowHistory(boundaryConditionFlowistoryReader.getNextDatum());
            }
          }

          // Set the restarted component flow:
          for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
          {
            (*component)->setRestartFlowFromHistory();
          }
        }


        if (mp_circuitData->m_numberOfVolumeTrackingComponenets > 0)
        {
          // Write the volumes of the volume tracking components, as netlistVolumes_surface_X.dat
          histFileReader boundaryConditionVolumeHistoryReader;
          boundaryConditionVolumeHistoryReader.setFileName(m_VolumeHistoryFileName);
          boundaryConditionVolumeHistoryReader.setNumColumns(mp_circuitData->m_numberOfVolumeTrackingComponenets + 1); // +1 for the timestep indexing column
          boundaryConditionVolumeHistoryReader.readAndSplitMultiSurfaceRestartFile();

          for (int stepToRead=0; stepToRead <= m_startingTimestepIndex; stepToRead++)
          {
            boundaryConditionVolumeHistoryReader.getNextDatum(); // ditch the timestep index
            for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
            {
              VolumeTrackingComponent* volumeTrackingComponent = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
              // If this component is actually a volume chamber, so it actually has a volume history we can write to the file:
              if (volumeTrackingComponent != NULL)
              {
                double readVolumeValue = boundaryConditionVolumeHistoryReader.getNextDatum();
                volumeTrackingComponent->setVolumeHistoryAtTimestep(readVolumeValue);
              }
            }
          }

          // Set the restarted stored volumes:
          for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
          {
            VolumeTrackingComponent* volumeTrackingComponent = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
            // If this component is actually a volume chamber, so it actually has a volume history we can write to the file:
            if (volumeTrackingComponent != NULL)
            {
              volumeTrackingComponent->setRestartVolumeFromHistory();
            }
          }
        }

        switchDiodeStatesIfNecessary();
    }
}

// nextTimestepWrite_start will be updated and returned to caller of this function.
void NetlistCircuit::writePressuresFlowsAndVolumes(int& nextTimestepWrite_start)
{
	// All the following writes can use the same nextTimestepWrite_end to determine how far to go when looking for data to write to the file:
      int nextTimestepWrite_end = mp_circuitData->getLengthOfHistoryData();

      try {
            {
              // Write the netlistPressures_zeroDDomainReplacement.dat
              basicFileWriter boundaryConditionPressureHistoryWriter;
              boundaryConditionPressureHistoryWriter.setFileName(m_PressureHistoryFileName);
      
              for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
              {
                boundaryConditionPressureHistoryWriter.writeStepIndex(stepToWrite);
                for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
                {
                  boundaryConditionPressureHistoryWriter.writeToFile(node->second->getFromPressureHistoryByTimestepIndex(stepToWrite));
                }
                boundaryConditionPressureHistoryWriter.writeEndLine();
              }
            }
      
            {
              // Write the netlistFlows_surface_X.dat
              basicFileWriter boundaryConditionFlowHistoryWriter;
              boundaryConditionFlowHistoryWriter.setFileName(m_FlowHistoryFileName);
      
              for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
              {
                boundaryConditionFlowHistoryWriter.writeStepIndex(stepToWrite);
                for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
                {
                  boundaryConditionFlowHistoryWriter.writeToFile((*component)->getFromFlowHistoryByTimestepIndex(stepToWrite));
                }
                boundaryConditionFlowHistoryWriter.writeEndLine();
              }
            }
      
      
            {
              // Write the volumes of the volume tracking components, as netlistVolumes_surface_X.dat
              basicFileWriter boundaryConditionVolumeHistoryWriter;
              boundaryConditionVolumeHistoryWriter.setFileName(m_VolumeHistoryFileName);
              for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
              {
                boundaryConditionVolumeHistoryWriter.writeStepIndex(stepToWrite);
                for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
                {
                  VolumeTrackingComponent* volumeTrackingComponent = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
                  // If this component is actually a volume chamber, so it actually has a volume history we can write to the file:
                  if (volumeTrackingComponent != NULL)
                  {
                    boundaryConditionVolumeHistoryWriter.writeToFile(volumeTrackingComponent->getVolumeHistoryAtTimestep(stepToWrite));
                  }
                }
                boundaryConditionVolumeHistoryWriter.writeEndLine();
              }
            }
      } catch (const std::exception& e) {
          std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
          throw;
      }

      // Set the starting place for the write on the next iteration (nextTimestepWrite_start is returned to the caller, by reference).
      nextTimestepWrite_start = nextTimestepWrite_end;
}

void NetlistCircuit::setPressureAndFlowPointers(double* pressurePointer, double* flowPointer)
{
    flow_n_ptrs.clear();
    flow_n_ptrs.push_back(flowPointer);

    pressure_n_ptrs.clear();
    pressure_n_ptrs.push_back(pressurePointer);
}

// void NetlistCircuit::identifyAtomicSubcircuits()
// {
//     createInitialCircuitDescriptionWithoutDiodes();

//     // The atomic subcircuits are those which cannot be broken down
//     // into further subcircuits by diodes. More complex ones (where some 
//     // valves are open, so some atomic subcircuits are joined) will 
//     // appear later.
//     assignComponentsToAtomicSubcircuits();
// }

// void NetlistCircuit::createInitialCircuitDescriptionWithoutDiodes()
// {
//     // Copy the data as a base for modification (i.e. the removal of the diode data)
//     mp_circuitDataWithoutDiodes = boost::make_shared<CircuitData> (*mp_circuitData);
//     // Prepare for diode removal by clearing the data that will be rebuilt:
//     // mp_circuitDataWithoutDiodes->deleteSubcircuitSpecificData();
//     // Rebuild the cleared vectors, but without the diodes:
//     int numberOfDiodes = 0;
//     for (auto component=mp_circuitDataWithoutDiodes->components.begin(); component!=mp_circuitDataWithoutDiodes->components.end(); component++)
//     {
//         if((*component)->getType() == Component_Diode)
//         {
//             component = mp_circuitDataWithoutDiodes->components.erase(component); // returns an iterator pointing to the new location of the element after the one that just got erased
//             component--; // decrement the returned iterator, so the for loop incrememnts it back to the correct (next) element.
//         }
//     }

//     // A tidier (and untested!) version of the above, without the ugly "component --" line:
//     // auto component=mp_circuitDataWithoutDiodes.components.begin();
//     // while (component!=mp_circuitDataWithoutDiodes.components.end())
//     // {
//     //     // both branches of this if guard result in an iterator to the next component after the one just checked (whether or not an erasing occurs)
//     //     bool componentNotADiode = !((*component)->type == Component_Diode);
//     //     if(componentNotADiode)
//     //     {
//     //         component++
//     //     }
//     //     else
//     //     {
//     //         // erase the diode:
//     //         component = mp_circuitDataWithoutDiodes.components.erase(component); // returns an iterator pointing to the new location of the element after the one that just got erased            
//     //     }
//     // }

//     mp_circuitDataWithoutDiodes->rebuildCircuitMetadata();

// }

// void NetlistCircuit::assignComponentsToAtomicSubcircuits()
// {
//     // This subroutine builds m_AtomicSubcircuitsComponentsBelongsTo, which is indexed by component, as they appear in mp_circuitDataWithoutDiodes.

//     // Group the components by partitioning the node indices into connected subcircuits (with diodes removed).
//     // We do this by parsing the circuit to discover the topology.
//     std::vector<bool> componentAssignedToASubcircuit(mp_circuitDataWithoutDiodes->numberOfComponents, false);
//     m_AtomicSubcircuitsComponentsBelongsTo.insert(m_AtomicSubcircuitsComponentsBelongsTo.begin(), mp_circuitDataWithoutDiodes->numberOfComponents, -1); // initialise with nonsense value
//     int subcircuitIndex = 0;
//     for (int startingComponent = 0; startingComponent < mp_circuitDataWithoutDiodes->numberOfComponents; startingComponent++)
//     {
//         // Only assign the startingComponent to a subcircuit if it doesn't yet belong to a subcircuit
//         if (!componentAssignedToASubcircuit.at(startingComponent))
//         {
//             // Because we immediately proceed to find the whole subcircuit, this is definitely a new
//             // subcircuit we're starting now.
//             m_AtomicSubcircuitsComponentsBelongsTo.at(startingComponent) = subcircuitIndex;
//             componentAssignedToASubcircuit.at(startingComponent) = true;

//             // search the rest of the list for the remaining components belonging to this subcircuit:
//             // Make a list to be populated with nodes belonging to the subcircuit, as they're found.
//             std::set<int> nodesInSubcircuit;
//             nodesInSubcircuit.insert(mp_circuitDataWithoutDiodes->components.at(startingComponent)->startNode->getIndex());
//             nodesInSubcircuit.insert(mp_circuitDataWithoutDiodes->components.at(startingComponent)->endNode->getIndex());
//             // Find all other components belonging to this subcircuit. We must iterate (while-loop) to keep
//             // finding new connections, until nothing changes on an iteration:
//             bool newComponentFoundOnLastIteration = true;
//             while (newComponentFoundOnLastIteration)
//             {
//                 newComponentFoundOnLastIteration = false;
//                 for (int potentialSubcircuitComponent = startingComponent+1; potentialSubcircuitComponent < mp_circuitDataWithoutDiodes->numberOfComponents; potentialSubcircuitComponent++)
//                 {
//                     if (componentAssignedToASubcircuit.at(potentialSubcircuitComponent) == false)
//                     {
//                         // get the nodes of the current component
//                         int componentStartNode = mp_circuitDataWithoutDiodes->components.at(potentialSubcircuitComponent)->startNode->getIndex();
//                         int componentEndNode = mp_circuitDataWithoutDiodes->components.at(potentialSubcircuitComponent)->endNode->getIndex();

//                         // Check whether either of these nodes connect potentialSubcircuitComponent to the current subcircuit.
//                         // If so, assign potentialSubcircuitComponent to the subcircuit.
//                         if(nodesInSubcircuit.count(componentStartNode) == 1 ||
//                             nodesInSubcircuit.count(componentEndNode) == 1) // the count() method is an unfortunate name; it can only return 0 or 1 because of std::set's unique element rule.
//                         {
//                             // Note that a re-parse is required, because the subcircuit grew:
//                             newComponentFoundOnLastIteration = true;

//                             m_AtomicSubcircuitsComponentsBelongsTo.at(potentialSubcircuitComponent) = subcircuitIndex;
//                             componentAssignedToASubcircuit.at(potentialSubcircuitComponent) = true;
                            
//                             // Record that these nodes belong to this subcircuit.
//                             // Note that the insert will simply do nothing in the case where the node
//                             // is already recorded as being part of the subcircuit. This is fine.
//                             nodesInSubcircuit.insert(componentStartNode);
//                             nodesInSubcircuit.insert(componentEndNode);
//                         }
//                     }
//                 }
//             }

//             subcircuitIndex++;
//         }
//     }
//     // We have now also counted the atomic subcircuits. Record that count!
//     m_NumberOfAtomicSubcircuits = subcircuitIndex;

//     // Sanity check, to make sure all the components got assigned to a subcircuit:
//     for (auto componentAssigned = componentAssignedToASubcircuit.begin(); componentAssigned!=componentAssignedToASubcircuit.end(); componentAssigned++)
//     {
//         assert(*componentAssigned);
//     }

// }

void NetlistCircuit::setPointersToBoundaryPressuresAndFlows(double* const interfacePressures, double* const interfaceFlows, const int& numberOfPointers)
{
    assert(pressure_n_ptrs.size()==0);
    assert(flow_n_ptrs.size()==0);
    for (int pointerIndex = 0; pointerIndex<numberOfPointers; pointerIndex++)
    {
        pressure_n_ptrs.push_back(&interfacePressures[pointerIndex]);
        flow_n_ptrs.push_back(&interfaceFlows[pointerIndex]);
    }
}

void NetlistCircuit::setInternalHistoryPressureFlowsAndVolumes()
{
    // for (auto node=mp_circuitDataWithoutDiodes.mapOfPressureNodes.begin(); node!=mp_circuitDataWithoutDiodes.mapOfPressureNodes.end(); node++)
    // Cycle and store the history pressures
    for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
    {
        if (node->second->hasHistoryPressure())
        {
            node->second->copyPressureToHistoryPressure();
        }
    }

    // for (auto component=mp_circuitDataWithoutDiodes.mapOfComponents.begin(); component!=mp_circuitDataWithoutDiodes.mapOfComponents.end(); component++)
    // Cycle and store the history flows
    for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    {
        if (component->second->hasHistoryFlow)
        {
            component->second->historyFlow = component->second->flow;
        }
    }

    // Store the volumes (currently just for VolumeTrackingPressureChambers. Make this more generic if new volume-tracking components are added later
    // - recommed using the hasHistoryVolume bool).
    for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    {
        VolumeTrackingComponent* pressureChamber = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
        // Ensure this actually is a VolumeTrackingComponent before going further:
        if (pressureChamber != NULL)
        {
            // Make the current volume into the new history volume:
            pressureChamber->cycleHistoryVolume();
        }
    }
}

void NetlistCircuit::recordPressuresFlowsAndVolumesInHistoryArrays()
{
    // for (auto node=mp_circuitDataWithoutDiodes.mapOfPressureNodes.begin(); node!=mp_circuitDataWithoutDiodes.mapOfPressureNodes.end(); node++)
    // Cycle and store the history pressures
    for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
    {
        // Store the pressure for writing to output file:
        node->second->appendToPressureHistory(node->second->getPressure());
    }

    // for (auto component=mp_circuitDataWithoutDiodes.mapOfComponents.begin(); component!=mp_circuitDataWithoutDiodes.mapOfComponents.end(); component++)
    // Cycle and store the history flows
    for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    {
        // Store the flow for writing to output file:
        component->second->appendToFlowHistory(component->second->flow);
    }

    // Store the volumes (currently just for VolumeTrackingPressureChambers. Make this more generic if new volume-tracking components are added later
    // - recommed using the hasHistoryVolume bool).
    for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    {
        VolumeTrackingComponent* pressureChamber = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
        // Ensure this actually is a VolumeTrackingComponent before going further:
        if (pressureChamber != NULL)
        {
            // Store the volume for writing to output file:
            pressureChamber->recordVolumeInHistory();
        }
    }

    setInternalHistoryPressureFlowsAndVolumes();
}

void NetlistCircuit::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    rebuildCircuitMetadata();
    detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    // cycleToSetHistoryPressuresFlowsAndVolumes();
}

void NetlistCircuit::finalizeLPNAtEndOfTimestep()
{
    recordPressureHistoryHistory();
    recordPressuresFlowsAndVolumesInHistoryArrays();
    switchDiodeStatesIfNecessary();
}

void NetlistCircuit::recordPressureHistoryHistory()
{
    for (auto node : mp_circuitData->mapOfPressureNodes)
    {
        // Store the pressure for writing to output file:
        node.second->copyHistoryPressureToHistoryHistoryPressure();
    }
}

void NetlistCircuit::recordPressureHistory()
{
    for (auto node : mp_circuitData->mapOfPressureNodes)
    {
        // Store the pressure for writing to output file:
        node.second->copyPressureToHistoryPressure();
    }
}

boost::shared_ptr<CircuitComponent> NetlistCircuit::getComponentByInputDataIndex(const int componentIndex)
{
    return mp_circuitData->getComponentByInputDataIndex(componentIndex);
}

boost::shared_ptr<CircuitPressureNode> NetlistCircuit::getNodeByInputDataIndex(const int componentIndex)
{
    return mp_circuitData->getNodeByInputDataIndex(componentIndex);
}

std::vector<std::pair<int,double*>> NetlistCircuit::getComponentInputDataIndicesAndFlows() const
{
    return mp_circuitData->getComponentInputDataIndicesAndFlows();
}

std::vector<std::pair<int,double*>> NetlistCircuit::getNodeInputDataIndicesAndPressures() const
{
    return mp_circuitData->getNodeInputDataIndicesAndPressures();
}

std::vector<std::pair<int,double*>> NetlistCircuit::getVolumeTrackingComponentInputDataIndicesAndVolumes() const
{
    return mp_circuitData->getVolumeTrackingComponentInputDataIndicesAndVolumes();
}

boost::shared_ptr<CircuitData> NetlistCircuit::getCircuitDescription()
{
    return mp_circuitData;
}

void NetlistCircuit::initialiseCircuit()
{
    // This function exists just so we can modify what initialiseCircuit does in subclasses without repeating code.
    mp_netlistFileReader = NetlistReader::Instance();
    mp_netlistXmlReader = NetlistXmlReader::Instance();
    initialiseCircuit_common();

    // The system is square in this case
    m_numberOfSystemRows = m_numberOfSystemColumns;

    createVectorsAndMatricesForCircuitLinearSystem();
}

void NetlistCircuit::initialiseCircuit_common()
{
  numberOfPrescribedPressuresAndFlows = mp_circuitData->numberOfPrescribedPressures + mp_circuitData->numberOfPrescribedFlows; // Just the sum of the previous two declared integers

  // Resize to contain the necessary flows and pressures, and initialise to zero:
  flowsInSubcircuit.resize(mp_circuitData->numberOfComponents,0.0);
  pressuresInSubcircuit.resize(mp_circuitData->numberOfPressureNodes,0.0);

  getMapOfPressHistoriesToCorrectPressNodes(); //initialises m_numberOfHistoryPressures
  getMapOfFlowHistoriesToCorrectComponents(); //initialises numberOfHistoryFlows
  getMapOfVolumeHistoriesToCorrectComponents(); // initialises numberOfHistoryVolumes
  getMapOfTrackedVolumesToCorrectComponents(); // initialises m_numberOfTrackedVolumes

  volumesInSubcircuit.resize(m_numberOfTrackedVolumes,0.0);

  m_numberOfSystemColumns = mp_circuitData->numberOfPressureNodes;
  m_numberOfSystemColumns += m_numberOfHistoryPressures;
  m_numberOfSystemColumns += mp_circuitData->numberOfComponents;
  m_numberOfSystemColumns += numberOfHistoryFlows;
  m_numberOfSystemColumns += m_numberOfTrackedVolumes;
  m_numberOfSystemColumns += numberOfHistoryVolumes;

  // columnMapSize = m_numberOfHistoryPressures + numberOfHistoryFlows + numberOfPrescribedPressures + numberOfPrescribedFlows;
  createListOfNodesWithMultipleIncidentCurrents();

  loadPressuresFlowsAndVolumesOnRestart();
  setInternalHistoryPressureFlowsAndVolumes();
  recordPressureHistoryHistory();
  if (m_startingTimestepIndex == 0)
  {
    m_oneshotIgnoreIncorrectFortranFlow = true;
    recordPressuresFlowsAndVolumesInHistoryArrays();
  }
}

int NetlistCircuit::getNumberOfHistoryPressures() const
{
    return m_numberOfHistoryPressures;
}

void NetlistCircuit::createVectorsAndMatricesForCircuitLinearSystem()
{
  PetscErrorCode errFlag;
  // Create m_systemMatrix and m_inverseOfSystemMatrix (to be filled later)
  errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_numberOfSystemRows,m_numberOfSystemColumns,NULL,&m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatZeroEntries(m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

  errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_numberOfSystemRows,m_numberOfSystemColumns,NULL,&m_inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatZeroEntries(m_inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  // To compute m_inverseOfSystemMatrix, we use Petsc to solve AX=B, with A,X and B all systemSize x systemSize matrices
  // B will just be m_identityMatrixForPetscInversionHack.
  errFlag = MatCreateSeqDense(PETSC_COMM_SELF,m_numberOfSystemRows,m_numberOfSystemColumns,NULL,&m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatZeroEntries(m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  // Fill the diagonal with ones:
  for (int ii=0; ii<m_numberOfSystemRows; ii++)
  {
      errFlag = MatSetValue(m_identityMatrixForPetscInversionHack,ii,ii,1.0,INSERT_VALUES);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  }
  errFlag = MatAssemblyBegin(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatAssemblyEnd(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);

  errFlag = VecCreate(PETSC_COMM_SELF,&m_RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecSetType(m_RHS,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_RHS a VECSEQ sequential vector
  errFlag = VecSetSizes(m_RHS,m_numberOfSystemRows,m_numberOfSystemRows); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecZeroEntries(m_RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyBegin(m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyEnd(m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

  errFlag = VecCreate(PETSC_COMM_SELF,&m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecSetType(m_solutionVector,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make m_solutionVector a VECSEQ sequential vector
  errFlag = VecSetSizes(m_solutionVector,m_numberOfSystemColumns,m_numberOfSystemColumns); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecZeroEntries(m_solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyBegin(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyEnd(m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistCircuit::createListOfNodesWithMultipleIncidentCurrents()
{
    // Note that this function also counts pressure nodes which are just
    // /between/ two components, eg. for the (two resistor) subcircuit:
    //        N0----[==R1==]----N1----[==R2==]----N2
    // This would count node N1 as appearing twice, and do a "Kirchoff" current
    // balance of the form "flow through R1 = flow through R2".
    //
    // It also catches and deals with true Kirchoff equations where a third (fourth, fifth,...)
    // component is connected to N1.

    int numberOfTimesNodeSeen;
    m_numberOfMultipleIncidentCurrentNodes = 0;

    // The node data from the input file is 1-indexed, so shift this to 1:numberOfPressureNodes, instead of 0:numberOfPressureNodes-1
    for(auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
    {
       int nodeIndex = node->second->getIndex();
       numberOfTimesNodeSeen = 0;
       for (int componentIndex = 0; componentIndex<mp_circuitData->numberOfComponents; componentIndex++)
       {
          if ((mp_circuitData->components.at(componentIndex)->startNode->getIndex() == nodeIndex) || 
              (mp_circuitData->components.at(componentIndex)->endNode->getIndex() == nodeIndex))
          {
             numberOfTimesNodeSeen++;
          }
       }
       if (numberOfTimesNodeSeen > 1)
       {
          if (kirchoffEquationAtNodeDeferredToInterfacingCircuit(nodeIndex))
          {
              m_nodesWithKirchoffEquationsDeferredToClosedLoop.push_back(nodeIndex);
          }
          else
          {
              listOfNodesWithMultipleIncidentCurrents.push_back(nodeIndex);
              m_numberOfMultipleIncidentCurrentNodes++;
          }
       }
    }

}

std::vector<int> NetlistCircuit::getNodesWithDeferredKirchoffEquations() const
{
    return m_nodesWithKirchoffEquationsDeferredToClosedLoop;
}

bool NetlistCircuit::kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const
{
    // In NetlistCircuit, there is no downstream circuit, so the return value is always false.
    return false;
}

void NetlistCircuit::getMapOfPressHistoriesToCorrectPressNodes()
{
    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
       // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
       if (mp_circuitData->components.at(ii)->getType() == Component_Capacitor)
       {
            listOfHistoryPressures.insert(mp_circuitData->components.at(ii)->startNode->getIndex());
            mp_circuitData->components.at(ii)->startNode->setHasHistoryPressure(true);
            
            listOfHistoryPressures.insert(mp_circuitData->components.at(ii)->endNode->getIndex());
            mp_circuitData->components.at(ii)->endNode->setHasHistoryPressure(true);
       }
    }

    m_numberOfHistoryPressures = listOfHistoryPressures.size();

    // Now do the actual generation of the pressure history node ordering map (for use when bulding the linear system matrix):
    int ii=0;
    for (auto iterator=listOfHistoryPressures.begin(); iterator != listOfHistoryPressures.end(); iterator++)
    {
       nodeIndexToPressureHistoryNodeOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
       ii++;
    }
}


void NetlistCircuit::getMapOfFlowHistoriesToCorrectComponents()
{

    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
       // Check for inductor, as these need flow "histories" (flow from the previous time-step) (for dQ/dt term).
     if(mp_circuitData->components.at(ii)->getType() == Component_Inductor)
       {
            listOfHistoryFlows.insert(ii);
        mp_circuitData->components.at(ii)->hasHistoryFlow = true;
       }
    }

    numberOfHistoryFlows = listOfHistoryFlows.size();

    // Now do the actual generation of the flow history component ordering map (for use when bulding the linear system matrix):
    int ii=0;
    for (auto iterator=listOfHistoryFlows.begin(); iterator != listOfHistoryFlows.end(); iterator++)
    {
       componentIndexToFlowHistoryComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
    }
}

void NetlistCircuit::getMapOfVolumeHistoriesToCorrectComponents()
{

  for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
  {
     // Check for VolumeTrackingPressureChambers, as these need volume "histories" (volume from the previous time-step) (for dVolume/dt term).
     if (boost::dynamic_pointer_cast<VolumeTrackingComponent> (mp_circuitData->components.at(ii)))
     {
        listOfHistoryVolumes.insert(ii);
        mp_circuitData->components.at(ii)->setHasHistoryVolume(true);
     }
  }

  numberOfHistoryVolumes = listOfHistoryVolumes.size();

  // Now do the actual generation of the pressure history node ordering map (for use when bulding the linear system matrix):
  int ii=0;
  for (auto iterator=listOfHistoryVolumes.begin(); iterator != listOfHistoryVolumes.end(); iterator++)
  {
     componentIndexToVolumeHistoryComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
  }
}

void NetlistCircuit::getMapOfTrackedVolumesToCorrectComponents()
{

  for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
  {
     // Check for VolumeTrackingPressureChambers, as these need volume tracking (for computing the current pressure, via the compliance/elastance).
     if (boost::dynamic_pointer_cast<VolumeTrackingComponent> (mp_circuitData->components.at(ii)))
     {
        listOfTrackedVolumes.insert(ii);
        // mp_circuitData->components.at(ii)->setHasTrackedVolume(true);
     }
  }

  m_numberOfTrackedVolumes = listOfTrackedVolumes.size();
  countVolumeTrackingPressureChambers();

  // Now do the actual generation of the volume node ordering map (for use when bulding the linear system matrix):
  int ii=0;
  for (auto iterator=listOfTrackedVolumes.begin(); iterator != listOfTrackedVolumes.end(); iterator++)
  {
     componentIndexToTrackedVolumeComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
  }
}

void NetlistCircuit::countVolumeTrackingPressureChambers()
{
    m_numberOfVolumeTrackingPressureChambers = 0;
    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
        // Check for VolumeTrackingPressureChambers, as these need volume tracking (for computing the current pressure, via the compliance/elastance).
        if (boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (mp_circuitData->components.at(ii)))
        {
            m_numberOfVolumeTrackingPressureChambers++;
        }
    }
}

void NetlistCircuit::generateLinearSystemFromPrescribedCircuit(const double alfi_delt)
{
    // Build m_systemMatrix
    generateLinearSystemWithoutFactorisation(alfi_delt);
    
    // LU factor m_systemMatrix
    PetscErrorCode errFlag = MatLUFactor(m_systemMatrix,NULL,NULL,NULL);CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

int NetlistCircuit::getSurfaceIndex() const
{
    return m_surfaceIndex;
}

void NetlistCircuit::generateLinearSystemWithoutFactorisation(const double alfi_delt)
{
    // This function assembles the system of (time-discretised) linear algebraic equations for the LPN.
    PetscErrorCode errFlag;

    m_locationsInRHSForUnstressedVolumesAndTheirValues.clear();

    errFlag = MatZeroEntries(m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    {
      int row = 0; // is the row in the matrix that we write to on each occasion
      for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
      {
          // bool componentIsOpenDiode = (mp_circuitData->components.at(ll)->type == Component_Diode &&
          //                              mp_circuitData->components.at(ll)->hasNonnegativePressureGradientAndNoBackflow());
          // open diodes are just implemented as zero-resistance resistors, closed diodes are zero-resistance resistors with prescribed zero flow
          if ((*component)->getType() == Component_Resistor)
          {
            // insert resistor relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            // std::cout << "current resistor value: " << currentParameterValue << std::endl;
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-currentParameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            row++;
          }
          else if  ((*component)->getType() == Component_Diode)
          {
            if ((*component)->permitsFlow()) // if the diode is open
            {
                // insert resistor(-eqsue) relationship into equation system for the open diode
                int startNode = (*component)->startNode->getIndex();
                errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                int endNode = (*component)->endNode->getIndex();
                errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

                double currentParameterValue = *((*component)->getParameterPointer());
                int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
                errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-currentParameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            else // the diode is closed
            {
                int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
                const double justAOneToImposeZeroFlow = 1.0;
                errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,justAOneToImposeZeroFlow,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);                
            }
            row++;
          }
          else if ((*component)->getType() == Component_Capacitor)
          {
            // insert capacitor relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            // std::cout << "current capacitor value: " << currentParameterValue << std::endl;
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-alfi_delt/currentParameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatSetValue(m_systemMatrix,row,nodeIndexToPressureHistoryNodeOrderingMap.at(startNode)+mp_circuitData->numberOfPressureNodes,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatSetValue(m_systemMatrix,row,nodeIndexToPressureHistoryNodeOrderingMap.at(endNode)+mp_circuitData->numberOfPressureNodes,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            row++;
          }
          else if ((*component)->getType() == Component_Inductor)
          {
            // insert inductor relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            
            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-currentParameterValue/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatSetValue(m_systemMatrix,row,componentIndexToFlowHistoryComponentOrderingMap.at(indexOfThisComponentsFlow)+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures+mp_circuitData->numberOfComponents,currentParameterValue/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            row++;
          }
          else if ((*component)->getType() == Component_VolumeTrackingPressureChamber)
          {
            // Two equations are needed for the VolumeTrackingPressureChamber:
            // 1) dVolume/dt = flow
            // 2) compliance * pressure = (volume - unstressed volume) ... unstressed vol will go on m_RHS of system, later.
            // Note that this means we increment row (row++) twice during this if-case
            // Do (1):
            // Insert the dt*flow term:
            int zeroIndexOfThisComponent = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,zeroIndexOfThisComponent+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // Insert the volume term (volume to-be-found during the next system solve):
            {
              int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            // Insert the volume history term:
            {
              int columnIndex = componentIndexToVolumeHistoryComponentOrderingMap.at(zeroIndexOfThisComponent) + m_numberOfTrackedVolumes + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            row++; // done twice in this if-case, because there are 2 equations to create for the VolumeTrackingPressureChamber
            
            // Now do (2) (see comment block above, within this if-case)
            // Do the compliance term:
            boost::shared_ptr<VolumeTrackingPressureChamber> volumeTrackingPressureChamber = boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (*component);
            if (!volumeTrackingPressureChamber->zeroVolumeShouldBePrescribed()) // test whether, on a previous attempt to solve the system, negative volumes were detected. If so, we'll do something different in the "else" below...
            {
              {
                int columnIndex = toZeroIndexing((*component)->startNode->getIndex());
                double valueToInsert = 1.0/(volumeTrackingPressureChamber->getElastance());
                errFlag = MatSetValue(m_systemMatrix,row,columnIndex,valueToInsert,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              }
              // Do the volume term:
              {
                int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
                errFlag = MatSetValue(m_systemMatrix,row,columnIndex,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              }
              // this is used to tell the RHS where to put the unstressed volume for this component:
              boost::shared_ptr<VolumeTrackingPressureChamber> downcastVolumeTrackingPressureChamber = boost::static_pointer_cast<VolumeTrackingPressureChamber> (*component);
              m_locationsInRHSForUnstressedVolumesAndTheirValues.push_back(std::make_pair(row, downcastVolumeTrackingPressureChamber->getUnstressedVolume() ) );
            }
            else //volumeTrackingPressureChamber->zeroVolumeShouldBePrescribed() == true, so instead of prescribing pressure based on volume, we allow the pressure to be a free variable, and prescribe zero volume for the chamber
            {
              // Prescribe zero volume:
              // VERY IMPORTANT: Note that we don't need to prescribe anything special on the m_RHS for this; the
              // prescribed volume is zero, we're doing a row for a component (which have zeros on the m_RHS anyway), and the
              // m_RHS is zeroed out before we start to fill it, so the zero will already be in place. But be aware of this
              // if you're making changes.
              int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // Reset the zero-volume marker on the component:p
              volumeTrackingPressureChamber->resetZeroVolumePrescription();
            }
            row++; // done twice in this if-case, because there are 2 equations to create for the VolumeTrackingPressureChamber
          }
          else if ((*component)->getType() == Component_VolumeTracking)
          {
            // Two equations are needed for the VolumeTrackingComponent:
            // 1) dVolume/dt = flow
            // Do (1):
            // Insert the dt*flow term:
            int zeroIndexOfThisComponent = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,zeroIndexOfThisComponent+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // Insert the volume term (volume to-be-found during the next system solve):
            {
              int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            // Insert the volume history term:
            {
              int columnIndex = componentIndexToVolumeHistoryComponentOrderingMap.at(zeroIndexOfThisComponent) + m_numberOfTrackedVolumes + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            row++;
          }
          else
          {
            throw std::runtime_error("EE: Unknown component type in netlist. Halting.");
          }
      }
    }

     // Do the equations for the nodes with multiple incident currents (Kirchoff):
     for (int mm = 0; mm < m_numberOfMultipleIncidentCurrentNodes; mm++)
     {
       for (int ll=0; ll<mp_circuitData->numberOfComponents; ll++)
       {
          bool foundMultipleIncidentCurrentsForEndNode = (mp_circuitData->components.at(ll)->endNode->getIndex() == listOfNodesWithMultipleIncidentCurrents.at(mm)); 
          if (foundMultipleIncidentCurrentsForEndNode)
          {
            int row = mm + mp_circuitData->numberOfComponents + m_numberOfVolumeTrackingPressureChambers;
            int column = ll + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures;
            errFlag = MatSetValue(m_systemMatrix,row,column,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }

          bool foundMultipleIncidentCurrentsForStartNode = (mp_circuitData->components.at(ll)->startNode->getIndex() == listOfNodesWithMultipleIncidentCurrents.at(mm));
          if (foundMultipleIncidentCurrentsForStartNode)
          {
            int row = mm + mp_circuitData->numberOfComponents + m_numberOfVolumeTrackingPressureChambers;
            int column = ll + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures;
            errFlag = MatSetValue(m_systemMatrix,row,column,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
       }
     }

     int rowsDoneSoFar = mp_circuitData->numberOfComponents + m_numberOfMultipleIncidentCurrentNodes + m_numberOfVolumeTrackingPressureChambers;

     // create the columnMap which tells us which system column each of the prescribed pressure, pressure-history, flow, or volume values belong to
     columnMap.clear();
     int tempUnknownVariableIndexWithinLinearSystem = 0; // just an indexing shift to keep track of where we need to write next
     // Do the prescribed pressures:
     for (auto prescribedPressureNode=mp_circuitData->mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=mp_circuitData->mapOfPrescribedPressureNodes.end(); prescribedPressureNode++)
     {
        columnMap.push_back(toZeroIndexing(prescribedPressureNode->second->getIndex()));
     }

     // Do the history pressures
     tempUnknownVariableIndexWithinLinearSystem += mp_circuitData->numberOfPressureNodes; // tempUnknownVariableIndexWithinLinearSystem is zero before this line; I'm doing it like this for clarity & consistency
     for (int ll=0; ll<m_numberOfHistoryPressures; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     // Do the prescribed flows
     tempUnknownVariableIndexWithinLinearSystem += m_numberOfHistoryPressures;
     for (auto prescribedFlowComponent=mp_circuitData->mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=mp_circuitData->mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
     {
        columnMap.push_back(toZeroIndexing(prescribedFlowComponent->second->getIndex())+tempUnknownVariableIndexWithinLinearSystem);
     }

     // Do the history flows
     tempUnknownVariableIndexWithinLinearSystem += mp_circuitData->numberOfComponents;
     for (int ll=0; ll<numberOfHistoryFlows; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     // Do the history volumes
     tempUnknownVariableIndexWithinLinearSystem += numberOfHistoryFlows + m_numberOfTrackedVolumes;
     for (int ll=0; ll<numberOfHistoryVolumes; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     if (columnMap.size() != m_numberOfSystemRows - rowsDoneSoFar)
     {
        throw std::runtime_error("EE: Netlist circuit appears to be malformed.");
     }

     // Set the prescribed-value equations (i.e. pressure_1 (LHS) = pressure_1 (m_RHS) - so really just a way of setting the prescribed values within the linear system)
     for (int ll = 0; ll < m_numberOfSystemRows - rowsDoneSoFar; ll++)
     {
       errFlag = MatSetValue(m_systemMatrix,rowsDoneSoFar + ll, columnMap.at(ll),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     }

     errFlag = MatAssemblyBegin(m_systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     errFlag = MatAssemblyEnd(m_systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    // std::cout << "System matrix for surface " << m_surfaceIndex << ":" << std::endl;
    // errFlag = MatView(m_systemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistCircuit::findLinearSystemIndicesOf3DInterfacePressureAndFlow()
{   
    // Find the location of the 3D interface pressure in the linear system:
    m_columnOf3DInterfacePrescribedPressureInLinearSystem.clear();
    m_locationOf3DInterfaceComputedPressureInSolutionVector.clear();
    if (mp_circuitData->hasPrescribedPressureAcrossInterface())
    {// m_numberOfTrackedVolumes swapped for m_numberOfVolumeTrackingPressureChambers
        const int tempIndexingShift = mp_circuitData->numberOfComponents + m_numberOfMultipleIncidentCurrentNodes + m_numberOfVolumeTrackingPressureChambers;
        int ll=0;
        for (auto prescribedPressureNode=mp_circuitData->mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=mp_circuitData->mapOfPrescribedPressureNodes.end(); prescribedPressureNode++ )
        {
            // Annotate the location of the 3D interface pressure in the linear system
            // (regardless of whether it is computed or prescribed)
            if (prescribedPressureNode->second->isAtBoundary())
            {
                m_columnOf3DInterfacePrescribedPressureInLinearSystem.push_back(ll + tempIndexingShift);
            }
            ll++;
        }
    }

    // Get the locations of the computed 3D interface pressure (which will just be a copy of an imposition in some cases):    
    for (auto pressureNode = mp_circuitData->mapOfPressureNodes.begin(); pressureNode != mp_circuitData->mapOfPressureNodes.end(); pressureNode++)
    {
        if (pressureNode->second->isAtBoundary())
        {
            m_locationOf3DInterfaceComputedPressureInSolutionVector.push_back(toZeroIndexing(pressureNode->second->getIndex()));
        }
    }
    

    // Find the location of the 3D interface flow in the linear system:
    m_columnOf3DInterfacePrescribedFlowInLinearSystem.clear();
    m_locationOf3DInterfaceComputedFlowInSolutionVector.clear();
    if (mp_circuitData->hasPrescribedFlowAcrossInterface())
    {//\friday m_numberOfTrackedVolumes swapped for m_numberOfVolumeTrackingPressureChambers
        const  int tempIndexingShift = mp_circuitData->numberOfComponents + m_numberOfMultipleIncidentCurrentNodes + m_numberOfVolumeTrackingPressureChambers +
                                        mp_circuitData->numberOfPrescribedPressures +
                                        m_numberOfHistoryPressures;

        int ll=0;
        for (auto prescribedFlowComponent=mp_circuitData->mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=mp_circuitData->mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
        {
            if (prescribedFlowComponent->second->connectsToNodeAtInterface())
            {
                m_columnOf3DInterfacePrescribedFlowInLinearSystem.push_back(ll + tempIndexingShift);
            }
            ll++;
        }
    }

    // Get the locations of the computed 3D interface flow (which will just be a copy of an imposition in some cases):    
    for (auto component = mp_circuitData->mapOfComponents.begin(); component != mp_circuitData->mapOfComponents.end(); component++)
    {
        if (component->second->connectsToNodeAtInterface())
        {
            const int index = toZeroIndexing(component->second->getIndex()) + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures;
            m_locationOf3DInterfaceComputedFlowInSolutionVector.push_back(index);
        }
    }

    // std::cout << "sizes are: " << columnIndexOf3DInterfacePressureInLinearSystem.size() << " " << columnIndexOf3DInterfaceFlowInLinearSystem.size() << std::endl;
}

// useHistoryHistoryPressure is for the Kalman filter, so we can update the Xn-1 variables to match
// the current Kalman particle
void NetlistCircuit::assembleRHS(const int timestepNumber, const bool useHistoryHistoryPressure)
{

    PetscErrorCode errFlag;
    errFlag = VecZeroEntries(m_RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // columnIndexOf3DInterfacePressureInLinearSystem.clear(); // dummy value, to be replaced!
    findLinearSystemIndicesOf3DInterfacePressureAndFlow();

    // Prescribed pressures
    int tempIndexingShift = mp_circuitData->numberOfComponents + m_numberOfMultipleIncidentCurrentNodes + m_numberOfVolumeTrackingPressureChambers;
    {
      int ll=0;
      for (auto prescribedPressureNode = mp_circuitData->mapOfPrescribedPressureNodes.begin(); prescribedPressureNode != mp_circuitData->mapOfPrescribedPressureNodes.end(); prescribedPressureNode++ )
      {
        // Coming from 'f' for 'fixed' in the input data:
        if (prescribedPressureNode->second->getPressurePrescriptionType() == Pressure_Fixed)
        {
            errFlag = VecSetValue(m_RHS, ll + tempIndexingShift, prescribedPressureNode->second->getPressure(), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);    
        }
        // Coming from 'l' for 'left-ventricular' in the input data:
        else if (prescribedPressureNode->second->getPressurePrescriptionType() == Pressure_LeftVentricular)
        {
            std::cerr << "this requires heart model. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
            std::exit(1);
        }
        else if (prescribedPressureNode->second->getPressurePrescriptionType() == Pressure_3DInterface)
        {
            // We only do this if the netlist is in Dirichlet BC mode:
            if (mp_circuitData->hasPrescribedPressureAcrossInterface())
            {
              try {
                double* pressurePointerToSet = pressure_n_ptrs.at(prescribedPressureNode->second->getPrescribedPressurePointerIndex());
                errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,*pressurePointerToSet,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              } catch (const std::exception& e) {
                std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
                throw;
              }
            }
        }
        else
        {
              throw std::runtime_error("Unknown pressure prescription value in Netlist.");
        }

        ll++;
      }
    }

    // History Pressures
    tempIndexingShift += mp_circuitData->numberOfPrescribedPressures;

    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
        int lll=0;
        for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
        {
            if (node->second->hasHistoryPressure())
            {
                double historyPressure;
                if (useHistoryHistoryPressure)
                {
                    historyPressure = node->second->getHistoryHistoryPressure();
                } else {
                    historyPressure = node->second->getHistoryPressure();
                }
                errFlag = VecSetValue(m_RHS,lll+tempIndexingShift,historyPressure,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
                lll++;
            }
        }
    }
    // Prescribed Flows:
    tempIndexingShift += m_numberOfHistoryPressures;
    // Scoping unit to include the second counter ll in the for loop, without having ll in-scope after the loop finishes:
    {
        int ll=0;
        for (auto prescribedFlowComponent=mp_circuitData->mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=mp_circuitData->mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
        {
            // std::cout << "surface " << m_IndexOfThisNetlistLPNInInputFile << "prescribed flow component " << prescribedFlowComponent->second->prescribedFlowType << std::endl;//++++
           if (prescribedFlowComponent->second->prescribedFlowType == Flow_3DInterface)
           {
            if (mp_circuitData->hasPrescribedFlowAcrossInterface())
            {
              double* flowPointerToSet;
              try {
                flowPointerToSet = flow_n_ptrs.at(prescribedFlowComponent->second->prescribedFlowPointerIndex);
              } catch (const std::exception& e) {
                std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
                throw;
              }
              // First, flip the sign of the flow, if necessary due to the orientation of the component at the 3D interface:
              double threeDFlowValue;
              if (m_oneshotIgnoreIncorrectFortranFlow)
              {
                // We know that we are in a restarted simulation, on the first
                // restarted step. In this case, the flow from Fortran will be
                // incorrect, so we use the value saved by the boundary condition
                // at the end of the last simulation.
                threeDFlowValue = prescribedFlowComponent->second->flow * mp_circuitData->getSignForPrescribed3DInterfaceFlow();
                m_oneshotIgnoreIncorrectFortranFlow = false;
              }
              else
              {
                threeDFlowValue = *flowPointerToSet * mp_circuitData->getSignForPrescribed3DInterfaceFlow();
              }
              assert(!isnan(threeDFlowValue));
              // Give the (possibly sign-corrected) flow to the linear system:
              errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,threeDFlowValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // columnIndexOf3DInterfaceFlowInLinearSystem.push_back(ll + tempIndexingShift);
              // nextFlowPointerIndex++;
            }
           }
           else if (prescribedFlowComponent->second->prescribedFlowType == Flow_Fixed)
           {
              errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,prescribedFlowComponent->second->getPrescribedFlow(), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
           }
         // else if (prescribedFlowComponent->second->prescribedFlowType == Flow_Diode_FixedWhenClosed)
         // {
         //    // Check whether the diode is closed; if so, prescribe the (zero) flow:
         //    if (!prescribedFlowComponent->second->hasNonnegativePressureGradientAndNoBackflow())
         //    {
         //      errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,prescribedFlowComponent->second->valueOfPrescribedFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
         //    }
         // }
           else
           {
                throw std::runtime_error("Unknown flow prescription value in Netlist.");
           }

         // if (prescribedFlowComponent->second->connectsToNodeAtInterface())
         // {
         //    columnIndexOf3DInterfaceFlowInLinearSystem.push_back(ll + tempIndexingShift);
         // }

           ll++;
        }
      }
    // for(int ll=0; ll<numberOfPrescribedFlows; ll++)
    // {
    //    if (subcircuitInputData.typeOfPrescribedFlows.at(ll) == Flow_3DInterface)
    //    {
    //       columnIndexOf3DInterfaceFlowInLinearSystem = ll + tempIndexingShift;
    //       errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,*flow_n_ptr,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //    }
    //    else if (subcircuitInputData.typeOfPrescribedFlows.at(ll) == Flow_Fixed)
    //    {
    //       errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,valueOfPrescribedFlows.at(ll), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //    }
    //    else
    //    {
    //         throw std::runtime_error("Unknown flow prescription value in Netlist.");
    //    }
    // }
    // History Flows
    tempIndexingShift = tempIndexingShift + mp_circuitData->numberOfPrescribedFlows;
    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
        int lll=0;
        for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
        {
            if (component->second->hasHistoryFlow)
            {
                errFlag = VecSetValue(m_RHS,lll + tempIndexingShift,component->second->historyFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
                lll++;      
            }
        }
        // for (auto iterator=listOfHistoryFlows.begin(); iterator!=listOfHistoryFlows.end(); iterator++)
        // {
        //    errFlag = VecSetValue(m_RHS,lll + tempIndexingShift,historyFlowsInSubcircuit.at(*iterator - 1), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // }
        // lll++;
      }

    // Give the linear system the history volumes, by putting them on the m_RHS:
    tempIndexingShift += numberOfHistoryFlows;
    {
      int ll=0;
      for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
      {
        if ((*component)->getHasHistoryVolume())
        {
          // currently, only VolumeTrackingComponents have history volumes. We might want to change this cast later, if new component types
          // with history volumes get added.
          boost::shared_ptr<VolumeTrackingComponent> volumeTrackingComponent = boost::dynamic_pointer_cast<VolumeTrackingComponent> (*component);
          double volume = volumeTrackingComponent->getHistoryVolume();
          int row = ll + tempIndexingShift;
          errFlag = VecSetValue(m_RHS,row,volume,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          ll++;
        }
      }
    }

    // Finally, add the unstressed volumes:
    for (auto rowIndexAndUnstressedVolumePair = m_locationsInRHSForUnstressedVolumesAndTheirValues.begin(); rowIndexAndUnstressedVolumePair != m_locationsInRHSForUnstressedVolumesAndTheirValues.end(); rowIndexAndUnstressedVolumePair++)
    {
        int row = rowIndexAndUnstressedVolumePair->first;
        double unstressedVolume = rowIndexAndUnstressedVolumePair->second;
        // minus unstressedVolume because compliance * pressure = (volume - unstressed volume).
        errFlag = VecSetValue(m_RHS, row, -unstressedVolume, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF, errFlag);
    }
 
    errFlag = VecAssemblyBegin(m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // if (m_surfaceIndex == -1){
    //     // std::cout << "m_RHS for surface " << m_surfaceIndex << ":" << std::endl;
    //     // errFlag = VecView(m_RHS,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //     // std::cout << "END m_RHS for surface " << m_surfaceIndex << std::endl;
    // }

}

void NetlistCircuit::updateLPN(const int timestepNumber)
{
    buildAndSolveLinearSystem(timestepNumber, m_delt);

    // Get the updated nodal pressures:
    giveNodesTheirPressuresFromSolutionVector();

    // Get the updated component flows:
    giveComponentsTheirFlowsFromSolutionVector();

    // Get the updated volumes:
    giveComponentsTheirVolumesFromSolutionVector();
}

void NetlistCircuit::giveNodesTheirPressuresFromSolutionVector()
{
  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;
  // std::cout << "called giveNodesTheirPressuresFromSolutionVector" << std::endl;

  // Look the nodes, handing them their new pressures from the circuit linear system solve:
  for (int pressureNodeIndex=0; pressureNodeIndex<mp_circuitData->numberOfPressureNodes; pressureNodeIndex++)
  {
      errFlag = VecGetValues(m_solutionVector,getSingleValue,&pressureNodeIndex,&pressuresInSubcircuit[pressureNodeIndex]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      // std::cout << "system matrix: " << m_surfaceIndex << std::endl;
      // MatView(m_systemMatrix,PETSC_VIEWER_STDOUT_WORLD);
      // std::cout << "m_RHS: " << std::endl;
      // VecView(m_RHS,PETSC_VIEWER_STDOUT_WORLD);
      // std::cout << "solution vec: " << std::endl;
      // VecView(m_solutionVector,PETSC_VIEWER_STDOUT_WORLD);
      assert(!isnan(pressuresInSubcircuit[pressureNodeIndex]));
      mp_circuitData->mapOfPressureNodes.at(toOneIndexing(pressureNodeIndex))->setPressure(pressuresInSubcircuit[pressureNodeIndex]);
      // std::cout << "just set pressure " << pressuresInSubcircuit[ll] << " at node " << ll << " of circuit " << m_surfaceIndex << std::endl;
  }
}

void NetlistCircuit::giveComponentsTheirFlowsFromSolutionVector()
{
  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;

  int firstFlowIndex = mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures;
  for (int ll=firstFlowIndex; ll<mp_circuitData->numberOfComponents+firstFlowIndex; ll++)
  {
      errFlag = VecGetValues(m_solutionVector,getSingleValue,&ll,&flowsInSubcircuit[ll-firstFlowIndex]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      assert(!isnan(flowsInSubcircuit[ll-firstFlowIndex]));
      mp_circuitData->mapOfComponents.at(toOneIndexing(ll-firstFlowIndex))->flow = flowsInSubcircuit[ll-firstFlowIndex];
  }
}

void NetlistCircuit::giveComponentsTheirVolumesFromSolutionVector()
{
  std::vector<double> volumes = getVolumesFromSolutionVector();
  // Reverse the volumes so we can pop off the back of it as we loop the mapOfVolumeTrackingComponents:
  std::reverse(volumes.begin(), volumes.end());
  for (auto component = mp_circuitData->mapOfVolumeTrackingComponents.begin(); component != mp_circuitData->mapOfVolumeTrackingComponents.end(); component++)
  {
      VolumeTrackingComponent* currentPressureChamber = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
      
      // Ensure we aren't dealing with negative volumes:
      //\todo REINSTATE!
      // assert(volumes.back() >= 0.0);
      assert(!isnan(volumes.back()));

      currentPressureChamber->setStoredVolume(volumes.back());
      volumes.pop_back();
  }

  // Ensure we've used all the retrieved volumes:
  assert(volumes.size() == 0);
}

// Compare with giveComponentsTheirVolumesFromSolutionVector. That function sets final, accepted volumes,
// whereas this function, giveComponentsTheirProposedVolumesFromSolutionVector, gives them /proposed/ volumes
// which are then checked for validity (i.e. non-negativity). Any proposed negative values trigger a re-computation
// of the solution to the circuit linear system, with an enforced zero-volume at the would-be negative volume locations.
void NetlistCircuit::giveComponentsTheirProposedVolumesFromSolutionVector()
{
  std::vector<double> volumes = getVolumesFromSolutionVector();
  // Reverse the volumes so we can pop off the back of it as we loop the mapOfVolumeTrackingComponents:
  std::reverse(volumes.begin(), volumes.end());
  for (auto component = mp_circuitData->mapOfVolumeTrackingComponents.begin(); component != mp_circuitData->mapOfVolumeTrackingComponents.end(); component++)
  {
      VolumeTrackingComponent* currentPressureChamber = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
      
      currentPressureChamber->setProposedVolume(volumes.back());
      volumes.pop_back();
  }

  // Ensure we've used all the retrieved volumes:
  assert(volumes.size() == 0);
}

std::vector<double> NetlistCircuit::getVolumesFromSolutionVector()
{
  std::vector<double> volumesToReturn;

  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;

  // The location of the first volume in the m_solutionVector:
  int firstVolumeIndex = mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->numberOfComponents + numberOfHistoryFlows;
  int volumeIndex = firstVolumeIndex;

  auto component = mp_circuitData->mapOfVolumeTrackingComponents.begin();
  while(component != mp_circuitData->mapOfVolumeTrackingComponents.end())
  {
    errFlag = VecGetValues(m_solutionVector,getSingleValue,&volumeIndex,&volumesInSubcircuit[volumeIndex-firstVolumeIndex]); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    VolumeTrackingComponent* currentPressureChamber = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
    volumesToReturn.push_back(volumesInSubcircuit[volumeIndex-firstVolumeIndex]);

    volumeIndex++;
    component++;
  }

  return volumesToReturn;
}

void NetlistCircuit::buildAndSolveLinearSystemForUpdatingHistoryVariablesToMatchCurrentKalmanParticle(const int timestepNumber, const double alfi_delt)
{
    generateLinearSystemFromPrescribedCircuit(alfi_delt);
    assembleRHS(timestepNumber, true);

    solveLinearSystem();
}

void NetlistCircuit::buildAndSolveLinearSystem(const int timestepNumber, const double alfi_delt)
{
  generateLinearSystemFromPrescribedCircuit(alfi_delt);
  assembleRHS(timestepNumber, false);

  solveLinearSystem();
}

void NetlistCircuit::solveLinearSystem()
{
  PetscErrorCode errFlag;
  // get the inverse of the system matrix:
  errFlag = MatMatSolve(m_systemMatrix,m_identityMatrixForPetscInversionHack,m_inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  // Release the m_systemMatrix so we can edit it again on the next iteration (we only need the just-computed m_inverseOfSystemMatrix for computations on this step now.)
  errFlag = MatSetUnfactored(m_systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

  // Solve the system
  errFlag = MatMult(m_inverseOfSystemMatrix,m_RHS,m_solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

std::pair<boundary_data_t,double> NetlistCircuit::computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber)
{
  buildAndSolveLinearSystem(timestepNumber,m_delt);

  PetscErrorCode errFlag;

  int numberOfValuesToGet=1;
  if (mp_circuitData->hasPrescribedFlowAcrossInterface()) // This boundary condition is receiving flow and returning pressure (Neumann mode if NetlistSubcircuit is a boundary condition)
  {
    int locationOfPressureInRHS = toZeroIndexing(mp_circuitData->getIndexOfNodeAtInterface());
    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,&locationOfPressureInRHS,&m_interfacePressure); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    return std::make_pair(Boundary_Pressure,m_interfacePressure);
  }
  else if (mp_circuitData->hasPrescribedPressureAcrossInterface()) // This boundary condition is receiving pressure and returning flow (Dirichlet mode if NetlistSubcircuit is a boundary condition)
  {
    int locationOfFlowInRHS = toOneIndexing(mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures + mp_circuitData->getIndexOfComponentConnectingToNodeAtInterface());

    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,&locationOfFlowInRHS,&m_interfaceFlow);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    return std::make_pair(Boundary_Flow,m_interfaceFlow);
  }
  else
  {
    std::stringstream errorMessage;
    errorMessage << "EE: Internal error when computing flow or pressure to pass to zero-D domain replacement." << std::endl;
    throw std::logic_error(errorMessage.str());
  }
}


std::pair<double,double> NetlistCircuit::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    assert(mp_circuitData->connectsTo3DDomain());

    buildAndSolveLinearSystem(timestepNumber,alfi_delt);

    // Extract the implicit coeffcients, for eventual passing to the FORTRAN
    // linear solve
    int rowToGet[] = {0};
    int numberOfValuesToGet=1;
    PetscScalar valueFromInverseOfSystemMatrix;
    
    std::pair<double,double> implicitCoefficientsToReturn;

    PetscErrorCode errFlag;
    errFlag = MatGetValues(m_inverseOfSystemMatrix,numberOfValuesToGet,rowToGet,numberOfValuesToGet,&m_columnOf3DInterfacePrescribedFlowInLinearSystem.at(0),&valueFromInverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    implicitCoefficientsToReturn.first = mp_circuitData->getSignForPrescribed3DInterfaceFlow() * valueFromInverseOfSystemMatrix;

    PetscScalar valueFromRHS;
    errFlag = VecGetValues(m_RHS,numberOfValuesToGet,&m_columnOf3DInterfacePrescribedFlowInLinearSystem.at(0),&valueFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    PetscScalar valueFromSolutionVector;
    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,rowToGet,&valueFromSolutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    implicitCoefficientsToReturn.second = valueFromSolutionVector - valueFromInverseOfSystemMatrix * valueFromRHS;//\todo make dynamic
    
    // std::cout << "m_inverseOfSystemMatrix: "<< std::endl;
    // MatView(m_inverseOfSystemMatrix,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "solution vector: "<< std::endl;
    // VecView(m_solutionVector,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "m_RHS vector: "<< std::endl;
    // VecView(m_RHS,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << std::setprecision(20) << "and just set " << implicitCoefficientsToReturn.first << " " <<implicitCoefficientsToReturn.second << std::endl;


    return implicitCoefficientsToReturn;
}

double NetlistCircuit::getInterfaceFlowSign() const
{
    return mp_circuitData->getSignForPrescribed3DInterfaceFlow();
}

void NetlistCircuit::computeHistoryVariablesToMatchCurrentKalmanFilterParticle(const int timestepNumber, const double alfi_delt)
{
    buildAndSolveLinearSystemForUpdatingHistoryVariablesToMatchCurrentKalmanParticle(timestepNumber, alfi_delt);
    giveNodesTheirPressuresFromSolutionVector();
    recordPressureHistory();
}

// This subroutine detects whether the last circuit linear system solve was invalid due to its producing negative
// volumes. The returned bool can be used to enforce a re-solve, with any negative pressures re-prescribed to be zero.
bool NetlistCircuit::areThereNegativeVolumes(const int timestepNumber, const double alfi_delt)
{
  buildAndSolveLinearSystem(timestepNumber, alfi_delt);
  // These volumes are "proposed", because if any are negative, we 
  // re-solve with zero-volume prescribed
  giveComponentsTheirProposedVolumesFromSolutionVector();

  bool thereAreNegativeVolumes = false;
  //\todo REINSTATE THIS LOOP!
  // for (auto component = mp_circuitData->mapOfVolumeTrackingComponents.begin(); component!=mp_circuitData->mapOfVolumeTrackingComponents.end(); component++)
  // {
  //   VolumeTrackingComponent* volumeTrackingPressureChamber = dynamic_cast<VolumeTrackingComponent*> (component->second.get());
  //   // Check for negative volumes:
  //   if (volumeTrackingPressureChamber->getProposedVolume() < 0.0)
  //   {
  //     thereAreNegativeVolumes = true;
  //     // Note that we're going to have to re-compute the circuit linear system solution, but this time
  //     // with the volume at this location enforced to zero (in order to avoid the negativity)
  //     volumeTrackingPressureChamber->enforceZeroVolumePrescription();
  //   }
  // }

  return thereAreNegativeVolumes;
}

bool NetlistCircuit::hasPrescribedPressureAcross3DInterface() const
{
    return mp_circuitData->hasPrescribedPressureAcrossInterface();
}

bool NetlistCircuit::hasPrescribedFlowAcross3DInterface() const
{
    return mp_circuitData->hasPrescribedFlowAcrossInterface();
}

std::vector<double*> NetlistCircuit::getCapacitorNodalHistoryPressurePointers() const
{
    return mp_circuitData->getCapacitorNodalHistoryPressurePointers();
}
