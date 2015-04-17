#include "NetlistCircuit.hxx"
#include "fileReaders.hxx"
#include "fileWriters.hxx"
#include <boost/make_shared.hpp>
#include "indexShifters.hxx"

// Statics:
int NetlistClosedLoopDownstreamCircuit::s_numberOfDownstreamCircuits = 0;

void NetlistCircuit::initialisePetscArrayNames()
{
    m_RHS = PETSC_NULL;
    m_solutionVector = PETSC_NULL;
    m_systemMatrix = PETSC_NULL;
    m_inverseOfSystemMatrix = PETSC_NULL;
    m_identityMatrixForPetscInversionHack = PETSC_NULL;
}

int NetlistCircuit::getNumberOfDegreesOfFreedom()
{
    return mp_circuitData->numberOfComponents + mp_circuitData->numberOfPressureNodes;
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
    createBasicCircuitDescription();
}

void NetlistCircuit::createBasicCircuitDescription()
{
    mp_circuitData->numberOfComponents = mp_netlistFileReader->getNumberOfComponents().at(m_IndexOfThisNetlistLPN);
    mp_circuitData->numberOfPressureNodes = mp_netlistFileReader->getNumberOfPressureNodes().at(m_IndexOfThisNetlistLPN);
    mp_circuitData->numberOfPrescribedPressures = mp_netlistFileReader->getNumberOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    mp_circuitData->numberOfPrescribedFlows = mp_netlistFileReader->getNumberOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);

    std::vector<circuit_component_t> retrievedComponentTypes = mp_netlistFileReader->getComponentTypes().at(m_IndexOfThisNetlistLPN);

    // Prepare space for the components in the circuit:
    assert(mp_circuitData->components.empty());
    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
        CircuitComponent* toPushBack;
        if (retrievedComponentTypes.at(ii) == Component_VolumeTrackingPressureChamber)
        {
            toPushBack = new VolumeTrackingPressureChamber(m_hstep,m_thisIsARestartedSimulation);
        }
        else
        {
            toPushBack = new CircuitComponent(m_hstep,m_thisIsARestartedSimulation);
        }

        mp_circuitData->components.push_back(boost::shared_ptr<CircuitComponent> (toPushBack));
        mp_circuitData->components.back()->setIndex(toOneIndexing(ii));
    }


    // Obtain the component- and node-level data for the circuit, for moving into the appropriate data structure CircuitData
    // We want to pop off the component types as we use them, but starting from the beginning of the vector. To do this, we reverse
    // the vector and then pop from the new end.
    std::reverse(retrievedComponentTypes.begin(), retrievedComponentTypes.end());
    std::vector<int> retrievedComponentStartNodes = mp_netlistFileReader->getComponentStartNodes().at(m_IndexOfThisNetlistLPN);
    std::reverse(retrievedComponentStartNodes.begin(), retrievedComponentStartNodes.end());
    std::vector<int> retrievedComponentEndNodes = mp_netlistFileReader->getComponentEndNodes().at(m_IndexOfThisNetlistLPN);
    std::reverse(retrievedComponentEndNodes.begin(), retrievedComponentEndNodes.end());
    std::vector<double> retrievedComponentParameterValues = mp_netlistFileReader->getComponentParameterValues().at(m_IndexOfThisNetlistLPN);
    std::reverse(retrievedComponentParameterValues.begin(), retrievedComponentParameterValues.end());

    std::vector<int> retrievedListOfPrescribedFlows = mp_netlistFileReader->getListOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::vector<circuit_component_flow_prescription_t> retrievedTypeOfPrescribedFlows = mp_netlistFileReader->getTypeOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::vector<double> retrievedValueOfPrescribedFlows = mp_netlistFileReader->getValueOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::map<int,double> retrievedInitialPressures = mp_netlistFileReader->getInitialPressures().at(m_IndexOfThisNetlistLPN);
    
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
                (*component)->valueOfPrescribedFlow = retrievedValueOfPrescribedFlows.at(prescribedFlow);
            }
        }

        (*component)->setParameterValue(retrievedComponentParameterValues.back());
        // make a copy of this value so we can reset it if necessary. Used for e.g. diode state changes.
        (*component)->parameterValueFromInputData = *((*component)->getParameterPointer());
        retrievedComponentParameterValues.pop_back();

        (*component)->startNode->setPressure(retrievedInitialPressures.at((*component)->startNode->getIndex()));
        (*component)->endNode->setPressure(retrievedInitialPressures.at((*component)->endNode->getIndex()));


    }

    // Some metadata is already set-up for this circuit, as a side-effect
    // of the above construction. This call completes the metadata; it's not
    // a problem that it also re-writes some of the existing metadata
    // (rewrites - but does not change - the values are identical!)
    mp_circuitData->rebuildCircuitMetadata();

    // Tell the node at the 3D interface that it connects to the 3D domain:
    {
        int threeDNodeIndex = mp_netlistFileReader->getIndicesOfNodesAt3DInterface().at(m_IndexOfThisNetlistLPN);
        mp_circuitData->initialiseNodeAndComponentAtInterface(threeDNodeIndex);
    }

    // mp_circuitData.switchDiodeStatesIfNecessary();
    // mp_circuitData.detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    // // Component indices are just consecutive integers by default, but sometimes non-consecutive numbering
    // // is needed; componentIndices allows for this.
    // // We initialise it now for the default case.
    // for (int ii=1; ii < mp_circuitData.numberOfComponents + 1; ii++)
    // {
    //     mp_circuitData.componentIndices.push_back(ii);
    // }
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
    mp_netlistFileReader = NetlistReader::Instance();
    std::vector<int> retrievedListOfPrescribedPressures = mp_netlistFileReader->getListOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    std::vector<circuit_nodal_pressure_prescription_t> retrievedTypeOfPrescribedPressures = mp_netlistFileReader->getTypeOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    std::vector<double> retrievedValueOfPrescribedPressures = mp_netlistFileReader->getValueOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);

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
    // node->getIndex() = indexOfNodeInInputData;
    // node->prescribedPressureType = typeOfPrescribedPressure;
    if (node->prescribedPressureType!=Pressure_NotPrescribed && node->prescribedPressureType!=Pressure_Null)
    {
        node->setPressure(retrievedValueOfPrescribedPressures.at(indexOfPrescribedPressure));
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

// nextTimestepWrite_start will be updated and returned to caller of this function.
void NetlistCircuit::writePressuresFlowsAndVolumes(int& nextTimestepWrite_start)
{
	// All the following writes can use the same nextTimestepWrite_end to determine how far to go when looking for data to write to the file:
      int nextTimestepWrite_end = mp_circuitData->getLengthOfHistoryData();
      
      {
        // Write the netlistPressures_zeroDDomainReplacement.dat
        basicFileWriter boundaryConditionPressureHistoryWriter;
        boundaryConditionPressureHistoryWriter.setFileName(m_PressureHistoryFileName);

        for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
        {
          boundaryConditionPressureHistoryWriter.writeStepIndex(stepToWrite);
          for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
          {
            boundaryConditionPressureHistoryWriter.writeToFile(node->second->m_entirePressureHistory.at(stepToWrite));
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
            boundaryConditionFlowHistoryWriter.writeToFile((*component)->m_entireFlowHistory.at(stepToWrite));
          }
          boundaryConditionFlowHistoryWriter.writeEndLine();
        }
      }


      {
        // Write the volumes of the pressure chambers, as netlistVolumes_surface_X.dat
        basicFileWriter boundaryConditionVolumeHistoryWriter;
        boundaryConditionVolumeHistoryWriter.setFileName(m_VolumeHistoryFileName);
        for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
        {
          boundaryConditionVolumeHistoryWriter.writeStepIndex(stepToWrite);
          for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
          {
            VolumeTrackingPressureChamber* pressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
            // If this component is actually a volume chamber, so it actually has a volume history we can write to the file:
            if (pressureChamber != NULL)
            {
              boundaryConditionVolumeHistoryWriter.writeToFile(pressureChamber->getVolumeHistoryAtTimestep(stepToWrite));
            }
          }
          boundaryConditionVolumeHistoryWriter.writeEndLine();
        }
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

void NetlistCircuit::cycleToSetHistoryPressuresFlowsAndVolumes()
{
    // for (auto node=mp_circuitDataWithoutDiodes.mapOfPressureNodes.begin(); node!=mp_circuitDataWithoutDiodes.mapOfPressureNodes.end(); node++)
    // Cycle and store the history pressures
    for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
    {
        // Store the pressure for writing to output file:
        node->second->m_entirePressureHistory.push_back(node->second->getPressure());

        if (node->second->hasHistoryPressure)
        {
            node->second->historyPressure = node->second->getPressure();
        }
    }

    // for (auto component=mp_circuitDataWithoutDiodes.mapOfComponents.begin(); component!=mp_circuitDataWithoutDiodes.mapOfComponents.end(); component++)
    // Cycle and store the history flows
    for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    {
        // Store the flow for writing to output file:
        component->second->m_entireFlowHistory.push_back(component->second->flow);

        if (component->second->hasHistoryFlow)
        {
            component->second->historyFlow = component->second->flow;
        }
    }

    // Store the volumes (currently just for VolumeTrackingPressureChambers. Make this more generic if new volume-tracking components are added later
    // - recommed using the hasHistoryVolume bool).
    for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    {
        VolumeTrackingPressureChamber* pressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
        // Ensure this actually is a VolumeTrackingPressureChamber before going further:
        if (pressureChamber != NULL)
        {
            // Store the volume for writing to output file:
            pressureChamber->recordVolumeInHistory();

            // Make the current volume into the new history volume:
            pressureChamber->cycleHistoryVolume();
        }
    }
}

void NetlistCircuit::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    rebuildCircuitMetadata();
    detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    cycleToSetHistoryPressuresFlowsAndVolumes();
}

void NetlistCircuit::finalizeLPNAtEndOfTimestep()
{
    switchDiodeStatesIfNecessary();
}

boost::shared_ptr<CircuitComponent> NetlistCircuit::getComponentByInputDataIndex(const int componentIndex)
{
    return mp_circuitData->getComponentByInputDataIndex(componentIndex);
}

void NetlistZeroDDomainCircuit::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    rebuildCircuitMetadata();
    cycleToSetHistoryPressuresFlowsAndVolumes();
}

void NetlistZeroDDomainCircuit::createCircuitDescription()
{
    // This function creates the internal CircuitData class format for a zero-D
    // Netlist, non-boundary-condition replacement for the 3D domain, for pure zero-D simulation.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    mp_netlistFileReader = NetlistReader::Instance();

    // Make the appropriate class to store the 3D domain replacement circuit data:
    mp_circuitData = boost::shared_ptr<Netlist3DDomainReplacementCircuitData> (new Netlist3DDomainReplacementCircuitData(m_hstep,m_numberOfNetlistsUsedAsBoundaryConditions));

    // we'll have a resistor for each netlist used as a boundary condition, plus one VolumeTrackingPressureChamber:
    mp_circuitData->numberOfComponents = 2*m_numberOfNetlistsUsedAsBoundaryConditions + 1;

    // The components will be arranged in a star (all start-nodes are connected together at a single point, all end-nodes only connect to one component)
    mp_circuitData->numberOfPressureNodes = 2*m_numberOfNetlistsUsedAsBoundaryConditions + 2;

    // The pressures will initially be received from the boundary conditions at the boundary condition interface nodes,
    // plus the zero-pressure prescription on the base of the VolumeTrackingPressureChamber:
    mp_circuitData->numberOfPrescribedPressures = 1;
    // None initially:
    mp_circuitData->numberOfPrescribedFlows = m_numberOfNetlistsUsedAsBoundaryConditions;

    // Create the component data:
    std::vector<circuit_component_t> componentTypes;
    for (int boundary=0; boundary<m_numberOfNetlistsUsedAsBoundaryConditions; boundary++)
    {
    	componentTypes.push_back(Component_Resistor);
        componentTypes.push_back(Component_Resistor);
    }
    // componentTypes.push_back(Component_VolumeTrackingPressureChamber);
    componentTypes.push_back(Component_Capacitor);

    // Prepare space for the components in the circuit:
    assert(mp_circuitData->components.empty());
    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
        CircuitComponent* toPushBack;
        if (componentTypes.at(ii) == Component_VolumeTrackingPressureChamber)
        {
            toPushBack = new VolumeTrackingPressureChamber(m_hstep,m_thisIsARestartedSimulation);
        }
        else
        {
            toPushBack = new CircuitComponent(m_hstep,m_thisIsARestartedSimulation);
        }

        mp_circuitData->components.push_back(boost::shared_ptr<CircuitComponent> (toPushBack));
        mp_circuitData->components.back()->setIndex(toOneIndexing(ii));
    }

    // Obtain the component- and node-level data for the circuit, for moving into the appropriate data structure CircuitData
    // We want to pop off the component types as we use them, but starting from the beginning of the vector. To do this, we reverse
    // the vector and then pop from the new end.
    std::reverse(componentTypes.begin(), componentTypes.end());
    // Create the component start node info (all are the same: we make it the highest-indexed node.)
    
    // int centralNodeIndex = mp_circuitData->numberOfComponents + 1;
    // std::vector<int> componentStartNodes(mp_circuitData->numberOfComponents, centralNodeIndex);

    // Deal out indices for the end nodes:
    std::vector<int> componentEndNodes;
    // First, do the end nodes for the dP/dQ resistors: those which are "past" the boundary, and have
    // resistance coming from the "implicit coefficient" dP/dQ passed by the
    // boundary condition to this zero-D domain. These are used to mirror the coupling
    // method that would be used if the domain were 3D.
    //
    // These end-nodes themselves will receive the second implicit coefficient as the 
    // pressure shift, on each time-step.
    //
    // (the domain is designed so that the
    // dP/dQ resistors happen to have the same indices as their end-nodes)
    for (int dpDqResistorIndex = 1; dpDqResistorIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; dpDqResistorIndex++)
    {
        componentEndNodes.push_back(dpDqResistorIndex);
    }
    // The next two node indices are used for the domain central node (it's basically star-shaped),
    // and the end-node singleton for the bottom of the compliance chamber.
    //
    // Skip these two, and now do the end nodes for the resistors which represent
    // the resistance of the 3D domain itself.
    for (int internalDomainResistorEndNodeIndex = m_numberOfNetlistsUsedAsBoundaryConditions + 3; internalDomainResistorEndNodeIndex <= mp_circuitData->numberOfPressureNodes; internalDomainResistorEndNodeIndex++)
    {
        componentEndNodes.push_back(internalDomainResistorEndNodeIndex);
    }
    // Finally, do the end-node at the base of the compliance chamber:
    componentEndNodes.push_back(m_numberOfNetlistsUsedAsBoundaryConditions + 1);


    std::vector<int> componentStartNodes;
    // All the components, except for the dP/dQ resistors, have the same
    // start-node (the centre-point of the star-shape of this domain).
    // We do the dP/dQ resistors first:
    for (int dpDqResistorStartNodeIndex = m_numberOfNetlistsUsedAsBoundaryConditions+3; dpDqResistorStartNodeIndex <= mp_circuitData->numberOfPressureNodes; dpDqResistorStartNodeIndex++)
    {
        componentStartNodes.push_back(dpDqResistorStartNodeIndex);
    }
    // All the other components share the same start node, at the centre of the star-shaped domain:
    for (int componentIndex = m_numberOfNetlistsUsedAsBoundaryConditions+1; componentIndex <= mp_circuitData->numberOfComponents; componentIndex++)
    {
        componentStartNodes.push_back(m_numberOfNetlistsUsedAsBoundaryConditions+2);
    }

    std::reverse(componentStartNodes.begin(), componentStartNodes.end()); // actually no point in this call, but it's tidier to leave it here for symmetry with the related calls below.

    std::reverse(componentEndNodes.begin(), componentEndNodes.end());

    // std::vector<double> componentParameterValues(m_numberOfNetlistsUsedAsBoundaryConditions, m_oneResistanceToGiveEachResistor);
    // componentParameterValues.push_back(m_elastanceToGiveVolumeTrackingPressureChamber);

    boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_circuitData);
    for (int dpDqResistorIndex = 1; dpDqResistorIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; dpDqResistorIndex++)
    {
        downcastDomainReplacementCircuitData->addToMapOfDpDqResistors(dpDqResistorIndex,downcastDomainReplacementCircuitData->components.at(toZeroIndexing(dpDqResistorIndex)));    
    }
    // downcastDomainReplacementCircuitData->addToMapOfDpDqResistors(1,downcastDomainReplacementCircuitData->components.at(0));
    // downcastDomainReplacementCircuitData->addToMapOfDpDqResistors(2,downcastDomainReplacementCircuitData->components.at(1));
    // downcastDomainReplacementCircuitData->addToMapOfDpDqResistors(3,downcastDomainReplacementCircuitData->components.at(2));
    
    std::vector<double> componentParameterValues;
    double initialResistanceForDpDqResistors = 0.0;
    for (int dpDqResistorIndex = 1; dpDqResistorIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; dpDqResistorIndex++)
    {
        componentParameterValues.push_back(initialResistanceForDpDqResistors);
    }
    for (int internalDomainResistorIndex = 1; internalDomainResistorIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; internalDomainResistorIndex++)
    {
        componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    }

    // componentParameterValues.push_back(initialResistanceForDpDqResistors);
    // componentParameterValues.push_back(initialResistanceForDpDqResistors);
    // componentParameterValues.push_back(initialResistanceForDpDqResistors);
    // componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    // componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    // componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    componentParameterValues.push_back(m_elastanceToGiveVolumeTrackingPressureChamber);

    std::reverse(componentParameterValues.begin(), componentParameterValues.end());

    

    std::vector<circuit_component_flow_prescription_t> typeOfPrescribedFlows;
    std::vector<double> valueOfPrescribedFlows;
    std::vector<int> listOfPrescribedFlows;
    for (int interfaceComponentIndex=0; interfaceComponentIndex < m_numberOfNetlistsUsedAsBoundaryConditions; interfaceComponentIndex++)
    {
        typeOfPrescribedFlows.push_back(Flow_3DInterface);
        valueOfPrescribedFlows.push_back(0.0);
        listOfPrescribedFlows.push_back(toOneIndexing(interfaceComponentIndex));
    }
    
    std::map<int,double> initialPressures;
    for (int node=1; node <= m_numberOfNetlistsUsedAsBoundaryConditions; node++)
    {
    	initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }
    // Do the node at the base of the pressure chamber:
    initialPressures.insert(std::make_pair(m_numberOfNetlistsUsedAsBoundaryConditions+1, 0.0));
    for (int node=m_numberOfNetlistsUsedAsBoundaryConditions+2; node <= mp_circuitData->numberOfPressureNodes; node++)
    {
        initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }

    // Do the final node (at the base of the pressure chamber):
    // initialPressures.insert(std::make_pair(mp_circuitData->numberOfPressureNodes, 0.0));
    
    // Loop over the components, assigning them (and their nodes) the appropriate properties to give the fully-described circuit:
    for (auto component = mp_circuitData->components.begin(); component != mp_circuitData->components.end(); component++)
    {
        (*component)->getType() = componentTypes.back();
        componentTypes.pop_back();

        int indexOfStartNodeInInputData = componentStartNodes.back();
        componentStartNodes.pop_back();
        
        // note that we're passing a boost::shared_ptr to the startNode here; if this is currently NULL, it will
        // be constructed during the call to setupPressureNode.
        setupPressureNode(indexOfStartNodeInInputData, (*component)->startNode, *component);

        int indexOfEndNodeInInputData = componentEndNodes.back();
        componentEndNodes.pop_back();

        // note that we're passing a boost::shared_ptr to the endNode here; if this is currently NULL, it will
        // be constructed during the call to setupPressureNode.
        setupPressureNode(indexOfEndNodeInInputData, (*component)->endNode, *component);

        (*component)->prescribedFlowType = Flow_NotPrescribed;  // initialise as a default, before replacing as necessary
        // for (int prescribedFlow=0; prescribedFlow<mp_circuitData->numberOfPrescribedFlows; prescribedFlow++)
        // {
        //     if (listOfPrescribedFlows.at(prescribedFlow) == (*component)->getIndex())
        //     {
        //         (*component)->prescribedFlowType = typeOfPrescribedFlows.at(prescribedFlow);
        //         (*component)->valueOfPrescribedFlow = valueOfPrescribedFlows.at(prescribedFlow);
        //     }
        // }

        (*component)->setParameterValue(componentParameterValues.back());
        // make a copy of this value so we can reset it if necessary. Used for e.g. diode state changes.
        (*component)->parameterValueFromInputData = *((*component)->getParameterPointer());
        componentParameterValues.pop_back();

        (*component)->startNode->setPressure(initialPressures.at((*component)->startNode->getIndex()));
        (*component)->endNode->setPressure(initialPressures.at((*component)->endNode->getIndex()));

        if ((*component)->getIndex() <= m_numberOfNetlistsUsedAsBoundaryConditions)
        {
            (*component)->prescribedFlowPointerIndex = (*component)->getIndex()-1;
        }


    }

    // Some metadata is already set-up for this circuit, as a side-effect
    // of the above construction. This call completes the metadata; it's not
    // a problem that it also re-writes some of the existing metadata
    // (rewrites - but does not change - the values are identical!)
    mp_circuitData->rebuildCircuitMetadata();

    // Tell the node at the 3D interface that it connects to the 3D domain:
    {
        std::vector<int> threeDNodeIndices;
        for (int nodeIndex=1; nodeIndex<=m_numberOfNetlistsUsedAsBoundaryConditions; nodeIndex++)
        {
        	threeDNodeIndices.push_back(nodeIndex);
        }
        boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_circuitData);
        downcastDomainReplacementCircuitData->initialiseNodesAndComponentsAtInterface_vector(threeDNodeIndices);
    }

    // mp_circuitData->switchDiodeStatesIfNecessary();
    // mp_circuitData->detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    // // Component indices are just consecutive integers by default, but sometimes non-consecutive numbering
    // // is needed; componentIndices allows for this.
    // // We initialise it now for the default case.
    // for (int ii=1; ii < mp_circuitData->numberOfComponents + 1; ii++)
    // {
    //     mp_circuitData->componentIndices.push_back(ii);
    // }
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
void NetlistZeroDDomainCircuit::setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode)
{
    std::vector<int> listOfPrescribedPressures;
    for (int pressureNodeIndex=0; pressureNodeIndex < m_numberOfNetlistsUsedAsBoundaryConditions+1; pressureNodeIndex++)
    {
    	listOfPrescribedPressures.push_back(pressureNodeIndex+1);
    }

    // All but one of the prescribed pressures are at the boundary interfaces:
    std::vector<circuit_nodal_pressure_prescription_t> typeOfPrescribedPressures;
    for (int pressureNodeIndex=0; pressureNodeIndex < m_numberOfNetlistsUsedAsBoundaryConditions; pressureNodeIndex++)
    {
    	typeOfPrescribedPressures.push_back(Pressure_3DInterface);
    }
    // The last prescribed pressure is the one at the base of the VolumeTrackingPressureChamber
    typeOfPrescribedPressures.push_back(Pressure_Fixed);

    std::vector<double> valueOfPrescribedPressures;
    // All but one of the prescribed pressures are at the boundary interfaces.
    for (int pressureNodeIndex=0; pressureNodeIndex < m_numberOfNetlistsUsedAsBoundaryConditions; pressureNodeIndex++)
    {
        //\todo remove \hardcoded to millimetres! units: 10^0 x Pa here
    	valueOfPrescribedPressures.push_back(1332.0); // ~10 mmHg
    }
    // The last prescribed pressure is the one at the base of the VolumeTrackingPressureChamber
    valueOfPrescribedPressures.push_back(0.0);


    // Discover whether this node has a prescribed pressure, and if so, what type:
    circuit_nodal_pressure_prescription_t typeOfPrescribedPressure = Pressure_NotPrescribed; // initialise, but chnage later if pressure is actually prescribed
    int indexOfPrescribedPressure = -1; // initialise to a nonsense value to detect errors
    for (int prescribedPressure=0; prescribedPressure<listOfPrescribedPressures.size(); prescribedPressure++)
    {
        if (listOfPrescribedPressures.at(prescribedPressure) == indexOfNodeInInputData)
        {
            typeOfPrescribedPressure = typeOfPrescribedPressures.at(prescribedPressure);
            indexOfPrescribedPressure = prescribedPressure;
        }
    }
    // Get the node (or create a new node if this one hasn't been made yet)
    node = mp_circuitData->ifExistsGetNodeOtherwiseConstructNode(indexOfNodeInInputData,typeOfPrescribedPressure,componentNeighbouringThisNode);    
    // node->getIndex() = indexOfNodeInInputData;
    // node->prescribedPressureType = typeOfPrescribedPressure;
    if (node->prescribedPressureType!=Pressure_NotPrescribed && node->prescribedPressureType!=Pressure_Null)
    {
        node->setPressure(valueOfPrescribedPressures.at(indexOfPrescribedPressure));
    }

    if (node->getIndex() <= m_numberOfNetlistsUsedAsBoundaryConditions)
    {
        node->prescribedPressurePointerIndex = node->getIndex()-1;
    }
}

boost::shared_ptr<CircuitData> NetlistCircuit::getCircuitDescription()
{
    return mp_circuitData;
}


void NetlistZeroDDomainCircuit::setBoundaryPrescriptionsAndBoundaryConditionTypes(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate)
{
    boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_circuitData);
    assert(downcastDomainReplacementCircuitData!=NULL);
	downcastDomainReplacementCircuitData->setBoundaryPrescriptionsAndBoundaryConditionTypes(boundaryFlowsOrPressuresAsAppropriate);
}

std::vector<double> NetlistZeroDDomainCircuit::getBoundaryPressures()
{
	std::vector<double> pressures;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		pressures.push_back(mp_circuitData->components.at(indexOfBoundaryInterfaceComponent)->startNode->getPressure());
	}
	return pressures;
}

std::vector<double> NetlistZeroDDomainCircuit::getBoundaryFlows()
{
	std::vector<double> flows;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		flows.push_back(mp_circuitData->components.at(indexOfBoundaryInterfaceComponent)->flow);
	}
	return flows;
}

void NetlistZeroDDomainCircuit::solveSystem(const int timestepNumber)
{
	buildAndSolveLinearSystem(timestepNumber,m_delt);
}

void NetlistZeroDDomainCircuit::setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients, std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAtBoundaries)
{
    int dpDqResistorIndex=0;
    for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
    {
        boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_circuitData);
        if (downcastDomainReplacementCircuitData->isADpDqResistor((*component)->getIndex())) // these are the extra resistors added for the zeroD domain replacement. Needs making generic (set for 3 outlets here!)
        {
            bool boundaryCurrentlyUsesDpDqResistance = (pressuresOrFlowsAtBoundaries.at(dpDqResistorIndex).first == Boundary_Pressure);
            if (boundaryCurrentlyUsesDpDqResistance)
            {
                // the component getIndex() needs to be converted to zero-indexing:
                double potentialResistance = allImplicitCoefficients.at((*component)->getIndex()-1).first;
                assert(!isnan(potentialResistance));

                // if (potentialResistance == 0.0)
                // {
                //     potentialResistance = 0.000001;
                // }
                (*component)->setParameterValue(potentialResistance);
            }
            dpDqResistorIndex++;
        }
    }
}

void NetlistCircuit::initialiseCircuit()
{
    // This function exists just so we can modify what initialiseCircuit does in subclasses without repeating code.
    initialiseCircuit_common();
    // The system is square in this case
    m_numberOfSystemRows = m_numberOfSystemColumns;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::initialiseCircuit()
{
    // Discover which pressure nodes don't need Kirchoff law applications yet, because such laws will be applied later
    // once this circuit is combined with the downstreamCircuit to make a (closed loop)-type boundary circuit.
    m_pressureNodesWhichConnectToDownstreamCircuits = NetlistDownstreamCircuitReader::Instance()->getSetOfNodesInBoundaryConditionWhichConnectToDownstreamCircuit(m_surfaceIndex);
    m_numberOfNodesConnectingToAnotherCircuit = m_pressureNodesWhichConnectToDownstreamCircuits.size();

    // This function exists just so we can modify what initialiseCircuit does in subclasses without repeating code.
    initialiseCircuit_common();

    m_numberOfSystemRows = m_numberOfSystemColumns - m_numberOfNodesConnectingToAnotherCircuit;
}

void NetlistClosedLoopDownstreamCircuit::initialiseCircuit()
{
    // Discover which pressure nodes don't need Kirchoff law applications yet, because such laws will be applied later
    // once this circuit is combined with the upstream boundary conditions to make a (closed loop)-type boundary circuit.
    std::vector<int> nodesConnectingToBoundaryCircuits = NetlistDownstreamCircuitReader::Instance()->getLocalBoundaryConditionInterfaceNodes(m_downstreamCircuitIndex);
    // Convert the vector to a set, which is more convenient for checking membership:
    for (auto node = nodesConnectingToBoundaryCircuits.begin(); node != nodesConnectingToBoundaryCircuits.end(); node++)
    {
        m_pressureNodesWhichConnectToBoundaryCircuits.insert(*node);
    }
    m_numberOfNodesConnectingToAnotherCircuit = m_pressureNodesWhichConnectToBoundaryCircuits.size();

    // This function exists just so we can modify what initialiseCircuit does in subclasses without repeating code.
    initialiseCircuit_common();

    m_numberOfSystemRows = m_numberOfSystemColumns - m_numberOfNodesConnectingToAnotherCircuit;
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

  createVectorsAndMatricesForCircuitLinearSystem();

  // columnMapSize = m_numberOfHistoryPressures + numberOfHistoryFlows + numberOfPrescribedPressures + numberOfPrescribedFlows;
  createListOfNodesWithMultipleIncidentCurrents();
}

int NetlistCircuit::getNumberOfHistoryPressures()
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
  errFlag = VecSetSizes(m_solutionVector,m_numberOfSystemRows,sym_numberOfSystemRowsstemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
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
       if (numberOfTimesNodeSeen > 1 && kirchoffEquationAtNodeNotDeferredToInterfacingCircuit(nodeIndex))
       {
          listOfNodesWithMultipleIncidentCurrents.push_back(nodeIndex);
          m_numberOfMultipleIncidentCurrentNodes++;
       }
    }

}

bool NetlistCircuit::kirchoffEquationAtNodeNotDeferredToInterfacingCircuit(const int nodeIndex) const
{
    // In NetlistCircuit, there is no downstream circuit, so the return value is always true.
    return true;
}

bool NetlistBoundaryCircuitWhenDownstreamCircuitsExist::kirchoffEquationAtNodeNotDeferredToInterfacingCircuit(const int nodeIndex) const
{
    bool nodeInterfacesWithDownstreamCircuit = false;
    if (m_pressureNodesWhichConnectToDownstreamCircuits.count(nodeIndex) == 1)
    {
        nodeInterfacesWithDownstreamCircuit = true;
    }
    return nodeInterfacesWithDownstreamCircuit;
}

bool NetlistClosedLoopDownstreamCircuit::kirchoffEquationAtNodeNotDeferredToInterfacingCircuit(const int nodeIndex) const
{
    bool nodeInterfacesWithDownstreamCircuit = false;
    if (m_pressureNodesWhichConnectToBoundaryCircuits.count(nodeIndex) == 1)
    {
        nodeInterfacesWithDownstreamCircuit = true;
    }
    return nodeInterfacesWithDownstreamCircuit;
}

void NetlistCircuit::getMapOfPressHistoriesToCorrectPressNodes()
{
    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
       // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
       if (mp_circuitData->components.at(ii)->getType() == Component_Capacitor)
       {
            listOfHistoryPressures.insert(mp_circuitData->components.at(ii)->startNode->getIndex());
          mp_circuitData->components.at(ii)->startNode->hasHistoryPressure = true;
            listOfHistoryPressures.insert(mp_circuitData->components.at(ii)->endNode->getIndex());
          mp_circuitData->components.at(ii)->endNode->hasHistoryPressure = true;
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
     if(mp_circuitData->components.at(ii)->getType() == Component_VolumeTrackingPressureChamber)
     {
        listOfHistoryVolumes.insert(ii);
        mp_circuitData->components.at(ii)->hasHistoryVolume = true;
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
     if(mp_circuitData->components.at(ii)->getType() == Component_VolumeTrackingPressureChamber)
     {
        listOfTrackedVolumes.insert(ii);
        mp_circuitData->components.at(ii)->hasTrackedVolume = true;
     }
  }

  m_numberOfTrackedVolumes = listOfTrackedVolumes.size();

  // Now do the actual generation of the volume node ordering map (for use when bulding the linear system matrix):
  int ii=0;
  for (auto iterator=listOfTrackedVolumes.begin(); iterator != listOfTrackedVolumes.end(); iterator++)
  {
     componentIndexToTrackedVolumeComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
  }
}

void NetlistCircuit::generateLinearSystemFromPrescribedCircuit(const double alfi_delt)
{
    // Build m_systemMatrix
    generateLinearSystemWithoutFactorisation(alfi_delt);
    
    // LU factor m_systemMatrix
    PetscErrorCode errFlag = MatLUFactor(m_systemMatrix,NULL,NULL,NULL);CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistCircuit::generateLinearSystemWithoutFactorisation(const double alfi_delt)
{
    // This function assembles the system of (time-discretised) linear algebraic equations for the LPN.
    PetscErrorCode errFlag;

    errFlag = MatZeroEntries(m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    {
      int row = 0; // is the row in the matrix that we write to on each occasion
      for (auto component=mp_circuitData->components.begin(); component!=mp_circuitData->components.end(); component++)
      {
          // bool componentIsOpenDiode = (mp_circuitData->components.at(ll)->type == Component_Diode &&
          //                              mp_circuitData->components.at(ll)->hasNonnegativePressureGradientAndNoBackflow());
          // open diodes are just implemented as zero-resistance resistors, closed diodes are zero-resistance resistors with prescribed zero flow
          if ((*component)->getType() == Component_Resistor || (*component)->getType() == Component_Diode)
          {
            // insert resistor(-eqsue) relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+mp_circuitData->numberOfPressureNodes+m_numberOfHistoryPressures,-currentParameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
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
              // std::cout << "setting in NetlistSubcircuit.cxx row and column: " << row << " " << columnIndex << std::endl;
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
              // Reset the zero-volume marker on the component:
              volumeTrackingPressureChamber->resetZeroVolumePrescription();
            }
            row++; // done twice in this if-case, because there are 2 equations to create for the VolumeTrackingPressureChamber

          }
          else
          {
            throw std::runtime_error("EE: Unknown component type in netlist. Halting.");
          }
      }
    }

     // Do the equations for the nodes with multiple incident currents:
     for (int mm = 0; mm < m_numberOfMultipleIncidentCurrentNodes; mm++)
     {
       for (int ll=0; ll<mp_circuitData->numberOfComponents; ll++)
       {
          bool foundMultipleIncidentCurrentsForEndNode = (mp_circuitData->components.at(ll)->endNode->getIndex() == listOfNodesWithMultipleIncidentCurrents.at(mm)); 
          if (foundMultipleIncidentCurrentsForEndNode)
          {
            int row = mm + mp_circuitData->numberOfComponents + m_numberOfTrackedVolumes;
            int column = ll + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures;
            errFlag = MatSetValue(m_systemMatrix,row,column,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }

          bool foundMultipleIncidentCurrentsForStartNode = (mp_circuitData->components.at(ll)->startNode->getIndex() == listOfNodesWithMultipleIncidentCurrents.at(mm));
          if (foundMultipleIncidentCurrentsForStartNode)
          {
            int row = mm + mp_circuitData->numberOfComponents + m_numberOfTrackedVolumes;
            int column = ll + mp_circuitData->numberOfPressureNodes + m_numberOfHistoryPressures;
            errFlag = MatSetValue(m_systemMatrix,row,column,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
       }
     }

     int rowsDoneSoFar = mp_circuitData->numberOfComponents + m_numberOfMultipleIncidentCurrentNodes + m_numberOfTrackedVolumes;

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

     // Set the prescribed-value equations (i.e. pressure_1 (LHS) = pressure_1 (m_RHS) - so really just a way of setting the prescribed values within the linear system)
     for (int ll = 0; ll < m_numberOfSystemRows - rowsDoneSoFar; ll++)
     {
       errFlag = MatSetValue(m_systemMatrix,rowsDoneSoFar + ll, columnMap.at(ll),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     }

     errFlag = MatAssemblyBegin(m_systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     errFlag = MatAssemblyEnd(m_systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    // std::cout << "System matrix for surface " << m_surfaceIndex << ":" << std::endl;
    //  errFlag = MatView(m_systemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistCircuit::assembleRHS(const int timestepNumber)
{

    PetscErrorCode errFlag;
    errFlag = VecZeroEntries(m_RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    historyPressuresInSubcircuit = pressuresInSubcircuit;
    historyFlowsInSubcircuit = flowsInSubcircuit;
    historyVolumesInSubcircuit = volumesInSubcircuit;

    columnIndexOf3DInterfacePressureInLinearSystem.clear(); // dummy value, to be replaced!

    // int nextPressurePointerIndex = 0; // for tracking which pressure pointer to use - useful when there are multiple pressure interfaces to other domains / subcircuits

    // Prescribed pressures
    int tempIndexingShift = mp_circuitData->numberOfComponents + m_numberOfMultipleIncidentCurrentNodes + m_numberOfTrackedVolumes;
    {
      int ll=0;
      for (auto prescribedPressureNode=mp_circuitData->mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=mp_circuitData->mapOfPrescribedPressureNodes.end(); prescribedPressureNode++ )
      {
        // Coming from 'f' for 'fixed' in the input data:
        if (prescribedPressureNode->second->prescribedPressureType == Pressure_Fixed)
        {
            errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,prescribedPressureNode->second->getPressure(),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);    
        }
        // Coming from 'l' for 'left-ventricular' in the input data:
        else if (prescribedPressureNode->second->prescribedPressureType == Pressure_LeftVentricular)
        {
            std::cerr << "this requires heart model. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
            std::exit(1);
        }
        else if (prescribedPressureNode->second->prescribedPressureType == Pressure_3DInterface)
        {
            // We only do this if the netlist is in Dirichlet BC mode:
            if (mp_circuitData->hasPrescribedPressureAcrossInterface())
            {
              columnIndexOf3DInterfacePressureInLinearSystem.push_back(ll + tempIndexingShift);
              double* pressurePointerToSet = pressure_n_ptrs.at(prescribedPressureNode->second->prescribedPressurePointerIndex);
              errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,*pressurePointerToSet,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // nextPressurePointerIndex++;
            }
        }
        else
        {
              throw std::runtime_error("Unknown pressure prescription value in Netlist.");
        }

        // // get the column index for the pressure in the linear system:
        // if (prescribedPressureNode->second->isAtBoundary())
        // {
        //   columnIndexOf3DInterfacePressureInLinearSystem.push_back(ll + tempIndexingShift);
        // }

        ll++;
      }
    }

    // for (int ll=0; ll<mp_circuitData->numberOfPrescribedPressures; ll++)
    // {
    //    // 'f' for 'fixed'
    //    if (mp_circuitData->typeOfPrescribedPressures.at(ll) == Pressure_Fixed)
    //    {
    //       errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,valueOfPrescribedPressures.at(ll),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //    }
    //    // 'l' for 'leftVentricular'
    //    else if (subcircuitInputData.typeOfPrescribedPressures.at(ll) == Pressure_LeftVentricular)
    //    {
    //       std::cout << "this requires heartmodel. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
    //       std::exit(1);
    //       // if ((timestepNumber == 0) || (timestepNumber == 1)) // treat case with no known IM pressure yet
    //       // {
    //       //    P_IM_mid_lasttimestep = 5000; // \todo find a better way of doing this; maybe input this value from file...
    //       //    P_IM_mid = 5000; // ... or set it based on the aortic valve state at simulation start
    //       // }
    //       // elseif (timestepNumber .eq. int(2)) then // treat case where only one IM pressure history point is known
    //       //    P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber-1)
    //       //    P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber)
    //       // else // get the previous intramyocardial pressure in the case where we have enough doata for this (see comment before start of "if" block)
    //       //    P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber-1)
    //       //    P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber)
    //       // end if

    //       // errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,P_IM_mid,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //       // int nn=0;
    //       // for (auto iterator=listOfHistoryPressures.begin(); iterator!=listOfHistoryPressures.end(); iterator++, nn++)
    //       // {
    //       //    if (*iterator == listOfPrescribedPressures.at(ll))
    //       //    {
    //       //       historyPressuresInSubcircuit.at(*iterator) = P_IM_mid_lasttimestep;
    //       //    }
    //       // }
    //    }
    //    else
    //    {
    //         throw std::runtime_error("Unknown pressure prescription value in Netlist.");
    //    }
    // }
    // History Pressures
    tempIndexingShift = tempIndexingShift + mp_circuitData->numberOfPrescribedPressures;
    // for(int ll=0; ll<m_numberOfHistoryPressures; ll++)

    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
        int lll=0;
        for (auto node=mp_circuitData->mapOfPressureNodes.begin(); node!=mp_circuitData->mapOfPressureNodes.end(); node++)
        {
            if (node->second->hasHistoryPressure)
            {
                errFlag = VecSetValue(m_RHS,lll+tempIndexingShift,node->second->historyPressure,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
                lll++;
            }
        }
        // for (auto iterator=listOfHistoryPressures.begin(); iterator!=listOfHistoryPressures.end(); iterator++)
        // {
        //     errFlag = VecSetValue(m_RHS,lll+tempIndexingShift,historyPressuresInSubcircuit.at(*iterator - 1),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        // }
        // lll++;
      }
    // Prescribed Flows:
    tempIndexingShift += m_numberOfHistoryPressures;
    columnIndexOf3DInterfaceFlowInLinearSystem.clear();
    // int nextFlowPointerIndex = 0;
    // Scoping unit to include the second counter ll in the for loop, without having ll in-scope after the loop finishes:
    {
        int ll=0;
        for (auto prescribedFlowComponent=mp_circuitData->mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=mp_circuitData->mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
        {
           if (prescribedFlowComponent->second->prescribedFlowType == Flow_3DInterface)
           {
            if (mp_circuitData->hasPrescribedFlowAcrossInterface())
            {
              columnIndexOf3DInterfaceFlowInLinearSystem.push_back(ll + tempIndexingShift);
              double* flowPointerToSet = flow_n_ptrs.at(prescribedFlowComponent->second->prescribedFlowPointerIndex);
              // First, flip the sign of the flow, if necessary due to the orientation of the component at the 3D interface:
              double threeDFlowValue = *flowPointerToSet * prescribedFlowComponent->second->m_signForPrescribed3DInterfaceFlow;
              assert(!isnan(threeDFlowValue));
              // Give the (possibly sign-corrected) flow to the linear system:
              errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,threeDFlowValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // nextFlowPointerIndex++;
            }
           }
           else if (prescribedFlowComponent->second->prescribedFlowType == Flow_Fixed)
           {
              errFlag = VecSetValue(m_RHS,ll + tempIndexingShift,prescribedFlowComponent->second->valueOfPrescribedFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
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
        if ((*component)->hasHistoryVolume)
        {
          // currently, only VolumeTrackingPresureChambers have history volumes. We might want to change this cast later, if new component types
          // with history volumes get added.
          boost::shared_ptr<VolumeTrackingPressureChamber> volumeTrackingPressureChamber = boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (*component);
          double volume = volumeTrackingPressureChamber->getHistoryVolume();
          int row = ll + tempIndexingShift;
          errFlag = VecSetValue(m_RHS,row,volume,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          ll++;
        }
      }
    }
 
    errFlag = VecAssemblyBegin(m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(m_RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // std::cout << "m_RHS for surface " << m_surfaceIndex << ":" << std::endl;
    // errFlag = VecView(m_RHS,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // std::cout << "END m_RHS for surface " << m_surfaceIndex << std::endl;

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

    // for (int volumeIndex = indexShift; volumeIndex < indexShift + m_numberOfTrackedVolumes; volumeIndex++)
    // {
    //     errFlag = VecGetValues(m_solutionVector,getSingleValue,&volumeIndex,&volumesInSubcircuit[volumeIndex-indexShift]); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    //     VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (mp_circuitData->mapOfVolumeTrackingComponents.at(toOneIndexing(volumeIndex-indexShift)).get());
    //     currentPressureChamber->setStoredVolume(volumesInSubcircuit[volumeIndex-indexShift]);
    //     currentPressureChamber->passPressureToStartNode();
    // }

    // // Update the volumes in each VolumeTrackingPressureChamber
    // for (auto component=mp_circuitData->mapOfComponents.begin(); component!=mp_circuitData->mapOfComponents.end(); component++)
    // {
    //   // detect the VolumeTrackingPressureChambers:
    //   if (component->second->type == Component_VolumeTrackingPressureChamber)
    //   {
    //     VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
    //     currentPressureChamber->updateStoredVolume(alfi_delt);
    //     currentPressureChamber->passPressureToStartNode();
    //   }
    // }

    // write(*,*) 'discrepancy:', (-this%P_a(1) - this%pressuresInSubcircuit(2))/1.2862d5 - this%flowsInSubcircuit(1)
}

void NetlistCircuit::giveNodesTheirPressuresFromSolutionVector()
{
  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;

  // Look the nodes, handing them their new pressures from the circuit linear system solve:
  for (int ll=0; ll<mp_circuitData->numberOfPressureNodes; ll++)
  {
      errFlag = VecGetValues(m_solutionVector,getSingleValue,&ll,&pressuresInSubcircuit[ll]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      // std::cout << "system matrix: " << m_surfaceIndex << std::endl;
      // MatView(m_systemMatrix,PETSC_VIEWER_STDOUT_WORLD);
      // std::cout << "m_RHS: " << std::endl;
      // VecView(m_RHS,PETSC_VIEWER_STDOUT_WORLD);
      // std::cout << "solution vec: " << std::endl;
      // VecView(m_solutionVector,PETSC_VIEWER_STDOUT_WORLD);
      assert(!isnan(pressuresInSubcircuit[ll]));
      mp_circuitData->mapOfPressureNodes.at(toOneIndexing(ll))->setPressure(pressuresInSubcircuit[ll]);
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
      VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
      
      // Ensure we aren't dealing with negative volumes:
      //\todo REINSTATE!
      // assert(volumes.back() >= 0.0);
      assert(!isnan(volumes.back()));

      currentPressureChamber->setStoredVolume(volumes.back());
      volumes.pop_back();
      currentPressureChamber->passPressureToStartNode();
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
      VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
      
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

    VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
    volumesToReturn.push_back(volumesInSubcircuit[volumeIndex-firstVolumeIndex]);

    volumeIndex++;
    component++;
  }

  return volumesToReturn;
}

void NetlistCircuit::buildAndSolveLinearSystem(const int timestepNumber, const double alfi_delt)
{
  generateLinearSystemFromPrescribedCircuit(alfi_delt);
  assembleRHS(timestepNumber);

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

std::pair<double,double> NetlistBoundaryCircuitWhenDownstreamCircuitsExist::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    // Call the downstream circuit(s?) and tell them to contact each NetlistBoundaryCondition to get
    // their contributions to the (closed loop)-type matrix, build the full matrix and solve it
    // 
    // Internally, that call will only do anything if this is the first boundary condition to try to computeImplicitCoefficients
    // this time. Otherwise, the system state has already been computed, and we don't need to do anything
    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
    {
        downstreamSubcircuit->lock()->buildAndSolveLinearSystemIfNotYetDone();
    }

    // Call the downstream circuits to get the implicit coefficients that they computed for this
    // boundary
    std::pair<double,double> returnValue;
    int counterToDetectErrors = 0;
    for (auto downstreamSubcircuit = m_netlistDownstreamLoopClosingSubcircuits.begin(); downstreamSubcircuit != m_netlistDownstreamLoopClosingSubcircuits.end(); downstreamSubcircuit++)
    {
        if (downstreamSubcircuit->lock()->boundaryConditionCircuitConnectsToThisDownstreamSubsection(m_IndexOfThisNetlistLPN))
        {
            returnValue = downstreamSubcircuit->lock()->getImplicitCoefficients(m_IndexOfThisNetlistLPN);
            counterToDetectErrors++;
        }
    }

    assert(counterToDetectErrors == 1);

    return returnValue;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getMatrixContribution(Mat& matrixFromThisBoundary)
{
    generateLinearSystemWithoutFactorisation(alfi_delt);
    matrixFromThisBoundary = m_systemMatrix;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getRHSContribution(Vec& rhsFromThisBoundary)
{
    assembleRHS(timestepNumber);
    rhsFromThisBoundary = m_RHS;
}

std::pair<double,double> NetlistCircuit::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    assert(mp_circuitData->connectsTo3DDomain());

    PetscErrorCode errFlag;
    {
      int safetyCounter = 0;
      bool solutionVectorMightHaveNegativeVolumes = true;
      // This loop keeps repeating the circuit linear system solve, until we have
      // ensured there are no negative volumes:
      while (solutionVectorMightHaveNegativeVolumes)
      {
        buildAndSolveLinearSystem(timestepNumber,alfi_delt);

        solutionVectorMightHaveNegativeVolumes = areThereNegativeVolumes(timestepNumber, alfi_delt);

        safetyCounter++;
        if (safetyCounter > 1 && safetyCounter < 5)
        {
          std::cout << "II: Redoing due to a detected negative-volume problem! ----------------------------------------------" << std::endl;
        }
        if (safetyCounter > safetyCounterLimit)
        {
          std::stringstream errorMessage;
          errorMessage << "EE: Took too long (" << safetyCounter << " repeated solves of the circuit linear system) to eradicate negative volumes at the domain boundary with index " << m_surfaceIndex << "." << std::endl;
          errorMessage << "This was probably caused by a bad (or an extremely large) Netlist circuit!" << std::endl;
          throw std::runtime_error(errorMessage.str());
        }
      }
    }

    // Extract the implicit coeffcients, for eventual passing to the FORTRAN
    // linear solve
    int rowToGet[] = {0};
    int numberOfValuesToGet=1;
    PetscScalar valueFromInverseOfSystemMatrix;
    
    std::pair<double,double> implicitCoefficientsToReturn;

    errFlag = MatGetValues(m_inverseOfSystemMatrix,numberOfValuesToGet,rowToGet,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem.at(0),&valueFromInverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    implicitCoefficientsToReturn.first = valueFromInverseOfSystemMatrix;

    PetscScalar valueFromRHS;
    errFlag = VecGetValues(m_RHS,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem.at(0),&valueFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    PetscScalar valueFromSolutionVector;
    errFlag = VecGetValues(m_solutionVector,numberOfValuesToGet,rowToGet,&valueFromSolutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    implicitCoefficientsToReturn.second = valueFromSolutionVector - valueFromInverseOfSystemMatrix * valueFromRHS;//\todo make dynamic
    
    // std::cout << "m_inverseOfSystemMatrix: "<< std::endl;
    // MatView(m_inverseOfSystemMatrix,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "solution vector: "<< std::endl;
    // VecView(m_solutionVector,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "m_RHS vector: "<< std::endl;
    // VecView(m_RHS,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "and just set " << implicitCoefficientsToReturn.first << " " <<implicitCoefficientsToReturn.second << std::endl;


    return implicitCoefficientsToReturn;
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
  //   VolumeTrackingPressureChamber* volumeTrackingPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
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


void NetlistClosedLoopDownstreamCircuit::createCircuitDescription()
{
    // This function takes the read-in netlist circuit description and converts it
    // to the internal CircuitData class format.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    mp_netlistFileReader = NetlistDownstreamCircuitReader::Instance();
    createBasicCircuitDescription();
    appendClosedLoopSpecificCircuitDescription();
}

void NetlistClosedLoopDownstreamCircuit::appendClosedLoopSpecificCircuitDescription()
{
    m_numberOfConnectedBoundaryConditions = getNumberOfBoundaryConditionsConnectedTo(m_downstreamCircuitIndex);
    m_connectedCircuitSurfaceIndices = getConnectedCircuitSurfaceIndices(m_downstreamCircuitIndex);
    m_localInterfacingNodes = getLocalBoundaryConditionInterfaceNodes(m_downstreamCircuitIndex);
    m_remoteInterfacingNodes = getRemoteBoundaryConditionInterfaceNodes(m_downstreamCircuitIndex);
}

void NetlistClosedLoopDownstreamCircuit::getMatrixContribution(Mat& matrixFromThisBoundary)
{
    generateLinearSystemWithoutFactorisation(alfi_delt);
    matrixFromThisBoundary = m_systemMatrix;
}

void NetlistClosedLoopDownstreamCircuit::getRHSContribution(Vec& rhsFromThisBoundary)
{
    assembleRHS(timestepNumber);
    rhsFromThisBoundary = m_RHS;
}

// Disable unwanted methods:
void NetlistClosedLoopDownstreamCircuit::detectWhetherClosedDiodesStopAllFlowAt3DInterface()
{
    throw std::logic_error("Method detectWhetherClosedDiodesStopAllFlowAt3DInterface() should not be called on the NetlistClosedLoopDownstreamCircuit, as it has no 3D interface.");
}