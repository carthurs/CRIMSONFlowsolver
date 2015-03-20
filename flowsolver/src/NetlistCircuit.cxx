#include "NetlistCircuit.hxx"
#include "fileReaders.hxx"
#include "fileWriters.hxx"
#include <boost/make_shared.hpp>
#include "indexShifters.hxx"

void NetlistCircuit::createCircuitDescription()
{
    // This function takes the read-in netlist circuit description and converts it
    // to the internal CircuitData class format.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    netlistReader* netlistReader_instance = netlistReader::Instance();
    mp_CircuitDescription->numberOfComponents = netlistReader_instance->getNumberOfComponents().at(m_IndexOfThisNetlistLPN);
    mp_CircuitDescription->numberOfPressureNodes = netlistReader_instance->getNumberOfPressureNodes().at(m_IndexOfThisNetlistLPN);
    mp_CircuitDescription->numberOfPrescribedPressures = netlistReader_instance->getNumberOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    mp_CircuitDescription->numberOfPrescribedFlows = netlistReader_instance->getNumberOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);

    std::vector<circuit_component_t> retrievedComponentTypes = netlistReader_instance->getComponentTypes().at(m_IndexOfThisNetlistLPN);

    // Prepare space for the components in the circuit:
    assert(mp_CircuitDescription->components.empty());
    for (int ii=0; ii<mp_CircuitDescription->numberOfComponents; ii++)
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

        mp_CircuitDescription->components.push_back(boost::shared_ptr<CircuitComponent> (toPushBack));
        mp_CircuitDescription->components.back()->setIndex(toOneIndexing(ii));
    }


    // Obtain the component- and node-level data for the circuit, for moving into the appropriate data structure CircuitData
    // We want to pop off the component types as we use them, but starting from the beginning of the vector. To do this, we reverse
    // the vector and then pop from the new end.
    std::reverse(retrievedComponentTypes.begin(), retrievedComponentTypes.end());
    std::vector<int> retrievedComponentStartNodes = netlistReader_instance->getComponentStartNodes().at(m_IndexOfThisNetlistLPN);
    std::reverse(retrievedComponentStartNodes.begin(), retrievedComponentStartNodes.end());
    std::vector<int> retrievedComponentEndNodes = netlistReader_instance->getComponentEndNodes().at(m_IndexOfThisNetlistLPN);
    std::reverse(retrievedComponentEndNodes.begin(), retrievedComponentEndNodes.end());
    std::vector<double> retrievedComponentParameterValues = netlistReader_instance->getComponentParameterValues().at(m_IndexOfThisNetlistLPN);
    std::reverse(retrievedComponentParameterValues.begin(), retrievedComponentParameterValues.end());

    std::vector<int> retrievedListOfPrescribedFlows = netlistReader_instance->getListOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::vector<circuit_component_flow_prescription_t> retrievedTypeOfPrescribedFlows = netlistReader_instance->getTypeOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::vector<double> retrievedValueOfPrescribedFlows = netlistReader_instance->getValueOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::map<int,double> retrievedInitialPressures = netlistReader_instance->getInitialPressures().at(m_IndexOfThisNetlistLPN);
    
    // Loop over the components, assigning them (and their nodes) the appropriate properties to give the fully-described circuit:
    for (auto component = mp_CircuitDescription->components.begin(); component != mp_CircuitDescription->components.end(); component++)
    {
        (*component)->getType() = retrievedComponentTypes.back();
        retrievedComponentTypes.pop_back();

        int indexOfStartNodeInInputData = retrievedComponentStartNodes.back();
        retrievedComponentStartNodes.pop_back();
        
        // note that we're passing a boost::shared_ptr to the startNode here; if this is currently NULL, it will
        // be constructed during the call to setupPressureNode.
        setupPressureNode(indexOfStartNodeInInputData, (*component)->startNode, *component);

        // (*component)->startNode = mp_CircuitDescription.ifExistsGetNodeOtherwiseConstructNode(indexOfStartNodeInInputData);
        // (*component)->startNode->getIndex() = indexOfStartNodeInInputData;

        // (*component)->startNode->prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        // for (int prescribedPressure=0; prescribedPressure<mp_CircuitDescription.numberOfPrescribedPressures; prescribedPressure++)
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
        // (*component)->endNode = mp_CircuitDescription.ifExistsGetNodeOtherwiseConstructNode(indexOfEndNodeInInputData);
        // (*component)->endNode->getIndex() = indexOfEndNodeInInputData;
        // (*component)->endNode->prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        // for (int prescribedPressure=0; prescribedPressure<mp_CircuitDescription.numberOfPrescribedPressures; prescribedPressure++)
        // {
        //     if (retrievedListOfPrescribedPressures.at(prescribedPressure) == (*component)->endNode->getIndex())
        //     {
        //         (*component)->endNode->prescribedPressureType = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->endNode->pressure = retrievedValueOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->endNode->pressure = (*component)->endNode->pressure;
        //     }
        // }

        (*component)->prescribedFlowType = Flow_NotPrescribed;  // initialise as a default, before replacing as necessary
        for (int prescribedFlow=0; prescribedFlow<mp_CircuitDescription->numberOfPrescribedFlows; prescribedFlow++)
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
    mp_CircuitDescription->rebuildCircuitMetadata();

    // Tell the node at the 3D interface that it connects to the 3D domain:
    {
        int threeDNodeIndex = netlistReader_instance->getIndicesOfNodesAt3DInterface().at(m_IndexOfThisNetlistLPN);
        mp_CircuitDescription->initialiseNodeAndComponentAtInterface(threeDNodeIndex);
    }

    // mp_CircuitDescription.switchDiodeStatesIfNecessary();
    // mp_CircuitDescription.detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    // // Component indices are just consecutive integers by default, but sometimes non-consecutive numbering
    // // is needed; componentIndices allows for this.
    // // We initialise it now for the default case.
    // for (int ii=1; ii < mp_CircuitDescription.numberOfComponents + 1; ii++)
    // {
    //     mp_CircuitDescription.componentIndices.push_back(ii);
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
    netlistReader* netlistReader_instance = netlistReader::Instance();
    std::vector<int> retrievedListOfPrescribedPressures = netlistReader_instance->getListOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    std::vector<circuit_nodal_pressure_prescription_t> retrievedTypeOfPrescribedPressures = netlistReader_instance->getTypeOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    std::vector<double> retrievedValueOfPrescribedPressures = netlistReader_instance->getValueOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);

    // Discover whether this node has a prescribed pressure, and if so, what type:
    circuit_nodal_pressure_prescription_t typeOfPrescribedPressure = Pressure_NotPrescribed; // initialise, but chnage later if pressure is actually prescribed
    int indexOfPrescribedPressure = -1; // initialise to a nonsense value to detect errors
    for (int prescribedPressure=0; prescribedPressure<mp_CircuitDescription->numberOfPrescribedPressures; prescribedPressure++)
    {
        if (retrievedListOfPrescribedPressures.at(prescribedPressure) == indexOfNodeInInputData)
        {
            typeOfPrescribedPressure = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
            indexOfPrescribedPressure = prescribedPressure;
        }
    }
    // Get the node (or create a new node if this one hasn't been made yet)
    node = mp_CircuitDescription->ifExistsGetNodeOtherwiseConstructNode(indexOfNodeInInputData,typeOfPrescribedPressure,componentNeighbouringThisNode);    
    // node->getIndex() = indexOfNodeInInputData;
    // node->prescribedPressureType = typeOfPrescribedPressure;
    if (node->prescribedPressureType!=Pressure_NotPrescribed && node->prescribedPressureType!=Pressure_Null)
    {
        node->setPressure(retrievedValueOfPrescribedPressures.at(indexOfPrescribedPressure));
    }
}

void NetlistCircuit::selectAndBuildActiveSubcircuits()
{
    rebuildCircuitMetadata();
    detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    // Actually build the active NetlistSubcircuit classes:
    mp_activeSubcircuits.clear();
    double alfi_delt_in = m_alfi_local*m_delt;
    boost::shared_ptr<NetlistSubcircuit> toPushBack(new NetlistSubcircuit(0, mp_CircuitDescription, flow_n_ptrs, pressure_n_ptrs, alfi_delt_in, m_surfaceIndex));
    mp_activeSubcircuits.push_back(toPushBack);
}

bool NetlistCircuit::flowPermittedAcross3DInterface() const
{
	return mp_CircuitDescription->flowPermittedAcross3DInterface();
}

bool NetlistCircuit::boundaryConditionTypeHasJustChanged()
{
	return mp_CircuitDescription->boundaryConditionTypeHasJustChanged();
}

void NetlistCircuit::closeAllDiodes()
{
	mp_CircuitDescription->closeAllDiodes();
}

void NetlistCircuit::detectWhetherClosedDiodesStopAllFlowAt3DInterface()
{
	mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
}

void NetlistCircuit::switchDiodeStatesIfNecessary()
{
	mp_CircuitDescription->switchDiodeStatesIfNecessary();
}

void NetlistCircuit::rebuildCircuitMetadata()
{
	mp_CircuitDescription->rebuildCircuitMetadata();
}

// nextTimestepWrite_start will be updated and returned to caller of this function.
void NetlistCircuit::writePressuresFlowsAndVolumes(int& nextTimestepWrite_start)
{
	// All the following writes can use the same nextTimestepWrite_end to determine how far to go when looking for data to write to the file:
      int nextTimestepWrite_end = mp_CircuitDescription->getLengthOfHistoryData();
      
      {
        // Write the netlistPressures_zeroDDomainReplacement.dat
        basicFileWriter boundaryConditionPressureHistoryWriter;
        boundaryConditionPressureHistoryWriter.setFileName(m_PressureHistoryFileName);

        for (int stepToWrite=nextTimestepWrite_start; stepToWrite<nextTimestepWrite_end; stepToWrite++)
        {
          boundaryConditionPressureHistoryWriter.writeStepIndex(stepToWrite);
          for (auto node=mp_CircuitDescription->mapOfPressureNodes.begin(); node!=mp_CircuitDescription->mapOfPressureNodes.end(); node++)
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
          for (auto component=mp_CircuitDescription->components.begin(); component!=mp_CircuitDescription->components.end(); component++)
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
          for (auto component=mp_CircuitDescription->mapOfComponents.begin(); component!=mp_CircuitDescription->mapOfComponents.end(); component++)
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

std::pair<double,double> NetlistCircuit::computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt)
{
    // Get the implicit coefficients from the identified subcircuit
    std::pair<double,double> implicitCoefficients;
    // if (mp_CircuitDescription.flowPermittedAcross3DInterface())
    // {
        int numberOfCircuitsClaimingToHaveA3DInterface = 0;
        for (auto subcircuit=mp_activeSubcircuits.begin(); subcircuit!=mp_activeSubcircuits.end(); subcircuit++)
        {
            if ((*subcircuit)->m_circuitData->connectsTo3DDomain() == true)
            {
                numberOfCircuitsClaimingToHaveA3DInterface++;
                implicitCoefficients = (*subcircuit)->computeImplicitCoefficients(timestepNumber,timeAtStepNplus1,alfi_delt);
            }
        }
        assert(numberOfCircuitsClaimingToHaveA3DInterface==1);
    // }
    // else
    // {
    //     implicitCoefficients.first=1.0;
    //     implicitCoefficients.second=0.0;
    // }

    return implicitCoefficients;
}

void NetlistCircuit::updateLPN(const int timestepNumber)
{
    for (auto activeSubcircuit = mp_activeSubcircuits.begin(); activeSubcircuit != mp_activeSubcircuits.end(); activeSubcircuit++)
    {
        (*activeSubcircuit)->updateInternalPressuresVolumesAndFlows(timestepNumber);
    }
}

void NetlistCircuit::setPressureAndFlowPointers(double* pressurePointer, double* flowPointer)
{
    flow_n_ptrs.clear();
    flow_n_ptrs.push_back(flowPointer);

    pressure_n_ptrs.clear();
    pressure_n_ptrs.push_back(pressurePointer);
}

void NetlistCircuit::identifyAtomicSubcircuits()
{
    createInitialCircuitDescriptionWithoutDiodes();

    // The atomic subcircuits are those which cannot be broken down
    // into further subcircuits by diodes. More complex ones (where some 
    // valves are open, so some atomic subcircuits are joined) will 
    // appear later.
    assignComponentsToAtomicSubcircuits();
}

void NetlistCircuit::createInitialCircuitDescriptionWithoutDiodes()
{
    // Copy the data as a base for modification (i.e. the removal of the diode data)
    mp_CircuitDescriptionWithoutDiodes = boost::make_shared<CircuitData> (*mp_CircuitDescription);
    // Prepare for diode removal by clearing the data that will be rebuilt:
    // mp_CircuitDescriptionWithoutDiodes->deleteSubcircuitSpecificData();
    // Rebuild the cleared vectors, but without the diodes:
    int numberOfDiodes = 0;
    for (auto component=mp_CircuitDescriptionWithoutDiodes->components.begin(); component!=mp_CircuitDescriptionWithoutDiodes->components.end(); component++)
    {
        if((*component)->getType() == Component_Diode)
        {
            component = mp_CircuitDescriptionWithoutDiodes->components.erase(component); // returns an iterator pointing to the new location of the element after the one that just got erased
            component--; // decrement the returned iterator, so the for loop incrememnts it back to the correct (next) element.
        }
    }

    // A tidier (and untested!) version of the above, without the ugly "component --" line:
    // auto component=mp_CircuitDescriptionWithoutDiodes.components.begin();
    // while (component!=mp_CircuitDescriptionWithoutDiodes.components.end())
    // {
    //     // both branches of this if guard result in an iterator to the next component after the one just checked (whether or not an erasing occurs)
    //     bool componentNotADiode = !((*component)->type == Component_Diode);
    //     if(componentNotADiode)
    //     {
    //         component++
    //     }
    //     else
    //     {
    //         // erase the diode:
    //         component = mp_CircuitDescriptionWithoutDiodes.components.erase(component); // returns an iterator pointing to the new location of the element after the one that just got erased            
    //     }
    // }

    mp_CircuitDescriptionWithoutDiodes->rebuildCircuitMetadata();

}

void NetlistCircuit::assignComponentsToAtomicSubcircuits()
{
    // This subroutine builds m_AtomicSubcircuitsComponentsBelongsTo, which is indexed by component, as they appear in mp_CircuitDescriptionWithoutDiodes.

    // Group the components by partitioning the node indices into connected subcircuits (with diodes removed).
    // We do this by parsing the circuit to discover the topology.
    std::vector<bool> componentAssignedToASubcircuit(mp_CircuitDescriptionWithoutDiodes->numberOfComponents, false);
    m_AtomicSubcircuitsComponentsBelongsTo.insert(m_AtomicSubcircuitsComponentsBelongsTo.begin(), mp_CircuitDescriptionWithoutDiodes->numberOfComponents, -1); // initialise with nonsense value
    int subcircuitIndex = 0;
    for (int startingComponent = 0; startingComponent < mp_CircuitDescriptionWithoutDiodes->numberOfComponents; startingComponent++)
    {
        // Only assign the startingComponent to a subcircuit if it doesn't yet belong to a subcircuit
        if (!componentAssignedToASubcircuit.at(startingComponent))
        {
            // Because we immediately proceed to find the whole subcircuit, this is definitely a new
            // subcircuit we're starting now.
            m_AtomicSubcircuitsComponentsBelongsTo.at(startingComponent) = subcircuitIndex;
            componentAssignedToASubcircuit.at(startingComponent) = true;

            // search the rest of the list for the remaining components belonging to this subcircuit:
            // Make a list to be populated with nodes belonging to the subcircuit, as they're found.
            std::set<int> nodesInSubcircuit;
            nodesInSubcircuit.insert(mp_CircuitDescriptionWithoutDiodes->components.at(startingComponent)->startNode->getIndex());
            nodesInSubcircuit.insert(mp_CircuitDescriptionWithoutDiodes->components.at(startingComponent)->endNode->getIndex());
            // Find all other components belonging to this subcircuit. We must iterate (while-loop) to keep
            // finding new connections, until nothing changes on an iteration:
            bool newComponentFoundOnLastIteration = true;
            while (newComponentFoundOnLastIteration)
            {
                newComponentFoundOnLastIteration = false;
                for (int potentialSubcircuitComponent = startingComponent+1; potentialSubcircuitComponent < mp_CircuitDescriptionWithoutDiodes->numberOfComponents; potentialSubcircuitComponent++)
                {
                    if (componentAssignedToASubcircuit.at(potentialSubcircuitComponent) == false)
                    {
                        // get the nodes of the current component
                        int componentStartNode = mp_CircuitDescriptionWithoutDiodes->components.at(potentialSubcircuitComponent)->startNode->getIndex();
                        int componentEndNode = mp_CircuitDescriptionWithoutDiodes->components.at(potentialSubcircuitComponent)->endNode->getIndex();

                        // Check whether either of these nodes connect potentialSubcircuitComponent to the current subcircuit.
                        // If so, assign potentialSubcircuitComponent to the subcircuit.
                        if(nodesInSubcircuit.count(componentStartNode) == 1 ||
                            nodesInSubcircuit.count(componentEndNode) == 1) // the count() method is an unfortunate name; it can only return 0 or 1 because of std::set's unique element rule.
                        {
                            // Note that a re-parse is required, because the subcircuit grew:
                            newComponentFoundOnLastIteration = true;

                            m_AtomicSubcircuitsComponentsBelongsTo.at(potentialSubcircuitComponent) = subcircuitIndex;
                            componentAssignedToASubcircuit.at(potentialSubcircuitComponent) = true;
                            
                            // Record that these nodes belong to this subcircuit.
                            // Note that the insert will simply do nothing in the case where the node
                            // is already recorded as being part of the subcircuit. This is fine.
                            nodesInSubcircuit.insert(componentStartNode);
                            nodesInSubcircuit.insert(componentEndNode);
                        }
                    }
                }
            }

            subcircuitIndex++;
        }
    }
    // We have now also counted the atomic subcircuits. Record that count!
    m_NumberOfAtomicSubcircuits = subcircuitIndex;

    // Sanity check, to make sure all the components got assigned to a subcircuit:
    for (auto componentAssigned = componentAssignedToASubcircuit.begin(); componentAssigned!=componentAssignedToASubcircuit.end(); componentAssigned++)
    {
        assert(*componentAssigned);
    }

}

void NetlistCircuit::createAtomicSubcircuitDescriptions()
{
    // Make the data for the atomic subcircuits:
    assert(m_CircuitDataForAtomicSubcircuits.empty());
    for(int currentSubcircuitIdx=0; currentSubcircuitIdx<m_NumberOfAtomicSubcircuits; currentSubcircuitIdx++)
    {
        boost::shared_ptr<CircuitData> toPushBack(new CircuitData(m_hstep));
        m_CircuitDataForAtomicSubcircuits.push_back(toPushBack);
        // Get a reference to the new CircuitData, just to avoid filling the screen with an infinite amount of code, below.
        CircuitData& currentSubcircuit = *(m_CircuitDataForAtomicSubcircuits.back());

        // Scoping unit for componentLocationIdx:
        {
            int componentLocationIdx = 0;
            for (auto componentsSubcircuitIdx=m_AtomicSubcircuitsComponentsBelongsTo.begin(); componentsSubcircuitIdx!=m_AtomicSubcircuitsComponentsBelongsTo.end(); componentsSubcircuitIdx++)
            {
                if (*componentsSubcircuitIdx == currentSubcircuitIdx)
                {
                    currentSubcircuit.components.push_back( mp_CircuitDescriptionWithoutDiodes->components.at(componentLocationIdx) );
                }
                componentLocationIdx++;
            }
        }

        currentSubcircuit.rebuildCircuitMetadata();
    }

    // Give each subcircuit its index in the parent vector:
    for (int ii=0; ii<m_CircuitDataForAtomicSubcircuits.size(); ii++)
    {
        m_CircuitDataForAtomicSubcircuits.at(ii)->index=ii;
    }

    // Identify the composite subcircuits (with legality assigned by user?)

    // Make the data for the composite subcircuits

    // Make sure that we have a sensible data structure so that the list of active subcircuits
    // can easily be mapped to the actual active circuits themselves.
}

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
    // for (auto node=mp_CircuitDescriptionWithoutDiodes.mapOfPressureNodes.begin(); node!=mp_CircuitDescriptionWithoutDiodes.mapOfPressureNodes.end(); node++)
    // Cycle and store the history pressures
    for (auto node=mp_CircuitDescription->mapOfPressureNodes.begin(); node!=mp_CircuitDescription->mapOfPressureNodes.end(); node++)
    {
        // Store the pressure for writing to output file:
        node->second->m_entirePressureHistory.push_back(node->second->getPressure());

        if (node->second->hasHistoryPressure)
        {
            node->second->historyPressure = node->second->getPressure();
        }
    }

    // for (auto component=mp_CircuitDescriptionWithoutDiodes.mapOfComponents.begin(); component!=mp_CircuitDescriptionWithoutDiodes.mapOfComponents.end(); component++)
    // Cycle and store the history flows
    for (auto component=mp_CircuitDescription->mapOfComponents.begin(); component!=mp_CircuitDescription->mapOfComponents.end(); component++)
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
    for (auto component=mp_CircuitDescription->mapOfComponents.begin(); component!=mp_CircuitDescription->mapOfComponents.end(); component++)
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

std::pair<boundary_data_t,double> NetlistCircuit::computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber)
{
    std::pair<boundary_data_t,double> pressureOrFlowToReturn;

    int counterToAvoidErrors = 0;
    for (auto subcircuit = mp_activeSubcircuits.begin(); subcircuit != mp_activeSubcircuits.end(); subcircuit++)
    {
        pressureOrFlowToReturn = (*subcircuit)->computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(timestepNumber);
        counterToAvoidErrors ++;
    }
    assert(counterToAvoidErrors == 1); // This is not built to deal with multiple subcircuits yet. If you go back to using multiple subcircuits, you'll need to fix this.

    return pressureOrFlowToReturn;
}

void NetlistCircuit::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    selectAndBuildActiveSubcircuits();
    cycleToSetHistoryPressuresFlowsAndVolumes();
}

void NetlistCircuit::finalizeLPNAtEndOfTimestep()
{
    switchDiodeStatesIfNecessary();
}

boost::shared_ptr<CircuitComponent> NetlistCircuit::getComponentByInputDataIndex(const int componentIndex)
{
    return mp_CircuitDescription->getComponentByInputDataIndex(componentIndex);
}

void NetlistZeroDDomainCircuit::createCircuitDescription()
{
    // This function creates the internal CircuitData class format for a zero-D
    // Netlist, non-boundary-condition replacement for the 3D domain, for pure zero-D simulation.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    netlistReader* netlistReader_instance = netlistReader::Instance();

    // Make the appropriate class to store the 3D domain replacement circuit data:
    mp_CircuitDescription = boost::shared_ptr<Netlist3DDomainReplacementCircuitData> (new Netlist3DDomainReplacementCircuitData(m_hstep,m_numberOfNetlistsUsedAsBoundaryConditions));

    // we'll have a resistor for each netlist used as a boundary condition, plus one VolumeTrackingPressureChamber:
    mp_CircuitDescription->numberOfComponents = 2*m_numberOfNetlistsUsedAsBoundaryConditions + 1;

    // The components will be arranged in a star (all start-nodes are connected together at a single point, all end-nodes only connect to one component)
    mp_CircuitDescription->numberOfPressureNodes = 2*m_numberOfNetlistsUsedAsBoundaryConditions + 2;

    // The pressures will initially be received from the boundary conditions at the boundary condition interface nodes,
    // plus the zero-pressure prescription on the base of the VolumeTrackingPressureChamber:
    mp_CircuitDescription->numberOfPrescribedPressures = 1;
    // None initially:
    mp_CircuitDescription->numberOfPrescribedFlows = m_numberOfNetlistsUsedAsBoundaryConditions;

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
    assert(mp_CircuitDescription->components.empty());
    for (int ii=0; ii<mp_CircuitDescription->numberOfComponents; ii++)
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

        mp_CircuitDescription->components.push_back(boost::shared_ptr<CircuitComponent> (toPushBack));
        mp_CircuitDescription->components.back()->setIndex(toOneIndexing(ii));
    }

    // Obtain the component- and node-level data for the circuit, for moving into the appropriate data structure CircuitData
    // We want to pop off the component types as we use them, but starting from the beginning of the vector. To do this, we reverse
    // the vector and then pop from the new end.
    std::reverse(componentTypes.begin(), componentTypes.end());
    // Create the component start node info (all are the same: we make it the highest-indexed node.)
    
    // int centralNodeIndex = mp_CircuitDescription->numberOfComponents + 1;
    // std::vector<int> componentStartNodes(mp_CircuitDescription->numberOfComponents, centralNodeIndex);

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
    for (int internalDomainResistorEndNodeIndex = m_numberOfNetlistsUsedAsBoundaryConditions + 3; internalDomainResistorEndNodeIndex <= mp_CircuitDescription->numberOfPressureNodes; internalDomainResistorEndNodeIndex++)
    {
        componentEndNodes.push_back(internalDomainResistorEndNodeIndex);
    }
    // Finally, do the end-node at the base of the compliance chamber:
    componentEndNodes.push_back(m_numberOfNetlistsUsedAsBoundaryConditions + 1);


    std::vector<int> componentStartNodes;
    // All the components, except for the dP/dQ resistors, have the same
    // start-node (the centre-point of the star-shape of this domain).
    // We do the dP/dQ resistors first:
    for (int dpDqResistorStartNodeIndex = m_numberOfNetlistsUsedAsBoundaryConditions+3; dpDqResistorStartNodeIndex <= mp_CircuitDescription->numberOfPressureNodes; dpDqResistorStartNodeIndex++)
    {
        componentStartNodes.push_back(dpDqResistorStartNodeIndex);
    }
    // All the other components share the same start node, at the centre of the star-shaped domain:
    for (int componentIndex = m_numberOfNetlistsUsedAsBoundaryConditions+1; componentIndex <= mp_CircuitDescription->numberOfComponents; componentIndex++)
    {
        componentStartNodes.push_back(m_numberOfNetlistsUsedAsBoundaryConditions+2);
    }

    std::reverse(componentStartNodes.begin(), componentStartNodes.end()); // actually no point in this call, but it's tidier to leave it here for symmetry with the related calls below.

    std::reverse(componentEndNodes.begin(), componentEndNodes.end());

    // std::vector<double> componentParameterValues(m_numberOfNetlistsUsedAsBoundaryConditions, m_oneResistanceToGiveEachResistor);
    // componentParameterValues.push_back(m_elastanceToGiveVolumeTrackingPressureChamber);

    boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_CircuitDescription);
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
    for (int node=m_numberOfNetlistsUsedAsBoundaryConditions+2; node <= mp_CircuitDescription->numberOfPressureNodes; node++)
    {
        initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }

    // Do the final node (at the base of the pressure chamber):
    // initialPressures.insert(std::make_pair(mp_CircuitDescription->numberOfPressureNodes, 0.0));
    
    // Loop over the components, assigning them (and their nodes) the appropriate properties to give the fully-described circuit:
    for (auto component = mp_CircuitDescription->components.begin(); component != mp_CircuitDescription->components.end(); component++)
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
        // for (int prescribedFlow=0; prescribedFlow<mp_CircuitDescription->numberOfPrescribedFlows; prescribedFlow++)
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
    mp_CircuitDescription->rebuildCircuitMetadata();

    // Tell the node at the 3D interface that it connects to the 3D domain:
    {
        std::vector<int> threeDNodeIndices;
        for (int nodeIndex=1; nodeIndex<=m_numberOfNetlistsUsedAsBoundaryConditions; nodeIndex++)
        {
        	threeDNodeIndices.push_back(nodeIndex);
        }
        boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_CircuitDescription);
        downcastDomainReplacementCircuitData->initialiseNodesAndComponentsAtInterface_vector(threeDNodeIndices);
    }

    // mp_CircuitDescription->switchDiodeStatesIfNecessary();
    // mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    // // Component indices are just consecutive integers by default, but sometimes non-consecutive numbering
    // // is needed; componentIndices allows for this.
    // // We initialise it now for the default case.
    // for (int ii=1; ii < mp_CircuitDescription->numberOfComponents + 1; ii++)
    // {
    //     mp_CircuitDescription->componentIndices.push_back(ii);
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
    node = mp_CircuitDescription->ifExistsGetNodeOtherwiseConstructNode(indexOfNodeInInputData,typeOfPrescribedPressure,componentNeighbouringThisNode);    
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
    return mp_CircuitDescription;
}

void NetlistZeroDDomainCircuit::selectAndBuildActiveSubcircuits()
{
    mp_CircuitDescription->rebuildCircuitMetadata();
    // mp_CircuitDescription.switchDiodeStatesIfNecessary();
    // mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    // Actually build the active NetlistSubcircuit classes:
    mp_activeSubcircuits.clear();
    double alfi_delt_in = m_alfi_local*m_delt;
    boost::shared_ptr<NetlistSubcircuit> newNetlistZeroDDomain(new NetlistSubcircuit(0, mp_CircuitDescription, flow_n_ptrs, pressure_n_ptrs,alfi_delt_in,m_surfaceIndex));
    newNetlistZeroDDomain->setThisisA3DDomainReplacement();
    mp_activeSubcircuits.push_back(newNetlistZeroDDomain);

}

void NetlistZeroDDomainCircuit::setBoundaryPrescriptionsAndBoundaryConditionTypes(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate)
{
    boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_CircuitDescription);
    assert(downcastDomainReplacementCircuitData!=NULL);
	downcastDomainReplacementCircuitData->setBoundaryPrescriptionsAndBoundaryConditionTypes(boundaryFlowsOrPressuresAsAppropriate);
}

std::vector<double> NetlistZeroDDomainCircuit::getBoundaryPressures()
{
	std::vector<double> pressures;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		pressures.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->startNode->getPressure());
	}
	return pressures;
}

std::vector<double> NetlistZeroDDomainCircuit::getBoundaryFlows()
{
	std::vector<double> flows;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		flows.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->flow);
	}
	return flows;
}

void NetlistZeroDDomainCircuit::solveSystem(const int timestepNumber)
{
	mp_activeSubcircuits.at(0)->buildAndSolveLinearSystem(timestepNumber);
}

void NetlistZeroDDomainCircuit::setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients, std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAtBoundaries)
{
    int dpDqResistorIndex=0;
    for (auto component=mp_CircuitDescription->components.begin(); component!=mp_CircuitDescription->components.end(); component++)
    {
        boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcastDomainReplacementCircuitData = boost::dynamic_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_CircuitDescription);
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