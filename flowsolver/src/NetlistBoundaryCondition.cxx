#include "NetlistBoundaryCondition.hxx"
#include "fileReaders.hxx"
#include "datatypesInCpp.hxx"
#include <cassert>
#include <algorithm>
#include <boost/make_shared.hpp>

// Statics:
int NetlistBoundaryCondition::numberOfInitialisedNetlistLPNs = 0;

void NetlistBoundaryCondition::initialiseModel()
{
    // Get the input data
    createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    identifyAtomicSubcircuits();

    // chop up the Circuitdata into subcircuits (including removing all mention of diodes from that data now)
    createAtomicSubcircuitDescriptions();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    mp_CircuitDescription->closeAllDiodes();
    mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    
    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    // m_atomicSubcircuitConnectionManager = boost::shared_ptr<AtomicSubcircuitConnectionManager>( toPassToSharedPtr );
}

void NetlistBoundaryCondition::finalizeLPNAtEndOfTimestep()
{
    mp_CircuitDescription->switchDiodeStatesIfNecessary();
}

void NetlistBoundaryCondition::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    selectAndBuildActiveSubcircuits();

    cycleToSetHistoryPressuresFlowsAndVolumes();
}

std::pair<double,double> NetlistBoundaryCondition::computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt)
{
    // Get the implicit coefficients from the identified subcircuit
    std::pair<double,double> implicitCoefficients;
    // if (mp_CircuitDescription.flowPermittedAcross3DInterface())
    // {
        int numberOfCircuitsClaimingToHaveA3DInterface = 0;
        for (auto subcircuit=m_activeSubcircuits.begin(); subcircuit!=m_activeSubcircuits.end(); subcircuit++)
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

    std::cout << "setting implicit coeffs in NetlistBoundaryCondition.cxx: " << implicitCoefficients.first << " " << implicitCoefficients.second << std::endl;
    return implicitCoefficients;
}

std::pair<boundary_data_t,double> NetlistBoundaryCondition::computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber)
{
    std::pair<boundary_data_t,double> pressureOrFlowToReturn;
    int counterToAvoidErrors = 0;
    for (auto subcircuit = m_activeSubcircuits.begin(); subcircuit != m_activeSubcircuits.end(); subcircuit++)
    {
        pressureOrFlowToReturn = (*subcircuit)->computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(timestepNumber);
        counterToAvoidErrors ++;
    }
    assert(counterToAvoidErrors == 1); // This is not built to deal with multiple subcircuits yet. If you go back to using multiple subcircuits, you'll need to fix this.
    return pressureOrFlowToReturn;
}

void NetlistBoundaryCondition::updateLPN()
{
    for (auto activeSubcircuit = m_activeSubcircuits.begin(); activeSubcircuit != m_activeSubcircuits.end(); activeSubcircuit++)
    {
        (*activeSubcircuit)->updateInternalPressuresVolumesAndFlows();
    }
}

void NetlistBoundaryCondition::selectAndBuildActiveSubcircuits()
{
    // // Identify which valves are open
    // for (int diodeIdx=0; diodeIdx< m_atomicSubcircuitConnectionManager->getNumberOfDiodes(); diodeIdx++)
    // {
    //     int diodeStartNodeIdx = m_atomicSubcircuitConnectionManager->m_diodeIndexingMap.at(diodeIdx)->startNode->indexInInputData;
    //     int diodeEndNodeIdx = m_atomicSubcircuitConnectionManager->m_diodeIndexingMap.at(diodeIdx)->endNode->indexInInputData;

    //     if(m_PressuresInLPN.at(diodeStartNodeIdx) >= m_PressuresInLPN.at(diodeEndNodeIdx))
    //     {
    //         m_atomicSubcircuitConnectionManager->setDiodeOpen(diodeIdx,true);
    //     }
    //     else
    //     {
    //         m_atomicSubcircuitConnectionManager->setDiodeOpen(diodeIdx,false);
    //     }
    // }

    // // populate ActiveSubcircuits
    // m_activeSubcircuitCircuitData.clear();
    // std::vector<bool> atomicSubcircuitAssignedToActiveCircuit(m_NumberOfAtomicSubcircuits,false);
    // if (m_atomicSubcircuitConnectionManager->getNumberOfDiodes() > 0)
    // {
    //     for (int diodeIdx=0; diodeIdx<m_atomicSubcircuitConnectionManager->getNumberOfDiodes(); diodeIdx++)
    //     {
    //         if (m_atomicSubcircuitConnectionManager->diodeIsOpen(diodeIdx))
    //         {
    //             int diodeStartNodeIdx = m_atomicSubcircuitConnectionManager->m_diodeIndexingMap.at(diodeIdx)->startNode->indexInInputData;
    //             int diodeEndNodeIdx = m_atomicSubcircuitConnectionManager->m_diodeIndexingMap.at(diodeIdx)->endNode->indexInInputData;
    //             for (auto activeSubcircuit = m_activeSubcircuitCircuitData.begin(); activeSubcircuit!=m_activeSubcircuitCircuitData.end(); activeSubcircuit++)
    //             {
    //                 // If startNode belongs to activeSubcircuit...
    //                 bool diodeStartNodeBelongsToActiveSubcircuit = ((*activeSubcircuit)->mapOfPressureNodes.count(diodeStartNodeIdx) == 1); 
    //                 bool atomicSubcircuitAtEndNodeOfDiodeNotAssignedToActiveSubcircuit = !(atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->index));
    //                 if (diodeStartNodeBelongsToActiveSubcircuit && atomicSubcircuitAtEndNodeOfDiodeNotAssignedToActiveSubcircuit)
    //                 {
    //                     if (diodeStartNodeConnectsCircuit(diodeIdx))
    //                     {
    //                         // Connect the end-node subcircuit to this activeSubcircuit:
    //                         (*activeSubcircuit)->components.reserve((*activeSubcircuit)->components.size() + m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->components.size());
    //                         for (auto componentToInsert = m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->components.begin();
    //                              componentToInsert != m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->components.end();
    //                              componentToInsert++)
    //                         {
    //                             boost::shared_ptr<CircuitComponent> componentToPushBack(&(**componentToInsert));
    //                             (*activeSubcircuit)->components.push_back(componentToPushBack);
    //                         }
    //                         atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->index) = true;
    //                     }
    //                     else // Diode start node is just a single node, connecting to no other componenet (but must have prescribed pressure!)
    //                     {
    //                         // In this case, we just re-index the diode end node to be the diode start node, throughout the whole subcircuit.
    //                         (*activeSubcircuit)->com
    //                     }
    //                 }
    //                 // If endNode belongs to activeSubcircuit...
    //                 bool diodeEndNodeBelongsToActiveSubcircuit = ((*activeSubcircuit)->mapOfPressureNodes.count(diodeEndNodeIdx) == 1);
    //                 bool atomicSubcircuitAtStartNodeOfDiodeNotAssignedToActiveSubcircuit = !(atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->index));
    //                 if (diodeEndNodeBelongsToActiveSubcircuit && atomicSubcircuitAtStartNodeOfDiodeNotAssignedToActiveSubcircuit)
    //                 {
    //                     // Add the atomic subcircuit to this activeSubcircuit:
    //                     (*activeSubcircuit)->components.reserve((*activeSubcircuit)->components.size() + m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->components.size());
    //                     for (auto componentToInsert = m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->components.begin();
    //                          componentToInsert != m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->components.end();
    //                          componentToInsert++)
    //                     {
    //                         boost::shared_ptr<CircuitComponent> componentToPushBack(&(**componentToInsert));
    //                         (*activeSubcircuit)->components.push_back(componentToPushBack);
    //                     }
    //                     atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->index) = true;
    //                 }
    //             }
    //             // If the atomic subcircuits at the start or end of the diode remain unassigned, add them as a new activeSubcircuit:
    //             if(atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->index) = false)
    //             {
    //                 m_activeSubcircuitCircuitData.push_back( m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx) );
    //                 atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToStartNode(diodeIdx)->index) = true;
    //             }
    //             if(atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->index) = false)
    //             {
    //                 m_activeSubcircuitCircuitData.push_back( m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx) );
    //                 atomicSubcircuitAssignedToActiveCircuit.at(m_atomicSubcircuitConnectionManager->getCircuitConnectedToEndNode(diodeIdx)->index) = true;
    //             }
    //         }
    //     }
    // }
    // else if (m_atomicSubcircuitConnectionManager->getNumberOfDiodes() == 0)
    // {
    //     // There's only one circuit (i.e. no diodes) in this case, so we just get the shared_ptr for its data and push it back.
    //     m_activeSubcircuitCircuitData.push_back(m_CircuitDataForAtomicSubcircuits.at(0));
    //     atomicSubcircuitAssignedToActiveCircuit.at(0) = true;
    // }
    // else
    // {
    //     throw std::runtime_error("EE: Invalid number of atomic subcircuits detected.");
    // }

    // // Build the active subcircuits' circuit metadata:
    // for (auto subcircuit=m_activeSubcircuitCircuitData.begin(); subcircuit!=m_activeSubcircuitCircuitData.end(); subcircuit++)
    // {
    //     (*subcircuit)->rebuildCircuitMetadata();
    // }

    // // Check all the atomic subcircuits got assigned:
    // for (auto assigned=atomicSubcircuitAssignedToActiveCircuit.begin(); assigned!=atomicSubcircuitAssignedToActiveCircuit.end(); assigned++)
    // {
    //     if(!*assigned)
    //     {
    //         throw std::runtime_error("EE: At least one atomic subcircuit not assigned during netlist assembly; likely a problem with the diodes.");
    //     }
    // }
    
    // // The atomic subcircuits are all now properly assigned to active subcircuits.
    // // Next, we go on to change the node indexing so that open diode start-nodes are 
    // // identified with their diode's end node.
    // // To do this, we replace all references to the diode's start-node with references
    // // to its end node.
    // //
    // // This completes the connection between the atomic subcircuits.
    // for (int diodeIdx=0; diodeIdx< m_atomicSubcircuitConnectionManager->getNumberOfDiodes(); diodeIdx++)
    // {
    //     if (m_atomicSubcircuitConnectionManager->diodeIsOpen(diodeIdx))
    //     {
    //         int diodeStartNodeIdx = m_atomicSubcircuitConnectionManager->m_diodeIndexingMap.at(diodeIdx)->startNode->indexInInputData;
    //         CircuitPressureNode diodeEndNode = *(m_atomicSubcircuitConnectionManager->m_diodeIndexingMap.at(diodeIdx)->endNode);
    //         for (auto activeSubcircuit = m_activeSubcircuitCircuitData.begin(); activeSubcircuit!=m_activeSubcircuitCircuitData.end(); activeSubcircuit++)
    //         {
    //             for (auto component=(*activeSubcircuit)->components.begin(); component!=(*activeSubcircuit)->components.end(); component++)
    //             {
    //                 //\todo consider cases at the edge of the domain here,
    //                 // where e.g. the diode start node has a prescribed pressure
    //                 // and is monopolar. The code currently here will over-write and
    //                 // break that pressure imposition!
    //                 if ((*component)->startNode->indexInInputData == diodeStartNodeIdx)
    //                 {
    //                     (*component)->startNode->indexInInputData = diodeEndNode.indexInInputData;
    //                 }
    //                 if ((*component)->endNode->indexInInputData == diodeStartNodeIdx)
    //                 {
    //                     (*component)->endNode->indexInInputData = diodeEndNode.indexInInputData;
    //                 }
    //             }
    //         }
    //     }
    // }

    // for (auto subcircuit=m_activeSubcircuitCircuitData.begin(); subcircuit!=m_activeSubcircuitCircuitData.end(); subcircuit++)
    // {
    //     (*subcircuit)->generateNodeAndComponentIndicesLocalToSubcircuit();
    // }

    // // Actually build the active NetlistSubcircuit classes:
    // m_activeSubcircuits.clear();
    // for (int subcircuitIdx=0; subcircuitIdx<m_activeSubcircuitCircuitData.size(); subcircuitIdx++)
    // {
    //     boost::shared_ptr<NetlistSubcircuit> toPushBack(new NetlistSubcircuit(subcircuitIdx, *(m_activeSubcircuitCircuitData.at(subcircuitIdx)), flow_n_ptr));
    //     m_activeSubcircuits.push_back(toPushBack);
    // }
    mp_CircuitDescription->rebuildCircuitMetadata();
    mp_CircuitDescription->generateNodeAndComponentIndicesLocalToSubcircuit();
    // mp_CircuitDescription.switchDiodeStatesIfNecessary();
    mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    // Actually build the active NetlistSubcircuit classes:
    m_activeSubcircuits.clear();
    double alfi_delt_in = alfi_local*delt;
    boost::shared_ptr<NetlistSubcircuit> toPushBack(new NetlistSubcircuit(0, mp_CircuitDescription, flow_n_ptrs, pressure_n_ptrs,alfi_delt_in,surfaceIndex));
    m_activeSubcircuits.push_back(toPushBack);

}

void NetlistBoundaryCondition::identifyAtomicSubcircuits()
{

    createInitialCircuitDescriptionWithoutDiodes();

    // The atomic subcircuits are those which cannot be broken down
    // into further subcircuits by diodes. More complex ones (where some 
    // valves are open, so some atomic subcircuits are joined) will 
    // appear later.
    assignComponentsToAtomicSubcircuits();

}

void NetlistBoundaryCondition::createCircuitDescription()
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
            toPushBack = new VolumeTrackingPressureChamber(hstep,thisIsARestartedSimulation);
        }
        else
        {
            toPushBack = new CircuitComponent(hstep,thisIsARestartedSimulation);
        }

        mp_CircuitDescription->components.push_back(boost::shared_ptr<CircuitComponent> (toPushBack));
        mp_CircuitDescription->components.back()->indexInInputData = ii+1; // This uses input data indexing, which is one-indexed. We add 1 to achieve this here.
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
        (*component)->type = retrievedComponentTypes.back();
        retrievedComponentTypes.pop_back();

        int indexOfStartNodeInInputData = retrievedComponentStartNodes.back();
        retrievedComponentStartNodes.pop_back();
        
        // note that we're passing a boost::shared_ptr to the startNode here; if this is currently NULL, it will
        // be constructed during the call to setupPressureNode.
        setupPressureNode(indexOfStartNodeInInputData, (*component)->startNode, *component);

        // (*component)->startNode = mp_CircuitDescription.ifExistsGetNodeOtherwiseConstructNode(indexOfStartNodeInInputData);
        // (*component)->startNode->indexInInputData = indexOfStartNodeInInputData;

        // (*component)->startNode->prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        // for (int prescribedPressure=0; prescribedPressure<mp_CircuitDescription.numberOfPrescribedPressures; prescribedPressure++)
        // {
        //     if (retrievedListOfPrescribedPressures.at(prescribedPressure) == (*component)->startNode->indexInInputData)
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
        // (*component)->endNode->indexInInputData = indexOfEndNodeInInputData;
        // (*component)->endNode->prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        // for (int prescribedPressure=0; prescribedPressure<mp_CircuitDescription.numberOfPrescribedPressures; prescribedPressure++)
        // {
        //     if (retrievedListOfPrescribedPressures.at(prescribedPressure) == (*component)->endNode->indexInInputData)
        //     {
        //         (*component)->endNode->prescribedPressureType = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->endNode->pressure = retrievedValueOfPrescribedPressures.at(prescribedPressure);
        //         (*component)->endNode->pressure = (*component)->endNode->pressure;
        //     }
        // }

        (*component)->prescribedFlowType = Flow_NotPrescribed;  // initialise as a default, before replacing as necessary
        for (int prescribedFlow=0; prescribedFlow<mp_CircuitDescription->numberOfPrescribedFlows; prescribedFlow++)
        {
            if (retrievedListOfPrescribedFlows.at(prescribedFlow) == (*component)->indexInInputData)
            {
                (*component)->prescribedFlowType = retrievedTypeOfPrescribedFlows.at(prescribedFlow);
                (*component)->valueOfPrescribedFlow = retrievedValueOfPrescribedFlows.at(prescribedFlow);
            }
        }

        (*component)->currentParameterValue = retrievedComponentParameterValues.back();
        // make a copy of this value so we can reset it if necessary. Used for e.g. diode state changes.
        (*component)->parameterValueFromInputData = (*component)->currentParameterValue;
        retrievedComponentParameterValues.pop_back();

        (*component)->startNode->setPressure(retrievedInitialPressures.at((*component)->startNode->indexInInputData));
        (*component)->endNode->setPressure(retrievedInitialPressures.at((*component)->endNode->indexInInputData));


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
void NetlistBoundaryCondition::setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode)
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
    // node->indexInInputData = indexOfNodeInInputData;
    // node->prescribedPressureType = typeOfPrescribedPressure;
    if (node->prescribedPressureType!=Pressure_NotPrescribed && node->prescribedPressureType!=Pressure_Null)
    {
        node->setPressure(retrievedValueOfPrescribedPressures.at(indexOfPrescribedPressure));
    }
}

void NetlistBoundaryCondition::createInitialCircuitDescriptionWithoutDiodes()
{
    // Copy the data as a base for modification (i.e. the removal of the diode data)
    mp_CircuitDescriptionWithoutDiodes = boost::make_shared<CircuitData> (*mp_CircuitDescription);
    // Prepare for diode removal by clearing the data that will be rebuilt:
    // mp_CircuitDescriptionWithoutDiodes->deleteSubcircuitSpecificData();
    // Rebuild the cleared vectors, but without the diodes:
    int numberOfDiodes = 0;
    for (auto component=mp_CircuitDescriptionWithoutDiodes->components.begin(); component!=mp_CircuitDescriptionWithoutDiodes->components.end(); component++)
    {
        if((*component)->type == Component_Diode)
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

void NetlistBoundaryCondition::assignComponentsToAtomicSubcircuits()
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
            nodesInSubcircuit.insert(mp_CircuitDescriptionWithoutDiodes->components.at(startingComponent)->startNode->indexInInputData);
            nodesInSubcircuit.insert(mp_CircuitDescriptionWithoutDiodes->components.at(startingComponent)->endNode->indexInInputData);
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
                        int componentStartNode = mp_CircuitDescriptionWithoutDiodes->components.at(potentialSubcircuitComponent)->startNode->indexInInputData;
                        int componentEndNode = mp_CircuitDescriptionWithoutDiodes->components.at(potentialSubcircuitComponent)->endNode->indexInInputData;

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

void NetlistBoundaryCondition::createAtomicSubcircuitDescriptions()
{
    // Make the data for the atomic subcircuits:
    assert(m_CircuitDataForAtomicSubcircuits.empty());
    for(int currentSubcircuitIdx=0; currentSubcircuitIdx<m_NumberOfAtomicSubcircuits; currentSubcircuitIdx++)
    {
        boost::shared_ptr<CircuitData> toPushBack(new CircuitData(hstep));
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

void NetlistBoundaryCondition::cycleToSetHistoryPressuresFlowsAndVolumes()
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

boost::shared_ptr<CircuitData> NetlistBoundaryCondition::getCircuitDescription()
{
    return mp_CircuitDescription;
}

// Processes the binaryMask for setting Dirichlet conditions.
// This boundary condition knows which mesh nodes lie at its surface (checked by the assert),
// and it sets 1 in binaryMask at the appropriate location for these nodes, if the boundary
// condition type is currently Dirichlet.
void NetlistBoundaryCondition::setDirichletConditionsIfNecessary(int* const binaryMask)
{
  if(mp_CircuitDescription->flowPermittedAcross3DInterface())
  {
    assert(hasListOfMeshNodesAtThisBoundary);
    // set ones in the binaryMask at the locations necessary to impose Dirichlet at this surface
    for (auto node=listOfMeshNodesAtThisBoundary.begin(); node!=listOfMeshNodesAtThisBoundary.end(); node++)
    {
      binaryMask[*node] = 0;
    }
  }
}

void NetlistBoundaryCondition::setPressureAndFlowPointers(double* pressurePointer, double* flowPointer)
{
    flow_n_ptrs.clear();
    flow_n_ptrs.push_back(flowPointer);

    pressure_n_ptrs.clear();
    pressure_n_ptrs.push_back(pressurePointer);
}