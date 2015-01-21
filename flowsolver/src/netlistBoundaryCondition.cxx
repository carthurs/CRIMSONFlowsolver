#include "netlistBoundaryCondition.hxx"
#include "fileReaders.hxx"
#include "datatypesInCpp.hxx"
#include <assert.h>
#include <algorithm>
#include <bitset>

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
    
    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    buildDiodeMetadata();

    // Make and store the actual subcircuits, where at least one valve is open
    createCompositeSubcircuitDescriptions();

}

std::pair<double,double> NetlistBoundaryCondition::computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt)
{
    // Idetify the appropriate subcircuit
    selectActiveSubcircuits();

    // Get the implicit coefficients from the identified subcircuit
    std::pair<double,double> implicitCoefficients;
    implicitCoefficients = m_SubcircuitsAsDelineatedByDiodes.at(m_ActiveSubcircuits.indexOfSubcircuitAt3DInterface)->computeImplicitCoefficients(timestepNumber,timeAtStepNplus1,alfi_delt);

    return implicitCoefficients;
}

void NetlistBoundaryCondition::updateLPN()
{
    // Identify which circuits are in play

    // Tell them to update

    // Recover their computed pressures and flows, store them in the NetlistBoundaryCondition

    passPressuresAndFlowsToAllSubcircuits();

}

void NetlistBoundaryCondition::passPressuresAndFlowsToAllSubcircuits()
{
    // loop the subcircuits, passing the data to them
}

void NetlistBoundaryCondition::selectActiveSubcircuits()
{
    // Identify which valves are open

    // populate ActiveSubcircuits
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
    m_CircuitDescription.numberOfComponents = netlistReader_instance->getNumberOfComponents().at(m_IndexOfThisNetlistLPN);
    m_CircuitDescription.numberOfPressureNodes = netlistReader_instance->getNumberOfPressureNodes().at(m_IndexOfThisNetlistLPN);
    m_CircuitDescription.numberOfPrescribedPressures = netlistReader_instance->getNumberOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    m_CircuitDescription.numberOfPrescribedFlows = netlistReader_instance->getNumberOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);

    // Prepare space for the components in the circuit:
    assert(m_CircuitDescription.components.empty());
    for (int ii=0; ii<m_CircuitDescription.numberOfComponents; ii++)
    {
        m_CircuitDescription.components.push_back(new CircuitComponent);
        m_CircuitDescription.components.back().indexInCircuit = ii;
    }


    // Obtain the component- and node-level data for the circuit, for moving into the appropriate data structure CircuitData
    std::vector<circuit_component_t> retrievedComponentTypes = netlistReader_instance->getComponentTypes().at(m_IndexOfThisNetlistLPN);
    // We want to pop off the component types as we use them, but starting from the beginning of the vector. To do this, we reverse
    // the vector and then pop from the new end.
    retrievedComponentTypes = std::reverse(retrievedComponentTypes.begin(), retrievedComponentTypes.end());
    std::vector<int> retrievedComponentStartNodes = netlistReader_instance->getComponentStartNodes().at(m_IndexOfThisNetlistLPN);
    retrievedComponentStartNodes = std::reverse(retrievedComponentStartNodes.begin(), retrievedComponentStartNodes.end());
    std::vector<int> retrievedComponentEndNodes = netlistReader_instance->getComponentEndNodes().at(m_IndexOfThisNetlistLPN);
    retrievedComponentEndNodes = std::reverse(retrievedComponentEndNodes.begin(), retrievedComponentEndNodes.end());
    std::vector<double> retrievedComponentParameterValues = netlistReader_instance->getComponentParameterValues().at(m_IndexOfThisNetlistLPN);
    retrievedComponentParameterValues = std::reverse(retrievedComponentParameterValues.begin(), retrievedComponentParameterValues.end());

    std::vector<int> retrievedListOfPrescribedPressures = netlistReader_instance->getListOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    std::vector<circuit_nodal_pressure_prescription_t> retrievedTypeOfPrescribedPressures = netlistReader_instance->getTypeOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);
    std::vector<double> retrievedValueOfPrescribedPressures = netlistReader_instance->getValueOfPrescribedPressures().at(m_IndexOfThisNetlistLPN);

    std::vector<int> retrievedListOfPrescribedFlows = netlistReader_instance->getListOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::vector<circuit_component_flow_prescription_t> retrievedTypeOfPrescribedFlows = netlistReader_instance->getTypeOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::vector<double> retrievedValueOfPrescribedFlows = netlistReader_instance->getValueOfPrescribedFlows().at(m_IndexOfThisNetlistLPN);
    std::map<int,double> retrievedInitialPressures = netlistReader_instance->getInitialPressures().at(m_IndexOfThisNetlistLPN);
    
    // Loop over the components, assigning them (and their nodes) the appropriate properties to give the fully-described circuit:
    for (int component = m_CircuitDescription.components.begin(); component != m_CircuitDescription.components.end(); component++)
    {
        component->type = retrievedComponentTypes.back();
        retrievedComponentTypes.pop_back();

        component->startNode.indexInCircuit = retrievedComponentStartNodes.back();
        retrievedComponentStartNodes.pop_back();
        component->startNode.prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        for (int prescribedPressure=0; prescribedPressure<m_CircuitDescription.numberOfPrescribedPressures; prescribedPressure++)
        {
            if (retrievedListOfPrescribedPressures.at(prescribedPressure) == component->startNode.indexInCircuit)
            {
                component->startNode.prescribedPressureType = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
                component->startNode.valueOfPrescribedPressure =retrievedValueOfPrescribedPressures.at(prescribedPressure);
            }
        }

        component->endNode.indexInCircuit = retrievedComponentEndNodes.back();
        retrievedComponentEndNodes.pop_back();
        component->endNode.prescribedPressureType = Pressure_NotPrescribed; // initialise as a default, before replacing as necessary
        for (int prescribedPressure=0; prescribedPressure<m_CircuitDescription.numberOfPrescribedPressures; prescribedPressure++)
        {
            if (retrievedListOfPrescribedPressures.at(prescribedPressure) == component->endNode.indexInCircuit)
            {
                component->endNode.prescribedPressureType = retrievedTypeOfPrescribedPressures.at(prescribedPressure);
                component->endNode.valueOfPrescribedPressure =retrievedValueOfPrescribedPressures.at(prescribedPressure);
            }
        }

        component->prescribedFlowType = Flow_NotPrescribed;  // initialise as a default, before replacing as necessary
        for (int prescribedFlow=0; prescribedFlow<m_CircuitDescription.numberOfPrescribedFlows; prescribedFlow++)
        {
            if (retrievedListOfPrescribedFlows.at(prescribedFlow) == component->indexInCircuit)
            {
                component->prescribedFlowType = retrievedTypeOfPrescribedFlows.at(prescribedFlow);
                component->valueOfPrescribedFlow = retrievedValueOfPrescribedFlows.at(prescribedFlow);
            }
        }

        component->parameterValue = retrievedComponentParameterValues.back();
        retrievedComponentParameterValues.pop_back();

        component->startNode.pressure = retrievedInitialPressures.find(component->startNode.indexInCircuit);
        component->endNode.pressure = retrievedInitialPressures.find(component->endNode.indexInCircuit);


    }

    // Some metadata is already set-up for this circuit, as a side-effect
    // of the above construction. This call completes the metadata; it's not
    // a problem that it also re-writes some of the existing metadata
    // (rewrites - but does not change - the values are identical!)
    m_CircuitDescription.rebuildCircuitMetadata();
    

    // // Component indices are just consecutive integers by default, but sometimes non-consecutive numbering
    // // is needed; componentIndices allows for this.
    // // We initialise it now for the default case.
    // for (int ii=1; ii < m_CircuitDescription.numberOfComponents + 1; ii++)
    // {
    //     m_CircuitDescription.componentIndices.push_back(ii);
    // }
}

void NetlistBoundaryCondition::createInitialCircuitDescriptionWithoutDiodes()
{
    // Copy the data as a base for modification (i.e. the removal of the diode data)
    m_CircuitDescriptionWithoutDiodes = m_CircuitDescription;
    // Prepare for diode removal by clearing the data that will be rebuilt:
    // m_CircuitDescriptionWithoutDiodes.deleteSubcircuitSpecificData();
    // Rebuild the cleared vectors, but without the diodes:
    int numberOfDiodes = 0;
    for (auto component=m_CircuitDescriptionWithoutDiodes.components.begin(); component!=m_CircuitDescriptionWithoutDiodes.components.end(); component++)
    {
        if(component->type == Component_Diode)
        {
            component = m_CircuitDescriptionWithoutDiodes.components.erase(component); // returns an iterator pointing to the new location of the element after the one that just got erased
            component--; // decrement the returned iterator, so the for loop incrememnts it back to the correct (next) element.
        }
    }

    m_CircuitDescriptionWithoutDiodes.rebuildCircuitMetadata();

}

void NetlistBoundaryCondition::assignComponentsToAtomicSubcircuits()
{
    // This subroutine builds m_AtomicSubcircuitsComponentsBelongsTo, which is indexed by component, as they appear in m_CircuitDescriptionWithoutDiodes.

    // Group the components by partitioning the node indices into connected subcircuits (with diodes removed).
    // We do this by parsing the circuit to discover the topology.
    std::vector<bool> componentAssignedToASubcircuit(m_CircuitDescriptionWithoutDiodes.numberOfComponents, false);
    m_AtomicSubcircuitsComponentsBelongsTo.insert(m_AtomicSubcircuitsComponentsBelongsTo.begin(), m_CircuitDescriptionWithoutDiodes.numberOfComponents, -1); // initialise with nonsense value
    int subcircuitIndex = 0;
    for (int startingComponent = 0; startingComponent < m_CircuitDescriptionWithoutDiodes.numberOfComponents; startingComponent++)
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
            nodesInSubcircuit.insert(m_CircuitDescriptionWithoutDiodes.components.at(startingComponent).startNode.indexInCircuit);
            nodesInSubcircuit.insert(m_CircuitDescriptionWithoutDiodes.components.at(startingComponent).endNode.indexInCircuit);
            // Find all other components belonging to this subcircuit. We must iterate (while-loop) to keep
            // finding new connections, until nothing changes on an iteration:
            bool newComponentFoundOnLastIteration = true;
            while (newComponentFoundOnLastIteration)
            {
                newComponentFoundOnLastIteration = false;
                for (int potentialSubcircuitComponent = startingComponent+1; potentialSubcircuitComponent < m_CircuitDescriptionWithoutDiodes.numberOfComponents; potentialSubcircuitComponent++)
                {
                    if (componentAssignedToASubcircuit.at(potentialSubcircuitComponent) = false)
                    {
                        // get the nodes of the current component
                        int componentStartNode = m_CircuitDescriptionWithoutDiodes.components.at(potentialSubcircuitComponent).startNode.indexInCircuit;
                        int componentEndNode = m_CircuitDescriptionWithoutDiodes.components.at(potentialSubcircuitComponent).endNode.indexInCircuit;

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
                            nodesInSubcircuit.insert(componentEndNode)

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
        m_CircuitDataForAtomicSubcircuits.push_back(new CircuitData);
        // Get a reference to the new CircuitData, just to avoid filling the screen with an infinite amount of code, below.
        CircuitData& currentSubcircuit = m_CircuitDataForAtomicSubcircuits.back();

        int componentLocationIdx = 0;
        for (auto componentsSubcircuitIdx=m_AtomicSubcircuitsComponentsBelongsTo.begin(); componentsSubcircuitIdx!=m_AtomicSubcircuitsComponentsBelongsTo.end(); componentsSubcircuitIdx++)
        {
            if (*componentsSubcircuitIdx == currentSubcircuitIdx)
            {
                currentSubcircuit.components.push_back(new CircuitComponent);
                currentSubcircuit.components.back() = m_CircuitDescriptionWithoutDiodes.components.at(componentLocationIdx);
            }
            componentLocationIdx++;
        }

        currentSubcircuit.rebuildCircuitMetadata();
    }
        

    //     // Initialise by just copying the data:
    //     // circDataToFill = m_CircuitDescription;
    //     // Delete the data that we don't want in the copy of the circuitDescription:
    //     // circDataToFill.deleteSubcircuitSpecificData();
    //     // Find the components which make up this atomic subcircuit, and build their Circuitdata:
    //     int componentIdxInNoDiodeData = 0;
    //     for (auto subcircuitMembershipIdx=m_AtomicSubcircuitsComponentsBelongsTo.begin(); subcircuitMembershipIdx!=m_AtomicSubcircuitsComponentsBelongsTo.end(); subcircuitMembershipIdx++)
    //     {
    //         if(*subcircuitMembershipIdx == currentSubcircuit)
    //         {
    //             circDataToFill.componentTypes.push_back(m_CircuitDescriptionWithoutDiodes.componentTypes.at(componentIdxInNoDiodeData));
    //             circDataToFill.componentStartNodes.push_back(m_CircuitDescriptionWithoutDiodes.componentStartNodes.at(componentIdxInNoDiodeData));
    //             circDataToFill.componentEndNodes.push_back(m_CircuitDescriptionWithoutDiodes.componentEndNodes.at(componentIdxInNoDiodeData));
    //             circDataToFill.componentParameterValues.push_back(m_CircuitDescriptionWithoutDiodes.componentParameterValues.at(componentIdxInNoDiodeData));
    //             circDataToFill.componentIndices.push_back(m_CircuitDescriptionWithoutDiodes.componentIndices.at(componentIdxInNoDiodeData));
                
    //             // Sets only allow one copy of each element, so this is a unique count
    //             circDataToFill.setOfPressureNodes.insert(circDataToFill.componentStartNodes.back());
    //             circDataToFill.setOfPressureNodes.insert(circDataToFill.componentEndNodes.back());
    //         }
    //         componentIdxInNoDiodeData++;
    //     }
    //     circDataToFill.numberOfComponents = circDataToFill.componentTypes.size();
    //     circDataToFill.numberOfPressureNodes = circDataToFill.setOfPressureNodes.size();
        
    //     // Copy the component indices into the convenience/utility data structure setOfComponentIndices.
    //     for (auto iterator=circDataToFill.componentStartNodes.begin(); iterator!=circDataToFill.componentStartNodes.end(); iterator++)
    //     {
    //         circDataToFill.setOfComponentIndices.insert(*iterator);
    //     }
    //     for (auto iterator=circDataToFill.componentEndNodes.begin(); iterator!=circDataToFill.componentEndNodes.end(); iterator++)
    //     {
    //         circDataToFill.setOfComponentIndices.insert(*iterator);
    //     }
    // }

    // // clean up the atomic subcircuit data, so that finally, only the necessary, relevant data is present.
    // for (auto subcircuit = m_CircuitDataForSubcircuits.begin(); subcircuit!=m_CircuitDataForSubcircuits.end; subcircuit++)
    // {
    //     // See if the prescribed pressure from the original circuit description belongs in this subcircuit:
    //     int typeOfPrescribedPressureLocationCounter = 0;
    //     for (auto prescribedPressure = m_CircuitDescriptionWithoutDiodes.listOfPrescribedPressures.begin(); prescribedPressure != m_CircuitDescriptionWithoutDiodes.listOfPrescribedPressures.end(); prescribedPressure++)
    //     {
    //         // If it does belong, add it to the listOfPrescribedPressures (add the appropriate typeOfPrescribedPressures to the appropriate vector, too)
    //         if (subcircuit->setOfPressureNodes.count(*prescribedPressure))
    //         {
    //             subcircuit->listOfPrescribedPressures.push_back(*prescribedPressure);
    //             subcircuit->typeOfPrescribedPressures.push_back(m_CircuitDescriptionWithoutDiodes.typeOfPrescribedPressures.at(typeOfPrescribedPressureLocationCounter));
    //         }

    //         typeOfPrescribedPressureLocationCounter++;
    //     }
    //     subcircuit->numberOfPrescribedPressures = listOfPrescribedPressures.size();

    //     // See if the prescribed flows from the original circuit description belongs in this subcircuit:
    //     int typeOfPrescribedFlowLocationCounter = 0;
    //     for (auto prescribedFlow = m_CircuitDescriptionWithoutDiodes.listOfPrescribedFlows.begin(); prescribedFlow!=m_CircuitDescriptionWithoutDiodes.listOfPrescribedFlows.end(); prescribedFlow++)
    //     {
    //         // if it does belong, add it to the listOfPrescribedFlows for the subcircuit (and add the appropriate typeOfPrescribedFlows to the appropriate vector, too)
    //         if (subcircuit->setOfComponentIndices.count(*prescribedFlow))
    //         {
    //             subcircuit->listOfPrescribedFlows.push_back(*prescribedFlow);
    //             subcircuit->typeOfPrescribedFlows.push_back(m_CircuitDescriptionWithoutDiodes.typeOfPrescribedFlows.at(typeOfPrescribedFlowLocationCounter));
    //         }
    //         typeOfPrescribedFlowLocationCounter++;
    //     }
    //     subcircuit->numberOfPrescribedFlows = listOfPrescribedFlows.size();
    // }

    // Identify the composite subcircuits (with legality assigned by user?)

    // Make the data for the composite subcircuits

    // Make sure that we have a sensible data structure so that the list of active subcircuits
    // can easily be mapped to the actual active circuits themselves.
}

void NetlistBoundaryCondition::createCompositeSubcircuitDescriptions()
{
    // Find the diodes in the origina input data:
    //- we index the diodes in the order that they appeared in the original netlist input file.
    std::map<int,
    for (component = m_CircuitDescription.components.begin(); component!=m_CircuitDescription.components.end(); component++)
    {
        if (*component->type == Component_Diode)
        {
            if ()
        }
    }

    // Identify all composite subcircuits by finding all possible connections between atomic subcircuits
    //- There are 2^(number of diodes) possible states. Index the states using integers,
    //- such that the binary representation of the index gives the pattern of diode open/closed states.
    //- the least significant bit will be the first valve, so with 5 valves, 00101 would
    //- indicate that valves 0 and 2 are open.
    //-
    //- we index the diodes in the order that they appeared in the original netlist input file.

    for (int diodeOpenPatternIdx=0; diodeOpenPatternIdx<pow(2,m_numberOfDiodes); diodeOpenPatternIdx++)
    {
        std::bitset<m_numberOfDiodes> diodeOpenPattern(diodeOpenPatternIdx);
        
    }

}

void NetlistBoundaryCondition::buildDiodeMetadata()
{
    // Count the diodes, get a map to them:
    m_numberOfDiodes = 0;
    for (auto component = m_CircuitDescription.components.begin(); component != m_CircuitDescription.components.end(); component++)
    {
        if (*component->type == Component_Diode)
        {
            m_diodeIndexingMap.insert(std::pair<int,boost::shared_ptr<CircuitComponent>> (m_numberOfDiodes, *component))
            m_numberOfDiodes++;
        }
    }
    
    m_atomicSubcircuitConnectionManager = new AtomicSubcircuitConnectionManager(m_numberOfDiodes);
    // Get a map from diode index to the two atomic subcircuits that it joins:
    for (int diodeIdx = 0; diodeIdx < m_numberOfDiodes; diodeIdx++)
    {
        boost::shared_ptr<CircuitComponent> currentDiode = m_diodeIndexingMap.find(diodeIdx);

        // Prepare a place to temporarily package the two adjoining circuits for this diode:
        std::pair<boost:shared_ptr<CircuitData>, boost::shared_ptr<CircuitData>> adjoiningAtomicSubcircuitsPair;

        // Loop the atomic subcircuits, looking for the subcircuit that the diode connects to:
        for (auto atomicSubcircuit = m_CircuitDataForAtomicSubcircuits.begin(); atomicSubcircuit != m_CircuitDataForAtomicSubcircuits.end(); atomicSubcircuit++)
        {
            if (*atomicSubcircuit->setOfPressureNodes.count(currentDiode->startNode) == 1)
            {
                adjoiningAtomicSubcircuitsPair.first = *atomicSubcircuit;
            }

            if (*atomicSubcircuit->setOfPressureNodes.count(currentDiode->endNode) == 1)
            {
                adjoiningAtomicSubcircuitsPair.second = *atomicSubcircuit;
            }
        }

        // Note that in the case where a diode is "hanging" (i.e. the diode has no further components on one side of it)
        // the relevant entry of adjoiningAtomicSubcircuitsPair will be NULL. Careful with this! The AtomicSubcircuitConnectionManager
        // is designed to help avoid making mistakes here.
        circuit_diode_node_t startNodeSubcircuitType;
        if (adjoiningAtomicSubcircuitsPair.first == NULL)
        {
            startNodeSubcircuitType = Node_IsMonopolar;
        }
        else
        {
            startNodeSubcircuitType = Node_ConnectsCircuit;
        }

        circuit_diode_node_t endNodeSubcircuitType;
        if (adjoiningAtomicSubcircuitsPair.second == NULL)
        {
            endNodeSubcircuitType = Node_IsMonopolar;
        }
        else
        {
            endNodeSubcircuitType = Node_ConnectsCircuit;
        }

        m_atomicSubcircuitConnectionManager->setDiodeInfo(diodeIdx,startNodeSubcircuitType,endNodeSubcircuitType,adjoiningAtomicSubcircuitsPair.first,adjoiningAtomicSubcircuitsPair.second);

    }
}