#include "NetlistZeroDDomainCircuit.hxx"
#include "indexShifters.hxx"

// NetlistZeroDDomainCircuit is the class which represents the LPN(s) which replaces the 3D domain(s) in pure zeroD mode.

void NetlistZeroDDomainCircuit::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    rebuildCircuitMetadata();
    // cycleToSetHistoryPressuresFlowsAndVolumes();
}

// This function takes a resistor which connects directly to one of the compliance units, and works
// out which connected component it belongs to (so it can be attached to the appropriate compliance unit)
int NetlistZeroDDomainCircuit::getConnectedTopologicalComponentIndexForInnerResistor(const int componentIndex) const
{
    assert(componentIndex > m_numberOfOutlets);
    assert(componentIndex <= 2*m_numberOfOutlets);
    int indexOfThisZeroDSurface = toZeroIndexing(componentIndex - m_numberOfOutlets);
    std::cout << "requested " << indexOfThisZeroDSurface << std::endl;
    for (auto& pair : m_mapFromZeroDSurfaceIndexToConnectedComponentIndex)
    {
        std::cout << "pair entry: " << pair.first << ", " << pair.second << std::endl;
    }
    return m_mapFromZeroDSurfaceIndexToConnectedComponentIndex.at(indexOfThisZeroDSurface);
}

void NetlistZeroDDomainCircuit::findNumberOfConnectedComponentsOf3DDomain()
{
    int maxConnectedComponentIndex = 0;
    for (auto netlistBoundaryConditionData = m_mapFromZeroDSurfaceIndexToConnectedComponentIndex.begin(); netlistBoundaryConditionData != m_mapFromZeroDSurfaceIndexToConnectedComponentIndex.end(); netlistBoundaryConditionData++)
    {
        if (netlistBoundaryConditionData->second > maxConnectedComponentIndex)
        {
            maxConnectedComponentIndex = netlistBoundaryConditionData->second;
        }
    }

    m_numberOfConnectedComponentsOf3DDomain = maxConnectedComponentIndex;
}

void NetlistZeroDDomainCircuit::createCircuitDescription()
{
    // This function creates the internal CircuitData class format for a zero-D
    // Netlist, non-boundary-condition replacement for the 3D domain, for pure zero-D simulation.

    const int totalNumberOfResistorsIn3DDomainReplacementCircuit = 2 * (m_numberOfOutlets);
    // we'll have a resistor for each netlist used as a boundary condition, plus one VolumeTrackingComponent for each connected component:
    mp_circuitData->numberOfComponents = totalNumberOfResistorsIn3DDomainReplacementCircuit + m_numberOfConnectedComponentsOf3DDomain;

    // The components will be arranged in a star (all start-nodes are connected together at a single point, all end-nodes only connect to one component)
    int numberOfAdditionalNodesDueToMultipleConnectedComponentsExisting = 2*(m_numberOfConnectedComponentsOf3DDomain - 1); // -1 because if there's only one connected component, there are no additional nodes created by splitting the 3D domain replacement circuit into disjoint pieces
    mp_circuitData->numberOfPressureNodes = totalNumberOfResistorsIn3DDomainReplacementCircuit + 2 + numberOfAdditionalNodesDueToMultipleConnectedComponentsExisting;

    // The pressures will initially be received from the boundary conditions at the boundary condition interface nodes,
    // plus the zero-pressure prescription on the base of the VolumeTrackingComponent (one for each connected component:
    mp_circuitData->numberOfPrescribedPressures = m_numberOfConnectedComponentsOf3DDomain;
    // None initially:
    mp_circuitData->numberOfPrescribedFlows = m_numberOfOutlets;

    // Create the component data:
    std::vector<circuit_component_t> componentTypes;
    for (int boundary=0; boundary < m_numberOfOutlets; boundary++)
    {
    	componentTypes.push_back(Component_Resistor);
        componentTypes.push_back(Component_Resistor);
    }

    for (int connectedComponentIndex = 0; connectedComponentIndex < m_numberOfConnectedComponentsOf3DDomain; connectedComponentIndex++)
    {
        componentTypes.push_back(Component_Capacitor);
    }

    // Prepare space for the components in the circuit:
    assert(mp_circuitData->components.empty());
    for (int ii=0; ii<mp_circuitData->numberOfComponents; ii++)
    {
        CircuitComponent* toPushBack;
        if (componentTypes.at(ii) == Component_VolumeTrackingPressureChamber)
        {
            double initialVolume = 130000.0; //\todo make adjustable
            double initialUnstressedVolume = 0.0; //\todo make adjustable
            toPushBack = new VolumeTrackingPressureChamber(m_hstep, m_thisIsARestartedSimulation, initialVolume, initialUnstressedVolume);
        }
        else if (componentTypes.at(ii) == Component_VolumeTracking)
        {
            double initialVolume = 130000.0; //\todo make adjustable
            toPushBack = new VolumeTrackingComponent(m_hstep, m_thisIsARestartedSimulation, initialVolume);
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
    //
    // As some of the outlets may not be Netlists (so far, the other possibility is BCT-type flow prescription)
    // we just assume that the netlists come first in the indexing, and the other types later.
    for (int dpDqResistorIndex = 1; dpDqResistorIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; dpDqResistorIndex++)
    {
        componentEndNodes.push_back(dpDqResistorIndex);
    }

    // Now do the non-netlist surface end-nodes (this could be done as part of the dpDqResistorIndex loop above),
    // but we split it here for human clarity.
    for (int prescribedFlowSurfaceResistorIndex = m_numberOfNetlistsUsedAsBoundaryConditions + 1; prescribedFlowSurfaceResistorIndex <= m_numberOfOutlets; prescribedFlowSurfaceResistorIndex++)
    {
        componentEndNodes.push_back(prescribedFlowSurfaceResistorIndex);
    }


    // The next 2*m_numberOfConnectedComponentsOf3DDomain node indices are used for the domain central node of each connected component
    // (it's basically star-shaped in each connected component), and the end-node singleton for the bottom of each compliance chamber.
    //
    // Skip these, and now do the end nodes for the resistors which represent
    // the resistance of the 3D domain itself.
    for (int internalDomainResistorEndNodeIndex = m_numberOfOutlets + 1 + 2*m_numberOfConnectedComponentsOf3DDomain; internalDomainResistorEndNodeIndex <= mp_circuitData->numberOfPressureNodes; internalDomainResistorEndNodeIndex++)
    {
        componentEndNodes.push_back(internalDomainResistorEndNodeIndex);
    }
    // Finally, do the end-node at the base of the compliance chambers:
    for (int complianceUnitIndex = totalNumberOfResistorsIn3DDomainReplacementCircuit+1; complianceUnitIndex <= mp_circuitData->numberOfComponents; complianceUnitIndex++)
    {
        int connectedComponentIndexOfThisComplianceUnit = complianceUnitIndex - totalNumberOfResistorsIn3DDomainReplacementCircuit;
        int endNodeIndexForThisComplianceUnit = 2*connectedComponentIndexOfThisComplianceUnit + m_numberOfOutlets - 1;
        componentEndNodes.push_back(endNodeIndexForThisComplianceUnit);
    }

    std::vector<int> componentStartNodes;
    // All the components, except for the dP/dQ resistors (or other BC type edge resistors), have the same
    // start-node (the centre-point of the star-shape of this domain).
    // We do the dP/dQ (and other edge-type) resistors first:
    //
    // the "+2*m_numberOfConnectedComponentsOf3D domain" is to make space for a compliance
    // chamber for each connected component
    for (int dpDqResistorStartNodeIndex = m_numberOfOutlets + 1 + 2*m_numberOfConnectedComponentsOf3DDomain; dpDqResistorStartNodeIndex <= mp_circuitData->numberOfPressureNodes; dpDqResistorStartNodeIndex++)
    {
        componentStartNodes.push_back(dpDqResistorStartNodeIndex);
    }
    // All the remaining components (compliance units & inner resistors) share the same start node, at the centre of the star-shaped domain (if there's
    // only one connected topological component of the 3D domain).
    //
    // If there's multiple topological connected components, all the remaining components /which are part of the same connected component/ have
    // the same start node. We set these indices now, as appropriate:
    //
    // resistors first....:
    for (int resistorIndex = m_numberOfOutlets+1; resistorIndex <= totalNumberOfResistorsIn3DDomainReplacementCircuit; resistorIndex++)
    {
        // componentStartNodes.push_back(m_numberOfNetlistsUsedAsBoundaryConditions+2);
        int connectedTopologicalComponentOfThisInnerResistor = getConnectedTopologicalComponentIndexForInnerResistor(resistorIndex);
        int startNodeIndexForThisInnerResistor = 2*connectedTopologicalComponentOfThisInnerResistor + m_numberOfOutlets;
        componentStartNodes.push_back(startNodeIndexForThisInnerResistor);
    }
    // ... then the compliance units:
    for (int complianceUnitIndex = totalNumberOfResistorsIn3DDomainReplacementCircuit+1; complianceUnitIndex <= mp_circuitData->numberOfComponents; complianceUnitIndex++)
    {
        int connectedComponentIndexOfThisComplianceUnit = complianceUnitIndex - totalNumberOfResistorsIn3DDomainReplacementCircuit;
        int startNodeIndexForThisComplianceUnit = 2*connectedComponentIndexOfThisComplianceUnit + m_numberOfOutlets;
        componentStartNodes.push_back(startNodeIndexForThisComplianceUnit);
    }

    std::reverse(componentStartNodes.begin(), componentStartNodes.end()); // actually no point in this call, but it's tidier to leave it here for symmetry with the related calls below.

    std::reverse(componentEndNodes.begin(), componentEndNodes.end());

    // std::vector<double> componentParameterValues(m_numberOfNetlistsUsedAsBoundaryConditions, m_oneResistanceToGiveEachResistor);
    // componentParameterValues.push_back(m_elastanceToGiveCentralCapacitor);

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
    for (int dpDqResistorIndex = 1; dpDqResistorIndex <= m_numberOfOutlets; dpDqResistorIndex++)
    {
        componentParameterValues.push_back(initialResistanceForDpDqResistors);
    }
    for (int internalDomainResistorIndex = 1; internalDomainResistorIndex <= m_numberOfOutlets; internalDomainResistorIndex++)
    {
        componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    }

    for (int complianceUnit = 0; complianceUnit < m_numberOfConnectedComponentsOf3DDomain; complianceUnit++)
    {
        componentParameterValues.push_back(m_elastanceToGiveCentralCapacitor);
    }

    std::reverse(componentParameterValues.begin(), componentParameterValues.end());

    std::vector<circuit_component_flow_prescription_t> typeOfPrescribedFlows;
    std::vector<double> valueOfPrescribedFlows;
    // std::vector<int> listOfPrescribedFlows;
    for (int interfaceComponentIndex=0; interfaceComponentIndex < m_numberOfNetlistsUsedAsBoundaryConditions; interfaceComponentIndex++)
    {
        typeOfPrescribedFlows.push_back(Flow_3DInterface);
        valueOfPrescribedFlows.push_back(0.0);
        // listOfPrescribedFlows.push_back(toOneIndexing(interfaceComponentIndex));
    }

    // Now do flow prescriptions for the prescribed flow type boundaries (these have no LPN BCM upstream of them, as Netlist BCs do)
    for (int interfaceComponentIndex = m_numberOfNetlistsUsedAsBoundaryConditions; interfaceComponentIndex < m_numberOfOutlets; interfaceComponentIndex++)
    {
        typeOfPrescribedFlows.push_back(Flow_Fixed);
        valueOfPrescribedFlows.push_back(1234.0);
        // listOfPrescribedFlows.push_back(toOneIndexing(interfaceComponentIndex));
    }
    
    std::map<int,double> initialPressures;
    // Do the nodes at the outside points of the domain which belong to resistors:
    for (int node=1; node <= m_numberOfOutlets; node++)
    {
    	initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }
    
    // Do the nodes at the bases of the pressure chambers:
    for (int complianceUnitIndex = 1; complianceUnitIndex <= m_numberOfConnectedComponentsOf3DDomain; complianceUnitIndex++)
    {
        int connectedComponentIndexOfThisComplianceUnit = complianceUnitIndex;
        int endNodeIndexForThisComplianceUnit = 2*connectedComponentIndexOfThisComplianceUnit + m_numberOfOutlets - 1;
        initialPressures.insert(std::make_pair(endNodeIndexForThisComplianceUnit, 0.0));
    }

    // Do the initial pressures at the centres of the star (top nodes of the compliance units):
    {
        int nodeIndex = m_numberOfOutlets + 2;
        for (int loopIndex = 0; loopIndex < m_numberOfConnectedComponentsOf3DDomain; loopIndex++)
        {
            initialPressures.insert(std::make_pair(nodeIndex,m_initialDomainPressure));
            nodeIndex += 2;
        }
    }

    // Do the initial pressures at the points between the pairs of resistors
    for (int nodeIndex = m_numberOfOutlets + 2 * m_numberOfConnectedComponentsOf3DDomain; nodeIndex <= mp_circuitData->numberOfPressureNodes; nodeIndex++)
    {
        initialPressures.insert(std::make_pair(nodeIndex, m_initialDomainPressure));
    }

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


        bool componentIsOuterResistorOnDirichletBranch = (*component)->getIndex() > m_numberOfNetlistsUsedAsBoundaryConditions && (*component)->getIndex() <= m_numberOfOutlets;
        if (componentIsOuterResistorOnDirichletBranch)
        {
            // Dirichlet / BCT surfaces have fixed flows:
            (*component)->prescribedFlowType = Flow_Fixed;
            (*component)->setPrescribedFlow(3.14159265); // recognisable nonsense initial value - in case it doesn't get properly replaced this value will be noticed in output data
            (*component)->addPythonControllerName(Controller_CustomPythonComponentFlowFile, "bctController");
        }
        else
        {
            (*component)->prescribedFlowType = Flow_NotPrescribed;
        }
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

    // Tell the node at the "3D" interface that it connects to the LPN domain:
    {
        std::vector<int> threeDNodeIndices;
        for (int nodeIndex = 1; nodeIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; nodeIndex++)
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
    // Do the Mop nodes (nodes at the edge of the NetlistZeroDDomainCircuit which attach to a dp/dq resistor)
    for (int pressureNodeIndex=1; pressureNodeIndex <= m_numberOfNetlistsUsedAsBoundaryConditions; pressureNodeIndex++)
    {
    	listOfPrescribedPressures.push_back(pressureNodeIndex);
    }

    // Now add the nodes at the base of the compliance units:
    {
        int highestIndexAmongstComplianceUnitBaseNodes = m_numberOfNetlistsUsedAsBoundaryConditions + 1 + 2*(m_numberOfConnectedComponentsOf3DDomain-1);
        for (int pressureNodeIndex=m_numberOfNetlistsUsedAsBoundaryConditions + 1; pressureNodeIndex <= highestIndexAmongstComplianceUnitBaseNodes; pressureNodeIndex+=2)
        {
            listOfPrescribedPressures.push_back(pressureNodeIndex);
        }
    }

    // All but one of the prescribed pressures are at the boundary interfaces:
    std::vector<circuit_nodal_pressure_prescription_t> typeOfPrescribedPressures;
    for (int pressureNodeIndex=0; pressureNodeIndex < m_numberOfNetlistsUsedAsBoundaryConditions; pressureNodeIndex++)
    {
    	typeOfPrescribedPressures.push_back(Pressure_3DInterface);
    }
    // The last prescribed pressures are all the ones at the bases of the VolumeTrackingComponent (or "compliance units" - they're now not always VolumeTrackingComponents)
    for (int complianceUnitIndex=0; complianceUnitIndex < m_numberOfConnectedComponentsOf3DDomain; complianceUnitIndex++)
    {
        typeOfPrescribedPressures.push_back(Pressure_Fixed);
    }

    std::vector<double> valueOfPrescribedPressures;
    // All but one of the prescribed pressures are at the boundary interfaces.
    for (int pressureNodeIndex=0; pressureNodeIndex < m_numberOfNetlistsUsedAsBoundaryConditions; pressureNodeIndex++)
    {
        //\todo remove \hardcoded to millimetres! units: 10^0 x Pa here
    	valueOfPrescribedPressures.push_back(1332.0); // ~10 mmHg
    }
    // The last prescribed pressures are all the ones at the bases of the VolumeTrackingComponent (or "compliance units" - they're now not always VolumeTrackingComponents)
    for (int complianceUnitIndex=0; complianceUnitIndex < m_numberOfConnectedComponentsOf3DDomain; complianceUnitIndex++)
    {
        valueOfPrescribedPressures.push_back(0.0);
    }


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

    if (node->getPressurePrescriptionType() != Pressure_NotPrescribed && node->getPressurePrescriptionType() != Pressure_Null)
    {
        node->setPrescribedPressure(valueOfPrescribedPressures.at(indexOfPrescribedPressure));
    }

    if (node->getIndex() <= m_numberOfNetlistsUsedAsBoundaryConditions)
    {
        node->setPrescribedPressurePointerIndex(node->getIndex()-1);
    }
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
                double potentialResistance = allImplicitCoefficients.at(toZeroIndexing((*component)->getIndex())).first;
                assert(!isnan(potentialResistance));

                (*component)->setParameterValue(potentialResistance);
                (*component)->setHasNoPrescribedFlow(); // defensive
            }

            bool boundaryCurrentlyHasPrescribedFlow = (pressuresOrFlowsAtBoundaries.at(dpDqResistorIndex).first == Boundary_Flow);
            if (boundaryCurrentlyHasPrescribedFlow)
            {
                double prescribedFlow = pressuresOrFlowsAtBoundaries.at(dpDqResistorIndex).second;
                (*component)->setPrescribedFlow(prescribedFlow);
            }

            dpDqResistorIndex++;
        }
    }
}

boost::shared_ptr<std::vector<std::pair<parameter_controller_t, int>>> NetlistZeroDDomainCircuit::getControlTypesAndComponentIndices() const
{
    return mp_circuitData->getControlTypesAndComponentIndices();
}