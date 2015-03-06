#include "Netlist3DDomainReplacement.hxx"
#include "fileReaders.hxx"

void Netlist3DDomainReplacement::createCircuitDescription()
{
    // This function creates the internal CircuitData class format for a zero-D
    // Netlist, non-boundary-condition replacement for the 3D domain, for pure zero-D simulation.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    netlistReader* netlistReader_instance = netlistReader::Instance();
    m_numberOfNetlistsUsedAsBoundaryConditions = netlistReader_instance->getNumberOfNetlistSurfaces();

    // Make the appropriate class to store the 3D domain replacement circuit data:
    mp_CircuitDescription = boost::shared_ptr<Netlist3DDomainReplacementCircuitData> (new Netlist3DDomainReplacementCircuitData(hstep,m_numberOfNetlistsUsedAsBoundaryConditions));

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
    std::reverse(componentTypes.begin(), componentTypes.end());
    // Create the component start node info (all are the same: we make it the highest-indexed node.)
    
    // int centralNodeIndex = mp_CircuitDescription->numberOfComponents + 1;
    // std::vector<int> componentStartNodes(mp_CircuitDescription->numberOfComponents, centralNodeIndex);

    std::vector<int> componentStartNodes;
    componentStartNodes.push_back(6);
    componentStartNodes.push_back(7);
    componentStartNodes.push_back(8);
    componentStartNodes.push_back(5);
    componentStartNodes.push_back(5);
    componentStartNodes.push_back(5);
    componentStartNodes.push_back(5);

    std::reverse(componentStartNodes.begin(), componentStartNodes.end()); // actually no point in this call, but it's tidier to leave it here for symmetry with the related calls below.

    // Deal out indices for the end nodes:
    std::vector<int> componentEndNodes;
    // for (int component=0; component<mp_CircuitDescription->numberOfComponents; component++)
    // {
    // 	componentEndNodes.push_back(component+1);
    // }

    componentEndNodes.push_back(1);
    componentEndNodes.push_back(2);
    componentEndNodes.push_back(3);
    componentEndNodes.push_back(6);
    componentEndNodes.push_back(7);
    componentEndNodes.push_back(8);
    componentEndNodes.push_back(4);

    std::reverse(componentEndNodes.begin(), componentEndNodes.end());

    // std::vector<double> componentParameterValues(m_numberOfNetlistsUsedAsBoundaryConditions, m_oneResistanceToGiveEachResistor);
    // componentParameterValues.push_back(m_elastanceToGiveVolumeTrackingPressureChamber);
    
    std::vector<double> componentParameterValues;
    componentParameterValues.push_back(NAN);
    componentParameterValues.push_back(NAN);
    componentParameterValues.push_back(NAN);
    componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
    componentParameterValues.push_back(m_oneResistanceToGiveEachResistor);
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
    for (int node=1; node < 4; node++)
    {
    	initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }
    // Do the node at the base of the pressure chamber:
    initialPressures.insert(std::make_pair(4, 0.0));
    for (int node=5; node <= mp_CircuitDescription->numberOfPressureNodes; node++)
    {
        initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }

    // Do the final node (at the base of the pressure chamber):
    // initialPressures.insert(std::make_pair(mp_CircuitDescription->numberOfPressureNodes, 0.0));
    
    // Loop over the components, assigning them (and their nodes) the appropriate properties to give the fully-described circuit:
    for (auto component = mp_CircuitDescription->components.begin(); component != mp_CircuitDescription->components.end(); component++)
    {
        (*component)->type = componentTypes.back();
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
        //     if (listOfPrescribedFlows.at(prescribedFlow) == (*component)->indexInInputData)
        //     {
        //         (*component)->prescribedFlowType = typeOfPrescribedFlows.at(prescribedFlow);
        //         (*component)->valueOfPrescribedFlow = valueOfPrescribedFlows.at(prescribedFlow);
        //     }
        // }

        (*component)->currentParameterValue = componentParameterValues.back();
        // make a copy of this value so we can reset it if necessary. Used for e.g. diode state changes.
        (*component)->parameterValueFromInputData = (*component)->currentParameterValue;
        componentParameterValues.pop_back();

        (*component)->startNode->setPressure(initialPressures.at((*component)->startNode->indexInInputData));
        (*component)->endNode->setPressure(initialPressures.at((*component)->endNode->indexInInputData));

        if ((*component)->indexInInputData <= m_numberOfNetlistsUsedAsBoundaryConditions)
        {
            (*component)->prescribedFlowPointerIndex = (*component)->indexInInputData-1;
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
void Netlist3DDomainReplacement::setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode)
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
    // node->indexInInputData = indexOfNodeInInputData;
    // node->prescribedPressureType = typeOfPrescribedPressure;
    if (node->prescribedPressureType!=Pressure_NotPrescribed && node->prescribedPressureType!=Pressure_Null)
    {
        node->setPressure(valueOfPrescribedPressures.at(indexOfPrescribedPressure));
    }

    if (node->indexInInputData <= m_numberOfNetlistsUsedAsBoundaryConditions)
    {
        node->prescribedPressurePointerIndex = node->indexInInputData-1;
    }
}

void Netlist3DDomainReplacement::setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate)
{
	m_boundaryFlowsOrPressuresAsAppropriate = boundaryFlowsOrPressuresAsAppropriate;
	boost::shared_ptr<Netlist3DDomainReplacementCircuitData> downcast3DDomainReplacementCircuitData = boost::static_pointer_cast<Netlist3DDomainReplacementCircuitData> (mp_CircuitDescription);
	downcast3DDomainReplacementCircuitData->setBoundaryPrescriptionsAndBoundaryConditionTypes(m_boundaryFlowsOrPressuresAsAppropriate);
}

std::vector<double> Netlist3DDomainReplacement::getBoundaryPressures()
{
	std::vector<double> pressures;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		pressures.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->startNode->getPressure());
	}
	return pressures;
}

std::vector<double> Netlist3DDomainReplacement::getBoundaryFlows()
{
	std::vector<double> flows;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		flows.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->flow);
	}
	return flows;
}

void Netlist3DDomainReplacement::solveSystem(const int timestepNumber)
{
	m_activeSubcircuits.at(0)->buildAndSolveLinearSystem(timestepNumber);
}

void Netlist3DDomainReplacement::initialiseModel()
{
    // Get the input data
    createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    identifyAtomicSubcircuits();

    // chop up the Circuitdata into subcircuits (including removing all mention of diodes from that data now)
    createAtomicSubcircuitDescriptions();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    // mp_CircuitDescription->closeAllDiodes();
    // mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    
    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    // m_atomicSubcircuitConnectionManager = boost::shared_ptr<AtomicSubcircuitConnectionManager>( toPassToSharedPtr );
}

void Netlist3DDomainReplacement::setPointersToBoundaryPressuresAndFlows(double* const interfacePressuresToBeReadBy3DDomainReplacement, double* const interfaceFlowsToBeReadBy3DDomainReplacement, const int& numberOfPointers)
{
    assert(pressure_n_ptrs.size()==0);
    assert(flow_n_ptrs.size()==0);
    for (int pointerIndex = 0; pointerIndex<numberOfPointers; pointerIndex++)
    {
        pressure_n_ptrs.push_back(&interfacePressuresToBeReadBy3DDomainReplacement[pointerIndex]);
        flow_n_ptrs.push_back(&interfaceFlowsToBeReadBy3DDomainReplacement[pointerIndex]);
    }
}

// void Netlist3DDomainReplacement::setPointersToBoundaryPressuresAndFlows(double* const mp_interfacePressures,double* const mp_interfaceFlows)
// {

// }

// std::vector<double> Netlist3DDomainReplacement::getPressureOrFlowPrescriptionsToReturnToNetlistBoundaryConditions()
// {
// 	std::vector<double> pressuresAndFlowsToGiveToBoundaryConditions;
// 	// we know the components at the boundary are just the first m_numberOfNetlistsUsedAsBoundaryConditions in the CircuitData. Loop over them!
// 	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
// 	{
// 		// Determine whether it's a Dirichlet or Neumann boundary condition at the interface (so pass a pressure or a flow to the boundary condition, respectively).
// 		if (mp_CircuitDescription->hasPrescribedPressureAcrossInterface())
// 		{
// 			// Gather the pressure or flow value into the return variable:
// 			pressuresAndFlowsToGiveToBoundaryConditions.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->flow);
// 		}
// 		else if (mp_CircuitDescription->hasPrescribedFlowAcrossInterface())
// 		{
// 			// Gather the pressure or flow value into the return variable:
// 			pressuresAndFlowsToGiveToBoundaryConditions.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->endNode->pressure);
// 		}
// 		else
// 		{
// 			std::stringstream errorMessage;
// 			errorMessage << "EE: Zero-D domain replacement has at least one bounday interface with neither a prescribed pressure nor a prescribed flow." << std::endl;
// 			throw std::logic_error(errorMessage);			
// 		}
// 	}
// 	return pressuresAndFlowsToGiveToBoundaryConditions;
// }

void Netlist3DDomainReplacement::selectAndBuildActiveSubcircuits()
{
    mp_CircuitDescription->rebuildCircuitMetadata();
    mp_CircuitDescription->generateNodeAndComponentIndicesLocalToSubcircuit();
    // mp_CircuitDescription.switchDiodeStatesIfNecessary();
    // mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    // Actually build the active NetlistSubcircuit classes:
    m_activeSubcircuits.clear();
    double alfi_delt_in = alfi_local*delt;
    boost::shared_ptr<NetlistSubcircuit> newNetlistZeroDDomain(new NetlistSubcircuit(0, mp_CircuitDescription, flow_n_ptrs, pressure_n_ptrs,alfi_delt_in,surfaceIndex));
    newNetlistZeroDDomain->setThisisA3DDomainReplacement();
    m_activeSubcircuits.push_back(newNetlistZeroDDomain);

}

void Netlist3DDomainReplacement::setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients)
{
    for (auto component=mp_CircuitDescription->components.begin(); component!=mp_CircuitDescription->components.end(); component++)
    {
        if ((*component)->indexInInputData <= 3)
        {
            // the component indexInInputData needs to be converted to zero-indexing:
            double potentialResistance = allImplicitCoefficients.at((*component)->indexInInputData-1).first;

            if (potentialResistance == 0.0)
            {
                potentialResistance = 0.000001;
            }
            (*component)->currentParameterValue = potentialResistance;
        }
    }
}

// void Netlist3DDomainReplacement::switchSurfaceToZeroFlow(const int& surfaceIndex)
// {

// }