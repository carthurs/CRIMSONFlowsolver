#include "Netlist3DDomainReplacement.hxx"
#include "fileReaders.hxx"

void Netlist3DDomainReplacement::createCircuitDescription()
{
    // This function creates the internal CircuitData class format for a zero-D
    // Netlist, non-boundary-condition replacement for the 3D domain, for pure zero-D simulation.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    netlistReader* netlistReader_instance = netlistReader::Instance();
    m_numberOfNetlistsUsedAsBoundaryConditions = netlistReader_instance->getNumberOfNetlistSurfaces();
    // we'll have a resistor for each netlist used as a boundary condition, plus one VolumeTrackingPressureChamber:
    mp_CircuitDescription->numberOfComponents = m_numberOfNetlistsUsedAsBoundaryConditions + 1;

    // The components will be arranged in a star (all start-nodes are connected together at a single point, all end-nodes only connect to one component)
    mp_CircuitDescription->numberOfPressureNodes = m_numberOfNetlistsUsedAsBoundaryConditions + 2;

    // The pressures will initially be received from the boundary conditions at the boundary condition interface nodes:
    mp_CircuitDescription->numberOfPrescribedPressures = m_numberOfNetlistsUsedAsBoundaryConditions;
    // None initially:
    mp_CircuitDescription->numberOfPrescribedFlows = 0;

    // Create the component data:
    std::vector<circuit_component_t> componentTypes;
    for (int boundary=0; boundary<m_numberOfNetlistsUsedAsBoundaryConditions; boundary++)
    {
    	componentTypes.push_back(Component_Resistor);
    }
    componentTypes.push_back(Component_VolumeTrackingPressureChamber);

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
    int centralNodeIndex = mp_CircuitDescription->numberOfComponents + 1;
    std::vector<int> componentStartNodes(mp_CircuitDescription->numberOfComponents, centralNodeIndex);
    std::reverse(componentStartNodes.begin(), componentStartNodes.end()); // actually no point in this call, but it's tidier to leave it here for symmetry with the related calls below.

    // Deal out indices for the end nodes:
    std::vector<int> componentEndNodes;
    for (int component=0; component<mp_CircuitDescription->numberOfComponents; component++)
    {
    	componentEndNodes.push_back(component+1);
    }
    std::reverse(componentEndNodes.begin(), componentEndNodes.end());

    std::vector<double> componentParameterValues(m_numberOfNetlistsUsedAsBoundaryConditions, m_oneResistanceToGiveEachResistor);
    componentParameterValues.push_back(m_elastanceToGiveVolumeTrackingPressureChamber);
    std::reverse(componentParameterValues.begin(), componentParameterValues.end());

    std::vector<int> listOfPrescribedFlows; // empty initially

    std::vector<circuit_component_flow_prescription_t> typeOfPrescribedFlows; // initially empty
    std::vector<double> valueOfPrescribedFlows; // empty
    std::map<int,double> initialPressures;
    for (int node=1; node < mp_CircuitDescription->numberOfPressureNodes; node++)
    {
    	initialPressures.insert(std::make_pair(node,m_initialDomainPressure));
    }
    // Do the final node (at the base of the pressure chamber):
    initialPressures.insert(std::make_pair(mp_CircuitDescription->numberOfPressureNodes, 0.0));
    
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
        for (int prescribedFlow=0; prescribedFlow<mp_CircuitDescription->numberOfPrescribedFlows; prescribedFlow++)
        {
            if (listOfPrescribedFlows.at(prescribedFlow) == (*component)->indexInInputData)
            {
                (*component)->prescribedFlowType = typeOfPrescribedFlows.at(prescribedFlow);
                (*component)->valueOfPrescribedFlow = valueOfPrescribedFlows.at(prescribedFlow);
            }
        }

        (*component)->currentParameterValue = componentParameterValues.back();
        // make a copy of this value so we can reset it if necessary. Used for e.g. diode state changes.
        (*component)->parameterValueFromInputData = (*component)->currentParameterValue;
        componentParameterValues.pop_back();

        (*component)->startNode->setPressure(initialPressures.at((*component)->startNode->indexInInputData));
        (*component)->endNode->setPressure(initialPressures.at((*component)->endNode->indexInInputData));


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

void Netlist3DDomainReplacement::setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<double> boundaryFlowsOrPressuresAsAppropriate)
{
	m_boundaryFlowsOrPressuresAsAppropriate = boundaryFlowsOrPressuresAsAppropriate;
}

std::vector<double> Netlist3DDomainReplacement::getBoundaryPressures()
{
	std::vector<double> pressures;
	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
	{
		pressures.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->endNode->getPressure());
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