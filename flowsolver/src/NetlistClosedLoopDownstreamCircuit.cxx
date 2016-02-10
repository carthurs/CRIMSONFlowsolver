#include "NetlistClosedLoopDownstreamCircuit.hxx"
#include "boost/pointer_cast.hpp"

// Statics:
int NetlistClosedLoopDownstreamCircuit::s_numberOfDownstreamCircuits = 0;

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

    createVectorsAndMatricesForCircuitLinearSystem();
}

bool NetlistClosedLoopDownstreamCircuit::kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const
{
    bool nodeInterfacesWithDownstreamCircuit = false;
    if (m_pressureNodesWhichConnectToBoundaryCircuits.count(nodeIndex) == 1)
    {
        nodeInterfacesWithDownstreamCircuit = true;
    }
    return nodeInterfacesWithDownstreamCircuit;
}

int NetlistClosedLoopDownstreamCircuit::convertInterfaceNodeIndexFromDownstreamToUpstreamCircuit(const int surfaceIndex, const int sharedNodeDownstreamIndex) const
{
    return m_circuitInterfaceNodeIndexMapDownstreamToUpstream.at(surfaceIndex, sharedNodeDownstreamIndex);
}

int NetlistClosedLoopDownstreamCircuit::convertInterfaceNodeIndexFromUpstreamToDownstreamCircuit(const int sharedNodeUpstreamIndex) const
{
    return m_circuitInterfaceNodeIndexMapUpstreamToDownstream.at(sharedNodeUpstreamIndex);
}

void NetlistClosedLoopDownstreamCircuit::createCircuitDescription()
{
    // This function takes the read-in netlist circuit description and converts it
    // to the internal CircuitData class format.

    // Get the reader class for the netlist data file, and ask it for the circuit description data:
    mp_netlistFileReader = NetlistDownstreamCircuitReader::Instance();
    mp_netlistXmlReader = NetlistDownstreamXmlReader::Instance();
    createBasicCircuitDescription();
    appendClosedLoopSpecificCircuitDescription();
}

void NetlistClosedLoopDownstreamCircuit::appendClosedLoopSpecificCircuitDescription()
{
    // NetlistDownstreamCircuitReader* downcastDownstreamCircuitReader = static_cast<NetlistDownstreamCircuitReader*> (mp_netlistFileReader);
	NetlistDownstreamXmlReader* downcastDownstreamCircuitReader = static_cast<NetlistDownstreamXmlReader*> (mp_netlistXmlReader);
    m_numberOfConnectedBoundaryConditions = downcastDownstreamCircuitReader->getNumberOfBoundaryConditionsConnectedTo(m_downstreamCircuitIndex);
    m_connectedCircuitSurfaceIndices = downcastDownstreamCircuitReader->getConnectedCircuitSurfaceIndices(m_downstreamCircuitIndex);

    for (auto upstreamSurfaceIndex = m_connectedCircuitSurfaceIndices.begin(); upstreamSurfaceIndex != m_connectedCircuitSurfaceIndices.end(); upstreamSurfaceIndex++)
    {
    	m_setOfAttachedBoundaryConditionIndices.insert(*upstreamSurfaceIndex);
    }

    m_localInterfacingNodes = downcastDownstreamCircuitReader->getLocalBoundaryConditionInterfaceNodes(m_downstreamCircuitIndex);
    m_remoteInterfacingNodes = downcastDownstreamCircuitReader->getRemoteBoundaryConditionInterfaceNodes(m_downstreamCircuitIndex);
    
    // Create a useful map which pairs up the two names for each interface node:
    // the node's index as seen by the downstream circuit of the closed loop,
    // and the nodes' index as seen by the upstream bounday condition circuit.
    assert(m_circuitInterfaceNodeIndexMapDownstreamToUpstream.size() == 0);
    for (int interfaceNodeIndex = 0; interfaceNodeIndex < m_localInterfacingNodes.size(); interfaceNodeIndex++)
    {
        m_circuitInterfaceNodeIndexMapDownstreamToUpstream.insert(m_connectedCircuitSurfaceIndices.at(interfaceNodeIndex),
                                                                      m_localInterfacingNodes.at(interfaceNodeIndex),
                                                                      m_remoteInterfacingNodes.at(interfaceNodeIndex)  );
    }
    // Also build the inverse of m_circuitInterfaceNodeIndexMapDownstreamToUpstream:
    assert(m_circuitInterfaceNodeIndexMapUpstreamToDownstream.size() == 0);
    for (int interfaceNodeIndex = 0; interfaceNodeIndex < m_localInterfacingNodes.size(); interfaceNodeIndex++)
    {
        m_circuitInterfaceNodeIndexMapUpstreamToDownstream.insert(std::make_pair(m_remoteInterfacingNodes.at(interfaceNodeIndex),
        																		m_localInterfacingNodes.at(interfaceNodeIndex)  ));
    }                                         
    
}

void NetlistClosedLoopDownstreamCircuit::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    rebuildCircuitMetadata();
    // cycleToSetHistoryPressuresFlowsAndVolumes();
}

void NetlistClosedLoopDownstreamCircuit::giveNodesAndComponentsTheirUpdatedValuesFromSolutionVector(const std::vector<PetscScalar> solutionEntriesForDownstreamCircuit)
{
    // This method receives the subvector from the closed loop m_solutionVector which corresponds
    // to just this downstream circuit.
    //
    // It then passes this data to its circuit components (pressures, flows, volumes).
    
    // Put the raw data we just received into m_solutionVector, so that
    // the function calls below can find it where they expect it to be:
    PetscErrorCode errFlag;
    for (int vectorEntry = 0; vectorEntry < solutionEntriesForDownstreamCircuit.size(); vectorEntry++)
    {
        PetscScalar valueToInsert = solutionEntriesForDownstreamCircuit.at(vectorEntry);
        errFlag = VecSetValue(m_solutionVector,vectorEntry,valueToInsert,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }

    giveNodesTheirPressuresFromSolutionVector();
    giveComponentsTheirFlowsFromSolutionVector();
    giveComponentsTheirVolumesFromSolutionVector();
    // giveComponentsTheirProposedVolumesFromSolutionVector();
}

bool NetlistClosedLoopDownstreamCircuit::boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const
{
	return (m_setOfAttachedBoundaryConditionIndices.count(boundaryConditionIndex) == 1);
}

void NetlistClosedLoopDownstreamCircuit::getMatrixContribution(const double alfi_delt, Mat& matrixFromThisDownstreamCircuit)
{
    generateLinearSystemWithoutFactorisation(alfi_delt);
    matrixFromThisDownstreamCircuit = m_systemMatrix;
}

void NetlistClosedLoopDownstreamCircuit::getRHSContribution(const int timestepNumber, Vec& rhsFromThisDownstreamCircuit)
{
    assembleRHS(timestepNumber);
    rhsFromThisDownstreamCircuit = m_RHS;
}


void NetlistClosedLoopDownstreamCircuit::getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(std::vector<int>& downstreamNodeIndices, std::vector<int>& upstreamNodeIndices, std::vector<int>& upstreamSurfaceIndices) const
{
    downstreamNodeIndices = m_localInterfacingNodes;
    upstreamNodeIndices = m_remoteInterfacingNodes;
    upstreamSurfaceIndices = m_connectedCircuitSurfaceIndices;
}

int NetlistClosedLoopDownstreamCircuit::getCircuitIndex() const
{
    return m_downstreamCircuitIndex;
}

// Disable unwanted methods:
void NetlistClosedLoopDownstreamCircuit::detectWhetherClosedDiodesStopAllFlowAt3DInterface()
{
    throw std::logic_error("Method detectWhetherClosedDiodesStopAllFlowAt3DInterface() should not be called on the NetlistClosedLoopDownstreamCircuit, as it has no 3D interface.");
}