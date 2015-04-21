#include "NetlistBoundaryCircuitWhenDownstreamCircuitsExist.hxx"
#include "ClosedLoopDownstreamSubsection.hxx"

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

bool NetlistBoundaryCircuitWhenDownstreamCircuitsExist::kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const
{
    bool nodeInterfacesWithDownstreamCircuit = false;
    if (m_pressureNodesWhichConnectToDownstreamCircuits.count(nodeIndex) == 1)
    {
        nodeInterfacesWithDownstreamCircuit = true;
    }
    return nodeInterfacesWithDownstreamCircuit;
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
        downstreamSubcircuit->lock()->buildAndSolveLinearSystemIfNotYetDone(timestepNumber, alfi_delt);
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

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getMatrixContribution(const double alfi_delt, Mat& matrixFromThisBoundary)
{
    generateLinearSystemWithoutFactorisation(alfi_delt);
    matrixFromThisBoundary = m_systemMatrix;
}

void NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getRHSContribution(const int timestepNumber, Vec& rhsFromThisBoundary)
{
    assembleRHS(timestepNumber);
    rhsFromThisBoundary = m_RHS;
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getLocationOf3DInterfaceFlowColumnInLinearSystem() const
{
    return columnIndexOf3DInterfaceFlowInLinearSystem.at(0);
}

int NetlistBoundaryCircuitWhenDownstreamCircuitsExist::getCircuitIndex() const
{
    return m_IndexOfThisNetlistLPN;
}