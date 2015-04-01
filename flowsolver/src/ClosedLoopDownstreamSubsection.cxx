#include "ClosedLoopDownstreamSubsection.hxx"
#include <pair>

bool ClosedLoopDownstreamSubsection::boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const
{
	return (m_setOfAttachedBoundaryConditionIndices.count(boundaryConditionIndex) == 1);
}

void ClosedLoopDownstreamSubsection::setPointerToNeighbouringBoundaryCondition(boost::shared_ptr<abstractBoundaryCondition> upstreamBC)
{
    m_upstreamBoundaryConditions.push_back(upstreamBC);
}

void ClosedLoopDownstreamSubsection::initialiseModel()
{
    // Get the input data
    mp_NetlistCircuit->createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    // mp_NetlistCircuit->identifyAtomicSubcircuits();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    mp_NetlistCircuit->closeAllDiodes();

    mp_NetlistCircuit->initialiseCircuit();

    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    //
}

void ClosedLoopDownstreamSubsection::buildAndSolveLinearSystemIfNotYetDone()
{
    // Check whether the linear system still needs to be built and solved; if not, do nothing.

    // Call the upstream boundary conditions to ask for their contributions to the (closed loop)-type
    // linear system:
    for (auto upstreamBoundaryCondition = m_upstreamBoundaryConditions.begin(); upstreamBoundaryCondition != m_upstreamBoundaryConditions.end(); upstreamBoundaryCondition++)
    {
        Mat matrixContribution;
        Vec rhsContribuiton;
        (*upstreamBoundaryCondition)->getMatrixContribuitons(matrixContribution, rhsContribuiton);
        m_matrixContributionsFromUpstreamBoundaryConditions.push_back(matrixContribution);
        m_rhsContributionsFromUpstreamBoundaryConditions.push_back(rhsContribuiton);
    }

    // Tile the matrices to make the big system matrix

    // Add the Kirchoff laws for the connecting nodes, and also set their pressures to be equal

    // Solve the matrix system

    // Set some "done" flag possibly so we know the system doesnt need building and solving - but where to reset this flag again?
}

std::pair<double,double> ClosedLoopDownstreamSubsection::getImplicitCoefficients(const int boundaryConditionIndex) const
{
    // assert the linear system has been solved:

    // Use the known offsets in the matrix of the boundary condition with index boundaryConditionIndex
    // to extract the implicit coefficients
}
