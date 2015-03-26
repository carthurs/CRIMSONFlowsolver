#include "ClosedLoopDownstreamSubsection.hxx"

bool ClosedLoopDownstreamSubsection::boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const
{
	return (m_setOfAttachedBoundaryConditionIndices.count(boundaryConditionIndex) == 1);
}

// void ClosedLoopDownstreamSubsection::initialiseModel()
// {
    // // Get the input data
    // mp_NetlistCircuit->createCircuitDescription();

    // // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    // // mp_NetlistCircuit->identifyAtomicSubcircuits();

    // // Initialise all diodes to their closed state, for stability
    // //\todo change this if you're restarting and the diodes need to be open at restart!
    // assert(false);
    // mp_NetlistCircuit->closeAllDiodes();
    // mp_NetlistCircuit->detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    // mp_NetlistCircuit->initialiseCircuit();

    // // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // // what connections should be made when a diode/valve opens.
    // // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    // //
// }