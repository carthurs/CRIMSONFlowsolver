#include "NetlistLoopClosingCircuit.hxx"

bool NetlistLoopClosingCircuit::boundaryConditionCircuitConnectsToThisLoopClosingCircuit(const int boundaryConditionIndex) const
{
	return (setOfAttachedBoundaryConditionIndices.count() == 1);
}