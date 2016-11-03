#ifndef BoundaryConditionFactory_HXX_
#define BoundaryConditionFactory_HXX_

#include <boost/shared_ptr.hpp>
#include "AbstractBoundaryCondition.hxx"
#include "datatypesInCpp.hxx"
#include "ClosedLoopDownstreamSubsection.hxx"

// Forward declarations:
class AbstractBoundaryCondition;

class BoundaryConditionFactory
{
 public:
 	BoundaryConditionFactory(const double hstep, const double delt, const double alfi, const int maxsurf, const int nstep, const int numLoopClosingCircuits, const bool simulationIsPurelyZeroD, const int startingTimestepIndex)
 	: m_hstep(hstep),
 	m_delt(delt),
 	m_alfi(alfi),
 	m_maxsurf(maxsurf),
 	m_nstep(nstep),
 	m_numLoopClosingNetlistCircuits(numLoopClosingCircuits),
 	m_simulationIsPurelyZeroD(simulationIsPurelyZeroD),
 	m_startingTimestepIndex(startingTimestepIndex)
 	{
 		m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt = false;
 	}

	boost::shared_ptr<AbstractBoundaryCondition> createBoundaryCondition(int surfaceIndex_in, boundary_condition_t boundaryType);

	void createNetlistLoopClosingCircuits(std::vector<boost::shared_ptr<ClosedLoopDownstreamSubsection>>& netlistDownstreamLoopClosingSubsections);

 private:
 	// We make boundaryConditionToReturn a member variable for safety: if 
 	// boundaryConditionToReturn throws anything durin construction, we dont
 	// want it being deleted by the shared pointer reference counter
 	// and causing segfaults before the throw error message has been 
 	// presented to the user.
 	boost::shared_ptr<AbstractBoundaryCondition> boundaryConditionToReturn;

	const double m_hstep;
	const double m_delt;
	const double m_alfi;
	const int m_maxsurf;
	const int m_nstep;
	const int m_numLoopClosingNetlistCircuits;
	const bool m_simulationIsPurelyZeroD;
	const int m_startingTimestepIndex;
	std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> mp_netlistDownstreamLoopClosingSubsectionsWeakPointers;

	bool m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt;
};

 #endif