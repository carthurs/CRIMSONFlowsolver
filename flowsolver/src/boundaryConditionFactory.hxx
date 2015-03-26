#ifndef BOUNDARYCONDITIONFACTORY_HXX_
#define BOUNDARYCONDITIONFACTORY_HXX_

#include <boost/shared_ptr.hpp>
#include "abstractBoundaryCondition.hxx"
#include "datatypesInCpp.hxx"

// Forward declarations:
class abstractBoundaryCondition;

class boundaryConditionFactory
{
 public:
 	boundaryConditionFactory(const double hstep, const double delt, const double alfi, const double lstep, const int maxsurf, const int nstep, const int numLoopClosingNetlistCircuits)
 	: m_hstep(hstep),
 	m_delt(delt),
 	m_alfi(alfi),
 	m_lstep(lstep),
 	m_maxsurf(maxsurf),
 	m_nstep(nstep),
 	m_numLoopClosingNetlistCircuits(numLoopClosingNetlistCircuits)
 	{
 		m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt = false;
 	}

	boost::shared_ptr<abstractBoundaryCondition> createBoundaryCondition(int surfaceIndex_in, boundary_condition_t boundaryType);

 private:
 	// We make boundaryConditionToReturn a member variable for safety: if 
 	// boundaryConditionToReturn throws anything durin construction, we dont
 	// want it being deleted by the shared pointer reference counter
 	// and causing segfaults before the throw error message has been 
 	// presented to the user.
 	boost::shared_ptr<abstractBoundaryCondition> boundaryConditionToReturn;

	const double m_hstep;
	const double m_delt;
	const double m_alfi;
	const double m_lstep;
	const int m_maxsurf;
	const int m_nstep;
	const int m_numLoopClosingNetlistCircuits;

	bool m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt;
};

 #endif