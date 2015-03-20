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
 	boundaryConditionFactory(const double hstep, const double delt, const double alfi, const double lstep, const int maxsurf, const int nstep)
 	: m_hstep(hstep),
 	m_delt(delt),
 	m_alfi(alfi),
 	m_lstep(lstep),
 	m_maxsurf(maxsurf),
 	m_nstep(nstep)
 	{
 	}

	boost::shared_ptr<abstractBoundaryCondition> createBoundaryCondition(int surfaceIndex_in, boundary_condition_t boundaryType);

 private:
	const double m_hstep;
	const double m_delt;
	const double m_alfi;
	const double m_lstep;
	const int m_maxsurf;
	const int m_nstep;
};

 #endif