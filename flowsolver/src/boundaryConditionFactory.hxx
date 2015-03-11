#ifndef BOUNDARYCONDITIONFACTORY_HXX_
#define BOUNDARYCONDITIONFACTORY_HXX_

#include <boost/shared_ptr.hpp>
#include "abstractBoundaryCondition.hxx"

// Forward declarations:
class abstractBoundaryCondition;

class boundaryConditionFactory
{
 public:
 	boundaryConditionFactory(const double hstep, const double delt, const double alfi, const double lstep)
 	: m_hstep(hstep),
 	m_delt(delt),
 	m_alfi(alfi),
 	m_lstep(lstep)
 	{
 	}

	boost::shared_ptr<abstractBoundaryCondition> createBoundaryCondition(int surfaceIndex_in, std::string boundaryType);

 private:
	const double m_hstep;
	const double m_delt;
	const double m_alfi;
	const double m_lstep;
};

 #endif