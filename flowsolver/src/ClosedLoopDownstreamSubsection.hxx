#ifndef CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_
#define CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_

#include <set>
#include <boost/shared_ptr.hpp>
#include "NetlistCircuit.hxx"
#include "abstractBoundaryCondition.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"

class ClosedLoopDownstreamSubsection
{
public:
	ClosedLoopDownstreamSubsection(const int index, const int hstep, const double delt, const double alfi, const int lstep)
	: m_index(index),
	m_hstep(hstep),
	m_delt(delt),
	m_alfi(alfi),
	m_lstep(lstep)
	{
		if (m_lstep > 0)
		{
			m_thisIsARestartedSimulation = true;
		}
		else
		{
			m_thisIsARestartedSimulation = false;
		}

		mp_NetlistCircuit = boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> (new NetlistClosedLoopDownstreamCircuit(m_hstep, m_thisIsARestartedSimulation, m_alfi, m_delt));
	}
	bool boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const;
	void setPointerToNeighbouringBoundaryCondition(boost::shared_ptr<abstractBoundaryCondition> upstreamBC);
	void buildAndSolveLinearSystemIfNotYetDone();
	std::pair<double,double> getImplicitCoefficients(const int boundaryConditionIndex) const;
private:
	std::set<int> m_setOfAttachedBoundaryConditionIndices;
	const int m_index;
	bool m_thisIsARestartedSimulation;

	std::vector<Mat> m_matrixContributionsFromUpstreamBoundaryConditions;
	std::vector<Vec> m_rhsContributionsFromUpstreamBoundaryConditions;

	std::vector<boost::shared_ptr<abstractBoundaryCondition>> m_upstreamBoundaryConditions;
	boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> mp_NetlistCircuit;

	void initialiseModel();
};

#endif