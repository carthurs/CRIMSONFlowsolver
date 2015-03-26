#ifndef CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_
#define CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_

#include <set>
#include <boost/shared_ptr.hpp>
#include "NetlistCircuit.hxx"

class ClosedLoopDownstreamSubsection
{
public:
	ClosedLoopDownstreamSubsection(const int index)
	: m_index(index)
	{
	}
	bool boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const;
private:
	std::set<int> m_setOfAttachedBoundaryConditionIndices;
	const int m_index;

	boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> mp_NetlistCircuit;
};

#endif