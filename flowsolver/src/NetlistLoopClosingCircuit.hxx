#ifndef NETLISTLOOPCLOSINGCIRCUIT_HXX_
#define NETLISTLOOPCLOSINGCIRCUIT_HXX_

#include <set>
#include <boost/shared_ptr.hpp>
#include "NetlistCircuit.hxx"

class NetlistLoopClosingCircuit
{
public:
	NetlistLoopClosingCircuit(index)
	: m_index(index)
	{
	}
	bool boundaryConditionCircuitConnectsToThisLoopClosingCircuit(const int boundaryConditionIndex) const;
private:
	std::set<int> setOfAttachedBoundaryConditionIndices;
	const int m_index;

	boost::shared_ptr<NetlistCircuit> mp_NetlistCircuit;
};

#endif