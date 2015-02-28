#ifndef PUREZERODDRIVER_HXX_
#define PUREZERODDRIVER_HXX_

#include "NetlistBoundaryCondition.hxx"
#include <boost/unique_ptr.hpp>

class PureZeroDDriver
{
public:
	void init();
	void iter_init();
	void iter_step();
	void iter_finalize();
	void finalize();
private:
	// this is not really a boundary condition here; we just use the machinery of the Netlist to make
	// a replacement for the 3D domain.
	boost::unique_ptr<NetlistBoundaryCondition> zeroDDomainLPN;
};

#endif