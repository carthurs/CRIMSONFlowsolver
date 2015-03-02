#ifndef PUREZERODDRIVER_HXX_
#define PUREZERODDRIVER_HXX_

#include "NetlistBoundaryCondition.hxx"
#include "Netlist3DDomainReplacement.hxx"
#include "boundaryConditionManager.hxx"
#include <boost/interprocess/smart_ptr/unique_ptr.hpp>

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
	boost::interprocess::unique_ptr<Netlist3DDomainReplacement> m_zeroDDomainLPN;
	boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
	std::vector<double> m_pressuresOrFlowsAtBoundaries;
	std::vector<double> m_FlowsOrPressuresToGiveToBoundaries;
	double* mp_interfaceFlows;
	double* mp_interfacePressures;
	int m_timestepNumber;

	void placePressuresAndFlowsInStorageArrays(std::vector<double> boundaryPressures, std::vector<double> boundaryFlows);
};

#endif