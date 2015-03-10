#ifndef NETLIST3DDOMAINREPLACEMENT_HXX_
#define NETLIST3DDOMAINREPLACEMENT_HXX_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "datatypesInCpp.hxx"
#include "NetlistCircuit.hxx"

class Netlist3DDomainReplacement
{
public:
	Netlist3DDomainReplacement(const int numberOfNetlistsUsedAsBoundaryConditions, const double oneResistanceToGiveEachResistor, const double elastanceToGiveVolumeTrackingPressureChamber, const double initialDomainPressure, const int hstep, const double alfi_local, const double delt)
	: m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions)
	{
		bool thisIsARestartedSimulation = false; //\todo fix this!
		mp_NetlistZeroDDomainCircuit = boost::shared_ptr<NetlistZeroDDomainCircuit> (new NetlistZeroDDomainCircuit(hstep, m_numberOfNetlistsUsedAsBoundaryConditions, thisIsARestartedSimulation, alfi_local, delt, oneResistanceToGiveEachResistor, elastanceToGiveVolumeTrackingPressureChamber, initialDomainPressure));
	}

	void setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate);
	std::vector<double> getBoundaryPressures();
	std::vector<double> getBoundaryFlows();

	void solveSystem(const int timestepNumber);

	void setPointersToBoundaryPressuresAndFlows(double* const interfacePressuresToBeReadBy3DDomainReplacement, double* const interfaceFlowsToBeReadBy3DDomainReplacement, const int& numberOfPointers);

	void initialiseModel();
	void initialiseAtStartOfTimestep();
	void updateLPN();
	void finalizeLPNAtEndOfTimestep();

	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_zeroDBoundaries_start);

	void setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients);

private:
	boost::shared_ptr<NetlistZeroDDomainCircuit> mp_NetlistZeroDDomainCircuit;
	const int m_numberOfNetlistsUsedAsBoundaryConditions;

	std::vector<std::pair<boundary_data_t,double>> m_boundaryFlowsOrPressuresAsAppropriate;

};

#endif