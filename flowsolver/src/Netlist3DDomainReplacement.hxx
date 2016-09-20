#ifndef NETLIST3DDOMAINREPLACEMENT_HXX_
#define NETLIST3DDOMAINREPLACEMENT_HXX_

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "datatypesInCpp.hxx"
#include "NetlistZeroDDomainCircuit.hxx"
#include "ControlSystemsManager.hxx"

class Netlist3DDomainReplacement
{
public:
	Netlist3DDomainReplacement(const int numberOfNetlistsUsedAsBoundaryConditions, const double oneResistanceToGiveEachResistor, const double elastanceToGiveVolumeTrackingPressureChamber, const double initialDomainPressure, const int hstep, const double alfi_local, const double delt, const std::map<int,int> mapFromZeroDSurfaceIndexToConnectedComponentIndex, const int startingTimestepIndex, const int numberOfDirichletBCTSurfaces, const int numberOfTimestepsBetweenRestarts)
	: m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions),
	m_startingTimestepIndex(startingTimestepIndex)
	{
		bool thisIsARestartedSimulation = false; //\todo fix this!
		mp_NetlistZeroDDomainCircuit = boost::shared_ptr<NetlistZeroDDomainCircuit> (new NetlistZeroDDomainCircuit(hstep, m_numberOfNetlistsUsedAsBoundaryConditions, thisIsARestartedSimulation, alfi_local, delt, oneResistanceToGiveEachResistor, elastanceToGiveVolumeTrackingPressureChamber, initialDomainPressure, mapFromZeroDSurfaceIndexToConnectedComponentIndex, startingTimestepIndex, numberOfDirichletBCTSurfaces));
		initialiseModel(delt, numberOfTimestepsBetweenRestarts);
	}

	void setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate);
	std::vector<double> getBoundaryPressures();
	std::vector<double> getBoundaryFlows();

	void solveSystem(const int timestepNumber);

	void setPointersToBoundaryPressuresAndFlows(double* const interfacePressuresToBeReadBy3DDomainReplacement, double* const interfaceFlowsToBeReadBy3DDomainReplacement, const int& numberOfPointers);

	void initialiseAtStartOfTimestep();
	void updateLPN(const int timestepNumber);
	void finalizeLPNAtEndOfTimestep();

	boost::shared_ptr<NetlistZeroDDomainCircuit> getCircuit();

	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_zeroDBoundaries_start);
	// void loadPressuresFlowsAndVolumesOnRestart(const int startingTimestepIndex);

	void setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients, std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAtBoundaries);

private:
	boost::shared_ptr<NetlistZeroDDomainCircuit> mp_NetlistZeroDDomainCircuit;
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	const int m_startingTimestepIndex;

	std::vector<std::pair<boundary_data_t,double>> m_boundaryFlowsOrPressuresAsAppropriate;
	boost::shared_ptr<ControlSystemsManager> mp_zeroDDomainControlSystemsManager;

	void initialiseModel(const double delt, const int numberOfTimestepsBetweenRestarts);

};

#endif