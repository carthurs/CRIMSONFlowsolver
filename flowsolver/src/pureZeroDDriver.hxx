#ifndef PUREZERODDRIVER_HXX_
#define PUREZERODDRIVER_HXX_

#include "NetlistBoundaryCondition.hxx"
#include "Netlist3DDomainReplacement.hxx"
#include "boundaryConditionManager.hxx"
#include <boost/shared_ptr.hpp>

class PureZeroDDriver
{
public:
	PureZeroDDriver(const int numberOfDirichletSurfaces, const double centralCapacitorCompliance) :
	m_numberOfDirichletSurfaces(numberOfDirichletSurfaces),
	m_complianceToGiveCentralCapacitor(centralCapacitorCompliance)
	{
		boundaryConditionManager_instance = boundaryConditionManager::Instance();
		checkIfThisIsARestartedSimulation();

		m_deltHasBeenSet = false;
		m_alfiHasBeenSet = false;
		m_hstepHasBeenSet = false;
		m_numberOfTimestepsBetweenRestartsHasBeenSet = false;
		m_mapFromNetlistIndexToConnectedComponentsHasBeenSet = false;
	}
	void init();
	void iter_init();
	void iter_step();
	void iter_finalize();
	void finalize();

	void setDelt(const double delt);
	void setAlfi(const double alfi);
	void setHstep(const int hstep);
	void setNtout(const int numberOfTimestepsBetweenRestarts);
	void setupConnectedComponents(const int num3DConnectedComponents, const int* const surfacesOfEachConnectedComponent, const int* const indicesOfNetlistSurfaces, const int* const indicesOfDirichletSurfaces);

private:
	// this is not really a boundary condition here; we just use the machinery of the Netlist to make
	// a replacement for the 3D domain.
	boost::shared_ptr<Netlist3DDomainReplacement> m_zeroDDomainLPN;
	boundaryConditionManager* boundaryConditionManager_instance;
	std::vector<std::pair<boundary_data_t,double>> m_pressuresOrFlowsAtBoundaries;
	double* mp_interfaceFlowsToBeReadByBoundaryConditions;
	double* mp_interfacePressuresToBeReadByBoundaryConditions;

	double* mp_interfaceFlowsToBeReadBy3DDomainReplacement;
	double* mp_interfacePressuresToBeReadBy3DDomainReplacement;

	int m_timestepNumber;
	int m_nextTimestepWrite_zeroDBoundaries_start;
	bool m_thisIsARestartedSimulation;

	double m_alfi;
	double m_delt;
	int m_hstep;
	int m_numberOfTimestepsBetweenRestarts;
	const int m_numberOfDirichletSurfaces;
	const double m_complianceToGiveCentralCapacitor;

	int m_numberOfNetlistsUsedAsBoundaryConditions;

	// the zero-D surfaces are just indexed by consecutive integers, staring from 1.
	// This is used for internal purposes only - the inlets/outlets still have the surface index
	// of their corresponding surface in the 3D model.
	std::map<int,int> m_mapFromZeroDSurfaceIndexToConnectedComponentIndex;

	bool m_deltHasBeenSet;
	bool m_alfiHasBeenSet;
	bool m_hstepHasBeenSet;
	bool m_numberOfTimestepsBetweenRestartsHasBeenSet;
	bool m_mapFromNetlistIndexToConnectedComponentsHasBeenSet;

	void placePressuresAndFlowsInStorageArrays_toGiveToBoundaryConditions(std::vector<double> boundaryPressures, std::vector<double> boundaryFlows);
	void placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement();
	void checkIfThisIsARestartedSimulation();
	void writeNumstartDotDat();
};

#endif