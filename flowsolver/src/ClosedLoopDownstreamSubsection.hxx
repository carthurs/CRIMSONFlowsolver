#ifndef CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_
#define CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_

#include <set>
#include <queue>
#include <boost/shared_ptr.hpp>
#include "abstractBoundaryCondition.hxx"
#include "NetlistClosedLoopDownstreamCircuit.hxx"
#include "NetlistBoundaryCircuitWhenDownstreamCircuitsExist.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"

// Forward declarations:
// class NetlistBoundaryCircuitWhenDownstreamCircuitsExist;

// Instances of this class manage whole closed-loop circuits currently, although they are designed
// to allow for extension to be only part of a single closed loop (in which case, they'll have to 
// be made to play nicely with each other when building the linear system for the whole closed loop)
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
		initialisePetscArrayNames();
		if (m_lstep > 0)
		{
			m_thisIsARestartedSimulation = true;
		}
		else
		{
			m_thisIsARestartedSimulation = false;
		}

		m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep = false;
		m_linearSystemAlreadyUpdatedOnThisTimestep = false;
		m_systemSize = 0;

		mp_NetlistCircuit = boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> (new NetlistClosedLoopDownstreamCircuit(m_hstep, m_thisIsARestartedSimulation, m_alfi, m_delt));
		initialiseModel();
	}

	~ClosedLoopDownstreamSubsection()
	{
		terminatePetscArrays();
	}

	void initialiseAtStartOfTimestep();
	bool boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const;
	void setPointerToNeighbouringBoundaryConditionCircuit(boost::shared_ptr<NetlistCircuit> upstreamBC);
	void buildAndSolveLinearSystemIfNotYetDone(const int timestepNumber, const double alfi_delt);
	void buildAndSolveLinearSystemForUpdateIfNotYetDone(const int timestepNumber, const double alfi_delt);
	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);
	std::pair<double,double> getImplicitCoefficients(const int boundaryConditionIndex) const;
	void markLinearSystemAsNeedingBuildingAgain();
	void markLinearSystemAsNeedingUpdatingAgain();
	std::vector<PetscScalar> getSolutionVectorEntriesCorrespondingToSurface(const int surfaceIndex) const;
	void giveNodesAndComponentsTheirUpdatedValues();
	double getComputedInterfacePressure(const int boundaryConditionIndex) const;
	double getComputedInterfaceFlow(const int boundaryConditionIndex) const;
	int getIndexOfClosedLoop_zeroIndexed() const;
	int getIndexOfClosedLoop_oneIndexed() const;
	boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> getNetlistCircuit() const;
private:
	bool m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep;
	bool m_linearSystemAlreadyUpdatedOnThisTimestep;
	const int m_index;
	const int m_hstep;
	const double m_delt;
	const double m_alfi;
	const int m_lstep;

	int m_systemSize;
	int m_numberOfUpstreamCircuits;
	bool m_thisIsARestartedSimulation;
	Mat m_closedLoopSystemMatrix;
	Mat m_inverseOfClosedLoopMatrix;
	Mat m_identityMatrixForPetscInversionHack;
	Vec m_closedLoopRHS;
	Vec m_solutionVector;
	int m_nextBlankSystemMatrixRow; // zero-indexed
	int m_nextBlankSystemMatrixColumn; // zero-indexed
	int m_nextBlankRhsRow; // zero-indexed
	std::pair<int,int> m_boundsOfDownstreamSolutionDataInSolutionVector;

	std::queue<std::pair<int, Mat>> m_matrixContributionsFromUpstreamBoundaryConditions;
	std::queue<Vec> m_rhsContributionsFromUpstreamBoundaryConditions;
	std::map<int,int> m_indicesOf3DInterfaceComputedFlowsInUpstreamSolutionVectors; // zero-indexed; by upstream index (i.e. order in netlist_surfaces.dat, not surfaceIndex from solver.inp).
	std::map<int,int> m_columnIndicesOf3DInterfacePrescribedFlowsInUpstreamLinearSystems; // zero-indexed; by upstream index (i.e. order in netlist_surfaces.dat, not surfaceIndex from solver.inp).
	std::map<int,int> m_indicesOf3DInterfaceComputedPressuresInUpstreamSolutionVectors; // zero-indexed; by upstream index (i.e. order in netlist_surfaces.dat, not surfaceIndex from solver.inp).
	std::map<int,int> m_columnIndicesOf3DInterfacePrescribedPressuresInUpstreamLinearSystems; // zero-indexed; by upstream index (i.e. order in netlist_surfaces.dat, not surfaceIndex from solver.inp).

	std::vector<boost::shared_ptr<NetlistCircuit>> m_upstreamBoundaryConditionCircuits;
	boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> mp_NetlistCircuit;

	std::map<int,int> m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix; // zero indexed
	std::map<int,int> m_indicesOfFirstColumnOfEachSubcircuitContributionInClosedLoopMatrix; // zero indexed

	std::map<int,std::set<int>> m_mapOfSurfaceIndicesConnectedToEachDownstreamInterfaceNode;
	std::map<int,std::pair<int,int>> m_mapOfSurfaceIndicesToRangeOfEntriesInSolutionVector;

	void terminatePetscArrays();
	void initialisePetscArrayNames();
	void createVectorsAndMatricesForCircuitLinearSystem();
	int getCircuitIndexFromSurfaceIndex(const int upstreamSurfaceIndex) const;
	void generateCircuitInterfaceNodeData();
	void buildAndSolveLinearSystem_internal(const int timestepNumber, const double alfi_delt);

	void initialiseModel();
	void createContiguousIntegerRange(const int startingInteger, const int numberOfIntegers, PetscInt* const arrayToFill);
	void appendKirchoffLawsAtInterfacesBetweenCircuits();
	void enforcePressureEqualityBetweenDuplicatedNodes();
	void writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(const boost::shared_ptr<const CircuitData> circuitData, const int multipleIncidentCurrentNode, const int row, const int numberOfHistoryPressures, const int columnOffset);
	std::vector<PetscScalar> extractContiguousRangeFromPetscVector(Vec vector, const int firstEntry, const int lastEntry) const;

	// std::queue has no clear() method, so we use this instead:
	template <typename Type>
	void clearQueue(std::queue<Type> queueToClear)
	{
		std::queue<Type> emptyQueue;
		std::swap(emptyQueue, queueToClear);
	}
};

#endif