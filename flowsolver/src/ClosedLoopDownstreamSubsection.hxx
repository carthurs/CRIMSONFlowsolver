#ifndef CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_
#define CLOSEDLOOPDOWNSTREAMSUBSECTION_HXX_

#include <set>
#include <queue>
#include <boost/shared_ptr.hpp>
#include "NetlistCircuit.hxx"
#include "abstractBoundaryCondition.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"

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
		m_systemSize = 0;

		mp_NetlistCircuit = boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> (new NetlistClosedLoopDownstreamCircuit(m_hstep, m_thisIsARestartedSimulation, m_alfi, m_delt));
	}

	~ClosedLoopDownstreamSubsection()
	{
		terminatePetscArrays();
	}

	bool boundaryConditionCircuitConnectsToThisDownstreamSubsection(const int boundaryConditionIndex) const;
	void setPointerToNeighbouringBoundaryConditionCircuit(boost::shared_ptr<NetlistCircuit> upstreamBC);
	void buildAndSolveLinearSystemIfNotYetDone();
	std::pair<double,double> getImplicitCoefficients(const int boundaryConditionIndex) const;
	void markLinearSystemAsNeedingBuildingAgain();
private:
	bool m_linearSystemAlreadyBuiltAndSolvedOnThisTimestep;
	std::set<int> m_setOfAttachedBoundaryConditionIndices;
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
	Vec m_closedLoopRHS;
	Vec m_solutionVector;
	int m_nextBlankSystemMatrixRow; // zero-indexed
	int m_nextBlankSystemMatrixColumn; // zero-indexed
	int m_nextBlankRhsRow; // zero-indexed

	std::queue<Mat> m_matrixContributionsFromUpstreamBoundaryConditions;
	std::queue<Vec> m_rhsContributionsFromUpstreamBoundaryConditions;

	std::vector<boost::shared_ptr<NetlistCircuit>> m_upstreamBoundaryConditionCircuits;
	boost::shared_ptr<NetlistClosedLoopDownstreamCircuit> mp_NetlistCircuit;

	std::vector<int> m_indicesOfFirstRowOfEachSubcircuitContributionInClosedLoopMatrix; // zero indexed

	void terminatePetscArrays();
	void initialisePetscArrayNames();
	void createVectorsAndMatricesForCircuitLinearSystem();

	void initialiseModel();
	void createContiguousIntegerRange(const int startingInteger, const int numberOfIntegers, PetscInt* const arrayToFill);
	void appendKirchoffLawsAtInterfacesBetweenCircuits();
	void enforcePressureEqualityBetweenDuplicatedNodes();
	void writePartOfKirchoffEquationIntoClosedLoopSysteMatrix(const boost::shared_ptr<const CircuitData> circuitData, const int multipleIncidentCurrentNode, const int row, const int numberOfHistoryPressures, const int columnOffset);
};

#endif