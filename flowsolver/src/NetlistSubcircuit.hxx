#ifndef NETLISTSUBCIRCUIT_HXX_
#define NETLISTSUBCIRCUIT_HXX_

#include "datatypesInCpp.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include <vector>
#include <set>
#include "CircuitData.hxx"
#include "indexShifters.hxx"

class NetlistSubcircuit
{
public:
	NetlistSubcircuit(const int indexOfThisSubcircuit, const CircuitData& circuitData, double* const flow_n_ptr, double* const pressure_n_ptr, double const alfi_delt_in, const int surfaceIndex_in)
	: indexOfThisSubcircuit(indexOfThisSubcircuit),
	 m_circuitData(circuitData),
	 flow_n_ptr(flow_n_ptr),
	 pressure_n_ptr(pressure_n_ptr),
	 alfi_delt(alfi_delt_in),
	 surfaceIndex(surfaceIndex_in)
	{
		initialiseSubcircuit();
		columnIndexOf3DInterfaceFlowInLinearSystem = 0;
		safetyCounterLimit = 1000;
	}

	~NetlistSubcircuit()
	{
		PetscErrorCode errFlag;
		errFlag = VecDestroy(&RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = VecDestroy(&solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&identityMatrixForPetscInversionHack); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	}

	void updateInternalPressuresVolumesAndFlows();
	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);
	
	const CircuitData m_circuitData;
private:
	void initialiseSubcircuit();
	void createVectorsAndMatricesForCircuitLinearSystem();
	void getListOfNodesWithMultipleIncidentCurrents();
	void getMapOfPressHistoriesToCorrectPressNodes();
	void getMapOfFlowHistoriesToCorrectComponents();
	void getMapOfVolumeHistoriesToCorrectComponents();
	void getMapOfTrackedVolumesToCorrectComponents();
	void generateLinearSystemFromPrescribedCircuit(const double alfi_delt);
	void assembleRHS(const int timestepNumber);
	void computeCircuitLinearSystemSolution();
	void giveNodesTheirPressuresFromSolutionVector();
	void giveComponentsTheirFlowsFromSolutionVector();
	void giveComponentsTheirVolumesFromSolutionVector();
	void giveComponentsTheirProposedVolumesFromSolutionVector();
	std::vector<double> getVolumesFromSolutionVector();
	bool areThereNegativeVolumes();

	const int indexOfThisSubcircuit;
	double* const flow_n_ptr;
	double* const pressure_n_ptr;
	double const alfi_delt;
	const int surfaceIndex;

	Mat systemMatrix;
	Mat inverseOfSystemMatrix;
	Mat identityMatrixForPetscInversionHack;
	Vec RHS;
	Vec solutionVector;

	std::vector<double> pressuresInSubcircuit;
	std::vector<double> historyPressuresInSubcircuit; // As pressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	std::vector<double> flowsInSubcircuit;            // Flow through each component in the LPN, in the order they appear in the netlist
	std::vector<double> historyFlowsInSubcircuit;	  // As flowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	std::vector<double> volumesInSubcircuit;
	std::vector<double> historyVolumesInSubcircuit;
	// circuitData subcircuitInputData;
	std::map<int,int> nodeIndexToPressureHistoryNodeOrderingMap;
	std::map<int,int> componentIndexToFlowHistoryComponentOrderingMap;
	std::map<int,int> componentIndexToVolumeHistoryComponentOrderingMap;
	std::map<int,int> componentIndexToTrackedVolumeComponentOrderingMap;
	PetscInt systemSize;
	std::vector<int> listOfNodesWithMultipleIncidentCurrents;
	int numberOfMultipleIncidentCurrentNodes;
	std::set<int> listOfHistoryPressures;            // generated from input data, listing pressure node indices and component flow indices where a history is needed (i.e. last time-step values for capacitors/inductors)
	std::set<int> listOfHistoryFlows;
	std::set<int> listOfHistoryVolumes;
	std::set<int> listOfTrackedVolumes;
	int numberOfPrescribedPressuresAndFlows;           // Just the sum of the previous two declared integers
	int numberOfHistoryPressures;
	int numberOfHistoryFlows;
	int numberOfHistoryVolumes;
	int numberOfTrackedVolumes;
	std::vector<int> columnMap;
	// int columnMapSize;//\todo check this is used
	int columnIndexOf3DInterfaceFlowInLinearSystem;
	int columnIndexOf3DInterfacePressureInLinearSystem;

	int safetyCounterLimit;

};

#endif