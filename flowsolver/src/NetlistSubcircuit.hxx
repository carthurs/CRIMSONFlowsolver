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
	NetlistSubcircuit(boost::shared_ptr<CircuitData> circuitData, const std::vector<double*> flow_n_ptrs, const std::vector<double*> pressure_n_ptrs, const int surfaceIndex_in, const double delt)
	: 
	 mp_circuitData(circuitData),
	 flow_n_ptrs(flow_n_ptrs),
	 pressure_n_ptrs(pressure_n_ptrs),
	 surfaceIndex(surfaceIndex_in),
	 m_delt(delt)
	{
		initialiseSubcircuit();
		// columnIndexOf3DInterfaceFlowInLinearSystem = 0;
		safetyCounterLimit = 1000;
	}

	~NetlistSubcircuit()
	{
		PetscErrorCode errFlag;
		errFlag = VecDestroy(&RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = VecDestroy(&solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&m_systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&m_inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&m_identityMatrixForPetscInversionHack); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	}

	void updateInternalPressuresVolumesAndFlows_internal(const int timestepNumber, const double alfi_delt);
	std::pair<double,double> computeImplicitCoefficients_internal(const int timestepNumber, const double timen_1, const double alfi_delt);
	std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement_internal(const int timestepNumber);
	
	boost::shared_ptr<CircuitData> mp_circuitData;
	void buildAndSolveLinearSystem_internal(const int timestepNumber, const double alfi_delt);
protected:

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
	void computeCircuitLinearSystemSolution(const int timestepNumber, const double alfi_delt);
	void giveNodesTheirPressuresFromSolutionVector();
	void giveComponentsTheirFlowsFromSolutionVector();
	void giveComponentsTheirVolumesFromSolutionVector();
	void giveComponentsTheirProposedVolumesFromSolutionVector();
	std::vector<double> getVolumesFromSolutionVector();
	bool areThereNegativeVolumes(const int timestepNumber, const double alfi_delt);

	const std::vector<double*> flow_n_ptrs;
	const std::vector<double*> pressure_n_ptrs;
	const int surfaceIndex;
	const double m_delt;

	Mat m_systemMatrix;
	Mat m_inverseOfSystemMatrix;
	Mat m_identityMatrixForPetscInversionHack;
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
	int m_numberOfTrackedVolumes;
	std::vector<int> columnMap;
	// int columnMapSize;//\todo check this is used
	std::vector<int> columnIndexOf3DInterfaceFlowInLinearSystem;
	std::vector<int> columnIndexOf3DInterfacePressureInLinearSystem;

	int safetyCounterLimit;

	PetscScalar m_interfaceFlow;
  	PetscScalar m_interfacePressure;

};

#endif