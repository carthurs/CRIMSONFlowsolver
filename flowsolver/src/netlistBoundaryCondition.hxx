#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include <set>
#include <vector>
#include "datatypesInCpp.hxx"

struct circuitData
{
	std::vector<circuit_component_t> componentTypes; // the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	std::vector<int> componentStartNodes;
	std::vector<int> componentEndNodes;
	std::vector<double> componentParameterValues;
	std::vector<circuit_nodal_pressure_prescription_t> typeOfPrescribedPressures;
	std::vector<circuit_component_flow_prescription_t> typeOfPrescribedFlows;
};

class netlistBoundaryCondition : public abstractBoundaryCondition
{
public:
	netlistBoundaryCondition(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		indexOfThisNetlistLPN = numberOfInitialisedNetlistLPNs;
		numberOfInitialisedNetlistLPNs++;
		initialiseModel();
	}

 	void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timen_1, double alfi_delt);
	void initialiseModel();

	void updateLPN();

	~netlistBoundaryCondition()
	{
		numberOfInitialisedNetlistLPNs--;
		PetscErrorCode errFlag;
		errFlag = VecDestroy(&RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = VecDestroy(&solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
		errFlag = MatDestroy(&identityMatrixForPetscInversionHack); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	}

private:
	static int numberOfInitialisedNetlistLPNs;
	int indexOfThisNetlistLPN;

	void getMapOfPressHistoriesToCorrectPressNodes();
	void getMapOfFlowHistoriesToCorrectComponents();
	void getListOfNodesWithMultipleIncidentCurrents();
	void generateLinearSystemFromPrescribedCircuit(double alfi_delt);
	void assembleRHS_netlistLPN(int timestepNumber);

	Mat systemMatrix;
	Mat inverseOfSystemMatrix;
	Mat identityMatrixForPetscInversionHack;
	Vec RHS;
	Vec solutionVector;
	std::vector<double> pressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	std::vector<double> historyPressuresInLPN;                // As pressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used
	std::vector<double> flowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	int numberOfComponents;
	int numberOfPressureNodes;
	// integer, allocatable :: localToGlobalSurfaceIndexMap(:)
	circuitData circuitInputData;
	std::map<int,int> nodeIndexToPressureHistoryNodeOrderingMap;
	std::map<int,int> componentIndexToFlowHistoryComponentOrderingMap;
	PetscInt systemSize;
	std::vector<int> listOfNodesWithMultipleIncidentCurrents;
	int numberOfMultipleIncidentCurrentNodes;
	std::vector<int> listOfPrescribedPressures;         // input data, listing node indices of prescribed pressures (e.g. LV pressure on a capacitor, venous pressure in open loop, etc.).
	std::vector<int> listOfPrescribedFlows;             // input data, listing node indices of prescribed flows, listing the component indices with prescribed flows (e.g. 3D outlet flow)
	std::set<int> listOfHistoryPressures;            // generated from input data, listing pressure node indices and component flow indices where a history is needed (i.e. last time-step values for capacitors/inductors)
	std::set<int> listOfHistoryFlows;
	int numberOfPrescribedPressures;
	std::vector<double> valueOfPrescribedPressures;
	std::vector<double> valueOfPrescribedFlows;
	int numberOfPrescribedFlows;
	int numberOfPrescribedPressuresAndFlows;           // Just the sum of the previous two declared integers
	int numberOfHistoryPressures;
	int numberOfHistoryFlows;
	std::vector<int> columnMap;
	double P_a;
	int columnIndexOf3DInterfaceFlowInLinearSystem;
	// std::vector<double> initialPressures;
	int columnMapSize;//\todo check this is used

};

#endif
