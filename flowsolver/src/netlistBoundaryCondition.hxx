#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"

class circuitData
{
public:
	std::vector<char> componentTypes; // the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.
	std::vector<int> componentStartNodes;
	std::vector<int> componentEndNodes;
	std::vector<double> componentParameterValues;
};

class netlistBoundaryCondition : public abstractBoundaryCondition
{
public:
	netlistBoundaryCondition(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		initialiseModel();
		indexOfThisNetlistLPN = numberOfInitialisedNetlistLPNs;
		numberOfInitialisedNetlistLPNs++;
	}
	void computeImplicitCoeff_solve(int timestepNumber)
	{

	}
 	void computeImplicitCoeff_update(int timestepNumber)
 	{

 	}
 	void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timen_1, double alfi_delt)
 	{
 		std::pair<double,double> dummyValue;
 		dummyValue.first=-3.14;
 		dummyValue.second=-2.718281828;
 		return dummyValue;
 	}
	void initialiseModel();

	~netlistBoundaryCondition()
	{
		numberOfInitialisedNetlistLPNs--;
	}

private:
	static int numberOfInitialisedNetlistLPNs;
	int indexOfThisNetlistLPN;

	Mat systemMatrix;
	Vec RHS;
	std::vector<double> pressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	std::vector<double> historyPressuresInLPN;                // As pressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used
	std::vector<double> flowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	Vec solutionVector;
	int numberOfComponents;
	int numberOfPressureNodes;
	// integer, allocatable :: localToGlobalSurfaceIndexMap(:)
	circuitData circuitInputData;
	// integer, allocatable :: nodeIndexToPressureHistoryNodeOrderingMap(:,:)
	// integer, allocatable :: componentIndexToFlowHistoryComponentOrderingMap(:,:)
	PetscInt systemSize;
	std::vector<int> listOfNodesWithMultipleIncidentCurrents;
	int numberOfMultipleIncidentCurrentNodes;
	std::vector<int> listOfPrescribedPressures;         // input data, listing node indices of prescribed pressures (e.g. LV pressure on a capacitor, venous pressure in open loop, etc.).
	std::vector<int> listOfPrescribedFlows;             // input data, listing node indices of prescribed flows, listing the component indices with prescribed flows (e.g. 3D outlet flow)
	std::vector<int> listOfHistoryPressures;            // generated from input data, listing pressure node indices and component flow indices where a history is needed (i.e. last time-step values for capacitors/inductors)
	std::vector<int> listOfHistoryFlows;
	int numberOfPrescribedPressures;
	std::vector<double> valueOfPrescribedPressures;
	std::vector<double> valueOfPrescribedFlows;
	std::vector<char> typeOfPrescribedPressures;
	int numberOfPrescribedFlows;
	std::vector<char> typeOfPrescribedFlows;
	int numberOfPrescribedPressuresAndFlows;           // Just the sum of the previous two declared integers
	int numberOfHistoryPressures;
	int numberOfHistoryFlows;
	std::vector<int> columnMap;
	Mat inverseOfSystemMatrix;
	double P_a;
	int columnIndexOf3DInterfaceFlowInLinearSystem;
	int columnMapSize;//\todo check this is used

};

#endif