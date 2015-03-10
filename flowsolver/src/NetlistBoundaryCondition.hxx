#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include <set>
#include <vector>
#include "datatypesInCpp.hxx"
#include "NetlistSubcircuit.hxx"
#include "CircuitData.hxx"
#include "AtomicSubcircuitConnectionManager.hxx"
#include "NetlistCircuit.hxx"

// The NetlistBoundaryCondition is really a manager class for a collection
// of subcircuits, divided by the diodes/valves in the input data.
// It determines which subcircuits should be used on any given time-step,
// and presents to the boundaryConditionManager an interface which is consistent
// with the other boundary condition types (such as the RCR), so boundaryConditionManager
// needn't concern itself with the fine details.
class NetlistBoundaryCondition : public abstractBoundaryCondition
{
	friend class testMultidom;
	FRIEND_TEST(testMultidom,checkNetlistComponentNeighbourPointers);
	FRIEND_TEST(testMultidom,checkClosedDiodeWithRemainingOpenPathDetected);
	FRIEND_TEST(testMultidom,checkClosedDiodeWithoutRemainingOpenPathDetected);
public:
	NetlistBoundaryCondition(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		m_IndexOfThisNetlistLPN = numberOfInitialisedNetlistLPNs;
		mp_NetlistCircuit = boost::shared_ptr<NetlistCircuit> (new NetlistCircuit(hstep,surfaceIndex_in, m_IndexOfThisNetlistLPN, thisIsARestartedSimulation, alfi_local, delt));
		// initialiseModel();
		numberOfInitialisedNetlistLPNs++;
	}

	int getIndexAmongstNetlists(){return m_IndexOfThisNetlistLPN;}

 	// void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);

	void updateLPN();
	void initialiseAtStartOfTimestep();

	bool flowPermittedAcross3DInterface();
	bool boundaryConditionTypeHasJustChanged();

	void setDirichletConditionsIfNecessary(int* const binaryMask);

	void finalizeLPNAtEndOfTimestep();

	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);

	~NetlistBoundaryCondition()
	{
		numberOfInitialisedNetlistLPNs--;
	}

	std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber);

	void setPressureAndFlowPointers(double* pressurePointer, double* flowPointer);
	void initialiseModel();
	int m_IndexOfThisNetlistLPN;

	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);
protected:
private:
	boost::shared_ptr<NetlistCircuit> mp_NetlistCircuit;
	static int numberOfInitialisedNetlistLPNs;

	int m_NumberOfAtomicSubcircuits; // These are what you get with all valves closed: no subcircuit divides an atomic subcircuit.
	int m_NumberOfSubcircuits;

	// boost::shared_ptr<AtomicSubcircuitConnectionManager> m_atomicSubcircuitConnectionManager;

	// double P_a;

};

#endif
