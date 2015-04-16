#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include <set>
#include <vector>
#include "datatypesInCpp.hxx"
#include "CircuitData.hxx"
#include "AtomicSubcircuitConnectionManager.hxx"
#include "NetlistCircuit.hxx"
#include <boost/weak_ptr.hpp>

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
	NetlistBoundaryCondition(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double lstep, const int maxsurf, const int nstep, const std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> downstreamSubcircuits)
	: abstractBoundaryCondition(surfaceIndex_in, hstep_in, delt_in, alfi_in, lstep, maxsurf, nstep),
	m_netlistDownstreamLoopClosingSubcircuits(downstreamSubcircuits)
	{
		m_IndexOfThisNetlistLPN = numberOfInitialisedNetlistLPNs;
		if (m_netlistDownstreamLoopClosingSubcircuits.size() > 0)
		{
			mp_NetlistCircuit = boost::shared_ptr<NetlistCircuit> (new NetlistBoundaryCircuitWhenDownstreamCircuitsExist(hstep,surfaceIndex_in, m_IndexOfThisNetlistLPN, thisIsARestartedSimulation, alfi_local, delt, m_netlistDownstreamLoopClosingSubcircuits));
		}
		else
		{
			mp_NetlistCircuit = boost::shared_ptr<NetlistCircuit> (new NetlistCircuit(hstep,surfaceIndex_in, m_IndexOfThisNetlistLPN, thisIsARestartedSimulation, alfi_local, delt));
		}
		// initialiseModel();
		numberOfInitialisedNetlistLPNs++;
	}

	int getIndexAmongstNetlists(){return m_IndexOfThisNetlistLPN;}
 	// void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);
	void updateLPN(const int timestepNumber);
	void initialiseAtStartOfTimestep();
	bool flowPermittedAcross3DInterface();
	bool boundaryConditionTypeHasJustChanged();
	void setDirichletConditionsIfNecessary(int* const binaryMask);
	void finalizeLPNAtEndOfTimestep();
	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);
	boost::shared_ptr<NetlistCircuit> getNetlistCircuit();

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

	std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> m_netlistDownstreamLoopClosingSubcircuits;

	// boost::shared_ptr<AtomicSubcircuitConnectionManager> m_atomicSubcircuitConnectionManager;

	// double P_a;

};

#endif
