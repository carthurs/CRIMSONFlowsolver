#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include <set>
#include <vector>
#include "datatypesInCpp.hxx"
#include "CircuitData.hxx"
#include "AtomicSubcircuitConnectionManager.hxx"
#include "NetlistCircuit.hxx"
#include "NetlistBoundaryCircuitWhenDownstreamCircuitsExist.hxx"
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
	NetlistBoundaryCondition(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double startingTimestepIndex, const int maxsurf, const int nstep, const std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> downstreamSubcircuits)
	: abstractBoundaryCondition(surfaceIndex_in, hstep_in, delt_in, alfi_in, startingTimestepIndex, maxsurf, nstep),
	m_netlistDownstreamLoopClosingSubcircuits(downstreamSubcircuits),
	m_startingTimestepIndex(startingTimestepIndex)
	{
		m_IndexOfThisNetlistLPNInInputFile = numberOfInitialisedNetlistLPNs;
		if (m_netlistDownstreamLoopClosingSubcircuits.size() > 0)
		{
			mp_NetlistCircuit = boost::shared_ptr<NetlistCircuit> (new NetlistBoundaryCircuitWhenDownstreamCircuitsExist(hstep,surfaceIndex_in, m_IndexOfThisNetlistLPNInInputFile, m_thisIsARestartedSimulation, alfi_local, delt, m_netlistDownstreamLoopClosingSubcircuits, startingTimestepIndex));
		}
		else
		{
			mp_NetlistCircuit = boost::shared_ptr<NetlistCircuit> (new NetlistCircuit(hstep,surfaceIndex_in, m_IndexOfThisNetlistLPNInInputFile, m_thisIsARestartedSimulation, alfi_local, delt, startingTimestepIndex));
		}
		// initialiseModel();
		numberOfInitialisedNetlistLPNs++;
	}

	int getIndexAmongstNetlists(){return m_IndexOfThisNetlistLPNInInputFile;}
 	// void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);
	virtual void updateLPN(const int timestepNumber);
	void initialiseAtStartOfTimestep();
	bool flowPermittedAcross3DInterface();
	bool boundaryConditionTypeHasJustChanged();
	void setDirichletConditionsIfNecessary(int* const binaryMask);
	void finalizeLPNAtEndOfTimestep();
	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);
	// void loadPressuresFlowsAndVolumesOnRestart(const int startingTimeStepIndex);
	boost::shared_ptr<NetlistCircuit> getNetlistCircuit();

	bool hasPrescribedPressureAcross3DInterface() const;
	bool hasPrescribedFlowAcross3DInterface() const;

	std::vector<double*> getCapacitorNodalHistoryPressurePointers() const;

	~NetlistBoundaryCondition()
	{
		numberOfInitialisedNetlistLPNs--;
	}

	std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement();

	void setPressureAndFlowPointers(double* pressurePointer, double* flowPointer);
	void initialiseModel();
	int m_IndexOfThisNetlistLPNInInputFile;

	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);
	void resetStateUsingKalmanFilteredEstimate(const double flow, const double pressure, const int timestepNumber);
protected:
private:
	boost::shared_ptr<NetlistCircuit> mp_NetlistCircuit;
	static int numberOfInitialisedNetlistLPNs;

	std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> m_netlistDownstreamLoopClosingSubcircuits;
	const int m_startingTimestepIndex;

	int m_NumberOfAtomicSubcircuits; // These are what you get with all valves closed: no subcircuit divides an atomic subcircuit.
	int m_NumberOfSubcircuits;

	// boost::shared_ptr<AtomicSubcircuitConnectionManager> m_atomicSubcircuitConnectionManager;

	// double P_a;

};

#endif
