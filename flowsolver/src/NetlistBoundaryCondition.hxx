#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"
#include <set>
#include <vector>
#include "datatypesInCpp.hxx"
#include "NetlistSubcircuit.hxx"
#include "CircuitData.hxx"
#include "AtomicSubcircuitConnectionManager.hxx"

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
	: abstractBoundaryCondition(surfaceIndex_in),
	  m_CircuitDescription(hstep),
	  m_CircuitDescriptionWithoutDiodes(hstep)
	{
		m_IndexOfThisNetlistLPN = numberOfInitialisedNetlistLPNs;
		initialiseModel();
		numberOfInitialisedNetlistLPNs++;
	}

 	void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt);

	void updateLPN();
	void initialiseAtStartOfTimestep();

	CircuitData& getCircuitDescription();

	void setDirichletConditionsIfNecessary(int* const binaryMask);

	~NetlistBoundaryCondition()
	{
		numberOfInitialisedNetlistLPNs--;
	}

private:
	static int numberOfInitialisedNetlistLPNs;
	int m_IndexOfThisNetlistLPN;

	int m_NumberOfAtomicSubcircuits; // These are what you get with all valves closed: no subcircuit divides an atomic subcircuit.
	int m_NumberOfSubcircuits;

	// boost::shared_ptr<AtomicSubcircuitConnectionManager> m_atomicSubcircuitConnectionManager;

	void initialiseModel();
	void createCircuitDescription();
	void identifySubciruits();
	void selectAndBuildActiveSubcircuits();
	void createInitialCircuitDescriptionWithoutDiodes();
	void assignComponentsToAtomicSubcircuits();
	void createSubcircuitDescriptions();
	void cycleToSetHistoryPressuresAndFlows();
	void netlistBoundaryCondition();
	void createAtomicSubcircuitDescriptions();
	void identifyAtomicSubcircuits();

	CircuitData m_CircuitDescription;
	CircuitData m_CircuitDescriptionWithoutDiodes;
	std::vector<boost::shared_ptr<CircuitData>> m_CircuitDataForAtomicSubcircuits;
	std::vector<boost::shared_ptr<CircuitData>> m_activeSubcircuitCircuitData;
	std::vector<boost::shared_ptr<NetlistSubcircuit>> m_activeSubcircuits;

	std::vector<int> m_AtomicSubcircuitsComponentsBelongsTo; // This is indexed by component, as they appear in m_CircuitDescriptionWithoutDiodes

	std::vector<double> m_PressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	std::vector<double> m_HistoryPressuresInLPN;                // As m_PressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	std::vector<double> m_FlowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	std::vector<double> m_HistoryFlowsInLPN;					  // As m_FlowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.

	// double P_a;

};

#endif
