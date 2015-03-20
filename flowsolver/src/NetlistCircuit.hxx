#ifndef NETLISTCIRCUIT_HXX_
#define NETLISTCIRCUIT_HXX_

#include "gtest/gtest_prod.h"
#include "CircuitData.hxx"
#include "NetlistSubcircuit.hxx"
#include <sstream>

class NetlistCircuit
{
	friend class testMultidom;
	FRIEND_TEST(testMultidom,checkNetlistComponentNeighbourPointers);
	FRIEND_TEST(testMultidom, checkClosedDiodeWithRemainingOpenPathDetected);
	FRIEND_TEST(testMultidom, checkClosedDiodeWithoutRemainingOpenPathDetected);
public:
	NetlistCircuit(const int hstep, const int surfaceIndex, const int indexOfThisNetlistLPN, const bool thisIsARestartedSimulation, const double alfi_local, const double delt)
	: m_surfaceIndex(surfaceIndex),
	m_IndexOfThisNetlistLPN(indexOfThisNetlistLPN),
	m_hstep(hstep),
	m_thisIsARestartedSimulation(thisIsARestartedSimulation),
	m_delt(delt),
	m_alfi_local(alfi_local)
	{
		mp_CircuitDescription = boost::shared_ptr<CircuitData> (new CircuitData(m_hstep));
		mp_CircuitDescriptionWithoutDiodes = boost::shared_ptr<CircuitData> (new CircuitData(m_hstep));

		std::stringstream pressureFileNameBuilder;
		pressureFileNameBuilder << "netlistPressures_surface_" << m_surfaceIndex << ".dat";
		m_PressureHistoryFileName = pressureFileNameBuilder.str();

		std::stringstream flowFileNameBuilder;
		flowFileNameBuilder << "netlistFlows_surface_" << m_surfaceIndex << ".dat";
		m_FlowHistoryFileName = flowFileNameBuilder.str();

		std::stringstream volumeFileNameBuilder;
		volumeFileNameBuilder << "netlistVolumes_surface_" << m_surfaceIndex << ".dat";
		m_VolumeHistoryFileName = volumeFileNameBuilder.str();
	}

	bool flowPermittedAcross3DInterface() const;
	bool boundaryConditionTypeHasJustChanged();
	void closeAllDiodes();
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
	void switchDiodeStatesIfNecessary();
	void rebuildCircuitMetadata();

	void setPressureAndFlowPointers(double* pressurePointer, double* flowPointer);
	void cycleToSetHistoryPressuresFlowsAndVolumes();

	void identifyAtomicSubcircuits();
	void initialiseAtStartOfTimestep();
	void finalizeLPNAtEndOfTimestep();
	boost::shared_ptr<CircuitData> getCircuitDescription();

	virtual void createCircuitDescription();
	virtual void buildSubcircuit();
	virtual ~NetlistCircuit() {}

	// This can be used to give more than one pressure and one flow pointer to the netlist. Useful if this Netlist
	// has multiple interfaces with other domains (e.g. if this is a Netlist replacement for the 3D domain.)
	void setPointersToBoundaryPressuresAndFlows(double* const interfacePressures, double* const interfaceFlows, const int& numberOfPointers);

	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);

	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt);
	void updateLPN(const int timestepNumber);

	std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber);
	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);
protected:
	// Overload constructor for subclasses to call:
	NetlistCircuit(const int hstep, const bool thisIsARestartedSimulation, const double alfi_local, const double delt)
	: m_hstep(hstep),
	m_surfaceIndex(-1),
	m_IndexOfThisNetlistLPN(-1),
	m_thisIsARestartedSimulation(thisIsARestartedSimulation),
	m_delt(delt),
	m_alfi_local(alfi_local)
	{
	}
	std::string m_PressureHistoryFileName;
	std::string m_FlowHistoryFileName;
	std::string m_VolumeHistoryFileName;
	boost::shared_ptr<CircuitData> mp_CircuitDescription;
	boost::shared_ptr<NetlistSubcircuit> mp_subcircuit;
	const int m_surfaceIndex;
	const bool m_thisIsARestartedSimulation;
	const double m_delt;
	const double m_alfi_local;
	const int m_hstep;
	std::vector<double*> pressure_n_ptrs;
	std::vector<double*> flow_n_ptrs;
	int m_NumberOfAtomicSubcircuits;
private:
	virtual void setupPressureNode(const int indexOfEndNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> component);
	void createInitialCircuitDescriptionWithoutDiodes();
	void assignComponentsToAtomicSubcircuits();

	boost::shared_ptr<CircuitData> mp_CircuitDescriptionWithoutDiodes;
	std::vector<boost::shared_ptr<CircuitData>> m_activeSubcircuitCircuitData;
	std::vector<int> m_AtomicSubcircuitsComponentsBelongsTo; // This is indexed by component, as they appear in mp_CircuitDescriptionWithoutDiodes

	const int m_IndexOfThisNetlistLPN;
	// std::vector<double> m_PressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	// std::vector<double> m_HistoryPressuresInLPN;                // As m_PressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	// std::vector<double> m_FlowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	// std::vector<double> m_HistoryFlowsInLPN;					  // As m_FlowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.

};

class NetlistZeroDDomainCircuit : public NetlistCircuit
{
public:
	NetlistZeroDDomainCircuit(int hstep, const int numberOfNetlistsUsedAsBoundaryConditions, const bool thisIsARestartedSimulation, const double alfi_local, const double delt, const double oneResistanceToGiveEachResistor, const double elastanceToGiveVolumeTrackingPressureChamber, const double initialDomainPressure)
	: NetlistCircuit(hstep, thisIsARestartedSimulation, delt, alfi_local),
	m_oneResistanceToGiveEachResistor(oneResistanceToGiveEachResistor),
	m_elastanceToGiveVolumeTrackingPressureChamber(elastanceToGiveVolumeTrackingPressureChamber),
	m_initialDomainPressure(initialDomainPressure),
	m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions)
	{
		mp_CircuitDescription = boost::shared_ptr<CircuitData> (new Netlist3DDomainReplacementCircuitData(hstep, numberOfNetlistsUsedAsBoundaryConditions));
		m_PressureHistoryFileName.clear(); // Defensive
		m_PressureHistoryFileName.append("netlistPressures_zeroDDomainReplacement.dat");
		m_FlowHistoryFileName.clear(); // Defensive
		m_FlowHistoryFileName.append("netlistFlows_zeroDDomainReplacement.dat");
		m_VolumeHistoryFileName.clear(); // Defensive
		m_VolumeHistoryFileName.append("netlistVolumes_zeroDDomainReplacement.dat");
	}

	void setBoundaryPrescriptionsAndBoundaryConditionTypes(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate);

	std::vector<double> getBoundaryPressures();
	std::vector<double> getBoundaryFlows();
	void solveSystem(const int timestepNumber);
	void setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients, std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAtBoundaries);
	void createCircuitDescription();
private:
	void setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode);
	void buildSubcircuit();
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	const double m_oneResistanceToGiveEachResistor;
	const double m_elastanceToGiveVolumeTrackingPressureChamber;
	const double m_initialDomainPressure;
};

#endif