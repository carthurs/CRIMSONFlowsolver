#ifndef NETLISTCIRCUIT_HXX_
#define NETLISTCIRCUIT_HXX_

#include "gtest/gtest_prod.h"
#include "CircuitData.hxx"
#include <sstream>
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include <set>
#include <boost/weak_ptr.hpp>
#include "ClosedLoopDownstreamSubsection.hxx"
#include "fileReaders.hxx"

class NetlistCircuit
{
	friend class testMultidom;
	FRIEND_TEST(testMultidom,checkNetlistComponentNeighbourPointers);
	FRIEND_TEST(testMultidom, checkClosedDiodeWithRemainingOpenPathDetected);
	FRIEND_TEST(testMultidom, checkClosedDiodeWithoutRemainingOpenPathDetected);
public:
	NetlistCircuit(const int hstep, const int surfaceIndex, const int indexOfThisNetlistLPN, const bool thisIsARestartedSimulation, const double alfi, const double delt)
	: m_surfaceIndex(surfaceIndex),
	m_IndexOfThisNetlistLPN(indexOfThisNetlistLPN),
	m_hstep(hstep),
	m_thisIsARestartedSimulation(thisIsARestartedSimulation),
	m_delt(delt),
	m_alfi(alfi)
	{
		initialisePetscArrayNames();

		safetyCounterLimit = 1000;
		mp_circuitData = boost::shared_ptr<CircuitData> (new CircuitData(m_hstep));
		// mp_circuitDataWithoutDiodes = boost::shared_ptr<CircuitData> (new CircuitData(m_hstep));

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

	virtual void initialiseCircuit();

	bool surfaceIndexMatches(const int surfaceIndexToTest) const;

	bool flowPermittedAcross3DInterface() const;
	bool boundaryConditionTypeHasJustChanged();
	void closeAllDiodes();
	virtual void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
	void switchDiodeStatesIfNecessary();
	void rebuildCircuitMetadata();

	void setPressureAndFlowPointers(double* pressurePointer, double* flowPointer);
	void cycleToSetHistoryPressuresFlowsAndVolumes();

	// void identifyAtomicSubcircuits();
	virtual void initialiseAtStartOfTimestep();
	void finalizeLPNAtEndOfTimestep();
	boost::shared_ptr<CircuitData> getCircuitDescription();
	int getNumberOfDegreesOfFreedom() const;
	std::vector<int> getNodesWithDeferredKirchoffEquations() const;

	virtual void createCircuitDescription();
	virtual ~NetlistCircuit()
	{
		terminatePetscArrays();
	}

	// This can be used to give more than one pressure and one flow pointer to the netlist. Useful if this Netlist
	// has multiple interfaces with other domains (e.g. if this is a Netlist replacement for the 3D domain.)
	void setPointersToBoundaryPressuresAndFlows(double* const interfacePressures, double* const interfaceFlows, const int& numberOfPointers);

	void writePressuresFlowsAndVolumes(int& nextTimestepWrite_start);

	virtual std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt);
	void updateLPN(const int timestepNumber);

	std::pair<boundary_data_t,double> computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber);
	boost::shared_ptr<CircuitComponent> getComponentByInputDataIndex(const int componentIndex);
protected:
	// Overload constructor for subclasses to call:
	NetlistCircuit(const int hstep, const bool thisIsARestartedSimulation, const double alfi, const double delt)
	: m_hstep(hstep),
	m_surfaceIndex(-1),
	m_IndexOfThisNetlistLPN(-1),
	m_thisIsARestartedSimulation(thisIsARestartedSimulation),
	m_delt(delt),
	m_alfi(alfi)
	{
		initialisePetscArrayNames();
	}
	std::string m_PressureHistoryFileName;
	std::string m_FlowHistoryFileName;
	std::string m_VolumeHistoryFileName;
	boost::shared_ptr<CircuitData> mp_circuitData;
	const int m_surfaceIndex;
	const bool m_thisIsARestartedSimulation;
	const double m_delt;
	const double m_alfi;
	const int m_hstep;
	std::vector<double*> pressure_n_ptrs;
	std::vector<double*> flow_n_ptrs;
	int m_NumberOfAtomicSubcircuits;

	NetlistReader* mp_netlistFileReader;

	void createBasicCircuitDescription();
	void createVectorsAndMatricesForCircuitLinearSystem();
	void createListOfNodesWithMultipleIncidentCurrents();
	void getMapOfPressHistoriesToCorrectPressNodes();
	void getMapOfFlowHistoriesToCorrectComponents();
	void getMapOfVolumeHistoriesToCorrectComponents();
	void getMapOfTrackedVolumesToCorrectComponents();
	void generateLinearSystemFromPrescribedCircuit(const double alfi_delt);
	void generateLinearSystemFromPrescribedCircuit(const double alfi_delt);
	void assembleRHS(const int timestepNumber);
	void giveNodesTheirPressuresFromSolutionVector();
	void giveComponentsTheirFlowsFromSolutionVector();
	void giveComponentsTheirVolumesFromSolutionVector();
	void giveComponentsTheirProposedVolumesFromSolutionVector();
	std::vector<double> getVolumesFromSolutionVector();
	bool areThereNegativeVolumes(const int timestepNumber, const double alfi_delt);
	void initialiseCircuit_common();

	Mat m_systemMatrix;
	Mat m_inverseOfSystemMatrix;
	Mat m_identityMatrixForPetscInversionHack;
	Vec m_RHS;
	Vec m_solutionVector;

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
	PetscInt m_numberOfSystemRows;
	PetscInt m_numberOfSystemColumns;
	std::vector<int> listOfNodesWithMultipleIncidentCurrents;
	int m_numberOfMultipleIncidentCurrentNodes;
	std::set<int> listOfHistoryPressures;            // generated from input data, listing pressure node indices and component flow indices where a history is needed (i.e. last time-step values for capacitors/inductors)
	std::set<int> listOfHistoryFlows;
	std::set<int> listOfHistoryVolumes;
	std::set<int> listOfTrackedVolumes;
	int numberOfPrescribedPressuresAndFlows;           // Just the sum of the previous two declared integers
	int m_numberOfHistoryPressures;
	int numberOfHistoryFlows;
	int numberOfHistoryVolumes;
	int m_numberOfTrackedVolumes;
	std::vector<int> columnMap;
	// int columnMapSize;//\todo check this is used
	std::vector<int> columnIndexOf3DInterfaceFlowInLinearSystem;
	std::vector<int> columnIndexOf3DInterfacePressureInLinearSystem;

	std::vector<int> m_nodesWithKirchoffEquationsDeferredToClosedLoop;

	int safetyCounterLimit;

	PetscScalar m_interfaceFlow;
  	PetscScalar m_interfacePressure;

  	const int m_IndexOfThisNetlistLPN;
	
	void buildAndSolveLinearSystem(const int timestepNumber, const double alfi_delt);
	const int getNumberOfHistoryPressures() const;

private:
	void initialisePetscArrayNames();
	void terminatePetscArrays();
	virtual void setupPressureNode(const int indexOfEndNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> component);
	virtual void kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;
	// void createInitialCircuitDescriptionWithoutDiodes();
	// void assignComponentsToAtomicSubcircuits();

	// boost::shared_ptr<CircuitData> mp_circuitDataWithoutDiodes;
	std::vector<boost::shared_ptr<CircuitData>> m_activeSubcircuitCircuitData;
	std::vector<int> m_AtomicSubcircuitsComponentsBelongsTo; // This is indexed by component, as they appear in mp_circuitDataWithoutDiodes

	// std::vector<double> m_PressuresInLPN;                       // Pressure at each LPN node, using the same node indexing as in the netlist
	// std::vector<double> m_HistoryPressuresInLPN;                // As m_PressuresInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.
	// std::vector<double> m_FlowsInLPN;                           // Flow through each component in the LPN, in the order they appear in the netlist
	// std::vector<double> m_HistoryFlowsInLPN;					  // As m_FlowsInLPN, but for any nodes with histories. /Most/ of the entries in this array will never be used.

};

// Forward declaration:
class ClosedLoopDownstreamSubsection;

class NetlistBoundaryCircuitWhenDownstreamCircuitsExist : public NetlistCircuit
{
public:
	NetlistBoundaryCircuitWhenDownstreamCircuitsExist(const int hstep, const int surfaceIndex, const int indexOfThisNetlistLPN, const bool thisIsARestartedSimulation, const double alfi, const double delt, const std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> downstreamSubcircuits)
	:NetlistCircuit(hstep, surfaceIndex, indexOfThisNetlistLPN, thisIsARestartedSimulation, alfi, delt),
	m_netlistDownstreamLoopClosingSubcircuits(downstreamSubcircuits)
	{}
	void initialiseCircuit();
	std::pair<double,double> computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt);
	void getMatrixContribution(Mat& matrixFromThisBoundary);
	void getRHSContribuiton(Vec& rhsFromThisBoundary);
protected:
private:
	std::set<int> m_pressureNodesWhichConnectToDownstreamCircuits;
	int m_numberOfNodesConnectingToAnotherCircuit;
	std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> m_netlistDownstreamLoopClosingSubcircuits;
	void kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;
};

class NetlistZeroDDomainCircuit : public NetlistCircuit
{
public:
	NetlistZeroDDomainCircuit(int hstep, const int numberOfNetlistsUsedAsBoundaryConditions, const bool thisIsARestartedSimulation, const double alfi, const double delt, const double oneResistanceToGiveEachResistor, const double elastanceToGiveVolumeTrackingPressureChamber, const double initialDomainPressure)
	: NetlistCircuit(hstep, thisIsARestartedSimulation, alfi, delt),
	m_oneResistanceToGiveEachResistor(oneResistanceToGiveEachResistor),
	m_elastanceToGiveVolumeTrackingPressureChamber(elastanceToGiveVolumeTrackingPressureChamber),
	m_initialDomainPressure(initialDomainPressure),
	m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions)
	{
		mp_circuitData = boost::shared_ptr<CircuitData> (new Netlist3DDomainReplacementCircuitData(hstep, numberOfNetlistsUsedAsBoundaryConditions));
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
	void initialiseAtStartOfTimestep();
private:
	void setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode);
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	const double m_oneResistanceToGiveEachResistor;
	const double m_elastanceToGiveVolumeTrackingPressureChamber;
	const double m_initialDomainPressure;
};

class NetlistClosedLoopDownstreamCircuit : public NetlistCircuit
{
public:
	NetlistClosedLoopDownstreamCircuit(const int hstep, const bool thisIsARestartedSimulation, const double alfi, const double delt)
	: NetlistCircuit(hstep, thisIsARestartedSimulation, alfi, delt)
	{
		m_downstreamCircuitIndex = s_numberOfDownstreamCircuits;
		s_numberOfDownstreamCircuits++;

		std::stringstream pressureFileNameBuilder;
		pressureFileNameBuilder << "netlistPressures_downstreamCircuit_" << m_downstreamCircuitIndex << ".dat";
		m_PressureHistoryFileName = pressureFileNameBuilder.str();

		std::stringstream flowFileNameBuilder;
		flowFileNameBuilder << "netlistFlows_downstreamCircuit_" << m_downstreamCircuitIndex << ".dat";
		m_FlowHistoryFileName = flowFileNameBuilder.str();

		std::stringstream volumeFileNameBuilder;
		volumeFileNameBuilder << "netlistVolumes_downstreamCircuit_" << m_downstreamCircuitIndex << ".dat";
		m_VolumeHistoryFileName = volumeFileNameBuilder.str();
	}

	void createCircuitDescription();
	void initialiseCircuit();

	void getMatrixContribution(Mat& matrixFromThisBoundary);
	void getRHSContribuiton(Vec& rhsFromThisBoundary);

	int getCircuitIndex() const;

	int convertInterfaceNodeIndexFromDownstreamToUpstreamCircuit(const int sharedNodeDownstreamIndex) const;

	void getSharedNodeDownstreamAndUpstreamAndCircuitUpstreamIndices(std::vector<int>& downstreamNodeIndices, std::vector<int>& upstreamNodeIndices, std::vector<int>& upstreamSurfaceIndices) const;

	~NetlistClosedLoopDownstreamCircuit()
	{
		s_numberOfDownstreamCircuits--;
	}
private:
	static int s_numberOfDownstreamCircuits;
	int m_downstreamCircuitIndex;
	int m_numberOfConnectedBoundaryConditions;
	
	std::vector<int> m_connectedCircuitSurfaceIndices;
	std::vector<int> m_localInterfacingNodes;
	std::vector<int> m_remoteInterfacingNodes;

	std::set<int> m_pressureNodesWhichConnectToBoundaryCircuits;
	std::map<int,int> m_circuitInterfaceNodeIndexMapDownstreamToUpstream;

	void appendClosedLoopSpecificCircuitDescription();
	bool kirchoffEquationAtNodeDeferredToInterfacingCircuit(const int nodeIndex) const;

	// Disabling methods that should never be called:
	void detectWhetherClosedDiodesStopAllFlowAt3DInterface();
};

#endif