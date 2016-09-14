#ifndef NETLISTZERODDOMAINCIRCUIT_HXX_
#define NETLISTZERODDOMAINCIRCUIT_HXX_

#include "NetlistCircuit.hxx"

class NetlistZeroDDomainCircuit : public NetlistCircuit
{
public:
	NetlistZeroDDomainCircuit(int hstep, const int numberOfNetlistsUsedAsBoundaryConditions, const bool thisIsARestartedSimulation, const double alfi, const double delt, const double oneResistanceToGiveEachResistor, const double complianceToGiveCentralCapacitor, const double initialDomainPressure, const std::map<int,int> mapFromNetlistIndexAmongstNetlistsToConnectedComponentIndex, const int startingTimestepIndex)
	: NetlistCircuit(hstep, -1, thisIsARestartedSimulation, alfi, delt, startingTimestepIndex),
	m_oneResistanceToGiveEachResistor(oneResistanceToGiveEachResistor),
	m_elastanceToGiveCentralCapacitor(complianceToGiveCentralCapacitor),
	m_initialDomainPressure(initialDomainPressure),
	m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions),
	m_mapFromNetlistIndexAmongstNetlistsToConnectedComponentIndex(mapFromNetlistIndexAmongstNetlistsToConnectedComponentIndex)
	{
		findNumberOfConnectedComponentsOf3DDomain();

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

	int getConnectedTopologicalComponentIndexForInnerResistor(const int componentIndex) const;
	void findNumberOfConnectedComponentsOf3DDomain();
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	const double m_oneResistanceToGiveEachResistor;
	const double m_elastanceToGiveCentralCapacitor;
	const double m_initialDomainPressure;
	int m_numberOfConnectedComponentsOf3DDomain;
	const std::map<int,int> m_mapFromNetlistIndexAmongstNetlistsToConnectedComponentIndex;
};

#endif