#ifndef NETLISTZERODDOMAINCIRCUIT_HXX_
#define NETLISTZERODDOMAINCIRCUIT_HXX_

#include "NetlistCircuit.hxx"

class NetlistZeroDDomainCircuit : public NetlistCircuit
{
public:
	NetlistZeroDDomainCircuit(int hstep, const int numberOfNetlistsUsedAsBoundaryConditions, const bool thisIsARestartedSimulation, const double alfi, const double delt, const double oneResistanceToGiveEachResistor, const std::vector<double>& compliancesToGiveCentralCapacitors, const double initialDomainPressure, const std::map<int,int> mapFromZeroDSurfaceIndexToConnectedComponentIndex, const int startingTimestepIndex, const int numberOfDirichletBCTSurfaces)
	: NetlistCircuit(hstep, -1, thisIsARestartedSimulation, alfi, delt, startingTimestepIndex),
	m_oneResistanceToGiveEachResistor(oneResistanceToGiveEachResistor),
	m_elastancesToGiveCentralCapacitors(compliancesToGiveCentralCapacitors),
	m_initialDomainPressure(initialDomainPressure),
	m_numberOfNetlistsUsedAsBoundaryConditions(numberOfNetlistsUsedAsBoundaryConditions),
	m_mapFromZeroDSurfaceIndexToConnectedComponentIndex(mapFromZeroDSurfaceIndexToConnectedComponentIndex),
	m_numberOfDirichletBCTSurfaces(numberOfDirichletBCTSurfaces)
	{
		findNumberOfConnectedComponentsOf3DDomain();

		mp_circuitData = boost::shared_ptr<CircuitData> (new Netlist3DDomainReplacementCircuitData(hstep, numberOfNetlistsUsedAsBoundaryConditions));
		m_PressureHistoryFileName.clear(); // Defensive
		m_PressureHistoryFileName.append("netlistPressures_zeroDDomainReplacement.dat");
		m_FlowHistoryFileName.clear(); // Defensive
		m_FlowHistoryFileName.append("netlistFlows_zeroDDomainReplacement.dat");
		m_VolumeHistoryFileName.clear(); // Defensive
		m_VolumeHistoryFileName.append("netlistVolumes_zeroDDomainReplacement.dat");

		m_numberOfOutlets = m_numberOfNetlistsUsedAsBoundaryConditions + m_numberOfDirichletBCTSurfaces;
	}

	void setBoundaryPrescriptionsAndBoundaryConditionTypes(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate);

	std::vector<double> getBoundaryPressures();
	std::vector<double> getBoundaryFlows();
	void solveSystem();
	void setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients, std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAtBoundaries);
	void createCircuitDescription() override;
	void initialiseAtStartOfTimestep() override;

	boost::shared_ptr<std::vector<std::pair<parameter_controller_t, int>>> getControlTypesAndComponentIndices() const;
private:
	void setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode) override;

	int getConnectedTopologicalComponentIndexForInnerResistor(const int componentIndex) const;
	void findNumberOfConnectedComponentsOf3DDomain();
	const int m_numberOfNetlistsUsedAsBoundaryConditions;
	const double m_oneResistanceToGiveEachResistor;
	const std::vector<double> m_elastancesToGiveCentralCapacitors;
	const double m_initialDomainPressure;
	const int m_numberOfDirichletBCTSurfaces;
    int m_numberOfOutlets;
	int m_numberOfConnectedComponentsOf3DDomain;
	const std::map<int,int> m_mapFromZeroDSurfaceIndexToConnectedComponentIndex;
};

#endif