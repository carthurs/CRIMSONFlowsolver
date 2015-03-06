#ifndef NETLIST3DDOMAINREPLACEMENT_HXX_
#define NETLIST3DDOMAINREPLACEMENT_HXX_
#include "NetlistBoundaryCondition.hxx"

class Netlist3DDomainReplacement : public NetlistBoundaryCondition
{
public:
	Netlist3DDomainReplacement(int surfaceIndex_in, const double m_oneResistanceToGiveEachResistor_in, const double m_elastanceToGiveVolumeTrackingPressureChamber_in, const double m_initialDomainPressure_in)
	: NetlistBoundaryCondition(surfaceIndex_in),
	m_oneResistanceToGiveEachResistor(m_oneResistanceToGiveEachResistor_in),
	m_elastanceToGiveVolumeTrackingPressureChamber(m_elastanceToGiveVolumeTrackingPressureChamber_in),
	m_initialDomainPressure(m_initialDomainPressure_in)
	{
	}

	void setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate);
	std::vector<double> getBoundaryPressures();
	std::vector<double> getBoundaryFlows();

	void solveSystem(const int timestepNumber);

	void setPointersToBoundaryPressuresAndFlows(double* const mp_interfacePressuresToBeReadBy3DDomainReplacement, double* const mp_interfaceFlowsToBeReadBy3DDomainReplacement, const int& numberOfPointers);

	void createCircuitDescription();

	void initialiseModel();

	void setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients);

private:
	int m_numberOfNetlistsUsedAsBoundaryConditions;

	void selectAndBuildActiveSubcircuits();

	const double m_oneResistanceToGiveEachResistor;
	const double m_elastanceToGiveVolumeTrackingPressureChamber;
	const double m_initialDomainPressure;
	std::vector<std::pair<boundary_data_t,double>> m_boundaryFlowsOrPressuresAsAppropriate;
	void setupPressureNode(const int indexOfNodeInInputData, boost::shared_ptr<CircuitPressureNode>& node, boost::shared_ptr<CircuitComponent> componentNeighbouringThisNode);

};

#endif