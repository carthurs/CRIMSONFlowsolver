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
		mp_CircuitDescription = boost::shared_ptr<Netlist3DDomainReplacementCircuitData> (new Netlist3DDomainReplacementCircuitData(hstep,m_numberOfNetlistsUsedAsBoundaryConditions));
	}

	void setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<double> boundaryFlowsOrPressuresAsAppropriate);
	std::vector<double> getBoundaryPressures();
	std::vector<double> getBoundaryFlows();

	void solveSystem(const int timestepNumber);

private:
	int m_numberOfNetlistsUsedAsBoundaryConditions;
	void createCircuitDescription();

	const double m_oneResistanceToGiveEachResistor;
	const double m_elastanceToGiveVolumeTrackingPressureChamber;
	const double m_initialDomainPressure;
	std::vector<double> m_boundaryFlowsOrPressuresAsAppropriate;

};

#endif