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
		mp_CircuitDescription = new Netlist3DDomainReplacementCircuitData(hstep,m_numberOfNetlistsUsedAsBoundaryConditions);
	}

private:
	int m_numberOfNetlistsUsedAsBoundaryConditions;
	void createCircuitDescription();

	const double m_oneResistanceToGiveEachResistor;
	const double m_elastanceToGiveVolumeTrackingPressureChamber;
	const double m_initialDomainPressure;


};

#endif