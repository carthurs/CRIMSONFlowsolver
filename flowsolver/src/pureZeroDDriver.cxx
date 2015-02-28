#include "pureZeroDDriver.hxx"

void PureZeroDDriver::init()
{
	double oneResistanceToGiveEachResistor = 0.001;
	double elastanceToGiveVolumeTrackingPressureChamber = 3000.0;
	double initialDomainPressure = 133.3 * 80; // 80 mmHg
	int negativeIndexForNetlistThatReplaces3DDomain = -1;
	zeroDDomainLPN = boost::unique_ptr<NetlistBoundaryCondition> (new Netlist3DDomainReplacement(negativeIndexForNetlistThatReplaces3DDomain, oneResistanceToGiveEachResistor, elastanceToGiveVolumeTrackingPressureChamber));
}

void PureZeroDDriver::iter_init()
{

}

void PureZeroDDriver::iter_step()
{

}

void PureZeroDDriver::iter_finalize()
{

}

void PureZeroDDriver::finalize()
{

}