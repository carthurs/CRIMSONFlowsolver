#include "abstractBoundaryCondition.hxx"


// Statics
int abstractBoundaryCondition::bcCount = 0;

// initialise the multidomain/LPN objects, this will need IFDEF for 3D and 1D codes
double abstractBoundaryCondition::getHop()
{
  return Hop;
}

double abstractBoundaryCondition::getdp_dq()
{
  return dp_dq;
}

void abstractBoundaryCondition::setLPNInflowPressure(double inflowPressure)
{
	// This will not be needed by all subclasses!
    LPNInflowPressure = inflowPressure;
}

void abstractBoundaryCondition::updatePressureAndFlowHistory()
{
  pressurehist[timdat.lstep] = pressure_n;
  flowhist[timdat.lstep] = *flow_n_ptr;
}