#ifndef CLOSEDLOOPBOUNDARYCONDITIONSUBSECTION_HXX_
#define CLOSEDLOOPBOUNDARYCONDITIONSUBSECTION_HXX_
#include "NetlistBoundaryCondition.hxx"

class ClosedLoopBoundaryConditionSubsection : public NetlistBoundaryCondition
{
public:
	ClosedLoopBoundaryConditionSubsection(const int surfaceIndex_in, const double hstep_in, const double delt_in, const double alfi_in, const double lstep, const int maxsurf, const int nstep)
	: NetlistBoundaryCondition(surfaceIndex_in, hstep_in, delt_in, alfi_in, lstep, maxsurf, nstep)
	{

	}
protected:
private:
};

#endif