#ifndef NETLISTSUBCIRCUIT_HXX_
#define NETLISTSUBCIRCUIT_HXX_

#include "datatypesInCpp.hxx"
#include "petscsys.h"
#include "petscmat.h"
#include "petscvec.h"
#include <vector>
#include <set>
#include "CircuitData.hxx"
#include "indexShifters.hxx"

class NetlistSubcircuit
{
public:
	NetlistSubcircuit(boost::shared_ptr<CircuitData> circuitData, const std::vector<double*> flow_n_ptrs, const std::vector<double*> pressure_n_ptrs, const int surfaceIndex_in, const double delt)
	: 
	 mp_circuitData(circuitData),
	 flow_n_ptrs(flow_n_ptrs),
	 pressure_n_ptrs(pressure_n_ptrs),
	 surfaceIndex(surfaceIndex_in),
	 m_delt(delt)
	{
		initialiseSubcircuit();
		// columnIndexOf3DInterfaceFlowInLinearSystem = 0;
		safetyCounterLimit = 1000;
	}

	~NetlistSubcircuit()
	{
		
	}

	
protected:

private:
	

};

#endif