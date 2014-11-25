#ifndef NETLISTBOUNDARYCONDITION_HXX_
#define NETLISTBOUNDARYCONDITION_HXX_

#include "abstractBoundaryCondition.hxx"

class netlist : public abstractBoundaryCondition
{
public:
	netlist(int surfaceIndex_in)
	: abstractBoundaryCondition(surfaceIndex_in)
	{
		initialiseModel();
	}
	void computeImplicitCoeff_solve(int timestepNumber)
	{

	}
 	void computeImplicitCoeff_update(int timestepNumber)
 	{

 	}
 	void updpressure_n1_withflow(){}
 	std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timen_1, double alfi_delt)
 	{
 		std::pair<double,double> dummyValue;
 		dummyValue.first=-3.14;
 		dummyValue.second=-2.718281828;
 		return dummyValue;
 	}
	void initialiseModel()
	{
		std::cout << "netlist Initialisation" << std::endl;
	}

};

#endif