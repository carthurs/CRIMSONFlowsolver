#ifndef ABSTRACTBOUNDARYCONDITION_HXX_
#define ABSTRACTBOUNDARYCONDITION_HXX_

#include "gtest/gtest_prod.h"
#include <string>
#include <stdexcept>
#include "common_c.h"

class abstractBoundaryCondition
 {
 	friend class boundaryConditionManager;
 	friend class testMultidom;
 	friend class basicFileWriter;
 	FRIEND_TEST(testMultidom, checkBoundaryConditionsMadeProperly);
 	FRIEND_TEST(testMultidom, checkRCRLinearInterpolators);
 	// FRIEND_TEST(testMultidom, checkDpDqAndHopFortranPasser)
 	FRIEND_TEST(testMultidom, checkImplicitConditionComputation_solve);
 	FRIEND_TEST(testMultidom, checkImplicitConditionComputation_update);
 	FRIEND_TEST(testMultidom, checkFlowAndPressureSetters);
 protected:
 	double dp_dq;
 	double Hop;
 	double dp_dq_n1;
 	double Hop_n1;
 	double pressure_n;
 	double flow_n;
 	int surfaceIndex;
 	int isactive;
 	double* flowhist;
 	double* pressurehist;
 	std::string flowfile;
    std::string pressurefile;
	double surfarea;
	double* flow_n_ptr;
    double flow_n1;
    double* pressure_n_ptr;
    double implicitcoeff;
    double implicitcoeff_n1; 
    int hstep;
    double delt;
    double alfi_local;
    int thisIsARestartedSimulation;

    virtual void computeImplicitCoeff_solve(int timestepNumber) = 0;
 	virtual void computeImplicitCoeff_update(int timestepNumber) = 0;
 	virtual void updpressure_n1_withflow() = 0;
	virtual std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timen_1, double alfi_delt) = 0;
	void updatePressureAndFlowHistory();
	virtual double linInterpolateTimeData(const double &currentTime, const int timeDataLength)
	{
		// std::cout << "Disallowed call to non-overridden (e.g. non-RCR) . Exiting.\n";
    	throw std::runtime_error("Disallowed call to non-overridden (e.g. non-RCR).");
    	return 0.0;
    };
public:
 	abstractBoundaryCondition(int surfaceIndex_in)
 	{
 	     hstep = inpdat.nstep[0] + timdat.lstep;
 	     delt = inpdat.Delt[0];
 	     // allocate arrays with +1 to size, in case hstep=0 (that would be undefined behaviour under new double)
 	     flowhist = new double [hstep+1];
         pressurehist = new double [hstep+1];
         
 	     surfaceIndex = surfaceIndex_in;
 	     dp_dq = 0.0;
    	 Hop = 0.0;
    	 bcCount++;
    	 index = bcCount;

    	 alfi_local = timdat.alfi;

		if (timdat.lstep > 0)
    	{
    		thisIsARestartedSimulation = 1;
    	}
    	else
    	{
    		thisIsARestartedSimulation = 0;
    	}
    	 
 	}

 	virtual double tempDataTestFunction()
 	{
 		return -3.14159265;
 	}

 	virtual ~abstractBoundaryCondition()
 	{
 	    delete[] flowhist;
 	    delete[] pressurehist;
 	    bcCount--;
 	}
 	virtual void initialiseModel() = 0;
 	double getdp_dq();
 	double getHop();
 	int index;
private:
 	static int bcCount;
 };

#endif