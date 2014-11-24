/*
 * multidom.h
 *
 *  Created on: Oct 7, 2014
 *      Author: klau, carthurs
 */

#ifndef MULTIDOM_H_
#define MULTIDOM_H_

// include the common block global variables
#include "common_c.h"
#include "fileReaders.hxx"
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <stdexcept>
#include "gtest/gtest_prod.h"
#include "debuggingToolsForCpp.hxx"
#include <boost/shared_ptr.hpp>
 
 
 class boundaryCondition
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
 	boundaryCondition(int surfaceIndex_in)
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

 	virtual ~boundaryCondition()
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
 

class boundaryConditionFactory
 {
 public:
	static boost::shared_ptr<boundaryCondition> createBoundaryCondition(int surfaceIndex_in, std::string boundaryType);
 };

class RCR : public boundaryCondition
{
public:
	RCR(int surfaceIndex_in)
	: boundaryCondition(surfaceIndex_in)
	{
		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		indexOfThisRCR = numberOfInitialisedRCRs;
		numberOfInitialisedRCRs++;

		initialiseModel();
		
		rcrtReader* rcrtReader_instance = rcrtReader::Instance();
		r1 = rcrtReader_instance->getR1()[indexOfThisRCR];
		c = rcrtReader_instance->getC()[indexOfThisRCR];
		r2 = rcrtReader_instance->getR2()[indexOfThisRCR];
		timeDataPdist = rcrtReader_instance->getTimeDataPdist()[indexOfThisRCR];
		lengthOftimeDataPdist = rcrtReader_instance->getNumDataRCR()[indexOfThisRCR];
	}
	
	void computeImplicitCoeff_solve(int timestepNumber);
 	void computeImplicitCoeff_update(int timestepNumber);
 	std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timeAtStepNplus1, double alfi_delt);
 	void updpressure_n1_withflow();

	
//  	procedure :: setimplicitcoeff_rcr => setimplicitcoeff_rcr
        // void setImplicitCoeff();
//     procedure :: updxvars_rcr => updxvars_rcr         
//     procedure :: writexvars_rcr => writexvars_rcr
//     procedure :: assign_ptrs_ext_rcr => assign_ptrs_ext_rcr

	double tempDataTestFunction()
	{
		return r1;
	}

	~RCR()
	{
		numberOfInitialisedRCRs--;
	}
protected:
	double linInterpolateTimeData(const double &currentTime, const int timeDataLength);

private:
	void initialiseModel();
	static int numberOfInitialisedRCRs;
	int indexOfThisRCR;
	int pdmax;
	double r1; // Proximal resistance
	double c; // Capacitance (compliance)
	double r2; // Distal Resistance
	std::vector<std::pair<double,double>> timeDataPdist; // Time-varying disal pressure data
	int lengthOftimeDataPdist;
};

class netlist : public boundaryCondition
{
public:
	netlist(int surfaceIndex_in)
	: boundaryCondition(surfaceIndex_in)
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

 class boundaryConditionManager
 {
 	friend class testMultidom;
 public:

 	static histFileReader* PHistReader;
 
    static boundaryConditionManager* Instance()
	{
		if (!instance) {
			instance = new boundaryConditionManager();
    		instance->ifRestartingLoadNecessaryData();
		}
		return instance;
	}

	static void Term()
	{
		if (thisIsARestartedSimulation)
		{
			delete boundaryConditionManager::PHistReader;
		}
		delete instance;
		instance = 0;
	}
	
    void setSurfaceList(std::vector<std::pair<int,std::string>> surfaceList);

    void ifRestartingLoadNecessaryData();
    
    void getImplicitCoeff_rcr(double* implicitCoeffs_toBeFilled);
    std::vector<boost::shared_ptr<boundaryCondition>>* getBoundaryConditions();

    void computeAllImplicitCoeff_solve(int timestepNumber);
    void computeAllImplicitCoeff_update(int timestepNumber);

    void updateAllRCRS_Pressure_n1_withflow();
    void updateAllRCRS_setflow_n(double* flows);
    void updateAllRCRS_setflow_n1(double* flows);

    void recordPressuresAndFlowsInHistoryArrays();

    void writePHistAndQHistRCR();

    ~boundaryConditionManager()
    {
    }
 
 private:
    boundaryConditionManager()
    {
    	numberOfRCRSurfaces = grcrbccom.numGRCRSrfs;
    	if (timdat.lstep > 0)
    	{
    		thisIsARestartedSimulation = 1;
    	}
    	else
    	{
    		thisIsARestartedSimulation = 0;
    	}
    }
    std::vector<boost::shared_ptr<boundaryCondition>> boundaryConditions;
    std::map<int,std::pair<double,double>> implicitCoefficientMap;
    boundaryConditionFactory factory;

    static int thisIsARestartedSimulation;
    static int numberOfRCRSurfaces;

    // For testing purposes, to clear the static class out before the next test begins
    // Note that you'll have to make the test class a friend in order to use this..
    // I've made it private on purpose!
    void tearDown()
    {
		boundaryConditionManager::Instance()->Term();
    }

    static boundaryConditionManager* instance;
 };

void multidom_initialise();
extern "C" void multidom_link(int);
void multidom_iter_initialise();
void multidom_iter_step();
void multidom_iter_finalise();
void multidom_finalise();

#endif /* MULTIDOM_H_ */
