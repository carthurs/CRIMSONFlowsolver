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
 
 
 class boundaryCondition
 {
 	friend class boundaryConditionManager;
 protected:
 	double dp_dq;
 	double Hop;
 	double dp_dq_n1;
 	double Hop_n1;
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
public:
 	boundaryCondition(int surfaceIndex_in)
 	{
 	     hstep = inpdat.nstep[0] + timdat.lstep;
 	     delt = inpdat.Delt[0];
 	     flowhist = new double [hstep];
         pressurehist = new double [hstep];
         
 	     surfaceIndex = surfaceIndex_in;
 	     dp_dq = 0.0;
    	 Hop = 0.0;
    	 bcCount++;
    	 index = bcCount;
    	 
 	}

 	virtual double tempDataTestFunction()
 	{
 		return -3.14159265;
 	}

 	~boundaryCondition()
 	{
 	    delete[] flowhist;
 	    delete[] pressurehist;
 	    bcCount--;
 	}
 	virtual void initialiseModel() = 0;
 	double getdp_dq();
 	double getHop();
 	int index;
protected:
	virtual void computeImplicitCoeff_solve(int timestepNumber) = 0;
 	virtual void computeImplicitCoeff_update(int timestepNumber) = 0;
	virtual std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timen_1, double alfi_delt) = 0;
	double linInterpolateTimeData(const std::vector<std::pair<double,double>> &timeData, const double &currentTime, const int timeDataLength);

private:
 	static int bcCount;
 };
 

class boundaryConditionFactory
 {
 public:
	static std::unique_ptr<boundaryCondition> createBoundaryCondition(int surfaceIndex, std::string boundaryType);
 };

class RCR : public boundaryCondition
{
public:
	RCR(int surfaceIndex)
	: boundaryCondition(surfaceIndex)
	{
		initialiseModel();
		Hop = 2.0;

		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		indexOfThisRCR = numberOfInitialisedRCRs;
		numberOfInitialisedRCRs++;
		
		rcrtReader* rcrtReader_instance = rcrtReader::Instance();
		r1 = rcrtReader_instance->getR1()[indexOfThisRCR];
		c = rcrtReader_instance->getC()[indexOfThisRCR];
		r2 = rcrtReader_instance->getR2()[indexOfThisRCR];
		timeDataPdist = rcrtReader_instance->rcrtReader::getTimeDataPdist()[indexOfThisRCR];
		lengthOftimeDataPdist = rcrtReader_instance->rcrtReader::getNumDataRCR()[indexOfThisRCR];
	}
	
	void computeImplicitCoeff_solve(int timestepNumber);
 	void computeImplicitCoeff_update(int timestepNumber);
 	std::pair<double,double> computeImplicitCoefficients(int timestepNumber, double timeAtStepNplus1, double alfi_delt);

	
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
	netlist(int surfaceIndex)
	: boundaryCondition(surfaceIndex)
	{
		initialiseModel();
		Hop = 3.0;
	}
	void computeImplicitCoeff_solve(int timestepNumber)
	{

	}
 	void computeImplicitCoeff_update(int timestepNumber)
 	{

 	}
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
 public:
 
    static boundaryConditionManager* Instance()
	{
		static boundaryConditionManager* instance = new boundaryConditionManager();
		return instance;
	}
	
    void setSurfaceList(std::vector<std::pair<int,std::string>> surfaceList);
    
    void getImplicitCoeff_rcr(double* implicitCoeffs_toBeFilled);
    std::vector<std::unique_ptr<boundaryCondition>>* getBoundaryConditions();

    void computeAllImplicitCoeff_solve(int timestepNumber);
    void computeAllImplicitCoeff_update(int timestepNumber);
 
 private:
    boundaryConditionManager()
    {
    }
    std::vector<std::unique_ptr<boundaryCondition>> boundaryConditions;
    std::map<int,std::pair<double,double>> implicitCoefficientMap;
    boundaryConditionFactory factory;
     
     
 };

void multidom_initialise();
extern "C" void multidom_link(int);
void multidom_iter_initialise();
void multidom_iter_step();
void multidom_iter_finalise();
void multidom_finalise();

#endif /* MULTIDOM_H_ */
