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
 protected:
 	double dp_dq;
 	double Hop;
 	int surfaceIndex;
 	int isactive;
 	double* flowhist;
 	double* pressurehist;
 	std::string flowfile;
    std::string pressurefile;
	double surfarea;
	double flow_n;
    double flow_n1;
    double pressure_n;
    double implicitcoeff;
    double implicitcoeff_n1; 
    int hstep;
public:
 	boundaryCondition(int surfaceIndex_in)
 	{
 	     hstep = inpdat.nstep[0] + timdat.lstep;
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
 	virtual void updateImplicitCoefficients() = 0;
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
	static std::unique_ptr<boundaryCondition> createBoundaryCondition(int surfaceIndex, std::string boundaryType);
 };

class RCR : public boundaryCondition
{
public:
	RCR(int surfaceIndex)
	: boundaryCondition(surfaceIndex)
	{
		initialiseModel();
		updateImplicitCoefficients();
		Hop = 2.0;

		// Note the index of this RCR (zero-indexed), and count its existance
		// (the order of these two lines is correct!)
		indexOfThisRCR = numberOfInitialisedRCRs;
		numberOfInitialisedRCRs++;
		
		rcrtReader* rcrtReader_instance = rcrtReader::Instance();
		r1=rcrtReader_instance->getR1()[indexOfThisRCR];
		c=rcrtReader_instance->getC()[indexOfThisRCR];
		r2=rcrtReader_instance->getR2()[indexOfThisRCR];
	}

	double tempDataTestFunction()
	{
		return r1;
	}

	~RCR()
	{
		numberOfInitialisedRCRs--;
	}

	void updateImplicitCoefficients()
	{
		// Do the RCR updating
		std::cout <<  "RCR updating..." << std::endl;
	}
private:
	void initialiseModel();
	static int numberOfInitialisedRCRs;
	int indexOfThisRCR;
	int pdmax;
	double r1; // Proximal resistance
	double c; // Capacitance (compliance)
	double r2; // Distal Resistance
	std::vector<std::pair<double,double>> pdval;
};

class netlist : public boundaryCondition
{
public:
	netlist(int surfaceIndex)
	: boundaryCondition(surfaceIndex)
	{
		initialiseModel();
		updateImplicitCoefficients();
		Hop = 3.0;
	}
	void updateImplicitCoefficients()
	{
		// Do the netlist updating
		std::cout <<  "netlist updating..." << std::endl;
	}
	void initialiseModel()
	{
		std::cout << "netlist Initialisation" << std::endl;
	}

};

void multidom_initialise();
extern "C" void multidom_link(int);
void multidom_iter_initialise();
void multidom_iter_step();
void multidom_iter_finalise();
void multidom_finalise();

#endif /* MULTIDOM_H_ */
