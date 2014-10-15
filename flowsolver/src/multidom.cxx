/*

 *
 *  Created on: Oct 7, 2014
 *      Author: klau
 */

#include "multidom.h"
#include <iostream>

// initialise the multidomain/LPN objects, this will need IFDEF for 3D and 1D codes

void multidom_initialise(){

	std::cout << "HELLO WORLD" << std::endl;

	if (grcrbccom.numGRCRSrfs > 0)
	{
		for (int i=0; i<grcrbccom.numGRCRSrfs; i++)
		{
			std::cout << grcrbccom.nsrflistGRCR[i+1] << std::endl;
		}
	}

}

//
extern "C" void multidom_link(int num)
{

	std::cout << "Index = " << num << std::endl;	
	    
}

void multidom_iter_initialise(){

}

void multidom_iter_step(){

}

void multidom_iter_finalise(){

}

void multidom_finalise(){

}


// set pointer to fortran arrays



ReducedOrderModel::ReducedOrderModel()
: surfid(0)
{

}

ReducedOrderModel::~ReducedOrderModel()
{
}

