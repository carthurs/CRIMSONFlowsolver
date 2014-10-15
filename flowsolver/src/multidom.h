/*
 * multidom.h
 *
 *  Created on: Oct 7, 2014
 *      Author: klau
 */

#ifndef MULTIDOM_H_
#define MULTIDOM_H_

// include the common block global variables
#include "common_c.h"

//

void multidom_initialise();
extern "C" void multidom_link(int);
void multidom_iter_initialise();
void multidom_iter_step();
void multidom_iter_finalise();
void multidom_finalise();



class ReducedOrderModel
{

public:

protected:

	int surfid;

public:

	// constructor
	ReducedOrderModel();

	// destructor
	~ReducedOrderModel();

	//
	int getSurfID();

};


#endif /* MULTIDOM_H_ */
