/*
 * multidom.hxx
 *
 *  Created on: Oct 7, 2014
 *      Author: klau, carthurs
 */

#ifndef MULTIDOM_H_
#define MULTIDOM_H_

// include the common block global variables
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <fstream>
#include "debuggingToolsForCpp.hxx"

void multidom_initialise();
extern "C" void multidom_link(int);
void multidom_iter_initialise();
void multidom_iter_step();
void multidom_iter_finalise();
void multidom_finalise();

#endif /* MULTIDOM_H_ */
