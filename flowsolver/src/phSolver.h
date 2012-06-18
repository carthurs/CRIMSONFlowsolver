/****************************************************************************** 

(c) 2009-2010 Scientific Computation Research Center, 
Rensselaer Polytechnic Institute. All rights reserved.

The LICENSE-SCOREC file included with this distribution describes the terms
of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/

/** \mainpage PHASTA Solver API
 *
 * @file      phSolver.h
 * @date      August 1, 2011
 * @version   0.1
 *
 * @brief     class for interfacing with phSolver via a functional interface
 *
 * @remark    class definition supporting interactions with phSolver via a functional interface and without files between calls solve() routine
 *
 */

#ifndef _PHSOLVER_H_
#define _PHSOLVER_H_

#include <errno.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <cstring>

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif
#include "common_c.h" //phSolver internal functions and global data - should this be in phSolver.h ???
#include "CInput.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "mpi.h"




using namespace std;
using namespace PHSOLVER;

/**
 * @class     phSolver
 * @brief     class for interfacing with phSolver via a functional interface
 * @remark    supports initialization of a single solver object, loading data into internal data structures from file or memory, and running the solver
 */
class phSolver {
public:

	/**
	 * @brief     initialize the phSolver instance
	 * @remark    initializes phSolver object, assumes solver.inp and input.config files are in the runtime current working directory, assumes MPI communicator is MPI_COMM_WORLD
	 *
	 * @param     none
	 * @return    none
	 */
	static phSolver* Instance();


	void setCommunicator(MPI_Comm comm);

	MPI_Comm getCommunicator();

	void readConfiguration();
	int readMeshAndSolution_fromFiles(char* pathToProcsCaseDir);

	int Solve();
	int SolverInit();

	void SolverForwardInit();
	void SolverForwardStep();
	void SolverForwardFinalize();

	void SolverFinalize();

	int getNumProcs();
	int getNumDofHolders();
	int getNumTimeSteps();

	int getTime();

	void setTime(int timeIn);

	int getSize();

	bool hasFinished();

    int input_fform();

	~phSolver();
private:

	char pc_pathToProcsCaseDir_[1024];

	/**
	 * @brief     default constructor
	 * @remark    initializes phSolver object, assumes solver.inp and input.config files are in the runtime current working directory, assumes MPI communicator is MPI_COMM_WORLD
	 *
	 * @param     none
	 * @return    none
	 */
	phSolver();
	static phSolver* p_Instance_;

	MPI_Comm iNewComm_C;

	int rank_;
	int numProcs_;

}; //end class phSolver

#endif

