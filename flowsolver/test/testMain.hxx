#ifndef TESTMAIN_HXX_
#define TESTMAIN_HXX_

#include "gtest/gtest.h"

/////////////////////////////////////////////////////////////////
//main function - builds flowsolver binary
/////////////////////////////////////////////////////////////////
#include "mpi.h"

#include <iostream>
#include <stdio.h>
#include "common_c.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "partition.h"
#include "input_fform.h"
#include "multidom.h"
#include "fortranPointerManager.hxx"
#include "fileReaders.hxx"
#include "boundaryConditionManager.hxx"

#include "debuggingToolsForCpp.hxx"

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif



// The fixture for testing class Foo.
	class testMain : public ::testing::Test {
	 protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.
	 	std::string dirBinaryCalledFrom;

	 	// This constructor should just do exactly what main.cxx does in estimation/src/main.cxx
		testMain() {
			dirBinaryCalledFrom = get_current_dir_name();
		}

		void setSimDirectory(std::string dir)
		{
			chdir(dir.c_str());
		}

		void clearOutOldFiles()
		{
			// Warning - this has the potential to delete multiple folders due to the *
			system("rm -rf *-procs-case");
			MPI_Barrier(MPI_COMM_WORLD);
		}

		void runSimulation()
		{
		   int rank;
		   int numProcsTotal,numProcs;
		   int ierr = 0;
		   char pathToProcsCaseDir[100];

		   // Moved this to the gtest_main.cc
		   // MPI_Init(&fake_argc,(char***)&fake_argv);

		   // save the communicator
		   MPI_Comm iNewComm_C = MPI_COMM_WORLD;
		   newcom.iNewComm = MPI_Comm_c2f(iNewComm_C); // modifies newcom in fortran common block

		   MPI_Comm_size(iNewComm_C, &numProcsTotal);
		   MPI_Comm_rank(iNewComm_C, &rank);

		   // Moved this to the gtest_main.cc
		   // if(fake_argc > 2 ){
			  //  static volatile int debuggerPresent =0;
			  //  while (!debuggerPresent ); // assign debuggerPresent=1
		   // }

		   // Dont need this during tests:
		   // if ( argc < 2 ) {
		   //    if ( 0 == rank ) {
		   //       std::cout << "Usage: " << fake_argv[0] << "<.inp file> <optional debug flag> \n";
		   //    }
		   //    // return 0;
		   // }

		   // read configuration file
		   int errFlag = input_fform();
		   if (errFlag != 0)
		   {
		   	   throw std::runtime_error("EE: Failed during parsing of input files.");
		   }

		   // Preprocess data and run the problem
		   // Partition the problem to the correct number of processors

		   if( rank == 0 )
		   {
		      //cout << "number of procs " << numprocs_perparticle << endl;
		      Partition_Problem( numProcsTotal );
		   }

		   MPI_Barrier(MPI_COMM_WORLD);

		   sprintf(pathToProcsCaseDir,"%d-procs-case",numProcsTotal);
		   chdir(pathToProcsCaseDir);
		   //sprintf(inpfilename,"%d-procs-case",numprocs_perparticle);

		   std::cout << "changing directory to " << pathToProcsCaseDir << std::endl;

		   input(&numProcsTotal, &rank);
		   proces();

		   itrdrv_init(); // initialize solver

		   // initialise reduced order boundary conditions
		   multidom_initialise();

		   fortranBoundaryDataPointerManager* pointerManager;
		   pointerManager = fortranBoundaryDataPointerManager::Get();
		   for (int kk = 1; kk <= inpdat.nstep[0]; kk++) {
			   itrdrv_iter_init();
			   itrdrv_iter_step();
			   itrdrv_iter_finalize();
		   }

		   itrdrv_finalize();
		   multidom_finalise();
		   // Moved this to the gtest_main.cc
		   // MPI_Finalize();


		   // return ierr;
		}

	  virtual ~testMain() {
	    // You can do clean-up work that doesn't throw exceptions here.
	    chdir(dirBinaryCalledFrom.c_str());
	  }

	  // If the constructor and destructor are not enough for setting up
	  // and cleaning up each test, you can define the following methods:

	  virtual void SetUp() {
	    // Code here will be called immediately after the constructor (right
	    // before each test).
	  }

	  virtual void TearDown() {
	    // Code here will be called immediately after each test (right
	    // before the destructor).
	    fortranBoundaryDataPointerManager::Get()->tearDown();
	    // SimvascularGlobalArrayTransfer::Get()->tearDown();
	    // boundaryConditionManager::Instance()->Term();
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibraryTestMain();

#endif