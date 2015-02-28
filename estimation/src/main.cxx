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
#include "multidom.hxx"
#include "fortranPointerManager.hxx"
#include "autoGeneratedVersionNumber.hxx"
#include "dateTools.hxx"
#include "petscsys.h"
#include "SimvascularGlobalArrayTransfer.h"
#include "itrPC.h"

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;

static char help[] = "Pure Flowsolver.\n\n";

int main(int argc, char * argv[]) {


   int rank;
   int numProcsTotal,numProcs;
   int ierr = 0;
   char pathToProcsCaseDir[100];

   // MPI_Init(&argc,&argv);
   PetscInitialize(&argc, &argv, (char *)0, help);

   // save the communicator
   MPI_Comm iNewComm_C = MPI_COMM_WORLD;
   newcom.iNewComm = MPI_Comm_c2f(iNewComm_C); // modifies newcom in fortran common block

   MPI_Comm_size(iNewComm_C, &numProcsTotal);
   MPI_Comm_rank(iNewComm_C, &rank);

   char buildNumber[100];
   char buildTime[100];
   getBuildNumber(buildNumber);
   getBuildTime(buildTime);
   std::cout << "This is Simvascular version " << buildNumber << ", built at " << buildTime << "." << std::endl;

   // Expiry date check (uncomment enableExpiryDate() call below to enable):
   expiryDate expiry = expiryDate();
   expiry.setExpiryDayOfMonth(14);
   expiry.setExpiryMonthOfYear(11);
   expiry.setExpiryYear(2014);
   // UNCOMMENT TO DO A BUILD WITH AN EXPIRY DATE!
   // expiry.enableExpiryDate();
   expiry.checkWhetherExpiryDatePassed();


   // Debugger snare:
   const char* debuggerFlag = "1";
   for(int ii=1; ii<argc; ii++)
   {
     // Look for a single "1" on the command line, indicating that we should
     // wait for the debugger...
     if(!strcmp(argv[ii], debuggerFlag))
     {
         static volatile int debuggerPresent =0;
         std::cout << "Debug flag spotted on the command line. Pausing to await debugger connection..." << std::endl;
         while (!debuggerPresent ); // assign debuggerPresent=1
     }
   }

   if ( argc < 2 ) {
      if ( 0 == rank ) {
         std::cout << "Usage: " << argv[0] << "<.inp file> <optional debug flag> \n";
      }
      return 0;
   }

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

   cout << "changing directory to " << pathToProcsCaseDir << endl;

   input(&numProcsTotal, &rank);
   proces(); // Includes the call to set up fortranBoundaryDataPointerManager pointerManager

   {
       // just initialise the time values that the abstractBoundaryCondition needs (when it's called in multidom_initialise).
       // This will be called again during itrdrv_init.
       // This is really ugly, but a proper fix will take days - it's a BIG refactor.
     int dummyInitialItseqValue=1;
     callFortranSetupTimeParameters(dummyInitialItseqValue);
   }

   // initialise reduced order boundary conditions
   multidom_initialise();

   if (nomodule.pureZeroDSimulation == 0)
   {
     itrdrv_init(); // initialize solver

     fortranBoundaryDataPointerManager* pointerManager;
     pointerManager = fortranBoundaryDataPointerManager::Get();

     for (int kk = 1; kk <= inpdat.nstep[0]; kk++) {
  	   itrdrv_iter_init();
        // std::cout << "3: C++ saw flow in dereferenced pointer: " << *(pointerManager->boundaryFlows.at(3)) << std::endl;
        //      // see elmgmr.f90 line 379 (approx) for the code I added (the call with "pressure") to make this update.....:
        //      std::cout << "3: C++ saw pressure in dereferenced pointer: " << *(pointerManager->boundaryPressures.at(3)) << std::endl;
  	   itrdrv_iter_step();
        // std::cout << "2: C++ saw flow in dereferenced pointer: " << *(pointerManager->boundaryFlows.at(3)) << std::endl;
        //      // see elmgmr.f90 line 379 (approx) for the code I added (the call with "pressure") to make this update.....:
        //      std::cout << "2: C++ saw pressure in dereferenced pointer: " << *(pointerManager->boundaryPressures.at(3)) << std::endl;
  	   itrdrv_iter_finalize();
             // std::cout << "1: C++ saw flow in dereferenced pointer: " << *(pointerManager->boundaryFlows.at(3)) << std::endl;
             // // see elmgmr.f90 line 379 (approx) for the code I added (the call with "pressure") to make this update.....:
             // std::cout << "1: C++ saw pressure in dereferenced pointer: " << *(pointerManager->boundaryPressures.at(3)) << std::endl;
     }

     itrdrv_finalize();
   }
   else // we disable to 3D domain, and just connect the boundary conditions to each other via a simple 0D circuit. This is for rapid prototyping and experimentation
   {
      assert(nomodule.pureZeroDSimulation==1);

      PureZeroDDriver pureZeroDDriver();

      pureZeroDDriver.init();

      // pointer manager?      

      for (int kk = 1; kk <= inpdat.nstep[0]; kk++)
      {
        pureZeroDDriver.iter_init();
        pureZeroDDriver.iter_step();
        pureZeroDDriver.iter_finalize();
      }

      pureZeroDDriver.finalize();

   }

   multidom_finalise();

   SimvascularGlobalArrayTransfer::Get()->tearDown();

   // MPI_Finalize();
   ierr = PetscFinalize();
   CHKERRQ(ierr);
   return ierr;
}

