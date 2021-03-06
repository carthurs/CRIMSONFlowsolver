/////////////////////////////////////////////////////////////////
//main function - builds flowsolver binary
/////////////////////////////////////////////////////////////////
#include <Python.h>
#include "CRIMSONPython.hxx"
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
#include "FortranBoundaryDataPointerManager.hxx"
#include "autoGeneratedVersionNumber.hxx"
#include "dateTools.hxx"
#include "petscsys.h"
#include "CrimsonGlobalArrayTransfer.h"
#include "itrPC.h"
#include "pureZeroDDriver.hxx"

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;

static char help[] = "Pure Flowsolver.\n\n";

int main(int argc, char **argv) {
   std::cout << "FLOWSOLVER STARTING" << std::endl;

   kalmanFilterActive.kalmanFilterOn = false;
   int rank;
   int numProcsTotal;
   int ierr = 0;
   char pathToProcsCaseDir[100];

   // MPI_Init(&argc,&argv);
   PetscInitialize(&argc, &argv, (char *)0, help);
   
   initialisePython();  

   // save the communicator
   MPI_Comm iNewComm_C = MPI_COMM_WORLD;
   newcom.iNewComm = MPI_Comm_c2f(iNewComm_C); // modifies newcom in fortran common block

   MPI_Comm_size(iNewComm_C, &numProcsTotal);
   MPI_Comm_rank(iNewComm_C, &rank);

   char buildNumber[100];
   char buildTime[100];
   getBuildNumber(buildNumber);
   getBuildTime(buildTime);
   if (rank == 0)
   {
    std::cout << "This is CRIMSON Flowsolver version " << buildNumber << ", built at " << buildTime << "." << std::endl;
   }

   // Expiry date check (uncomment enableExpiryDate() call below to enable):
   expiryDate expiry = expiryDate();
   expiry.setExpiryDayOfMonth(31);
   expiry.setExpiryMonthOfYear(12);
   expiry.setExpiryYear(2021);
   // UNCOMMENT TO DO A BUILD WITH AN EXPIRY DATE!
   // expiry.enableExpiryDate();
   expiry.checkWhetherExpiryDatePassed();

   if ( argc < 2 ) {
      if ( 0 == rank ) {
         std::cout << "Usage: " << argv[0] << "<.inp file> <optional debug flag> \n";
      }
      return 0;
   }

   // Debugger snare:
   if (rank==0)
   {
     const char* debuggerFlag = "1";
     for(int ii=1; ii<argc; ii++)
     {
       // Look for a single "1" on the command line, indicating that we should
       // wait for the debugger...
       if(!strcmp(argv[ii], debuggerFlag))
       {
           static volatile int debuggerPresent =0;
           std::cout << "Debug flag spotted on the command line. Pausing to await debugger connection..." << std::endl;
           while (!debuggerPresent )
           {
             sleep(2);
             std::cout << "Awaiting debugger connection..." << std::endl;
           }; // assign debuggerPresent=1
       }
     }
   }

   // wait for the debugger if to release MPI process with rank 0, if a debug flag was present:
   MPI_Barrier(MPI_COMM_WORLD);

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
   proces(); // Includes the call to set up FortranBoundaryDataPointerManager pointerManager

   {
       // just initialise the time values that the AbstractBoundaryCondition needs (when it's called in multidom_initialise).
       // This will be called again during itrdrv_init.
       // This is really ugly, but a proper fix will take days - it's a BIG refactor.
     int dummyInitialItseqValue=1;
     callFortranSetupTimeParameters(dummyInitialItseqValue);
   }

   // initialise reduced order boundary conditions
   multidom_initialise();

   if (nomodule.pureZeroDSimulation == 0)
   {
     multidomSetupControlSystems();
     itrdrv_init(); // initialize solver

     // FortranBoundaryDataPointerManager* pointerManager;
     // pointerManager = FortranBoundaryDataPointerManager::Get();

     for (int kk = 1; kk <= inpdat.nstep[0]; kk++) {
      multidom_iter_initialise();
  	   itrdrv_iter_init();

  	   itrdrv_iter_step();

  	   itrdrv_iter_finalize();
       multidom_iter_finalise();
     }

     itrdrv_finalize();
   }
   else // we disable to 3D domain, and just connect the boundary conditions to each other via a simple 0D circuit. This is for rapid prototyping and experimentation
   {
      assert(nomodule.pureZeroDSimulation==1);

      // numDirCalcSrfs = number of prescribed flow (i.e. bct; dirichlet) surfaces
      std::vector<double> zeroDDomainCompliancesForEachConnectedComponent;
      double* const zeroDDomainCompliancesEndPointer = nomodule.zeroDDomainCompliances + nomodule.num3DConnectedComponents + 1; // Fortran indexing
      zeroDDomainCompliancesForEachConnectedComponent.assign(nomodule.zeroDDomainCompliances + 1, zeroDDomainCompliancesEndPointer);

      PureZeroDDriver pureZeroDDriver(nomodule.numDirCalcSrfs, zeroDDomainCompliancesForEachConnectedComponent);

      pureZeroDDriver.setDelt(inpdat.Delt[0]);
      pureZeroDDriver.setAlfi(timdat.alfi);
      pureZeroDDriver.setHstep(inpdat.nstep[0] + timdat.currentTimestepIndex);
      pureZeroDDriver.setNtout(outpar.ntout);
      // nsrflistDirCalc = list of prescribed flow surface indices
      pureZeroDDriver.setupConnectedComponents(nomodule.num3DConnectedComponents, nomodule.surfacesOfEachConnectedComponent, nomodule.indicesOfNetlistSurfaces, nomodule.nsrflistDirCalc);

      pureZeroDDriver.init();
      multidomSetupControlSystems();

      // pointer manager?      

      for (int kk = 1; kk <= inpdat.nstep[0]; kk++)
      {
        pureZeroDDriver.iter_init();
        pureZeroDDriver.iter_step();
        pureZeroDDriver.iter_finalize();
        multidom_iter_finalise();
      }

      pureZeroDDriver.finalize();

   }

   multidom_finalise();

   CrimsonGlobalArrayTransfer::Get()->tearDown();

   // MPI_Finalize();
   Py_Finalize();
   ierr = PetscFinalize();
   CHKERRQ(ierr);
   return ierr;
}
