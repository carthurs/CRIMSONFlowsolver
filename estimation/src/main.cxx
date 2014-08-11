/////////////////////////////////////////////////////////////////
//main function
/////////////////////////////////////////////////////////////////
#include "mpi.h"

#include <iostream>
#include "common_c.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "partition.h"
#include "input_fform.h"

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif

using namespace std;

int main(int argc, char * argv[]) {

   int rank;
   int numProcsTotal,numProcs;
   int ierr = 0;
   char pathToProcsCaseDir[100];

   MPI_Init(&argc,&argv);

   // save the communicator
   MPI_Comm iNewComm_C = MPI_COMM_WORLD;
   newcom.iNewComm = MPI_Comm_c2f(iNewComm_C); // modifies newcom in fortran common block

   MPI_Comm_size(iNewComm_C, &numProcsTotal);
   MPI_Comm_rank(iNewComm_C, &rank);

   if(argc > 2 ){
	   static volatile int debuggerPresent =0;
	   while (!debuggerPresent ); // assign debuggerPresent=1
   }

   if ( argc < 2 ) {
      if ( 0 == rank ) {
         std::cout << "Usage: " << argv[0] << "<.inp file> <optional debug flag> \n";
      }
      return 0;
   }

   // read configuration file
   input_fform();

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
   proces();
   itrdrv_init(); // initialize solver

   for (int kk = 1; kk <= inpdat.nstep[0]; kk++) {
	   itrdrv_iter_init();
	   itrdrv_iter_step();
	   itrdrv_iter_finalize();
   }

   itrdrv_finalize();

   MPI_Finalize();
   return ierr;
}

