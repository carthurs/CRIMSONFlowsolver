/////////////////////////////////////////////////////////////////
//main function
/////////////////////////////////////////////////////////////////
#include "mpi.h"

#include <iostream>

#include "partition.h"
#include "phSolver.h"
#include "SCField.h"

using namespace std;

int main(int argc, char * argv[]) {

   int rank;
   int numProcsTotal,numProcs;
   int ierr = 0;
   char inpfilename[100];

//   MPI_Comm newcomm;
//   int color,key;
//   int numparticles = 1;
//   int numprocs_perparticle;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &numProcsTotal);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


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

   // create MPI communicator for one simulation
   // by splitting MPI_COMM_WORLD
//   numprocs_perparticle = numProcsTotal / numparticles;
//   color = rank / numprocs_perparticle;
//   key = 0;
//   MPI_Comm_split(MPI_COMM_WORLD,color,key,&newcomm);

   // Initialize phSolver 
   phSolver* phS = phSolver::Instance();

   // save the communicator
   //phS->setCommunicator(newcomm);
   phS->setCommunicator(MPI_COMM_WORLD);

   // read configuration file
   phS->readConfiguration();

   // update numProcsTotal and rank
//   MPI_Comm_size(newcomm, &numProcs);
//   MPI_Comm_rank(newcomm, &rank);

   // Preprocess data and run the problem
   // Partition the problem to the correct number of processors

   if( rank == 0 )
   {
      //cout << "number of procs " << numprocs_perparticle << endl;
      Partition_Problem( numProcsTotal );
   }

   MPI_Barrier(MPI_COMM_WORLD);

   sprintf(inpfilename,"%d-procs-case",numProcsTotal);
   //sprintf(inpfilename,"%d-procs-case",numprocs_perparticle);

   cout << "changing directory to " << inpfilename << endl;

   phS->readMeshAndSolution_fromFiles(inpfilename);

   int solveReturn = phS->SolverInit(); // initialize solver
   if ( 0 == solveReturn ) {

      for (int kk = 1; kk <= phS->getNumTimeSteps(); kk++) {
         phS->SolverForwardInit();
         phS->SolverForwardStep();
         phS->SolverForwardFinalize();
      }

      phS->SolverFinalize();
   }
   else {
      if (rank == 0)
         fprintf(stderr, "Solve failed ... exiting\n");
      ierr = 1;
   }

   //MPI_Comm_free(&newcomm);
   MPI_Finalize();
   return ierr;
}

