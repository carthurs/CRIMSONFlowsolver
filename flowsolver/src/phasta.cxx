/*********************************************************************

Copyright (c) 2000-2007, Stanford University, 
    Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
    Charles A. Taylor (see SimVascular Acknowledgements file 
    for additional contributors to the source code).

All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions 
are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution. 
Neither the name of the Stanford University or Rensselaer Polytechnic
Institute nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior 
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

**********************************************************************/
#include "mpi.h"

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include <iostream>

#include "partition.h"
#include "input.h"
#include "proces.h"
#include "itrdrv.h"
#include "input_fform.h"

#ifdef intel
#include <direct.h>
//#define input INPUT
//#define proces PROCES
#define chdir _chdir
#else
#include <unistd.h>
#ifndef ibm6000
//#define input input_
//#define proces proces_
#define timer timer_
#endif
#endif

extern "C" char phasta_iotype[80];
char phasta_iotype[80];

extern int SONFATH;

//extern "C" void proces();
//extern "C" void itrdrv();
//extern "C" void input(int* size,int* myrank);
//extern int input_fform();

int myrank; /* made file global for ease in debugging */

void
catchDebugger() {
    static volatile int debuggerPresent =0;
    while (!debuggerPresent ); // assign debuggerPresent=1
}

// some useful debugging functions

void
pdarray( void* darray , int start, int end ) {
    for( int i=start; i < end; i++ ){
        std::cout << ((double*)darray)[i] << std::endl;
    }
}

void
piarray( void* iarray , int start, int end ) {
    for( int i=start; i < end; i++ ){
        std::cout << ((int*)iarray)[i] << std::endl;
    }
}

extern "C" int 
phasta( int argc,   
        char *argv[] ) {
  
    int size,ierr;
    char inpfilename[100];

    MPI_Init(&argc,&argv);
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
#if (defined WIN32)
    if(argc > 2 ){
      catchDebugger();
    }
#endif
#if ( defined LAUNCH_GDB ) && !( defined WIN32 )

    if ( getenv( "catchDebugger" ) ) {

        int parent_pid = getpid();
        int gdb_child = fork();

        if( gdb_child == 0 ) {
     
            std::cout << "Debugger Process initiating" << std::endl;
            strstream exec_string;
         
#if ( defined decalp )
            exec_string <<"xterm -e idb " 
                        << " -pid "<< parent_pid <<" "<< argv[0] << std::endl;
#endif
#if ( defined LINUX )
            exec_string <<"idb -gui" 
                        << " -pid "<< parent_pid <<" "<< argv[0] << std::endl;
#endif
#if ( defined SUN4 )
            exec_string <<"xterm -e dbx " 
                        << " - "<< parent_pid <<" "<< argv[0] << std::endl;
#endif
#if ( defined IRIX )
            exec_string <<"xterm -e dbx " 
                        << " -p "<< parent_pid <<" "<< argv[0] << std::endl;
#endif
            system( exec_string.str() );
            exit(0);
        }
        catchDebugger();
    }

#endif

    /* Input data  */
    //if(argc > 1 ){
    //    strcpy(inpfilename,argv[1]);
    //} else {
    //    strcpy(inpfilename,"solver.inp");
    //}
    ierr = input_fform();
    if(!ierr){

        /* Preprocess data and run the problem  */
        /* Partition the problem to the correct number of processors */
        if( size > 1 ) {
            if( myrank == 0 ) {
                 Partition_Problem( size, phasta_iotype, 
                                    phasta_iotype, SONFATH );
                 MPI_Barrier(MPI_COMM_WORLD);
            } else { 
                 MPI_Barrier(MPI_COMM_WORLD);
                 sprintf(inpfilename,"%d-procs_case/",size);
                 if( !chdir( inpfilename ) ) {
                    std::cout << " changed directory to " 
                              << inpfilename << std::endl;
                 }else {
                    std::cerr << "could not change to the problem directory "
                              << inpfilename << std::endl;
                    return 1;
                 }
            }
        }

        input(&size,&myrank);
        /* now we can start the solver */
        proces();
		itrdrv();
    }
    else{
        printf("error during reading ascii input \n");
    }
    
    MPI_Finalize();
    return 0;
}
