// Copyright 2006, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "CRIMSONPython.hxx"
#include <cstdlib>
#include <stdio.h>
#include "mpi.h"
#include "petscsys.h"
#include "gtest/gtest.h"
#include "boost/filesystem.hpp"


// Hack to force the compiler to link this main to the test libraries
int PullInMyLibrary();
static int dummy = PullInMyLibrary();
int PullInMyLibraryTestMultidom();
static int dummy2 = PullInMyLibraryTestMultidom();
int PullInMyLibraryTestMain();
static int dummy3 = PullInMyLibraryTestMain();
int PullInMyLibraryTestOrphans();
static int dummy4 = PullInMyLibraryTestOrphans();
int PullInMyLibraryTestMainWithZeroDDomain();
static int dummy5 =PullInMyLibraryTestMainWithZeroDDomain();

static char help[] = "Google test, modified to work with CRIMSON.\n\n";

GTEST_API_ int main(int argc, char **argv) {

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

  // Fake command line input for the flowsolver, so we don't have to pass it:
  int fake_argc = 2;
  char *fake_argv_temp[4] = {"testMain","1","solver.inp",NULL};
  char** fake_argv= (char**) fake_argv_temp;

  // MPI_Init(&fake_argc,(char***)&fake_argv);
  PetscInitialize(&fake_argc, &fake_argv, (char *)0, help);

  initialisePython();

  int testSuccess = 0;

  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);

  // We don't actually use testSuccess, but I'll leave it here in case
  // we find a use for it later .
  testSuccess = RUN_ALL_TESTS();

  Py_Finalize();
  int ierr;
  ierr = PetscFinalize();
  CHKERRQ(ierr);
  return 0;
}
