#ifndef TESTMAINWITHZERODDOMAIN_HXX_
#define TESTMAINWITHZERODDOMAIN_HXX_

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
#include "pureZeroDDriver.hxx"
#include "itrdrv.h"
#include "itrPC.h"
#include "partition.h"
#include "input_fform.h"
#include "multidom.hxx"
#include "FortranBoundaryDataPointerManager.hxx"
#include "fileReaders.hxx"
#include "BoundaryConditionManager.hxx"
#include "CrimsonGlobalArrayTransfer.h"

#include "debuggingToolsForCpp.hxx"
#include <boost/filesystem/path.hpp>

#ifdef intel
#include <direct.h>
#else
#include <unistd.h>
#endif



// The fixture for testing class Foo.
	class testMainWithZeroDDomain : public ::testing::Test {
	 public:
	 protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.
	 	boost::filesystem::path dirBinaryCalledFrom;

	 	// This test environment should setup exactly as main.cxx does in estimation/src/main.cxx
		testMainWithZeroDDomain() {
			MPI_Barrier(MPI_COMM_WORLD);
			dirBinaryCalledFrom = boost::filesystem::current_path();
			getRank();
		}

		void setSimDirectory(std::string dir)
		{
			int success = chdir(dir.c_str());
			if (success==0)
			{
				std::cout << "II: Changed to dir: " << dir << " rank was: " << m_rank << std::endl;
			}
			else
			{
				std::cout << "EE: Failed to change to dir: " << dir << " rank was: " << m_rank << " success was: " << success << std::endl;
				perror("EEEE: Failed to change dir in test");
			}
		}

		void clearOutOldFiles()
		{
			MPI_Barrier(m_iNewComm_C);
			if (m_rank == 0)
			{
				// Warning - this has the potential to delete multiple folders due to the *
				system("rm -rf *-procs-case");
			}
			MPI_Barrier(m_iNewComm_C);
		}

		void runSimulation()
		{
		   char pathToProcsCaseDir[100];

		   // The zero-D domain should just be run single-threaded
		   // if (m_rank==0)
		   {
				// read configuration file
				int errFlag = input_fform();
				if (errFlag != 0)
				{
					throw std::runtime_error("EE: Failed during parsing of input files.");
				}

				// Preprocess data and run the problem
				// Partition the problem to the correct number of processors

				if( m_rank == 0 )
				{
					//cout << "number of procs " << numprocs_perparticle << endl;
					try {
						Partition_Problem( m_numProcsTotal );
					} catch (const std::exception& e) {
					    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
					    throw;
					}
				}

				MPI_Barrier(m_iNewComm_C);

				sprintf(pathToProcsCaseDir,"%d-procs-case",m_numProcsTotal);
				boost::filesystem::path thisDir = boost::filesystem::current_path();
				int errStat = chdir(pathToProcsCaseDir);

				if (errStat != 0)
				{
					std::cerr << "Failed to change to directory " << pathToProcsCaseDir << ". Rank is: " << m_rank << std::endl;
					perror("EE: Failed to change to test directory.");
					throw std::runtime_error("EE: Failed to change to test directory.");
			  	}
				else
				{
				   std::cout << "changing directory to " << pathToProcsCaseDir << std::endl;
				}

				// input(&m_numProcsTotal, &m_rank);
				// proces();

				try {
					// just initialise the time values that the AbstractBoundaryCondition needs (when it's called in multidom_initialise).
					// This will be called again during itrdrv_init.
					// This is really ugly, but a proper fix will take days - it's a BIG refactor.
					int dummyInitialItseqValue=1;
					callFortranSetupTimeParameters(dummyInitialItseqValue);
				} catch (const std::exception& e) {
				    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
				    throw;
				}

				// initialise reduced order boundary conditions
				try {
					multidom_initialise();
				} catch (const std::exception& e) {
				    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
				    throw;
				}

				std::vector<double> zeroDDomainCompliancesForEachConnectedComponent;
                double* const zeroDDomainCompliancesEndPointer = nomodule.zeroDDomainCompliances + nomodule.num3DConnectedComponents + 1; // Fortran indexing
			    zeroDDomainCompliancesForEachConnectedComponent.assign(nomodule.zeroDDomainCompliances + 1, zeroDDomainCompliancesEndPointer);

			    std::cout << std::endl;

				PureZeroDDriver pureZeroDDriver(nomodule.numDirCalcSrfs, zeroDDomainCompliancesForEachConnectedComponent);
				try {
					pureZeroDDriver.setDelt(inpdat.Delt[0]);
				    pureZeroDDriver.setAlfi(timdat.alfi);
				    pureZeroDDriver.setHstep(inpdat.nstep[0] + timdat.currentTimestepIndex);
				    pureZeroDDriver.setNtout(1);
	
				    pureZeroDDriver.setupConnectedComponents(nomodule.num3DConnectedComponents, nomodule.surfacesOfEachConnectedComponent, nomodule.indicesOfNetlistSurfaces, nomodule.nsrflistDirCalc);
	
					pureZeroDDriver.init();
	
					multidomSetupControlSystems();
				} catch (const std::exception& e) {
				    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
				    throw;
				}

				for (int kk = 1; kk <= inpdat.nstep[0]; kk++)
				{
					try {
						pureZeroDDriver.iter_init();
					} catch (const std::exception& e) {
					    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
					    throw;
					}
					try {
						pureZeroDDriver.iter_step();
					} catch (const std::exception& e) {
					    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
					    throw;
					}
					try {
						pureZeroDDriver.iter_finalize();
					} catch (const std::exception& e) {
					    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
					    throw;
					}
					try {
						multidom_iter_finalise();
					} catch (const std::exception& e) {
					    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
					    throw;
					}
				}

				pureZeroDDriver.finalize();
			   
				multidom_finalise();
				deallocate_arrays();
		           
		    }
		    MPI_Barrier(m_iNewComm_C);
		}

	  virtual ~testMainWithZeroDDomain() {
	    // You can do clean-up work that doesn't throw exceptions here.
	    boost::filesystem::current_path(dirBinaryCalledFrom);
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
	    FortranBoundaryDataPointerManager::Get()->tearDown();
	    CrimsonGlobalArrayTransfer::Get()->tearDown();
	    // BoundaryConditionManager::Instance()->Term();
	  }
	private:
		int m_rank;
	    int m_numProcsTotal;
        MPI_Comm m_iNewComm_C;

		void getRank()
		{
		   // save the communicator
		   m_iNewComm_C = MPI_COMM_WORLD;
		   newcom.iNewComm = MPI_Comm_c2f(m_iNewComm_C); // modifies newcom in fortran common block
	           MPI_Barrier(m_iNewComm_C);

		   MPI_Comm_size(m_iNewComm_C, &m_numProcsTotal);
		   MPI_Comm_rank(m_iNewComm_C, &m_rank);
		}

	  // Objects declared here can be used by all tests in the test case for Foo.
	};


int PullInMyLibraryTestMainWithZeroDDomain();
#endif