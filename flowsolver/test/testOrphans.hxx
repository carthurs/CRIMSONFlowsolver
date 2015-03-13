#ifndef TESTORPHANS_HPP_
#define TESTORPHANS_HPP_

#include "multidom.hxx"
#include "fileReaders.hxx"
#include "boundaryConditionManager.hxx"
#include "gtest/gtest.h"
#include "common_c.h"
#include <typeinfo>
#include "debuggingToolsForCpp.hxx"
#include <boost/shared_ptr.hpp>

	// The fixture for testing class Foo.
	class testOrphans : public ::testing::Test {

	protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.
	  netlistReader* netlistReader_instance;
	  boundaryConditionManager* boundaryConditionManager_instance;
	  fortranBoundaryDataPointerManager* fortranPointerManager_instance;

	  testOrphans() {
	  	// get the boundary condition manager
	  	timdat.lstep = 0; // this is needed because the boundaryConditionManager constructor looks for it explicitly!
  		boundaryConditionManager_instance = boundaryConditionManager::Instance();
  		overrideMissingDataForTesting();
	  }

	  virtual ~testOrphans() {
	    // You can do clean-up work that doesn't throw exceptions here.
	    // (*fortranPointerManager_instance).~fortranBoundaryDataPointerManager();
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
	    // multidom_finalise();
	    // fortranPointerManager_instance->tearDown();
	    // retrievedBoundaryConditions = 0;
	    netlistReader_instance->Term();
	    boundaryConditionManager_instance->Term();
	    fortranPointerManager_instance->tearDown();
	  }

	  void overrideMissingDataForTesting() {
	  	// This function is to be used to create fake data to avoid testing problems
	  	// due to no solver.inp being read prior to testing some of the smaller features.
	  	//
	  	// Please document each override very carefully, and write a note to the console!

	  	// This is because timdat.lstep is used to determine whether we're restarting a simulation
	  	// and we're not doing that during tests.
	  	std::cout << "Information -- I'm overriding the following variables for this test..." << std::endl;
	  	std::cout << "Information -- The overrides are in testOrphans.hxx." << std::endl;
	  	std::cout << "Information -- delt, hstep, alfi, lstep, ntout." << std::endl;
	  	boundaryConditionManager_instance->setDelt(0.01);
		boundaryConditionManager_instance->setHstep(5);
		boundaryConditionManager_instance->setAlfi(0.5);
		boundaryConditionManager_instance->setLstep(0);
		boundaryConditionManager_instance->setNtout(1);
		boundaryConditionManager_instance->setMaxsurf(199);
		boundaryConditionManager_instance->setNstep(5);
	  	// grcrbccom.numGRCRSrfs;
	  	// std::cout << "nomodule.numNetlistLPNSrfs" << std::endl;
    //     nomodule.numNetlistLPNSrfs = int(1);
        // nomodule.numControlledCoronarySrfs;
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibraryTestOrphans();

#endif
	
