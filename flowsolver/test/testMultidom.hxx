#ifndef TESTMULTIDOM_HPP_
#define TESTMULTIDOM_HPP_

#include "multidom.hxx"
#include "fileReaders.hxx"
#include "abstractBoundaryCondition.hxx"
#include "fortranPointerManager.hxx"
#include "boundaryConditionManager.hxx"
#include "gtest/gtest.h"
#include "common_c.h"
#include <typeinfo>
#include "debuggingToolsForCpp.hxx"
#include <boost/shared_ptr.hpp>

	// The fixture for testing class Foo.
	class testMultidom : public ::testing::Test {

	protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.

	   rcrtReader* rcrtReader_instance;
	   controlledCoronaryReader* controlledCoronaryReader_instance;
	   std::vector<boost::shared_ptr<abstractBoundaryCondition>>* retrievedBoundaryConditions;
	   boundaryConditionManager* boundaryConditionManager_instance;
	   fortranBoundaryDataPointerManager* fortranPointerManager_instance;
	   netlistReader* netlistReader_instance;

	   double press1;
	   double press2;
	   double press3;
	   double press4;

	   double flow1;
	   double flow2;
	   double flow3;
	   double flow4;

	   // A fake timdat for testing
	   double alfi_local;
	   // Fake delt
	   double delt;

	  testMultidom() {
	  }

	  virtual ~testMultidom() {
	    // You can do clean-up work that doesn't throw exceptions here.
	    // (*fortranPointerManager_instance).~fortranBoundaryDataPointerManager();
	  }

	  // If the constructor and destructor are not enough for setting up
	  // and cleaning up each test, you can define the following methods:

	  virtual void SetUp() {
	    // Code here will be called immediately after the constructor (right
	    // before each test).
	  	overrideMissingDataForTesting();

	  	boundaryConditionManager_instance = boundaryConditionManager::Instance();

	  	boundaryConditionManager_instance->setDelt(0.001);
		boundaryConditionManager_instance->setHstep(5);
		boundaryConditionManager_instance->setAlfi(0.5);
		boundaryConditionManager_instance->setLstep(0);
		boundaryConditionManager_instance->setNtout(1);
		boundaryConditionManager_instance->setMaxsurf(MAXSURF);
		boundaryConditionManager_instance->setNstep(5);
	  	
	  	press1 = 1000;
	    press2 = 2000;
	    press3 = 3000;
	    press4 = 4000;

	    flow1 = 33;
	    flow2 = 66;
	    flow3 = 99;
	    flow4 = 132;

	    // Create fake (i.e. non-FORTRAN) pointer manager
	    // fortranBoundaryDataPointerManager* fortranPointerManager_instance;
	    fortranPointerManager_instance = fortranBoundaryDataPointerManager::Get();
	    
	    fortranPointerManager_instance->boundaryFlows.clear();
	    fortranPointerManager_instance->boundaryPressures.clear();
	    
	    // Insert fake pointer data:
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(3,&flow1));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(7,&flow2));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(9,&flow3));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(11,&flow4));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(12,&flow4));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(13,&flow4));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(14,&flow4));
	    
		fortranPointerManager_instance->hasBoundaryFlows = true;

	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(3,&press1));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(7,&press2));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(9,&press3));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(11,&press4));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(12,&press4));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(13,&press4));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(14,&press4));        

	    fortranPointerManager_instance->hasBoundaryPressures = true;

	  	// Setup the file reader for the RCRTs
		rcrtReader_instance = rcrtReader::Instance();
		rcrtReader_instance->setFileName("rcrt_test.dat");
		rcrtReader_instance->readAndSplitMultiSurfaceInputFile();

		// Setup the file reader for the coronaries:
		controlledCoronaryReader_instance = controlledCoronaryReader::Instance();
		controlledCoronaryReader_instance->setFileName("controlled_coronaries_test.dat");
		controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile();

		// Setup the netlist reader:
		netlistReader_instance = netlistReader::Instance();
	    netlistReader_instance->setFileName("netlist_surfaces.dat");
	    netlistReader_instance->readAndSplitMultiSurfaceInputFile();

		
		std::vector<std::pair<int,boundary_condition_t>> surfaceList;
		surfaceList.clear();
		

	    // Describe 4 test BCs that we want to construct:
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (3,BoundaryCondition_RCR));
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (7,BoundaryCondition_Netlist));
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (9,BoundaryCondition_RCR));
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (11,BoundaryCondition_ControlledCoronary));
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (12,BoundaryCondition_Netlist));
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (13,BoundaryCondition_Netlist));
	    surfaceList.push_back(std::pair <int,boundary_condition_t> (14,BoundaryCondition_Netlist));

	    // get the boundary condition manager
		// boundaryConditionManager_instance->boundaryConditions.clear();
		boundaryConditionManager_instance->setSurfaceList(surfaceList);
		
		retrievedBoundaryConditions = boundaryConditionManager::Instance()->getBoundaryConditions();
		

		(*retrievedBoundaryConditions)[0]->dp_dq = 1.234;
		(*retrievedBoundaryConditions)[1]->dp_dq = 2.234;
		(*retrievedBoundaryConditions)[2]->dp_dq = 3.234;
		(*retrievedBoundaryConditions)[3]->dp_dq = 4.234;
		

		(*retrievedBoundaryConditions)[0]->Hop = 1.567;
		(*retrievedBoundaryConditions)[1]->Hop = 2.567;
		(*retrievedBoundaryConditions)[2]->Hop = 3.567;
		(*retrievedBoundaryConditions)[3]->Hop = 4.567;

		// Force an update of the pressure_n in the boundary conditions
		// (this was introduced as a least-bad hack to deal with order or initialisation,
		// which has a knock-on effect that we must use it here, too...
		// as a bonus, this acts as a test of setPressureFromFortran.)
		boundaryConditionManager_instance->setPressureFromFortran();

		// (*retrievedBoundaryConditions)[0]->pressure_n = *((*retrievedBoundaryConditions)[0]->pressure_n_ptr);
		// (*retrievedBoundaryConditions)[1]->pressure_n = *((*retrievedBoundaryConditions)[1]->pressure_n_ptr);
		// (*retrievedBoundaryConditions)[2]->pressure_n = *((*retrievedBoundaryConditions)[2]->pressure_n_ptr);
		// (*retrievedBoundaryConditions)[3]->pressure_n = *((*retrievedBoundaryConditions)[3]->pressure_n_ptr);
		// (*retrievedBoundaryConditions)[4]->pressure_n = *((*retrievedBoundaryConditions)[4]->pressure_n_ptr);
		// (*retrievedBoundaryConditions)[5]->pressure_n = *((*retrievedBoundaryConditions)[5]->pressure_n_ptr);
		// (*retrievedBoundaryConditions)[6]->pressure_n = *((*retrievedBoundaryConditions)[6]->pressure_n_ptr);

	    // alfi_local = 0.5;
	    // delt = 0.001;
	    
	 //    (*retrievedBoundaryConditions)[0]->delt = delt;
		// (*retrievedBoundaryConditions)[1]->delt = delt;
		// (*retrievedBoundaryConditions)[2]->delt = delt;
		// (*retrievedBoundaryConditions)[3]->delt = delt;
		// (*retrievedBoundaryConditions)[4]->delt = delt;
		// (*retrievedBoundaryConditions)[5]->delt = delt;
		// (*retrievedBoundaryConditions)[6]->delt = delt;
		

		// (*retrievedBoundaryConditions)[0]->alfi_local = alfi_local;
		// (*retrievedBoundaryConditions)[1]->alfi_local = alfi_local;
		// (*retrievedBoundaryConditions)[2]->alfi_local = alfi_local;
		// (*retrievedBoundaryConditions)[3]->alfi_local = alfi_local;
		// (*retrievedBoundaryConditions)[4]->alfi_local = alfi_local;
		// (*retrievedBoundaryConditions)[5]->alfi_local = alfi_local;
		// (*retrievedBoundaryConditions)[6]->alfi_local = alfi_local;
		
	  }

	  virtual void TearDown() {
	    // Code here will be called immediately after each test (right
	    // before the destructor).
	    multidom_finalise();
	    fortranPointerManager_instance->tearDown();
	    // retrievedBoundaryConditions = 0;
	  }

	  void overrideMissingDataForTesting() {
	  	// This function is to be used to create fake data to avoid testing problems
	  	// due to no solver.inp being read prior to testing some of the smaller features.
	  	//
	  	// Please document each override very carefully, and write a note to the console!

	  	// This is because timdat.lstep is used to determine whether we're restarting a simulation
	  	// and we're not doing that during tests.
	  	// std::cout << "Information -- I'm overriding the following variables for this test..." << std::endl;
	  	// std::cout << "Information -- The overrides are in testMultidom.hxx." << std::endl;
	  	// std::cout << " timdat.lstep" << std::endl;
	  	// timdat.lstep = int(0);
	  	// std::cout << " inpdat.Delt[0]" << std::endl;
	  	// inpdat.Delt[0] = 0.01;
	  	// grcrbccom.numGRCRSrfs;
	  	// std::cout << "nomodule.numNetlistLPNSrfs" << std::endl;
    //     nomodule.numNetlistLPNSrfs = int(1);
        // nomodule.numControlledCoronarySrfs;
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibraryTestMultidom();

#endif
	
