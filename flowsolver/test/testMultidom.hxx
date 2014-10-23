#ifndef TESTMULTIDOM_HPP_
#define TESTMULTIDOM_HPP_

#include "multidom.cxx"
#include "gtest/gtest.h"
#include "common_c.h"
#include <typeinfo>

	// The fixture for testing class Foo.
	class testMultidom : public ::testing::Test {
	 protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.

	   rcrtReader* rcrtReader_instance;
	   std::vector<std::unique_ptr<boundaryCondition>>* retrievedBoundaryConditions;
	   fortranBoundaryDataPointerManager* fortranPointerManager_instance;
	   double press1;
	   double press2;
	   double press3;

	   double flow1;
	   double flow2;
	   double flow3;

	  testMultidom() {

	  	press1 = 1000;
	    press2 = 2000;
	    press3 = 3000;

	    flow1 = 33;
	    flow2 = 66;
	    flow3 = 99;

	    // Create fake (i.e. non-FORTRAN) pointer manager
	    fortranPointerManager_instance = fortranBoundaryDataPointerManager::Get();
	    
	    // Insert fake pointer data:
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(3,&press1));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(7,&press2));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(2341,&press3));

	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(3,&flow1));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(7,&flow2));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(2341,&flow3));



	  	// Setup the file reader for the RCRTs
		rcrtReader_instance = rcrtReader::Instance();
		rcrtReader_instance->setFileName("rcrt_test.dat");
		rcrtReader_instance->readAndSplitMultiSurfaceInputFile();

		std::vector<std::pair<int,std::string>> surfaceList;

	    // Describe 3 test BCs that we want to construct:
	    surfaceList.push_back(std::pair <int,std::string> (3,"rcr"));
	    surfaceList.push_back(std::pair <int,std::string> (7,"netlist"));
	    surfaceList.push_back(std::pair <int,std::string> (2341,"rcr"));

	    // get the boundary condition manager
		boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
		boundaryConditionManager_instance->setSurfaceList(surfaceList);
		retrievedBoundaryConditions = boundaryConditionManager_instance->getBoundaryConditions();
	  }

	  virtual ~testMultidom() {
	    // You can do clean-up work that doesn't throw exceptions here.
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
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibraryTestMultidom();

#endif
	
