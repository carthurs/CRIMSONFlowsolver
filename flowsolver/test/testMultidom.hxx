#ifndef TESTMULTIDOM_HPP_
#define TESTMULTIDOM_HPP_

#include "multidom.cxx"
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
	   std::vector<boost::shared_ptr<boundaryCondition>>* retrievedBoundaryConditions;
	   boundaryConditionManager* boundaryConditionManager_instance;
	   fortranBoundaryDataPointerManager* fortranPointerManager_instance;

	   double press1;
	   double press2;
	   double press3;

	   double flow1;
	   double flow2;
	   double flow3;

	   // A fake timdat for testing
	   double alfi_local;
	   // Fake delt
	   double delt;

	  testMultidom() {
	  	boundaryConditionManager_instance = boundaryConditionManager::Instance();
	  	
	  	press1 = 1000;
	    press2 = 2000;
	    press3 = 3000;

	    flow1 = 33;
	    flow2 = 66;
	    flow3 = 99;

	    // Create fake (i.e. non-FORTRAN) pointer manager
	    // fortranBoundaryDataPointerManager* fortranPointerManager_instance;
	    fortranPointerManager_instance = fortranBoundaryDataPointerManager::Get();
	    
	    fortranPointerManager_instance->boundaryFlows.clear();
	    fortranPointerManager_instance->boundaryPressures.clear();
	    
	    // Insert fake pointer data:
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(3,&flow1));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(7,&flow2));
	    fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(9,&flow3));

	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(3,&press1));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(7,&press2));
	    fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(9,&press3));
        

	  	// Setup the file reader for the RCRTs
		rcrtReader_instance = rcrtReader::Instance();
		rcrtReader_instance->setFileName("rcrt_test.dat");
		rcrtReader_instance->readAndSplitMultiSurfaceInputFile();

		
		std::vector<std::pair<int,std::string>> surfaceList;
		surfaceList.clear();
		

	    // Describe 3 test BCs that we want to construct:
	    surfaceList.push_back(std::pair <int,std::string> (3,"rcr"));
	    surfaceList.push_back(std::pair <int,std::string> (7,"netlist"));
	    surfaceList.push_back(std::pair <int,std::string> (9,"rcr"));

	    

	    // get the boundary condition manager
		// boundaryConditionManager_instance->boundaryConditions.clear();
		boundaryConditionManager_instance->setSurfaceList(surfaceList);
		
		retrievedBoundaryConditions = boundaryConditionManager::Instance()->getBoundaryConditions();
		

		(*retrievedBoundaryConditions)[0]->dp_dq = 1.234;
		(*retrievedBoundaryConditions)[1]->dp_dq = 2.234;
		(*retrievedBoundaryConditions)[2]->dp_dq = 3.234;
		

		(*retrievedBoundaryConditions)[0]->Hop = 1.567;
		(*retrievedBoundaryConditions)[1]->Hop = 2.567;
		(*retrievedBoundaryConditions)[2]->Hop = 3.567;
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
	    alfi_local = 0.5;
	    delt = 0.001;
	    (*retrievedBoundaryConditions)[0]->delt = delt;
		(*retrievedBoundaryConditions)[1]->delt = delt;
		(*retrievedBoundaryConditions)[2]->delt = delt;
		

		(*retrievedBoundaryConditions)[0]->alfi_local = alfi_local;
		(*retrievedBoundaryConditions)[1]->alfi_local = alfi_local;
		(*retrievedBoundaryConditions)[2]->alfi_local = alfi_local;
		
	  }

	  virtual void TearDown() {
	    // Code here will be called immediately after each test (right
	    // before the destructor).
	  	std::cout << "here1" <<std::endl;
	    
	    std::cout << "here2" <<std::endl;

	    boundaryConditionManager_instance->Term();
	    fortranPointerManager_instance->tearDown();
	    rcrtReader_instance->tearDown();
	    retrievedBoundaryConditions = 0;
	    
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	};

// Hack to force the compiler to link this test to the relevant main() for testing
	int PullInMyLibraryTestMultidom();

#endif
	
