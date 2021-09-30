#ifndef TESTMULTIDOM_HPP_
#define TESTMULTIDOM_HPP_

#include "multidom.hxx"
#include "fileReaders.hxx"
#include "AbstractBoundaryCondition.hxx"
#include "FortranBoundaryDataPointerManager.hxx"
#include "BoundaryConditionManager.hxx"
#include "gtest/gtest.h"
#include "common_c.h"
#include <typeinfo>
#include "debuggingToolsForCpp.hxx"
#include <boost/shared_ptr.hpp>
#include "CrimsonGlobalArrayTransfer.h"

	// The fixture for testing class Foo.
	class testMultidom : public ::testing::Test {

	protected:
	  // You can remove any or all of the following functions if its body
	  // is empty.

	   RcrtReader* rcrtReader_instance;
	   ControlledCoronaryReader* controlledCoronaryReader_instance;
	   std::vector<boost::shared_ptr<AbstractBoundaryCondition>>* retrievedBoundaryConditions;
	   BoundaryConditionManager* boundaryConditionManager_instance;
	   FortranBoundaryDataPointerManager* fortranPointerManager_instance;
	   NetlistReader* netlistReader_instance;

	   double press1;
	   double press2;
	   double press3;
	   double press4;

	   double flow1;
	   double flow2;
	   double flow3;
	   double flow4;

	   // A fake timdat for testing
	   double m_generalizedAlphaMethodAlpha;
	   // Fake delt
	   double delt;

	  testMultidom() {
	  }

	  virtual ~testMultidom() {
	    // You can do clean-up work that doesn't throw exceptions here.
	    // (*fortranPointerManager_instance).~FortranBoundaryDataPointerManager();
	  }

	  // If the constructor and destructor are not enough for setting up
	  // and cleaning up each test, you can define the following methods:

	  virtual void SetUp() {
	    // Code here will be called immediately after the constructor (right
	    // before each test).
	  	overrideMissingDataForTesting();

	  	boundaryConditionManager_instance = BoundaryConditionManager::Instance();

	  	boundaryConditionManager_instance->setDelt(0.001);
		boundaryConditionManager_instance->setHstep(5);
		boundaryConditionManager_instance->setAlfi(0.5);
		// boundaryConditionManager_instance->setStartingTimestepIndex(0);
		boundaryConditionManager_instance->setNtout(1);
		boundaryConditionManager_instance->setMaxsurf(MAXSURF);
		boundaryConditionManager_instance->setNstep(5);
		boundaryConditionManager_instance->setNumLoopClosingnetlistCircuits(0);
		boundaryConditionManager_instance->setMasterControlScriptPresent(0);
	  	
	  	press1 = 1000;
	    press2 = 2000;
	    press3 = 3000;
	    press4 = 4000;

	    flow1 = 33;
	    flow2 = 66;
	    flow3 = 99;
	    flow4 = 132;

	    // Create fake (i.e. non-FORTRAN) pointer manager
	    // FortranBoundaryDataPointerManager* fortranPointerManager_instance;
	    fortranPointerManager_instance = FortranBoundaryDataPointerManager::Get();
	    
	    fortranPointerManager_instance->m_boundaryFlows.clear();
	    fortranPointerManager_instance->m_boundaryPressures.clear();
	    
	    // Insert fake pointer data:
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(3,&flow1));
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(7,&flow2));
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(9,&flow3));
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(11,&flow4));
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(12,&flow4));
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(13,&flow4));
	    fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(14,&flow4));
	    
		fortranPointerManager_instance->m_hasBoundaryFlows = true;

	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(3,&press1));
	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(7,&press2));
	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(9,&press3));
	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(11,&press4));
	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(12,&press4));
	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(13,&press4));
	    fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(14,&press4));        

	    fortranPointerManager_instance->m_hasBoundaryPressures = true;

	  	// Setup the file reader for the RCRTs
		rcrtReader_instance = RcrtReader::Instance();
		rcrtReader_instance->setFileName("rcrt_test.dat");
		rcrtReader_instance->readAndSplitMultiSurfaceInputFile();

		// Setup the file reader for the coronaries:
		controlledCoronaryReader_instance = ControlledCoronaryReader::Instance();
		controlledCoronaryReader_instance->setFileName("controlled_coronaries_test.dat");
		controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile();

		// Setup the netlist reader:
		// boost::filesystem::current_path(boost::filesystem::path("basicTestFiles"));
		netlistReader_instance = NetlistReader::Instance();
		if (boost::filesystem::exists(boost::filesystem::path("netlist_surfaces.dat")))
  		{
		    netlistReader_instance->setFileName("netlist_surfaces.dat");
		    netlistReader_instance->readAndSplitMultiSurfaceInputFile();
		    // for converting old netlist specification file format to new (generally not important for actual simulations)
    		netlistReader_instance->writeCircuitSpecificationInXmlFormat();
		}
		// boost::filesystem::current_path(boost::filesystem::path(".."));

		
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
		CrimsonGlobalArrayTransfer::Get()->initialiseForRCRFiltering(2);
		boundaryConditionManager_instance->setSimulationModePurelyZeroD(0);
		boundaryConditionManager_instance->setSurfaceList(surfaceList);
		
		retrievedBoundaryConditions = BoundaryConditionManager::Instance()->getBoundaryConditions();
		

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
		
	  }

	  virtual void TearDown() {
	    // Code here will be called immediately after each test (right
	    // before the destructor).
	    multidom_finalise();
	    fortranPointerManager_instance->tearDown();
	    boundaryConditionManager_instance->tearDown();
	    CrimsonGlobalArrayTransfer::Get()->tearDown();
	    // retrievedBoundaryConditions = 0;

	    // we do this because the new NetlistXmlReader does not allow the use of 
	  	// non-standard file names, and the tests here work with e.g. 
	  	// 
	  	// netlist_surfaces_bad3DInterfaceComponentOrientation.dat
	  	//
	  	// which becomes just
	  	//
	  	// netlist_surfaces.xml
	  	//
	  	// This wil clash with subsequent tests, so we move it safely to
	  	// a storage directory
	  	MPI_Barrier(MPI_COMM_WORLD);
	  	int rank;
	  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  	if (rank==0)
	  	{
		  	boost::filesystem::path oldXmlFilesDirectory("oldXmlFiles");
		  	if (!boost::filesystem::exists(oldXmlFilesDirectory))
		  	{
				boost::filesystem::create_directory(oldXmlFilesDirectory);
			}

		  	boost::filesystem::path basePath("oldXmlFiles/netlist_surfaces.xml");
		  	boost::filesystem::path pathToTry = basePath;
		  	int fileSuffix = 0;
		  	while (boost::filesystem::exists(pathToTry))
		  	{
		  		pathToTry = basePath;
		  		pathToTry += boost::filesystem::path("."+boost::lexical_cast<std::string>(fileSuffix));
		  		fileSuffix++;
		  	}
		    boost::filesystem::rename(boost::filesystem::path("netlist_surfaces.xml"), pathToTry);
		    std::cout << "Moving the just-generated netlist_surfaces.xml to ./oldXmlFiles/" << std::endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	  }

	  void overrideMissingDataForTesting() {
	  	// This function is to be used to create fake data to avoid testing problems
	  	// due to no solver.inp being read prior to testing some of the smaller features.
	  	//
	  	// Please document each override very carefully, and write a note to the console!

	  	// This is because timdat.currentTimestepIndex is used to determine whether we're restarting a simulation
	  	// and we're not doing that during tests.
	  	// std::cout << "Information -- I'm overriding the following variables for this test..." << std::endl;
	  	// std::cout << "Information -- The overrides are in testMultidom.hxx." << std::endl;
	  	// std::cout << " timdat.currentTimestepIndex" << std::endl;
	  	// timdat.currentTimestepIndex = int(0);
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
	
