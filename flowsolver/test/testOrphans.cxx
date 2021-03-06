#include "gtest/gtest.h"
#include "testOrphans.hxx"
#include "fileReaders.hxx"
#include "FortranBoundaryDataPointerManager.hxx"
#include "NetlistBoundaryCondition.hxx"
#include <stdexcept>

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestOrphans() { return 0; }

// The component at the 3D interface must have a node which connects to no other components (i.e. just to the 3D domain)
TEST_F(testOrphans, checkNetlistDetectsBadComponentAt3DInterface) {
  // Create fake (i.e. non-FORTRAN) pointer manager
  fortranPointerManager_instance = FortranBoundaryDataPointerManager::Get();

  fortranPointerManager_instance->m_boundaryFlows.clear();
  fortranPointerManager_instance->m_boundaryPressures.clear();

  double fakeFlow = 0.0;
  double fakePressure = 0.0;
  // Insert fake pointer data:
  fortranPointerManager_instance->m_boundaryFlows.insert(std::pair<int,double*>(2,&fakeFlow));
  fortranPointerManager_instance->m_hasBoundaryFlows = true;

  fortranPointerManager_instance->m_boundaryPressures.insert(std::pair<int,double*>(2,&fakePressure));
  fortranPointerManager_instance->m_hasBoundaryPressures = true;

  // Setup the netlist reader:
  netlistReader_instance = NetlistReader::Instance();
  // if (boost::filesystem::exists(boost::filesystem::path("netlist_surfaces_badComponentAt3DInterface.dat")))
  // {
    netlistReader_instance->setFileName("netlist_surfaces_badComponentAt3DInterface.dat");
    netlistReader_instance->readAndSplitMultiSurfaceInputFile();
  // }
  // for converting old netlist specification file format to new (generally not important for actual simulations)
  netlistReader_instance->writeCircuitSpecificationInXmlFormat();

  std::vector<std::pair<int,boundary_condition_t>> surfaceList;
  surfaceList.clear();

  // Describe test BC that we want to construct:
  surfaceList.push_back(std::pair <int,boundary_condition_t> (2,BoundaryCondition_Netlist));

  // get the boundary condition manager
  boundaryConditionManager_instance = BoundaryConditionManager::Instance();
  // Testing for "EE: Only one component may be directly connected to the 3D interface in the Netlist."
  EXPECT_THROW(boundaryConditionManager_instance->setSurfaceList(surfaceList), std::runtime_error);
  // boundaryConditionManager_instance->setSurfaceList(surfaceList);
}