#include "gtest/gtest.h"
#include "testOrphans.hxx"
#include "fileReaders.hxx"
#include "fortranPointerManager.hxx"
#include "NetlistBoundaryCondition.hxx"
#include <stdexcept>

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestOrphans() { return 0; }

// The component at the 3D interface must be oriented so that the 3D domain is connected
// to its startNode (not its endNode). Which is which is determined by the order of the
// two node indices given in netlist_surfaces.dat forthe component.
TEST_F(testOrphans, checkNetlistDetectsBad3DInterfaceComponentOrientation) {
  // Create fake (i.e. non-FORTRAN) pointer manager
  fortranPointerManager_instance = fortranBoundaryDataPointerManager::Get();

  fortranPointerManager_instance->boundaryFlows.clear();
  fortranPointerManager_instance->boundaryPressures.clear();

  double fakeFlow = 0.0;
  double fakePressure = 0.0;
  // Insert fake pointer data:
  fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(2,&fakeFlow));
  fortranPointerManager_instance->hasBoundaryFlows = true;

  fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(2,&fakePressure));
  fortranPointerManager_instance->hasBoundaryPressures = true;

  // Setup the netlist reader:
  netlistReader_instance = NetlistReader::Instance();
  netlistReader_instance->setFileName("netlist_surfaces_bad3DInterfaceComponentOrientation.dat");
  netlistReader_instance->readAndSplitMultiSurfaceInputFile();

  std::vector<std::pair<int,boundary_condition_t>> surfaceList;
  surfaceList.clear();

  // Describe test BC that we want to construct:
  surfaceList.push_back(std::pair <int,boundary_condition_t> (2,BoundaryCondition_Netlist));

  // Testing for "EE: The netlist component at the 3D interface must be connected to it via its start node, not its end node."
  EXPECT_THROW(boundaryConditionManager_instance->setSurfaceList(surfaceList), std::runtime_error);
  // boundaryConditionManager_instance->setSurfaceList(surfaceList);
}

// The component at the 3D interface must have a node which connects to no other components (i.e. just to the 3D domain)
TEST_F(testOrphans, checkNetlistDetectsBadComponentAt3DInterface) {
  // Create fake (i.e. non-FORTRAN) pointer manager
  fortranPointerManager_instance = fortranBoundaryDataPointerManager::Get();

  fortranPointerManager_instance->boundaryFlows.clear();
  fortranPointerManager_instance->boundaryPressures.clear();

  double fakeFlow = 0.0;
  double fakePressure = 0.0;
  // Insert fake pointer data:
  fortranPointerManager_instance->boundaryFlows.insert(std::pair<int,double*>(2,&fakeFlow));
  fortranPointerManager_instance->hasBoundaryFlows = true;

  fortranPointerManager_instance->boundaryPressures.insert(std::pair<int,double*>(2,&fakePressure));
  fortranPointerManager_instance->hasBoundaryPressures = true;

  // Setup the netlist reader:
  netlistReader_instance = NetlistReader::Instance();
  netlistReader_instance->setFileName("netlist_surfaces_badComponentAt3DInterface.dat");
  netlistReader_instance->readAndSplitMultiSurfaceInputFile();

  std::vector<std::pair<int,boundary_condition_t>> surfaceList;
  surfaceList.clear();

  // Describe test BC that we want to construct:
  surfaceList.push_back(std::pair <int,boundary_condition_t> (2,BoundaryCondition_Netlist));

  // get the boundary condition manager
  boundaryConditionManager_instance = boundaryConditionManager::Instance();
  // Testing for "EE: Only one component may be directly connected to the 3D interface in the Netlist."
  EXPECT_THROW(boundaryConditionManager_instance->setSurfaceList(surfaceList), std::runtime_error);
  // boundaryConditionManager_instance->setSurfaceList(surfaceList);
}