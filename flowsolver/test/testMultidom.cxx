#include "testMultidom.hxx"
#include "gtest/gtest.h"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMultidom() { return 0; }

// Tests that the fileReader reads the first entry of the test-case rcrt.dat
TEST_F(testMultidom, boundaryConditionsMadeProperly) {
  // Check we got the right boundary conditions
	// EXPECT_TRUE(1==1);
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[0])==typeid(RCR));
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[1])==typeid(netlist));
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[2])==typeid(RCR));

  EXPECT_EQ((*retrievedBoundaryConditions)[0]->surfaceIndex,3);
  EXPECT_EQ((*retrievedBoundaryConditions)[1]->surfaceIndex,7);
  EXPECT_EQ((*retrievedBoundaryConditions)[2]->surfaceIndex,2342);
}
