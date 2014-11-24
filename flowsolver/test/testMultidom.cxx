#include "gtest/gtest.h"
#include "testMultidom.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMultidom() { return 0; }

// Tests that the correct types of boundary conditions were made, with the correct given surface indices.
TEST_F(testMultidom, checkBoundaryConditionsMadeProperly) {
  // Check we got the right boundary conditions
	// EXPECT_TRUE(1==1);
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[0])==typeid(RCR));
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[1])==typeid(netlist));
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[2])==typeid(RCR));

  EXPECT_EQ((*retrievedBoundaryConditions)[0]->surfaceIndex,3);
  EXPECT_EQ((*retrievedBoundaryConditions)[1]->surfaceIndex,7);
  EXPECT_EQ((*retrievedBoundaryConditions)[2]->surfaceIndex,9);

  EXPECT_EQ(*((*retrievedBoundaryConditions)[0]->flow_n_ptr),flow1);
  // EXPECT_EQ(*((*retrievedBoundaryConditions)[1]->flow_n_ptr),flow2); // uncomment once Netlist initialisation sets this pointer
  EXPECT_EQ(*((*retrievedBoundaryConditions)[2]->flow_n_ptr),flow3);

  EXPECT_EQ(*((*retrievedBoundaryConditions)[0]->pressure_n_ptr),press1);
  // EXPECT_EQ(*((*retrievedBoundaryConditions)[1]->pressure_n_ptr),press2); // uncomment once Netlist initialisation sets this pointer
  EXPECT_EQ(*((*retrievedBoundaryConditions)[2]->pressure_n_ptr),press3);
  
}

TEST_F(testMultidom, checkRCRLinearInterpolators)
{
	// Test using the values read from rcrt_test.dat	
	double testValue;

	testValue = (*retrievedBoundaryConditions)[2]->linInterpolateTimeData(-1.0,3);
	EXPECT_DOUBLE_EQ(testValue,0.0);
	testValue = ((*retrievedBoundaryConditions)[2])->linInterpolateTimeData(0.25,3);
	EXPECT_DOUBLE_EQ(testValue,0.85);
	testValue = (*retrievedBoundaryConditions)[2]->linInterpolateTimeData(0.5,3);
	EXPECT_DOUBLE_EQ(testValue,1.7);
	testValue = (*retrievedBoundaryConditions)[2]->linInterpolateTimeData(0.6,3);
	EXPECT_DOUBLE_EQ(testValue,1.85);
	testValue = (*retrievedBoundaryConditions)[2]->linInterpolateTimeData(10.0,3);
	EXPECT_DOUBLE_EQ(testValue,2.3);

	// Check the netlists throw an exception when we try to call their 
	// non-existent linInterpolateTimeData (should hit the virtual method 
	// in the base class, and throw.)
	EXPECT_ANY_THROW((*retrievedBoundaryConditions)[1]->linInterpolateTimeData(10.0,3));
}

TEST_F(testMultidom, checkDpDqAndHopFortranPasser)
{
	// The braces initialise the array to zero
	double implicitCoeffs_toBeFilled[400] = {};
	

	boundaryConditionManager_instance->getImplicitCoeff_rcr(implicitCoeffs_toBeFilled);
	

	EXPECT_DOUBLE_EQ(implicitCoeffs_toBeFilled[0],1.234);
	EXPECT_DOUBLE_EQ(implicitCoeffs_toBeFilled[1],3.234);
	EXPECT_DOUBLE_EQ(implicitCoeffs_toBeFilled[2],0.0);
	EXPECT_DOUBLE_EQ(implicitCoeffs_toBeFilled[200],1.567);
	EXPECT_DOUBLE_EQ(implicitCoeffs_toBeFilled[201],3.567);
	EXPECT_DOUBLE_EQ(implicitCoeffs_toBeFilled[202],0.0);
}

TEST_F(testMultidom, checkImplicitConditionComputation_solve)
{
	(*retrievedBoundaryConditions)[0]->computeImplicitCoeff_solve(int(1));
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[0]->getdp_dq(),0.38423669121132564);
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[0]->getHop(),9.859056500341553e+02);

	(*retrievedBoundaryConditions)[2]->computeImplicitCoeff_solve(int(1));
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[2]->getdp_dq(),44.123009016155727);
	// Just trying out a different sort of EXPECT, this time with a tolerance option:
	EXPECT_NEAR((*retrievedBoundaryConditions)[2]->getHop(),-1.368175114687044e+03,1e-12);
}

TEST_F(testMultidom, checkImplicitConditionComputation_update)
{
	(*retrievedBoundaryConditions)[2]->computeImplicitCoeff_update(int(1));
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[2]->dp_dq_n1,44.123018032309020);
	// Just trying out a different sort of EXPECT, this time with a tolerance option:
	EXPECT_NEAR((*retrievedBoundaryConditions)[2]->Hop_n1,-1.368173229374138e+03,1e-12);
}

TEST_F(testMultidom, checkFlowAndPressureSetters)
{
	
	(*retrievedBoundaryConditions)[0]->dp_dq_n1 = 534.1;
	(*retrievedBoundaryConditions)[0]->Hop_n1 = 0.001;
	double fakeFlow1 = 10.0;
	(*retrievedBoundaryConditions)[0]->flow_n_ptr = &fakeFlow1;
	
	(*retrievedBoundaryConditions)[2]->dp_dq_n1 = 534.1;
	(*retrievedBoundaryConditions)[2]->Hop_n1 = 0.001;
	double fakeFlow2 = 15.0;
	(*retrievedBoundaryConditions)[2]->flow_n_ptr = &fakeFlow2;
	
	boundaryConditionManager_instance->updateAllRCRS_Pressure_n1_withflow();
	

	// Check the setting and computation worked on the RCRs (these have indices 0 and 2 in the tests)
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[0]->pressure_n, 5341.001);
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[2]->pressure_n, 8011.501);


	// Check the n1 setflow works:
	// Set some fake flow data ( this may be unnecessary now.. )
	double flows[5] = {10.0, 15.0, 0.0, 0.0, 0.0};
	
	boundaryConditionManager_instance->updateAllRCRS_setflow_n(flows);
	
	boundaryConditionManager_instance->updateAllRCRS_setflow_n1(flows);
	
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[0]->flow_n1, 10.0);
	EXPECT_DOUBLE_EQ((*retrievedBoundaryConditions)[2]->flow_n1, 15.0);
}
