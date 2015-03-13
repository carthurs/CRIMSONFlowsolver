#include "gtest/gtest.h"
#include "testMultidom.hxx"
#include "RCR.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "indexShifters.hxx"
#include <stack>
#include <boost/weak_ptr.hpp>

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMultidom() { return 0; }

// Tests that the correct types of boundary conditions were made, with the correct given surface indices.
TEST_F(testMultidom, checkBoundaryConditionsMadeProperly) {
  // Check we got the right boundary conditions
	// EXPECT_TRUE(1==1);
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[0])==typeid(RCR));
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[1])==typeid(NetlistBoundaryCondition));
  EXPECT_TRUE(typeid(*(*retrievedBoundaryConditions)[2])==typeid(RCR));

  EXPECT_EQ((*retrievedBoundaryConditions)[0]->surfaceIndex,3);
  EXPECT_EQ((*retrievedBoundaryConditions)[1]->surfaceIndex,7);
  EXPECT_EQ((*retrievedBoundaryConditions)[2]->surfaceIndex,9);

  EXPECT_EQ(*((*retrievedBoundaryConditions)[0]->flow_n_ptrs.at(0)),flow1);
  // EXPECT_EQ(*((*retrievedBoundaryConditions)[1]->flow_n_ptr),flow2); // uncomment once Netlist initialisation sets this pointer
  EXPECT_EQ(*((*retrievedBoundaryConditions)[2]->flow_n_ptrs.at(0)),flow3);

  EXPECT_EQ(*((*retrievedBoundaryConditions)[0]->pressure_n_ptrs.at(0)),press1);
  // EXPECT_EQ(*((*retrievedBoundaryConditions)[1]->pressure_n_ptr),press2); // uncomment once Netlist initialisation sets this pointer
  EXPECT_EQ(*((*retrievedBoundaryConditions)[2]->pressure_n_ptrs.at(0)),press3);
  
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
	(*retrievedBoundaryConditions)[0]->flow_n_ptrs.at(0) = &fakeFlow1;
	
	(*retrievedBoundaryConditions)[2]->dp_dq_n1 = 534.1;
	(*retrievedBoundaryConditions)[2]->Hop_n1 = 0.001;
	double fakeFlow2 = 15.0;
	(*retrievedBoundaryConditions)[2]->flow_n_ptrs.at(0) = &fakeFlow2;
	
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

TEST_F(testMultidom, checkNetlistComponentNeighbourPointers)
{
	// Get the netlist boundary condition
	NetlistBoundaryCondition* netlistBC = dynamic_cast<NetlistBoundaryCondition*>((*retrievedBoundaryConditions).at(1).get());

	// Find the index of the component connected to the 3D domain, so we can use it as a starting point for walking the component tree::
	int indexOfComponentAt3DInterface;
	for (auto component = netlistBC->mp_NetlistCircuit->mp_CircuitDescription->components.begin(); component!= netlistBC->mp_NetlistCircuit->mp_CircuitDescription->components.end(); component++)
	{
		if ((*component)->prescribedFlowType == Flow_3DInterface)
		{
			indexOfComponentAt3DInterface = (*component)->indexInInputData;
		}
	}

	// Beginning from the component at the 3D interface, walk the component tree.
	// First, some data to check against:
	std::stack<circuit_component_t> expectedComponentTypes;
	expectedComponentTypes.push(Component_Resistor);
	expectedComponentTypes.push(Component_Capacitor);
	expectedComponentTypes.push(Component_Resistor);
	std::stack<int> expectedEndNodeNeighbourCounts;
	expectedEndNodeNeighbourCounts.push(0);
	expectedEndNodeNeighbourCounts.push(0);
	expectedEndNodeNeighbourCounts.push(2);
	std::stack<double> expectedParameterValues;
	expectedParameterValues.push(1600.0);
	expectedParameterValues.push(0.001278473);
	expectedParameterValues.push(58.43089);

	// Put the 3D interface component on a stack of components that need walking from:
	std::stack<boost::weak_ptr<CircuitComponent>> componentsNeedingChecking;
	boost::weak_ptr<CircuitComponent> toPushOntoStack(netlistBC->mp_NetlistCircuit->mp_CircuitDescription->mapOfComponents.at(indexOfComponentAt3DInterface));
	componentsNeedingChecking.push(toPushOntoStack);
	// To keep track of which components have been already checked:
	std::vector<bool> componentsWhichHaveBeenChecked(netlistBC->mp_NetlistCircuit->mp_CircuitDescription->components.size(),false);
	while(!componentsNeedingChecking.empty())
	{
		// pop the stack to set the current component:
		boost::weak_ptr<CircuitComponent> currentComponent(componentsNeedingChecking.top());
		componentsNeedingChecking.pop();
		// For the current component:
		// Ensure we've not done this component yet (avoids circular problems)
		if (componentsWhichHaveBeenChecked.at(toZeroIndexing(currentComponent.lock()->indexInInputData)) == false)
		{
			// note that we're checking this component:
			componentsWhichHaveBeenChecked.at(toZeroIndexing(currentComponent.lock()->indexInInputData)) = true;
			// check the type is as expected
			EXPECT_TRUE(currentComponent.lock()->getType() == expectedComponentTypes.top());
			expectedComponentTypes.pop();
			// check the number of neighbours is as expected
			EXPECT_EQ(currentComponent.lock()->neighbouringComponentsAtEndNode.size(), expectedEndNodeNeighbourCounts.top());
			expectedEndNodeNeighbourCounts.pop();
			// Put all the neighbours on a stack
			for (auto neighbouringComponent=currentComponent.lock()->neighbouringComponentsAtEndNode.begin(); neighbouringComponent!=currentComponent.lock()->neighbouringComponentsAtEndNode.end(); neighbouringComponent++)
			{
				componentsNeedingChecking.push(*neighbouringComponent);
			}
		}
	}
}

TEST_F(testMultidom, checkClosedDiodeWithRemainingOpenPathDetected)
{
	// Get the netlist boundary condition
	NetlistBoundaryCondition* netlistBC_thirdNetlist = dynamic_cast<NetlistBoundaryCondition*>((*retrievedBoundaryConditions).at(5).get());
	EXPECT_TRUE(netlistBC_thirdNetlist->mp_NetlistCircuit->mp_CircuitDescription->m_flowPermittedAcross3DInterface);
}

TEST_F(testMultidom, checkClosedDiodeWithoutRemainingOpenPathDetected)
{
	// Get the netlist boundary condition
	NetlistBoundaryCondition* netlistBC_fourthNetlist = dynamic_cast<NetlistBoundaryCondition*>((*retrievedBoundaryConditions).at(6).get());
	EXPECT_FALSE(netlistBC_fourthNetlist->mp_NetlistCircuit->mp_CircuitDescription->m_flowPermittedAcross3DInterface);
}