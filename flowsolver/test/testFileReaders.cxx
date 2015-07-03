#include "testFileReaders.hxx"
#include "gtest/gtest.h"
#include "datatypesInCpp.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibrary() { return 0; }


// Tests that the fileReader reads the first entry of the test-case rcrt.dat
TEST_F(testFileReaders, checkFileReaderGotFirstRCR) {
  rcrtReader* rcrtReader_instance = rcrtReader::Instance();

  int pdmax = rcrtReader_instance->getPdmax();
  EXPECT_EQ(3,pdmax);

  // First RCR Params
  int numberOfRCRDatapoints1 = rcrtReader_instance->getNumDataRCR()[0];
  EXPECT_EQ(2,numberOfRCRDatapoints1);
  
  double r1 = rcrtReader_instance->getR1()[0];
  EXPECT_DOUBLE_EQ(3.815889462107935e-01, r1);

  double c =  rcrtReader_instance->getC()[0];
  EXPECT_DOUBLE_EQ(1.885526988247284e-01, c);

  double r2 = rcrtReader_instance->getR2()[0];
  EXPECT_DOUBLE_EQ(1.740713618386789e+00, r2);

  std::vector<std::pair<double,double>> timeData1 = rcrtReader_instance->getTimeDataPdist()[0];
  EXPECT_DOUBLE_EQ(0.0,timeData1[0].first);
  EXPECT_DOUBLE_EQ(0.0,timeData1[0].second);
  EXPECT_DOUBLE_EQ(1.1,timeData1[1].first);
  EXPECT_DOUBLE_EQ(0.0,timeData1[1].second);

}

// Tests that the fileReader reads the second entry of the test-case rcrt.dat
TEST_F(testFileReaders, checkFileReaderGotSecondRCR) {
  rcrtReader* rcrtReader_instance = rcrtReader::Instance();

  // Second RCR Params
  int numberOfRCRDatapoints2 = rcrtReader_instance->getNumDataRCR()[1];
  EXPECT_EQ(3,numberOfRCRDatapoints2);
  
  double r1 = rcrtReader_instance->getR1()[1];
  EXPECT_DOUBLE_EQ(44.123, r1);

  double c =  rcrtReader_instance->getC()[1];
  EXPECT_DOUBLE_EQ(55.456, c);

  double r2 = rcrtReader_instance->getR2()[1];
  EXPECT_DOUBLE_EQ(66.567, r2);

  std::vector<std::vector<std::pair<double,double>>> timeData2v = rcrtReader_instance->getTimeDataPdist();
  std::vector<std::pair<double,double>> timeData2 = timeData2v[1];
  EXPECT_DOUBLE_EQ(0.0,timeData2[0].first);
  EXPECT_DOUBLE_EQ(0.0,timeData2[0].second);
  EXPECT_DOUBLE_EQ(0.5,timeData2[1].first);
  EXPECT_DOUBLE_EQ(1.7,timeData2[1].second);
  EXPECT_DOUBLE_EQ(0.9,timeData2[2].first);
  EXPECT_DOUBLE_EQ(2.3,timeData2[2].second);
}

TEST_F(testFileReaders, checkControlledCoronaryReader) {
  controlledCoronaryReader* controlledCoronaryReader_instance = controlledCoronaryReader::Instance();
  controlledCoronaryReader_instance->setFileName("controlled_coronaries_test.dat");
  controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile();

  std::vector<double> returnedVector;

  returnedVector = controlledCoronaryReader_instance->getResistanceNearAorta();
  EXPECT_DOUBLE_EQ(1.286208e5,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(2.286208e5,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getComplianceNearAorta();
  EXPECT_DOUBLE_EQ(4.5e-7,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(5.5e-7,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getMidResistance();
  EXPECT_DOUBLE_EQ(1.286208e5,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(2.286208e5,returnedVector.at(1));
  
  returnedVector = controlledCoronaryReader_instance->getIntramyocardialCompliance();
  EXPECT_DOUBLE_EQ(2.7e-6,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(3.7e-6,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getDistalResistance();
  EXPECT_DOUBLE_EQ(7.5e5,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(8.5e5,returnedVector.at(1));
  
  returnedVector = controlledCoronaryReader_instance->getMinimumAllowedResistance();
  EXPECT_DOUBLE_EQ(10.0,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(20.0,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getMaximumAllowedResistance();
  EXPECT_DOUBLE_EQ(1e8,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(2e8,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getPerfusionBedMVO2_previous();
  EXPECT_DOUBLE_EQ(0.0272,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.0272,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getPerfusionBedMVO2_current();
  EXPECT_DOUBLE_EQ(0.0273,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.0273,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getProportionOfMyocardiumPerfusedByThisSurface();
  EXPECT_DOUBLE_EQ(0.1,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.1,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getMetabolicFeedbackGain();
  EXPECT_DOUBLE_EQ(0.31,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.31,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getAlphaAdrenergicFeedforwardGain();
  EXPECT_DOUBLE_EQ(0.0,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.0,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getBetaAdrenergicFeedforwardGain();
  EXPECT_DOUBLE_EQ(0.0,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.0,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getFeedbackDamping();
  EXPECT_DOUBLE_EQ(0.26,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(1.26,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getO2DemandIntegrationWindow();
  EXPECT_DOUBLE_EQ(5.0,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(6.0,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getCapacitorNearAortaTopPressure();
  EXPECT_DOUBLE_EQ(10000.0,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(20000.0,returnedVector.at(1));

  returnedVector = controlledCoronaryReader_instance->getIntramyocardialCapacitorTopPressure();
  EXPECT_DOUBLE_EQ(11000.0,returnedVector.at(0));
  EXPECT_DOUBLE_EQ(21000.0,returnedVector.at(1));


  controlledCoronaryReader_instance->Term();
}

TEST_F(testFileReaders, checkControlledCoronaryReaderMalformedFileDetection) {
  controlledCoronaryReader* controlledCoronaryReader_instance = controlledCoronaryReader::Instance();
  controlledCoronaryReader_instance->setFileName("controlled_coronaries_test_malformed.dat");
  
  EXPECT_ANY_THROW(controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile());

  controlledCoronaryReader_instance->Term();
}

TEST_F(testFileReaders, checkHistFileReader) {
  
  histFileReader histFileReader_instance = histFileReader();

  histFileReader_instance.setFileName("PressHist.dat");
  histFileReader_instance.setNumColumns(2);
  histFileReader_instance.readFileInternalMetadata();
  histFileReader_instance.readAndSplitMultiSurfaceRestartFile();

  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double readResult = histFileReader_instance.getReadFileData(1,5);

  EXPECT_EQ(10645.5858581080,readResult);
}

TEST_F(testFileReaders, checkNetlistReader) {
  NetlistReader* netlistReader_instance = NetlistReader::Instance();
  netlistReader_instance->setFileName("netlist_surfaces.dat");
  netlistReader_instance->readAndSplitMultiSurfaceInputFile();

  {
    std::vector<double> returnedVectorOfDoubles;
    returnedVectorOfDoubles = netlistReader_instance->getComponentParameterValues(0);
    EXPECT_EQ(returnedVectorOfDoubles.at(0),58.43089e0);
    EXPECT_EQ(returnedVectorOfDoubles.at(1),0.001278473e0);
    EXPECT_EQ(returnedVectorOfDoubles.at(2),1600.0e0);
  }

  {
    std::vector<double> returnedVectorOfDoubles;
    returnedVectorOfDoubles = netlistReader_instance->getComponentParameterValues(1);
    EXPECT_EQ(returnedVectorOfDoubles.at(0),0.001278473e1);
    EXPECT_EQ(returnedVectorOfDoubles.at(1),60.43089e1);
    EXPECT_EQ(returnedVectorOfDoubles.at(2),1600.0e1);
    EXPECT_EQ(returnedVectorOfDoubles.at(3),1e-4);
    EXPECT_EQ(returnedVectorOfDoubles.at(4),1e2);
  }

  {
    double initialVolume = netlistReader_instance->getComponentInitialVolume(1,2);
    EXPECT_EQ(initialVolume, 130000.0);
  }

  {
    double initialVolume = netlistReader_instance->getComponentInitialVolume(1,4);
    EXPECT_EQ(initialVolume, 140000.0);
  }

  {
    EXPECT_THROW(netlistReader_instance->getComponentInitialVolume(1,3), std::runtime_error);
  }

  std::vector<std::vector<double>> returnedVectorOfDoubleVectors = netlistReader_instance->getValueOfPrescribedPressures();
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(0).at(0),0.01e0);
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(0).at(1),0.11e0);
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(1).at(0),1.1e0);
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(1).at(1),1.2e0);
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(1).at(2),1.3e0);
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(1).at(3),1.4e0);

  returnedVectorOfDoubleVectors = netlistReader_instance->getValueOfPrescribedFlows();
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(0).at(0),-1.0e0);
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(1).at(0),-1.1e0);

  std::vector<std::map<int,double>> returnedVectorOfIntToDoubleMaps = netlistReader_instance->getInitialPressures();
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(1),206663.0e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(2),206663.1e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(3),0.0e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(4),0.1e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(1).at(1),106663.0e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(1).at(2),106663.1e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(1).at(3),106663.2e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(1).at(4),1.0e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(1).at(5),1.1e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(1).at(6),1.2e0);


  std::vector<int> returnedVectorOfInts;
  returnedVectorOfInts = netlistReader_instance->getNumberOfPressureNodes();
  EXPECT_EQ(returnedVectorOfInts.at(0),4);
  EXPECT_EQ(returnedVectorOfInts.at(1),6);

  returnedVectorOfInts = netlistReader_instance->getNumberOfComponents();
  EXPECT_EQ(returnedVectorOfInts.at(0),3);
  EXPECT_EQ(returnedVectorOfInts.at(1),5);

  returnedVectorOfInts = netlistReader_instance->getNumberOfPrescribedPressures();
  EXPECT_EQ(returnedVectorOfInts.at(0),2);
  EXPECT_EQ(returnedVectorOfInts.at(1),4);

  returnedVectorOfInts = netlistReader_instance->getNumberOfPrescribedFlows();
  EXPECT_EQ(returnedVectorOfInts.at(0),1);
  EXPECT_EQ(returnedVectorOfInts.at(1),1);
  
  std::vector<std::vector<int>> returnedVectorOfIntVectors;
  returnedVectorOfIntVectors = netlistReader_instance->getComponentStartNodes();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),1);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(1),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(2),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(0),3);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(1),1);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(2),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(3),4);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(4),4);

  returnedVectorOfIntVectors = netlistReader_instance->getComponentEndNodes();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(1),3);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(2),4);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(0),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(1),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(2),4);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(3),5);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(4),6);

  returnedVectorOfIntVectors = netlistReader_instance->getListOfPrescribedPressures();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),3);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(1),4);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(0),3);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(1),4);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(2),5);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(3),6);

  returnedVectorOfIntVectors = netlistReader_instance->getListOfPrescribedFlows();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),1);
  EXPECT_EQ(returnedVectorOfIntVectors.at(1).at(0),1);


  std::vector<std::vector<circuit_nodal_pressure_prescription_t>> returnedVectorOfPressurePrescriptions;
  returnedVectorOfPressurePrescriptions = netlistReader_instance->getTypeOfPrescribedPressures();
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(0).at(0) == Pressure_Fixed);
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(0).at(1) == Pressure_Fixed);
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(1).at(0) == Pressure_Fixed);
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(1).at(1) == Pressure_Fixed);
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(1).at(2) == Pressure_Fixed);
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(1).at(3) == Pressure_Fixed);

  std::vector<std::vector<circuit_component_flow_prescription_t>> returnedVectorOfFlowPrescriptions;
  returnedVectorOfFlowPrescriptions = netlistReader_instance->getTypeOfPrescribedFlows();
  EXPECT_TRUE(returnedVectorOfFlowPrescriptions.at(0).at(0) == Flow_3DInterface);
  EXPECT_TRUE(returnedVectorOfFlowPrescriptions.at(1).at(0) == Flow_3DInterface);

  std::vector<std::vector<circuit_component_t>> returnedVectorOfComponentTypes;
  returnedVectorOfComponentTypes = netlistReader_instance->getComponentTypes();
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(0).at(0) == Component_Resistor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(0).at(1) == Component_Capacitor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(0).at(2) == Component_Resistor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(1).at(0) == Component_Resistor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(1).at(1) == Component_Capacitor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(1).at(2) == Component_VolumeTracking);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(1).at(3) == Component_Inductor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(1).at(4) == Component_VolumeTrackingPressureChamber);

  netlistReader_instance->Term();
}


TEST_F(testFileReaders, checkNetlistDownstreamCircuitReader)
{
  NetlistDownstreamCircuitReader* downstreamReader_instance = NetlistDownstreamCircuitReader::Instance();
  downstreamReader_instance->setFileName("netlist_closed_loop_downstream.dat");
  downstreamReader_instance->readAndSplitMultiSurfaceInputFile();

  {
    std::vector<double> returnedVectorOfDoubles;
    returnedVectorOfDoubles = downstreamReader_instance->getComponentParameterValues(0);
    EXPECT_EQ(returnedVectorOfDoubles.at(0),58.43089e0);
    EXPECT_EQ(returnedVectorOfDoubles.at(1),0.001278473e0);
    EXPECT_EQ(returnedVectorOfDoubles.at(2),1400.0e0);
  }

  std::vector<std::vector<double>> returnedVectorOfDoubleVectors = downstreamReader_instance->getValueOfPrescribedPressures();
  EXPECT_EQ(returnedVectorOfDoubleVectors.at(0).at(0),0.01e0);

  // returnedVectorOfDoubleVectors = downstreamReader_instance->getValueOfPrescribedFlows();
  // EXPECT_EQ(returnedVectorOfDoubleVectors.at(0).at(0),-1.0e0);

  std::vector<std::map<int,double>> returnedVectorOfIntToDoubleMaps = downstreamReader_instance->getInitialPressures();
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(1),133.0e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(2),133.1e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(3),0.1e0);
  EXPECT_EQ(returnedVectorOfIntToDoubleMaps.at(0).at(4),133.2e0);

  std::vector<int> returnedVectorOfInts;
  returnedVectorOfInts = downstreamReader_instance->getNumberOfPressureNodes();
  EXPECT_EQ(returnedVectorOfInts.at(0),4);

  returnedVectorOfInts = downstreamReader_instance->getNumberOfComponents();
  EXPECT_EQ(returnedVectorOfInts.at(0),3);

  returnedVectorOfInts = downstreamReader_instance->getNumberOfPrescribedPressures();
  EXPECT_EQ(returnedVectorOfInts.at(0),1);

  // returnedVectorOfInts = downstreamReader_instance->getNumberOfPrescribedFlows();
  // EXPECT_EQ(returnedVectorOfInts.at(0),1);
  // EXPECT_EQ(returnedVectorOfInts.at(1),1);
  
  std::vector<std::vector<int>> returnedVectorOfIntVectors;
  returnedVectorOfIntVectors = downstreamReader_instance->getComponentStartNodes();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),1);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(1),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(2),2);

  returnedVectorOfIntVectors = downstreamReader_instance->getComponentEndNodes();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),2);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(1),3);
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(2),4);

  returnedVectorOfIntVectors = downstreamReader_instance->getListOfPrescribedPressures();
  EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),3);

  // returnedVectorOfIntVectors = downstreamReader_instance->getListOfPrescribedFlows();
  // EXPECT_EQ(returnedVectorOfIntVectors.at(0).at(0),1);


  std::vector<std::vector<circuit_nodal_pressure_prescription_t>> returnedVectorOfPressurePrescriptions;
  returnedVectorOfPressurePrescriptions = downstreamReader_instance->getTypeOfPrescribedPressures();
  EXPECT_TRUE(returnedVectorOfPressurePrescriptions.at(0).at(0) == Pressure_Fixed);

  // std::vector<std::vector<circuit_component_flow_prescription_t>> returnedVectorOfFlowPrescriptions;
  // returnedVectorOfFlowPrescriptions = downstreamReader_instance->getTypeOfPrescribedFlows();
  // EXPECT_TRUE(returnedVectorOfFlowPrescriptions.at(0).at(0) == Flow_3DInterface);
  // EXPECT_TRUE(returnedVectorOfFlowPrescriptions.at(1).at(0) == Flow_3DInterface);

  std::vector<std::vector<circuit_component_t>> returnedVectorOfComponentTypes;
  returnedVectorOfComponentTypes = downstreamReader_instance->getComponentTypes();
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(0).at(0) == Component_Resistor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(0).at(1) == Component_Capacitor);
  EXPECT_TRUE(returnedVectorOfComponentTypes.at(0).at(2) == Component_Resistor);

  int downstreamCircuitIndex = 0;
  int returnedInteger = downstreamReader_instance->getNumberOfBoundaryConditionsConnectedTo(downstreamCircuitIndex);
  EXPECT_EQ(returnedInteger,3);

  returnedVectorOfInts = downstreamReader_instance->getConnectedCircuitSurfaceIndices(downstreamCircuitIndex);
  EXPECT_EQ(returnedVectorOfInts.at(0),5);
  EXPECT_EQ(returnedVectorOfInts.at(1),6);
  EXPECT_EQ(returnedVectorOfInts.at(2),7);

  returnedVectorOfInts = downstreamReader_instance->getLocalBoundaryConditionInterfaceNodes(downstreamCircuitIndex);
  EXPECT_EQ(returnedVectorOfInts.at(0),1);
  EXPECT_EQ(returnedVectorOfInts.at(1),4);
  EXPECT_EQ(returnedVectorOfInts.at(2),4);

  returnedVectorOfInts = downstreamReader_instance->getRemoteBoundaryConditionInterfaceNodes(downstreamCircuitIndex);
  EXPECT_EQ(returnedVectorOfInts.at(0),5);
  EXPECT_EQ(returnedVectorOfInts.at(1),4);
  EXPECT_EQ(returnedVectorOfInts.at(2),3);

  downstreamReader_instance->Term();
}