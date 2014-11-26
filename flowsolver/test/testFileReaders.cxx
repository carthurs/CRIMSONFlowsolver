#include "testFileReaders.hxx"
#include "gtest/gtest.h"

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


  controlledCoronaryReader_instance->Term();
}

TEST_F(testFileReaders, checkControlledCoronaryReaderMalformedFileDetection) {
  controlledCoronaryReader* controlledCoronaryReader_instance = controlledCoronaryReader::Instance();
  controlledCoronaryReader_instance->setFileName("controlled_coronaries_test_malformed.dat");
  
  std::cout << "This test expects an error message on the next line:" << std::endl;
  EXPECT_ANY_THROW(controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile());

  controlledCoronaryReader_instance->Term();
}