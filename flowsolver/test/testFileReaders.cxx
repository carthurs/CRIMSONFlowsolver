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
