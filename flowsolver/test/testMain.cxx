#include "testMain.hxx"
#include "debuggingToolsForCpp.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMain() { return 0; }

TEST_F(testMain, checkSimpleShortSimulationWithRCRs) {
  
  std::string simDir = "mainTests/basic";
  setSimDirectory(simDir);
  clearOutOldFiles();

  runSimulation();

  histFileReader* QHistReader = new histFileReader();
  QHistReader->setFileName("QHistRCR.dat");
  QHistReader->setNumColumns(2);
  QHistReader->readAndSplitMultiSurfaceRestartFile();

  double finalQHistRCRValue = ((QHistReader->dataReadFromFile).at(5))[1];
  EXPECT_DOUBLE_EQ(finalQHistRCRValue,714.921079082528);

  delete QHistReader;
}


TEST_F(testMain, checkRestartWorks_simpleShortSimulationWithRCRs) {

  std::string simDir = "mainTests/restart";
  setSimDirectory(simDir);
  clearOutOldFiles();
  
  runSimulation();
  histFileReader* PHistReader = new histFileReader();
  PHistReader->setFileName("PHistRCR.dat");
  PHistReader->setNumColumns(2);
  PHistReader->readAndSplitMultiSurfaceRestartFile();

  double finalPHistRCRValue = ((PHistReader->dataReadFromFile).at(10))[1];

  EXPECT_NEAR(finalPHistRCRValue,0.1162541088E+05,10e-5);

  delete PHistReader;
}