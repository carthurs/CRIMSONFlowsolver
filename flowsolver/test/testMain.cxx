#include "testMain.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMain() { return 0; }

TEST_F(testMain, checkSimpleShortSimulationWithRCRs) {
  
  std::string simDir = "mainTests/basic";
  setSimDirectoryAndClearoutOldFiles(simDir);

  runSimulation();

  histFileReader* QHistReader = new histFileReader();
  QHistReader->setFileName("QHistRCR.dat");
  QHistReader->setNumColumns(2);
  QHistReader->readAndSplitMultiSurfaceRestartFile();

  double finalQhistRCRValue = ((QHistReader->dataReadFromFile).at(5))[1];
  EXPECT_NEAR(finalQhistRCRValue,0.7149210791E+03,1e-7);

  delete QHistReader;
}