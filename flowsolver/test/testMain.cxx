#include "testMain.hxx"
#include "debuggingToolsForCpp.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMain() { return 0; }

TEST_F(testMain, checkRCRSimpleShortSimulation) {
  
  std::string simDir = "mainTests/basic";
  setSimDirectory(simDir);
  clearOutOldFiles();

  runSimulation();

  histFileReader* QHistReader = new histFileReader();
  QHistReader->setFileName("QHistRCR.dat");
  QHistReader->setNumColumns(2);
  QHistReader->readAndSplitMultiSurfaceRestartFile();

  double finalQHistRCRValue = ((QHistReader->dataReadFromFile).at(5))[1];
  EXPECT_NEAR(finalQHistRCRValue,714.921079082528,1e-7);

  delete QHistReader;
}


TEST_F(testMain, checkRestartWorks_RCRSimpleShortSimulation) {

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

TEST_F(testMain, checkCoronarySimpleShortSimulation) {
  std::string simDir = "mainTests/coronary/completeLPN";
  setSimDirectory(simDir);
  clearOutOldFiles();

  runSimulation();

  // Check a value in FlowHist.dat
  histFileReader FlowHistReader = histFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(2);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double finalFlowHistValue = FlowHistReader.getReadFileData(1,5);
  EXPECT_NEAR(finalFlowHistValue,3.628656482154824E-003,1e-10);

  // Check a value in PressHist.dat
  histFileReader PressHistReader = histFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(2);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double finalPressHistValue = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(finalPressHistValue,10772.5198465269,1e-6);
}

TEST_F(testMain, checkCoronaryCanEmulateKnownRCRResults) {
  std::string simDir = "mainTests/coronary/emulateRCR";
  setSimDirectory(simDir);
  clearOutOldFiles();

  runSimulation();

  // Check PressHist.dat
  histFileReader PressHistReader = histFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(2);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double readResult = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(10645.5858581080,readResult,1e-7);

  // Check FlowHist.dat
  histFileReader FlowHistReader = histFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(2);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  readResult = FlowHistReader.getReadFileData(1,5);
  EXPECT_NEAR(478.116982120136,readResult,1e-7);
}

TEST_F(testMain, checkNetlistCanEmulateKnownRCRResults) {
  std::string simDir = "mainTests/netlist/emulateRCR";
  setSimDirectory(simDir);
  clearOutOldFiles();

  runSimulation();

  // Check PressHist.dat
  histFileReader PressHistReader = histFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(2);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double readResult = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(10645.5858581080,readResult,1e-8);

  // Check FlowHist.dat
  histFileReader FlowHistReader = histFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(2);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  readResult = FlowHistReader.getReadFileData(1,5);
  EXPECT_NEAR(478.116982120136,readResult,1e-7);
}