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

  double finalQHistRCRValue = ((QHistReader->m_dataReadFromFile).at(5))[1];
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

  double finalPHistRCRValue = ((PHistReader->m_dataReadFromFile).at(10))[1];

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
  EXPECT_NEAR(finalPressHistValue,10772.5198465269,1e-5);
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

//\todo reinstate or remove:
// TEST_F(testMain, checkNetlistSimpleDiode) {
//   std::string simDir = "mainTests/netlist/simpleDiodeTest";
//   setSimDirectory(simDir);
//   clearOutOldFiles();

//   // runSimulation();

//   // Check PressHist.dat
//   // histFileReader PressHistReader = histFileReader();
//   // PressHistReader.setFileName("PressHist.dat");
//   // PressHistReader.setNumColumns(2);
//   // PressHistReader.readFileInternalMetadata();
//   // PressHistReader.readAndSplitMultiSurfaceRestartFile();
//   // // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
//   // double readResult = PressHistReader.getReadFileData(1,5);
//   // EXPECT_NEAR(10645.5858581080,readResult,1e-8);

//   // // Check FlowHist.dat
//   // histFileReader FlowHistReader = histFileReader();
//   // FlowHistReader.setFileName("FlowHist.dat");
//   // FlowHistReader.setNumColumns(2);
//   // FlowHistReader.readFileInternalMetadata();
//   // FlowHistReader.readAndSplitMultiSurfaceRestartFile();
//   // // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
//   // readResult = FlowHistReader.getReadFileData(1,5);
//   // EXPECT_NEAR(478.116982120136,readResult,1e-7);
// }

TEST_F(testMain, checkPreKalmanPreGlobalNodeNumberingGeombcRuns) {
  std::string simDir = "mainTests/legacy/preKalmanPreGlobalNodenumbering";
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
  double pressHistResult = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(106420.076723820,pressHistResult,1e-8);

  // Check FlowHist.dat
  histFileReader FlowHistReader = histFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(2);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double flowHistResult = FlowHistReader.getReadFileData(1,5);
  EXPECT_NEAR(1.077519395498929e-2,flowHistResult,1e-7);
}

TEST_F(testMain, checkNetlistHeartModel) {
  // This test uses a solver.inp which (on purpose) does not take enough
  // iterations (Step Construction 0 1 0 1...) for decent convergence.
  // The reason for this is that the inaccuracy causes the aortic valve to flap
  // open and closed during a very short test - this is fine for us as we
  // just want to test that the valve is doing its job!
  setSimDirectory("mainTests/netlist/heart");
  clearOutOldFiles();

  runSimulation();

  // Check PressHist.dat
  histFileReader PressHistReader = histFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(3);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  
  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double pressHistResult = PressHistReader.getReadFileData(0,6);
  EXPECT_NEAR(1003.49069102218,pressHistResult,1e-5);
  // ...second column
  pressHistResult = PressHistReader.getReadFileData(1,6);
  EXPECT_NEAR(968.660015172370,pressHistResult,1e-6);
  // ... third column
  pressHistResult = PressHistReader.getReadFileData(2,6);
  EXPECT_NEAR(964.764355883842,pressHistResult,1e-6);

  // Check FlowHist.dat
  histFileReader FlowHistReader = histFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(3);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double flowHistResult = FlowHistReader.getReadFileData(0,6);
  EXPECT_NEAR(0.0,flowHistResult,1e-7);
  // ... 2nd column:
  flowHistResult = FlowHistReader.getReadFileData(1,6);
  EXPECT_NEAR(73.4892153364009,flowHistResult,1e-5);
  // ...third column:
  flowHistResult = FlowHistReader.getReadFileData(2,6);
  EXPECT_NEAR(-62.2259850634413,flowHistResult,1e-5);
}

TEST_F(testMain, checkClosedLoopWithHeart) {
  setSimDirectory("mainTests/netlist/closedLoopHeart");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check PressHist.dat
  {
    histFileReader zeroDDomainPressures = histFileReader();
    zeroDDomainPressures.setFileName("PressHist.dat");
    zeroDDomainPressures.setNumColumns(3);
    zeroDDomainPressures.readFileInternalMetadata();
    zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();

    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = zeroDDomainPressures.getReadFileData(0,6);
    EXPECT_NEAR(10646.3794331755,pressureResult,1e-5);
    // ...second column
    pressureResult = zeroDDomainPressures.getReadFileData(1,6);
    EXPECT_NEAR(10646.3789516346,pressureResult,1e-6);
    // ... third column
    pressureResult = zeroDDomainPressures.getReadFileData(2,6);
    EXPECT_NEAR(10646.3708851326,pressureResult,1e-6);
  }

  // Check FlowHist.dat
  {
    histFileReader zeroDDomainFlows = histFileReader();
    zeroDDomainFlows.setFileName("FlowHist.dat");
    zeroDDomainFlows.setNumColumns(3);
    zeroDDomainFlows.readFileInternalMetadata();
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double flowResult = zeroDDomainFlows.getReadFileData(0,6);
    EXPECT_NEAR(0.000000000000000E+000,flowResult,1e-8);
    // ...second column
    flowResult = zeroDDomainFlows.getReadFileData(1,6);
    EXPECT_NEAR(0.328453445161543,flowResult,1e-4);
    // ... third column
    flowResult = zeroDDomainFlows.getReadFileData(2,6);
    EXPECT_NEAR(6.117932527654634E-002,flowResult,1e-5);
  }

  // Check netlistFlows_surface_5.dat
  {
    histFileReader zeroDDomainFlows = histFileReader();
    zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
    zeroDDomainFlows.setNumColumns(6);
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = zeroDDomainFlows.getReadFileData(1,5);
    EXPECT_NEAR(5.62122509164454e-305,pressureResult,1e-8);
  }

  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
  {
    histFileReader closedLoopDownstreamPressures = histFileReader();
    closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
    closedLoopDownstreamPressures.setNumColumns(5);
    closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

    double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,5);
    EXPECT_NEAR(533.423900130907,pressureResult,1e-8);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(2,5);
    EXPECT_NEAR(533.423900130907,pressureResult,1e-8);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(3,5);
    EXPECT_NEAR(0.0,pressureResult,1e-8);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(4,5);
    EXPECT_NEAR(4580.01154285495,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
    histFileReader closedLoopDownstreamPressures = histFileReader();
    closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
    closedLoopDownstreamPressures.setNumColumns(4);
    closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

    double flowResult = closedLoopDownstreamPressures.getReadFileData(1,5);
    EXPECT_NEAR(6.25264122964437e-308,flowResult,1e-8);

    flowResult = closedLoopDownstreamPressures.getReadFileData(2,5);
    EXPECT_NEAR(40465.8764272404,flowResult,1e-8);

    flowResult = closedLoopDownstreamPressures.getReadFileData(3,5);
    EXPECT_NEAR(-40465.8764272404,flowResult,1e-8);
  }
}

TEST_F(testMain, checkMixedNetlistAndRCRT)
{
  setSimDirectory("mainTests/legacy/netlistWithRCRTs");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

   // Check PressHist.dat
  {
    histFileReader pressHistData = histFileReader();
    pressHistData.setFileName("PressHist.dat");
    pressHistData.setNumColumns(3);
    pressHistData.readFileInternalMetadata();
    pressHistData.readAndSplitMultiSurfaceRestartFile();

    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = pressHistData.getReadFileData(0,6);
    EXPECT_NEAR(790.498746283900,pressureResult,1e-5);
    // ...second column
    pressureResult = pressHistData.getReadFileData(1,6);
    EXPECT_NEAR(7653.31439101382,pressureResult,1e-6);
    // ... third column
    pressureResult = pressHistData.getReadFileData(2,6);
    EXPECT_NEAR(9376.67945078140,pressureResult,1e-6);
  }

  // Check FlowHist.dat
  {
    histFileReader flowHistData = histFileReader();
    flowHistData.setFileName("FlowHist.dat");
    flowHistData.setNumColumns(3);
    flowHistData.readFileInternalMetadata();
    flowHistData.readAndSplitMultiSurfaceRestartFile();
    
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double flowResult = flowHistData.getReadFileData(0,6);
    EXPECT_NEAR(137208.420657606,flowResult,1e-8);
    // ...second column
    flowResult = flowHistData.getReadFileData(1,6);
    EXPECT_NEAR(-97144.4580441544,flowResult,1e-4);
    // ... third column
    flowResult = flowHistData.getReadFileData(2,6);
    EXPECT_NEAR(-40084.8691875509,flowResult,1e-5);
  }

  // Check netlistFlows_surface_5.dat
  {
    histFileReader zeroDDomainFlows = histFileReader();
    zeroDDomainFlows.setFileName("netlistPressures_surface_5.dat");
    zeroDDomainFlows.setNumColumns(7);
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = zeroDDomainFlows.getReadFileData(1,5);
    EXPECT_NEAR(798.293386775323,pressureResult,1e-8);
  }

  {
    histFileReader* PHistReader = new histFileReader();
    PHistReader->setFileName("PHistRCR.dat");
    PHistReader->setNumColumns(3);
    PHistReader->readAndSplitMultiSurfaceRestartFile();

    double finalPHistRCRValue = PHistReader->getReadFileData(1,5);//((PHistReader->m_dataReadFromFile).at(5))[1];
    EXPECT_NEAR(finalPHistRCRValue,7664.81545745205,10e-5);

    finalPHistRCRValue = PHistReader->getReadFileData(2,5);//((PHistReader->m_dataReadFromFile).at(5))[2];
    EXPECT_NEAR(finalPHistRCRValue,9409.53605840001,10e-5);

    delete PHistReader;
  }  

  {
    histFileReader* QHistReader = new histFileReader();
    QHistReader->setFileName("QHistRCR.dat");
    QHistReader->setNumColumns(3);
    QHistReader->readAndSplitMultiSurfaceRestartFile();

    double finalPHistRCRValue = QHistReader->getReadFileData(1,5);//((QHistReader->m_dataReadFromFile).at(5))[1];
    EXPECT_NEAR(finalPHistRCRValue,-97144.4580441544,10e-5);

    finalPHistRCRValue = QHistReader->getReadFileData(2,5);//((QHistReader->m_dataReadFromFile).at(5))[2];
    EXPECT_NEAR(finalPHistRCRValue,-40084.8691875509,10e-5);

    delete QHistReader;
  }
}
