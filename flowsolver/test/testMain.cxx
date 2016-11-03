#include "testMain.hxx"
#include "debuggingToolsForCpp.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMain() { return 0; }

TEST_F(testMain, checkRCRSimpleShortSimulation) {
  std::string simDir = "mainTests/basic";
  setSimDirectory(simDir);
  clearOutOldFiles();
  runSimulation();

  HistFileReader* QHistReader = new HistFileReader();
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
  HistFileReader* PHistReader = new HistFileReader();
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
  HistFileReader FlowHistReader = HistFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(2);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double finalFlowHistValue = FlowHistReader.getReadFileData(1,5);
  EXPECT_NEAR(finalFlowHistValue,3.628656482154824E-003,1e-10);

  // Check a value in PressHist.dat
  HistFileReader PressHistReader = HistFileReader();
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
  HistFileReader PressHistReader = HistFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(2);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double readResult = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(10645.5858581080,readResult,1e-7);

  // Check FlowHist.dat
  HistFileReader FlowHistReader = HistFileReader();
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
  HistFileReader PressHistReader = HistFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(2);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double readResult = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(10645.5858581080,readResult,1e-7);

  // Check FlowHist.dat
  HistFileReader FlowHistReader = HistFileReader();
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
//   // HistFileReader PressHistReader = HistFileReader();
//   // PressHistReader.setFileName("PressHist.dat");
//   // PressHistReader.setNumColumns(2);
//   // PressHistReader.readFileInternalMetadata();
//   // PressHistReader.readAndSplitMultiSurfaceRestartFile();
//   // // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
//   // double readResult = PressHistReader.getReadFileData(1,5);
//   // EXPECT_NEAR(10645.5858581080,readResult,1e-8);

//   // // Check FlowHist.dat
//   // HistFileReader FlowHistReader = HistFileReader();
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
  HistFileReader PressHistReader = HistFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(2);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double pressHistResult = PressHistReader.getReadFileData(1,5);
  EXPECT_NEAR(106420.076723820,pressHistResult,1e-8);

  // Check FlowHist.dat
  HistFileReader FlowHistReader = HistFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(2);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 2nd column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double flowHistResult = FlowHistReader.getReadFileData(1,5);
  EXPECT_NEAR(1.077519395498929e-2,flowHistResult,1e-7);
}

#ifndef DISABLE_ACUSIM_TESTS
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
  HistFileReader PressHistReader = HistFileReader();
  PressHistReader.setFileName("PressHist.dat");
  PressHistReader.setNumColumns(3);
  PressHistReader.readFileInternalMetadata();
  PressHistReader.readAndSplitMultiSurfaceRestartFile();
  
  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double pressHistResult = PressHistReader.getReadFileData(0,11);
  EXPECT_NEAR(632.670675940888,pressHistResult,1e-1);
  // ...second column
  pressHistResult = PressHistReader.getReadFileData(1,11);
  EXPECT_NEAR(556.262903111148,pressHistResult,1e-1);
  // ... third column
  pressHistResult = PressHistReader.getReadFileData(2,11);
  EXPECT_NEAR(544.651629587945,pressHistResult,1e-1);

  // Check FlowHist.dat
  HistFileReader FlowHistReader = HistFileReader();
  FlowHistReader.setFileName("FlowHist.dat");
  FlowHistReader.setNumColumns(3);
  FlowHistReader.readFileInternalMetadata();
  FlowHistReader.readAndSplitMultiSurfaceRestartFile();
  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
  double flowHistResult = FlowHistReader.getReadFileData(0,11);
  EXPECT_NEAR(-889.383031270501,flowHistResult,1e-1);
  // ... 2nd column:
  flowHistResult = FlowHistReader.getReadFileData(1,11);
  EXPECT_NEAR(636.889350513117,flowHistResult,1e-1);
  // ...third column:
  flowHistResult = FlowHistReader.getReadFileData(2,11);
  EXPECT_NEAR(252.486584230077,flowHistResult,1e-1);

  // valve closed case:
  flowHistResult = FlowHistReader.getReadFileData(0,5);
  EXPECT_EQ(0,flowHistResult);
}
#endif

#ifndef DISABLE_ACUSIM_TESTS
TEST_F(testMain, checkClosedLoopWithHeart) {
  setSimDirectory("mainTests/netlist/closedLoopHeart");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check PressHist.dat
  {
    HistFileReader zeroDDomainPressures = HistFileReader();
    zeroDDomainPressures.setFileName("PressHist.dat");
    zeroDDomainPressures.setNumColumns(3);
    zeroDDomainPressures.readFileInternalMetadata();
    zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();

    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = zeroDDomainPressures.getReadFileData(0,6);
    EXPECT_NEAR(553.277617338762,pressureResult,1e-5);
    // ...second column
    pressureResult = zeroDDomainPressures.getReadFileData(1,6);
    EXPECT_NEAR(548.211423767513,pressureResult,1e-6);
    // ... third column
    pressureResult = zeroDDomainPressures.getReadFileData(2,6);
    EXPECT_NEAR(541.282726604241,pressureResult,1e-6);
  }

  // Check FlowHist.dat
  {
    HistFileReader zeroDDomainFlows = HistFileReader();
    zeroDDomainFlows.setFileName("FlowHist.dat");
    zeroDDomainFlows.setNumColumns(3);
    zeroDDomainFlows.readFileInternalMetadata();
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double flowResult = zeroDDomainFlows.getReadFileData(0,6);
    EXPECT_NEAR(-307.419802901504,flowResult,1e-8);
    // ...second column
    flowResult = zeroDDomainFlows.getReadFileData(1,6);
    EXPECT_NEAR(270.813353397073,flowResult,1e-4);
    // ... third column
    flowResult = zeroDDomainFlows.getReadFileData(2,6);
    EXPECT_NEAR(37.8070383450010,flowResult,1e-5);
  }

  // Check netlistFlows_surface_5.dat
  {
    HistFileReader zeroDDomainFlows = HistFileReader();
    zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
    zeroDDomainFlows.setNumColumns(6);
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double flowResult = zeroDDomainFlows.getReadFileData(1,5);
    EXPECT_NEAR(-287.802631760359,flowResult,1e-5);
  }

  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
  {
    HistFileReader closedLoopDownstreamPressures = HistFileReader();
    closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
    closedLoopDownstreamPressures.setNumColumns(5);
    closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

    double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,5);
    EXPECT_NEAR(533.100216837156,pressureResult,1e-8);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(2,5);
    EXPECT_NEAR(533.100216837156,pressureResult,1e-8);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(3,5);
    EXPECT_NEAR(0.0,pressureResult,1e-8);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(4,5);
    EXPECT_NEAR(535.821025669528,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
    HistFileReader closedLoopDownstreamPressures = HistFileReader();
    closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
    closedLoopDownstreamPressures.setNumColumns(4);
    closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

    double flowResult = closedLoopDownstreamPressures.getReadFileData(1,5);
    EXPECT_EQ(0.0,flowResult);

    flowResult = closedLoopDownstreamPressures.getReadFileData(2,5);
    EXPECT_NEAR(27.2080883237159,flowResult,1e-8);

    flowResult = closedLoopDownstreamPressures.getReadFileData(3,5);
    EXPECT_NEAR(-27.2080883237159,flowResult,1e-8);
  }
}
#endif

#ifndef DISABLE_ACUSIM_TESTS
TEST_F(testMain, checkMixedNetlistAndRCRT)
{
  setSimDirectory("mainTests/legacy/netlistWithRCRTs");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

   // Check PressHist.dat
  {
    HistFileReader pressHistData = HistFileReader();
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
    HistFileReader flowHistData = HistFileReader();
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
    HistFileReader zeroDDomainFlows = HistFileReader();
    zeroDDomainFlows.setFileName("netlistPressures_surface_5.dat");
    zeroDDomainFlows.setNumColumns(7);
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = zeroDDomainFlows.getReadFileData(1,5);
    EXPECT_NEAR(798.293386775323,pressureResult,1e-8);
  }

  {
    HistFileReader* PHistReader = new HistFileReader();
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
    HistFileReader* QHistReader = new HistFileReader();
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
#endif

#ifndef DISABLE_ACUSIM_TESTS
TEST_F(testMain, checkClosedLoopWithHeartRestart) {
  setSimDirectory("mainTests/netlist/closedLoopHeartRestart");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check PressHist.dat
  {
    HistFileReader zeroDDomainPressures = HistFileReader();
    zeroDDomainPressures.setFileName("PressHist.dat");
    zeroDDomainPressures.setNumColumns(3);
    zeroDDomainPressures.readFileInternalMetadata();
    zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();

    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double pressureResult = zeroDDomainPressures.getReadFileData(0,11);
    EXPECT_NEAR(632.747782002238,pressureResult,1e-1);
    // ...second column
    pressureResult = zeroDDomainPressures.getReadFileData(1,11);
    EXPECT_NEAR(558.547714500901,pressureResult,1e-1);
    // ... third column
    pressureResult = zeroDDomainPressures.getReadFileData(2,11);
    EXPECT_NEAR(547.411051345682,pressureResult,1e-1);
  }

  // Check FlowHist.dat
  {
    HistFileReader zeroDDomainFlows = HistFileReader();
    zeroDDomainFlows.setFileName("FlowHist.dat");
    zeroDDomainFlows.setNumColumns(3);
    zeroDDomainFlows.readFileInternalMetadata();
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double flowResult = zeroDDomainFlows.getReadFileData(0,11);
    EXPECT_NEAR(-851.496522781367,flowResult,1e-1);
    // ...second column
    flowResult = zeroDDomainFlows.getReadFileData(1,11);
    EXPECT_NEAR(610.138952390568,flowResult,1e-1);
    // ... third column
    flowResult = zeroDDomainFlows.getReadFileData(2,11);
    EXPECT_NEAR(241.349979602039,flowResult,1e-1);
  }

  // Check netlistFlows_surface_5.dat
  {
    HistFileReader zeroDDomainFlows = HistFileReader();
    zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
    zeroDDomainFlows.setNumColumns(6);
    zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
    // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
    double flowResult = zeroDDomainFlows.getReadFileData(1,11);
    EXPECT_NEAR(-851.496522781367,flowResult,1e-1);
  }

  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
  {
    HistFileReader closedLoopDownstreamPressures = HistFileReader();
    closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
    closedLoopDownstreamPressures.setNumColumns(5);
    closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

    double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,11);
    EXPECT_NEAR(533.100554393682,pressureResult,1e-1);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(2,11);
    EXPECT_NEAR(533.100554393682,pressureResult,1e-1);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(3,11);
    EXPECT_EQ(0.0,pressureResult);

    pressureResult = closedLoopDownstreamPressures.getReadFileData(4,11);
    EXPECT_NEAR(535.928307373969,pressureResult,1e-1);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
    HistFileReader closedLoopDownstreamPressures = HistFileReader();
    closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
    closedLoopDownstreamPressures.setNumColumns(4);
    closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

    double flowResult = closedLoopDownstreamPressures.getReadFileData(1,11);
    EXPECT_EQ(0.0, flowResult);

    flowResult = closedLoopDownstreamPressures.getReadFileData(2,11);
    EXPECT_NEAR(28.2775298028741,flowResult,1e-1);

    flowResult = closedLoopDownstreamPressures.getReadFileData(3,11);
    EXPECT_NEAR(-28.2775298028741,flowResult,1e-1);
  }
}
#endif