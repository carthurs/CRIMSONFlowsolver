#include "testMainWithZeroDDomain.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMainWithZeroDDomain() { return 0; }

TEST_F(testMainWithZeroDDomain, checkWithNetlistRCRs) {
  setSimDirectory("mainTests/zeroDDomain/netlistRCRs");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_-1.dat.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(9);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,5);
	  EXPECT_NEAR(10000.7824501887,pressureResult,1e-8);
	  // ...second column
	  pressureResult = zeroDDomainPressures.getReadFileData(2,5);
	  EXPECT_NEAR(1349.01093286719,pressureResult,1e-8);
	  // ... third column
	  pressureResult = zeroDDomainPressures.getReadFileData(3,5);
	  EXPECT_NEAR(1349.01093286719,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(4,5);
	  EXPECT_NEAR(-4.09407660591749e-12,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(5,5);
	  EXPECT_NEAR(9963.38206912058,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(6,5);
	  EXPECT_NEAR(9963.41939485764,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(7,5);
	  EXPECT_NEAR(9954.78489062223,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(8,5);
	  EXPECT_NEAR(9954.78489062223,pressureResult,1e-8);
  }

  // Check netlistFlows_surface_-1.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double flowResult = zeroDDomainFlows.getReadFileData(1,5);
	  EXPECT_NEAR(-37.3257370576594,flowResult,1e-8);
	  // ... 2nd column:
	  flowResult = zeroDDomainFlows.getReadFileData(2,5);
	  EXPECT_NEAR(8597.17849834857,flowResult,1e-8);
	  // ...third column:
	  flowResult = zeroDDomainFlows.getReadFileData(3,5);
	  EXPECT_NEAR(8597.17849834857,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(4,5);
	  EXPECT_NEAR(-37.3257370576612,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(5,5);
	  EXPECT_NEAR(8597.17849834857,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(6,5);
	  EXPECT_NEAR(8597.17849834857,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(7,5);
	  EXPECT_NEAR(-17157.0312596395,flowResult,1e-8);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
		zeroDDomainFlows.setNumColumns(4);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = zeroDDomainFlows.getReadFileData(1,5);
		EXPECT_NEAR(130.091368993213,flowResult,1e-8);
  }

  // Check netlistFlows_surface_6.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_6.dat");
		zeroDDomainFlows.setNumColumns(4);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = zeroDDomainFlows.getReadFileData(1,5);
		EXPECT_NEAR(8772.95873611765,flowResult,1e-8);
  }
  
  // Check netlistFlows_surface_7.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_7.dat");
		zeroDDomainFlows.setNumColumns(4);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = zeroDDomainFlows.getReadFileData(1,5);
		EXPECT_NEAR(8772.95873611765,flowResult,1e-8);
  }
}

TEST_F(testMainWithZeroDDomain, checkClosedLoopWithHeart) {
  setSimDirectory("mainTests/zeroDDomain/closedLoopHeart");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(9);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,1000);
	  EXPECT_NEAR(19268.692618963,pressureResult,2e-8);
	  // ...second column
	  pressureResult = zeroDDomainPressures.getReadFileData(2,1000);
	  EXPECT_NEAR(12714.0753689733,pressureResult,1e-8);
	  // ... third column
	  pressureResult = zeroDDomainPressures.getReadFileData(3,1000);
	  EXPECT_NEAR(12714.0753689733,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(4,1000);
	  EXPECT_NEAR(4.54590964207829e-10,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(5,1000);
	  EXPECT_NEAR(16437.5433223572,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(6,1000);
	  EXPECT_NEAR(16691.0938886357,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(7,1000);
	  EXPECT_NEAR(16318.102383974,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(8,1000);
	  EXPECT_NEAR(16318.102383974,pressureResult,1e-8);
  }

  // Check netlistFlows_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double flowResult = zeroDDomainFlows.getReadFileData(1,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,5e-7);
	  // ... 2nd column:
	  flowResult = zeroDDomainFlows.getReadFileData(2,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,3e-7);
	  // ...third column:
	  flowResult = zeroDDomainFlows.getReadFileData(3,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(4,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,5e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(5,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(6,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(7,1000);
	  EXPECT_NEAR(14668.6895121776,flowResult,2e-7);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
		zeroDDomainFlows.setNumColumns(6);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double pressureResult = zeroDDomainFlows.getReadFileData(1,1000);
		EXPECT_NEAR(-245488.583296014,pressureResult,5e-7);
  }

  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
  	closedLoopDownstreamPressures.setNumColumns(5);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(528.30042814354,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(528.30042814354,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(0.0,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,1000);
  	EXPECT_NEAR(5410.64387151683,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
  	closedLoopDownstreamPressures.setNumColumns(4);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double flowResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(9.08301652082719e-305,flowResult,1e-8);

  	flowResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(48823.4344337329,flowResult,1e-8);

  	flowResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(-48823.4344337329,flowResult,1e-8);
  }
}

TEST_F(testMainWithZeroDDomain, checkClosedLoopWithHeartAndEndNodeOfComponentConnectedTo3D) {
  setSimDirectory("mainTests/zeroDDomain/closedLoopHeartReversedInterfaceComponent");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(9);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,1000);
	  EXPECT_NEAR(19268.692618963,pressureResult,2e-8);
	  // ...second column
	  pressureResult = zeroDDomainPressures.getReadFileData(2,1000);
	  EXPECT_NEAR(12714.0753689733,pressureResult,1e-8);
	  // ... third column
	  pressureResult = zeroDDomainPressures.getReadFileData(3,1000);
	  EXPECT_NEAR(12714.0753689733,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(4,1000);
	  EXPECT_NEAR(4.54590964207829e-10,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(5,1000);
	  EXPECT_NEAR(16437.5433223572,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(6,1000);
	  EXPECT_NEAR(16691.0938886357,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(7,1000);
	  EXPECT_NEAR(16318.102383974,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(8,1000);
	  EXPECT_NEAR(16318.102383974,pressureResult,1e-8);
  }

  // Check netlistFlows_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double flowResult = zeroDDomainFlows.getReadFileData(1,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,5e-7);
	  // ... 2nd column:
	  flowResult = zeroDDomainFlows.getReadFileData(2,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,3e-7);
	  // ...third column:
	  flowResult = zeroDDomainFlows.getReadFileData(3,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(4,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,5e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(5,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(6,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(7,1000);
	  EXPECT_NEAR(14668.6895121776,flowResult,2e-7);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
		zeroDDomainFlows.setNumColumns(6);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = zeroDDomainFlows.getReadFileData(1,1000);
		EXPECT_NEAR(245488.583296014,flowResult,5e-7);
  }

  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
  	closedLoopDownstreamPressures.setNumColumns(5);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(528.30042814354,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(528.30042814354,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(0.0,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,1000);
  	EXPECT_NEAR(5410.64387151683,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
  	closedLoopDownstreamPressures.setNumColumns(4);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double flowResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(9.08301652082719e-305,flowResult,1e-8);

  	flowResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(48823.4344337329,flowResult,1e-8);

  	flowResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(-48823.4344337329,flowResult,1e-8);
  }
}

TEST_F(testMainWithZeroDDomain, checkPythonElastanceController) {
  setSimDirectory("mainTests/zeroDDomain/pythonParameterControllers");
  clearOutOldFiles();

  try {
  	runSimulation();
  } catch (const std::exception& e) {
      std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
      throw;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(9);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,1000);
	  EXPECT_NEAR(19268.692618963,pressureResult,2e-8);
	  // ...second column
	  pressureResult = zeroDDomainPressures.getReadFileData(2,1000);
	  EXPECT_NEAR(12714.0753689733,pressureResult,1e-8);
	  // ... third column
	  pressureResult = zeroDDomainPressures.getReadFileData(3,1000);
	  EXPECT_NEAR(12714.0753689733,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(4,1000);
	  EXPECT_NEAR(4.54590964207829e-10,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(5,1000);
	  EXPECT_NEAR(16437.5433223572,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(6,1000);
	  EXPECT_NEAR(16691.0938886357,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(7,1000);
	  EXPECT_NEAR(16318.102383974,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(8,1000);
	  EXPECT_NEAR(16318.102383974,pressureResult,1e-8);
  }

  // Check netlistFlows_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double flowResult = zeroDDomainFlows.getReadFileData(1,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,5e-7);
	  // ... 2nd column:
	  flowResult = zeroDDomainFlows.getReadFileData(2,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,3e-7);
	  // ...third column:
	  flowResult = zeroDDomainFlows.getReadFileData(3,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(4,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,5e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(5,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(6,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,3e-7);

	  flowResult = zeroDDomainFlows.getReadFileData(7,1000);
	  EXPECT_NEAR(14668.6895121776,flowResult,2e-7);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
		zeroDDomainFlows.setNumColumns(6);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = zeroDDomainFlows.getReadFileData(1,1000);
		EXPECT_NEAR(-245488.583296014,flowResult,5e-7);
  }

  // Check netlistPressures_surface_5.dat - this is really a very minimal
  // test of the nodal python parameter controller - it's set up to prescribe
  // cosine-in-time pressure on a node which does not affect the results
  // (the base of the VolumeTrackingComponent - whose pressure is not relevant
  // to any calculations). We just read the pressure that it has on the last
  // timestep and confirm that it is what would be expected after 1000 steps
  // of cosine prescription.
  {
		histFileReader heartModelPressures = histFileReader();
		heartModelPressures.setFileName("netlistPressures_surface_5.dat");
		heartModelPressures.setNumColumns(7);
		heartModelPressures.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double pressureResult = heartModelPressures.getReadFileData(6,1000);
		EXPECT_NEAR(0.990355044196067,pressureResult,1e-8);
  }

  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
  	closedLoopDownstreamPressures.setNumColumns(5);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(528.30042814354,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(528.30042814354,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(0.0,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,1000);
  	EXPECT_NEAR(5410.64387151683,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
  	closedLoopDownstreamPressures.setNumColumns(4);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double flowResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(9.08301652082719e-305,flowResult,1e-8);

  	flowResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(48823.4344337329,flowResult,1e-8);

  	flowResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(-48823.4344337329,flowResult,1e-8);
  }
}

// This test includes Python controllers on both pressure nodes
// and resistances in the downstream closed loop.
TEST_F(testMainWithZeroDDomain, checkDownstreamPythonElastanceController) {
	  setSimDirectory("mainTests/zeroDDomain/downstreamPythonParameterControllers");
	  clearOutOldFiles();

	  runSimulation();
	  MPI_Barrier(MPI_COMM_WORLD);

	  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
	  {
	  	histFileReader closedLoopDownstreamPressures = histFileReader();
	  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
	  	closedLoopDownstreamPressures.setNumColumns(5);
	  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

	  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
	  	EXPECT_NEAR(498.519721797159,pressureResult,1e-8);

	  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
	  	EXPECT_NEAR(498.519721797159,pressureResult,1e-8);

	  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
	  	EXPECT_NEAR(0.990355044196067,pressureResult,1e-8);

	  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,1000);
	  	EXPECT_NEAR(5489.82349624579,pressureResult,1e-8);
	  }

	  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
	  {
	  	histFileReader closedLoopDownstreamPressures = histFileReader();
	  	closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
	  	closedLoopDownstreamPressures.setNumColumns(4);
	  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

	  	double flowResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
	  	EXPECT_NEAR(0.0,flowResult,1e-8);

	  	flowResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
	  	EXPECT_NEAR(49913.0377444863,flowResult,1e-8);

	  	flowResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
	  	EXPECT_NEAR(-49913.0377444863,flowResult,1e-8);
	  }
  }

// This test includes a Python prescribed flow controller on one of the Windkessels.
// There's no closed loop.
TEST_F(testMainWithZeroDDomain, checkPythonFlowDatFilePrescriber) 
{
  setSimDirectory("mainTests/zeroDDomain/pythonFlowPrescriptions");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_7.dat (the Windkessel that has a flow controller on its distal resistor)
  {
  	histFileReader flowControlledWindkesselPressures = histFileReader();
  	flowControlledWindkesselPressures.setFileName("netlistPressures_surface_7.dat");
  	flowControlledWindkesselPressures.setNumColumns(5);
  	flowControlledWindkesselPressures.readAndSplitMultiSurfaceRestartFile();

  	double pressureResult = flowControlledWindkesselPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(26060.1257730391,pressureResult,1e-8);

  	pressureResult = flowControlledWindkesselPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(15738.7041709838,pressureResult,1e-8);

  	pressureResult = flowControlledWindkesselPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(15441.7041709838,pressureResult,1e-8);

  	pressureResult = flowControlledWindkesselPressures.getReadFileData(4,1000);
  	EXPECT_EQ(0.0,pressureResult);
  }

  // Check netlistFlows_surface_5.dat (the heart model circuit)
  {
  	histFileReader heartModelFlows = histFileReader();
  	heartModelFlows.setFileName("netlistFlows_surface_5.dat");
  	heartModelFlows.setNumColumns(6);
  	heartModelFlows.readAndSplitMultiSurfaceRestartFile();

  	double flowResult = heartModelFlows.getReadFileData(1,1000);
  	EXPECT_NEAR(-749259.15077706,flowResult,3e-8);

  	flowResult = heartModelFlows.getReadFileData(2,1000);
  	EXPECT_NEAR(749259.15077706,flowResult,3e-8);

  	flowResult = heartModelFlows.getReadFileData(3,1000);
  	EXPECT_EQ(0.0,flowResult);

  	flowResult = heartModelFlows.getReadFileData(4,1000);
  	EXPECT_EQ(0.0,flowResult);

  	flowResult = heartModelFlows.getReadFileData(5,1000);
  	EXPECT_NEAR(-749259.15077706,flowResult,3e-8);
  }
}

// This test includes a Python prescribed pressure controller on one of the Windkessels.
// There's no closed loop.
TEST_F(testMainWithZeroDDomain, checkPythonPressureDatFilePrescriber) 
{
  setSimDirectory("mainTests/zeroDDomain/pythonPressurePrescriptions");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_5.dat (the heart model circuit)
  {
  	histFileReader heartPressureFile = histFileReader();
  	heartPressureFile.setFileName("netlistPressures_surface_5.dat");
  	heartPressureFile.setNumColumns(7);
  	heartPressureFile.readAndSplitMultiSurfaceRestartFile();

  	double heartModelPressure = heartPressureFile.getReadFileData(1,1000);
  	EXPECT_NEAR(26463.9215456089,heartModelPressure,1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(2,1000);
  	EXPECT_NEAR(26526.8020490521,heartModelPressure,1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(3,1000);
  	EXPECT_NEAR(26535.1245604121,heartModelPressure,1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(4,1000);
  	EXPECT_NEAR(26535.1245604121,heartModelPressure, 1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(5,1000);
  	EXPECT_EQ(533.2,heartModelPressure);

  	heartModelPressure = heartPressureFile.getReadFileData(6,1000);
  	EXPECT_NEAR(0.990355044196067,heartModelPressure, 1e-8);
  }

  // Check netlistFlows_surface_7.dat (the Windkessel that has a Python pressures controller on its downstream pressure node)
  {
  	histFileReader windkesselFlows = histFileReader();
  	windkesselFlows.setFileName("netlistFlows_surface_7.dat");
  	windkesselFlows.setNumColumns(4);
  	windkesselFlows.readAndSplitMultiSurfaceRestartFile();

  	double flowInWindkesselThatHasPythonPrescribedPressure = windkesselFlows.getReadFileData(1,1000);
  	EXPECT_NEAR(407584.54276605,flowInWindkesselThatHasPythonPrescribedPressure,2e-8);

  	flowInWindkesselThatHasPythonPrescribedPressure = windkesselFlows.getReadFileData(2,1000);
  	EXPECT_NEAR(39466.4517277277,flowInWindkesselThatHasPythonPrescribedPressure,1e-8);

  	flowInWindkesselThatHasPythonPrescribedPressure = windkesselFlows.getReadFileData(3,1000);
  	EXPECT_NEAR(368118.091038322,flowInWindkesselThatHasPythonPrescribedPressure,2e-8);

  }
}

// This test has a masterController.py, which is not associated
// to any Netlist surface, but it sends a master control signal to
// all other surfaces, which use it to adjust their control output.
//
// We have a closed loop (downstream) circuit, and there is a 
// controller in there too, recieving the masterController broadcast.
TEST_F(testMainWithZeroDDomain, checkPythonControlBroadcastsAndMasterController) 
{
  setSimDirectory("mainTests/zeroDDomain/pythonBroadcastCommunications");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_5.dat (the heart model circuit)
  {
  	histFileReader heartPressureFile = histFileReader();
  	heartPressureFile.setFileName("netlistPressures_surface_5.dat");
  	heartPressureFile.setNumColumns(7);
  	heartPressureFile.readAndSplitMultiSurfaceRestartFile();

  	double heartModelPressure = heartPressureFile.getReadFileData(1,1000);
  	EXPECT_NEAR(20044.9904142384,heartModelPressure,1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(2,1000);
  	EXPECT_NEAR(20114.4151803054,heartModelPressure,1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(3,1000);
  	EXPECT_NEAR(20118.0514049299,heartModelPressure,1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(4,1000);
  	EXPECT_NEAR(20118.0514049299,heartModelPressure, 1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(5,1000);
  	EXPECT_NEAR(540.802034459378,heartModelPressure, 1e-8);

  	heartModelPressure = heartPressureFile.getReadFileData(6,1000);
  	EXPECT_NEAR(0.990355044196067,heartModelPressure, 1e-8);
  }

  // Check netlistFlows_surface_7.dat (the Windkessel that has a Python pressures controller on its downstream pressure node)
  {
  	histFileReader windkesselFlows = histFileReader();
  	windkesselFlows.setFileName("netlistFlows_surface_7.dat");
  	windkesselFlows.setNumColumns(4);
  	windkesselFlows.readAndSplitMultiSurfaceRestartFile();

  	double flowInWindkesselThatHasPythonPrescribedPressure = windkesselFlows.getReadFileData(1,1000);
  	EXPECT_NEAR(174767.004330894,flowInWindkesselThatHasPythonPrescribedPressure,1e-7);

  	flowInWindkesselThatHasPythonPrescribedPressure = windkesselFlows.getReadFileData(2,1000);
  	EXPECT_NEAR(27077.5494999395,flowInWindkesselThatHasPythonPrescribedPressure,1e-8);

  	flowInWindkesselThatHasPythonPrescribedPressure = windkesselFlows.getReadFileData(3,1000);
  	EXPECT_NEAR(147689.454830954,flowInWindkesselThatHasPythonPrescribedPressure,1e-7);
  }

  // Check netlistPressures_zeroDDomainReplacement.dat
  {
  	histFileReader zeroDDomainReplacementPressureFile = histFileReader();
  	zeroDDomainReplacementPressureFile.setFileName("netlistPressures_zeroDDomainReplacement.dat");
  	zeroDDomainReplacementPressureFile.setNumColumns(9);
  	zeroDDomainReplacementPressureFile.readAndSplitMultiSurfaceRestartFile();

  	double zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(1,1000);
  	EXPECT_NEAR(23769.3996600964,zeroDDomainNodalPressure,1e-8);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(2,1000);
  	EXPECT_NEAR(14049.1436358498,zeroDDomainNodalPressure,1e-8);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(3,1000);
  	EXPECT_NEAR(14049.1436358498,zeroDDomainNodalPressure,1e-8);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(4,1000);
  	EXPECT_EQ(0.0,zeroDDomainNodalPressure);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(5,1000);
  	EXPECT_NEAR(19606.4633221664,zeroDDomainNodalPressure, 1e-8);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(6,1000);
  	EXPECT_NEAR(19976.7484056073,zeroDDomainNodalPressure, 1e-8);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(7,1000);
  	EXPECT_NEAR(19428.1963177128,zeroDDomainNodalPressure, 1e-8);

  	zeroDDomainNodalPressure = zeroDDomainReplacementPressureFile.getReadFileData(8,1000);
  	EXPECT_NEAR(19428.1963177128,zeroDDomainNodalPressure, 1e-8);
  }

   // Check netlistPressures_downstreamCircuit_0.dat
  {
  	histFileReader downstreamCircuitNodalPressureFile = histFileReader();
  	downstreamCircuitNodalPressureFile.setFileName("netlistPressures_downstreamCircuit_0.dat");
  	downstreamCircuitNodalPressureFile.setNumColumns(5);
  	downstreamCircuitNodalPressureFile.readAndSplitMultiSurfaceRestartFile();

  	double downstreamCircuitNodalPressure = downstreamCircuitNodalPressureFile.getReadFileData(1,1000);
  	EXPECT_NEAR(540.802034459378,downstreamCircuitNodalPressure, 1e-8);

  	downstreamCircuitNodalPressure = downstreamCircuitNodalPressureFile.getReadFileData(2,1000);
  	EXPECT_NEAR(540.802034459378,downstreamCircuitNodalPressure, 1e-8);

  	downstreamCircuitNodalPressure = downstreamCircuitNodalPressureFile.getReadFileData(3,1000);
  	EXPECT_NEAR(0.990355044196067,downstreamCircuitNodalPressure, 1e-8);

  	downstreamCircuitNodalPressure = downstreamCircuitNodalPressureFile.getReadFileData(4,1000);
  	EXPECT_NEAR(5956.31193444729,downstreamCircuitNodalPressure, 1e-8);
  }
}

TEST_F(testMainWithZeroDDomain, checkTwoDomainsWithClosedLoop) {
  setSimDirectory("mainTests/zeroDDomain/twoDomainsClosedLoop");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(13);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,1000);
	  EXPECT_NEAR(18997.5594514993,pressureResult,1e-8);
	  // ...second column
	  pressureResult = zeroDDomainPressures.getReadFileData(2,1000);
	  EXPECT_NEAR(14052.5402246331,pressureResult,1e-8);
	  // ... third column
	  pressureResult = zeroDDomainPressures.getReadFileData(3,1000);
	  EXPECT_NEAR(19122.6055974038,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(4,1000);
	  EXPECT_NEAR(14050.6840783447,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(5,1000);
	  EXPECT_NEAR(0.0,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(6,1000);
	  EXPECT_NEAR(17548.0750926851,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(7,1000);
	  EXPECT_NEAR(0.0,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(8,1000);
	  EXPECT_NEAR(17641.193233404,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(9,1000);
	  EXPECT_NEAR(17677.8872436832,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(10,1000);
	  EXPECT_NEAR(17435.9457814192,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(11,1000);
	  EXPECT_NEAR(17773.9717975567,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(12,1000);
	  EXPECT_NEAR(17526.1528519115,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamFlows = histFileReader();
  	closedLoopDownstreamFlows.setFileName("netlistFlows_downstreamCircuit_0.dat");
  	closedLoopDownstreamFlows.setNumColumns(4);
  	closedLoopDownstreamFlows.readAndSplitMultiSurfaceRestartFile();

  	double flowResult = closedLoopDownstreamFlows.getReadFileData(1,1000);
  	EXPECT_NEAR(0.0,flowResult,1e-8);

  	flowResult = closedLoopDownstreamFlows.getReadFileData(2,1000);
  	EXPECT_NEAR(33867.6988549767,flowResult,1e-8);

  	flowResult = closedLoopDownstreamFlows.getReadFileData(3,1000);
  	EXPECT_NEAR(-33867.6988549767,flowResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_1.dat");
  	closedLoopDownstreamPressures.setNumColumns(5);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,1000);
  	EXPECT_NEAR(528.058272620996,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,1000);
  	EXPECT_NEAR(528.058272620996,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,1000);
  	EXPECT_NEAR(0.0,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,1000);
  	EXPECT_NEAR(3913.83581629739,pressureResult,1e-8);
  }
}

TEST_F(testMainWithZeroDDomain, checkPythonControlledCoronary) {
	setSimDirectory("mainTests/netlist/coronaryControlPython");
	clearOutOldFiles();
	runSimulation();
	MPI_Barrier(MPI_COMM_WORLD);

	{
		histFileReader mvo2HistoryReader = histFileReader();
		mvo2HistoryReader.setFileName("MVO2History.dat");
		mvo2HistoryReader.setNumColumns(1);
		mvo2HistoryReader.readAndSplitMultiSurfaceRestartFile();

		double mvo2Result = mvo2HistoryReader.getReadFileData(0,13);
		EXPECT_NEAR(4.525274383785653498e+01, mvo2Result, 1e-8);
	}

	{
		histFileReader r_d_historyReader = histFileReader();
		r_d_historyReader.setFileName("R_d_history.dat");
		r_d_historyReader.setNumColumns(1);
		r_d_historyReader.readAndSplitMultiSurfaceRestartFile();

		double r_d_result = r_d_historyReader.getReadFileData(0,6640);
		EXPECT_NEAR(4.350306627092590617e+01, r_d_result, 1e-8);
	}
	
	{
		histFileReader r_p_historyReader = histFileReader();
		r_p_historyReader.setFileName("R_p_history.dat");
		r_p_historyReader.setNumColumns(1);
		r_p_historyReader.readAndSplitMultiSurfaceRestartFile();

		double r_p_result = r_p_historyReader.getReadFileData(0,6640);
		EXPECT_NEAR(6.475766148479123352, r_p_result, 1e-8);
	}

	{
		histFileReader flowHistReader = histFileReader();
		flowHistReader.setFileName("netlistFlows_surface_6.dat");
		flowHistReader.setNumColumns(6);
		flowHistReader.readAndSplitMultiSurfaceRestartFile();

		double lastTimestepCoronaryFlow = flowHistReader.getReadFileData(1,7000);
		EXPECT_NEAR(879.323421682284, lastTimestepCoronaryFlow, 1e-8);
	}	
}

TEST_F(testMainWithZeroDDomain, checkRestartWithNetlistRCRs) {
  setSimDirectory("mainTests/zeroDDomain/restartNetlistRCRs");
  clearOutOldFiles();

  try {
  	runSimulation();
  } catch (const std::exception& e) {
      std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
      throw;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_-1.dat.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(9);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,10);
	  EXPECT_NEAR(9999.80845698841,pressureResult,1e-8);
  }

  // Check netlistFlows_surface_-1.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  // ...third column:
	  double flowResult = zeroDDomainFlows.getReadFileData(3,10);
	  EXPECT_NEAR(7791.68348266783,flowResult,1e-8);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader netlist1Flow = histFileReader();
		netlist1Flow.setFileName("netlistFlows_surface_5.dat");
		netlist1Flow.setNumColumns(4);
		netlist1Flow.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = netlist1Flow.getReadFileData(1,10);
		EXPECT_NEAR(-657.680077994277,flowResult,1e-8);
  }

  // Check netlistFlows_surface_6.dat
  {
		histFileReader netlist2Flow = histFileReader();
		netlist2Flow.setFileName("netlistFlows_surface_6.dat");
		netlist2Flow.setNumColumns(4);
		netlist2Flow.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = netlist2Flow.getReadFileData(1,10);
		EXPECT_NEAR(7943.47173011129,flowResult,1e-8);
  }
  
  // Check netlistFlows_surface_7.dat
  {
		histFileReader netlist3Flow = histFileReader();
		netlist3Flow.setFileName("netlistFlows_surface_7.dat");
		netlist3Flow.setNumColumns(4);
		netlist3Flow.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = netlist3Flow.getReadFileData(1,10);
		EXPECT_NEAR(7943.47173011129,flowResult,1e-8);
  }

  // Check netlistPressures_surface_5.dat
  {
		histFileReader netlist1Pressure = histFileReader();
		netlist1Pressure.setFileName("netlistPressures_surface_5.dat");
		netlist1Pressure.setNumColumns(5);
		netlist1Pressure.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = netlist1Pressure.getReadFileData(1,10);
		EXPECT_NEAR(9341.47083042585,flowResult,1e-8);
  }
}

TEST_F(testMainWithZeroDDomain, checkRestartTwoDomainsWithClosedLoop) {
  setSimDirectory("mainTests/zeroDDomain/restartTwoDomainsClosedLoop");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_zeroDDomainReplacement.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_zeroDDomainReplacement.dat");
	  zeroDDomainPressures.setNumColumns(13);
	  zeroDDomainPressures.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double pressureResult = zeroDDomainPressures.getReadFileData(1,10);
	  EXPECT_NEAR(10637.4012376802,pressureResult,1e-8);
	  // ...second column
	  pressureResult = zeroDDomainPressures.getReadFileData(2,10);
	  EXPECT_NEAR(10624.7099151321,pressureResult,1e-8);
	  // ... third column
	  pressureResult = zeroDDomainPressures.getReadFileData(3,10);
	  EXPECT_NEAR(10631.8631238922,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(4,10);
	  EXPECT_NEAR(10616.5271591846,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(5,10);
	  EXPECT_NEAR(0.0,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(6,10);
	  EXPECT_NEAR(10637.4012376802,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(7,10);
	  EXPECT_NEAR(0.0,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(8,10);
	  EXPECT_NEAR(10631.8631238922,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(9,10);
	  EXPECT_NEAR(10637.4012376802,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(10,10);
	  EXPECT_NEAR(10636.9941270919,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(11,10);
	  EXPECT_NEAR(10631.8631238922,pressureResult,1e-8);

	  pressureResult = zeroDDomainPressures.getReadFileData(12,10);
	  EXPECT_NEAR(10631.3717575833,pressureResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamFlows = histFileReader();
  	closedLoopDownstreamFlows.setFileName("netlistFlows_downstreamCircuit_0.dat");
  	closedLoopDownstreamFlows.setNumColumns(4);
  	closedLoopDownstreamFlows.readAndSplitMultiSurfaceRestartFile();

  	double flowResult = closedLoopDownstreamFlows.getReadFileData(1,10);
  	EXPECT_NEAR(0.0,flowResult,1e-8);

  	flowResult = closedLoopDownstreamFlows.getReadFileData(2,10);
  	EXPECT_NEAR(25207.683008632,flowResult,1e-8);

  	flowResult = closedLoopDownstreamFlows.getReadFileData(3,10);
  	EXPECT_NEAR(-25207.683008632,flowResult,1e-8);
  }

  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
  {
  	histFileReader closedLoopDownstreamPressures = histFileReader();
  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_1.dat");
  	closedLoopDownstreamPressures.setNumColumns(5);
  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,10);
  	EXPECT_NEAR(533.554895556009,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,10);
  	EXPECT_NEAR(533.554895556009,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,10);
  	EXPECT_NEAR(0.0,pressureResult,1e-8);

  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,10);
  	EXPECT_NEAR(3056.36087755523,pressureResult,1e-8);
  }
}

// Note that this test is imperfect for 2 reasons:
// 1) Nobody checks the pickling, this is just a test of unpickling
// 2) The results are checked against an initial run of this test,
//    which gave values very close (within 1mmHg) to the non-restarted
//    version of this test. So there is a slight discrepancy in the
//    restart, and this should one day be investigated.
TEST_F(testMainWithZeroDDomain, checkPythonUnpicklers) {
	  setSimDirectory("mainTests/zeroDDomain/pythonPicklers");
	  clearOutOldFiles();

	  runSimulation();
	  MPI_Barrier(MPI_COMM_WORLD);

	  // Check netlistPressures_downstreamCircuit_0.dat (the loop-closing circuit)
	  {
	  	histFileReader closedLoopDownstreamPressures = histFileReader();
	  	closedLoopDownstreamPressures.setFileName("netlistPressures_downstreamCircuit_0.dat");
	  	closedLoopDownstreamPressures.setNumColumns(5);
	  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();
                       
	  	double pressureResult = closedLoopDownstreamPressures.getReadFileData(1,201);
	  	EXPECT_NEAR(551.944803344236,pressureResult,1e-8);

	  	pressureResult = closedLoopDownstreamPressures.getReadFileData(2,201);
	  	EXPECT_NEAR(551.944803344236,pressureResult,1e-8);

	  	pressureResult = closedLoopDownstreamPressures.getReadFileData(3,201);
	  	EXPECT_NEAR(0.980066577841242,pressureResult,1e-8);

	  	pressureResult = closedLoopDownstreamPressures.getReadFileData(4,201);
	  	EXPECT_NEAR(5779.50484489326,pressureResult,1e-8);
	  }

	  // Check netlistFlows_downstreamCircuit_0.dat (the loop-closing circuit)
	  {
	  	histFileReader closedLoopDownstreamPressures = histFileReader();
	  	closedLoopDownstreamPressures.setFileName("netlistFlows_downstreamCircuit_0.dat");
	  	closedLoopDownstreamPressures.setNumColumns(4);
	  	closedLoopDownstreamPressures.readAndSplitMultiSurfaceRestartFile();

	  	double flowResult = closedLoopDownstreamPressures.getReadFileData(1,201);
	  	EXPECT_NEAR(0.0,flowResult,1e-8);

	  	flowResult = closedLoopDownstreamPressures.getReadFileData(2,201);
	  	EXPECT_NEAR(52275.6004154903,flowResult,1e-8);

	  	flowResult = closedLoopDownstreamPressures.getReadFileData(3,201);
	  	EXPECT_NEAR(-52275.6004154903,flowResult,1e-8);
	  }
  }


TEST_F(testMainWithZeroDDomain, checkBctFlowPrescriber) {
  setSimDirectory("mainTests/zeroDDomain/bctFlowPrescriber");
  clearOutOldFiles();

  try {
  	runSimulation();
  } catch (const std::exception& e) {
      std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
      throw;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_-1.dat.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  double flowResult = zeroDDomainFlows.getReadFileData(5,1501);
	  EXPECT_NEAR(259.543214447641,flowResult,1e-8);
	  flowResult = zeroDDomainFlows.getReadFileData(6,1501);
	  EXPECT_NEAR(84.6680132335589,flowResult,1e-8);
	  flowResult = zeroDDomainFlows.getReadFileData(7,1501);
	  EXPECT_NEAR(-24.7370038486727,flowResult,1e-8);
  }

  // Check netlistFlows_surface_-1.dat
  {
	  histFileReader zeroDDomainFlows = histFileReader();
	  zeroDDomainFlows.setFileName("netlistFlows_zeroDDomainReplacement.dat");
	  zeroDDomainFlows.setNumColumns(8);
	  zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
	  // Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
	  // ...third column:
	  double flowResult = zeroDDomainFlows.getReadFileData(2,1501);
	  EXPECT_NEAR(259.543214447641,flowResult,1e-8);
	  
	  flowResult = zeroDDomainFlows.getReadFileData(3,1501);
	  EXPECT_NEAR(84.6680132335589,flowResult,1e-8);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader netlist1Flow = histFileReader();
		netlist1Flow.setFileName("netlistFlows_surface_5.dat");
		netlist1Flow.setNumColumns(4);
		netlist1Flow.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = netlist1Flow.getReadFileData(1,1501);
		EXPECT_NEAR(-319.619773209747,flowResult,1e-8);
  }

  // Check netlistFlows_surface_6.dat
  {
		histFileReader netlist2Flow = histFileReader();
		netlist2Flow.setFileName("netlistFlows_surface_6.dat");
		netlist2Flow.setNumColumns(4);
		netlist2Flow.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double flowResult = netlist2Flow.getReadFileData(1,1501);
		EXPECT_NEAR(259.579506526332,flowResult,1e-8);
  }
  
  // Check netlistPressures_surface_5.dat
  {
		histFileReader netlist1Pressure = histFileReader();
		netlist1Pressure.setFileName("netlistPressures_surface_5.dat");
		netlist1Pressure.setNumColumns(5);
		netlist1Pressure.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double pressureResult = netlist1Pressure.getReadFileData(1,1501);
		EXPECT_NEAR(770.012379257636,pressureResult,1e-8);
  }

  // Check netlistPressures_surface_6.dat
  {
		histFileReader netlist1Pressure = histFileReader();
		netlist1Pressure.setFileName("netlistPressures_surface_6.dat");
		netlist1Pressure.setNumColumns(5);
		netlist1Pressure.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double pressureResult = netlist1Pressure.getReadFileData(1,1501);
		EXPECT_NEAR(769.442541551951,pressureResult,1e-8);
  }
}