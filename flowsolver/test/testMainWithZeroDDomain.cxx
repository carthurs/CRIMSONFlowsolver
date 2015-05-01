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
	  EXPECT_NEAR(19268.692618963,pressureResult,1e-8);
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
	  EXPECT_NEAR(-253550.566278505,flowResult,1e-8);
	  // ... 2nd column:
	  flowResult = zeroDDomainFlows.getReadFileData(2,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,1e-8);
	  // ...third column:
	  flowResult = zeroDDomainFlows.getReadFileData(3,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(4,1000);
	  EXPECT_NEAR(-253550.566278505,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(5,1000);
	  EXPECT_NEAR(119440.938383164,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(6,1000);
	  EXPECT_NEAR(119440.938383163,flowResult,1e-8);

	  flowResult = zeroDDomainFlows.getReadFileData(7,1000);
	  EXPECT_NEAR(14668.6895121776,flowResult,1e-8);
  }

  // Check netlistFlows_surface_5.dat
  {
		histFileReader zeroDDomainFlows = histFileReader();
		zeroDDomainFlows.setFileName("netlistFlows_surface_5.dat");
		zeroDDomainFlows.setNumColumns(6);
		zeroDDomainFlows.readAndSplitMultiSurfaceRestartFile();
		// Get the data from timestep 5, 1st column (this method searches for the timestep by value, whereas the columns are zero-indexed)
		double pressureResult = zeroDDomainFlows.getReadFileData(1,1000);
		EXPECT_NEAR(-245488.583296014,pressureResult,1e-8);
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
