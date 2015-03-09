#include "testMainWithZeroDDomain.hxx"

// Hack to force the compiler to link this test to the relevant main() for testing
int PullInMyLibraryTestMainWithZeroDDomain() { return 0; }

TEST_F(testMainWithZeroDDomain, checkWithRCRs) {
  // This test uses a solver.inp which (on purpose) does not take enough
  // iterations (Step Construction 0 1 0 1...) for decent convergence.
  // The reason for this is that the inaccuracy causes the aortic valve to flap
  // open and closed during a very short test - this is fine for us as we
  // just want to test that the valve is doing its job!
  setSimDirectory("mainTests/zeroDDomain/RCRs");
  clearOutOldFiles();

  runSimulation();
  MPI_Barrier(MPI_COMM_WORLD);

  // Check netlistPressures_surface_-1.dat.dat
  {
	  histFileReader zeroDDomainPressures = histFileReader();
	  zeroDDomainPressures.setFileName("netlistPressures_surface_-1.dat");
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
	  zeroDDomainFlows.setFileName("netlistFlows_surface_-1.dat");
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
