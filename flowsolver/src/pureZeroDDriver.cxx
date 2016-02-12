#include "pureZeroDDriver.hxx"
#include "fileIOHelpers.hxx"
#include "NetlistXmlReader.hxx"
#include <cmath>
#include <sstream>
#include <stdexcept>
#include "mpi.h"

void PureZeroDDriver::init()
{
	m_numberOfNetlistsUsedAsBoundaryConditions = boundaryConditionManager_instance->getNumberOfNetlistSurfaces();
	assert(m_numberOfNetlistsUsedAsBoundaryConditions > 0);

	double oneResistanceToGiveEachResistor = 0.001;
	double complianceToGiveCentralCapacitor = 0.0102461E1; // to give ~10mmHg pressure at the initial vol
	double initialDomainPressure = 10664.0;//133.3 * 10; // 80 mmHg
	assert(m_hstepHasBeenSet);
	assert(m_mapFromNetlistIndexToConnectedComponentsHasBeenSet);

	m_zeroDDomainLPN = boost::shared_ptr<Netlist3DDomainReplacement> (new Netlist3DDomainReplacement(m_numberOfNetlistsUsedAsBoundaryConditions, oneResistanceToGiveEachResistor, complianceToGiveCentralCapacitor, initialDomainPressure, m_hstep, m_alfi, m_delt, m_mapFromNetlistIndexAmongstNetlistsToConnectedComponentIndex, m_timestepNumber));

	// Arrays to park the boundary pressures and flows generated by the zero-D replacement for the 3D domain
	// We do this so we can recycle the code that expects boundary flow/pressure pointers from Fortran from the 3D domain case.
	mp_interfaceFlowsToBeReadByBoundaryConditions = new double[m_numberOfNetlistsUsedAsBoundaryConditions];
	mp_interfacePressuresToBeReadByBoundaryConditions = new double[m_numberOfNetlistsUsedAsBoundaryConditions];

	// Arrays to park boundary pressures and flows generated by the boundary conditions, to be read by the 3D domain replacement Netlist:
	mp_interfaceFlowsToBeReadBy3DDomainReplacement = new double[m_numberOfNetlistsUsedAsBoundaryConditions];
	mp_interfacePressuresToBeReadBy3DDomainReplacement = new double[m_numberOfNetlistsUsedAsBoundaryConditions];

	// Initialise to zero:
	for (int ii=0; ii<m_numberOfNetlistsUsedAsBoundaryConditions; ii++)
	{
		mp_interfaceFlowsToBeReadByBoundaryConditions[ii] = 0.0;
		mp_interfacePressuresToBeReadByBoundaryConditions[ii] = 0.0;
		
		mp_interfaceFlowsToBeReadBy3DDomainReplacement[ii] = 0.0;
		mp_interfacePressuresToBeReadBy3DDomainReplacement[ii] = 0.0;
	}

	// Give the boundary conditions the pointers to the boundary pressures and flows.
	// We do this like this because it mirrors the way the Fortran 3D domain would give the boundary
	// conditions their pointers to boundary pressures and flows.
	boundaryConditionManager_instance->setZeroDDomainReplacementPressuresAndFlows(mp_interfacePressuresToBeReadByBoundaryConditions,mp_interfaceFlowsToBeReadByBoundaryConditions);

	// m_zeroDDomainLPN->m_IndexOfThisNetlistLPNInInputFile = 0; // \todo remove, this is a hack to get the fake "file reader" to read from the right place
	m_zeroDDomainLPN->initialiseModel();
	m_zeroDDomainLPN->setPointersToBoundaryPressuresAndFlows(mp_interfacePressuresToBeReadBy3DDomainReplacement, mp_interfaceFlowsToBeReadBy3DDomainReplacement, m_numberOfNetlistsUsedAsBoundaryConditions);
	// m_zeroDDomainLPN->loadPressuresFlowsAndVolumesOnRestart(m_timestepNumber);
	// boundaryConditionManager_instance->loadAllNetlistComponentFlowsAndNodalPressures();
}

void PureZeroDDriver::checkIfThisIsARestartedSimulation()
{
  SimpleFileReader numstartReader("numstart.dat");

  bool successfullyReadNumstartDotDat = false;
  std::string numstartString = numstartReader.getNextDataSplitBySpacesOrEndOfLine(successfullyReadNumstartDotDat);
  assert(successfullyReadNumstartDotDat);

  m_timestepNumber = boost::lexical_cast<int>(numstartString);

  if (m_timestepNumber > 0)
  {
    m_thisIsARestartedSimulation = true;
    m_nextTimestepWrite_zeroDBoundaries_start = m_timestepNumber + 1;
  }
  else
  {
    m_thisIsARestartedSimulation = false;
    m_nextTimestepWrite_zeroDBoundaries_start = 0;
  }
}

void PureZeroDDriver::setDelt(const double delt)
{
	m_delt = delt;
	m_deltHasBeenSet = true;
}

void PureZeroDDriver::setAlfi(const double alfi)
{
	m_alfi = alfi;
	m_alfiHasBeenSet = true;
}

void PureZeroDDriver::setHstep(const int hstep)
{
	m_hstep = hstep;
	m_hstepHasBeenSet = true;
}

void PureZeroDDriver::setNtout(const int ntout)
{
	m_ntout = ntout;
	m_ntoutHasBeenSet = true;
}

void PureZeroDDriver::iter_init()
{
	std::cout << "============ Doing timestep number " << m_timestepNumber << " ============" << std::endl;
	assert(m_alfiHasBeenSet);
	assert(m_deltHasBeenSet);
	assert(m_ntoutHasBeenSet);
	m_zeroDDomainLPN->initialiseAtStartOfTimestep();
	boundaryConditionManager_instance->initialiseLPNAtStartOfTimestep_netlist();
	boundaryConditionManager_instance->updateBoundaryConditionControlSystems();
}

void PureZeroDDriver::iter_step()
{
	// Also need to actually solve the boundary conditions etc. here!

	m_pressuresOrFlowsAtBoundaries = boundaryConditionManager_instance->getBoundaryPressuresOrFlows_zeroDDomainReplacement(m_timestepNumber);
	boundaryConditionManager_instance->computeAllImplicitCoeff_solve(m_timestepNumber);
	std::map<int,std::pair<double,double>> allNetlistBoundaryImplicitCoeffs = boundaryConditionManager_instance->getImplicitCoeff_netlistLPNs_toPassTo3DDomainReplacement();
	for (int ii =0; ii<m_numberOfNetlistsUsedAsBoundaryConditions; ii++)
	{
		try {
			if (m_pressuresOrFlowsAtBoundaries.at(ii).first == Boundary_Pressure)
			{
				assert(!isnan(allNetlistBoundaryImplicitCoeffs.at(ii).second));
				m_pressuresOrFlowsAtBoundaries.at(ii).second = allNetlistBoundaryImplicitCoeffs.at(ii).second;
			}
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw e;
		}
	}
	m_zeroDDomainLPN->setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(m_pressuresOrFlowsAtBoundaries);
	m_zeroDDomainLPN->setDpDqResistances(allNetlistBoundaryImplicitCoeffs,m_pressuresOrFlowsAtBoundaries);
	placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement();
	m_zeroDDomainLPN->solveSystem(m_timestepNumber);

	std::vector<double> boundaryPressures = m_zeroDDomainLPN->getBoundaryPressures();
	std::vector<double> boundaryFlows = m_zeroDDomainLPN->getBoundaryFlows();
	
	placePressuresAndFlowsInStorageArrays_toGiveToBoundaryConditions(boundaryPressures, boundaryFlows);

}

void PureZeroDDriver::iter_finalize()
{
	m_zeroDDomainLPN->updateLPN(m_timestepNumber);
	m_zeroDDomainLPN->finalizeLPNAtEndOfTimestep();
	m_zeroDDomainLPN->writePressuresFlowsAndVolumes(m_nextTimestepWrite_zeroDBoundaries_start);
	
	boundaryConditionManager_instance->updateAllNetlistLPNs(m_timestepNumber);
	boundaryConditionManager_instance->finalizeLPNAtEndOfTimestep_netlists();
	bool thisIsAWritingStep = ( m_timestepNumber % m_ntout == 0);
	if (thisIsAWritingStep)
	{
		writeNumstartDotDat();
		boundaryConditionManager_instance->writeAllNetlistComponentFlowsAndNodalPressures();
	}

	m_timestepNumber++;
}

void PureZeroDDriver::finalize()
{
	NetlistXmlReader::Term();
	NetlistDownstreamXmlReader::Term();
	delete[] mp_interfaceFlowsToBeReadByBoundaryConditions;
	delete[] mp_interfacePressuresToBeReadByBoundaryConditions;

	delete[] mp_interfaceFlowsToBeReadBy3DDomainReplacement;
	delete[] mp_interfacePressuresToBeReadBy3DDomainReplacement;
}

void PureZeroDDriver::placePressuresAndFlowsInStorageArrays_toGiveToBoundaryConditions(std::vector<double> boundaryPressures, std::vector<double> boundaryFlows)
{
	for (int boundaryConditionIndex = 0; boundaryConditionIndex < boundaryConditionManager_instance->getNumberOfNetlistSurfaces(); boundaryConditionIndex++)
	{
		try {
			mp_interfaceFlowsToBeReadByBoundaryConditions[boundaryConditionIndex] = boundaryFlows.at(boundaryConditionIndex);
			mp_interfacePressuresToBeReadByBoundaryConditions[boundaryConditionIndex] = boundaryPressures.at(boundaryConditionIndex);
		} catch (const std::exception& e) {
		    std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
		    throw e;
		}
	}
}

void PureZeroDDriver::placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement()
{
	// std::cout << "mp_interfaceFlowsToBeReadBy3DDomainReplacement & mp_interfacePressuresToBeReadBy3DDomainReplacement:" << std::endl;
	{
		auto boundaryPressureOrFlow = m_pressuresOrFlowsAtBoundaries.begin();
		int boundaryConditionIndex = 0;
		while(boundaryPressureOrFlow != m_pressuresOrFlowsAtBoundaries.end())
		{
			if (boundaryPressureOrFlow->first == Boundary_Flow)
			{
				assert(!isnan(boundaryPressureOrFlow->second));
				mp_interfaceFlowsToBeReadBy3DDomainReplacement[boundaryConditionIndex] = boundaryPressureOrFlow->second;
				// Set NaN in the pressure storage for this point, so that we notice if we accidentally read it for this surface!
				// (we should only ever read the flow at location boundaryConditionIndex)
				mp_interfacePressuresToBeReadBy3DDomainReplacement[boundaryConditionIndex] = NAN;
			}
			else if (boundaryPressureOrFlow->first == Boundary_Pressure)
			{
				assert(!isnan(boundaryPressureOrFlow->second));
				mp_interfacePressuresToBeReadBy3DDomainReplacement[boundaryConditionIndex] = boundaryPressureOrFlow->second;
				// Set NaN in the flow storage for this point, so that we notice if we accidentally read it for this surface!
				mp_interfaceFlowsToBeReadBy3DDomainReplacement[boundaryConditionIndex] = NAN;
				// (we should only ever read the pressure at location boundaryConditionIndex)
			}
			else
			{
				std::stringstream errorMessage;
				errorMessage << "EE: Internal error in PureZeroDDriver." << std::endl;
				throw std::logic_error(errorMessage.str());
			}
			// std::cout << mp_interfaceFlowsToBeReadBy3DDomainReplacement[boundaryConditionIndex] << " " << mp_interfacePressuresToBeReadBy3DDomainReplacement[boundaryConditionIndex] << std::endl;


			// Loop housekeeping:
			boundaryPressureOrFlow++;
			boundaryConditionIndex++;
		}
	}
	
}


void PureZeroDDriver::setupConnectedComponents(const int num3DConnectedComponents, const int* const surfacesOfEachConnectedComponent, const int* const indicesOfNetlistSurfaces)
{
	int numberOfNetlistSurfaces = boundaryConditionManager_instance->getNumberOfNetlistSurfaces();
	// Split the zero D domain replacement, in the case
	// where there are multiple connected components of the 3D domain:
	// As a preliminary, build this map. Not strictly necessary here, but it makes things clearer later:
	std::map<int,int> mapFromNetlistSurfaceIndexToIndexAmongstNetlistsInInputData;
	for (int netlistSurfaceIndexInInputData = 0; netlistSurfaceIndexInInputData < numberOfNetlistSurfaces; netlistSurfaceIndexInInputData++)
	{
	  mapFromNetlistSurfaceIndexToIndexAmongstNetlistsInInputData.insert(std::make_pair(indicesOfNetlistSurfaces[netlistSurfaceIndexInInputData+1], netlistSurfaceIndexInInputData)); // +1 because the input_fform.cxx arrays start at their 1th entry (!)
	}

	int connectedComponentIndex = 1;
	// Loop the input data from solver.inp's "List of Surfaces In Each Connected Component Separated by -1s":
	for (int index = 0; index < numberOfNetlistSurfaces + num3DConnectedComponents - 1; index++)
	{
	  // If it's not a separator symbol "-1", used to mark the end of a connected component, 
	  int surfaceIndexOrNextConnectedComponentFlag = surfacesOfEachConnectedComponent[index+1]; // +1 because the input_fform.cxx arrays start at their 1th entry (!)
	  bool isANextConnectedComponentFlag = (surfaceIndexOrNextConnectedComponentFlag == -1);
	  if (!isANextConnectedComponentFlag)
	  {
	    int surfaceIndex = surfaceIndexOrNextConnectedComponentFlag;
	    int indexAmongstNetlists;
	    try {
	    	indexAmongstNetlists = mapFromNetlistSurfaceIndexToIndexAmongstNetlistsInInputData.at(surfaceIndex);
	    } catch (const std::exception& e) {
	        std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
	        throw e;
	    }
	    m_mapFromNetlistIndexAmongstNetlistsToConnectedComponentIndex.insert(std::make_pair(indexAmongstNetlists, connectedComponentIndex));
	    indexAmongstNetlists++;
	  }
	  else
	  {
	    connectedComponentIndex++;
	  }
	}

	m_mapFromNetlistIndexToConnectedComponentsHasBeenSet = true;
}

void PureZeroDDriver::writeNumstartDotDat()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0)
	{
		std::string numstartFileName = std::string("numstart.dat");

		std::ofstream fileHandle = std::ofstream(numstartFileName.c_str(), std::ios::trunc);
		
		if (!fileHandle.is_open())
		{
			std::stringstream errorBuilder;
			errorBuilder << "EE: Failed to open " << numstartFileName << " for writing!" << std::endl;
			throw std::runtime_error(errorBuilder.str());
		}

		fileHandle << m_timestepNumber+1 << std::endl;

		fileHandle.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
}