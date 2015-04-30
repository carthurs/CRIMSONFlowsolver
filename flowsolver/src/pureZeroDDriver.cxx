#include "pureZeroDDriver.hxx"
#include "fileIOHelpers.hxx"
#include <cmath>
#include <sstream>
#include <stdexcept>

void PureZeroDDriver::init()
{
	m_timestepNumber = 0;

	m_numberOfNetlistsUsedAsBoundaryConditions = boundaryConditionManager_instance->getNumberOfNetlistSurfaces();
	assert(m_numberOfNetlistsUsedAsBoundaryConditions > 0);

	double oneResistanceToGiveEachResistor = 0.001;
	double elastanceToGiveVolumeTrackingPressureChamber = 0.0102461E1; // to give ~10mmHg pressure at the initial vol
	double initialDomainPressure = 10664.0;//133.3 * 10; // 80 mmHg
	assert(m_hstepHasBeenSet);
	m_zeroDDomainLPN = boost::shared_ptr<Netlist3DDomainReplacement> (new Netlist3DDomainReplacement(m_numberOfNetlistsUsedAsBoundaryConditions, oneResistanceToGiveEachResistor, elastanceToGiveVolumeTrackingPressureChamber, initialDomainPressure, m_hstep, m_alfi, m_delt));

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
	// m_zeroDDomainLPN->setPointersToBoundaryPressuresAndFlows(mp_interfacePressuresToBeReadByBoundaryConditions,mp_interfaceFlowsToBeReadByBoundaryConditions);
	// for (int ii=0; ii<m_zeroDDomainLPN->getCircuit()->getCircuitDescription()->numberOfComponents; ii++)
	// {
	// 	std::cout << "index component: " << m_zeroDDomainLPN->getCircuit()->getCircuitDescription()->components.at(ii)->getIndex() << std::endl;
	// 	std::cout << "type: " << m_zeroDDomainLPN->getCircuit()->getCircuitDescription()->components.at(ii)->getType() << std::endl;
	// 	std::cout << "parametervalue: " << *(m_zeroDDomainLPN->getCircuit()->getCircuitDescription()->components.at(ii)->getParameterPointer()) << std::endl;
	// 	std::cout << "start node index: " << m_zeroDDomainLPN->getCircuit()->getCircuitDescription()->components.at(ii)->startNode->getIndex() << std::endl;
	// 	std::cout << "end node index: " << m_zeroDDomainLPN->getCircuit()->getCircuitDescription()->components.at(ii)->endNode->getIndex() << std::endl;
	// }
}

void PureZeroDDriver::checkIfThisIsARestartedSimulation()
{
  SimpleFileReader numstartReader("numstart.dat");

  bool success = false;
  std::string numstartString = numstartReader.getNextDataSplitBySpacesOrEndOfLine(success);
  assert(success);

  int valueFromNumstartDotDat = boost::lexical_cast<int>(numstartString);

  if (valueFromNumstartDotDat > 0)
  {
    m_thisIsARestartedSimulation = true;
    m_nextTimestepWrite_zeroDBoundaries_start = valueFromNumstartDotDat + 1; // +1 because numstart should contain the step just written before the program last terminated. So we need to start writing on the next (+1 th) time-step.
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
	boundaryConditionManager_instance->updateAllControlSystems();
}

void PureZeroDDriver::iter_step()
{
	// Also need to actually solve the boundary conditions etc. here!

	m_pressuresOrFlowsAtBoundaries = boundaryConditionManager_instance->getBoundaryPressuresOrFlows_zeroDDomainReplacement(m_timestepNumber);
	boundaryConditionManager_instance->computeAllImplicitCoeff_solve(m_timestepNumber);
	std::map<int,std::pair<double,double>> allNetlistBoundaryImplicitCoeffs = boundaryConditionManager_instance->getImplicitCoeff_netlistLPNs_toPassTo3DDomainReplacement();
	for (int ii =0; ii<m_numberOfNetlistsUsedAsBoundaryConditions; ii++)
	{
		if (m_pressuresOrFlowsAtBoundaries.at(ii).first == Boundary_Pressure)
		{
			assert(!isnan(allNetlistBoundaryImplicitCoeffs.at(ii).second));
			m_pressuresOrFlowsAtBoundaries.at(ii).second = allNetlistBoundaryImplicitCoeffs.at(ii).second;
		}
		// std::cout << "Gave 0D domain: " << ii << " " << m_pressuresOrFlowsAtBoundaries.at(ii).first << " " << m_pressuresOrFlowsAtBoundaries.at(ii).second << std::endl;
		// std::cout << "Computed implicit coeff at surface " << ii << ": " << allNetlistBoundaryImplicitCoeffs.at(ii).first << " " << allNetlistBoundaryImplicitCoeffs.at(ii).second << std::endl;
	}
	m_zeroDDomainLPN->setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(m_pressuresOrFlowsAtBoundaries);
	m_zeroDDomainLPN->setDpDqResistances(allNetlistBoundaryImplicitCoeffs,m_pressuresOrFlowsAtBoundaries);
	placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement();
	m_zeroDDomainLPN->solveSystem(m_timestepNumber);

	std::vector<double> boundaryPressures = m_zeroDDomainLPN->getBoundaryPressures();
	std::vector<double> boundaryFlows = m_zeroDDomainLPN->getBoundaryFlows();
	// for (int ii =0; ii<3; ii++)
	// {
	// 	std::cout << "Gave Boundary Conditions: " << ii << " " << boundaryPressures.at(ii) << " " << boundaryFlows.at(ii) << std::endl;
	// }
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
		boundaryConditionManager_instance->writeAllNetlistComponentFlowsAndNodalPressures();
	}

	m_timestepNumber++;
}

void PureZeroDDriver::finalize()
{
	delete[] mp_interfaceFlowsToBeReadByBoundaryConditions;
	delete[] mp_interfacePressuresToBeReadByBoundaryConditions;

	delete[] mp_interfaceFlowsToBeReadBy3DDomainReplacement;
	delete[] mp_interfacePressuresToBeReadBy3DDomainReplacement;
}

void PureZeroDDriver::placePressuresAndFlowsInStorageArrays_toGiveToBoundaryConditions(std::vector<double> boundaryPressures, std::vector<double> boundaryFlows)
{
	for (int boundaryConditionIndex = 0; boundaryConditionIndex < boundaryConditionManager_instance->getNumberOfNetlistSurfaces(); boundaryConditionIndex++)
	{
		mp_interfaceFlowsToBeReadByBoundaryConditions[boundaryConditionIndex] = boundaryFlows.at(boundaryConditionIndex);
		// std::cout << "set flow: " << boundaryFlows.at(boundaryConditionIndex) << std::endl;
		mp_interfacePressuresToBeReadByBoundaryConditions[boundaryConditionIndex] = boundaryPressures.at(boundaryConditionIndex);
		// std::cout << "set pressure: " << boundaryPressures.at(boundaryConditionIndex) << std::endl;
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