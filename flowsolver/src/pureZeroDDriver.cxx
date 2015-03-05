#include "pureZeroDDriver.hxx"
#include "fileIOHelpers.hxx"
#include <cmath>
#include <sstream>
#include <stdexcept>

void PureZeroDDriver::init()
{
	m_timestepNumber = 0;

	double oneResistanceToGiveEachResistor = 0.001;
	double elastanceToGiveVolumeTrackingPressureChamber = 0.0102461; // to give ~10mmHg pressure at the initial vol
	double initialDomainPressure = 133.3 * 10; // 80 mmHg
	int negativeIndexForNetlistThatReplaces3DDomain = -1;
	m_zeroDDomainLPN = boost::shared_ptr<Netlist3DDomainReplacement> (new Netlist3DDomainReplacement(negativeIndexForNetlistThatReplaces3DDomain, oneResistanceToGiveEachResistor, elastanceToGiveVolumeTrackingPressureChamber, initialDomainPressure));

	assert(boundaryConditionManager_instance->getNumberOfNetlistSurfaces() > 0);
	// Arrays to park the boundary pressures and flows generated by the zero-D replacement for the 3D domain
	// We do this so we can recycle the code that expects boundary flow/pressure pointers from Fortran from the 3D domain case.
	mp_interfaceFlowsToBeReadByBoundaryConditions = new double[boundaryConditionManager_instance->getNumberOfNetlistSurfaces()];
	mp_interfacePressuresToBeReadByBoundaryConditions = new double[boundaryConditionManager_instance->getNumberOfNetlistSurfaces()];

	// Arrays to park boundary pressures and flows generated by the boundary conditions, to be read by the 3D domain replacement Netlist:
	mp_interfaceFlowsToBeReadBy3DDomainReplacement = new double[boundaryConditionManager_instance->getNumberOfNetlistSurfaces()];
	mp_interfacePressuresToBeReadBy3DDomainReplacement = new double[boundaryConditionManager_instance->getNumberOfNetlistSurfaces()];

	// Initialise to zero:
	for (int ii=0; ii<boundaryConditionManager_instance->getNumberOfNetlistSurfaces(); ii++)
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

	// m_zeroDDomainLPN->m_IndexOfThisNetlistLPN = 0; // \todo remove, this is a hack to get the fake "file reader" to read from the right place
	m_zeroDDomainLPN->initialiseModel();
	m_zeroDDomainLPN->setPointersToBoundaryPressuresAndFlows(mp_interfacePressuresToBeReadBy3DDomainReplacement, mp_interfaceFlowsToBeReadBy3DDomainReplacement, boundaryConditionManager_instance->getNumberOfNetlistSurfaces());
	// m_zeroDDomainLPN->setPointersToBoundaryPressuresAndFlows(mp_interfacePressuresToBeReadByBoundaryConditions,mp_interfaceFlowsToBeReadByBoundaryConditions);
}

void PureZeroDDriver::iter_init()
{
	std::cout << "============ Doing timestep number " << m_timestepNumber << " ============" << std::endl;
	m_zeroDDomainLPN->initialiseAtStartOfTimestep();
	boundaryConditionManager_instance->initialiseLPNAtStartOfTimestep_netlist();
}

void PureZeroDDriver::iter_step()
{
	// Also need to actually solve the boundary conditions etc. here!

	m_pressuresOrFlowsAtBoundaries = boundaryConditionManager_instance->getBoundaryPressuresOrFlows_zeroDDomainReplacement(m_timestepNumber);
	for (int ii =0; ii<3; ii++)
	{
		std::cout << "Gave 0D domain: " << ii << " " << m_pressuresOrFlowsAtBoundaries.at(ii).first << " " << m_pressuresOrFlowsAtBoundaries.at(ii).second << std::endl;
	}
	m_zeroDDomainLPN->setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(m_pressuresOrFlowsAtBoundaries);
	placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement();
	m_zeroDDomainLPN->solveSystem(m_timestepNumber);

	std::vector<double> boundaryPressures = m_zeroDDomainLPN->getBoundaryPressures();
	std::vector<double> boundaryFlows = m_zeroDDomainLPN->getBoundaryFlows();
	for (int ii =0; ii<3; ii++)
	{
		std::cout << "Gave Boundary Conditions: " << ii << " " << boundaryPressures.at(ii) << " " << boundaryFlows.at(ii) << std::endl;
	}
	placePressuresAndFlowsInStorageArrays_toGiveToBoundaryConditions(boundaryPressures, boundaryFlows);

}

void PureZeroDDriver::iter_finalize()
{
	m_zeroDDomainLPN->updateLPN();
	m_zeroDDomainLPN->finalizeLPNAtEndOfTimestep();
	std::vector<boost::shared_ptr<abstractBoundaryCondition>> oneElementVectorToPass;
	// boost::shared_ptr<Netlist3DDomainReplacement> castToAbstractBC(m_zeroDDomainLPN.get());
	oneElementVectorToPass.push_back(m_zeroDDomainLPN);
	writeNetlistFlowsPressuresAndVolumes(oneElementVectorToPass, m_nextTimestepWrite_zeroDBoundaries_start);
	
	boundaryConditionManager_instance->updateAllNetlistLPNs();
	boundaryConditionManager_instance->finalizeLPNAtEndOfTimestep_netlists();
	boundaryConditionManager_instance->writeAllNetlistComponentFlowsAndNodalPressures();

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
		mp_interfacePressuresToBeReadByBoundaryConditions[boundaryConditionIndex] = boundaryPressures.at(boundaryConditionIndex);
	}
}

void PureZeroDDriver::placePressuresAndFlowsInStorageArrays_toGiveTo3DDomainReplacement()
{
	{
		auto boundaryPressureOrFlow = m_pressuresOrFlowsAtBoundaries.begin();
		int boundaryConditionIndex = 0;
		while(boundaryPressureOrFlow != m_pressuresOrFlowsAtBoundaries.end())
		{
			if (boundaryPressureOrFlow->first == Boundary_Flow)
			{
				mp_interfaceFlowsToBeReadBy3DDomainReplacement[boundaryConditionIndex] = boundaryPressureOrFlow->second;
				// Set NaN in the pressure storage for this point, so that we notice if we accidentally read it for this surface!
				// (we should only ever read the flow at location boundaryConditionIndex)
				mp_interfacePressuresToBeReadBy3DDomainReplacement[boundaryConditionIndex] = NAN;
			}
			else if (boundaryPressureOrFlow->first == Boundary_Pressure)
			{
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

			// Loop housekeeping:
			boundaryPressureOrFlow++;
			boundaryConditionIndex++;
		}
	}
}