#include "Netlist3DDomainReplacement.hxx"
#include "fileReaders.hxx"

void Netlist3DDomainReplacement::setFlowOrPressurePrescriptionsFromNetlistBoundaryConditions(std::vector<std::pair<boundary_data_t,double>> boundaryFlowsOrPressuresAsAppropriate)
{
	m_boundaryFlowsOrPressuresAsAppropriate = boundaryFlowsOrPressuresAsAppropriate;
    mp_NetlistZeroDDomainCircuit->setBoundaryPrescriptionsAndBoundaryConditionTypes(m_boundaryFlowsOrPressuresAsAppropriate);
}

std::vector<double> Netlist3DDomainReplacement::getBoundaryPressures()
{
	std::vector<double> pressures = mp_NetlistZeroDDomainCircuit->getBoundaryPressures();
	return pressures;
}

std::vector<double> Netlist3DDomainReplacement::getBoundaryFlows()
{
	std::vector<double> flows = mp_NetlistZeroDDomainCircuit->getBoundaryFlows();
	return flows;
}

void Netlist3DDomainReplacement::solveSystem(const int timestepNumber)
{
    mp_NetlistZeroDDomainCircuit->solveSystem(timestepNumber);
}

void Netlist3DDomainReplacement::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    mp_NetlistZeroDDomainCircuit->initialiseAtStartOfTimestep();
}

void Netlist3DDomainReplacement::updateLPN(const int timestepNumber)
{
    // Idetify and construct the appropriate subcircuits for this timestep
    mp_NetlistZeroDDomainCircuit->updateLPN(timestepNumber);
}

void Netlist3DDomainReplacement::finalizeLPNAtEndOfTimestep()
{
    mp_NetlistZeroDDomainCircuit->finalizeLPNAtEndOfTimestep();
}

void Netlist3DDomainReplacement::writePressuresFlowsAndVolumes(int& nextTimestepWrite_zeroDBoundaries_start)
{
    mp_NetlistZeroDDomainCircuit->writePressuresFlowsAndVolumes(nextTimestepWrite_zeroDBoundaries_start);
}

void Netlist3DDomainReplacement::initialiseModel()
{
    // Get the input data
    mp_NetlistZeroDDomainCircuit->createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    mp_NetlistZeroDDomainCircuit->identifyAtomicSubcircuits();

    // chop up the Circuitdata into subcircuits (including removing all mention of diodes from that data now)
    mp_NetlistZeroDDomainCircuit->createAtomicSubcircuitDescriptions();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    // mp_CircuitDescription->closeAllDiodes();
    // mp_CircuitDescription->detectWhetherClosedDiodesStopAllFlowAt3DInterface();
    
    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    // m_atomicSubcircuitConnectionManager = boost::shared_ptr<AtomicSubcircuitConnectionManager>( toPassToSharedPtr );
}

void Netlist3DDomainReplacement::setPointersToBoundaryPressuresAndFlows(double* const interfacePressuresToBeReadBy3DDomainReplacement, double* const interfaceFlowsToBeReadBy3DDomainReplacement, const int& numberOfPointers)
{
    mp_NetlistZeroDDomainCircuit->setPointersToBoundaryPressuresAndFlows(interfacePressuresToBeReadBy3DDomainReplacement, interfaceFlowsToBeReadBy3DDomainReplacement, numberOfPointers);
}

// void Netlist3DDomainReplacement::setPointersToBoundaryPressuresAndFlows(double* const mp_interfacePressures,double* const mp_interfaceFlows)
// {

// }

// std::vector<double> Netlist3DDomainReplacement::getPressureOrFlowPrescriptionsToReturnToNetlistBoundaryConditions()
// {
// 	std::vector<double> pressuresAndFlowsToGiveToBoundaryConditions;
// 	// we know the components at the boundary are just the first m_numberOfNetlistsUsedAsBoundaryConditions in the CircuitData. Loop over them!
// 	for (int indexOfBoundaryInterfaceComponent = 0; indexOfBoundaryInterfaceComponent < m_numberOfNetlistsUsedAsBoundaryConditions; indexOfBoundaryInterfaceComponent++)
// 	{
// 		// Determine whether it's a Dirichlet or Neumann boundary condition at the interface (so pass a pressure or a flow to the boundary condition, respectively).
// 		if (mp_CircuitDescription->hasPrescribedPressureAcrossInterface())
// 		{
// 			// Gather the pressure or flow value into the return variable:
// 			pressuresAndFlowsToGiveToBoundaryConditions.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->flow);
// 		}
// 		else if (mp_CircuitDescription->hasPrescribedFlowAcrossInterface())
// 		{
// 			// Gather the pressure or flow value into the return variable:
// 			pressuresAndFlowsToGiveToBoundaryConditions.push_back(mp_CircuitDescription->components.at(indexOfBoundaryInterfaceComponent)->endNode->pressure);
// 		}
// 		else
// 		{
// 			std::stringstream errorMessage;
// 			errorMessage << "EE: Zero-D domain replacement has at least one bounday interface with neither a prescribed pressure nor a prescribed flow." << std::endl;
// 			throw std::logic_error(errorMessage);			
// 		}
// 	}
// 	return pressuresAndFlowsToGiveToBoundaryConditions;
// }

void Netlist3DDomainReplacement::setDpDqResistances(std::map<int,std::pair<double,double>> allImplicitCoefficients, std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAtBoundaries)
{
    mp_NetlistZeroDDomainCircuit->setDpDqResistances(allImplicitCoefficients,pressuresOrFlowsAtBoundaries);
}

// void Netlist3DDomainReplacement::switchSurfaceToZeroFlow(const int& surfaceIndex)
// {

// }