#include "NetlistBoundaryCondition.hxx"
#include "fileReaders.hxx"
#include "datatypesInCpp.hxx"
#include <cassert>
#include <algorithm>

// Statics:
int NetlistBoundaryCondition::numberOfInitialisedNetlistLPNs = 0;

void NetlistBoundaryCondition::initialiseModel()
{
    {
        int numberOfPointers = 1;
        try {
            mp_NetlistCircuit->setPointersToBoundaryPressuresAndFlows(pressure_n_ptrs.at(0), flow_n_ptrs.at(0), numberOfPointers);
        } catch (const std::exception& e) {
            std::cout << e.what() << " observed at line " << __LINE__ << " of " << __FILE__ << std::endl;
            throw e;
        }
    }
    // Get the input data
    mp_NetlistCircuit->createCircuitDescription();

    // Determine how many subcircuits are needed, and note which components belong to each subcircuit
    // mp_NetlistCircuit->identifyAtomicSubcircuits();

    // Initialise all diodes to their closed state, for stability
    //\todo change this if you're restarting and the diodes need to be open at restart!
    mp_NetlistCircuit->closeAllDiodes();
    mp_NetlistCircuit->detectWhetherClosedDiodesStopAllFlowAt3DInterface();

    mp_NetlistCircuit->initialiseCircuit();

    // count the diodes, and set up the AtomicSubcircuitConnectionManager, which is used it working out
    // what connections should be made when a diode/valve opens.
    // AtomicSubcircuitConnectionManager* toPassToSharedPtr = new AtomicSubcircuitConnectionManager(mp_CircuitDescription,m_CircuitDataForAtomicSubcircuits);
    // m_atomicSubcircuitConnectionManager = boost::shared_ptr<AtomicSubcircuitConnectionManager>( toPassToSharedPtr );
}

void NetlistBoundaryCondition::finalizeLPNAtEndOfTimestep()
{
    mp_NetlistCircuit->finalizeLPNAtEndOfTimestep();
}

void NetlistBoundaryCondition::initialiseAtStartOfTimestep()
{
    // Idetify and construct the appropriate subcircuits for this timestep
    mp_NetlistCircuit->initialiseAtStartOfTimestep();
}

std::pair<double,double> NetlistBoundaryCondition::computeImplicitCoefficients(const int timestepNumber, const double timeAtStepNplus1, const double alfi_delt)
{
    // Get the implicit coefficients from the identified subcircuit
    std::pair<double,double> implicitCoefficients;
    implicitCoefficients = mp_NetlistCircuit->computeImplicitCoefficients(timestepNumber,timeAtStepNplus1,alfi_delt);

    return implicitCoefficients;
}

std::pair<boundary_data_t,double> NetlistBoundaryCondition::computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber)
{
    std::pair<boundary_data_t,double> pressureOrFlowToReturn;
    pressureOrFlowToReturn = mp_NetlistCircuit->computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(timestepNumber);
    return pressureOrFlowToReturn;
}

void NetlistBoundaryCondition::updateLPN(const int timestepNumber)
{
    mp_NetlistCircuit->updateLPN(timestepNumber);
}

void NetlistBoundaryCondition::writePressuresFlowsAndVolumes(int& nextTimestepWrite_start)
{
    mp_NetlistCircuit->writePressuresFlowsAndVolumes(nextTimestepWrite_start);
}

// void NetlistBoundaryCondition::loadPressuresFlowsAndVolumesOnRestart(const int startingTimeStepIndex)
// {
//     mp_NetlistCircuit->loadPressuresFlowsAndVolumesOnRestart(startingTimeStepIndex);
// }

// Processes the binaryMask for setting Dirichlet conditions.
// This boundary condition knows which mesh nodes lie at its surface (checked by the assert),
// and it sets 0 in binaryMask at the appropriate location for these nodes, if the boundary
// condition type is currently Dirichlet.
void NetlistBoundaryCondition::setDirichletConditionsIfNecessary(int* const binaryMask)
{
  if(flowPermittedAcross3DInterface())
  {
    assert(hasListOfMeshNodesAtThisBoundary);
    // set zero in the binaryMask at the locations necessary to impose Dirichlet at this surface
    for (auto node=listOfMeshNodesAtThisBoundary.begin(); node!=listOfMeshNodesAtThisBoundary.end(); node++)
    {
      binaryMask[*node] = 0;
    }
  }
}

void NetlistBoundaryCondition::setPressureAndFlowPointers(double* pressurePointer, double* flowPointer)
{
    // This needs tidying up; currently the pressure_n_ptrs and flow_n_ptrs
    // are duplicated here and in the NetlistCircuit.
    //
    // We may not even need vectors here any more...
    // mp_NetlistCircuit->setPressureAndFlowPointers(pressurePointer,flowPointer);
    flow_n_ptrs.clear();
    flow_n_ptrs.push_back(flowPointer);

    pressure_n_ptrs.clear();
    pressure_n_ptrs.push_back(pressurePointer);
}

bool NetlistBoundaryCondition::flowPermittedAcross3DInterface()
{
    return mp_NetlistCircuit->flowPermittedAcross3DInterface();
}

bool NetlistBoundaryCondition::boundaryConditionTypeHasJustChanged()
{
    return mp_NetlistCircuit->boundaryConditionTypeHasJustChanged();
}

boost::shared_ptr<CircuitComponent> NetlistBoundaryCondition::getComponentByInputDataIndex(const int componentIndex)
{
    return mp_NetlistCircuit->getComponentByInputDataIndex(componentIndex);
}

boost::shared_ptr<NetlistCircuit> NetlistBoundaryCondition::getNetlistCircuit()
{
    return mp_NetlistCircuit;
}

bool NetlistBoundaryCondition::hasPrescribedPressureAcross3DInterface() const
{
    return mp_NetlistCircuit->hasPrescribedPressureAcross3DInterface();
}

bool NetlistBoundaryCondition::hasPrescribedFlowAcross3DInterface() const
{
    return mp_NetlistCircuit->hasPrescribedFlowAcross3DInterface();
}