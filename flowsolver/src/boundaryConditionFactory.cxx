#include "boundaryConditionFactory.hxx"
#include "RCR.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "controlledCoronary.hxx"
#include "datatypesInCpp.hxx"
#include "ClosedLoopDownstreamSubsection.hxx"
#include <boost/weak_ptr.hpp>

boost::shared_ptr<abstractBoundaryCondition> boundaryConditionFactory::createBoundaryCondition (int surfaceIndex, boundary_condition_t boundaryType)
{
  assert(m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt);

  if (boundaryType == BoundaryCondition_RCR)
  {
    boundaryConditionToReturn = boost::shared_ptr<abstractBoundaryCondition> (new RCR(surfaceIndex, m_hstep, m_delt, m_alfi, m_startingTimestepIndex, m_maxsurf, m_nstep));
    
    if (!m_simulationIsPurelyZeroD) {
      boundaryConditionToReturn->getPressureAndFlowPointersFromFortran();
    }

    return boundaryConditionToReturn;
  }
  else if (boundaryType == BoundaryCondition_Netlist)
  {
    // Identify and gather the weak pointers to the downstream loop-closing subcircuits which connect to this boundary condition.
    // These will be passed to the NetlistBoundaryCondition constructor so that it knows which, if any, circuits are downstream of it
    std::vector<boost::weak_ptr<ClosedLoopDownstreamSubsection>> gatheredDownstreamSubcircuits;
    for (auto downstreamSubcircuit = mp_netlistDownstreamLoopClosingSubsectionsWeakPointers.begin(); downstreamSubcircuit != mp_netlistDownstreamLoopClosingSubsectionsWeakPointers.end(); downstreamSubcircuit++)
    {
      if(downstreamSubcircuit->lock()->boundaryConditionCircuitConnectsToThisDownstreamSubsection(surfaceIndex))
      {
        gatheredDownstreamSubcircuits.push_back(*downstreamSubcircuit);
      }
    }

    boundaryConditionToReturn = boost::shared_ptr<abstractBoundaryCondition> (new NetlistBoundaryCondition(surfaceIndex, m_hstep, m_delt, m_alfi, m_startingTimestepIndex, m_maxsurf, m_nstep, gatheredDownstreamSubcircuits));

    if (!m_simulationIsPurelyZeroD) {
      boundaryConditionToReturn->getPressureAndFlowPointersFromFortran();
      // We can't yet initialise the model if this is a pure zero-d simulation,
      // because the boundary pressure/flow pointers aren't ready yet.
      // (because normally they're provided by Fortran, and zeroD doesn't use the Fortran code)
      // Responsibility for this will be taken by the pureZeroDDriver instead, in the call to
      // boundaryConditionManager::setZeroDDomainReplacementPressuresAndFlows().
      boundaryConditionToReturn->initialiseModel();
    }

    boost::shared_ptr<NetlistBoundaryCondition> downcastNetlistBoundaryCondition = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (boundaryConditionToReturn);

    // We finish off by giving the ClosedLoopdownstreamSubcircuits pointers to the just-constructed NetlistBoundaryCondition, if the two are connected directly.
    // This completes the allowing of each to access each other using smart pointers.
    for (auto downstreamSubcircuit = gatheredDownstreamSubcircuits.begin(); downstreamSubcircuit != gatheredDownstreamSubcircuits.end(); downstreamSubcircuit++)
    {
      downstreamSubcircuit->lock()->setPointerToNeighbouringBoundaryConditionCircuit(downcastNetlistBoundaryCondition->getNetlistCircuit());
    }
    
    return boundaryConditionToReturn;
  }
  else if (boundaryType == BoundaryCondition_ControlledCoronary)
  {
    boundaryConditionToReturn = boost::shared_ptr<abstractBoundaryCondition> (new controlledCoronary(surfaceIndex, m_hstep, m_delt, m_alfi, m_startingTimestepIndex, m_maxsurf, m_nstep));

    if (!m_simulationIsPurelyZeroD) {
      boundaryConditionToReturn->getPressureAndFlowPointersFromFortran();
    }

    return boundaryConditionToReturn;
  }
  else
  {
    std::cout << "Unknown boundary type. Exiting.\n";
    std::exit(1);
  }
}

// The "venous system" type circuits - which we generically name "netlist loop closing circuits"
// (because we may use these in non-venous contexts) must be built before we start
// constructing the boundary conditions that connect to them.
void boundaryConditionFactory::createNetlistLoopClosingCircuits(std::vector<boost::shared_ptr<ClosedLoopDownstreamSubsection>>& netlistDownstreamLoopClosingSubsections)
{
  assert(!m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt);

  for (int loopClosingCircuitIndex=1; loopClosingCircuitIndex <= m_numLoopClosingNetlistCircuits; loopClosingCircuitIndex++)
  {
    boost::shared_ptr<ClosedLoopDownstreamSubsection> sharedPtrToPushBack(new ClosedLoopDownstreamSubsection(loopClosingCircuitIndex, m_hstep, m_delt, m_alfi, m_startingTimestepIndex));
    netlistDownstreamLoopClosingSubsections.push_back(sharedPtrToPushBack);

    // These weak pointers to the loop closing circuits will be given to the boundary conditions that connect to them.
    boost::weak_ptr<ClosedLoopDownstreamSubsection> weakPtrToPushBack(sharedPtrToPushBack);
    mp_netlistDownstreamLoopClosingSubsectionsWeakPointers.push_back(weakPtrToPushBack);
  }

  m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt = true;
}