#include "boundaryConditionFactory.hxx"
#include "RCR.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "controlledCoronary.hxx"
#include "datatypesInCpp.hxx"
#include "NetlistLoopClosingCircuit.hxx"

boost::shared_ptr<abstractBoundaryCondition> boundaryConditionFactory::createBoundaryCondition (int surfaceIndex, boundary_condition_t boundaryType)
{
  assert(m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt);

  if (boundaryType == BoundaryCondition_RCR)
  {
    boundaryConditionToReturn = boost::shared_ptr<abstractBoundaryCondition> (new RCR(surfaceIndex, m_hstep, m_delt, m_alfi, m_lstep, m_maxsurf, m_nstep));
    
    boundaryConditionToReturn->getPressureAndFlowPointersFromFortran();

    return boundaryConditionToReturn;
  }
  else if (boundaryType == BoundaryCondition_Netlist)
  {
    boundaryConditionToReturn = boost::shared_ptr<abstractBoundaryCondition> (new NetlistBoundaryCondition(surfaceIndex, m_hstep, m_delt, m_alfi, m_lstep, m_maxsurf, m_nstep));

    boundaryConditionToReturn->getPressureAndFlowPointersFromFortran();
    
    boundaryConditionToReturn->initialiseModel();
    
    return boundaryConditionToReturn;
  }
  else if (boundaryType == BoundaryCondition_ControlledCoronary)
  {
    boundaryConditionToReturn = boost::shared_ptr<abstractBoundaryCondition> (new controlledCoronary(surfaceIndex, m_hstep, m_delt, m_alfi, m_lstep, m_maxsurf, m_nstep));

    boundaryConditionToReturn->getPressureAndFlowPointersFromFortran();

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
void boundaryConditionFactory::createNetlistLoopClosingCircuits(std::vector<boost::shared_ptr<NetlistLoopClosingCircuit>>& netlistLoopClosingCircuits)
{
  assert(!m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt);

  for (int loopClosingCircuitIndex=1; loopClosingCircuitIndex <= m_numLoopClosingNetlistCircuits; loopClosingCircuitIndex++)
  {
    boost::shared_ptr<NetlistLoopClosingCircuit> toPushBack(new NetlistLoopClosingCircuit(loopClosingCircuitIndex));
    netlistLoopClosingCircuits.push_back(toPushBack);
  }

  m_anyNeededNetlistLoopClosingCircuitsHaveBeenBuilt = true;
}