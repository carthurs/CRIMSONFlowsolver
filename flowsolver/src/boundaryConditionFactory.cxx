#include "boundaryConditionFactory.hxx"
#include "RCR.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "controlledCoronary.hxx"
#include "datatypesInCpp.hxx"

boost::shared_ptr<abstractBoundaryCondition> boundaryConditionFactory::createBoundaryCondition (int surfaceIndex, boundary_condition_t boundaryType)
{
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