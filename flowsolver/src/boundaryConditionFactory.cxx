#include "boundaryConditionFactory.hxx"
#include "RCR.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "controlledCoronary.hxx"

boost::shared_ptr<abstractBoundaryCondition> boundaryConditionFactory::createBoundaryCondition (int surfaceIndex, std::string boundaryType)
{
  if (boundaryType.compare("rcr") == 0)
  {
    boost::shared_ptr<abstractBoundaryCondition> newRCR(new RCR(surfaceIndex, m_hstep, m_delt, m_alfi, m_lstep, m_maxsurf, m_nstep));
    
    newRCR->getPressureAndFlowPointersFromFortran();

    return newRCR;
  }
  else if (boundaryType.compare("netlist") == 0)
  {
    boost::shared_ptr<abstractBoundaryCondition> newNetlist(new NetlistBoundaryCondition(surfaceIndex, m_hstep, m_delt, m_alfi, m_lstep, m_maxsurf, m_nstep));

    newNetlist->getPressureAndFlowPointersFromFortran();
    
    newNetlist->initialiseModel();
    
    return newNetlist;
  }
  else if (boundaryType.compare("controlledCoronary") == 0)
  {
    boost::shared_ptr<abstractBoundaryCondition> newControlledCoronary(new controlledCoronary(surfaceIndex, m_hstep, m_delt, m_alfi, m_lstep, m_maxsurf, m_nstep));

    newControlledCoronary->getPressureAndFlowPointersFromFortran();

    return newControlledCoronary;
  }
  else
  {
    std::cout << "Unknown boundary type. Exiting.\n";
    std::exit(1);
  }
}