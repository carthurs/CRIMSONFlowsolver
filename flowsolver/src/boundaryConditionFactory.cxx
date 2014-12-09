#include "boundaryConditionFactory.hxx"
#include "RCR.hxx"
#include "netlistBoundaryCondition.hxx"
#include "controlledCoronary.hxx"

boost::shared_ptr<abstractBoundaryCondition> boundaryConditionFactory::createBoundaryCondition (int surfaceIndex, std::string boundaryType)
 {

  if (boundaryType.compare("rcr") == int(0))
  {
    
    return boost::shared_ptr<abstractBoundaryCondition> (new RCR(surfaceIndex));
  }
  else if (boundaryType.compare("netlist") == int(0))
  {
    return boost::shared_ptr<abstractBoundaryCondition> (new netlistBoundaryCondition(surfaceIndex));
  }
  else if (boundaryType.compare("controlledCoronary") == int(0))
  {
    return boost::shared_ptr<abstractBoundaryCondition> (new controlledCoronary(surfaceIndex));
  }
  else
  {
    std::cout << "Unknown boundary type. Exiting.\n";
    std::exit(1);
  }
}