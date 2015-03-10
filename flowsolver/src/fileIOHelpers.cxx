#include "fileIOHelpers.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "fileWriters.hxx"
#include <sstream>

// nextTimestepWrite_start will be updated and returned to caller of this function.
void writeNetlistFlowsPressuresAndVolumes(const std::vector<boost::shared_ptr<abstractBoundaryCondition>>& boundaryConditions, int& nextTimestepWrite_start)
{
  int nextTimestepWrite_inout;
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    NetlistBoundaryCondition* netlistBoundaryCondition = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
    if (netlistBoundaryCondition != NULL) // if the boundaryCondition was a NetlistBoundaryCondition (or subclass thereof)
    {
      nextTimestepWrite_inout = nextTimestepWrite_start;
      // nextTimestepWrite_inout will be updated and returned to caller of this function.
      // but we need to pass the same value each time, hence the use of nextTimestepWrite_inout.
      //\todo clean up this mess!
      netlistBoundaryCondition->writePressuresFlowsAndVolumes(nextTimestepWrite_inout);

    }
  }

  // Set the last value of nextTimestepWrite_inout, modified by writePressuresFlowsAndVolumes, to be nextTimestepWrite_start,
  // so that we know which timestep to start writing data fromon the next call to this function
  nextTimestepWrite_start = nextTimestepWrite_inout;
}