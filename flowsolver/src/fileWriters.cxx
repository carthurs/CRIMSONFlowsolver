#include "fileWriters.hxx"

void basicFileWriter::writeHistories_rcr(std::vector<std::unique_ptr<boundaryCondition>> &boundaryConditions)
{


  (*fileHandle) << timdat.lstep;

  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {

      // Write the data for the surface
      (*fileHandle) << " " << (*iterator)->pressurehist[timdat.lstep];
    }
  }

  // Write the end-of-line info
  (*fileHandle) << std::endl;
}