/*

 *
 *  Created on: Oct 7, 2014
 *      Author: klau, carthurs
 */

#include "common_c.h"
#include "multidom.h"
#include "fortranPointerManager.hxx"
#include "fileWriters.hxx"
#include <typeinfo>
#include "fileReaders.hxx"
#include "boundaryConditionManager.hxx"

void multidom_initialise(){

  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();

  // Make the file readers for the classes of surface present in this simulation,
  // and make them read their files:
  if (boundaryConditionManager_instance->getNumberOfRCRSurfaces() > 0)
  {
    rcrtReader* rcrtReader_instance = rcrtReader::Instance();
    rcrtReader_instance->setFileName("rcrt.dat");
    rcrtReader_instance->readAndSplitMultiSurfaceInputFile();
  }

  if (boundaryConditionManager_instance->getNumberOfControlledCoronarySurfaces() > 0)
  {
    controlledCoronaryReader* controlledCoronaryReader_instance = controlledCoronaryReader::Instance();
    controlledCoronaryReader_instance->setFileName("controlled_coronaries.dat");
    controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile();
  }

  if (boundaryConditionManager_instance->getNumberOfNetlistSurfaces() > 0)
  {
    netlistReader* netlistReader_instance = netlistReader::Instance();
    netlistReader_instance->setFileName("netlist_surfaces.dat");
    netlistReader_instance->readAndSplitMultiSurfaceInputFile();
  }

  // Assemble the list of global surface numbers and types. This will be used
  // by the boundaryConditionFactory to build the boundary conditions.
  std::vector<std::pair<int,std::string>> surfaceList;
  for (int ii = 0; ii < boundaryConditionManager_instance->getNumberOfRCRSurfaces(); ii++)
  {
    surfaceList.push_back(std::pair <int,std::string> (grcrbccom.nsrflistGRCR[ii+1],"rcr"));
  }
  for (int ii = 0; ii < boundaryConditionManager_instance->getNumberOfControlledCoronarySurfaces(); ii++)
  {
    surfaceList.push_back(std::pair <int,std::string> (nomodule.indicesOfCoronarySurfaces[ii+1],"controlledCoronary"));
  }
  for (int ii = 0; ii < boundaryConditionManager_instance->getNumberOfNetlistSurfaces() ; ii++)
  {
    surfaceList.push_back(std::pair<int,std::string> (nomodule.indicesOfNetlistSurfaces[ii+1],"netlist"));
  }
  // Write loops here for all the other surface types!

  boundaryConditionManager_instance->setSurfaceList(surfaceList);
  
  // std::vector<boost::shared_ptr<abstractBoundaryCondition>>* retrievedBoundaryConditions;
  // retrievedBoundaryConditions = boundaryConditionManager_instance->getBoundaryConditions();


  // if (grcrbccom.numGRCRSrfs > 0)
  // {
  //  for (int i=0; i < grcrbccom.numGRCRSrfs; i++)
  //  {
  //    std::cout << "writing here" << std::endl;
  //    std::cout << grcrbccom.nsrflistGRCR[i+1] << std::endl;
  //  }
  // }

}


// set pointer to fortran arrays
void multidom_iter_initialise(){

}

void multidom_iter_step(){

}

void multidom_iter_finalise(){

}

void multidom_finalise(){
  boundaryConditionManager::Instance()->Term();
  rcrtReader::Instance()->Term();
  controlledCoronaryReader::Instance()->Term();
  netlistReader::Instance()->Term();
}

