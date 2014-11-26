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

  // Make the file readers for the classes of surface present in this simulation:
  if (grcrbccom.numGRCRSrfs > 0)
  {
    rcrtReader* rcrtReader_instance = rcrtReader::Instance();
    rcrtReader_instance->setFileName("rcrt.dat");
    rcrtReader_instance->readAndSplitMultiSurfaceInputFile();
  }

  if (nomodule.numControlledCoronarySrfs > 0)
  {
    controlledCoronaryReader* controlledCoronaryReader_instance = controlledCoronaryReader::Instance();
    controlledCoronaryReader_instance->setFileName("controlled_coronaries.dat");
    controlledCoronaryReader_instance->readAndSplitMultiSurfaceInputFile();
  }

  std::vector<std::pair<int,std::string>> surfaceList;

  // loop through rcr boundaries listed in the input file, surface numbers read from the common_c.h
  for (int i = 0; i < grcrbccom.numGRCRSrfs; i++)
  {
    surfaceList.push_back(std::pair <int,std::string> (grcrbccom.nsrflistGRCR[i+1],"rcr"));
  }
  for (int ii = 0; ii < nomodule.numControlledCoronarySrfs; ii++)
  {
    surfaceList.push_back(std::pair <int,std::string> (nomodule.indicesOfCoronarySurfaces[ii+1],"controlledCoronary"));
  }
  // Write loops here for all the other surface types!
  
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->setSurfaceList(surfaceList);
  
  std::vector<boost::shared_ptr<abstractBoundaryCondition>>* retrievedBoundaryConditions;
  retrievedBoundaryConditions = boundaryConditionManager_instance->getBoundaryConditions();


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
}

