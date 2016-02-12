/*

 *
 *  Created on: Oct 7, 2014
 *      Author: klau, carthurs
 */

#include "common_c.h"
#include "multidom.hxx"
#include "fortranPointerManager.hxx"
#include "fileWriters.hxx"
#include <typeinfo>
#include "fileReaders.hxx"
#include "boundaryConditionManager.hxx"
#include <boost/filesystem.hpp>


void multidom_initialise(){

  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();

  // Global data is a terrible idea, so we turn int into class data instead.
  // Pass the necessary variables in to the boundary condition manager
  // ...please try not to access anything in the common block at any level
  // lower than the boundaryConditionManager - if you need something from 
  // the common block to use in a boundary condition, make a set method here
  // and pass it in explicitly to the BC. It's much easier to keep track of this way.
  boundaryConditionManager_instance->setSimulationModePurelyZeroD(nomodule.pureZeroDSimulation);
  boundaryConditionManager_instance->setDelt(inpdat.Delt[0]);
  boundaryConditionManager_instance->setHstep(inpdat.nstep[0] + timdat.lstep);
  boundaryConditionManager_instance->setAlfi(timdat.alfi);
  // boundaryConditionManager_instance->setLstep(timdat.lstep);
  boundaryConditionManager_instance->setNtout(outpar.ntout);
  boundaryConditionManager_instance->setMaxsurf(MAXSURF);
  boundaryConditionManager_instance->setNstep(inpdat.nstep[0]);
  boundaryConditionManager_instance->setNumberOfRCRSurfaces(grcrbccom.numGRCRSrfs);
  boundaryConditionManager_instance->setNumberOfControlledCoronarySurfaces(nomodule.numControlledCoronarySrfs);
  boundaryConditionManager_instance->setNumberOfNetlistSurfaces(nomodule.numNetlistLPNSrfs);
  boundaryConditionManager_instance->setNumLoopClosingnetlistCircuits(nomodule.numLoopClosingCircuits);
  boundaryConditionManager_instance->setMasterControlScriptPresent(nomodule.hasMasterPythonControlScript);

  boundaryConditionManager_instance->ifRestartingLoadNecessaryData();

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
    NetlistReader* netlistReader_instance = NetlistReader::Instance();
    if (boost::filesystem::exists(boost::filesystem::path("netlist_surfaces.dat")))
    {
      netlistReader_instance->setFileName("netlist_surfaces.dat");
      netlistReader_instance->readAndSplitMultiSurfaceInputFile();
    }

    // for converting old netlist specification file format to new (generally not important for actual simulations)
    netlistReader_instance->writeCircuitSpecificationInXmlFormat();
  }

  if (boundaryConditionManager_instance->getNumberOfDownsreamClosedLoopCircuits() > 0)
  {
    NetlistReader* closedLoopDownstreamReader_instance = NetlistDownstreamCircuitReader::Instance();
    if (boost::filesystem::exists(boost::filesystem::path("netlist_closed_loop_downstream.dat")))
    {
      closedLoopDownstreamReader_instance->setFileName("netlist_closed_loop_downstream.dat");
      closedLoopDownstreamReader_instance->readAndSplitMultiSurfaceInputFile();
      // for converting old netlist specification file format to new (generally not important for actual simulations)
      NetlistDownstreamCircuitReader* downcastDownstreamCircuitReader = static_cast<NetlistDownstreamCircuitReader*> (closedLoopDownstreamReader_instance);
      downcastDownstreamCircuitReader->writeDownstreamCircuitSpecificationInXmlFormat();
    }
  }

  // Assemble the list of global surface numbers and types. This will be used
  // by the boundaryConditionFactory to build the boundary conditions.
  std::vector<std::pair<int,boundary_condition_t>> surfaceList;
  for (int ii = 0; ii < boundaryConditionManager_instance->getNumberOfRCRSurfaces(); ii++)
  {
    surfaceList.push_back(std::pair <int,boundary_condition_t> (grcrbccom.nsrflistGRCR[ii+1],BoundaryCondition_RCR));
  }
  for (int ii = 0; ii < boundaryConditionManager_instance->getNumberOfControlledCoronarySurfaces(); ii++)
  {
    surfaceList.push_back(std::pair <int,boundary_condition_t> (nomodule.indicesOfCoronarySurfaces[ii+1],BoundaryCondition_ControlledCoronary));
  }
  for (int ii = 0; ii < boundaryConditionManager_instance->getNumberOfNetlistSurfaces() ; ii++)
  {
    surfaceList.push_back(std::pair<int,boundary_condition_t> (nomodule.indicesOfNetlistSurfaces[ii+1],BoundaryCondition_Netlist));
  }
  // Write loops here for all the other surface types!

  boundaryConditionManager_instance->setSurfaceList(surfaceList);
}

void multidomSetupControlSystems()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();

  if (boundaryConditionManager_instance->getNumberOfNetlistSurfaces() > 0)
  {
    boundaryConditionManager_instance->createControlSystems();
  }
}



void multidom_iter_initialise()
{
}

void multidom_iter_step()
{
}

void multidom_iter_finalise()
{
  boundaryConditionManager::Instance()->markClosedLoopLinearSystemsForRebuilding();
  boundaryConditionManager::Instance()->incrementTimestepIndex();
  // boundaryConditionManager::Instance()->storeAllBoundaryConditionFlowsAndPressuresAtStartOfTimestep();
}

void multidom_finalise(){
  boundaryConditionManager::Term();
  rcrtReader::Term();
  controlledCoronaryReader::Term();
  NetlistReader::Term();
  NetlistDownstreamCircuitReader::Term();
  NetlistXmlReader::Term();
}

