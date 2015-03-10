#include "boundaryConditionManager.hxx"
#include "RCR.hxx"
#include "controlledCoronary.hxx"
#include "NetlistBoundaryCondition.hxx"
#include "fortranPointerManager.hxx"
#include "fileWriters.hxx"
#include "fileIOHelpers.hxx"

// This file contains (and should continue to contain) all the tools needed to control the boundary conditions.
//
// This includes functions which can be called from Fortran, and should be the sole point of interface between Fortran and C++
// for the boundary conditions, as far as is possible.
//
// One thing that might be though of as an exception to this rule is the fortranBoundaryDataPointerManager class, but this
// is really a way of setting up the link between Fortran and C++, not so much a way of allowing Fortran to control the BCs.
//
// Another exception is that some of the global data is accessed directly by the C++ classes, such as in the constructor
// for the abstractBoundaryCondition. This is not ideal, and should be phased out slowly so that we have fewer points
// of interface between the two languages.

// Static class static member variables:
boundaryConditionManager* boundaryConditionManager::instance = 0;
histFileReader* boundaryConditionManager::PHistReader = NULL;
int boundaryConditionManager::thisIsARestartedSimulation = 0;

// Functions which affect features of the abstract class:
void boundaryConditionManager::setNumberOfRCRSurfaces(const int numGRCRSrfs)
{
  m_NumberOfRCRSurfaces = numGRCRSrfs;
}
void boundaryConditionManager::setNumberOfControlledCoronarySurfaces(const int numControlledCoronarySrfs)
{
  m_NumberOfControlledCoronarySurfaces = numControlledCoronarySrfs;
}
void boundaryConditionManager::setNumberOfNetlistSurfaces(const int numNetlistLPNSrfs)
{
  m_NumberOfNetlistSurfaces = numNetlistLPNSrfs;
}
void boundaryConditionManager::setDelt(const double delt)
{
  m_delt = delt;
}

void boundaryConditionManager::giveBoundaryConditionsListsOfTheirAssociatedMeshNodes(const int* ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray)
{
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    (*boundaryCondition)->setListOfMeshNodesAtThisBoundary(ndsurf_nodeToBoundaryAssociationArray, lengthOfNodeToBoundaryAssociationArray);
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGiveBoundaryConditionsListsOfTheirAssociatedMeshNodes(const int*& ndsurf_nodeToBoundaryAssociationArray, const int& lengthOfNodeToBoundaryAssociationArray)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->giveBoundaryConditionsListsOfTheirAssociatedMeshNodes(ndsurf_nodeToBoundaryAssociationArray, lengthOfNodeToBoundaryAssociationArray);
}


// RCR Boundary condition specific functions
void boundaryConditionManager::setPressureFromFortran()
{
  // see the called funciton setPressureFromFortran comments for details of what this does.
  for (auto boundaryCondition = boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(RCR))
    {
      RCR* downcastRCR = dynamic_cast<RCR*>(boundaryCondition->get());
      downcastRCR->setPressureFromFortran();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPSetPressureFromFortran()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->setPressureFromFortran();
}

void boundaryConditionManager::getImplicitCoeff_rcr(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // +MAXSURF+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      
      implicitCoeffs_toBeFilled[writeLocation+MAXSURF+1] = (*iterator)->getHop();
      
      writeLocation++;
    }
  }
  
}
// ---WRAPPED BY--->
extern "C" void callCppGetImplicitCoeff_rcr(double*& implicitCoeffs_toBeFilled)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_rcr(implicitCoeffs_toBeFilled);
}

void boundaryConditionManager::updateAllRCRS_Pressure_n1_withflow()
{
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      (*iterator)->updpressure_n1_withflow();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllRCRS_Pressure_n1_withflow()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllRCRS_Pressure_n1_withflow();
}

void boundaryConditionManager::setSurfaceList(const std::vector<std::pair<int,std::string>> surfaceList)
{
  // Build a factory
  boundaryConditionFactory factory;
  
  for (auto iterator = surfaceList.begin(); iterator != surfaceList.end(); iterator++)
  {
    boundaryConditions.push_back(factory.createBoundaryCondition(iterator->first,iterator->second));
  }
}

void boundaryConditionManager::setZeroDDomainReplacementPressuresAndFlows(double* zeroDDomainPressures, double* zeroDDomainFlows)
{
  for (auto boundaryCondition = boundaryConditions.begin(); boundaryCondition != boundaryConditions.end(); boundaryCondition++)
  {
    boost::shared_ptr<NetlistBoundaryCondition> downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
    assert(downcastNetlist != NULL);
    
    // Get the (zero-indexed) Netlist index; this gives us the appropriate location of the pointers in the input variables
    int netlistIndex = downcastNetlist->getIndexAmongstNetlists();
    // Give the appropriate memory addresses of the pressures and flows to this NetlistBC:
    downcastNetlist->setPressureAndFlowPointers(&zeroDDomainPressures[netlistIndex], &zeroDDomainFlows[netlistIndex]);
  }
}

void boundaryConditionManager::ifRestartingLoadNecessaryData()
{
  if (thisIsARestartedSimulation)
  {
    // Load PHistRCR.dat, necessary for setting the pressure data in the 
    // LPN at the boundary when restarting
    if (m_NumberOfRCRSurfaces > 0)
    {
      PHistReader = new histFileReader();
      PHistReader->setFileName("PHistRCR.dat");
      PHistReader->setNumColumns(m_NumberOfRCRSurfaces+1);
      PHistReader->readAndSplitMultiSurfaceRestartFile();
    }
  }
}

std::vector<boost::shared_ptr<abstractBoundaryCondition>>* boundaryConditionManager::getBoundaryConditions()
{
    return &boundaryConditions;
}

void boundaryConditionManager::computeAllImplicitCoeff_solve(const int timestepNumber)
{
  fortranBoundaryDataPointerManager* fortranBoundaryDataPointerManager_instance = fortranBoundaryDataPointerManager::Get();
  for (auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    (*iterator)->computeImplicitCoeff_solve(timestepNumber);
  }
}
// ---WRAPPED BY--->
extern "C" void callCppComputeAllImplicitCoeff_solve(int& timestepNumber)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeAllImplicitCoeff_solve(timestepNumber);
}

void boundaryConditionManager::computeAllImplicitCoeff_update(const int timestepNumber)
{
  for (auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    (*iterator)->computeImplicitCoeff_update(timestepNumber);
  }
}
// ---WRAPPED BY--->
extern "C" void callCppComputeAllImplicitCoeff_update(int& timestepNumber)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->computeAllImplicitCoeff_update(timestepNumber);
}


void boundaryConditionManager::updateAllRCRS_setflow_n(const double* const flows)
{
  int readLocation = 0;
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      (*iterator)->flow_n = flows[readLocation];
      readLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllRCRS_setflow_n(double*& flows)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllRCRS_setflow_n(flows); 
}


void boundaryConditionManager::updateAllRCRS_setflow_n1(const double* const flows)
{
  int readLocation = 0;
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(RCR))
    {
      (*iterator)->flow_n1 = flows[readLocation];
      readLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllRCRS_setflow_n1(double*& flows)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllRCRS_setflow_n1(flows); 
}

void boundaryConditionManager::recordPressuresAndFlowsInHistoryArrays()
{
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    (*iterator)->updatePressureAndFlowHistory();
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPRecordPressuresAndFlowsInHistoryArrays()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->recordPressuresAndFlowsInHistoryArrays();
}

void boundaryConditionManager::writePHistAndQHistRCR()
{
  // Open a file writer to append to Phist
  basicFileWriter phistrcr_writer;
  phistrcr_writer.setFileName("PHistRCR.dat");

  // Open a file writer to append to Qhist
  basicFileWriter qhistrcr_writer;
  qhistrcr_writer.setFileName("QHistRCR.dat");

  // Loop over all the updates since the last restart was written:
  for (int i=timdat.lstep-outpar.ntout+int(1); i<timdat.lstep+int(1); i++)
  {
    phistrcr_writer.writeStepIndex(i);
    qhistrcr_writer.writeStepIndex(i);

    // Loop the boundary conditions looking for the RCRs
    for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
    {
      if (typeid(**iterator)==typeid(RCR))
      {
        // Write the pressure and flow for this timestep (indexed i)
        phistrcr_writer.writeToFile((*iterator)->pressurehist[i]);
        qhistrcr_writer.writeToFile((*iterator)->flowhist[i]);
      }
    }
    phistrcr_writer.writeEndLine();
    qhistrcr_writer.writeEndLine();
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPWritePHistAndQHistRCR()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->writePHistAndQHistRCR();
}


// =========== Controlled Coronary Block ===========

// void boundaryConditionManager::setSurfacePressure_controlledCoronary(double* coronarySurfacePressures)
// {
//   int readLocation = int(0);
//   for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
//   {
//     if (typeid(**iterator)==typeid(controlledCoronary))
//     {
//      (*iterator)->setLPNInflowPressure(coronarySurfacePressures[readLocation]);
//      readLocation++;
//     }
//   }
// }
// // ---WRAPPED BY--->
// extern "C" void callCppSetSurfacePressure_controlledCoronary(double*& coronarySurfacePressures)
// {
//   boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
//   boundaryConditionManager_instance->setSurfacePressure_controlledCoronary(coronarySurfacePressures);
// }

void boundaryConditionManager::getImplicitCoeff_controlledCoronary(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(controlledCoronary))
    {
      
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // +MAXSURF+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      
      implicitCoeffs_toBeFilled[writeLocation+MAXSURF+1] = (*iterator)->getHop();
      
      writeLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppGetImplicitCoeff_controlledCoronary(double*& implicitCoeffs_toBeFilled) 
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_controlledCoronary(implicitCoeffs_toBeFilled);
}

void boundaryConditionManager::updateAllControlledCoronaryLPNs()
{
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(controlledCoronary))
    {
      (*iterator)->updateLPN();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppUpdateAllControlledCoronaryLPNs()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllControlledCoronaryLPNs();
}


void boundaryConditionManager::finalizeLPNAtEndOfTimestep_controlledCoronary()
{
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(controlledCoronary))
    {
      (*iterator)->finalizeLPNAtEndOfTimestep();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppfinalizeLPNAtEndOfTimestep_controlledCoronary()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->finalizeLPNAtEndOfTimestep_controlledCoronary();
}

void boundaryConditionManager::finalizeLPNAtEndOfTimestep_netlists()
{
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      (*iterator)->finalizeLPNAtEndOfTimestep();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCppfinalizeLPNAtEndOfTimestep_netlists()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->finalizeLPNAtEndOfTimestep_netlists();
}


// void boundaryConditionManager::updateAllControlledCoronaryLPNs_Pressure_n1_withflow()
// {
//   for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
//   {
//     if (typeid(**iterator)==typeid(controlledCoronary))
//     {
//       (*iterator)->updpressure_n1_withflow();
//     }
//   }
// }
// // ---WRAPPED BY--->
// extern "C" void callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow()
// {
//   boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
//   boundaryConditionManager_instance->updateAllRCRS_Pressure_n1_withflow();
// }

// ========== Controlled Coronary Block End =========

// ========== Netlist LPN Block Start =========
void boundaryConditionManager::initialiseLPNAtStartOfTimestep_netlist()
{
  for (auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      (*iterator)->initialiseAtStartOfTimestep();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPInitialiseLPNAtStartOfTimestep_netlist()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->initialiseLPNAtStartOfTimestep_netlist();
}


void boundaryConditionManager::updateAllNetlistLPNs()
{
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      (*iterator)->updateLPN();
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllNetlistLPNs()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllNetlistLPNs();
}

std::map<int,std::pair<double,double>> boundaryConditionManager::getImplicitCoeff_netlistLPNs_toPassTo3DDomainReplacement()
{
  std::map<int,std::pair<double,double>> allNetlistImplicitCoefficients;
  
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      std::pair<double,double> thisNetlistsImplicitCoefficients;
      
      thisNetlistsImplicitCoefficients.first = (*iterator)->getdp_dq();
      thisNetlistsImplicitCoefficients.second = (*iterator)->getHop();

      boost::shared_ptr<NetlistBoundaryCondition> downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*iterator);

      allNetlistImplicitCoefficients.insert(std::make_pair(downcastNetlist->getIndexAmongstNetlists(), thisNetlistsImplicitCoefficients));
    }
  }

  return allNetlistImplicitCoefficients;
}

void boundaryConditionManager::getImplicitCoeff_netlistLPNs(double* const implicitCoeffs_toBeFilled)
{
  // This code is a bit tricky, becase FORTRAN/C++ interfacing doesn't yet support passing arrays which are sized
  // at run-time to C++ from FORTRAN. Therefore, I've had to just pass a pointer to the first entry, and then manage
  // dereferencing of that pointer manually to fill the whole array, but with the FORTRAN column-major array structure,
  // as opposed to the C++ row-major standard.
  int writeLocation = 0;
  
  for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      
      implicitCoeffs_toBeFilled[writeLocation] = (*iterator)->getdp_dq();
      // +MAXSURF+1 here to move to the next column of the array (the +1 is annoying, and is because of weird design decisions in old FORTRAN code)
      
      implicitCoeffs_toBeFilled[writeLocation+MAXSURF+1] = (*iterator)->getHop();
      
      writeLocation++;
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetImplicitCoeff_netlistLPNs(double*& implicitCoeffs_toBeFilled) 
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getImplicitCoeff_netlistLPNs(implicitCoeffs_toBeFilled);
}

// The purpose of this function is to detect when flow across the 3D interface is blocked due
// to diode closure (or due to flows being prescribed -- still \todo).
// It provides the Fortran code with an array of zeros and ones, one for each boundary node in the mesh;
// a 1 indicates that the boundary should be left as-is from the initial input data state (meaning flow is allowed, and the BC is Neumann),
// whereas a 0 indicates that it should be switched to Dirichlet.
void boundaryConditionManager::getBinaryMaskToAdjustNodalBoundaryConditions(int* const binaryMask, const int binaryMaskLength)
{
  // Begin by setting the binary mask to all ones (i.e. flagging 1 to set Dirichlet boundary conditions at all nodes)
  // ... we will set zeros where we want Neumann conditions in a moment...
  for (int maskLocation=0; maskLocation<binaryMaskLength; maskLocation++)
  {
    binaryMask[maskLocation] = 1;
  }
  // Ask the boundary conditions to set zeros where they want Neumann conditions
  for (auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
  {
    if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
    {
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(iterator->get());
      downcastNetlist->setDirichletConditionsIfNecessary(binaryMask);
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetBinaryMaskToAdjustNodalBoundaryConditions(int*& binaryMask, const int& binaryMaskLength)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getBinaryMaskToAdjustNodalBoundaryConditions(binaryMask, binaryMaskLength);
}

void boundaryConditionManager::getNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(int& numBCsWhichDisallowFlow)
{
  // Ensure we start from zero, before we count the surfaces which disallow flow due to closed valves (so we have to switch to Dirichlet)
  numBCsWhichDisallowFlow = 0;
  // ...do the counting (currently only netlists have valves...)
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition))
    {
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (!downcastNetlist->flowPermittedAcross3DInterface())
      {
        numBCsWhichDisallowFlow++;
      }
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(int& numBCsWhichDisallowFlow)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(numBCsWhichDisallowFlow);
}

void boundaryConditionManager::getNumberOfNetlistBoundaryConditionsWhichCurrentlyAllowFlow(int& numBCsWhichAllowFlow)
{
  // Ensure we start from zero, before we count the surfaces which allow flow due to closed valves (so we have to switch to Dirichlet)
  numBCsWhichAllowFlow = 0;
  // ...do the counting (currently only netlists have valves...)
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition))
    {
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (downcastNetlist->flowPermittedAcross3DInterface())
      {
        numBCsWhichAllowFlow++;
      }
    }
  }
}
// ---WRAPPED BY--->
extern "C" void callCPPGetNumberOfNetlistsWhichCurrentlyAllowFlow(int& numBCsWhichAllowFlow)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->getNumberOfNetlistBoundaryConditionsWhichCurrentlyAllowFlow(numBCsWhichAllowFlow);
}

// Takes a surface index and a reference to an int - if flow is allowed across this surface,
// returns 1 in the referenced int, otherwise returns zero in that int.
void boundaryConditionManager::discoverWhetherFlowPermittedAcrossSurface(const int& queriedSurfaceIndex, int& flowIsPermitted)
{
  // Begin by assuming flow is permitted; this will be changed below if flow is not permitted.
  flowIsPermitted = 1;
  // find the queried surface:
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    bool thisIsANetlist = typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition);
    if ((*boundaryCondition)->surfaceIndex == queriedSurfaceIndex && thisIsANetlist)
    {
      // Discover whether we should report that flow is permitted or not:
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (!(downcastNetlist->flowPermittedAcross3DInterface()))
      {
        flowIsPermitted = 0;
      }
    }
  }
}
//---WRAPPED BY--->
extern "C" void callCPPDiscoverWhetherFlowPermittedAcrossSurface(const int& queriedSurfaceIndex, int& flowIsPermitted)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->discoverWhetherFlowPermittedAcrossSurface(queriedSurfaceIndex, flowIsPermitted);
}

void boundaryConditionManager::haveBoundaryConditionTypesChanged(int& boundaryConditionTypesHaveChanged)
{
  // Warning: This thing is broken: the booleans it checks for are not reset correctly currently. Do not use without fixing first!
  assert(false);
  // Begin by assuming boundary conditions are as they were on the previous time-step; this will be changed below if the assumption is false.
  boundaryConditionTypesHaveChanged = 0;
  // find netlists
  for (auto boundaryCondition=boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition)==typeid(NetlistBoundaryCondition))
    {
      // Discover whether we should report a change in boundary condition type (Neumann/Dirichlet):
      NetlistBoundaryCondition* downcastNetlist = dynamic_cast<NetlistBoundaryCondition*>(boundaryCondition->get());
      if (downcastNetlist->boundaryConditionTypeHasJustChanged())
      {
        boundaryConditionTypesHaveChanged = 1;
      }
    }
  }
}
//---WRAPPED BY--->
extern "C" void callCPPHaveBoundaryConditionTypesChanged(int& boundaryConditionTypesHaveChanged)
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->haveBoundaryConditionTypesChanged(boundaryConditionTypesHaveChanged);
}

// void boundaryConditionManager::setSurfacePressure_netlistLPNs(double* netlistSurfacePressures)
// {
//   int readLocation = int(0);
//   for(auto iterator=boundaryConditions.begin(); iterator!=boundaryConditions.end(); iterator++)
//   {
//     if (typeid(**iterator)==typeid(NetlistBoundaryCondition))
//     {
//      (*iterator)->setLPNInflowPressure(netlistSurfacePressures[readLocation]);
//      readLocation++;
//     }
//   }
// }
// // ---WRAPPED BY--->
// extern "C" void callCppSetSurfacePressure_netlistLPNs(double*& netlistSurfacePressures)
// {
//   boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
//   boundaryConditionManager_instance->setSurfacePressure_netlistLPNs(netlistSurfacePressures);
// }

void boundaryConditionManager::writeAllNetlistComponentFlowsAndNodalPressures()
{
  writeNetlistFlowsPressuresAndVolumes(boundaryConditions, m_nextTimestepWrite_netlistBoundaries_start);
}
// ---WRAPPED BY--->
extern "C" void callCPPWriteAllNetlistComponentFlowsAndNodalPressures()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->writeAllNetlistComponentFlowsAndNodalPressures();
}

// Control systems specific functions
void boundaryConditionManager::updateAllControlSystems()
{
  mp_controlSystemsManager->updateAllControlSystems();
}
// ---WRAPPED BY--->
extern "C" void callCPPUpdateAllControlSystems()
{
  boundaryConditionManager* boundaryConditionManager_instance = boundaryConditionManager::Instance();
  boundaryConditionManager_instance->updateAllControlSystems();
}

void boundaryConditionManager::createControlSystems()
{
  // Instantiate the manager
  mp_controlSystemsManager = std::unique_ptr<ControlSystemsManager>(new ControlSystemsManager(m_delt));
  
  // Get the reader class for the netlist data file, and ask it for the control description data:
  netlistReader* netlistReader_instance = netlistReader::Instance();
  
  // Get info for the components that need control (number of these, the component indices in the netlist, and the control types for each)
  // std::vector<int> numberOfComponentsWithControl = getNumberOfComponentsWithControl();
  std::vector<std::map<int,parameter_controller_t>> mapsOfComponentControlTypes = netlistReader_instance->getMapsOfComponentControlTypesForEachSurface();

  // Get info for the nodes that need control (number of these, the nodes indices in the netlist, and the control types for each)
  // std::vector<int> numberOfNodesWithControl = getNumberOfNodesWithControl();
  std::vector<std::map<int,parameter_controller_t>> mapsOfNodeControlTypes = netlistReader_instance->getMapsOfNodalControlTypesForEachSurface();


  // Check for the existence of netlists with input data setting up control of any of 
  // the components. If any are found, initialise the control appropriately.
  for (auto boundaryCondition = boundaryConditions.begin(); boundaryCondition!=boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition) == typeid(NetlistBoundaryCondition))
    {
      // Downcast to a shared_ptr to a Netlist:
      boost::shared_ptr<NetlistBoundaryCondition> currentNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
      // We now initialise all the controls which affect this netlist...
      int netlistIndex = currentNetlist->getIndexAmongstNetlists();
      // Create the controls for components by looping over the pairs which give the component index in the netlist, together with its prescribed control type from netlist_surfaces.dat:
      for (auto componentIndexAndControlType = mapsOfComponentControlTypes.at(netlistIndex).begin(); componentIndexAndControlType != mapsOfComponentControlTypes.at(netlistIndex).end(); componentIndexAndControlType++)
      {
        mp_controlSystemsManager->createParameterController(componentIndexAndControlType->second, currentNetlist, componentIndexAndControlType->first);
      }
      // Create the controls for nodes by looping over the pairs which give the component index in the netlist, together with its prescribed control type from netlist_surfaces.dat:
      for (auto nodeIndexAndControlType = mapsOfNodeControlTypes.at(netlistIndex).begin(); nodeIndexAndControlType != mapsOfNodeControlTypes.at(netlistIndex).end(); nodeIndexAndControlType++)
      {
        mp_controlSystemsManager->createParameterController(nodeIndexAndControlType->second, currentNetlist, nodeIndexAndControlType->first);
      }
    }
  }

}

std::vector<std::pair<boundary_data_t,double>> boundaryConditionManager::getBoundaryPressuresOrFlows_zeroDDomainReplacement(const int timestepNumber)
{
  std::vector<std::pair<boundary_data_t,double>> pressuresOrFlowsAsAppropriate;
  for (auto boundaryCondition = boundaryConditions.begin(); boundaryCondition != boundaryConditions.end(); boundaryCondition++)
  {
    if (typeid(**boundaryCondition) == typeid(NetlistBoundaryCondition))
    {
      boost::shared_ptr<NetlistBoundaryCondition> downcastNetlist = boost::dynamic_pointer_cast<NetlistBoundaryCondition> (*boundaryCondition);
      pressuresOrFlowsAsAppropriate.push_back(downcastNetlist->computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(timestepNumber));
    }
    else
    {
      std::stringstream errorMessage;
      errorMessage << "EE: You can only use a zero-D replacement for the 3D domain if all the boundary conditions are Netlists." << std::endl;
      throw std::runtime_error(errorMessage.str());
    }
  }
  assert(pressuresOrFlowsAsAppropriate.size() == m_NumberOfNetlistSurfaces);
  return pressuresOrFlowsAsAppropriate;
}