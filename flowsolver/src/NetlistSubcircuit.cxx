#include "NetlistSubcircuit.hxx"
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include "debuggingToolsForCpp.hxx"
#include <cfloat>

void NetlistSubcircuit::initialiseSubcircuit()
{
  numberOfPrescribedPressuresAndFlows = m_circuitData->numberOfPrescribedPressures + m_circuitData->numberOfPrescribedFlows; // Just the sum of the previous two declared integers

  // Resize to contain the necessary flows and pressures, and initialise to zero:
  flowsInSubcircuit.resize(m_circuitData->numberOfComponents,0.0);
  pressuresInSubcircuit.resize(m_circuitData->numberOfPressureNodes,0.0);

  getMapOfPressHistoriesToCorrectPressNodes(); //initialises numberOfHistoryPressures
  getMapOfFlowHistoriesToCorrectComponents(); //initialises numberOfHistoryFlows
  getMapOfVolumeHistoriesToCorrectComponents(); // initialises numberOfHistoryVolumes
  getMapOfTrackedVolumesToCorrectComponents(); // initialises m_numberOfTrackedVolumes

  volumesInSubcircuit.resize(m_numberOfTrackedVolumes,0.0);

  systemSize = m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->numberOfComponents + numberOfHistoryFlows + m_numberOfTrackedVolumes + numberOfHistoryVolumes;

  createVectorsAndMatricesForCircuitLinearSystem();

  // columnMapSize = numberOfHistoryPressures + numberOfHistoryFlows + numberOfPrescribedPressures + numberOfPrescribedFlows;
  getListOfNodesWithMultipleIncidentCurrents();
}

void NetlistSubcircuit::createVectorsAndMatricesForCircuitLinearSystem()
{
  PetscErrorCode errFlag;
  // Create m_systemMatrix and m_inverseOfSystemMatrix (to be filled later)
  errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatZeroEntries(m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

  errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&m_inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatZeroEntries(m_inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  // To compute m_inverseOfSystemMatrix, we use Petsc to solve AX=B, with A,X and B all systemSize x systemSize matrices
  // B will just be m_identityMatrixForPetscInversionHack.
  errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatZeroEntries(m_identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  // Fill the diagonal with ones:
  for (int ii=0; ii<systemSize; ii++)
  {
      errFlag = MatSetValue(m_identityMatrixForPetscInversionHack,ii,ii,1.0,INSERT_VALUES);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  }
  errFlag = MatAssemblyBegin(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = MatAssemblyEnd(m_identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);

  errFlag = VecCreate(PETSC_COMM_SELF,&RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecSetType(RHS,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make RHS a VECSEQ sequential vector
  errFlag = VecSetSizes(RHS,systemSize,systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecZeroEntries(RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyBegin(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyEnd(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

  errFlag = VecCreate(PETSC_COMM_SELF,&solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecSetType(solutionVector,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make solutionVector a VECSEQ sequential vector
  errFlag = VecSetSizes(solutionVector,systemSize,systemSize); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecZeroEntries(solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyBegin(solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  errFlag = VecAssemblyEnd(solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistSubcircuit::getListOfNodesWithMultipleIncidentCurrents()
{
    // Note that this function also counts pressure nodes which are just
    // /between/ two components, eg. for the (two resistor) subcircuit:
    //        N0----[==R1==]----N1----[==R2==]----N2
    // This would count node N1 as appearing twice, and do a "Kirchoff" current
    // balance of the form "flow through R1 = flow through R2".
    //
    // It also catches and deals with true Kirchoff equations wherea third (fourth, fifth,...)
    // component is connected to N1.

    int numberOfTimesNodeSeen;

    numberOfMultipleIncidentCurrentNodes = int(0);

    // The node data from the input file is 1-indexed, so shift this to 1:numberOfPressureNodes, instead of 0:numberOfPressureNodes-1
    for(auto node=m_circuitData->mapOfPressureNodes.begin(); node!=m_circuitData->mapOfPressureNodes.end(); node++)
    {
       int nodeIndex = node->second->getIndex();
       numberOfTimesNodeSeen = int(0);
       for (int ii = 0; ii<m_circuitData->numberOfComponents; ii++)
       {
          if ((m_circuitData->components.at(ii)->startNode->getIndex() == nodeIndex) || 
          	  (m_circuitData->components.at(ii)->endNode->getIndex() == nodeIndex))
          {
             numberOfTimesNodeSeen++;
          }
       }
       if (numberOfTimesNodeSeen > int(1))
       {
          listOfNodesWithMultipleIncidentCurrents.push_back(nodeIndex);
          numberOfMultipleIncidentCurrentNodes++;
       }
    }

}

void NetlistSubcircuit::getMapOfPressHistoriesToCorrectPressNodes()
{
    for (int ii=0; ii<m_circuitData->numberOfComponents; ii++)
    {
       // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
       if (m_circuitData->components.at(ii)->getType() == Component_Capacitor)
       {
       		listOfHistoryPressures.insert(m_circuitData->components.at(ii)->startNode->getIndex());
          m_circuitData->components.at(ii)->startNode->hasHistoryPressure = true;
       		listOfHistoryPressures.insert(m_circuitData->components.at(ii)->endNode->getIndex());
          m_circuitData->components.at(ii)->endNode->hasHistoryPressure = true;
       }
    }

    numberOfHistoryPressures = listOfHistoryPressures.size();

    // Now do the actual generation of the pressure history node ordering map (for use when bulding the linear system matrix):
    int ii=0;
    for (auto iterator=listOfHistoryPressures.begin(); iterator != listOfHistoryPressures.end(); iterator++)
    {
       nodeIndexToPressureHistoryNodeOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
       ii++;
    }
}


void NetlistSubcircuit::getMapOfFlowHistoriesToCorrectComponents()
{

	for (int ii=0; ii<m_circuitData->numberOfComponents; ii++)
	{
	   // Check for inductor, as these need flow "histories" (flow from the previous time-step) (for dQ/dt term).
     if(m_circuitData->components.at(ii)->getType() == Component_Inductor)
	   {
	   		listOfHistoryFlows.insert(ii);
        m_circuitData->components.at(ii)->hasHistoryFlow = true;
	   }
	}

	numberOfHistoryFlows = listOfHistoryFlows.size();

	// Now do the actual generation of the flow history component ordering map (for use when bulding the linear system matrix):
	int ii=0;
	for (auto iterator=listOfHistoryFlows.begin(); iterator != listOfHistoryFlows.end(); iterator++)
	{
	   componentIndexToFlowHistoryComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
	}
}

void NetlistSubcircuit::getMapOfVolumeHistoriesToCorrectComponents()
{

  for (int ii=0; ii<m_circuitData->numberOfComponents; ii++)
  {
     // Check for VolumeTrackingPressureChambers, as these need volume "histories" (volume from the previous time-step) (for dVolume/dt term).
     if(m_circuitData->components.at(ii)->getType() == Component_VolumeTrackingPressureChamber)
     {
        listOfHistoryVolumes.insert(ii);
        m_circuitData->components.at(ii)->hasHistoryVolume = true;
     }
  }

  numberOfHistoryVolumes = listOfHistoryVolumes.size();

  // Now do the actual generation of the pressure history node ordering map (for use when bulding the linear system matrix):
  int ii=0;
  for (auto iterator=listOfHistoryVolumes.begin(); iterator != listOfHistoryVolumes.end(); iterator++)
  {
     componentIndexToVolumeHistoryComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
  }
}

void NetlistSubcircuit::getMapOfTrackedVolumesToCorrectComponents()
{

  for (int ii=0; ii<m_circuitData->numberOfComponents; ii++)
  {
     // Check for VolumeTrackingPressureChambers, as these need volume tracking (for computing the current pressure, via the compliance/elastance).
     if(m_circuitData->components.at(ii)->getType() == Component_VolumeTrackingPressureChamber)
     {
        listOfTrackedVolumes.insert(ii);
        m_circuitData->components.at(ii)->hasTrackedVolume = true;
     }
  }

  m_numberOfTrackedVolumes = listOfTrackedVolumes.size();

  // Now do the actual generation of the volume node ordering map (for use when bulding the linear system matrix):
  int ii=0;
  for (auto iterator=listOfTrackedVolumes.begin(); iterator != listOfTrackedVolumes.end(); iterator++)
  {
     componentIndexToTrackedVolumeComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
     ii++;
  }
}

void NetlistSubcircuit::generateLinearSystemFromPrescribedCircuit(const double alfi_delt)
{
    // This function assembles the system of (time-discretised) linear algebraic equations for the LPN.
    PetscErrorCode errFlag;

    errFlag = MatZeroEntries(m_systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    {
      int row = 0; // is the row in the matrix that we write to on each occasion
      for (auto component=m_circuitData->components.begin(); component!=m_circuitData->components.end(); component++)
      {
          // bool componentIsOpenDiode = (m_circuitData->components.at(ll)->type == Component_Diode &&
          //                              m_circuitData->components.at(ll)->hasNonnegativePressureGradientAndNoBackflow());
          // open diodes are just implemented as zero-resistance resistors, closed diodes are zero-resistance resistors with prescribed zero flow
          if ((*component)->getType() == Component_Resistor || (*component)->getType() == Component_Diode)
          {
            // insert resistor(-eqsue) relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures,-currentParameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            row++;
          }
          else if ((*component)->getType() == Component_Capacitor)
          {
            // insert capacitor relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures,-alfi_delt/currentParameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatSetValue(m_systemMatrix,row,nodeIndexToPressureHistoryNodeOrderingMap.at(startNode)+m_circuitData->numberOfPressureNodes,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatSetValue(m_systemMatrix,row,nodeIndexToPressureHistoryNodeOrderingMap.at(endNode)+m_circuitData->numberOfPressureNodes,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            row++;
          }
          else if ((*component)->getType() == Component_Inductor)
          {
            // insert inductor relationship into equation system
            int startNode = (*component)->startNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            
            int endNode = (*component)->endNode->getIndex();
            errFlag = MatSetValue(m_systemMatrix,row,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);

            double currentParameterValue = *((*component)->getParameterPointer());
            int indexOfThisComponentsFlow = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,indexOfThisComponentsFlow+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures,-currentParameterValue/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            errFlag = MatSetValue(m_systemMatrix,row,componentIndexToFlowHistoryComponentOrderingMap.at(indexOfThisComponentsFlow)+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures+m_circuitData->numberOfComponents,currentParameterValue/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            row++;
          }
          else if ((*component)->getType() == Component_VolumeTrackingPressureChamber)
          {
            // Two equations are needed for the VolumeTrackingPressureChamber:
            // 1) dVolume/dt = flow
            // 2) compliance * pressure = (volume - unstressed volume) ... unstressed vol will go on RHS of system, later.
            // Note that this means we increment row (row++) twice during this if-case
            // Do (1):
            // Insert the dt*flow term:
            int zeroIndexOfThisComponent = toZeroIndexing((*component)->getIndex());
            errFlag = MatSetValue(m_systemMatrix,row,zeroIndexOfThisComponent+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures,-alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            // Insert the volume term (volume to-be-found during the next system solve):
            {
              int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            // Insert the volume history term:
            {
              int columnIndex = componentIndexToVolumeHistoryComponentOrderingMap.at(zeroIndexOfThisComponent) + m_numberOfTrackedVolumes + m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->numberOfComponents + numberOfHistoryFlows;
              // std::cout << "setting in NetlistSubcircuit.cxx row and column: " << row << " " << columnIndex << std::endl;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
            }
            row++; // done twice in this if-case, because there are 2 equations to create for the VolumeTrackingPressureChamber

            // Now do (2) (see comment block above, within this if-case)
            // Do the compliance term:
            boost::shared_ptr<VolumeTrackingPressureChamber> volumeTrackingPressureChamber = boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (*component);
            if (!volumeTrackingPressureChamber->zeroVolumeShouldBePrescribed()) // test whether, on a previous attempt to solve the system, negative volumes were detected. If so, we'll do something different in the "else" below...
            {
              {
                int columnIndex = toZeroIndexing((*component)->startNode->getIndex());
                double valueToInsert = 1.0/(volumeTrackingPressureChamber->getElastance());
                errFlag = MatSetValue(m_systemMatrix,row,columnIndex,valueToInsert,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              }
              // Do the volume term:
              {
                int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->numberOfComponents + numberOfHistoryFlows;
                errFlag = MatSetValue(m_systemMatrix,row,columnIndex,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              }
            }
            else //volumeTrackingPressureChamber->zeroVolumeShouldBePrescribed() == true, so instead of prescribing pressure based on volume, we allow the pressure to be a free variable, and prescribe zero volume for the chamber
            {
              // Prescribe zero volume:
              // VERY IMPORTANT: Note that we don't need to prescribe anything special on the RHS for this; the
              // prescribed volume is zero, we're doing a row for a component (which have zeros on the RHS anyway), and the
              // RHS is zeroed out before we start to fill it, so the zero will already be in place. But be aware of this
              // if you're making changes.
              int columnIndex = componentIndexToTrackedVolumeComponentOrderingMap.at(zeroIndexOfThisComponent) + m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->numberOfComponents + numberOfHistoryFlows;
              errFlag = MatSetValue(m_systemMatrix,row,columnIndex,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // Reset the zero-volume marker on the component:
              volumeTrackingPressureChamber->resetZeroVolumePrescription();
            }
            row++; // done twice in this if-case, because there are 2 equations to create for the VolumeTrackingPressureChamber

          }
          else
          {
            throw std::runtime_error("EE: Unknown component type in netlist. Halting.");
          }
      }
    }

     // Do the equations for the nodes with multiple incident currents:
     for (int mm=0; mm<numberOfMultipleIncidentCurrentNodes; mm++)
     {
       for (int ll=0; ll<m_circuitData->numberOfComponents; ll++)
       {
          bool foundMultipleIncidentCurrentsForEndNode = (m_circuitData->components.at(ll)->endNode->getIndex() == listOfNodesWithMultipleIncidentCurrents.at(mm)); 
          if (foundMultipleIncidentCurrentsForEndNode)
          {
            int row = mm + m_circuitData->numberOfComponents + m_numberOfTrackedVolumes;
            errFlag = MatSetValue(m_systemMatrix,row,ll+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }

          bool foundMultipleIncidentCurrentsForStartNode = (m_circuitData->components.at(ll)->startNode->getIndex() == listOfNodesWithMultipleIncidentCurrents.at(mm));
          if (foundMultipleIncidentCurrentsForStartNode)
          {
            int row = mm + m_circuitData->numberOfComponents + m_numberOfTrackedVolumes;
            errFlag = MatSetValue(m_systemMatrix,row,ll+m_circuitData->numberOfPressureNodes+numberOfHistoryPressures,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
       }
     }

     int rowsDoneSoFar = m_circuitData->numberOfComponents + numberOfMultipleIncidentCurrentNodes + m_numberOfTrackedVolumes;

     // create the columnMap which tells us which system column each of the prescribed pressure, pressure-history, flow, or volume values belong to
     columnMap.clear();
     int tempUnknownVariableIndexWithinLinearSystem = 0; // just an indexing shift to keep track of where we need to write next
     // Do the prescribed pressures:
     for (auto prescribedPressureNode=m_circuitData->mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=m_circuitData->mapOfPrescribedPressureNodes.end(); prescribedPressureNode++)
     {
     	columnMap.push_back(toZeroIndexing(prescribedPressureNode->second->getIndex()));
     }

     // Do the history pressures
     tempUnknownVariableIndexWithinLinearSystem += m_circuitData->numberOfPressureNodes; // tempUnknownVariableIndexWithinLinearSystem is zero before this line; I'm doing it like this for clarity & consistency
     for (int ll=0; ll<numberOfHistoryPressures; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     // Do the prescribed flows
     tempUnknownVariableIndexWithinLinearSystem += numberOfHistoryPressures;
     for (auto prescribedFlowComponent=m_circuitData->mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=m_circuitData->mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
     {
     	columnMap.push_back(toZeroIndexing(prescribedFlowComponent->second->getIndex())+tempUnknownVariableIndexWithinLinearSystem);
     }

     // Do the history flows
     tempUnknownVariableIndexWithinLinearSystem += m_circuitData->numberOfComponents;
     for (int ll=0; ll<numberOfHistoryFlows; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     // Do the history volumes
     tempUnknownVariableIndexWithinLinearSystem += numberOfHistoryFlows + m_numberOfTrackedVolumes;
     for (int ll=0; ll<numberOfHistoryVolumes; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     // Set the prescribed-value equations (i.e. pressure_1 (LHS) = pressure_1 (RHS) - so really just a way of setting the prescribed values within the linear system)
     for (int ll = 0; ll<systemSize - rowsDoneSoFar; ll++)
     {
       errFlag = MatSetValue(m_systemMatrix,rowsDoneSoFar + ll, columnMap.at(ll),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     }

     errFlag = MatAssemblyBegin(m_systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     errFlag = MatAssemblyEnd(m_systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    // std::cout << "System matrix for surface " << surfaceIndex << ":" << std::endl;
    //  errFlag = MatView(m_systemMatrix,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
     errFlag = MatLUFactor(m_systemMatrix,NULL,NULL,NULL);CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistSubcircuit::assembleRHS(const int timestepNumber)
{

    PetscErrorCode errFlag;
    errFlag = VecZeroEntries(RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    historyPressuresInSubcircuit = pressuresInSubcircuit;
    historyFlowsInSubcircuit = flowsInSubcircuit;
    historyVolumesInSubcircuit = volumesInSubcircuit;

    columnIndexOf3DInterfacePressureInLinearSystem.clear(); // dummy value, to be replaced!

    // int nextPressurePointerIndex = 0; // for tracking which pressure pointer to use - useful when there are multiple pressure interfaces to other domains / subcircuits

    // Prescribed pressures
    int tempIndexingShift = m_circuitData->numberOfComponents + numberOfMultipleIncidentCurrentNodes + m_numberOfTrackedVolumes;
    {
      int ll=0;
      for (auto prescribedPressureNode=m_circuitData->mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=m_circuitData->mapOfPrescribedPressureNodes.end(); prescribedPressureNode++ )
      {
        // Coming from 'f' for 'fixed' in the input data:
        if (prescribedPressureNode->second->prescribedPressureType == Pressure_Fixed)
        {
            errFlag = VecSetValue(RHS,ll + tempIndexingShift,prescribedPressureNode->second->getPressure(),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);	
        }
        // Coming from 'l' for 'left-ventricular' in the input data:
        else if (prescribedPressureNode->second->prescribedPressureType == Pressure_LeftVentricular)
        {
          	std::cerr << "this requires heart model. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
            std::exit(1);
        }
        else if (prescribedPressureNode->second->prescribedPressureType == Pressure_3DInterface)
        {
            // We only do this if the netlist is in Dirichlet BC mode:
            if (m_circuitData->hasPrescribedPressureAcrossInterface())
            {
              columnIndexOf3DInterfacePressureInLinearSystem.push_back(ll + tempIndexingShift);
              double* pressurePointerToSet = pressure_n_ptrs.at(prescribedPressureNode->second->prescribedPressurePointerIndex);
              errFlag = VecSetValue(RHS,ll + tempIndexingShift,*pressurePointerToSet,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // nextPressurePointerIndex++;
            }
        }
        else
        {
        	  throw std::runtime_error("Unknown pressure prescription value in Netlist.");
        }

        // // get the column index for the pressure in the linear system:
        // if (prescribedPressureNode->second->isAtBoundary())
        // {
        //   columnIndexOf3DInterfacePressureInLinearSystem.push_back(ll + tempIndexingShift);
        // }

        ll++;
      }
    }

    // for (int ll=0; ll<m_circuitData->numberOfPrescribedPressures; ll++)
    // {
    //    // 'f' for 'fixed'
    //    if (m_circuitData->typeOfPrescribedPressures.at(ll) == Pressure_Fixed)
    //    {
    //       errFlag = VecSetValue(RHS,ll + tempIndexingShift,valueOfPrescribedPressures.at(ll),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //    }
    //    // 'l' for 'leftVentricular'
    //    else if (subcircuitInputData.typeOfPrescribedPressures.at(ll) == Pressure_LeftVentricular)
    //    {
    //       std::cout << "this requires heartmodel. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
    //       std::exit(1);
    //       // if ((timestepNumber == 0) || (timestepNumber == 1)) // treat case with no known IM pressure yet
    //       // {
    //       //    P_IM_mid_lasttimestep = 5000; // \todo find a better way of doing this; maybe input this value from file...
    //       //    P_IM_mid = 5000; // ... or set it based on the aortic valve state at simulation start
    //       // }
    //       // elseif (timestepNumber .eq. int(2)) then // treat case where only one IM pressure history point is known
    //       //    P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber-1)
    //       //    P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber)
    //       // else // get the previous intramyocardial pressure in the case where we have enough doata for this (see comment before start of "if" block)
    //       //    P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber-1)
    //       //    P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber)
    //       // end if

    //       // errFlag = VecSetValue(RHS,ll + tempIndexingShift,P_IM_mid,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //       // int nn=0;
    //       // for (auto iterator=listOfHistoryPressures.begin(); iterator!=listOfHistoryPressures.end(); iterator++, nn++)
    //       // {
    //       //    if (*iterator == listOfPrescribedPressures.at(ll))
    //       //    {
    //       //       historyPressuresInSubcircuit.at(*iterator) = P_IM_mid_lasttimestep;
    //       //    }
    //       // }
    //    }
    //    else
    //    {
    //         throw std::runtime_error("Unknown pressure prescription value in Netlist.");
    //    }
    // }
    // History Pressures
    tempIndexingShift = tempIndexingShift + m_circuitData->numberOfPrescribedPressures;
    // for(int ll=0; ll<numberOfHistoryPressures; ll++)

    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
	    int lll=0;
	    for (auto node=m_circuitData->mapOfPressureNodes.begin(); node!=m_circuitData->mapOfPressureNodes.end(); node++)
	    {
	    	if (node->second->hasHistoryPressure)
	    	{
	    		errFlag = VecSetValue(RHS,lll+tempIndexingShift,node->second->historyPressure,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	    		lll++;
	    	}
	    }
	    // for (auto iterator=listOfHistoryPressures.begin(); iterator!=listOfHistoryPressures.end(); iterator++)
	    // {
	    //     errFlag = VecSetValue(RHS,lll+tempIndexingShift,historyPressuresInSubcircuit.at(*iterator - 1),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	    // }
	    // lll++;
	  }
    // Prescribed Flows:
    tempIndexingShift += numberOfHistoryPressures;
    columnIndexOf3DInterfaceFlowInLinearSystem.clear();
    // int nextFlowPointerIndex = 0;
    // Scoping unit to include the second counter ll in the for loop, without having ll in-scope after the loop finishes:
    {
    	int ll=0;
	    for (auto prescribedFlowComponent=m_circuitData->mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=m_circuitData->mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
	    {
	       if (prescribedFlowComponent->second->prescribedFlowType == Flow_3DInterface)
	       {
            if (m_circuitData->hasPrescribedFlowAcrossInterface())
            {
  	          columnIndexOf3DInterfaceFlowInLinearSystem.push_back(ll + tempIndexingShift);
              double* flowPointerToSet = flow_n_ptrs.at(prescribedFlowComponent->second->prescribedFlowPointerIndex);
              // First, flip the sign of the flow, if necessary due to the orientation of the component at the 3D interface:
              double threeDFlowValue = *flowPointerToSet * prescribedFlowComponent->second->m_signForPrescribed3DInterfaceFlow;
              assert(!isnan(threeDFlowValue));
              // Give the (possibly sign-corrected) flow to the linear system:
  	          errFlag = VecSetValue(RHS,ll + tempIndexingShift,threeDFlowValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
              // nextFlowPointerIndex++;
            }
	       }
	       else if (prescribedFlowComponent->second->prescribedFlowType == Flow_Fixed)
	       {
	          errFlag = VecSetValue(RHS,ll + tempIndexingShift,prescribedFlowComponent->second->valueOfPrescribedFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	       }
         // else if (prescribedFlowComponent->second->prescribedFlowType == Flow_Diode_FixedWhenClosed)
         // {
         //    // Check whether the diode is closed; if so, prescribe the (zero) flow:
         //    if (!prescribedFlowComponent->second->hasNonnegativePressureGradientAndNoBackflow())
         //    {
         //      errFlag = VecSetValue(RHS,ll + tempIndexingShift,prescribedFlowComponent->second->valueOfPrescribedFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
         //    }
         // }
	       else
	       {
	            throw std::runtime_error("Unknown flow prescription value in Netlist.");
	       }

         // if (prescribedFlowComponent->second->connectsToNodeAtInterface())
         // {
         //    columnIndexOf3DInterfaceFlowInLinearSystem.push_back(ll + tempIndexingShift);
         // }

	       ll++;
	    }
	  }
    // for(int ll=0; ll<numberOfPrescribedFlows; ll++)
    // {
    //    if (subcircuitInputData.typeOfPrescribedFlows.at(ll) == Flow_3DInterface)
    //    {
    //       columnIndexOf3DInterfaceFlowInLinearSystem = ll + tempIndexingShift;
    //       errFlag = VecSetValue(RHS,ll + tempIndexingShift,*flow_n_ptr,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //    }
    //    else if (subcircuitInputData.typeOfPrescribedFlows.at(ll) == Flow_Fixed)
    //    {
    //       errFlag = VecSetValue(RHS,ll + tempIndexingShift,valueOfPrescribedFlows.at(ll), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    //    }
    //    else
    //    {
    //         throw std::runtime_error("Unknown flow prescription value in Netlist.");
    //    }
    // }
    // History Flows
    tempIndexingShift = tempIndexingShift + m_circuitData->numberOfPrescribedFlows;
    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
	    int lll=0;
	    for (auto component=m_circuitData->mapOfComponents.begin(); component!=m_circuitData->mapOfComponents.end(); component++)
	    {
	    	if (component->second->hasHistoryFlow)
	    	{
	    		errFlag = VecSetValue(RHS,lll + tempIndexingShift,component->second->historyFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	    		lll++;		
	    	}
	    }
	    // for (auto iterator=listOfHistoryFlows.begin(); iterator!=listOfHistoryFlows.end(); iterator++)
	    // {
	    //    errFlag = VecSetValue(RHS,lll + tempIndexingShift,historyFlowsInSubcircuit.at(*iterator - 1), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	    // }
	    // lll++;
	  }

    // Give the linear system the history volumes, by putting them on the RHS:
    tempIndexingShift += numberOfHistoryFlows;
    {
      int ll=0;
      for (auto component=m_circuitData->components.begin(); component!=m_circuitData->components.end(); component++)
      {
        if ((*component)->hasHistoryVolume)
        {
          // currently, only VolumeTrackingPresureChambers have history volumes. We might want to change this cast later, if new component types
          // with history volumes get added.
          boost::shared_ptr<VolumeTrackingPressureChamber> volumeTrackingPressureChamber = boost::dynamic_pointer_cast<VolumeTrackingPressureChamber> (*component);
          double volume = volumeTrackingPressureChamber->getHistoryVolume();
          int row = ll + tempIndexingShift;
          errFlag = VecSetValue(RHS,row,volume,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          ll++;
        }
      }
    }
 
    errFlag = VecAssemblyBegin(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // std::cout << "RHS for surface " << surfaceIndex << ":" << std::endl;
    // errFlag = VecView(RHS,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // std::cout << "END RHS for surface " << surfaceIndex << std::endl;

}

void NetlistSubcircuit::updateInternalPressuresVolumesAndFlows(const int timestepNumber)
{
    computeCircuitLinearSystemSolution(timestepNumber);

    // Get the updated nodal pressures:
    giveNodesTheirPressuresFromSolutionVector();

    // Get the updated component flows:
    giveComponentsTheirFlowsFromSolutionVector();

    // Get the updated volumes:
    giveComponentsTheirVolumesFromSolutionVector();

    // for (int volumeIndex = indexShift; volumeIndex < indexShift + m_numberOfTrackedVolumes; volumeIndex++)
    // {
    //     errFlag = VecGetValues(solutionVector,getSingleValue,&volumeIndex,&volumesInSubcircuit[volumeIndex-indexShift]); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    //     VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (m_circuitData->mapOfVolumeTrackingComponents.at(toOneIndexing(volumeIndex-indexShift)).get());
    //     currentPressureChamber->setStoredVolume(volumesInSubcircuit[volumeIndex-indexShift]);
    //     currentPressureChamber->passPressureToStartNode();
    // }

    // // Update the volumes in each VolumeTrackingPressureChamber
    // for (auto component=m_circuitData->mapOfComponents.begin(); component!=m_circuitData->mapOfComponents.end(); component++)
    // {
    //   // detect the VolumeTrackingPressureChambers:
    //   if (component->second->type == Component_VolumeTrackingPressureChamber)
    //   {
    //     VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
    //     currentPressureChamber->updateStoredVolume(alfi_delt);
    //     currentPressureChamber->passPressureToStartNode();
    //   }
    // }

    // write(*,*) 'discrepancy:', (-this%P_a(1) - this%pressuresInSubcircuit(2))/1.2862d5 - this%flowsInSubcircuit(1)
}

void NetlistSubcircuit::computeCircuitLinearSystemSolution(const int timestepNumber)
{
  buildAndSolveLinearSystem(timestepNumber);
}

void NetlistSubcircuit::giveNodesTheirPressuresFromSolutionVector()
{
  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;

  // Look the nodes, handing them their new pressures from the circuit linear system solve:
  for (int ll=0; ll<m_circuitData->numberOfPressureNodes; ll++)
  {
      errFlag = VecGetValues(solutionVector,getSingleValue,&ll,&pressuresInSubcircuit[ll]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      // std::cout << "system matrix: " << surfaceIndex << std::endl;
      // MatView(m_systemMatrix,PETSC_VIEWER_STDOUT_WORLD);
      // std::cout << "RHS: " << std::endl;
      // VecView(RHS,PETSC_VIEWER_STDOUT_WORLD);
      // std::cout << "solution vec: " << std::endl;
      // VecView(solutionVector,PETSC_VIEWER_STDOUT_WORLD);
      assert(!isnan(pressuresInSubcircuit[ll]));
      m_circuitData->mapOfPressureNodes.at(toOneIndexing(ll))->setPressure(pressuresInSubcircuit[ll]);
  }
}

void NetlistSubcircuit::giveComponentsTheirFlowsFromSolutionVector()
{
  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;

  int firstFlowIndex = m_circuitData->numberOfPressureNodes + numberOfHistoryPressures;
  for (int ll=firstFlowIndex; ll<m_circuitData->numberOfComponents+firstFlowIndex; ll++)
  {
      errFlag = VecGetValues(solutionVector,getSingleValue,&ll,&flowsInSubcircuit[ll-firstFlowIndex]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
      assert(!isnan(flowsInSubcircuit[ll-firstFlowIndex]));
      m_circuitData->mapOfComponents.at(toOneIndexing(ll-firstFlowIndex))->flow = flowsInSubcircuit[ll-firstFlowIndex];
  }
}

void NetlistSubcircuit::giveComponentsTheirVolumesFromSolutionVector()
{
  std::vector<double> volumes = getVolumesFromSolutionVector();
  // Reverse the volumes so we can pop off the back of it as we loop the mapOfVolumeTrackingComponents:
  std::reverse(volumes.begin(), volumes.end());
  for (auto component = m_circuitData->mapOfVolumeTrackingComponents.begin(); component != m_circuitData->mapOfVolumeTrackingComponents.end(); component++)
  {
      VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
      
      // Ensure we aren't dealing with negative volumes:
      //\todo REINSTATE!
      // assert(volumes.back() >= 0.0);
      assert(!isnan(volumes.back()));

      currentPressureChamber->setStoredVolume(volumes.back());
      volumes.pop_back();
      currentPressureChamber->passPressureToStartNode();
  }

  // Ensure we've used all the retrieved volumes:
  assert(volumes.size() == 0);
}

// Compare with giveComponentsTheirVolumesFromSolutionVector. That function sets final, accepted volumes,
// whereas this function, giveComponentsTheirProposedVolumesFromSolutionVector, gives them /proposed/ volumes
// which are then checked for validity (i.e. non-negativity). Any proposed negative values trigger a re-computation
// of the solution to the circuit linear system, with an enforced zero-volume at the would-be negative volume locations.
void NetlistSubcircuit::giveComponentsTheirProposedVolumesFromSolutionVector()
{
  std::vector<double> volumes = getVolumesFromSolutionVector();
  // Reverse the volumes so we can pop off the back of it as we loop the mapOfVolumeTrackingComponents:
  std::reverse(volumes.begin(), volumes.end());
  for (auto component = m_circuitData->mapOfVolumeTrackingComponents.begin(); component != m_circuitData->mapOfVolumeTrackingComponents.end(); component++)
  {
      VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
      
      currentPressureChamber->setProposedVolume(volumes.back());
      volumes.pop_back();
  }

  // Ensure we've used all the retrieved volumes:
  assert(volumes.size() == 0);
}

std::vector<double> NetlistSubcircuit::getVolumesFromSolutionVector()
{
  std::vector<double> volumesToReturn;

  PetscErrorCode errFlag;

  // A self-documenting name for the request given to VecGetValues():
  int getSingleValue=1;

  // The location of the first volume in the solutionVector:
  int firstVolumeIndex = m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->numberOfComponents + numberOfHistoryFlows;
  int volumeIndex = firstVolumeIndex;

  auto component = m_circuitData->mapOfVolumeTrackingComponents.begin();
  while(component != m_circuitData->mapOfVolumeTrackingComponents.end())
  {
    errFlag = VecGetValues(solutionVector,getSingleValue,&volumeIndex,&volumesInSubcircuit[volumeIndex-firstVolumeIndex]); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    VolumeTrackingPressureChamber* currentPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
    volumesToReturn.push_back(volumesInSubcircuit[volumeIndex-firstVolumeIndex]);

    volumeIndex++;
    component++;
  }

  return volumesToReturn;
}

void NetlistSubcircuit::buildAndSolveLinearSystem(const int timestepNumber)
{
  generateLinearSystemFromPrescribedCircuit(m_alfi_delt);
  assembleRHS(timestepNumber);

  PetscErrorCode errFlag;
  // get the inverse of the system matrix:
  errFlag = MatMatSolve(m_systemMatrix,m_identityMatrixForPetscInversionHack,m_inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
  // Release the m_systemMatrix so we can edit it again on the next iteration (we only need the just-computed m_inverseOfSystemMatrix for computations on this step now.)
  errFlag = MatSetUnfactored(m_systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);

  // Solve the system
  errFlag = MatMult(m_inverseOfSystemMatrix,RHS,solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

std::pair<boundary_data_t,double> NetlistSubcircuit::computeAndGetFlowOrPressureToGiveToZeroDDomainReplacement(const int timestepNumber)
{
  buildAndSolveLinearSystem(timestepNumber);

  PetscErrorCode errFlag;

  int numberOfValuesToGet=1;
  if (m_circuitData->hasPrescribedFlowAcrossInterface()) // This boundary condition is receiving flow and returning pressure (Neumann mode if NetlistSubcircuit is a boundary condition)
  {
    int locationOfPressureInRHS = toZeroIndexing(m_circuitData->getIndexOfNodeAtInterface());
    errFlag = VecGetValues(solutionVector,numberOfValuesToGet,&locationOfPressureInRHS,&m_interfacePressure); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    return std::make_pair(Boundary_Pressure,m_interfacePressure);
  }
  else if (m_circuitData->hasPrescribedPressureAcrossInterface()) // This boundary condition is receiving pressure and returning flow (Dirichlet mode if NetlistSubcircuit is a boundary condition)
  {
    int locationOfFlowInRHS = toOneIndexing(m_circuitData->numberOfPressureNodes + numberOfHistoryPressures + m_circuitData->getIndexOfComponentConnectingToNodeAtInterface());

    errFlag = VecGetValues(solutionVector,numberOfValuesToGet,&locationOfFlowInRHS,&m_interfaceFlow);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    return std::make_pair(Boundary_Flow,m_interfaceFlow);
  }
  else
  {
    std::stringstream errorMessage;
    errorMessage << "EE: Internal error when computing flow or pressure to pass to zero-D domain replacement." << std::endl;
    throw std::logic_error(errorMessage.str());
  }
}

std::pair<double,double> NetlistSubcircuit::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
    assert(m_circuitData->connectsTo3DDomain());

    PetscErrorCode errFlag;
    {
      int safetyCounter = 0;
      bool solutionVectorMightHaveNegativeVolumes = true;
      // This loop keeps repeating the circuit linear system solve, until we have
      // ensured there are no negative volumes:
      while (solutionVectorMightHaveNegativeVolumes)
      {
        buildAndSolveLinearSystem(timestepNumber);

        solutionVectorMightHaveNegativeVolumes = areThereNegativeVolumes(timestepNumber);

        safetyCounter++;
        if (safetyCounter > 1 && safetyCounter < 5)
        {
          std::cout << "II: Redoing due to a detected zero-volume problem! ----------------------------------------------" << std::endl;
        }
        if (safetyCounter > safetyCounterLimit)
        {
          std::stringstream errorMessage;
          errorMessage << "EE: Took too long (" << safetyCounter << " repeated solves of the circuit linear system) to eradicate negative volumes at the domain boundary with index " << surfaceIndex << "." << std::endl;
          errorMessage << "This was probably caused by a bad (or an extremely large) Netlist circuit!" << std::endl;
          throw std::runtime_error(errorMessage.str());
        }
      }
    }

    // Extract the implicit coeffcients, for eventual passing to the FORTRAN
    // linear solve
    int rowToGet[] = {0};
    int numberOfValuesToGet=1;
    PetscScalar valueFromInverseOfSystemMatrix;
    std::pair<double,double> returnValue;
    

    PetscScalar valueFromRHS;
    // if (columnIndexOf3DInterfaceFlowInLinearSystem.size()!=0) //\todo remove this hack and have the mat/vecgetvalues unguarded here!!
    // {
      errFlag = MatGetValues(m_inverseOfSystemMatrix,numberOfValuesToGet,rowToGet,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem.at(0),&valueFromInverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
      returnValue.first = valueFromInverseOfSystemMatrix;

      errFlag = VecGetValues(RHS,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem.at(0),&valueFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
      PetscScalar valueFromSolutionVector;
      errFlag = VecGetValues(solutionVector,numberOfValuesToGet,rowToGet,&valueFromSolutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
      returnValue.second = valueFromSolutionVector - valueFromInverseOfSystemMatrix * valueFromRHS;//\todo make dynamic
    // }
    // else
    // {
    //   returnValue.first = DBL_MAX;
    //   returnValue.second = 0.0;
    // }
    
    // std::cout << "m_inverseOfSystemMatrix: "<< std::endl;
    // MatView(m_inverseOfSystemMatrix,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "solution vector: "<< std::endl;
    // VecView(solutionVector,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "RHS vector: "<< std::endl;
    // VecView(RHS,PETSC_VIEWER_STDOUT_WORLD);
    // std::cout << "and just set " << returnValue.first << " " <<returnValue.second << std::endl;


    return returnValue;
}

// This subroutine detects whether the last circuit linear system solve was invalid due to its producing negative
// volumes. The returned bool can be used to enforce a re-solve, with any negative pressures re-prescribed to be zero.
bool NetlistSubcircuit::areThereNegativeVolumes(const int timestepNumber)
{
  computeCircuitLinearSystemSolution(timestepNumber);
  // These volumes are "proposed", because if any are negative, we 
  // re-solve with zero-volume prescribed
  giveComponentsTheirProposedVolumesFromSolutionVector();

  bool thereAreNegativeVolumes = false;
  //\todo REINSTATE THIS LOOP!
  // for (auto component = m_circuitData->mapOfVolumeTrackingComponents.begin(); component!=m_circuitData->mapOfVolumeTrackingComponents.end(); component++)
  // {
  //   VolumeTrackingPressureChamber* volumeTrackingPressureChamber = dynamic_cast<VolumeTrackingPressureChamber*> (component->second.get());
  //   // Check for negative volumes:
  //   if (volumeTrackingPressureChamber->getProposedVolume() < 0.0)
  //   {
  //     thereAreNegativeVolumes = true;
  //     // Note that we're going to have to re-compute the circuit linear system solution, but this time
  //     // with the volume at this location enforced to zero (in order to avoid the negativity)
  //     volumeTrackingPressureChamber->enforceZeroVolumePrescription();
  //   }
  // }

  return thereAreNegativeVolumes;
}
