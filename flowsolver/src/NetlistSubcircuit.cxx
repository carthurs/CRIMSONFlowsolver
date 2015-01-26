#include "NetlistSubcircuit.hxx"
#include <cassert>
#include <stdexcept>
#include <iostream>

void NetlistSubcircuit::initialiseSubcircuit()
{
    numberOfPrescribedPressuresAndFlows = m_circuitData.numberOfPrescribedPressures + m_circuitData.numberOfPrescribedFlows; // Just the sum of the previous two declared integers

    // Resize to contain the necessary flows and pressures, and initialise to zero:
	flowsInSubcircuit.resize(m_circuitData.numberOfComponents,0.0);
	pressuresInSubcircuit.resize(m_circuitData.numberOfPressureNodes,0.0);

	getMapOfPressHistoriesToCorrectPressNodes(); //initialises numberOfHistoryPressures
	getMapOfFlowHistoriesToCorrectComponents(); //initialises numberOfHistoryFlows


	systemSize = m_circuitData.numberOfPressureNodes + numberOfHistoryPressures + m_circuitData.numberOfComponents + numberOfHistoryFlows;

    PetscErrorCode errFlag;
    // Create systemMatrix and inverseOfSystemMatrix (to be filled later)
	errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

	errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // To compute inverseOfSystemMatrix, we use Petsc to solve AX=B, with A,X and B all systemSize x systemSize matrices
    // B will just be identityMatrixForPetscInversionHack.
    errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatZeroEntries(identityMatrixForPetscInversionHack);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // Fill the diagonal with ones:
    for (int ii=0; ii<systemSize; ii++)
    {
        errFlag = MatSetValue(identityMatrixForPetscInversionHack,ii,ii,1.0,INSERT_VALUES);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    errFlag = MatAssemblyBegin(identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatAssemblyEnd(identityMatrixForPetscInversionHack,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);

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

    // columnMapSize = numberOfHistoryPressures + numberOfHistoryFlows + numberOfPrescribedPressures + numberOfPrescribedFlows;
    getListOfNodesWithMultipleIncidentCurrents();
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
    for(auto node=m_circuitData.mapOfPressureNodes.begin(); node!=m_circuitData.mapOfPressureNodes.end(); node++)
    {
       int nodeIndex = node->second->indexInInputData;
       numberOfTimesNodeSeen = int(0);
       for (int ii = 0; ii<m_circuitData.numberOfComponents; ii++)
       {
          if ((m_circuitData.components.at(ii)->startNode->indexInInputData == nodeIndex) || 
          	  (m_circuitData.components.at(ii)->endNode->indexInInputData == nodeIndex))
          {
             numberOfTimesNodeSeen++;
          }
       }
       if (numberOfTimesNodeSeen > int(1))
       {
          // this acts as a flag to the next loop, which will make the final listOfNodesWithMultipleIncidentCurrents.
          // listOfNodesWithMultipleIncidentCurrents_temp(node,kk) = int(1)
          listOfNodesWithMultipleIncidentCurrents.push_back(nodeIndex);
          numberOfMultipleIncidentCurrentNodes++;
       }
    }

}

void NetlistSubcircuit::getMapOfPressHistoriesToCorrectPressNodes()
{
    for (int ii=0; ii<m_circuitData.numberOfComponents; ii++)
    {
       // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
       if (m_circuitData.components.at(ii)->type == Component_Capacitor)
       {
       		listOfHistoryPressures.insert(m_circuitData.components.at(ii)->startNode->indexInInputData);
       		listOfHistoryPressures.insert(m_circuitData.components.at(ii)->endNode->indexInInputData);
       }
    }

    numberOfHistoryPressures = listOfHistoryPressures.size();

    // Now do the actual generation of the pressure history node ordering map:
    int ii=0;
    for (auto iterator=listOfHistoryPressures.begin(); iterator != listOfHistoryPressures.end(); iterator++, ii++)
    {
       nodeIndexToPressureHistoryNodeOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
    }
}


void NetlistSubcircuit::getMapOfFlowHistoriesToCorrectComponents()
{

	for (int ii=0; ii<m_circuitData.numberOfComponents; ii++)
	{
	   // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
        if(m_circuitData.components.at(ii)->type == Component_Inductor)
	   {
	   		listOfHistoryFlows.insert(ii);
	   }
	}

	numberOfHistoryFlows = listOfHistoryFlows.size();

	// Now do the actual generation of the pressure history node ordering map:
	int ii=0;
	for (auto iterator=listOfHistoryFlows.begin(); iterator != listOfHistoryFlows.end(); iterator++, ii++)
	{
	   componentIndexToFlowHistoryComponentOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
	}
}

void NetlistSubcircuit::generateLinearSystemFromPrescribedCircuit(const double alfi_delt)
{
    // This function assembles the system of (time-discretised) linear algebraic equations for the LPN.

    PetscErrorCode errFlag;

    errFlag = MatZeroEntries(systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    for(int ll=0; ll<m_circuitData.numberOfComponents; ll++)
    {
        if (m_circuitData.components.at(ll)->type == Component_Resistor)
        {
          // insert resistor relationship into equation system
          int startNode = m_circuitData.components.at(ll)->startNode->indexInInputData;
          errFlag = MatSetValue(systemMatrix,ll,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
          int endNode = m_circuitData.components.at(ll)->endNode->indexInInputData;
          errFlag = MatSetValue(systemMatrix,ll,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -this%circuitData(ll,3,kk)
          double parameterValue = m_circuitData.components.at(ll)->parameterValue;
          errFlag = MatSetValue(systemMatrix,ll,ll+m_circuitData.numberOfPressureNodes+numberOfHistoryPressures,-parameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }
        else if (m_circuitData.components.at(ll)->type == Component_Capacitor)
        {
          // insert capacitor relationship into equation system
          // this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
          int startNode = m_circuitData.components.at(ll)->startNode->indexInInputData;
          errFlag = MatSetValue(systemMatrix,ll,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
          int endNode = m_circuitData.components.at(ll)->endNode->indexInInputData;
          errFlag = MatSetValue(systemMatrix,ll,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -alfi_delt/this%circuitData(ll,3,kk)
          double parameterValue = m_circuitData.components.at(ll)->parameterValue;
          errFlag = MatSetValue(systemMatrix,ll,ll+m_circuitData.numberOfPressureNodes+numberOfHistoryPressures,-alfi_delt/parameterValue,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,this%nodeIndexToPressureHistoryNodeOrderingMap(int(this%circuitData(ll,1,kk)),kk) + this%numberOfPressureNodes(kk),kk) = -1.0d0
          errFlag = MatSetValue(systemMatrix,ll,nodeIndexToPressureHistoryNodeOrderingMap.at(startNode)+m_circuitData.numberOfPressureNodes,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,this%nodeIndexToPressureHistoryNodeOrderingMap(int(this%circuitData(ll,2,kk)),kk) + this%numberOfPressureNodes(kk),kk) = 1.0d0
          errFlag = MatSetValue(systemMatrix,ll,nodeIndexToPressureHistoryNodeOrderingMap.at(endNode)+m_circuitData.numberOfPressureNodes,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }
        else if (m_circuitData.components.at(ll)->type == Component_Inductor)
        {
          // insert inductor relationship into equation system
          // this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
          int startNode = m_circuitData.components.at(ll)->startNode->indexInInputData;
          errFlag = MatSetValue(systemMatrix,ll,toZeroIndexing(startNode),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
          int endNode = m_circuitData.components.at(ll)->endNode->indexInInputData;
          errFlag = MatSetValue(systemMatrix,ll,toZeroIndexing(endNode),-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -this%circuitData(ll,3,kk)/alfi_delt
          double parameterValue = m_circuitData.components.at(ll)->parameterValue;
          errFlag = MatSetValue(systemMatrix,ll,ll+m_circuitData.numberOfPressureNodes+numberOfHistoryPressures,-parameterValue/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,this%componentIndexToFlowHistoryComponentOrderingMap(ll,kk) + this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + this%numberOfComponents(kk),kk) = this%circuitData(ll,3,kk)/alfi_delt
          errFlag = MatSetValue(systemMatrix,ll,componentIndexToFlowHistoryComponentOrderingMap.at(ll)+m_circuitData.numberOfPressureNodes+numberOfHistoryPressures+m_circuitData.numberOfComponents,parameterValue/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }
        else
        {
          throw std::runtime_error("EE: Unknown component type in netlist. Halting.");
        }
    }

     // Do the equations for the nodes with multiple incident currents:
     for (int mm=0; mm<numberOfMultipleIncidentCurrentNodes; mm++)
     {
       for (int ll=0; ll<m_circuitData.numberOfComponents; ll++)
       {
          if (m_circuitData.components.at(ll)->endNode->indexInInputData == listOfNodesWithMultipleIncidentCurrents.at(mm))
          {
             // this%systemMatrix(mm+this%numberOfComponents(kk), this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + ll,kk) = 1.0d0
            errFlag = MatSetValue(systemMatrix,mm+m_circuitData.numberOfComponents,ll+m_circuitData.numberOfPressureNodes+numberOfHistoryPressures,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
          if (m_circuitData.components.at(ll)->startNode->indexInInputData == listOfNodesWithMultipleIncidentCurrents.at(mm))
          {
             // this%systemMatrix(mm+this%numberOfComponents(kk), this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + ll,kk) = -1.0d0
            errFlag = MatSetValue(systemMatrix,mm+m_circuitData.numberOfComponents,ll+m_circuitData.numberOfPressureNodes+numberOfHistoryPressures,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
       }
     }

     int rowsDoneSoFar = m_circuitData.numberOfComponents + numberOfMultipleIncidentCurrentNodes;

     // create the columnMap which tells us which system column each of the prescribed pressure, pressure-history or flow values belong to
     int tempUnknownVariableIndexWithinLinearSystem = 0; // just an indexing shift to keep track of where we need to write next
     // for (int ll=0; ll<m_circuitData.numberOfPrescribedPressures; ll++)
     // {	
     //   columnMap.push_back(listOfPrescribedPressures.at(ll) - 1 + tempUnknownVariableIndexWithinLinearSystem); // -1 to convert to zero-indexing (listOfPrescribedPressures is 1-indexed)
     // }
     for (auto prescribedPressureNode=m_circuitData.mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=m_circuitData.mapOfPrescribedPressureNodes.end(); prescribedPressureNode++)
     {
     	columnMap.push_back(prescribedPressureNode->second->indexLocalToSubcircuit);
     }

     tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + m_circuitData.numberOfPressureNodes; // tempUnknownVariableIndexWithinLinearSystem is zero before this line; I'm doing it like this for clarity & consistency
     for (int ll=0; ll<numberOfHistoryPressures; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + numberOfHistoryPressures;
     // for (int ll=0; ll<numberOfPrescribedFlows; ll++)
     // {
     //   columnMap.push_back(listOfPrescribedFlows.at(ll) - 1 + tempUnknownVariableIndexWithinLinearSystem); // -1 to convert listOfPrescribedFlows to zero-indexing
     // }
     for (auto prescribedFlowComponent=m_circuitData.mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=m_circuitData.mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
     {
     	columnMap.push_back(prescribedFlowComponent->second->indexLocalToSubcircuit);
     }

     tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + m_circuitData.numberOfComponents;
     for (int ll=0; ll<numberOfHistoryFlows; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     // Set the prescribed-value equations (i.e. pressure_1 (LHS) = pressure_1 (RHS) - so really just a way of setting the prescribed values within the linear system)
     for (int ll = 0; ll<systemSize - rowsDoneSoFar; ll++)
     {
       errFlag = MatSetValue(systemMatrix,rowsDoneSoFar + ll, columnMap.at(ll),1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
     }

    errFlag = MatAssemblyBegin(systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatAssemblyEnd(systemMatrix,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    errFlag = MatLUFactor(systemMatrix,NULL,NULL,NULL);CHKERRABORT(PETSC_COMM_SELF,errFlag);
}

void NetlistSubcircuit::assembleRHS(const int timestepNumber)
{

    PetscErrorCode errFlag;
    errFlag = VecZeroEntries(RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    historyPressuresInSubcircuit = pressuresInSubcircuit;
    historyFlowsInSubcircuit = flowsInSubcircuit;

    // Prescribed pressures
    int tempIndexingShift = m_circuitData.numberOfComponents + numberOfMultipleIncidentCurrentNodes;
    {
      int ll=0;
      for (auto prescribedPressureNode=m_circuitData.mapOfPrescribedPressureNodes.begin(); prescribedPressureNode!=m_circuitData.mapOfPrescribedPressureNodes.end(); prescribedPressureNode++ )
      {
        // Coming from 'f' for 'fixed' in the input data:
        if (prescribedPressureNode->second->prescribedPressureType == Pressure_Fixed)
        {
        	errFlag = VecSetValue(RHS,ll + tempIndexingShift,prescribedPressureNode->second->valueOfPrescribedPressure,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);	
        }
        // Coming from 'l' for 'left-ventricular' in the input data:
        else if (prescribedPressureNode->second->prescribedPressureType == Pressure_LeftVentricular)
        {
        	std::cerr << "this requires heart model. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
            	std::exit(1);
        }
        else
        {
        	throw std::runtime_error("Unknown pressure prescription value in Netlist.");
        }
        ll++;
      }
    }
    // for (int ll=0; ll<m_circuitData.numberOfPrescribedPressures; ll++)
    // {
    //    // 'f' for 'fixed'
    //    if (m_circuitData.typeOfPrescribedPressures.at(ll) == Pressure_Fixed)
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
    tempIndexingShift = tempIndexingShift + m_circuitData.numberOfPrescribedPressures;
    // for(int ll=0; ll<numberOfHistoryPressures; ll++)

    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
	    int lll=0;
	    for (auto node=m_circuitData.mapOfPressureNodes.begin(); node!=m_circuitData.mapOfPressureNodes.end(); node++)
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
    tempIndexingShift = tempIndexingShift + numberOfHistoryPressures;
    // Scoping unit to include the second counter ll in the for loop, without having ll in-scope after the loop finishes:
    {
    	int ll=0;
	    for (auto prescribedFlowComponent=m_circuitData.mapOfPrescribedFlowComponents.begin(); prescribedFlowComponent!=m_circuitData.mapOfPrescribedFlowComponents.end(); prescribedFlowComponent++)
	    {
	       if (prescribedFlowComponent->second->prescribedFlowType == Flow_3DInterface)
	       {
	          columnIndexOf3DInterfaceFlowInLinearSystem = ll + tempIndexingShift;
	          errFlag = VecSetValue(RHS,ll + tempIndexingShift,*flow_n_ptr,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	       }
	       else if (prescribedFlowComponent->second->prescribedFlowType == Flow_Fixed)
	       {
	          errFlag = VecSetValue(RHS,ll + tempIndexingShift,prescribedFlowComponent->second->valueOfPrescribedFlow, INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
	       }
	       else
	       {
	            throw std::runtime_error("Unknown flow prescription value in Netlist.");
	       }
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
    tempIndexingShift = tempIndexingShift + m_circuitData.numberOfPrescribedFlows;
    // Scoping unit to include the second counter lll in the for loop, without having lll in-scope after the loop finishes:
    {
	    int lll=0;
	    for (auto component=m_circuitData.mapOfComponents.begin(); component!=m_circuitData.mapOfComponents.end(); component++)
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
 
    errFlag = VecAssemblyBegin(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

}

void NetlistSubcircuit::updateInternalPressuresAndFlows()
{
    PetscErrorCode errFlag;

    if (m_circuitData.connectsTo3DDomain())
    {
    	errFlag = VecSetValue(RHS,columnIndexOf3DInterfaceFlowInLinearSystem,*flow_n_ptr,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag); //\todo make this write to the correct entry of RHS, dynamically, and read the correct pointer when there are multiple netlist LPNs
    }
    errFlag = MatMult(inverseOfSystemMatrix,RHS,solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // A self-documenting name for the request given to VecGetValues():
    int getSingleValue=1;

    // Get the updated nodal pressures:
    for (int ll=0; ll<m_circuitData.numberOfPressureNodes; ll++)
    {
        errFlag = VecGetValues(solutionVector,getSingleValue,&ll,&pressuresInSubcircuit[ll]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        m_circuitData.mapOfPressureNodes.at(toOneIndexing(ll))->pressure = pressuresInSubcircuit[ll];
    }

    // Get the updated component flows:
    int indexShift = m_circuitData.numberOfPressureNodes + numberOfHistoryPressures;
    for (int ll=indexShift; ll<m_circuitData.numberOfComponents+indexShift; ll++)
    {
        errFlag = VecGetValues(solutionVector,getSingleValue,&ll,&flowsInSubcircuit[ll-indexShift]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        m_circuitData.mapOfComponents.at(toOneIndexing(ll-indexShift))->flow = flowsInSubcircuit[ll-indexShift];
    }

    // write(*,*) 'discrepancy:', (-this%P_a(1) - this%pressuresInSubcircuit(2))/1.2862d5 - this%flowsInSubcircuit(1)
}

std::pair<double,double> NetlistSubcircuit::computeImplicitCoefficients(const int timestepNumber, const double timen_1, const double alfi_delt)
{
	assert(m_circuitData.connectsTo3DDomain());

	generateLinearSystemFromPrescribedCircuit(alfi_delt);
    assembleRHS(timestepNumber);

    PetscErrorCode errFlag;
    // get the inverse of the system matrix:
    errFlag = MatMatSolve(systemMatrix,identityMatrixForPetscInversionHack,inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatSetUnfactored(systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    // Solve the system
    errFlag = MatMult(inverseOfSystemMatrix,RHS,solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // Extract the implicit coeffcients, for eventual passing to the FORTRAN
    // linear solve
    int rowToGet[] = {0};
    int numberOfValuesToGet=1;
    PetscScalar valueFromInverseOfSystemMatrix;
    errFlag = MatGetValues(inverseOfSystemMatrix,numberOfValuesToGet,rowToGet,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem,&valueFromInverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    std::pair<double,double> returnValue;
    returnValue.first = valueFromInverseOfSystemMatrix;

    PetscScalar valueFromRHS;
    errFlag = VecGetValues(RHS,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem,&valueFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    PetscScalar valueFromSolutionVector;
    errFlag = VecGetValues(solutionVector,numberOfValuesToGet,rowToGet,&valueFromSolutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    returnValue.second = valueFromSolutionVector - valueFromInverseOfSystemMatrix * valueFromRHS;//\todo make dynamic

    return returnValue;
}

inline int NetlistSubcircuit::toZeroIndexing(const int oneIndexedValue)
{
	int zeroIndexedValue = oneIndexedValue - 1;
	return zeroIndexedValue;
}

inline int NetlistSubcircuit::toOneIndexing(const int zeroIndexedValue)
{
	int oneIndexedValue = zeroIndexedValue + 1;
	return oneIndexedValue;
}