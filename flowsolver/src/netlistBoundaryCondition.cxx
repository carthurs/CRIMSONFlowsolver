#include "netlistBoundaryCondition.hxx"
#include "fileReaders.hxx"
#include "datatypesInCpp.hxx"
#include <assert.h>

// Statics:
int netlistBoundaryCondition::numberOfInitialisedNetlistLPNs = 0;

void netlistBoundaryCondition::initialiseModel()
{
	netlistReader* netlistReader_instance = netlistReader::Instance();
	numberOfComponents = netlistReader_instance->getNumberOfComponents().at(indexOfThisNetlistLPN);

	numberOfPressureNodes = netlistReader_instance->getNumberOfPressureNodes().at(indexOfThisNetlistLPN);

	numberOfPrescribedPressures = netlistReader_instance->getNumberOfPrescribedPressures().at(indexOfThisNetlistLPN);
	numberOfPrescribedFlows = netlistReader_instance->getNumberOfPrescribedFlows().at(indexOfThisNetlistLPN);

	circuitInputData.componentTypes = netlistReader_instance->getComponentTypes().at(indexOfThisNetlistLPN);
	circuitInputData.componentStartNodes = netlistReader_instance->getComponentStartNodes().at(indexOfThisNetlistLPN);
	circuitInputData.componentEndNodes = netlistReader_instance->getComponentEndNodes().at(indexOfThisNetlistLPN);

	circuitInputData.componentParameterValues = netlistReader_instance->getComponentParameterValues().at(indexOfThisNetlistLPN);
	listOfPrescribedPressures = netlistReader_instance->getListOfPrescribedPressures().at(indexOfThisNetlistLPN);
	listOfPrescribedFlows = netlistReader_instance->getListOfPrescribedFlows().at(indexOfThisNetlistLPN);
	valueOfPrescribedPressures = netlistReader_instance->getValueOfPrescribedPressures().at(indexOfThisNetlistLPN);
	valueOfPrescribedFlows = netlistReader_instance->getValueOfPrescribedFlows().at(indexOfThisNetlistLPN);
	circuitInputData.typeOfPrescribedPressures = netlistReader_instance->getTypeOfPrescribedPressures().at(indexOfThisNetlistLPN);
	circuitInputData.typeOfPrescribedFlows = netlistReader_instance->getTypeOfPrescribedFlows().at(indexOfThisNetlistLPN);
	pressuresInLPN = netlistReader_instance->getInitialPressures().at(indexOfThisNetlistLPN);


	// integer :: rank
	// integer :: err

     // integer :: expectedNumberOfLines
     // integer :: nextReadComponentCount
     // integer :: newval
     // integer :: numberOfLines
     // integer :: iostatus
     // logical :: file_exists

     // allocate(this%localToGlobalSurfaceIndexMap(this%numberOfLPNSurfaces))

     // call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
     
     // inquire(file='netlist_surfaces_restart.dat', exist = file_exists)
     // if (file_exists) then
     //    open(unit=73,file='netlist_surfaces_restart.dat',status='old')
     //    if (rank .eq. int(0)) then
     //       write(*,*) '====> Loading netlist restart status from a previous simulation (netlist_surfaces_restart.dat)'
     //    end if
     // else
     //    inquire(file="netlist_surfaces.dat", exist=file_exists)
     //    if (file_exists) then
     //       open(unit=73,file='netlist_surfaces.dat',status='old')
     //    else
     //       if (rank .eq. int(0)) then
     //          write(*,*) 'ERROR: netlist_surfaces.dat not found. Exiting.'
     //       end if
     //       stop
     //    end if
     // end if

     // read(73,*) // Read & ignore comment line: "# List of components in a format similar to that for netlist. All comment lines must be present, but the comment content is itself irrelevant."
     
     // read(73,*) // Read & ignore comment line: "# Maximum number of components that any of the netlist boundary LPNs have:"
     // read(73,*) this%maxComponents
     // allocate(this%circuitData(this%maxComponents, 3, this%numberOfLPNSurfaces))
     // allocate(this%circuitData_componentTypes(this%maxComponents,this%numberOfLPNSurfaces))
     
     // read(73,*) // Read & ignore comment line: "# Maximum number of prescribed pressure nodes that the LPNs have:"
     // read(73,*) this%maxPrescribedPressureNodes
     // allocate(this%listOfPrescribedPressures(this%maxPrescribedPressureNodes,this%numberOfLPNSurfaces))
     // allocate(this%valueOfPrescribedPressures(this%maxPrescribedPressureNodes,this%numberOfLPNSurfaces))
     // allocate(this%typeOfPrescribedPressures(this%maxPrescribedPressureNodes,this%numberOfLPNSurfaces))

     // read(73,*) // Read & ignore comment line: "# Maximum number of prescribed flows that the LPNs have:"
     // read(73,*) this%maxPrescribedFlows
     // allocate(this%listOfPrescribedFlows(this%maxPrescribedFlows,this%numberOfLPNSurfaces))
     // allocate(this%valueOfPrescribedFlows(this%maxPrescribedFlows,this%numberOfLPNSurfaces))
     // allocate(this%typeOfPrescribedFlows(this%maxPrescribedFlows,this%numberOfLPNSurfaces))

     // read(73,*) // Read & ignore comment line: "# Maximum number of nodes that any LPN has:"
     // read(73,*) this%maxPressureNodes
     // allocate(this%pressuresInLPN(this%maxPressureNodes,this%numberOfLPNSurfaces))   // Pressure at each LPN node, using the same node indexing as in the netlist
     // allocate(this%historyPressuresInLPN(this%maxPressureNodes,this%numberOfLPNSurfaces))

     // do ii=1, this%numberOfLPNSurfaces
     //    read(73,*) // Read & ignore comment line: "### Begin first netlist boundary condition model"
     //    read(73,*) // Read & ignore comment line: "# Number of Components"
     //    read(73,*) this%numberOfComponents(ii)
     //    do jj=1, this%numberOfComponents(ii)
     //       read(73,*) // Read & ignore comment line: "# Component jj type:"
     //       read(73,*) this%circuitData_componentTypes(jj,ii)
     //       read(73,*) // Read & ignore comment line: "# Component jj details (start-node index, end-node index, associated parameter (resistance for resistors, capacitance for capacitors):"
     //       read(73,*) this%circuitData(jj,1,ii)
     //       read(73,*) this%circuitData(jj,2,ii)
     //       read(73,*) this%circuitData(jj,3,ii)
     //    end do
     //    read(73,*) // Read & ignore comment line: "# Number of prescribed pressure nodes:"
     //    read(73,*) this%numberOfPrescribedPressures(ii)
     //    if (this%numberOfPrescribedPressures(ii) .gt. int(0)) then
     //       read(73,*) // Read & ignore comment line: "# Indices of nodes with prescribed pressures:" !\todo make this work with not-three columns
     //       do jj=1, this%numberOfPrescribedPressures(ii)
     //          read(73,*) this%listOfPrescribedPressures(jj,ii)
     //       end do
     //       read(73,*) // Read & ignore comment line: "# Prescribed pressure values / scalings (dependent on types, given by next component):"
     //       do jj=1, this%numberOfPrescribedPressures(ii)
     //          read(73,*) this%valueOfPrescribedPressures(jj,ii)
     //       end do
     //       read(73,*) // Read & ignore comment line: "# Prescribed pressure types (f=fixed to value given in previous line, l=left ventricular pressure, scaled by value given in previous line):"
     //       do jj=1, this%numberOfPrescribedPressures(ii)
     //          read(73,*) this%typeOfPrescribedPressures(jj,ii)
     //       end do
     //    end if
     //    read(73,*) // Read & ignore comment line: "# Number of prescribed flows:"
     //    read(73,*) this%numberOfPrescribedFlows(ii)
     //    if (this%numberOfPrescribedFlows(ii) .gt. int(0)) then
     //       read(73,*) // Read & ignore comment line: "# Indices of components with prescribed flows"
     //       do jj=1, this%numberOfPrescribedFlows(ii)
     //          read(73,*) this%listOfPrescribedFlows(jj,ii)
     //       end do
     //       read(73,*) // Read & ignore comment line: "# Values of prescribed flows (3D interface set to -1; this value is irrelevant and unused by the code):"
     //       do jj=1, this%numberOfPrescribedFlows(ii)
     //          read(73,*) this%valueOfPrescribedFlows(jj,ii)
     //       end do
     //       read(73,*) // Read & ignore comment line: "# Types of prescribed flows (t=threeD domain interface)"
     //       do jj=1, this%numberOfPrescribedFlows(ii)
     //          read(73,*) this%typeOfPrescribedFlows(jj,ii)
     //       end do
     //    end if
     //    read(73,*) // Read & ignore comment line: "# Number of pressure nodes (including everything- 3D interface, zero-pressure points, internal nodes, etc.):"
     //    read(73,*) this%numberOfPressureNodes(ii)
     //    read(73,*) // Read & ignore comment line: "# Initial pressures at the pressure nodes:"
     //    do jj=1, this%numberOfPressureNodes(ii)
     //       read(73,*) this%pressuresInLPN(jj,ii) //\todo this will break restarts - you must load these from an alternative place in that case!
     //    end do
     //    read(73,*) // Read & ignore comment line: "### End first netlist boundary condition model"
     //  enddo
     // close(73)

     numberOfPrescribedPressuresAndFlows = numberOfPrescribedPressures + numberOfPrescribedFlows; // Just the sum of the previous two declared integers

     //\todo make these dynamic
     // allocate(this%circuitData(5,3))
     // allocate(this%circuitData_componentTypes(5))

     // this%circuitData = reshape((/1.0d0,2.0d0,2.0d0,4.0d0,4.0d0,  2.0d0,3.0d0,4.0d0,5.0d0,6.0d0,  1.286208d5,4.5d-7,1.286208d5,2.7d-6,7.5d5/), (/5,3/)) //column-major content, then the array shape
     // this%circuitData_componentTypes = (/'r','c','r','c','r'/) ! the data in here will be the stripped first column of the netilst, identifying each line of circuitData as being r=resistor, c=capacitor, etc.

     // this%circuitData(1,3) = 1.286208000000000e5
     // this%circuitData(2,3) = real(4.5000000000e-7,8)

     // open(unit=73,file='netlist.dat',status='old')
     // read(73,*) this%circuitData(1,3)
     // read(73,*) this%circuitData(2,3)
     // read(73,*) this%circuitData(3,3)
     // read(73,*) this%circuitData(4,3)
     // read(73,*) this%circuitData(5,3)
     // close(73)

     // this%circuitData(3,3) = 1.286208000000000e5
     // this%circuitData(4,3) = 2.7000000000e-6
     // this%circuitData(5,3) = real(7.5000000000e5,8)

     // this%numberOfPrescribedPressures = 3 // \todo set this dynamically
     // this%numberOfPrescribedFlows = 1 // \todo set this dynamically (when inductor added) - the 1 is for the inflow face (will need to hard-code this for now!)
     
     // allocate(this%listOfPrescribedFlows(this%numberOfPrescribedFlows))

     // this%listOfPrescribedPressures = (/3,5,6/) //\todo make dynamic
     // this%listOfPrescribedFlows = (/1/) //\todo needed? \todo make dynamic

     //\todo make these an input
     //\note that for non-fixed values (non-'f'), such as 'l' for LV pressure, the valueOfPrescribedPressure gives a scaling (so here we're using 1x LV pressure for the 2nd prescribed node)
     // this%valueOfPrescribedPressures = (/0d0,1d0,0d0/)
     // this%typeOfPrescribedPressures = (/'f','l','f'/)
     // this%typeOfPrescribedFlows = (/'3'/) //'3' for 3D domain

     // this%numberOfComponents = 5 // \todo set this dynamically
     // this%numberOfPressureNodes = 6 // \todo set this dynamically

     // Resize to contain the necessary flows, and zero it out
	flowsInLPN.resize(numberOfComponents,0.0);

	getMapOfPressHistoriesToCorrectPressNodes(); //initialises numberOfHistoryPressures
	getMapOfFlowHistoriesToCorrectComponents(); //initialises numberOfHistoryFlows


	systemSize = numberOfPressureNodes + numberOfHistoryPressures + numberOfComponents + numberOfHistoryFlows;

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

    columnMapSize = numberOfHistoryPressures + numberOfHistoryFlows + numberOfPrescribedPressures + numberOfPrescribedFlows;
    getListOfNodesWithMultipleIncidentCurrents();
//        Populate the map which converts LPN surface indices (local) to surface indices (global) - for FlowHist etc.
     // do ii=1, this%numberOfLPNSurfaces
     //    do jj=1, numCalcSrfs
     //       if(this%surfids(ii).eq.nsrflistCalc(jj)) then
     //          this%localToGlobalSurfaceIndexMap(ii) = jj
     //       endif
     //    enddo
     // enddo
}

void netlistBoundaryCondition::getMapOfPressHistoriesToCorrectPressNodes()
{
	// std::vector<int>::iterator locationOfSmallestValueSoFar;
	// int tempValue;
	// int maxHistoryPressures;
	// int writeCounter;

    for (int ii=0; ii<numberOfComponents; ii++)
    {
       // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
       if (circuitInputData.componentTypes.at(ii) == Component_Capacitor)
       {
       		listOfHistoryPressures.insert(circuitInputData.componentStartNodes.at(ii));
       		listOfHistoryPressures.insert(circuitInputData.componentEndNodes.at(ii));
       }
    }

    // First, ensure that there are any history pressures at all:
    // if (listOfHistoryPressures_temp.size() != int(0))
    // {
	   //  // Sort the history nodes into ascending order.
	   //  // This is an implementation of the simple, inefficient "Selection Sort" algorithm.
	   //  // 
	   //  // Get a pointer to the last element of the iterator:
	   //  auto justBeforeTheVectorEnd = --listOfHistoryPressures_temp.end();
	   //  for(auto iterator=listOfHistoryPressures_temp.begin(); iterator<justBeforeTheVectorEnd; iterator++)
	   //  {
	   //     locationOfSmallestValueSoFar = iterator;
	   //     for (auto iterator2 = iterator++; iterator2<justBeforeTheVectorEnd; iterator2++)
	   //     {
	   //        if(*locationOfSmallestValueSoFar > *iterator2)
	   //        {
	   //           locationOfSmallestValueSoFar = iterator2;
	   //        }
	   //     }
	   //     // Swap the values in the vector, using just the iterators which point to them:
	   //     std::iter_swap(iterator,locationOfSmallestValueSoFar);
	   //  }
	    
	   //  // Remove duplicates from the history nodes array:
	   //  // Begin by counting the number of unique entries:
    //    tempValue = listOfHistoryPressures_temp(1,kk)
    //    this%numberOfHistoryPressures(kk) = int(1)
    //    do ii = 2, 2*this%numberOfComponents(kk)
    //       if ((listOfHistoryPressures_temp(ii,kk) .ne. tempValue) .and. (listOfHistoryPressures_temp(ii,kk) .ne. int(0))) then
    //          tempValue = listOfHistoryPressures_temp(ii,kk)
    //          this%numberOfHistoryPressures(kk) = this%numberOfHistoryPressures(kk) + 1
    //       else
    //          listOfHistoryPressures_temp(ii,kk) = int(0) // we just do this for ease when removing the duplicates in the next loop... read ahead and you'll see what I mean
    //       end if
    //    end do
    // }
    // else
    // {
    //    numberOfHistoryPressures = int(0);
    // }

    numberOfHistoryPressures = listOfHistoryPressures.size();

     // Now get the unique listOfHistoryPressures array, using the number of unique entries, and the fact that listOfHistoryPressures_temp now contains unique values, interspersed with zeros
     // maxHistoryPressures = maxval(this%numberOfHistoryPressures)

     // allocate(this%nodeIndexToPressureHistoryNodeOrderingMap(this%maxPressureNodes,this%numberOfLPNSurfaces))
     // this%nodeIndexToPressureHistoryNodeOrderingMap = int(-1) // this initialisation to a clearly-invalid index should help catch indexing bugs

     // Copy the sorted values into the member array:
     // allocate(this%listOfHistoryPressures(maxHistoryPressures,this%numberOfLPNSurfaces))

     // do kk=1, this%numberOfLPNSurfaces
     //    writeCounter = 1
     //    do ii = 1, 2*this%numberOfComponents(kk)
     //       if (listOfHistoryPressures_temp(ii,kk) .ne. int(0)) then
     //          this%listOfHistoryPressures(writeCounter,kk) = listOfHistoryPressures_temp(ii,kk)
     //          writeCounter = writeCounter + 1
     //       end if
     //    end do

        // // Now do the actual generation of the pressure history node ordering map:
        // do ii=1, this%numberOfHistoryPressures(kk)
        //    this%nodeIndexToPressureHistoryNodeOrderingMap(this%listOfHistoryPressures(ii,kk),kk) = ii
        // end do
     // end do

     // Now do the actual generation of the pressure history node ordering map:
    int ii=0;
    for (auto iterator=listOfHistoryPressures.begin(); iterator != listOfHistoryPressures.end(); iterator++, ii++)
    {
       nodeIndexToPressureHistoryNodeOrderingMap.insert( std::pair<int,int> ( *iterator, ii ) );
    }
}


void netlistBoundaryCondition::getMapOfFlowHistoriesToCorrectComponents()
{

	for (int ii=0; ii<numberOfComponents; ii++)
	{
	   // Check for capacitor, as these need pressure "histories" (pressure from the previous time-step) at their end-nodes (for dP/dt term).
        if(circuitInputData.componentTypes.at(ii) == Component_Inductor)
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

      //    allocate(listOfHistoryFlows_temp(2*this%maxComponents,this%numberOfLPNSurfaces))
      //    allocate(this%numberOfHistoryFlows(this%numberOfLPNSurfaces))
      //    listOfHistoryFlows_temp = int(0)
      //    do kk = 1, this%numberOfLPNSurfaces
      //       cursor = int(0)
      //       do ii=1, this%numberOfComponents(kk)
      //          ! Check for inductor.
      //          if (this%circuitData_componentTypes(ii,kk) == 'i') then
      //             cursor = cursor + int(1)
      //             listOfHistoryFlows_temp(cursor,kk) = ii
      //          end if
      //       end do

      //       ! Sort the history flows into ascending order (will have lots of zeros at beginning of array after sort)
      //       ! This is an implementation of the simple, inefficient "Selection Sort" algorithm.
      //       do ii=1, cursor - int(1)
      //          locationOfSmallestValueSoFar = ii
      //          do jj=ii+1, cursor - 1
      //             if(listOfHistoryFlows_temp(locationOfSmallestValueSoFar,kk) > listOfHistoryFlows_temp(jj,kk)) then
      //                locationOfSmallestValueSoFar = jj
      //             end if
      //          end do
      //          tempValue = listOfHistoryFlows_temp(ii,kk)
      //          listOfHistoryFlows_temp(ii,kk) = listOfHistoryFlows_temp(locationOfSmallestValueSoFar,kk)
      //          listOfHistoryFlows_temp(locationOfSmallestValueSoFar,kk) = tempValue
      //       end do

            
      //       ! Remove duplicates from the history flows array:
      //       ! Begin by counting the number of unique entries:

      //       ! First, ensure that there are any history flows at all:
      //       if (maxval(listOfHistoryFlows_temp(:,kk)) .ne. int(0)) then
      //          tempValue = listOfHistoryFlows_temp(1,kk)
      //          this%numberOfHistoryFlows(kk) = int(1)
      //          do ii = 1, 2*this%numberOfComponents(kk)
      //             if ((listOfHistoryFlows_temp(ii,kk) .ne. tempValue) .and. (listOfHistoryFlows_temp(ii,kk) .ne. int(0))) then
      //                tempValue = listOfHistoryFlows_temp(ii,kk)
      //                this%numberOfHistoryFlows(kk) = this%numberOfHistoryFlows(kk) + 1
      //             ! else
      //             !    listOfHistoryFlows_temp(ii,kk) = int(0) ! we just do this for ease when removing the duplicates in the next loop... read ahead and you'll see what I mean
      //             end if
      //          end do
      //       else
      //          this%numberOfHistoryFlows(kk) = int(0)
      //       end if
      //       ! Now get the unique listOfHistoryFlows array, using the number of unique entries, and the fact that listOfHistoryFlows_temp now contains unique values, interspersed with zeros
      //    end do

      //    maxHistoryFlows = maxval(this%numberOfHistoryFlows)

      //    allocate(this%componentIndexToFlowHistoryComponentOrderingMap(this%maxComponents,this%numberOfLPNSurfaces))
      //    this%componentIndexToFlowHistoryComponentOrderingMap = int(-1) ! this initialisation to a clearly-invalid index should help catch indexing bugs

      //    ! Copy the sorted values into the member array:
      //    allocate(this%listOfHistoryFlows(maxHistoryFlows,this%numberOfLPNSurfaces))
      //    this%listOfHistoryFlows = int(0)
      //    do kk=1, this%numberOfLPNSurfaces
      //       do ii = 1, 2*this%numberOfComponents(kk)
      //          if (listOfHistoryFlows_temp(ii,kk) .ne. int(0)) then
      //             this%listOfHistoryFlows(ii,kk) = listOfHistoryFlows_temp(ii,kk)
      //          end if
      //       end do

      //       ! Now do the actual generation of the pressure history node ordering map:
      //       do ii=1, this%numberOfHistoryFlows(kk)
      //          this%componentIndexToFlowHistoryComponentOrderingMap(this%listOfHistoryFlows(ii,kk),kk) = ii
      //       end do
      //    end do

      // end subroutine getMapOfFlowHistoriesToCorrectComponents
}


void netlistBoundaryCondition::getListOfNodesWithMultipleIncidentCurrents()
{
    // Note that this function also counts pressure nodes which are just
    // /between/ two components, eg. for the (two resistor) subcircuit:
    //        N0--[==R1==]--N1--[==R2==]--N2
    // This would count node N1 as appearing twice, and do a "Kirchoff" current
    // balance of the form "flow through R1 = flow through R2".
    //
    // It also catches and deals with true Kirchoff equations wherea third (fourth, fifth,...)
    // component is connected to N1.

    int numberOfTimesNodeSeen;

    numberOfMultipleIncidentCurrentNodes = int(0);

    // The node data from the input file is 1-indexed, so shift this to 1:numberOfPressureNodes, instead of 0:numberOfPressureNodes-1
    for(int node=1; node < numberOfPressureNodes+1; node++)
    {
       numberOfTimesNodeSeen = int(0);
       for (int ii = 0; ii<numberOfComponents; ii++)
       {
          if ((circuitInputData.componentStartNodes.at(ii) == node) || (circuitInputData.componentEndNodes.at(ii) == node))
          {
             numberOfTimesNodeSeen++;
          }
       }
       if (numberOfTimesNodeSeen > int(1))
       {
          // this acts as a flag to the next loop, which will make the final listOfNodesWithMultipleIncidentCurrents.
          // listOfNodesWithMultipleIncidentCurrents_temp(node,kk) = int(1)
          listOfNodesWithMultipleIncidentCurrents.push_back(node);
          numberOfMultipleIncidentCurrentNodes++;
       }
    }

         // allocate(this%listOfNodesWithMultipleIncidentCurrents(int(maxMultipleIncidentCurrents),this%numberOfLPNSurfaces))
         
         // // I'm using the allocatable array nextWriteLocation2 instead of the single nextWriteLocation counter (and the =int(0) resetting inside the loop below)
         // // because we get weird crashes otherwise. This seems to be related to a bug with ifort's loop collapse optimisation -http://software.intel.com/en-us/forums/topic/505605
         // // I can avoid it by just having the nextWriteLocation2 with allocatable attribute above, or with the write(*,*) below, or with no compiler optimisations, or with -O1, or
         // // with -O2 -no-simd -no-vec. I think this hack fixes it; try again when the bug is fixed in ifort to confirm (probably ifort > 14; follow the link: http://software.intel.com/en-us/forums/topic/505605)
         // allocate(nextWriteLocation2(this%numberOfLPNSurfaces))
         // nextWriteLocation2 = int(0)

         // do kk=1, this%numberOfLPNSurfaces
         //    // Write the final listOfNodesWithMultipleIncidentCurrents by collapsing the flagged locations from the previous loop into their corresponding node indices
         //    // nextWriteLocation = int(0)
         //    do node = 1, this%numberOfPressureNodes(kk)
         //       if (listOfNodesWithMultipleIncidentCurrents_temp(node,kk) .eq. int(1)) then
         //          nextWriteLocation2(kk) = nextWriteLocation2(kk) + int(1)
         //          this%listOfNodesWithMultipleIncidentCurrents(nextWriteLocation2(kk),kk) = node
         //       end if
         //    end do
         //    // write(*,*) 'nrlis',nextWriteLocation2(kk)
         // end do
}

std::pair<double,double> netlistBoundaryCondition::computeImplicitCoefficients(int timestepNumber, double timeAtStepNplus1, double alfi_delt)
{
    generateLinearSystemFromPrescribedCircuit(alfi_delt);
    assembleRHS_netlistLPN(timestepNumber);

    PetscErrorCode errFlag;

    // inverseOfSystemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk) =
    //                invertSquareMatrix(this%systemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk));

    // get the inverse of the system matrix:
    errFlag = MatMatSolve(systemMatrix,identityMatrixForPetscInversionHack,inverseOfSystemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = MatSetUnfactored(systemMatrix); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    
    // solutionVector(1:this%systemSize(kk),kk) =                              &
    //   matmul(this%inverseOfSystemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk), this%RHS(1:this%systemSize(kk),kk))
    
    // Solve the system
    errFlag = MatMult(inverseOfSystemMatrix,RHS,solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);

    // !\todo make this generic!
    // coeff(kk,1) = this%inverseOfSystemMatrix(1,this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk)
    // coeff(kk,2) = this%solutionVector(1,kk) -                                                             &
    //                 this%inverseOfSystemMatrix(1,this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk)* &
    //                 this%RHS(this%columnIndexOf3DInterfaceFlowInLinearSystem(kk),kk) //\todo make dynamic

    std::pair<double,double> returnValue;
    PetscScalar valueFromInverseOfSystemMatrix;
    PetscScalar valueFromRHS;
    PetscScalar valueFromSolutionVector;

    int rowIndex[] = {0};
    int numberOfValuesToGet=1;
    errFlag = MatGetValues(inverseOfSystemMatrix,numberOfValuesToGet,rowIndex,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem,&valueFromInverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    returnValue.first = valueFromInverseOfSystemMatrix;

    errFlag = VecGetValues(RHS,numberOfValuesToGet,&columnIndexOf3DInterfaceFlowInLinearSystem,&valueFromRHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecGetValues(solutionVector,numberOfValuesToGet,rowIndex,&valueFromSolutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
    returnValue.second = valueFromSolutionVector - valueFromInverseOfSystemMatrix * valueFromRHS;//\todo make dynamic

    return returnValue;
}

void netlistBoundaryCondition::generateLinearSystemFromPrescribedCircuit(double alfi_delt)
{
    // This function assembles the system of (time-discretised) linear algebraic equations for the LPN.

    PetscErrorCode errFlag;

    errFlag = MatZeroEntries(systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    for(int ll=0; ll<numberOfComponents; ll++)
    {
        if (circuitInputData.componentTypes.at(ll) == Component_Resistor)
        {
          // insert resistor relationship into equation system
          errFlag = MatSetValue(systemMatrix,ll,circuitInputData.componentStartNodes.at(ll)-1,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
          errFlag = MatSetValue(systemMatrix,ll,circuitInputData.componentEndNodes.at(ll)-1,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -this%circuitData(ll,3,kk)
          errFlag = MatSetValue(systemMatrix,ll,ll+numberOfPressureNodes+numberOfHistoryPressures,-circuitInputData.componentParameterValues.at(ll),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }
        else if (circuitInputData.componentTypes.at(ll) == Component_Capacitor)
        {
          // insert capacitor relationship into equation system
          // this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
          errFlag = MatSetValue(systemMatrix,ll,circuitInputData.componentStartNodes.at(ll)-1,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
          errFlag = MatSetValue(systemMatrix,ll,circuitInputData.componentEndNodes.at(ll)-1,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -alfi_delt/this%circuitData(ll,3,kk)
          errFlag = MatSetValue(systemMatrix,ll,ll+numberOfPressureNodes+numberOfHistoryPressures,-alfi_delt/circuitInputData.componentParameterValues.at(ll),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,this%nodeIndexToPressureHistoryNodeOrderingMap(int(this%circuitData(ll,1,kk)),kk) + this%numberOfPressureNodes(kk),kk) = -1.0d0
          errFlag = MatSetValue(systemMatrix,ll,nodeIndexToPressureHistoryNodeOrderingMap[circuitInputData.componentStartNodes.at(ll)]+numberOfPressureNodes,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,this%nodeIndexToPressureHistoryNodeOrderingMap(int(this%circuitData(ll,2,kk)),kk) + this%numberOfPressureNodes(kk),kk) = 1.0d0
          errFlag = MatSetValue(systemMatrix,ll,nodeIndexToPressureHistoryNodeOrderingMap[circuitInputData.componentEndNodes.at(ll)]+numberOfPressureNodes,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }
        else if (circuitInputData.componentTypes.at(ll) == Component_Inductor)
        {
          // insert inductor relationship into equation system
          // this%systemMatrix(ll,int(this%circuitData(ll,1,kk)),kk) = 1.0d0
          errFlag = MatSetValue(systemMatrix,ll,circuitInputData.componentStartNodes.at(ll)-1,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,int(this%circuitData(ll,2,kk)),kk) = -1.0d0
          errFlag = MatSetValue(systemMatrix,ll,circuitInputData.componentEndNodes.at(ll)-1,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,ll+this%numberOfPressureNodes(kk)+this%numberOfHistoryPressures(kk),kk) = -this%circuitData(ll,3,kk)/alfi_delt
          errFlag = MatSetValue(systemMatrix,ll,ll+numberOfPressureNodes+numberOfHistoryPressures,-circuitInputData.componentParameterValues.at(ll)/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // this%systemMatrix(ll,this%componentIndexToFlowHistoryComponentOrderingMap(ll,kk) + this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + this%numberOfComponents(kk),kk) = this%circuitData(ll,3,kk)/alfi_delt
          errFlag = MatSetValue(systemMatrix,ll,componentIndexToFlowHistoryComponentOrderingMap[ll]+numberOfPressureNodes+numberOfHistoryPressures+numberOfComponents,circuitInputData.componentParameterValues.at(ll)/alfi_delt,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
        }
        else
        {
          throw std::runtime_error("EE: Unknown component type in netlist. Halting.");
        }
    }

     // Do the equations for the nodes with multiple incident currents:
     for (int mm=0; mm<numberOfMultipleIncidentCurrentNodes; mm++)
     {
       for (int ll=0; ll<numberOfComponents; ll++)
       {
          if (circuitInputData.componentEndNodes.at(ll) == listOfNodesWithMultipleIncidentCurrents.at(mm))
          {
             // this%systemMatrix(mm+this%numberOfComponents(kk), this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + ll,kk) = 1.0d0
            errFlag = MatSetValue(systemMatrix,mm+numberOfComponents,ll+numberOfPressureNodes+numberOfHistoryPressures,1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
          if (circuitInputData.componentStartNodes.at(ll) == listOfNodesWithMultipleIncidentCurrents.at(mm))
          {
             // this%systemMatrix(mm+this%numberOfComponents(kk), this%numberOfPressureNodes(kk) + this%numberOfHistoryPressures(kk) + ll,kk) = -1.0d0
            errFlag = MatSetValue(systemMatrix,mm+numberOfComponents,ll+numberOfPressureNodes+numberOfHistoryPressures,-1.0,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          }
       }
     }

     int rowsDoneSoFar = numberOfComponents + numberOfMultipleIncidentCurrentNodes;

     // create the columnMap which tells us which system column each of the prescribed pressure, pressure-history or flow values belong to
     int tempUnknownVariableIndexWithinLinearSystem = 0; // just an indexing shift to keep track of where we need to write next
     for (int ll=0; ll<numberOfPrescribedPressures; ll++)
     {
       columnMap.push_back(listOfPrescribedPressures.at(ll) - 1 + tempUnknownVariableIndexWithinLinearSystem); // -1 to convert to zero-indexing (listOfPrescribedPressures is 1-indexed)
     }

     tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + numberOfPressureNodes; // tempUnknownVariableIndexWithinLinearSystem is zero before this line; I'm doing it like this for clarity & consistency
     for (int ll=0; ll<numberOfHistoryPressures; ll++)
     {
       columnMap.push_back(ll + tempUnknownVariableIndexWithinLinearSystem);
     }

     tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + numberOfHistoryPressures;
     for (int ll=0; ll<numberOfPrescribedFlows; ll++)
     {
       columnMap.push_back(listOfPrescribedFlows.at(ll) - 1 + tempUnknownVariableIndexWithinLinearSystem); // -1 to convert listOfPrescribedFlows to zero-indexing
     }

     tempUnknownVariableIndexWithinLinearSystem = tempUnknownVariableIndexWithinLinearSystem + numberOfComponents;
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

void netlistBoundaryCondition::assembleRHS_netlistLPN(int timestepNumber)
{

    PetscErrorCode errFlag;
    errFlag = VecZeroEntries(RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);

    historyPressuresInLPN = pressuresInLPN;

    // Prescribed pressures
    int tempIndexingShift = numberOfComponents + numberOfMultipleIncidentCurrentNodes;
    for (int ll=0; ll<numberOfPrescribedPressures; ll++)
    {
       // 'f' for 'fixed'
       if (circuitInputData.typeOfPrescribedPressures.at(ll) == Pressure_Fixed)
       {
          errFlag = VecSetValue(RHS,ll + tempIndexingShift,valueOfPrescribedPressures.at(ll),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
       }
       // 'l' for 'leftVentricular'
       else if (circuitInputData.typeOfPrescribedPressures.at(ll) == Pressure_LeftVentricular)
       {
          std::cout << "this requires heartmodel. Also should make boundaryConditionManager able to provide P_IM..whatevers." << std::endl;
          std::exit(1);
          // if ((timestepNumber == 0) || (timestepNumber == 1)) // treat case with no known IM pressure yet
          // {
          //    P_IM_mid_lasttimestep = 5000; // \todo find a better way of doing this; maybe input this value from file...
          //    P_IM_mid = 5000; // ... or set it based on the aortic valve state at simulation start
          // }
          // elseif (timestepNumber .eq. int(2)) then // treat case where only one IM pressure history point is known
          //    P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber-1)
          //    P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber)
          // else // get the previous intramyocardial pressure in the case where we have enough doata for this (see comment before start of "if" block)
          //    P_IM_mid_lasttimestep = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber-1)
          //    P_IM_mid = this%valueOfPrescribedPressures(ll,kk) * hrt%plv_hist(timestepNumber)
          // end if

          // errFlag = VecSetValue(RHS,ll + tempIndexingShift,P_IM_mid,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
          // int nn=0;
          // for (auto iterator=listOfHistoryPressures.begin(); iterator!=listOfHistoryPressures.end(); iterator++, nn++)
          // {
          //    if (*iterator == listOfPrescribedPressures.at(ll))
          //    {
          //       historyPressuresInLPN.at(*iterator) = P_IM_mid_lasttimestep;
          //    }
          // }
       }
       else
       {
            throw std::runtime_error("Unknown pressure prescription value in Netlist.");
       }
    }
    // History Pressures
    tempIndexingShift = tempIndexingShift + numberOfPrescribedPressures;
    // for(int ll=0; ll<numberOfHistoryPressures; ll++)
    int lll=0;
    for (auto iterator=listOfHistoryPressures.begin(); iterator!=listOfHistoryPressures.end(); iterator++, lll++)
    {
       // this%RHS(ll + tempIndexingShift,kk) = this%historyPressuresInLPN(this%listOfHistoryPressures(ll,kk),kk)
        errFlag = VecSetValue(RHS,lll+tempIndexingShift,historyPressuresInLPN.at(*iterator - 1),INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
    // Prescribed Flows
    tempIndexingShift = tempIndexingShift + numberOfHistoryPressures;
    for(int ll=0; ll<numberOfPrescribedFlows; ll++)
    {
       if (circuitInputData.typeOfPrescribedFlows.at(ll) == Flow_3DInterface)
       {
          columnIndexOf3DInterfaceFlowInLinearSystem = ll + tempIndexingShift;
          errFlag = VecSetValue(RHS,ll + tempIndexingShift,*flow_n_ptr,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
       }
       else if (circuitInputData.typeOfPrescribedFlows.at(ll) == Flow_Fixed)
       {
          errFlag = VecSetValue(RHS,ll + tempIndexingShift,valueOfPrescribedFlows.at(ll), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
       }
       else
       {
            throw std::runtime_error("Unknown flow prescription value in Netlist.");
       }
    }
    // History Flows
    tempIndexingShift = tempIndexingShift + numberOfPrescribedFlows;
    lll=0;
    // for(int ll=0; ll<numberOfHistoryFlows; ll++)
    for (auto iterator=listOfHistoryFlows.begin(); iterator!=listOfHistoryFlows.end(); iterator++, lll++)
    {
       // this%RHS(ll + tempIndexingShift,kk) = this%flowsInLPN(this%listOfHistoryFlows(ll,kk),kk)
       errFlag = VecSetValue(RHS,lll + tempIndexingShift,flowsInLPN.at(*iterator - 1), INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }
 
    errFlag = VecAssemblyBegin(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    errFlag = VecAssemblyEnd(RHS); CHKERRABORT(PETSC_COMM_SELF,errFlag);

}

// void netlistBoundaryCondition::updateLPN_netlistLPN()
void netlistBoundaryCondition::updateLPN()
{
    PetscErrorCode errFlag;

    errFlag = VecSetValue(RHS,columnIndexOf3DInterfaceFlowInLinearSystem,*flow_n_ptr,INSERT_VALUES); CHKERRABORT(PETSC_COMM_SELF,errFlag); //\todo make this write to the correct entry of RHS, dynamically, and read the correct pointer when there are multiple netlist LPNs
    errFlag = MatMult(inverseOfSystemMatrix,RHS,solutionVector); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    // this%solutionVector(1:this%systemSize(kk),kk) = matmul(this%inverseOfSystemMatrix(1:this%systemSize(kk),1:this%systemSize(kk),kk), this%RHS(1:this%systemSize(kk),kk))

    // A self-documenting name for the request given to VecGetValues():
    int getSingleValue=1;

    // Get the updated nodal pressures:
    for (int ll=0; ll<numberOfPressureNodes; ll++)
    {
        errFlag = VecGetValues(solutionVector,getSingleValue,&ll,&pressuresInLPN[ll]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }

    // Get the updated component flows:
    int indexShift = numberOfPressureNodes + numberOfHistoryPressures;
    for (int ll=indexShift; ll<numberOfComponents+indexShift; ll++)
    {
        errFlag = VecGetValues(solutionVector,getSingleValue,&ll,&flowsInLPN[ll-indexShift]); CHKERRABORT(PETSC_COMM_SELF,errFlag);
    }

    // write(*,*) 'discrepancy:', (-this%P_a(1) - this%pressuresInLPN(2))/1.2862d5 - this%flowsInLPN(1)

}