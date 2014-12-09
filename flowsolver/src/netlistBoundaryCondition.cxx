#include "netlistBoundaryCondition.hxx"

// Statics:
int netlistBoundaryCondition::numberOfInitialisedNetlistLPNs = 0;

void netlistBoundaryCondition::initialiseModel()
{
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

     // Resize to contai
	flowsInLPN.resize(numberOfComponents,0.0);

	// call this%getMapOfPressHistoriesToCorrectPressNodes() //initialises this%numberOfHistoryPressures
	// call this%getMapOfFlowHistoriesToCorrectComponents() //initialises this%numberOfHistoryFlows

	systemSize = numberOfPressureNodes + numberOfHistoryPressures + numberOfComponents + numberOfHistoryFlows;

    PetscErrorCode errFlag;
	errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&systemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
	errFlag = MatCreateSeqDense(PETSC_COMM_SELF,systemSize,systemSize,NULL,&inverseOfSystemMatrix);CHKERRABORT(PETSC_COMM_SELF,errFlag);
	
	errFlag = VecCreate(PETSC_COMM_SELF,&RHS);CHKERRABORT(PETSC_COMM_SELF,errFlag);
	errFlag = VecSetType(RHS,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make RHS a VECSEQ sequential vector

	errFlag = VecCreate(PETSC_COMM_SELF,&solutionVector);CHKERRABORT(PETSC_COMM_SELF,errFlag);
	errFlag = VecSetType(solutionVector,VECSEQ);CHKERRABORT(PETSC_COMM_SELF,errFlag); // Make RHS a VECSEQ sequential vector

    columnMapSize = numberOfHistoryPressures + numberOfHistoryFlows + numberOfPrescribedPressures + numberOfPrescribedFlows;

    // call this%getListOfNodesWithMultipleIncidentCurrents()

//        Populate the map which converts LPN surface indices (local) to surface indices (global) - for FlowHist etc.
     // do ii=1, this%numberOfLPNSurfaces
     //    do jj=1, numCalcSrfs
     //       if(this%surfids(ii).eq.nsrflistCalc(jj)) then
     //          this%localToGlobalSurfaceIndexMap(ii) = jj
     //       endif
     //    enddo
     // enddo
}

