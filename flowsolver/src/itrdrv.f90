module itrDrvVars ! this needs to be cleaned up

    use memLS

    integer npermdims !perhaps these can go in global or common
    integer ntmpdims

    integer npermdimss
    integer ntmpdimss

    integer nstp
    real*8 rmub
    real*8 rmue

    integer    ifuncs(6)

    real*8  almit, alfit, gamit

    !integer istp
    real*8 xi
    integer ilss
    integer istepc
    integer isolve
    real*8 deltt, dtglt
    integer iupdate
    real*8 soln, asoln, sgn, asolng
    integer ierr
    integer ifile
    integer ntoutv

    real*8 tfact

    !for memLS
    INTEGER memLS_nFaces, gnNo, nNo, faIn, facenNo
    INTEGER, ALLOCATABLE :: ltg(:), gNodes(:)
    REAL*8, ALLOCATABLE :: sV(:,:)

    CHARACTER*128 fileName
    TYPE(memLS_commuType) communicator
    TYPE(memLS_lhsType) memLS_lhs
    TYPE(memLS_lsType) memLS_ls

    ! memLS copies for the heart model - systolic/diastolic 
    INTEGER memLS_nFaces_s
    INTEGER, ALLOCATABLE :: gNodes_s(:)
    REAL*8, ALLOCATABLE :: sV_s(:,:)
    TYPE(memLS_lhsType) memLS_lhs_d      
    TYPE(memLS_lhsType) memLS_lhs_s   
    TYPE(memLS_commuType) communicator_s
    TYPE(memLS_lsType) memLS_ls_d
    TYPE(memLS_lsType) memLS_ls_s

    ! iBC copies for heart only, the original iBC is now in globalArrays.f90
    real*8, allocatable, dimension(:,:) :: iBCs
    real*8, allocatable, dimension(:)   :: iBCd    

end module

!! This iterative driver is the semi-discrete, predictor multi-corrector
!! algorithm. It contains the Hulbert Generalized Alpha method which
!! is 2nd order accurate for Rho_inf from 0 to 1.  The method can be
!! made  first-order accurate by setting Rho_inf=-1. It uses CGP and
!! GMRES iterative solvers.
!!
!! working arrays:
!!  y      (nshg,ndof)           : Y variables
!!  x      (nshg,nsd)            : node coordinates
!!  iBC    (nshg)                : BC codes
!!  BC     (nshg,ndofBC)         : BC constraint parameters
!!  iper   (nshg)                : periodicity table
!!
!!
!! Zdenek Johan,  Winter 1991.  (Fortran 90)
!! Alberto Figueroa, Winter 2004.  CMM-FSI
!! Irene Vignon, Fall 2004. Impedance BC

!
!.... ---------------> initialization and pre-processing <---------------
!
subroutine itrdrv_init() bind(C, name="itrdrv_init")

    use iso_c_binding
    use cpp_interface
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl
    use phcommonvars
    use itrDrvVars
    use multidomain

    implicit none
    !IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    integer i,j,jj,k

    integer lesid
    integer nkvecs, nsclrsol, nsolflow
    integer icnt
    integer ndofs

    integer ntempdims
    integer isolsc
    integer indx

    integer lstep0

    integer :: surfids(0:MAXSURF)

    ! added target for gfortran
    integer, dimension(nshg), target :: binaryMask

    !
    !.... For linear solver Library
    !
    integer eqnType, prjFlag, presPrjFlag, verbose  ! init
    character*1024    servername ! init
    character*20        license_f_name ! init

    integer numberOfCppManagedBoundaryConditions

!--------------------------------------------------------------------
!   Setting up memLS

    IF (memLSFlag .EQ. 1) THEN
        CALL memLS_LS_CREATE(memLS_ls, LS_TYPE_NS, dimKry=Kspace,relTol=epstol(8), &
                             relTolIn=(/epstol(1),epstol(7)/), maxItr=nPrjs, maxItrIn=(/nGMRES,maxIters/))

        CALL memLS_COMMU_CREATE(communicator, MPI_COMM_WORLD)

        ! make heart model copies of the memLS objects
        IF (iheart .GT. int(0)) THEN
            
            ! default memLS object has prescribed_velocity on the inflow face, in the heart model velocity = 0, i.e. diastole
            memLS_ls_d = memLS_ls
            
            ! recreate memLS_ls for systole, we setup the reduced order boundary condition later
            CALL memLS_LS_CREATE(memLS_ls_s, LS_TYPE_NS, dimKry=Kspace, relTol=epstol(8), & 
                                 relTolIn=(/epstol(1),epstol(7)/), maxItr=nPrjs, maxItrIn=(/nGMRES,maxIters/))

            ! recreate communicator for the systolic copy
            CALL memLS_COMMU_CREATE(communicator_s, MPI_COMM_WORLD)

        END IF 

        IF (numpe .GT. 1) THEN
            WRITE(fileName,*) myrank
            fileName = "ltg.dat."//ADJUSTL(TRIM(fileName))
            OPEN(1,FILE=fileName)
            READ(1,*) gnNo
            READ(1,*) nNo
            if (allocated(ltg)) then
              DEALLOCATE(ltg)
            endif
            ALLOCATE(ltg(nNo))
            READ(1,*) ltg
            CLOSE(1)
        ELSE
            gnNo = nshg
            nNo = nshg
            if (allocated(ltg)) then
              DEALLOCATE(ltg)
            endif
            ALLOCATE(ltg(nNo))
            DO i=1, nNo
                ltg(i) = i
            END DO
        END IF
    ELSE

    ! find the machine name so that we set the license key properly

        call MPI_BARRIER(INEWCOMM,ierr)

        license_f_name='license.dat'

        call SolverLicenseServer(servername)

    END IF

    !
    ! only master should be verbose
    !

    if(numpe.gt.0 .and. myrank.ne.master)iverbose=0
    !
    inquire(file='xyzts.dat',exist=exts)

    if(exts) then

        open(unit=626,file='xyzts.dat',status='old')
        read(626,*) ntspts, freq, tolpt, iterat, varcod
        call sTD             ! sets data structures

        do jj=1,ntspts       ! read coordinate data where solution desired
            read(626,*) ptts(jj,1),ptts(jj,2),ptts(jj,3)
        enddo

        varts = zero

    endif
    !
    !.... open history and aerodynamic forces files
    !
    if (myrank .eq. master) then
        open (unit=ihist,  file=fhist,  status='unknown')
        open (unit=iforce, file=fforce, status='unknown')
        open (unit=76, file="fort.76", status='unknown')
    endif
    !
    !.... initialize
    !
    ifuncs(:)  = 0              ! func. evaluation counter
    istep  = 0
    !yold   = y
    !acold  = ac

    rerr = zero
    ybar = yold

    if (ideformwall.eq.1) then
        ubar = uold
    endif

    xdist = zero
    xdnv = zero
    df_fem = zero

    !
    !.... ---------------> initialize LesLib Library <---------------
    !
    !.... assign parameter values
    !
    do i = 1, 100
        numeqns(i) = i
    enddo
    nKvecs       = Kspace
    prjFlag      = iprjFlag
    presPrjFlag  = ipresPrjFlag
    verbose      = iverbose
    !
    !.... determine how many scalar equations we are going to need to solve
    !
    nsolt=mod(impl(1),2)      ! 1 if solving temperature
    nsclrsol=nsolt+nsclr      ! total number of scalars solved At
                              ! some point we probably want to create
                              ! a map, considering stepseq(), to find
                              ! what is actually solved and only
                              ! dimension lhs to the appropriate
                              ! size. (see 1.6.1 and earlier for a
                              ! "failed" attempt at this).


    nsolflow=mod(impl(1),100)/10  ! 1 if solving flow

    !
    !.... Now, call lesNew routine to initialize
    !     memory space
    !
    call genadj(colm, rowp, icnt )  ! preprocess the adjacency list

    numberOfCppManagedBoundaryConditions = 0
    call callCPPGetNumberOfCppManagedBoundaryConditions(numberOfCppManagedBoundaryConditions)

    ! If there are netlist boundary conditions, we will get CPP to provide the
    ! details on how iBC should be edited in order to implement Dirichlet
    ! or Neumann boundary conditions on a surface, dependent on whether
    ! flow is permitted (or prescribed!) across that surface.
    !
    ! The major reason for this functionality is for netlists which have valves;
    ! in this case there will be a switching of the boundary condition type
    ! if the valves block the flow across the 3D interface.

    ! Tell the boundary condition objects which boundary nodes belong to their surface:
    call callCPPGiveBoundaryConditionsListsOfTheirAssociatedMeshNodes(c_loc(ndsurf), nshg)
    ! Get the correct boundary condition type flag array iBC for the current valve state
    ! in the netlists:
    if (numberOfCppManagedBoundaryConditions .gt. int(0)) then
        if (.not. allocated(iBC_original)) then
            allocate(iBC_original(nshg))
            iBC_original = iBC
        endif

        ! Fill out the array with ones
        binaryMask = 1
        ! Set the apprporiate zeros for the nodes where the Dirichlet conditions should be applied:
        call callCPPGetBinaryMaskToAdjustNodalBoundaryConditions(c_loc(binaryMask), nshg)
        ! Reset iBC, so we work with a clean copy
        iBC = iBC_original
        ! zero out the entries of iBC where the boundary condition type has become
        ! Neumann at that node, as annotated by binaryMask from the CPP boundary
        ! condition objects.
        where(binaryMask .eq. int(0))
            iBC = int(0)
        end where
    endif

    nnz_tot=icnt ! this is exactly the number of non-zero blocks on
                 ! this proc

    if (nsolflow.eq.1) then
        lesId   = numeqns(1)
        eqnType = 1
        nDofs   = 4

        !--------------------------------------------------------------------
        !     Rest of configuration of memLS is added here, where we have LHS
        !     pointers
        IF (memLSFlag .EQ. 1) THEN

            ! This big commented block has been moved to the function rebuildMemLS_lhs, so that we can call it whenever
            ! the state of the netlist valves changes.
            !
            ! This approach replaces the "two memLS" approach previously used for the heart model; instead of
            ! pre-constructing a systole and a diastole version of memLS, we'll just poll the boundary conditions
            ! to ask them if their valve state has changed since the last time-step, in such a way to cause the 
            ! flow to be newly permitted / banned across one of the netlist surfaces. In that case,
            ! we rebuild memLS for the new state.

            call rebuildMemLS_lhs()


            ! IF  (ipvsq .GE. 2) THEN
            !     memLS_nFaces = 1 + numResistSrfs + numImpSrfs + numRCRSrfs + numGRCRSrfs + numControlledCoronarySrfs + numNetlistLPNSrfs
            ! ELSE
            !     memLS_nFaces = 1
            ! END IF

            ! CALL memLS_LHS_CREATE(memLS_lhs, communicator, gnNo, nNo, nnz_tot, ltg, colm, rowp, memLS_nFaces)

            ! faIn = 1
            ! facenNo = 0
            ! DO i=1, nshg
            !     IF (IBITS(iBC(i),3,3) .NE. 0)  facenNo = facenNo + 1
            ! END DO
            ! if ((.not.allocated(gNodes)).and.(.not.allocated(sV))) then
            !   ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
            ! endif
            ! sV = 0D0
            ! j = 0
            ! DO i=1, nshg
            !     IF (IBITS(iBC(i),3,3) .NE. 0) THEN
            !         j = j + 1
            !         gNodes(j) = i
            !         IF (.NOT.BTEST(iBC(i),3)) sV(1,j) = 1D0
            !         IF (.NOT.BTEST(iBC(i),4)) sV(2,j) = 1D0
            !         IF (.NOT.BTEST(iBC(i),5)) sV(3,j) = 1D0
            !     END IF
            ! END DO
            ! CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo, nsd, BC_TYPE_Dir, gNodes, sV)

            ! IF  (ipvsq .GE. 2) THEN
            !     DO k = 1, numResistSrfs
            !         faIn = faIn + 1
            !         CALL AddNeumannBCTomemLS(nsrflistResist(k), faIn, memLS_lhs)
            !     END DO
            !     DO k = 1, numImpSrfs
            !         faIn = faIn + 1
            !         CALL AddNeumannBCTomemLS(nsrflistImp(k), faIn, memLS_lhs)
            !     END DO
            !     DO k = 1, numRCRSrfs
            !         faIn = faIn + 1
            !         CALL AddNeumannBCTomemLS(nsrflistRCR(k), faIn, memLS_lhs)
            !     END DO
            !     DO k = 1, numGRCRSrfs
            !         faIn = faIn + 1
            !         CALL AddNeumannBCTomemLS(nsrflistGRCR(k), faIn, memLS_lhs)
            !     END DO
            !     do k=1, numControlledCoronarySrfs
            !         faIn = faIn + 1
            !         call AddNeumannBCTomemLS(indicesOfCoronarySurfaces(k), faIn, memLS_lhs)
            !     end do
            !     do k=1, numNetlistLPNSrfs
            !         faIn = faIn + 1
            !         call AddNeumannBCTomemLS(indicesOfNetlistSurfaces(k),faIn,memLS_lhs)
            !     end do
            ! END IF

            ! create systolic memLS_lhs
            IF (iheart .gt. int(0)) THEN

                ! add 1 for the heart model surface
                memLS_nFaces_s = memLS_nFaces + int(1)
                ! create diastolic memLS_lhs
                memLS_lhs_d = memLS_lhs                
                ! get surface id of the heart model
                surfids = hrt%getsurfids()
                ! here 1 is added to the number of surfaces in memLS_nFaces_s
                CALL memLS_LHS_CREATE(memLS_lhs_s, communicator_s, gnNo, nNo, nnz_tot, ltg, colm, rowp, memLS_nFaces_s) 


                ! here we are counting the number of dirchlet nodes excluding the heart model face
                faIn = 1
                facenNo = 0
                DO i=1, nshg
                    IF (IBITS(iBC(i),3,3) .NE. 0 .AND. ndsurf(i) .NE. surfids(1)) THEN
                        facenNo = facenNo + 1
                    END IF
                END DO
                if ((.not.allocated(gNodes_s)).and.(.not.allocated(sV_s))) then
                  ALLOCATE(gNodes_s(facenNo), sV_s(nsd,facenNo))
                endif
                sV_s = 0D0

                j = 0

                DO i = 1, nshg               
                    IF (IBITS(iBC(i),3,3) .NE. 0 .AND. ndsurf(i) .NE. surfids(1)) THEN
                        j = j + 1
                        gNodes_s(j) = i
                        IF (.NOT.BTEST(iBC(i),3)) sV_s(1,j) = 1D0
                        IF (.NOT.BTEST(iBC(i),4)) sV_s(2,j) = 1D0
                        IF (.NOT.BTEST(iBC(i),5)) sV_s(3,j) = 1D0
                    END IF
                END DO

                CALL memLS_BC_CREATE(memLS_lhs_s, faIn, facenNo, nsd, BC_TYPE_Dir, gNodes_s, sV_s)

                ! same as above we add the reduced order surfaces
                IF  (ipvsq .GE. 2) THEN
                    DO k = 1, numResistSrfs
                        faIn = faIn + 1
                        CALL AddNeumannBCTomemLS(nsrflistResist(k), faIn, memLS_lhs_s)
                    END DO
                    DO k = 1, numImpSrfs
                        faIn = faIn + 1
                        CALL AddNeumannBCTomemLS(nsrflistImp(k), faIn, memLS_lhs_s)
                    END DO
                    DO k = 1, numRCRSrfs
                        faIn = faIn + 1
                        CALL AddNeumannBCTomemLS(nsrflistRCR(k), faIn, memLS_lhs_s)
                    END DO
                    DO k = 1, numGRCRSrfs
                        faIn = faIn + 1
                        CALL AddNeumannBCTomemLS(nsrflistGRCR(k), faIn, memLS_lhs_s)
                    END DO
                    do k=1, numControlledCoronarySrfs
                        faIn = faIn + 1
                        call AddNeumannBCTomemLS(indicesOfCoronarySurfaces(k), faIn, memLS_lhs_s)
                    end do
                    do k=1, numNetlistLPNSrfs
                        faIn = faIn + 1
                        call AddNeumannBCTomemLS(indicesOfNetlistSurfaces(k), faIn, memLS_lhs_s)
                    end do

                    ! here we add the heart model surface
                    faIn = faIn + 1 ! add one to skip to next surface
                    CALL AddNeumannBCTomemLS(surfids(1), faIn, memLS_lhs_s)  
                END IF

            END IF 

        ELSE
!--------------------------------------------------------------------

#ifndef NO_ACUSIM
            call myfLesNew( lesId,          41994, &
                            eqnType, &
                            nDofs,          minIters,       maxIters, &
                            nKvecs,         prjFlag,        nPrjs, &
                            presPrjFlag,    nPresPrjs,      epstol(1), &
                            prestol,        verbose,        statsflow, &
                            nPermDims,      nTmpDims,      servername  )
#endif
        END IF

        if (.not.allocated(aperm)) then
          allocate (aperm(nshg,nPermDims))
        endif
        if (.not.allocated(atemp)) then
          allocate (atemp(nshg,nTmpDims))
        endif
        if (.not.allocated(lhsP)) then
          allocate (lhsP(4,nnz_tot))
        endif
        if (.not.allocated(lhsK)) then
          allocate (lhsK(9,nnz_tot))
        endif

        call readLesRestart( lesId,  aperm, nshg, myrank, currentTimestepIndex, nPermDims )

    else
        nPermDims = 0
        nTempDims = 0
    endif


    if(nsclrsol.gt.0) then
        do isolsc=1,nsclrsol
            lesId       = numeqns(isolsc+1)
            eqnType     = 2
            nDofs       = 1
            presPrjflag = 0
            nPresPrjs   = 0
            prjFlag     = 1
            indx=isolsc+2-nsolt ! complicated to keep epstol(2) for
                             ! temperature followed by scalars
            call myfLesNew( lesId,          41994, &
            eqnType, &
            nDofs,          minIters,       maxIters, &
            nKvecs,         prjFlag,        nPrjs, &
            presPrjFlag,    nPresPrjs,      epstol(indx), &
            prestol,        verbose,        statssclr, &
            nPermDimsS,     nTmpDimsS,      servername )
        enddo
        !
        !  Assume all scalars have the same size needs
        !
        if (.not.allocated(apermS)) then
          allocate (apermS(nshg,nPermDimsS,nsclrsol))
        endif
        if (.not.allocated(atempS)) then
          allocate (atempS(nshg,nTmpDimsS))  !they can all share this
        endif
        if (.not.allocated(lhsS)) then
          allocate (lhsS(nnz_tot,nsclrsol))
        endif
    !
    ! actually they could even share with atemp but leave that for later
    !
    else
        nPermDimsS = 0
        nTmpDimsS  = 0
    endif
    !
    !...  prepare lumped mass if needed
    !
    if((flmpr.ne.0).or.(flmpl.ne.0)) call genlmass(x, shp,shgl)
    !
    !.... -----------------> End of initialization <-----------------
    !
    !.....open the necessary files to gather time series
    !
    lstep0 = currentTimestepIndex+1
        
    ! ******************************************************** ! 
    ! *** initialize the initial condition for heart model *** !
    ! ******************************************************** !

    ! make copies of iBC        
    if (iheart .gt. int(0)) then
        if (.not.allocated(iBCd)) then
          allocate(iBCd(nshg))
        endif
        iBCd = iBC
        if (.not.allocated(iBCs)) then
          allocate(iBCs(1,nshg))
        endif
        iBCs(1,:) = iBC
        surfids = hrt%getsurfids()
        do i = 1,nshg
            if(ndsurf(i) .eq. surfids(1)) then
                iBCs(1,i)=0
            end if
        end do
    end if

  
    if(iheart .gt. int(0)) then
        
        ! set pressures and flows at 3D/0D interface
        call initreducedordermodel(y,hrt,'multidomain') 
        
        ! assign pointers for the filter
        call hrt%assign_ptrs_ext()

        ! set internal variables for the left heart     
        ! if (isystemic .ne. int(1)) then
            ! call hrt%initxvars(currentTimestepIndex)                   
        ! end if 

        call hrt%initxvars(currentTimestepIndex)                   

        ! if valve open, switch iBC
        if (hrt%isavopen() .eq. 1) then
            iBC = iBCs(1,:)
        else
            iBC = iBCd
        end if        

    end if

    ! ******************************************************** !
    ! ***                                                  *** !
    ! ******************************************************** !
    
    !
    !.... satisfy the boundary conditions
    !
    if (currentTimestepIndex .eq. 0 ) then
        call itrBC (y, ac,  iBC, BC, iper, ilwork)
    endif
    yold = y
    acold = ac

    !
    !.... loop through the time sequences
    !

    itseq = 1

    !AD         tcorecp1 = second(0)
    !AD         tcorewc1 = second(-1)
    !
    !.... set up the time integration parameters
    !
    nstp   = nstep(itseq)
    nitr   = niter(itseq)
    LCtime = loctim(itseq)
    dtol(:)= deltol(itseq,:)

    call itrSetup ( y, acold ) ! sets up alfi

   
    !
    ! *** set \alpha_{i} in multidomain module
    !
    call setsimv_alfi(alfi)

    !...initialize the coefficients for the impedance convolution,
    !   which are functions of alphaf so need to do it after itrSetup
    if(numImpSrfs.gt.zero) then
        call calcImpConvCoef (numImpSrfs, ntimeptpT)
    endif
    !
    !...initialize the initial condition P(0)-RQ(0)-Pd(0) for RCR BC
    !   need ndsurf so should be after initNABI
    if(numRCRSrfs.gt.zero) then
        call calcRCRic(y,nsrflistRCR,numRCRSrfs)
    endif
    !
    !...compute area and initial flow and pressure of TRCR boundary surfaces
    !
    if(numTRCRSrfs.gt.zero) then
        call calcTRCRic(y,nsrflistTRCR,numTRCRSrfs)
    endif


    !---------------------------------- Nan rcr
    if (numGRCRSrfs .gt. 0 ) then
        !call grcrbc_Initialize()
        !call grcrbc_SetInternalState(y)
        
        if (nrcractive .eq. 1) then
            nrcr = nrcrconstructor(numGRCRSrfs,nsrflistGRCR)
            
            !call initreducedordermodel(y, nrcr, 'legacy')
            call initreducedordermodel(y, nrcr, 'multidomain')

            ! assign pointers for the filter
            call nrcr%assign_ptrs_ext()

        end if

    endif
   !---------------------------------- Nan rcr
   !
    !...calculate area and initial pressure and flow for CalcSurfaces
    !
    if(numCalcSrfs.gt.zero) then
        call calcCalcic(y,nsrflistCalc,numCalcSrfs)
    endif
    !
    !...initialize the initial condition for Coronary BC
    !   need ndsurf so should be after initNABI
    if(numCORSrfs.gt.zero) then
        call calcCORic(y,nsrflistCOR,numCORSrfs)
    endif
    !
    !.... allocate LHS and RHS arrays required for the constrained surfaces
    !
    if(Lagrange.gt.zero) then
        call calcLagrangeic(nsrflistLagrange,numLagrangeSrfs)
    endif

    ! ************************************************************ !
    ! *** intialise multidomain container and systemic circuit *** !
    ! ************************************************************ !

    if (multidomainactive .eq. 1) then
        call initmultidomaincontainer(y,multidom) ! set using flows and pressures at t = t_{n}
        !if (sysactive) then
        !  call sys%initxvars(currentTimestepIndex)
        !  call initreducedordermodel(y,sys,'multidomain')
        !end if
    end if

    ! ************************************************************ !
    ! ************************************************************ !

    !
    !.... deformable wall initialization
    !
    if(ideformwall.eq.1) then
        call vlmwStTri(uref,nodetagfield,x)
        ! call solveWallProb(rowp,colm,ilwork,iBC,BC,iper)
    end if


    if (idistancenudge.eq.1) then
        if (myrank.eq.zero) then
            write(*,*) "computing distance to wall data surfaces (initial)"
        end if
        call ElmDist(uold,x,xdist,xdnv,df_fem)
    end if

    !
    !  find the last solve of the flow in the step sequence so that we will
    !         know when we are at/near end of step
    !
    !         ilast=0
    nitr=0  ! count number of flow solves in a step (# of iterations)
    do i=1,seqsize
        if(stepseq(i).eq.0) nitr=nitr+1
    enddo
    !
    !.... loop through the time steps
    !
    istop=0
    rmub=datmat(1,2,1)
    !if(rmutarget.gt.0) then
    !    rmue=rmutarget
    !else
        rmue=datmat(1,2,1) ! keep constant
    !endif

    ! This function is needed to poke the RCRs to get them to re-read
    ! the pressure pointer (pressure_n_ptr) to initialise with the 3D
    ! domain pressure at the surface.
    ! 
    ! This is only done if this is not a restarted simulation (its internal
    ! logic guards this).
    !
    ! It's a least-bad hack to get around the chicken-and-egg initialisation
    ! order issues for the 3D and multidomain. If we want to avoid this,
    ! the Fortran code needs refactoring.
    call callCPPSetPressureFromFortran()
    ! call callCPPLoadAllNetlistComponentFlowsAndNodalPressures()

    CONTAINS

    ! SUBROUTINE AddNeumannBCTomemLS(srfID, faIn)

    ! INTEGER, INTENT(IN) :: srfID, faIn

    ! INTEGER facenNo, i, j

    ! facenNo = 0
    ! DO i = 1, nshg
    !  IF (srfID .EQ. ndsurf(i)) THEN
    !     facenNo = facenNo + 1
    !  END IF
    ! END DO
    ! IF (ALLOCATED(gNodes)) DEALLOCATE(gNodes, sV)
    ! ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
    ! sV = 0D0
    ! j = 0
    ! DO i = 1, nshg
    !  IF (srfID .EQ. ndsurf(i)) THEN
    !     j = j + 1
    !     gNodes(j) = i
    !     sV(:,j) = NABI(i,1:3)
    !  END IF
    ! END DO
    ! CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo, nsd, BC_TYPE_Neu, gNodes, sV)

    ! RETURN
    ! END SUBROUTINE AddNeumannBCTomemLS

end subroutine itrdrv_init
!
! replacement subroutine for systolic/diastolic copies
!

SUBROUTINE AddNeumannBCTomemLS(srfID, faIn, memLS_lhs_inout)
    use iso_c_binding
    use cpp_interface
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl
    use phcommonvars
    use itrDrvVars, only : facenNo
    use memLS
    use multidomain

    implicit none

    INTEGER, INTENT(IN) :: srfID, faIn
      
    INTEGER i, j

    TYPE(memLS_lhsType), INTENT(INOUT) :: memLS_lhs_inout
    INTEGER, ALLOCATABLE :: gNodes_tmp(:)
    REAL*8, ALLOCATABLE :: sV_tmp(:,:)

    facenNo = 0
    DO i = 1, nshg
        IF (srfID .EQ. ndsurf(i)) THEN
            facenNo = facenNo + 1
        END IF
    END DO
    IF (ALLOCATED(gNodes_tmp)) THEN
        DEALLOCATE(gNodes_tmp)
    END IF
    IF (ALLOCATED(sV_tmp)) THEN      
        DEALLOCATE(sV_tmp)
    END IF

    if ((.not.allocated(gNodes_tmp)).and.(.not.allocated(sV_tmp))) then
      ALLOCATE(gNodes_tmp(facenNo), sV_tmp(nsd,facenNo))
    endif

    sV_tmp = 0D0
    gNodes_tmp = int(0)
    j = 0
    DO i = 1, nshg
        IF (srfID .EQ. ndsurf(i)) THEN
            j = j + 1
            gNodes_tmp(j) = i
            sV_tmp(:,j) = NABI(i,1:3)
        END IF
    END DO

    CALL memLS_BC_CREATE(memLS_lhs_inout, faIn, facenNo, nsd, BC_TYPE_Neu, gNodes_tmp, sV_tmp)      


    RETURN
END SUBROUTINE AddNeumannBCTomemLS

! This subroutine is the place for initialisation calls which require that the CPP
! boundary conditions have already been initialised
subroutine initialiseInfoForCPPBoundaryConditions()
    implicit none
end subroutine initialiseInfoForCPPBoundaryConditions

subroutine rebuildMemLS_lhs()
        use iso_c_binding
        use cpp_interface
        use shapeTable
        use globalArrays
        use pvsQbi     !gives us splag (the spmass at the end of this run
        use specialBC !gives us itvn
        use timedata   !allows collection of time series
        use convolImpFlow !for Imp bc
        use convolRCRFlow !for RCR bc
        use convolTRCRFlow !for time-varying RCR bc

        !use grcrbc ! Nan rcr

        use convolCORFlow !for Coronary bc
        use incpBC        !for INCP bc
        use calcFlowPressure !to save history of flow and pressure of bc surfaces
        use LagrangeMultipliers
        use deformableWall
        use ResidualControl
        use phcommonvars
        use itrDrvVars
        use multidomain

        implicit none

        integer numBCsWhichDisallowFlow
        integer flowIsPermitted
        integer i, j, k
        integer thisIsANetlistSurface
        integer thisSurfacePermitsFlow
        integer surface

        IF  (ipvsq .GE. 2) THEN
            call callCPPGetNumberOfBoundaryConditionsWhichCurrentlyDisallowFlow(numBCsWhichDisallowFlow)
            memLS_nFaces = 1 + numResistSrfs + numImpSrfs + numRCRSrfs + numGRCRSrfs + numControlledCoronarySrfs + numNetlistLPNSrfs - numBCsWhichDisallowFlow
            ! write(*,*) "numBCsWhichDisallowFlow:", numBCsWhichDisallowFlow
        ELSE
            memLS_nFaces = 1
        END IF

        ! Free the LHS, but first:
        ! check whether LHS has been created (so we don't attempt to free a non-existent LHS):
        if (memLS_lhs%foC) then
            call memLS_LHS_FREE(memLS_lhs)
        endif

        CALL memLS_LHS_CREATE(memLS_lhs, communicator, gnNo, nNo, nnz_tot, ltg, colm, rowp, memLS_nFaces)

        ! here we are counting the number of dirchlet nodes excluding the Netlists
        faIn = 1
        facenNo = 0
        DO i=1, nshg
            IF (IBITS(iBC(i),3,3) .NE. 0)  then
                ! Check whether this node is on one a surface belonging to a Netlist boundary condition:
                thisIsANetlistSurface = int(0)
                thisSurfacePermitsFlow = int(1)
                do surface = 1, numNetlistLPNSrfs
                    if (indicesOfNetlistSurfaces(surface) .eq. ndsurf(i)) then
                        thisIsANetlistSurface = int(1)
                        call callCPPDiscoverWhetherFlowPermittedAcrossSurface(indicesOfNetlistSurfaces(surface),thisSurfacePermitsFlow)
                    end if
                end do

                ! If this is not a netlist surface:
                if((thisIsANetlistSurface .eq. int(0)) .or. (thisSurfacePermitsFlow .eq. int(0))) then
                    facenNo = facenNo + 1
                endif
            ENDIF
        END DO
        if (allocated(gNodes)) then
            deallocate(gNodes)
        endif
        ALLOCATE(gNodes(facenNo))

        if (allocated(sV)) then
            deallocate(sV)
        endif
        allocate(sV(nsd,facenNo))


        sV = 0D0
        j = 0
        DO i=1, nshg
            IF (IBITS(iBC(i),3,3) .NE. 0) THEN
                ! Check whether this node is on one a surface belonging to a Netlist boundary condition:
                thisIsANetlistSurface = int(0)
                thisSurfacePermitsFlow = int(1)
                do surface = 1, numNetlistLPNSrfs
                    if (indicesOfNetlistSurfaces(surface) .eq. ndsurf(i)) then
                        thisIsANetlistSurface = int(1)
                        call callCPPDiscoverWhetherFlowPermittedAcrossSurface(indicesOfNetlistSurfaces(surface),thisSurfacePermitsFlow)
                    end if
                end do

                ! If this is not a netlist surface:
                if((thisIsANetlistSurface .eq. int(0)) .or. (thisSurfacePermitsFlow .eq. int(0))) then
                    j = j + 1
                    gNodes(j) = i
                    IF (.NOT.BTEST(iBC(i),3)) sV(1,j) = 1D0
                    IF (.NOT.BTEST(iBC(i),4)) sV(2,j) = 1D0
                    IF (.NOT.BTEST(iBC(i),5)) sV(3,j) = 1D0
                endif
            END IF
        END DO
        CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo, nsd, BC_TYPE_Dir, gNodes, sV)

        IF  (ipvsq .GE. 2) THEN
            DO k = 1, numResistSrfs
                faIn = faIn + 1
                CALL AddNeumannBCTomemLS(nsrflistResist(k), faIn, memLS_lhs)
            END DO
            DO k = 1, numImpSrfs
                faIn = faIn + 1
                CALL AddNeumannBCTomemLS(nsrflistImp(k), faIn, memLS_lhs)
            END DO
            DO k = 1, numRCRSrfs
                faIn = faIn + 1
                CALL AddNeumannBCTomemLS(nsrflistRCR(k), faIn, memLS_lhs)
            END DO
            DO k = 1, numGRCRSrfs
                faIn = faIn + 1
                CALL AddNeumannBCTomemLS(nsrflistGRCR(k), faIn, memLS_lhs)
            END DO
            do k=1, numControlledCoronarySrfs
                faIn = faIn + 1
                call AddNeumannBCTomemLS(indicesOfCoronarySurfaces(k), faIn, memLS_lhs)
            end do
            do k=1, numNetlistLPNSrfs
                call callCPPDiscoverWhetherFlowPermittedAcrossSurface(indicesOfNetlistSurfaces(k),flowIsPermitted)
                if (flowIsPermitted .eq. int(1)) then
                    faIn = faIn + 1
                    call AddNeumannBCTomemLS(indicesOfNetlistSurfaces(k),faIn,memLS_lhs)
                endif
            end do
        END IF
    end subroutine rebuildMemLS_lhs


!
!.... ---------------> initialize the time step <---------------
!
subroutine itrdrv_iter_init() bind(C, name="itrdrv_iter_init")

    use iso_c_binding
    use cpp_interface
    use errorManagement, only: write_to_stderr
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    use multidomain, only: nrcr, nrcractive, hrt
    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl
    use phcommonvars
    use itrDrvVars

    implicit none
    !IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    ! added target for gfortran
    integer, dimension(nshg), target :: binaryMask
    integer numberOfCppManagedBoundaryConditions, ii

    ! ! Ensure that the CurrentIter counter has been reset (detects e.g. problems with
    ! ! solver.inp requesting a MinNumIter which exceeds the number of steps in the
    ! ! step construction)
    ! if(rescontrol .gt. 0) then
    !     if (CurrentIter .ne. 0) then
    !         call write_to_stderr("WARNING: CurrentIter not reset to zero. Does Minimum Required Iterations exceed the number of steps in Step Construction in solver.inp?")
    !     end if
    ! end if

    call callCPPInitialiseLPNAtStartOfTimestep_netlist()

    ! write(*,*) "setting pressure for C++ RCRs 2"
    ! call callCPPSetPressureFromFortran() ! enable for kalman filter

! ********************************************* !
! *** heart model boundary condition switch *** !
! ********************************************* !

          if (iheart .eq. 1) then  
            if (hrt%isavopen() .eq. 1) then
              iBC = iBCs(1,:)
            else
              iBC = iBCd
            end if
          end if

        ! Get the correct boundary condition type flag array iBC for the current valve state
        ! in the netlists:
        numberOfCppManagedBoundaryConditions = 0
        call callCPPGetNumberOfCppManagedBoundaryConditions(numberOfCppManagedBoundaryConditions)
        if (numberOfCppManagedBoundaryConditions .gt. int(0)) then
            ! Fill out the array with zeros
            binaryMask = 1
            ! Set the apprporiate zeros for the nodes where the Neumann conditions should be applied:
            call callCPPGetBinaryMaskToAdjustNodalBoundaryConditions(c_loc(binaryMask), nshg)
            ! Reset iBC, so we work with a clean copy
            iBC = iBC_original
            ! zero out the entries of iBC where the boundary condition type is
            ! Neumann at that node, as annotated by binaryMask from the CPP boundary
            ! condition objects.
            where(binaryMask .eq. int(0))
                iBC = int(0)
            end where
        endif


! ********************************************* c
! ********************************************* c

    xi=(istep+1)*1.0/nstp
    datmat(1,2,1)=rmub*(1.0-xi)+xi*rmue
    !            write(*,*) "current mol. visc = ", datmat(1,2,1)
    
!.... if we have time varying boundary conditions update the values of BC.
!     these will be for time step n+1 so use currentTimestepIndex+1
    
    if (itvn .gt. 0) then
       call BCint((currentTimestepIndex+1)*Delt(1), x, BC, iBC)
    end if 

    !
    ! ... calc the pressure contribution that depends on the history for the imp BC
    !
    if(numImpSrfs.gt.0) call pHist(poldImp,QHistImp,ImpConvCoef, &
    ntimeptpT,numImpSrfs)
    !
    ! ... calc the pressure contribution that depends on the history for the RCR BC
    !
    if(numRCRSrfs.gt.0) then
        call CalcHopRCR (Delt(itseq), currentTimestepIndex, numRCRSrfs)
        call CalcRCRConvCoef(currentTimestepIndex,numRCRSrfs)
        call pHist(poldRCR,QHistRCR,RCRConvCoef, &
        nstep+nptsRCR,numRCRSrfs)
    endif
    !
    ! ... compute coefficients required for time-varying RCR BCs
    !
    if(numTRCRSrfs.gt.0) then
        call CalcTRCRConvCoef(currentTimestepIndex,numTRCRSrfs)
        call CalcHopTRCR (Delt(itseq), currentTimestepIndex, numTRCRSrfs)
        call pHist(poldTRCR,QHistTRCR,TRCRConvCoef, &
        nstep(1)+nptsTRCR,numTRCRSrfs)
    endif

    ! Kalman filter reset status of boundary conditions after the estimation step:
    if (kalmanFilterOn) then
        call resetBoundaryConditionStateForKalmanFilter(yold, numGRCRSrfs, nsrflistGRCR)
        call resetBoundaryConditionStateForKalmanFilter(yold, numControlledCoronarySrfs, indicesOfCoronarySurfaces)
        call resetBoundaryConditionStateForKalmanFilter(yold, numNetlistLPNSrfs, indicesOfNetlistSurfaces)
    endif

    ! ! Nan rcr ----------------------------------
    ! if(numGRCRSrfs.gt.0) then
    !     !call grcrbc_ComputeImplicitCoefficients(currentTimestepIndex, yold)

    !     if (nrcractive) then
    !         ! first reset flow for the filter
    !         call reset_flow_n(yold, nrcr)
    !         ! calculate the implicit coefficients
    !         ! call nrcr%setimplicitcoeff(currentTimestepIndex) !\cppHook

    !         ! call callCppComputeAllImplicitCoeff_solve(currentTimestepIndex)
    !         ! call callCppComputeAllImplicitCoeff_update(currentTimestepIndex)
    !     end if

    ! endif

    ! Moved here from above - it's a generic update for everything,
    ! so shouldn't need guarding
    call callCppComputeAllImplicitCoeff_solve(currentTimestepIndex)
    call callCppComputeAllImplicitCoeff_update(currentTimestepIndex)
    ! ------------------------------------------

    !
    ! ... calc the pressure contribution that depends on the history for the Coronary BC
    !
    if(numCORSrfs.gt.0) then
        call CalcCORConvCoef(currentTimestepIndex,numCORSrfs)
        call pHist(poldCOR, QHistCOR, CORConvCoef, &
        nstep+nptsCOR,numCORSrfs)
        call CalcHopCOR (Delt(itseq), currentTimestepIndex, nsrflistCOR,  &
        numCORSrfs,yold)
    endif
    !
    !.... calculate the coefficients
    !
    if (numINCPSrfs .gt. 0) then
        call CalcINCPConvCoef(currentTimestepIndex, numINCPSrfs)
        call pHist(poldINCP, QHistINCP, &
        INCPConvCoef, nstep+nptsINCP, numINCPSrfs)
        call CalcINCPCoef(Delt(itseq), currentTimestepIndex, &
        nsrflistINCP, numINCPSrfs, yold)
    endif
    !
    ! Decay of scalars
    !
    if(nsclr.gt.0 .and. tdecay.ne.1) then
        yold(:,6:ndof)=y(:,6:ndof)*tdecay
        BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*tdecay
    endif

    if(nosource.eq.1) BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*0.8


!    if(iLES.gt.0) then  !complicated stuff has moved to
!                                       !routine below
!        call lesmodels(yold,  acold,     shgl,      shp,  &
!        iper,  ilwork,    rowp,      colm, &
!        nsons, ifath,     x,    &
!        iBC,   BC)
!
!
!    endif

    !.... set traction BCs for modeled walls
    !
!    if (itwmod.ne.0) then
!        call asbwmod(yold,   acold,   x,      BC,     iBC, &
!        iper,   ilwork,  ifath,  velbar)
!    endif

end subroutine itrdrv_iter_init

!
!.... ---------------> a single time step <---------------
!
subroutine itrdrv_iter_step() bind(C, name="itrdrv_iter_step")

    use iso_c_binding
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    use multidomain, only: nrcr, nrcractive, hrt, multidom, multidomainactive, newCoronaryActive
    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl
    use phcommonvars
    use itrDrvVars

    use cpp_interface

    implicit none
    !IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    integer j
    integer boundaryConditionRebuildNeeded

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !.... -----------------------> predictor phase <-----------------------
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call itrPredict(yold, y,   acold,  ac ,  uold,  u)

    call itrBC (y,  ac,  iBC,  BC,  iper,ilwork)

    if(nsolt.eq.1) then
        isclr=0
        call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
    endif
    do isclr=1,nsclr
        call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
    enddo


    iter=0
    ilss=0  ! this is a switch thrown on first solve of LS redistance


    do istepc=1,seqsize
        icode=stepseq(istepc)
        if(mod(icode,5).eq.0) then ! this is a solve
            isolve=icode/10
            if(icode.eq.0) then ! flow solve (encoded as 0)
                !
                iter   = iter+1
                ifuncs(1)  = ifuncs(1) + 1
                !
                Force(1) = zero
                Force(2) = zero
                Force(3) = zero
                HFlux    = zero
                lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1)))

                ! ************************************ !
                ! heart model switch for memLS_lhs *** !
                ! ************************************ !

                if (iheart .gt. int(0)) then
                    if (hrt%isavopen() .eq. 1) then
                        memLS_lhs = memLS_lhs_s
                        memLS_ls = memLS_ls_s
                    else
                        memLS_lhs = memLS_lhs_d
                        memLS_ls = memLS_ls_d
                    end if
                end if

                ! If we have netlists, check if a change in valve state within the netlists
                ! means that we need to change boundary condition type (i.e. rebuild memLS_lhs)
                ! ...
                ! This replaces the old method used for the stand-alone heart model, which worked
                ! with two memLS systems, one for systole, and one for diastole.
                if ( (memLSFlag .EQ. 1) .and. (numNetlistLPNSrfs .gt. int(0))) then
                    ! if a rebuild is needed:
                    ! call callCPPHaveBoundaryConditionTypesChanged(boundaryConditionRebuildNeeded)
                    ! write(*,*) "BC state change flag: ", boundaryConditionRebuildNeeded
                    ! if (boundaryConditionRebuildNeeded .eq. int(1)) then
                        ! write(*,*) "rebuilding memLS linear system..."
                        call rebuildMemLS_lhs()
                    ! endif
                endif

                ! ************************************ !                    
                ! ************************************ !

                call SolFlow(y,             ac,        u, &
                yold,          acold,     uold, &
                x,             xdist,     xdnv, &
                iBC,           BC,        res, &
                nPermDims,     nTmpDims,  aperm, &
                atemp,         iper,           &
                ilwork,        shp,       shgl, &
                shpb,          shglb,     rowp,      &
                colm,          lhsK,      lhsP, &
                solinc,        rerr, &
                memLS_lhs,     memLS_ls,  memLS_nFaces)


            else          ! scalar type solve
                if (icode.eq.5) then ! Solve for Temperature
                            ! (encoded as (nsclr+1)*10)
                    isclr=0
                    ifuncs(2)  = ifuncs(2) + 1
                    j=1
                else       ! solve a scalar  (encoded at isclr*10)
                    isclr=isolve
                    ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                    j=isclr+nsolt
                    if((iLSet.eq.2).and.(ilss.eq.0) &
                    .and.(isclr.eq.2)) then
                        ilss=1 ! throw switch (once per step)
                        y(:,7)=y(:,6) ! redistance field initialized
                        ac(:,7)   = zero
                        call itrBCSclr ( y,  ac,  iBC,  BC, iper, &
                        ilwork)
                        !
                        !....store the flow alpha, gamma parameter values and assigm them the
                        !....Backward Euler parameters to solve the second levelset scalar
                        !
                        alfit=alfi
                        gamit=gami
                        almit=almi
                        Deltt=Delt(1)
                        Dtglt=Dtgl
                        alfi = 1
                        gami = 1
                        almi = 1
                        !     Delt(1)= Deltt ! Give a pseudo time step
                        Dtgl = one / Delt(1)
                    endif  ! level set eq. 2
                endif ! deciding between temperature and scalar

                lhs = 1 - min(1,mod(ifuncs(isclr+2)-1, &
                LHSupd(isclr+2)))

                call SolSclr(y,             ac,        u, &
                yold,          acold,     uold, &
                x,             iBC, &
                BC,            nPermDimsS,nTmpDimsS,   &
                apermS(1,1,j), atempS,    iper, &
                ilwork,        shp,       shgl, &
                shpb,          shglb,     rowp,      &
                colm,          lhsS(1,j),  &
                solinc(1,isclr+5))


            endif         ! end of scalar type solve

        else ! this is an update  (mod did not equal zero)

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !.... -----------------------> corrector phase <-----------------------
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            iupdate=icode/10  ! what to update

            if(icode.eq.1) then !update flow

                call itrCorrect ( y,    ac,    u,   solinc)
                call itrBC (y,  ac,  iBC,  BC, iper, ilwork)

            else  ! update scalar

                isclr=iupdate  !unless
                if(icode.eq.6) isclr=0
                !if(iRANS.lt.0)then  ! RANS
                !    call itrCorrectSclrPos(y,ac,solinc(1,isclr+5))
                !else
                    call itrCorrectSclr (y, ac, solinc(1,isclr+5))
                !endif
                if (ilset.eq.2 .and. isclr.eq.2)  then
                    if (ivconstraint .eq. 1) then
                        call itrBCSclr (  y,  ac,  iBC,  BC, iper, &
                        ilwork)
                        !
                        ! ... applying the volume constraint on second level set scalar
                        !
                        call solvecon (y,    x,      iBC,  BC,  &
                        iper, ilwork, shp,  shgl)
                    !
                    endif   ! end of volume constraint calculations
                endif      ! end of redistance calculations
                !
                call itrBCSclr (  y,  ac,  iBC,  BC, iper, &
                ilwork)
            endif
        endif         !end of switch between solve or update

        if(rescontrol .gt. 0) then
            if (controlResidual .lt. ResCriteria .and.  &
            CurrentIter .ge. MinNumIter) then
                CurrentIter = 0

                exit
            endif
        endif

    enddo            ! loop over sequence in step

    if (ioform .eq. 2) then

        call stsGetStats( y,      yold,     ac,     acold, &
        u,      uold,      &
        x,      xdist,    xdnv, &
        shp,    shgl,     shpb,   shglb, &
        iBC,    BC,       iper,   ilwork, &
        rowp,   colm,     lhsK,   lhsP )

    endif

    !
    !  Find the solution at the end of the timestep and move it to old
    !
    !
    ! ...First to reassign the parameters for the original time integrator scheme
    !
    if((iLSet.eq.2).and.(ilss.eq.1)) then
        alfi =alfit
        gami =gamit
        almi =almit
        Delt(1)=Deltt
        Dtgl =Dtglt
    endif

    call itrUpdate( yold,  acold,   uold,  y,    ac,   u)

    !
    ! once all non-linear iterations finished, update the container with the final y 
    ! note this is being done for each particle, it must also be done with the final 
    ! y in the finalise step after the filter has averaged y
    !
    if (multidomainactive .eq. 1) then
        call updmultidomaincontainer(y,multidom,'velocity')
        call updmultidomaincontainer(y,multidom,'pressure')
    end if 


    ! Nan rcr ----------------------------------
    ! update P,Q variables
    if(numGRCRSrfs.gt.0) then
        !call grcrbc_UpdateInternalState(y)

        ! this is here because it is particle dependent as the WK pressure is part of the state
        ! the pressure_n is a set via the pointer set previously 
        if (nrcractive .eq. 1) then
            ! call updreducedordermodel(y,nrcr,'update') !\todo add this back in maybe when you do the Kalman filter?
            !
            ! At least when the kalman filter is off, this only needs to be done
            ! at the end of the time-step, so I moved it to itrdrv_iter_finalize()
            ! call callCPPUpdateAllRCRS_Pressure_n1_withflow()
        end if

    endif

    
    ! ------------------------------------------

    ! if(numControlledCoronarySrfs .gt. 0) then
    !     call callCppUpdateAllControlledCoronaryLPNs()
    ! endif

    ! if(numNetlistLPNSrfs .gt. 0) then
    !     call callCPPUpdateAllNetlistLPNs()
    ! endif

    ! interface to compute distances to observed wall motion
!    if (imeasdist.eq.1) then
!        if (myrank.eq.zero) then
!            write(*,*) "computing distance to wall data surfaces (end of step)"
!        end if
!        call ElmDist(u,x,xdist,xdnv,df_fem)
!    end if

end subroutine itrdrv_iter_step

!
!.... ---------------> finalize the time step <---------------
!

subroutine itrdrv_iter_finalize() bind(C, name="itrdrv_iter_finalize")

    use iso_c_binding
    use cpp_interface
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl

    use multidomain, only: nrcr, hrt, multidomainactive, multidom, newCoronaryActive, nrcractive

    use phcommonvars
    use itrDrvVars
    use ale

    implicit none
    !IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    integer jj
    integer ifail

    real*8 uMesh1(nshg), uMesh2(nshg), uMesh3(nshg)
    real*8 relativeVelocity(nshg,3)
    real*8 updatedMeshCoordinates(nshg,3)

    ! ! Update boundary conditions to the final pressure, conforming to the final flow:
    ! if (nrcractive) then
    !     call callCPPUpdateAllRCRS_Pressure_n1_withflow()
    ! end if

    if (newCoronaryActive .eq. 1) then
        ! call callCPPUpdateAllControlledCoronaryLPNs_Pressure_n1_withflow()

        ! One final update of the internal pressures in the LPN to conform to the 
        ! end-of-timestep flow, and
        ! save the LPN internal pressures.
        call callCppfinalizeLPNAtEndOfTimestep_controlledCoronary()
    end if

    if(numNetlistLPNSrfs .gt. 0) then
        call callCPPUpdateAllNetlistLPNs(currentTimestepIndex)
        call callCppfinalizeLPNAtEndOfTimestep_netlists()
    endif

    call itrBC (yold, acold,  iBC,  BC,  iper,ilwork)

    istep = istep + 1
    currentTimestepIndex = currentTimestepIndex + 1

    ! interface to compute distances to observed wall motion
    if (imeasdist.eq.1) then
        if (myrank.eq.zero) then
            write(*,*) "computing distance to wall data surfaces (end of step -- finalize)"
        end if
        call ElmDist(u,x,xdist,xdnv,df_fem)
    end if

    !
    !.... in case yold,acold,uold were changed externally
    !.... sync y,ac,u with yold,acold,uold
    !
    y = yold
    ac = acold
    u = uold

    !
    !.... update average displacement
    !
    if (iupdateprestress.eq.1) then
         tfact = one/istep

         ubar = tfact*uold + (one-tfact)*ubar

    endif

    !
    !.... update reference displacement
    !.... at the very last time step
    !
    !
    if (nstep(1).eq.istep .and. iupdateprestress.eq.1) then
        uref = uref + ubar
        u = u - ubar
        uold = uold - ubar
        write(*,*) 'updating uref with avg'
    end if

    !
    ! ... calculate relative velocity KDL, MA
    ! 

    call getMeshVelocities(uMesh1, uMesh2, uMesh3, nshg)
    relativeVelocity(:,1) = y(:,1) - uMesh1(:)
    relativeVelocity(:,2) = y(:,2) - uMesh2(:)
    relativeVelocity(:,3) = y(:,3) - uMesh3(:)

    ! calculate updated mesh coordinates
    ! x is stored in the global arrays module
    updatedMeshCoordinates(:,1) = x(:,1)
    updatedMeshCoordinates(:,2) = x(:,2)
    updatedMeshCoordinates(:,3) = x(:,3)

    !
    ! ... write out the solution
    !
    if ((irs .ge. 1) .and. (mod(currentTimestepIndex, ntout) .eq. 0)) then
        call restar ('out ', yold, ac)
        if(ideformwall.eq.1) then
            call write_displ(myrank, currentTimestepIndex, nshg, 3, uold, uref)
            if (imeasdist.eq.1) then
                call write_distl(myrank, currentTimestepIndex, nshg, 1, xdist) ! should use nshg or numnp?
            end if            
        end if

        if (aleOn .eq. int(1)) then
#if DEBUG_ALE == 1
            call Write_Residual(myrank, currentTimestepIndex, nshg, 4, res) 
#endif
            ! call Write_Relative_Velocity(myrank, currentTimestepIndex, nshg, 3, relativeVelocity) 
            call appendDoubleFieldToRestart(myrank, currentTimestepIndex, nshg, 3, relativeVelocity, "relative velocity") 
            call appendDoubleFieldToRestart(myrank, currentTimestepIndex, nshg, 3, updatedMeshCoordinates, "updated mesh coordinates") 
        end if 
        
    endif
    
    ! ************************** !
    ! *** update heart model *** !
    ! ************************** !
    
    !if (iheart .gt. int(0) .and. isystemic .ne. int(1)) then
    if (iheart .gt. int(0)) then

        
        call hrt%iterate_hrt(currentTimestepIndex,'update')
        call updreducedordermodel(y,hrt,'update') 
        
        ! update internal variables in the heart with new flow at t = t_{n+1}
        call hrt%updxvars(currentTimestepIndex)        
        call hrt%writexvars(currentTimestepIndex)
        
        ! call hrt%write_activation_history_hrt(currentTimestepIndex) !moved this to avoid missing latest activation value...

    end if

    ! ************************** !
    ! ************************** !

    !
    ! ... update the flow history for the impedance convolution, filter it and write it out
    !
    if(numImpSrfs.gt.zero) then
        call UpdHistConv(y,nsrflistImp,numImpSrfs) !uses Delt(1)
    endif

    !
    ! ... update the flow history for the RCR convolution
    !
    if(numRCRSrfs.gt.zero) then
        call UpdHistConv(y,nsrflistRCR,numRCRSrfs) !uses currentTimestepIndex
        call UpdRCR(y,nsrflistRCR,numRCRSrfs)
    endif

    !
    ! UPDATE CONTAINER WITH AVERAGE FROM THE PARTICLES
    ! THIS SHOULD BE MOVED UP, BUT WILL THAT PUSH THE HEART MODEL OUT OF SYNC?!!
    !
    if (multidomainactive .eq. 1) then
        call updmultidomaincontainer(y,multidom,'velocity')
        call updmultidomaincontainer(y,multidom,'pressure')
    end if 

    if (nrcractive .eq. 1) then
        call callCPPUpdateAllRCRS_Pressure_n1_withflow()
    end if
    !call callCPPUpdateAllRCRS_Pressure_n1_withflow()
    !
    ! *** update flow and pressure history in the numerical RCR
    !
    if(numGRCRSrfs.gt.zero) then        
        
        ! here we re-update the RCR with the innovated y (for the flow only)
        ! the corresponding innovated pressure_n is externally set by the filter directly (see the subroutine updreducedordermodel)
        call updreducedordermodel(y,nrcr,'update-flow-only')
         ! Update boundary conditions to the final pressure, conforming to the final flow:

        ! update pressure and flow history arrays
        call nrcr%updxvars(currentTimestepIndex)

        ! write pressure and flow history arrays
        ! call nrcr%writexvars(currentTimestepIndex)


    endif

    ! wrt fcns 
    call callCPPRecordPressuresAndFlowsInHistoryArrays()

    ! wrt dta
    if ((irs .ge. 1) .and. (mod(currentTimestepIndex, ntout) .eq. 0)) then
        if(myrank.eq.zero) then
            call callCPPWritePHistAndQHistRCR()
            call callCPPWriteAllNetlistComponentFlowsAndNodalPressures()
        endif
    end if 
    !Now the files are in the computer!



    !
    ! ... update flow and pressure history of time-varying RCR BCs
    !
    if(numTRCRSrfs.gt.zero) then
        call UpdHistConv(y, nsrflistTRCR, numTRCRSrfs)
        call UpdTRCR(y, nsrflistTRCR, numTRCRSrfs)
    endif

!    ! Nan rcr ----------------------------------
!    if(numGRCRSrfs.gt.0) then
!        call grcrbc_UpdateInternalState(y)
!    endif
!    ! ------------------------------------------

    !
    ! ... update the flow history for the Coronary convolution
    !
    if(numCORSrfs.gt.zero) then
        call UpdHistConv(y,nsrflistCOR,numCORSrfs) !uses currentTimestepIndex
        call UpdHistPlvConv(y,Delt(itseq),currentTimestepIndex,nsrflistCOR, numCORSrfs)
    endif
    !
    ! ... update the flow history for the CalcSurfaces
    !
    if(numCalcSrfs.gt.zero) then
        call Updcalc(y,nsrflistCalc,numCalcSrfs)
    endif
    !
    !.... calculate the values of constraint functions and write Lagrange Multipliers
    !
    if (Lagrange .gt. 0) then
        call UpdateLagrangeCoef(y,colm,rowp,nsrflistLagrange,numLagrangeSrfs)
    endif

    !
    !.... compute the consistent boundary flux
    !
    !if(abs(itwmod).ne.1) &
    call Bflux ( yold,      acold,      uold,      &
    x,         xdist,      xdnv, &
    shp,       shgl,       shpb,    &
    shglb,     ilwork,     iBC, &
    BC,        iper)


    !...  dump TIME SERIES

    if (exts) then

        if (mod(currentTimestepIndex-1,freq).eq.0) then

            do jj = 1, ntspts

                if (numpe > 1) then

                    soln = varts(jj)
                    asoln = abs(soln)

                    !     if(jj.eq.24) then
                    !     write(*,*) soln
                    !     write(*,*)"and..."
                    !     endif

                    if (asoln.ne.zero) then
                        sgn = soln/asoln
                    else
                        sgn = 1
                    endif

                    call MPI_ALLREDUCE ( asoln, asolng, 1,  &
                    MPI_DOUBLE_PRECISION, MPI_MAX, &
                    INEWCOMM,ierr)
                    varts(jj) = sgn * asolng

                endif

                if (myrank.eq.zero) then
                    ifile = 1000+jj
                    write(ifile,555) varts(jj)
                    call flush(ifile)
                endif

            enddo


            varts = zero  ! reset the array for next step

555         format(e18.11)

        endif

    endif


    !
    !.... update and the aerodynamic forces
    !
    !call forces ( yold,  ilwork )

!    if((irscale.ge.0).or.(itwmod.gt.0))  &
!    call getvel (yold,     ilwork, iBC, &
!    nsons,    ifath, velbar)

!    if((irscale.ge.0).and.(myrank.eq.master)) then
!        call genscale(yold,       x,       iper,  &
!        iBC,     ifath,   velbar, &
!        nsons)
!    endif
    !
    !....  print out results.
    !
!    ntoutv=max(ntout,100)   ! velb is not needed so often
!    if ((irs .ge. 1) .and. (mod(currentTimestepIndex, ntout) .eq. 0)) then
!        if( (mod(currentTimestepIndex, ntoutv) .eq. 0) .and. &
!        ((irscale.ge.0).or.(itwmod.gt.0) .or.  &
!        ((nsonmax.eq.1).and.(iLES.gt.0)))) &
!        call rwvelb  ('out ',  velbar  ,ifail)
!    endif
    !
    !.... end of the NSTEP loop
    !
    !
    !.... -------------------> error calculation  <-----------------
    !
    if(ierrcalc.eq.1 .or. ioybar.eq.1) then
        !$$$c
        !$$$c compute average
        !$$$c
        !$$$               tfact=one/istep
        !$$$               ybar =tfact*yold + (one-tfact)*ybar

        ! compute average
        ! ybar(:,1) - ybar(:,3) is average velocity components
        ! ybar(:,4) is average pressure
        ! ybar(:,5) is average speed
        ! averaging procedure justified only for identical time step sizes
        ! istep is number of time step
        !
        tfact=one/istep

        ! ybar to contain the averaged ((u,v,w),p)-field
        ! and speed average, i.e sqrt(u^2+v^2+w^2)

        ybar(:,1) = tfact*yold(:,1) + (one-tfact)*ybar(:,1)
        ybar(:,2) = tfact*yold(:,2) + (one-tfact)*ybar(:,2)
        ybar(:,3) = tfact*yold(:,3) + (one-tfact)*ybar(:,3)
        ybar(:,4) = tfact*yold(:,4) + (one-tfact)*ybar(:,4)
        !
        dummyVar  = sqrt(yold(:,1)**2+yold(:,2)**2+yold(:,3)**2)

        if (istep .eq. 1) then
            ybar(:,5) = dummyVar
        else
            ybar(:,5) = tfact*dummyVar + (one-tfact)*ybar(:,5)
        endif
        !
        ! compute rms
        !
        rerr(:, 7)=rerr(:, 7)+(yold(:,1)-ybar(:,1))**2
        rerr(:, 8)=rerr(:, 8)+(yold(:,2)-ybar(:,2))**2
        rerr(:, 9)=rerr(:, 9)+(yold(:,3)-ybar(:,3))**2
        rerr(:,10)=rerr(:,10)+(yold(:,4)-ybar(:,4))**2
    endif

    call callCPPUpdateBoundaryConditionControlSystems()

 !if(istop.eq.1000) exit ! stop when delta small (see rstatic)

!2000  continue



end subroutine itrdrv_iter_finalize


!
!.... ---------------------->  Post Processing  <----------------------
!
subroutine itrdrv_finalize() bind(C, name="itrdrv_finalize")

    use iso_c_binding
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl
    use phcommonvars
    use itrDrvVars
    use debuggingTools

    use multidomain

    ! For deallocation:
    use readarrays, only: nBC, BCinp, iBCtmp
    use turbSA, only: wnrm, otwn
    use periodicity, only: rcount

    implicit none
    !IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    integer i,ku

    integer iqoldsiz
    integer itmp
    integer isize
    integer nitems
    integer lesid
    integer ifail

    real*8    tf(nshg,ndof)
    character*5  cname

    character*20    fname2,fmt2,fnamer2 ! for appending ybar
    integer         iarray(50) ! integers for headers

    !.... print out the last step
    !
    if ((irs .ge. 1) .and. ((mod(currentTimestepIndex, ntout) .ne. 0) .or. &
    (nstp .eq. 0))) then
!        if( &
!        ((irscale.ge.0).or.(itwmod.gt.0) .or.  &
!        ((nsonmax.eq.1).and.iLES.gt.0))) &
!        call rwvelb  ('out ',  velbar  ,ifail)
        write(*,*) "OUT2 max value of y: ", maxval(y)
        write(*,*) "OUT2 min value of y: ", minval(y)
        call restar ('out ',  yold  ,ac)
        if(ideformwall.eq.1) then
            call write_displ(myrank, currentTimestepIndex, nshg, 3, u, uref )
            if (imeasdist.eq.1) then
                call ElmDist(u,x,xdist,xdnv,df_fem)
                call write_distl(myrank, currentTimestepIndex, nshg,1,xdist)
            end if
        end if
    endif


    lesId   = numeqns(1)
    IF (memLSflag .NE. 1) THEN
        call saveLesRestart( lesId,  aperm , nshg, myrank, currentTimestepIndex, nPermDims )
    ENDIF


    if(ierrcalc.eq.1) then

        !
        !.....smooth the error indicators
        !
        do i=1,ierrsmooth
            call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
        end do
        !
        !.... open the output file
        !
        iqoldsiz=nshg*ndof*2
        call write_error(myrank, currentTimestepIndex, nshg, 10, rerr )


    endif

    if(ioybar.eq.1) then

        itmp = 1
        if (currentTimestepIndex .gt. 0) itmp = int(log10(float(currentTimestepIndex)))+1
        write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
        write (fname2,fmt2) currentTimestepIndex

        fname2 = trim(fname2) // cname(myrank+1)
        !
        !.... open  files
        !
        call openfile(  fname2//c_null_char,  c_char_"append?"//c_null_char, irstin )

        fnamer2 = 'ybar'
        isize = nshg*5
        nitems = 3
        iarray(1) = nshg
        iarray(2) = 5
        iarray(3) = currentTimestepIndex
        call writeheader(irstin, fnamer2//c_null_char,iarray, nitems, isize, &
        c_char_"double"//c_null_char, iotype )

        nitems = nshg*5
        call writedatablock(irstin, fnamer2//c_null_char,ybar, nitems, &
        c_char_"double"//c_null_char, iotype)

        call closefile( irstin, c_char_"append" )

    endif

    if ( ( ihessian .eq. 1 ) .and. ( numpe < 2 )  )then

        uhess = zero
        gradu = zero
        tf = zero

        do ku=1,nshg
            !           tf(ku,1) = x(ku,1)**2+2*x(ku,1)*x(ku,2)
            tf(ku,1) = x(ku,1)**3
        end do

        call hessian( yold, x,     shp,  shgl,   iBC,  &
        shpb, shglb, iper, ilwork, uhess, gradu )

        call write_hessian( uhess, gradu, nshg )

    endif
    !
    !.... close history and aerodynamic forces files
    !
    if (myrank .eq. master) then
        close (ihist)
        close (iforce)
    endif

    IF (memLSflag .eq. 1) THEN
        call memLS_LS_FREE(memLS_ls)
        call memLS_COMMU_FREE(communicator)
        call memLS_LHS_FREE(memLS_lhs)
        if (iheart .gt. 0) then
            call memLS_COMMU_FREE(communicator_s)
            call memLS_LHS_FREE(memLS_lhs_s)
        endif
    endif

    if (myrank .eq. master) then
        close(76)
    endif


    call deallocate_arrays_fortran()


5   format(1X,F15.10,3X,F15.10,3X,F15.10,3X,F15.10)
444 format(6(2x,e14.7))
    
    return

end subroutine itrdrv_finalize

subroutine deallocate_arrays() bind(C,name="deallocate_arrays")
    implicit none
    call deallocate_arrays_fortran()
end subroutine deallocate_arrays

subroutine deallocate_arrays_fortran()
    use iso_c_binding
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    !use grcrbc ! Nan rcr

    use convolCORFlow !for Coronary bc
    use incpBC        !for INCP bc
    use calcFlowPressure !to save history of flow and pressure of bc surfaces
    use LagrangeMultipliers
    use deformableWall
    use ResidualControl
    use phcommonvars
    use itrDrvVars
    use debuggingTools

    use multidomain

    ! For deallocation:
    use readarrays, only: nBC, BCinp, iBCtmp
    use turbSA, only: wnrm, otwn
    use periodicity, only: rcount

    implicit none

    if(allocated(FlowHist)) then
        deallocate(FlowHist)
    endif

    if(allocated(PressHist)) then
        deallocate(PressHist)
    endif

    ! if (associated(mBET%p)) then
    !     nullify(mBET%p)
    ! endif
    if (allocated(nBC)) then
        deallocate(nBC)
    endif
    if (allocated(wnrm)) then
        deallocate(wnrm)
    endif
    if (allocated(otwn)) then
        deallocate(otwn)
    endif
    if (allocated(iper)) then
        deallocate(iper)
    endif
    if (associated(mWNodes%p)) then
        deallocate(mWNodes%p)
        nullify(mWNodes%p)
    endif
    if (associated(mWNodes_gtlmap%p)) then
        deallocate(mWNodes_gtlmap%p)
        nullify(mWNodes_gtlmap%p)
    endif
    if(allocated(gNodes)) then
        deallocate(gNodes)
    endif
    if(allocated(sV)) then
        deallocate(sV)
    endif
    if(allocated(gNodes_s)) then
        deallocate(gNodes_s)
    endif
    if(allocated(sV_s)) then
        deallocate(sV_s)
    endif
    if (allocated(aperm)) then
        deallocate(aperm)
    endif
    if (allocated(atemp)) then
        deallocate(atemp)
    endif
    if (allocated(lhsP)) then
        deallocate(lhsP)
    endif
    if (allocated(lhsK)) then
        deallocate(lhsK)
    endif
    if(allocated(apermS)) then
        deallocate(apermS)
    endif
    if(allocated(atempS)) then
        deallocate(atempS)
    endif    
    if(allocated(lhsS)) then
        deallocate(lhsS)
    endif
    if(allocated(iBCd)) then
        deallocate(iBCd)
    endif
    if(allocated(iBCs)) then
        deallocate(iBCs)
    endif    
    ! if(allocated(gNodes_tmp)) then
    !     deallocate(gNodes_tmp)
    ! endif
    ! if(allocated(sV_tmp)) then
    !     deallocate(sV_tmp)
    ! endif
    if(allocated(rcount)) then
        deallocate(rcount)
    endif
    if(allocated(NABI)) then
        deallocate(NABI)
    endif
    if(allocated(NASC)) then
        deallocate(NASC)
    endif
    if(allocated(ndsurf)) then
        deallocate(ndsurf)
    endif
    if(Lagrange .gt. 0) then
        if(allocated(PNABI)) then
            deallocate(PNABI)
        endif
        if(allocated(NANBIJ)) then
            deallocate(NANBIJ)
        endif
    endif
    ! if(allocated(tmpshpb)) then
    !     deallocate(tmpshpb)
    ! endif
    ! if(allocated(Res4)) then
    !     deallocate(Res4)
    ! endif
    if(allocated(ilwork)) then
        deallocate(ilwork)
    endif
    ! if(allocated(ilworkread)) then
    !     deallocate(ilworkread)
    ! endif
    ! if(allocated(xread)) then
    !     deallocate(xread)
    ! endif
    ! if(allocated(ntagread)) then
    !     deallocate(ntagread)
    ! endif
    ! if(allocated(nBC)) then
    !     deallocate(nBC)
    ! endif
    ! if(allocated(nBCread)) then
    !     deallocate(nBCread)
    ! endif
    if(allocated(iBCtmp)) then
        deallocate(iBCtmp)
    endif
    ! if(allocated(iBCtmpread)) then
    !     deallocate(iBCtmpread)
    ! endif
    if(allocated(BCinp)) then
        deallocate(BCinp)
    endif
    ! if(allocated(BCinpread)) then
    !     deallocate(BCinpread)
    ! endif
    ! if(allocated(iperread)) then
    !     deallocate(iperread)
    ! endif
    if(allocated(inodesuniq)) then
        deallocate(inodesuniq) ! Valgrind thinks this is already free'd by now
    endif
    if(allocated(ilinobsfunc_sol)) then
        deallocate(ilinobsfunc_sol)
    endif
    if(allocated(ilinobsfunc_acc)) then
        deallocate(ilinobsfunc_acc)
    endif
    if(allocated(ilinobsfunc_disp)) then
        deallocate(ilinobsfunc_disp)
    endif
    if(allocated(obsfunc_dist)) then
        deallocate(obsfunc_dist)
    endif
    !
    ! if(allocated(acread)) then
    !     deallocate(acread)
    ! endif
    ! if(allocated(uread)) then
    !     deallocate(uread)
    ! endif
    if (allocated(iBC_original)) then
        deallocate(iBC_original)
    endif
    call destroyGlobalArrays()

    call tearDownMultidomain()

end subroutine deallocate_arrays_fortran


!
! *****************************************************************************
! *** initialise reduced order model with flows and pressures from 3d model ***
! *****************************************************************************
!

subroutine initreducedordermodel(y, rom, asciiformat)

    use multidomain
    use boundarymodule, only: area, integrScalar, GetFlowQ
    use phcommonvars

    real*8 :: y(nshg, ndof)
    type(reducedorder) :: rom
    integer :: nsurf, i
    integer :: srflist(0:MAXSURF)
    real*8  :: integpress(0:MAXSURF)
    real*8  :: currpress(0:MAXSURF)
    real*8  :: surfarea(0:MAXSURF)      
    real*8  :: currflow(0:MAXSURF)
    real*8  :: unity(nshg)
    character(len=*) :: asciiformat      

!   get number of surfaces 
    nsurf = rom%getsurfnum()      

!   get surface list in 0:MAXSURF array
    srflist = rom%getsurfids()

! !   calculate area by integrating unity over the surface 
! !   unity has dimensions of nshg = total number of nodes
!     unity(:) = one 
!     call integrScalar(surfarea,unity,srflist,nsurf) 
!     call rom%setarea(nsurf,surfarea)

    ! get surface area from boundarymodule 
    do i = 1, nsurf
        surfarea(i) = area(srflist(i))
        write(*,*) i,' ',srflist(i),' ',surfarea(i)
    end do

    ! set reduced order model
    call rom%setarea(nsurf,surfarea)


    if (currentTimestepIndex .eq. zero) then

!       integrate flow field on surface in normal direction
        call GetFlowQ(currflow,y(:,1:3),srflist,nsurf)

!       integrate pressure
        call integrScalar(integpress,y(:,4),srflist,nsurf)

!       get area and divide integrate pressure
        currpress(1:nsurf) = integpress(1:nsurf)/surfarea(1:nsurf)

!       set flows and pressure in reduced order model at step n
        call rom%setflow_n(nsurf,currflow) 
        call rom%setpressure_n(nsurf,currpress)           

    elseif (currentTimestepIndex .gt. zero) then
       
!       hack for heart model !!!
        call rom%loadflowfile(currentTimestepIndex,asciiformat)
        call rom%loadpressurefile(currentTimestepIndex,asciiformat) !! currentTimestepIndex index starts at zero  

    end if

    return

end subroutine initreducedordermodel

! ***********************************************************************
! *** initialise multidomain container with areas, pressure and flows ***
! ***********************************************************************


subroutine initmultidomaincontainer(y,mdc)

    use multidomain, only: multidomaincontainer
    use boundarymodule, only: area, integrScalar, GetFlowQ
    use phcommonvars

    real*8 :: y(nshg, ndof)
    type(multidomaincontainer) :: mdc
    integer :: nsurf 
    integer :: srflist(0:MAXSURF)
    real*8 :: integpress(0:MAXSURF)
    real*8 :: currpress(0:MAXSURF)
    real*8 :: surfarea(0:MAXSURF)      
    real*8 :: currflow(0:MAXSURF)
    real*8 :: unity(nshg)

    ! get number of surfaces 
    nsurf = mdc%getsurfnum()      

    if (nsurf .gt. int(0)) then
    
        ! get surface list in 0:MAXSURF array
        srflist = mdc%getsurfids()

        ! get surface area from boundarymodule 
        do i = 1, nsurf
           surfarea(i) = area(srflist(i))
           write(*,*) 'mdc area: ', i,' ',srflist(i),' ',surfarea(i)
        end do


        ! calculate area by integrating unity over the surface 
        ! unity has dimensions of nshg = total number of nodes
        unity(:) = one 
        call integrScalar(surfarea,unity,srflist,nsurf) 

        do i = 1, nsurf
           write(*,*) 'mdc area: ', i,' ',srflist(i),' ',surfarea(i)
        end do


        call mdc%setarea(nsurf,surfarea)
         
        ! integrate flow field on surface in normal direction
        call GetFlowQ(currflow,y(:,1:3),srflist,nsurf)

        ! integrate pressure
        call integrScalar(integpress,y(:,4),srflist,nsurf)
        ! get area and divide integrate pressure
        currpress(1:nsurf) = integpress(1:nsurf)/surfarea(1:nsurf)

        ! set flows and pressure in reduced order model at step n
        call mdc%setflow(nsurf,currflow)          
        call mdc%setpressure(nsurf,currpress)

      end if 

      return
      end subroutine initmultidomaincontainer

! ********************************************************
! *** subroutine to reset flow at time step n using the fluid solution y
!     added because of the filter modifies the state at the timestep n
!     KDL NAN 22/08/14
!
! *** upgraded and renamed to support a full step-back of all
!     types of boundary conditions
!     CA 2/2/16
! ********************************************************

subroutine resetBoundaryConditionStateForKalmanFilter(y, numberOfSurfaces, listOfSurfaces)

    use multidomain
    use boundarymodule, only: GetFlowQ, integrScalar
    use phcommonvars
    use cpp_interface

    real*8, intent(in) :: y(nshg, ndof)
    integer, intent(in) :: numberOfSurfaces
    integer, intent(in) :: listOfSurfaces(0:MAXSURF)
    real*8 :: surfaceFlows(0:MAXSURF)
    real*8 :: POnly(nshg)
    real*8 :: surfacePressures(0:MAXSURF)
    real*8 :: flowSurfaceArea(0:MAXSURF)
    integer :: indexAmongstCurrentSurfaceType
    integer :: surfaceIndex

    ! integrate flow field on surface in normal direction
    call GetFlowQ(surfaceFlows, y(:,1:3), listOfSurfaces, numberOfSurfaces)

    POnly(:)=y(:,4) ! pressure
    call integrScalar(surfacePressures, POnly, listOfSurfaces, numberOfSurfaces) !get surface pressure integral
    call integrScalar(flowSurfaceArea, POnly*0 + 1.0, listOfSurfaces, numberOfSurfaces)

    ! compute the area-averaged pressure over the surface
    do indexAmongstCurrentSurfaceType = 1, numberOfSurfaces
        surfacePressures(indexAmongstCurrentSurfaceType) = surfacePressures(indexAmongstCurrentSurfaceType) / flowSurfaceArea(indexAmongstCurrentSurfaceType) !(50.265) ! just a hack until i sort out the areas for netlists... one model only, and assuming radius is 4 units
    end do

    ! tell the boundary condition in C++ to do the actual resetting to the state the Kalman filter wants us to begin from:
    do indexAmongstCurrentSurfaceType = 1, numberOfSurfaces
        surfaceIndex = listOfSurfaces(indexAmongstCurrentSurfaceType)
        call callCPPResetStateUsingKalmanFilteredEstimate(surfaceFlows(indexAmongstCurrentSurfaceType), surfacePressures(indexAmongstCurrentSurfaceType), surfaceIndex, currentTimestepIndex)
    end do

end subroutine resetBoundaryConditionStateForKalmanFilter


!      
! ********************************************************
! *** update pressure and flows in reduced order model ***
! ********************************************************
!
subroutine updreducedordermodel(y,rom,varchar)

    use multidomain
    use iso_c_binding
    use cpp_interface
    use boundarymodule, only: GetFlowQ, integrScalar
    use phcommonvars

    real*8 :: y(nshg, ndof)
    type(reducedorder) :: rom
    integer :: nsurf 
    integer :: srflist(0:MAXSURF)
    real*8 :: integpress(0:MAXSURF)
    real*8 :: currpress(0:MAXSURF)
    real*8 :: surfarea(0:MAXSURF)      
    real*8 :: currflow(0:MAXSURF)
    character(len=*) :: varchar 
    character(len=*), parameter :: updchar = 'update' ! n+1
    character(len=*), parameter :: updchar_flow = 'update-flow-only' ! n+1    
    character(len=*), parameter :: solchar = 'solve'  ! n+alf

    ! get number of surfaces 
    nsurf = rom%getsurfnum()

    ! get surface list in 0:MAXSURF array
    srflist = rom%getsurfids()

    ! integrate flow field on surface in normal direction
    call GetFlowQ(currflow,y(:,1:3),srflist,nsurf)

    ! if update at t = t_{n+1}
    if (varchar .eq. updchar) then
  
        ! set flows and pressure at n+1 step to step n in reduced order
        ! model before moving onto next step
        call rom%setflow_n(nsurf,currflow)
        call rom%setflow_n1(nsurf,currflow)

        ! if(numGRCRSrfs.gt.0) then
            !call grcrbc_UpdateInternalState(y)

            ! this is here because it is particle dependent as the WK pressure is part of the state
            ! the pressure_n is a set via the pointer set previously 
            ! if (nrcractive) then ! commented this IF whilst doing the C++ additions - check it's necessary! \todo
                ! call callCPPUpdateAllRCRS_setflow_n(c_loc(currflow(1)))
                ! call callCPPUpdateAllRCRS_setflow_n1(c_loc(currflow(1)))
            ! end if

        ! endif      

        ! if flowsolver only then update pressure the normal way
        ! if estimator, do nothing as pressure is set by the filter from the average of all the particles
        
        ! check if pressure is required to be passed to the reduced order model
        ! for e.g. the heart model when the aortic valve is shut
              
        if (rom%ispressureupdate() .eq. 1) then
      
            ! integrate pressure
            call integrScalar(integpress,y(:,4),srflist,nsurf)
    
            ! get area and divide integrate pressure
            surfarea = rom%getarea()
            currpress(1:nsurf) = integpress(1:nsurf)/surfarea(1:nsurf)
    
            ! update model with integrated value
            call rom%updpressure_n1_withvalue(currpress)


            ! else update pressure from pressure flow relationship
        else

            ! update pressure from pressure/flow from saved flow in reduced model
            ! sets it to pressure_n
            call rom%updpressure_n1_withflow()     
    
        end if

        ! These old calls are on specific classes, so for the c++ implementation
        ! I'm going to move them outside this function.

        ! if (rom%classNameString .eq. 'controlledCoronaryModel') then
        !   !Update coronary internal pressures, now we're done with this timestep:
        !   call controlledCoronarySurfaces%updateLPN_coronary(currentTimestepIndex)
        ! else if (rom%classNameString .eq. 'netlistLPN') then
        !   call netlistLPNSurfaces%updateLPN_netlistLPN()
        ! end if

    else if (varchar .eq. updchar_flow) then

        ! here we only update the flows, this is relevant for the filter 
        ! the filter sets the pressure in the reduced order model via the pointer
        call rom%setflow_n(nsurf,currflow) 

    else if (varchar .eq. solchar) then

        ! set flows at current n+alf step 
        call rom%setflow_n1(nsurf,currflow) 

    end if


end subroutine updreducedordermodel
