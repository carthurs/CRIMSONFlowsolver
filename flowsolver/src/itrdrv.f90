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
    !
    !.... For linear solver Library
    !
    integer eqnType, prjFlag, presPrjFlag, verbose  ! init
    character*1024    servername ! init
    character*20        license_f_name ! init


!--------------------------------------------------------------------
!   Setting up memLS

    IF (memLSFlag .EQ. 1) THEN
        CALL memLS_LS_CREATE(memLS_ls, LS_TYPE_NS, dimKry=Kspace,relTol=epstol(8), &
                             relTolIn=(/epstol(1),epstol(7)/), maxItr=nPrjs, maxItrIn=(/nGMRES,maxIters/))

        CALL memLS_COMMU_CREATE(communicator, MPI_COMM_WORLD)

        IF (numpe .GT. 1) THEN
            WRITE(fileName,*) myrank
            fileName = "ltg.dat."//ADJUSTL(TRIM(fileName))
            OPEN(1,FILE=fileName)
            READ(1,*) gnNo
            READ(1,*) nNo
            ALLOCATE(ltg(nNo))
            READ(1,*) ltg
            CLOSE(1)
        ELSE
            gnNo = nshg
            nNo = nshg
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
            IF  (ipvsq .GE. 2) THEN
                memLS_nFaces = 1 + numResistSrfs + numImpSrfs + numRCRSrfs + numGRCRSrfs
            ELSE
                memLS_nFaces = 1
            END IF

            CALL memLS_LHS_CREATE(memLS_lhs, communicator, gnNo, nNo, nnz_tot, ltg, colm, rowp, memLS_nFaces)

            faIn = 1
            facenNo = 0
            DO i=1, nshg
                IF (IBITS(iBC(i),3,3) .NE. 0)  facenNo = facenNo + 1
            END DO
            ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
            sV = 0D0
            j = 0
            DO i=1, nshg
                IF (IBITS(iBC(i),3,3) .NE. 0) THEN
                    j = j + 1
                    gNodes(j) = i
                    IF (.NOT.BTEST(iBC(i),3)) sV(1,j) = 1D0
                    IF (.NOT.BTEST(iBC(i),4)) sV(2,j) = 1D0
                    IF (.NOT.BTEST(iBC(i),5)) sV(3,j) = 1D0
                END IF
            END DO
            CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo, nsd, BC_TYPE_Dir, gNodes, sV)

            IF  (ipvsq .GE. 2) THEN
                DO k = 1, numResistSrfs
                    faIn = faIn + 1
                    CALL AddNeumannBCTomemLS(nsrflistResist(k), faIn)
                END DO
                DO k = 1, numImpSrfs
                    faIn = faIn + 1
                    CALL AddNeumannBCTomemLS(nsrflistImp(k), faIn)
                END DO
                DO k = 1, numRCRSrfs
                    faIn = faIn + 1
                    CALL AddNeumannBCTomemLS(nsrflistRCR(k), faIn)
                END DO
                DO k = 1, numGRCRSrfs
                    faIn = faIn + 1
                    CALL AddNeumannBCTomemLS(nsrflistGRCR(k), faIn)
                END DO
            END IF
        ELSE
!--------------------------------------------------------------------


            call myfLesNew( lesId,          41994, &
                            eqnType, &
                            nDofs,          minIters,       maxIters, &
                            nKvecs,         prjFlag,        nPrjs, &
                            presPrjFlag,    nPresPrjs,      epstol(1), &
                            prestol,        verbose,        statsflow, &
                            nPermDims,      nTmpDims,      servername  )

        END IF

        allocate (aperm(nshg,nPermDims))
        allocate (atemp(nshg,nTmpDims))
        allocate (lhsP(4,nnz_tot))
        allocate (lhsK(9,nnz_tot))

        call readLesRestart( lesId,  aperm, nshg, myrank, lstep, nPermDims )

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
        allocate (apermS(nshg,nPermDimsS,nsclrsol))
        allocate (atempS(nshg,nTmpDimsS))  !they can all share this
        allocate (lhsS(nnz_tot,nsclrsol))
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
    lstep0 = lstep+1
    !
    !.... make two copies of iBC array to use a switch for INCP boundary condition
    !
    if (incp .gt. zero) then
        allocate(iBCd(nshg))
        iBCd = iBC
        if (numINCPSrfs .eq. one) then
            allocate(iBCs(1,nshg))
            iBCs(1,:) = iBC
            do i=1, nshg
                if(ndsurf(i).eq.nsrflistINCP(1)) then
                    iBCs(1,i)=0
                endif
            enddo
        elseif (numINCPSrfs .eq. two) then
            allocate(iBCs(3,nshg))
            do k=1, 3
                iBCs(k,:) = iBC
            enddo

            do i=1, nshg
                if(ndsurf(i).eq.nsrflistINCP(1)) then
                    iBCs(1,i)=0
                    iBCs(2,i)=0
                elseif(ndsurf(i).eq.nsrflistINCP(2)) then
                    iBCs(1,i)=0
                    iBCs(3,i)=0
                endif
            enddo
        endif
    endif
    !
    !...initialize the initial condition for INCP BC
    !
    if(numINCPSrfs.gt.zero) then
        call calcINCPic(Delt(1), y, nsrflistINCP, numINCPSrfs)
    endif

    if (incp.gt.zero) then
        if (numINCPSrfs .eq. one) then
            if(Qaorta(lstep+1,1) .lt. zero) then
                iBC = iBCs(1,:)   !systole
                inactive(1)=1
                INCPSwitch = 1
            endif
        elseif (numINCPSrfs .eq. two) then
            if(Qaorta(lstep+1,1) .lt. zero .and.  &
            Qaorta(lstep+1,2) .lt. zero) then
                iBC = iBCs(1,:)   !systole
                inactive(1:2)=1
                INCPSwitch = 1
            endif
        endif
    endif

    !
    !
    !.... satisfy the boundary conditions
    !
    call itrBC (y, ac,  iBC, BC, iper, ilwork)
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

    call itrSetup ( y, acold )
    
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
        
        if (nrcractive) then
            nrcr = nrcrconstructor(numGRCRSrfs,nsrflistGRCR)
            call initreducedordermodel(y, nrcr, 'legacy')
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

    CONTAINS

    SUBROUTINE AddNeumannBCTomemLS(srfID, faIn)

    INTEGER, INTENT(IN) :: srfID, faIn

    INTEGER facenNo, i, j

    facenNo = 0
    DO i = 1, nshg
     IF (srfID .EQ. ndsurf(i)) THEN
        facenNo = facenNo + 1
     END IF
    END DO
    IF (ALLOCATED(gNodes)) DEALLOCATE(gNodes, sV)
    ALLOCATE(gNodes(facenNo), sV(nsd,facenNo))
    sV = 0D0
    j = 0
    DO i = 1, nshg
     IF (srfID .EQ. ndsurf(i)) THEN
        j = j + 1
        gNodes(j) = i
        sV(:,j) = NABI(i,1:3)
     END IF
    END DO
    CALL memLS_BC_CREATE(memLS_lhs, faIn, facenNo, nsd, BC_TYPE_Neu, gNodes, sV)

    RETURN
    END SUBROUTINE AddNeumannBCTomemLS





end subroutine


!
!.... ---------------> initialize the time step <---------------
!
subroutine itrdrv_iter_init() bind(C, name="itrdrv_iter_init")

    use iso_c_binding
    use shapeTable
    use globalArrays
    use pvsQbi     !gives us splag (the spmass at the end of this run
    use specialBC !gives us itvn
    use timedata   !allows collection of time series
    use convolImpFlow !for Imp bc
    use convolRCRFlow !for RCR bc
    use convolTRCRFlow !for time-varying RCR bc

    use multidomain, only: nrcr, nrcractive
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

    if (incp.gt.zero) then  ! works only when there is one "INCP" srf
        if (numINCPSrfs .eq. one) then
            if(PLV(lstep+1,1) .gt. Paorta(lstep+1,1)) then
                iBC = iBCs(1,:)   !systole
                inactive(1)=1
                INCPSwitch = 1
            elseif (INCPSwitch .gt. 0 .and.  &
            Qaorta(lstep+1,1) .le. zero) then
                iBC = iBCs(1,:)   !systole
                inactive(1)=1
                INCPSwitch = 1
            else
                iBC=iBCd   !diastole
                inactive(1)=nsrflistINCP(1)
                INCPSwitch = 0
            endif
        elseif (numINCPSrfs .eq. two) then
            if(PLV(lstep+1,1) .gt. Paorta(lstep+1,1) .and.  &
            PLV(lstep+1,2) .gt. Paorta(lstep+1,2)) then
                iBC = iBCs(1,:)   !systole
                inactive(1:2)=1
                INCPSwitch = 1
            elseif(PLV(lstep+1,1) .gt. Paorta(lstep+1,1) .and.  &
            PLV(lstep+1,2) .le. Paorta(lstep+1,2)) then
                iBC = iBCs(2,:)   !systole
                inactive(1)=1
                inactive(2)=nsrflistINCP(2)
                INCPSwitch = 1
            elseif(PLV(lstep+1,1) .le. Paorta(lstep+1,1) .and.  &
            PLV(lstep+1,2) .gt. Paorta(lstep+1,2)) then
                iBC = iBCs(3,:)   !systole
                inactive(2)=1
                inactive(1)=nsrflistINCP(1)
                INCPSwitch = 1
            elseif (INCPSwitch.gt.0 .and. Qaorta(lstep+1,1).le.zero  &
            .and. Qaorta(lstep+1,2) .le. zero) then
                iBC = iBCs(1,:)   !systole
                inactive(1:2)=1
                INCPSwitch = 1
            elseif (INCPSwitch.gt.0 .and. Qaorta(lstep+1,1).le.zero  &
            .and. Qaorta(lstep+1,2) .gt. zero) then
                iBC = iBCs(2,:)   !systole
                inactive(1)=1
                inactive(2)=nsrflistINCP(2)
                INCPSwitch = 1
            elseif (INCPSwitch.gt.0 .and. Qaorta(lstep+1,1).gt.zero  &
            .and. Qaorta(lstep+1,2) .le. zero) then
                iBC = iBCs(3,:)   !systole
                inactive(2)=1
                inactive(1)=nsrflistINCP(1)
                INCPSwitch = 1
            else
                iBC=iBCd   !diastole
                inactive(1)=nsrflistINCP(1)
                INCPSwitch = 0
            endif
        endif
    endif

    xi=(istep+1)*1.0/nstp
    datmat(1,2,1)=rmub*(1.0-xi)+xi*rmue
    !            write(*,*) "current mol. visc = ", datmat(1,2,1)
    !.... if we have time varying boundary conditions update the values of BC.
    !     these will be for time step n+1 so use lstep+1
    !
    if(itvn.gt.0) call BCint((lstep+1)*Delt(1), shp, shgl,  &
    shpb, shglb, x, BC, iBC)

    !
    ! ... calc the pressure contribution that depends on the history for the imp BC
    !
    if(numImpSrfs.gt.0) call pHist(poldImp,QHistImp,ImpConvCoef, &
    ntimeptpT,numImpSrfs)
    !
    ! ... calc the pressure contribution that depends on the history for the RCR BC
    !
    if(numRCRSrfs.gt.0) then
        call CalcHopRCR (Delt(itseq), lstep, numRCRSrfs)
        call CalcRCRConvCoef(lstep,numRCRSrfs)
        call pHist(poldRCR,QHistRCR,RCRConvCoef, &
        nstep+nptsRCR,numRCRSrfs)
    endif
    !
    ! ... compute coefficients required for time-varying RCR BCs
    !
    if(numTRCRSrfs.gt.0) then
        call CalcTRCRConvCoef(lstep,numTRCRSrfs)
        call CalcHopTRCR (Delt(itseq), lstep, numTRCRSrfs)
        call pHist(poldTRCR,QHistTRCR,TRCRConvCoef, &
        nstep(1)+nptsTRCR,numTRCRSrfs)
    endif

    ! Nan rcr ----------------------------------
    if(numGRCRSrfs.gt.0) then
        !call grcrbc_ComputeImplicitCoefficients(lstep, yold)

        if (nrcractive) then
            ! first reset flow for the filter
            call reset_flow_n(yold, nrcr)
            ! calculate the implicit coefficients
            call nrcr%setimplicitcoeff(lstep)
        end if 

    endif
    ! ------------------------------------------

    !
    ! ... calc the pressure contribution that depends on the history for the Coronary BC
    !
    if(numCORSrfs.gt.0) then
        call CalcCORConvCoef(lstep,numCORSrfs)
        call pHist(poldCOR, QHistCOR, CORConvCoef, &
        nstep+nptsCOR,numCORSrfs)
        call CalcHopCOR (Delt(itseq), lstep, nsrflistCOR,  &
        numCORSrfs,yold)
    endif
    !
    !.... calculate the coefficients
    !
    if (numINCPSrfs .gt. 0) then
        call CalcINCPConvCoef(lstep, numINCPSrfs)
        call pHist(poldINCP, QHistINCP, &
        INCPConvCoef, nstep+nptsINCP, numINCPSrfs)
        call CalcINCPCoef(Delt(itseq), lstep, &
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

end subroutine

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

    use multidomain, only: nrcr, nrcractive
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

    integer j

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

    ! Nan rcr ----------------------------------
    ! update P,Q variables
    if(numGRCRSrfs.gt.0) then
        !call grcrbc_UpdateInternalState(y)

        if (nrcractive) then
            call updreducedordermodel(y,nrcr,'update')
        end if

    endif
    ! ------------------------------------------

    ! interface to compute distances to observed wall motion
!    if (imeasdist.eq.1) then
!        if (myrank.eq.zero) then
!            write(*,*) "computing distance to wall data surfaces (end of step)"
!        end if
!        call ElmDist(u,x,xdist,xdnv,df_fem)
!    end if

end subroutine

!
!.... ---------------> finalize the time step <---------------
!

subroutine itrdrv_iter_finalize() bind(C, name="itrdrv_iter_finalize")

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

    implicit none
    !IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    integer jj
    integer ifail

    call itrBC (yold, acold,  iBC,  BC,  iper,ilwork)

    istep = istep + 1
    lstep = lstep + 1

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
    ! ... write out the solution
    !
    if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
        call restar ('out ',  yold  ,ac)
        if(ideformwall.eq.1) then
            call write_displ(myrank, lstep, nshg, 3, uold, uref )
            if (imeasdist.eq.1) then
                call write_distl(myrank, lstep, nshg, 1, xdist ) ! should use nshg or numnp?
            end if
        end if
    endif
    !
    ! ... update the flow history for the INCP convolution
    !
    if(numINCPSrfs.gt.zero) then
        call UpdHistConv(y,nsrflistINCP,numINCPSrfs)
        call UpdHeartModel(Delt(itseq),y,nsrflistINCP,numINCPSrfs,lstep)
    endif
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
        call UpdHistConv(y,nsrflistRCR,numRCRSrfs) !uses lstep
        call UpdRCR(y,nsrflistRCR,numRCRSrfs)
    endif
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
        call UpdHistConv(y,nsrflistCOR,numCORSrfs) !uses lstep
        call UpdHistPlvConv(y,Delt(itseq),lstep,nsrflistCOR, numCORSrfs)
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

        if (mod(lstep-1,freq).eq.0) then

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
!    if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
!        if( (mod(lstep, ntoutv) .eq. 0) .and. &
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

 !if(istop.eq.1000) exit ! stop when delta small (see rstatic)

!2000  continue



end subroutine


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
    if ((irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or. &
    (nstp .eq. 0))) then
!        if( &
!        ((irscale.ge.0).or.(itwmod.gt.0) .or.  &
!        ((nsonmax.eq.1).and.iLES.gt.0))) &
!        call rwvelb  ('out ',  velbar  ,ifail)
        call restar ('out ',  yold  ,ac)
        if(ideformwall.eq.1) then
            call write_displ(myrank, lstep, nshg, 3, u, uref )
            if (imeasdist.eq.1) then
                call ElmDist(u,x,xdist,xdnv,df_fem)
                call write_distl(myrank, lstep, nshg,1,xdist)
            end if
        end if
    endif


    lesId   = numeqns(1)
    IF (memLSflag .NE. 1) THEN
        call saveLesRestart( lesId,  aperm , nshg, myrank, lstep, nPermDims )
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
        call write_error(myrank, lstep, nshg, 10, rerr )


    endif

    if(ioybar.eq.1) then

        itmp = 1
        if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
        write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
        write (fname2,fmt2) lstep

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
        iarray(3) = lstep
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
5   format(1X,F15.10,3X,F15.10,3X,F15.10,3X,F15.10)
444 format(6(2x,e14.7))

    return

end subroutine


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


    if (lstep .eq. zero) then

!       integrate flow field on surface in normal direction
        call GetFlowQ(currflow,y(:,1:3),srflist,nsurf)

!       integrate pressure
        call integrScalar(integpress,y(:,4),srflist,nsurf)

!       get area and divide integrate pressure
        currpress(1:nsurf) = integpress(1:nsurf)/surfarea(1:nsurf)

!       set flows and pressure in reduced order model at step n
        call rom%setflow_n(nsurf,currflow) 
        call rom%setpressure_n(nsurf,currpress)              

    elseif (lstep .gt. zero) then
       
!       hack for heart model !!!
        call rom%loadflowfile(lstep,asciiformat)
        call rom%loadpressurefile(lstep,asciiformat) !! lstep index starts at zero  

    end if

    return

end subroutine initreducedordermodel


! ********************************************************
! *** subroutine to reset flow at time step n using the fluid solution y
!     added because of the filter modifies the state at the timestep n
!     KDL NAN 22/08/14
! ********************************************************

subroutine reset_flow_n(y,rom)

    use multidomain
    use boundarymodule, only: GetFlowQ
    use phcommonvars

    real*8 :: y(nshg, ndof)
    type(reducedorder) :: rom
    integer :: nsurf 
    integer :: srflist(0:MAXSURF)
    real*8 :: currflow(0:MAXSURF)

    ! get number of surfaces 
    nsurf = rom%getsurfnum()

    ! get surface list in 0:MAXSURF array
    srflist = rom%getsurfids()

    ! integrate flow field on surface in normal direction
    call GetFlowQ(currflow,y(:,1:3),srflist,nsurf)

    ! reset flow at step n in reduced order model
    call rom%setflow_n(nsurf,currflow) 

end subroutine 


!      
! ********************************************************
! *** update pressure and flows in reduced order model ***
! ********************************************************
!
subroutine updreducedordermodel(y,rom,varchar)

    use multidomain
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

        ! check if pressure is required to be passed to the reduced order model
        ! for e.g. the heart model when the aortic valve is shut
          
        if (rom%ispressureupdate()) then
  
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
            call rom%updpressure_n1_withflow()
    
        end if

        ! if (rom%classNameString .eq. 'controlledCoronaryModel') then
        !   !Update coronary internal pressures, now we're done with this timestep:
        !   call controlledCoronarySurfaces%updateLPN_coronary(lstep)
        ! else if (rom%classNameString .eq. 'netlistLPN') then
        !   call netlistLPNSurfaces%updateLPN_netlistLPN()
        ! end if

    elseif (varchar .eq. solchar) then

        ! set flows at current n+alf step 
        call rom%setflow_n1(nsurf,currflow) 
    
    end if 

end subroutine 
