module grcrbc_internal

    implicit none

    real*8, allocatable :: parameters_RCR(:,:), parameters_Pdist(:,:,:) !inputs

    real*8, allocatable :: coeff_1_implicit(:),coeff_2_implicit(:)
    real*8, allocatable :: RCoverDt(:) ! (C * Rd / dt)
!    real*8, allocatable :: PDist_alpha(:)     ! Pd_n+1 distal pressure
!    real*8, allocatable :: PDist_current(:)  ! Pd_n distal pressure

    real*8, allocatable :: P_current(:)
    real*8, allocatable :: Q_current(:)

    real*8, allocatable :: SrfArea(:) ! area of the outlet face

    integer numtimepoints_Pdistmax
    integer, allocatable :: numtimepoints_Pdist(:)
    real*8  PDist ! one constant distal pressure everywhere

contains

    subroutine Initialize()

        use iso_c_binding
        use phcommonvars

        implicit none

        integer i,j,k,n
        integer numDataRCR

        ! open rcr parameter file to read values and allocate arrays
        open(unit=818, file='rcrt.dat',status='old')
        read (818,*) numtimepoints_Pdistmax

        allocate (numtimepoints_Pdist(numGRCRSrfs))
        allocate (parameters_Pdist(numtimepoints_Pdistmax,2,numGRCRSrfs))
        allocate (parameters_RCR(3,numGRCRSrfs))

        allocate (coeff_1_implicit(numGRCRSrfs))
        allocate (coeff_2_implicit(numGRCRSrfs))
        allocate (RCoverDt(numGRCRSrfs))
!        allocate (PDist_current(numGRCRSrfs))
!        allocate (PDist_alpha(numGRCRSrfs))

        allocate (P_current(numGRCRSrfs))
        allocate (Q_current(numGRCRSrfs))

        allocate (SrfArea(numGRCRSrfs))

        parameters_Pdist = 0
        parameters_RCR = 0

        do k=1,numGRCRSrfs
            read (818,*) numDataRCR
            numtimepoints_Pdist(k) = numDataRCR
            do j=1,3
                read(818,*) parameters_RCR(j,k) ! reads Rp,C,Rd in that order
            enddo
            do j=1,numtimepoints_Pdist(k)
                read(818,*) (parameters_Pdist(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
        enddo
        close(818)

        PDist = parameters_Pdist(1,2,1) ! read just the first value
        !write(*,*) 'PDIST: ', PDist

        call phgloballumpedparameterarrayassignpointer(c_loc(P_current), &
            c_loc(Q_current), c_loc(parameters_RCR), c_loc(PDist))

    end subroutine

    subroutine SetInternalState(y)

        use phcommonvars

        implicit none

        real*8 y(nshg,ndof)
        real*8 P_temp(nshg)        ! temporary array for nodal pressures

        real*8 PSrf_temp(0:MAXSURF)  ! temporary array for srf pressures
        real*8 QSrf_temp(0:MAXSURF)  ! temporary array for srf flows
        real*8 ASrf_temp(0:MAXSURF)  ! temporary array for srf area

        ! compute the initial integrated pressure at each face (initial P_n)

        P_temp = one ! temporary array set to ones in order to compute face area

        call integrScalar(ASrf_temp,P_temp,nsrflistGRCR,numGRCRsrfs) ! first compute face area
        SrfArea(:) = ASrf_temp(1:numGRCRSrfs)

        call integrScalar(PSrf_temp,y(:,4),nsrflistGRCR,numGRCRsrfs) ! get pressure integral
        P_current(:) = PSrf_temp(1:numGRCRSrfs) / SrfArea(:) ! compute average pressure over the pace

        call GetFlowQ(QSrf_temp,y(:,1:3),nsrflistGRCR,numGRCRsrfs) !get flow integral
        Q_current(:) = QSrf_temp(1:numGRCRSrfs)

        !write(*,*) "initial distal P ", P_current
        !write(*,*) "initial distal Q ", Q_current

    end subroutine

    ! interpolate distal pressure from values in the rcrt.dat file
    ! input: ctime
    ! output: Pdist

    ! this routine and others like it should be replaced with a more general interpolate routine

!    subroutine PdistInterpolate(ctime,Pdist)
!
!        use phcommonvars
!        implicit none
!
!        real*8  ctime, ptime, wr
!        integer j, k
!        integer nlast, nper
!        real*8  Pdist(0:MAXSURF)
!
!        do k=1,numGRCRSrfs
!            nlast=numtimepoints_Pdist(k)     ! number of time series to interpolate from
!            nper=ctime/parameters_Pdist(nlast,1,k) !number of periods completed to shift off
!            ptime = ctime-nper*parameters_Pdist(nlast,1,k)  ! now time in periodic domain
!
!            do j=2,nlast   !loop to find the interval that we are in
!
!                if(parameters_Pdist(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
!
!                    wr=(ptime-parameters_Pdist(j-1,1,k)) / ( parameters_Pdist(j,1,k)-parameters_Pdist(j-1,1,k) )
!
!                    Pdist(k)= parameters_Pdist(j-1,2,k)*(one-wr) + parameters_Pdist(j,2,k)*wr
!
!                    exit
!
!                endif
!
!            enddo
!        enddo
!        return
!
!    end subroutine

    ! update integrated pressure
    ! input: stepn, stepsize

    ! will need to re examine use of global variables

    subroutine ComputeImplicitCoefficients(stepn,y)

        use phcommonvars
        implicit none

        real*8 y(nshg,ndof)

        real*8 QSrf_temp(0:MAXSURF)
        real*8 Pdist_temp(0:MAXSURF)

        integer stepn

        real*8 Dt

        Dt = alfi*Delt(1)

!        call PdistInterpolate(Delt(1)*stepn,PDist_current(:))
!        call PdistInterpolate(Delt(1)*(stepn+alfi),Pdist_temp(:))

!        PDist_alpha = Pdist_temp(1:numGRCRSrfs);
        !write(*,*) Delt(1)*(stepn+alfi),PDist_alpha

        call GetFlowQ(QSrf_temp,y(:,1:3),nsrflistGRCR,numGRCRsrfs) !get flow integral
        Q_current(:) = QSrf_temp(1:numGRCRSrfs)

        RCoverDt(:) = parameters_RCR(2,:)*parameters_RCR(3,:) / ( Dt )

        coeff_1_implicit(:) = one / (one + RCoverDt(:)) * ( parameters_RCR(3,:)+parameters_RCR(1,:)*(one+RCoverDt(:)) )

!        coeff_2_implicit(:) = one / (one + RCoverDt(:)) * ( RCoverDt(:)*( P_current(:)+PDist_alpha(:)-PDist_current(:)-parameters_RCR(1,:)*Q_current(:) ) + PDist_alpha(:) )

        coeff_2_implicit(:) = one / (one + RCoverDt(:)) * ( RCoverDt(:)*( P_current(:)-parameters_RCR(1,:)*Q_current(:) ) + PDist )
        !coeff_2_implicit(:) = one / (one + RCoverDt(:)) * ( RCoverDt(:)*( P_current(:)-parameters_RCR(1,:)*Q_current(:) ) )

!        write(*,*) "rank ", myrank, " PDist_current ", PDist_current
!        write(*,*) "rank ", myrank, " PDist_current ", PDist_current

        !write(*,*) "rank ", myrank, " coeff_1_implicit ", coeff_1_implicit
        !write(*,*) "rank ", myrank, " coeff_2_implicit ", coeff_2_implicit

    end subroutine

    ! update pressure and flow variables
    ! input: stepn, stepsize

    subroutine UpdateInternalState(y)

        use phcommonvars

        implicit none

        real*8 y(nshg,ndof)

        real*8 QSrf_temp(0:MAXSURF)

        call GetFlowQ(QSrf_temp,y(:,1:3),nsrflistGRCR,numGRCRsrfs) !get flow integral
        Q_current(:) = QSrf_temp(1:numGRCRSrfs)

        !write(*,*) "rank ", myrank, " numGRCRsrfs ", numGRCRsrfs
        !write(*,*) "rank ", myrank, " distal Q ", QSrf_temp(1:numGRCRSrfs)

        P_current(:) =  coeff_1_implicit(:) * Q_current(:) + coeff_2_implicit(:)

        if (myrank .eq. 0) then
        !write(*,*) "rank ", myrank, " distal P ", P_current

        end if

    end subroutine

end module

module grcrbc

    use grcrbc_internal, only: &
      grcrbc_Initialize => Initialize, &
      grcrbc_SetInternalState => SetInternalState, &
!      grcrbc_PdistInterpolate => PdistInterpolate, &
      grcrbc_ComputeImplicitCoefficients => ComputeImplicitCoefficients, &
      grcrbc_UpdateInternalState => UpdateInternalState, &
      grcrbc_parameters_RCR => parameters_RCR, &
      grcrbc_parameters_Pdist => parameters_Pdist, &
!      grcrbc_PDist_current => PDist_current, &
      grcrbc_P_current => P_current, &
      grcrbc_Q_current => Q_current, &
      grcrbc_coeff_1_implicit => coeff_1_implicit, &
      grcrbc_coeff_2_implicit => coeff_2_implicit

end module

subroutine Dgrcrbc

    use grcrbc_internal
    use phcommonvars
    implicit none

    if (igrcrfile .gt. 0) then
        if (allocated(numtimepoints_Pdist)) deallocate (numtimepoints_Pdist)
        if (allocated(parameters_Pdist)) deallocate (parameters_Pdist)
        if (allocated(parameters_RCR)) deallocate (parameters_RCR)
        if (allocated(coeff_1_implicit)) deallocate (coeff_1_implicit)
        if (allocated(coeff_2_implicit)) deallocate (coeff_2_implicit)
        if (allocated(RCoverDt)) deallocate (RCoverDt)
!        if (allocated(PDist_current)) deallocate (PDist_current)
!        if (allocated(PDist_alpha)) deallocate (PDist_alpha)
        if (allocated(P_current)) deallocate (P_current)
        if (allocated(Q_current)) deallocate (Q_current)
        if (allocated(SrfArea)) deallocate (SrfArea)
    end if

end subroutine
