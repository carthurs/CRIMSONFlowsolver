!    module globalArrays contains the arrays that formerly were passed
!    into (or initialized in) itrdrv but now must be reformed on every 
!    adaptivity step by a call
!      to process (note some will appear to be missing because they are
!      already part of other modules that will also be included in itrdrv


module globalArrays
    use, intrinsic :: iso_c_binding
    integer, target, allocatable :: nodetagfield(:)
    real (c_double), target, allocatable :: y(:,:)
    real (c_double), target, allocatable :: ac(:,:)
    real (c_double), target, allocatable :: u(:,:)
    real (c_double), target, allocatable :: yold(:,:)
    real (c_double), target, allocatable :: acold(:,:)
    real (c_double), target, allocatable :: uold(:,:)
    real*8, allocatable :: uref(:,:)
    real*8, allocatable :: ubar(:,:)

    real*8, allocatable :: xdist(:)
    real*8, allocatable :: xdnv(:,:)
    real*8, allocatable :: df_fem(:)

    real (c_double), allocatable :: temporary_array(:,:)

    integer (c_int), target, allocatable :: iBC(:)
    real*8, allocatable :: BC(:,:)
    real*8, allocatable :: solinc(:,:)
    real*8, allocatable :: res(:,:)
    integer, allocatable :: rowp(:,:)
    integer, allocatable :: colm(:)
    real*8, allocatable :: rerr(:,:)
    real*8, allocatable :: ybar(:,:)
    real*8, allocatable :: dummyVar(:)
    real*8, allocatable :: uhess(:,:)
    real*8, allocatable :: gradu(:,:)


    ! the ones listed below were formerly in readarrays but since they are
    ! unchanged in phasta preprocessing they have been moved here to make
    ! clear that everything in readarrays can be destroyed after reading
    ! while the arrays found here live as long as the mesh is unchanged
    ! A consequence is that now readnblk must include this data module so
    ! that  these arrays can be filled as before

    real (c_double), target, allocatable :: x(:,:)
    integer (c_int), target, allocatable :: ilwork(:)
    integer (c_int), target, allocatable :: iper(:)

    integer, allocatable :: ifath(:)
    integer, allocatable :: nsons(:)

    real*8, allocatable :: velbar(:,:)

    !
    ! new arrays that are used specifically for data assimilation
    !

    integer (c_int), target, allocatable :: inodesuniq(:) ! local unique index of nodes
    integer (c_int), target, allocatable :: ilinobsfunc_sol(:,:) ! simple linear observation functions
    integer (c_int), target, allocatable :: ilinobsfunc_acc(:,:)
    integer (c_int), target, allocatable :: ilinobsfunc_disp(:,:)
    integer (c_int), target, allocatable :: obsfunc_dist(:) ! distance observation flag
 
    !linear algebra matrices (actually allocated in itrdrv.f but destroyed
    !with global destructor)

    real*8, allocatable, dimension(:,:) :: aperm,  atemp, atempS
    real*8, allocatable, dimension(:,:,:) :: apermS

    real*8, allocatable, dimension(:,:) :: lhsP, lhsK, lhsS

end module


!-----------------------------------------------------------------------
!
!     Initialize:
!
!-----------------------------------------------------------------------
subroutine initGlobalArrays
    !
    use globalArrays
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    if (.not. allocated(nodetagfield)) allocate (nodetagfield(numnp))
    if (.not. allocated(y)) allocate (y(nshg,ndof))
    if (.not. allocated(ac)) allocate (ac(nshg,ndof))
    if (.not. allocated(yold)) allocate (yold(nshg,ndof))
    if (.not. allocated(acold)) allocate (acold(nshg,ndof))

    if (.not. allocated(temporary_array)) allocate(temporary_array(nshg,ndof))

    if (.not. allocated(ibc)) allocate (ibc(nshg))
    if (.not. allocated(BC)) allocate (BC(nshg,ndofBC))
    if (.not. allocated(solinc)) allocate (solinc(nshg,ndof))
    if (.not. allocated(rowp)) allocate (rowp(nshg,nnz))
    if (.not. allocated(colm)) allocate (colm(nshg+1))
    !allocate (x(nshg,nsd))
    if (.not. allocated(x)) allocate (x(numnp,nsd))
    if (.not. allocated(iper)) allocate (iper(nshg))
    if (.not. allocated(res)) allocate (res(nshg,nflow))
    !
    ! note that ilwork, ifath, nsons and velbar are all allocated in readnblk.f
    ! destructor will need to kill these.  Same is true for

    !
    ! conditional allocates
    !
    ! these are currently passed but should be moved to a
    ! module so that they can be conditionally allocated to save
    ! memory when not used
    if (.not. allocated(u)) allocate (u(nshg,nsd))
    if (.not. allocated(uold)) allocate (uold(nshg,nsd))
    if (.not. allocated(uref)) allocate (uref(nshg,nsd))
    if (.not. allocated(ubar)) allocate (ubar(nshg,nsd))
    if (.not. allocated(xdist)) allocate (xdist(nshg))
    if (.not. allocated(xdnv)) allocate (xdnv(nshg,nsd))
    if (.not. allocated(df_fem)) allocate (df_fem(nshg))
      

    if (.not. allocated(rerr)) allocate (rerr(nshg,numerr)) ! this one was pulled out because it
                                                            ! is passed down.  Soon, the error
                                                            ! fields should be put in a module and
                                                            ! "used" only where needed to save
                                                            ! this substantial memory not being
                                                            ! allocated when not needed
    if(ierrcalc.ne.0) then       ! these only live at itrdrv level and thus
                                  ! safe to conditional allocate
        if (.not. allocated(ybar)) allocate (ybar(nshg,5))
        if (.not. allocated(dummyVar)) allocate (dummyVar(nshg))
    endif
    if(ihessian.ne.0) then
        if (.not. allocated(uhess)) allocate (uhess(nshg,27))
        if (.not. allocated(gradu)) allocate (gradu(nshg,9))
    endif

end

!c-----------------------------------------------------------------------
!
!     GlobalDestruction Deallocate all allocatables in phasta that
!      depend on the mesh:
!
!-----------------------------------------------------------------------
subroutine GlobalDestruction
    !
    use globalArrays
    use pointer_data
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    !deallocate (y) ! this is a produced array that will be deallocated from the driver
    if (allocated(nodetagfield)) deallocate (nodetagfield)
    if (allocated(ac)) deallocate (ac)
    if (allocated(yold)) deallocate (yold)
    if (allocated(acold)) deallocate (acold)
    if (allocated(ibc)) deallocate (ibc)
    if (allocated(BC)) deallocate (BC)
    if (allocated(solinc)) deallocate (solinc)
    if (allocated(rowp)) deallocate (rowp)
    if (allocated(colm)) deallocate (colm)
    if (allocated(x)) deallocate (x)
    if (allocated(iper)) deallocate (iper)
    if (allocated(res)) deallocate (res)
    !
    ! conditional allocates
    !
    ! the next four are currently passed but should be moved to a
    ! module so that they can be conditionally allocated to save
    ! memory when not used

    if (allocated(u)) deallocate (u)
    if (allocated(uold)) deallocate (uold)
    if (allocated(xdist)) deallocate(xdist)
    if (allocated(xdnv)) deallocate(xdnv)
    if (allocated(df_fem)) deallocate(df_fem)


    if (numpe>1) then
        call Dctypes(ilwork) ! free the MPI derived datatypes
    endif

    if (allocated(ilwork)) deallocate (ilwork)
    if (allocated(ifath)) deallocate (ifath)
    if (allocated(nsons)) deallocate (nsons)
    if (allocated(velbar)) deallocate (velbar)
    !
    !
    call DShapeTable  ! destroy the shape functions
      
    call DdeformableWall
      
    call DmeasureWallDistance

    ! rerr one was pulled out because it
    ! is passed down.  Soon, the error
    ! fields should be put in a module and
    ! "used" only where needed to save
    ! this substantial memory not being
    ! allocated when not needed
    if (allocated(rerr)) deallocate (rerr)

    ! these only live at itrdrv level and thus safe to conditional allocate
    if(ierrcalc.ne.0) then
        if (allocated(ybar)) deallocate (ybar)
        if (allocated(dummyVar)) deallocate (dummyVar)
    endif
    if(ihessian.ne.0) then
        if (allocated(uhess)) deallocate (uhess)
        if (allocated(gradu)) deallocate (gradu)
    endif
      
    !
    !   Now destroy all of the other allocatables
    !
    !if(ilset.gt.0) call DspatVarVol
!    if (iLES/10 .eq. 2) call Daveall
!    if (iLES/10 .eq. 1) then
!        call DlhsGkeep
!        call Drlssave
!    endif
    call DspecialBC           ! conditionals inside
    call DincpBC              ! conditionals inside
    call DconvolImpFlow       ! conditionals inside
    call DconvolRCRFlow       ! conditionals inside
    call DconvolCORFlow       ! conditionals inside
    call DcalcFlowPressure    ! conditionals inside
    call DLagrangeMultipliers ! conditionals inside
    if(idtn.eq.1) call Ddtnmod
!    if(myrank.eq. master .and. irscale.ge.0) call DSPEBC
    call Dturbsa          ! always allocated???
    !      if((flmpr.ne.0) .or. (flmpl.ne.0))
    call Dgmass
    call DpvsQbi ! No conditional??

    call Dgrcrbc ! Nan rcr
    !
    !.... determine how many scalar equations we are going to need to solve
    !
    nsolflow=mod(impl(1),100)/10  ! 1 if solving flow
    nsolt=mod(impl(1),2)      ! 1 if solving temperature
    nsclrsol=nsolt+nsclr

    if (nsolflow.eq.1) then
        if (allocated(aperm)) deallocate (aperm)
        if (allocated(atemp)) deallocate (atemp)
        if (allocated(lhsP)) deallocate (lhsP)
        if (allocated(lhsK)) deallocate (lhsK)
    endif
    if(nsclrsol.gt.0) then
        if (allocated(apermS)) deallocate (apermS)
        if (allocated(atempS)) deallocate (atempS)
        if (allocated(lhsS)) deallocate (lhsS)
    endif

    !  I got tired of looking at SPEBC  and Andres stuff so it probably
    !               needs a more careful look when we try to get this stuff
    !               going but hopefully we can do some testing without that
    !               stuff first
    !
    !.... loop over the element-blocks
    !
    do iblk = 1, nelblk  ! I am not sure if this is the correct way to
                          ! deallocate memory at a pointer
        if (associated(mien(iblk)%p)) deallocate(mien(iblk)%p)
        if (associated(mmat(iblk)%p)) deallocate(mmat(iblk)%p)
        if (associated(mxmudmi(iblk)%p)) deallocate(mxmudmi(iblk)%p)
    enddo

    !
    !.... loop over the boundary-blocks
    !
    do iblk = 1, nelblb
        if (associated(mienb(iblk)%p)) deallocate(mienb(iblk)%p)
        if (associated(miBCB(iblk)%p)) deallocate(miBCB(iblk)%p)
        if (associated(mBCB(iblk)%p)) deallocate(mBCB(iblk)%p)
        if (associated(mmatb(iblk)%p)) deallocate(mmatb(iblk)%p)
    enddo
      

    return
      
end
