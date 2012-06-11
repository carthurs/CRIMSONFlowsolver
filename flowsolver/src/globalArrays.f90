!    module globalArrays contains the arrays that formerly were passed
!    into (or initialized in) itrdrv but now must be reformed on every 
!    adaptivity step by a call
!      to process (note some will appear to be missing because they are
!      already part of other modules that will also be included in itrdrv


      module globalArrays
      use, intrinsic :: iso_c_binding
      real (c_double), target, allocatable :: y(:,:)
      real*8, allocatable :: ac(:,:)
      real (c_double), target, allocatable :: u(:,:)
      real*8, allocatable :: yold(:,:)
      real*8, allocatable :: acold(:,:)
      real*8, allocatable :: uold(:,:)

      integer (c_int), target, allocatable :: iBC(:)
      real*8, allocatable :: BC(:,:)
      real*8, allocatable :: solinc(:,:)
      real*8, allocatable :: res(:,:)
      real*8, allocatable :: rowp(:,:)
      real*8, allocatable :: colm(:)
      real*8, allocatable :: rerr(:,:)
      real*8, allocatable :: ybar(:,:)
      real*8, allocatable :: dummyVar(:)
      real*8, allocatable :: uhess(:,:)
      real*8, allocatable :: gradu(:,:)
      
      real*8, allocatable :: xdist(:)
      real*8, allocatable :: xdnv(:,:)


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

      integer, allocatable :: inodesuniq(:) ! local unique index of nodes
      integer, allocatable :: ilinobsfunc_sol(:,:) ! simple linear observation functions
      integer, allocatable :: ilinobsfunc_acc(:,:)
      integer, allocatable :: ilinobsfunc_disp(:,:)
 
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
      if (.not. allocated(y)) allocate (y(nshg,ndof))
      if (.not. allocated(ac)) allocate (ac(nshg,ndof))
      if (.not. allocated(yold)) allocate (yold(nshg,ndof))
      if (.not. allocated(acold)) allocate (acold(nshg,ndof))
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
      ! the next two are currently passed but should be moved to a
      ! module so that they can be conditionally allocated to save
      ! memory when not used
      if (.not. allocated(u)) allocate (u(nshg,nsd))
      if (.not. allocated(uold)) allocate (uold(nshg,nsd))
      if (.not. allocated(xdist)) allocate (xdist(numnp))
      if (.not. allocated(xdnv)) allocate (xdnv(numnp,nsd))
      

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

!
!

!!     Temporary subroutine for proces_adapt
!      subroutine initGlobalArrays_adapt
!!
!      use globalArrays
!  use phcommonvars  
!  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!      !allocate (y(nshg,ndof))
!      !allocate (ac(nshg,ndof))
!      allocate (yold(nshg,ndof))
!      allocate (acold(nshg,ndof))
!      !allocate (ibc(nshg))
!      !allocate (BC(nshg,ndofBC))
!      allocate (solinc(nshg,ndof))
!      allocate (rowp(nshg,nnz))
!      allocate (colm(nshg+1))
!      !allocate (x(nshg,nsd))
!      !allocate (iper(nshg))
!      allocate (res(nshg,nflow))
!!
!! note that ilwork, ifath, nsons and velbar are all allocated in readnblk.f
!! destructor will need to kill these.  Same is true for
!
!!
!! conditional allocates
!!
!      ! the next two are currently passed but should be moved to a
!      ! module so that they can be conditionally allocated to save
!      ! memory when not used
!      allocate (u(nshg,nsd))
!      allocate (uold(nshg,nsd))
!
!      allocate (rerr(nshg,numerr)) ! this one was pulled out because it
!                                   ! is passed down.  Soon, the error
!                                   ! fields should be put in a module and
!                                   ! "used" only where needed to save
!                                   ! this substantial memory not being
!                                   ! allocated when not needed
!      if(ierrcalc.ne.0) then       ! these only live at itrdrv level and thus
!                                   ! safe to conditional allocate
!         allocate (ybar(nshg,5))
!         allocate (dummyVar(nshg))
!      endif
!      if(ihessian.ne.0) then
!         allocate (uhess(nshg,27))
!         allocate (gradu(nshg,9))
!      endif
!
!      end


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
!
! conditional allocates
!
      ! the next two are currently passed but should be moved to a
      ! module so that they can be conditionally allocated to save
      ! memory when not used
      
      if (allocated(u)) deallocate (u)
      if (allocated(uold)) deallocate (uold)
      if (allocated(xdist)) deallocate (xdist)
      if (allocated(xdnv)) deallocate (xdnv)
      
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
      if (iLES/10 .eq. 2) call Daveall
      if (iLES/10 .eq. 1) then
         call DlhsGkeep 
         call Drlssave 
      endif
      call DspecialBC           ! conditionals inside
      call DincpBC              ! conditionals inside
      call DconvolImpFlow       ! conditionals inside
      call DconvolRCRFlow       ! conditionals inside
      call DconvolCORFlow       ! conditionals inside
      call DcalcFlowPressure    ! conditionals inside
      call DLagrangeMultipliers ! conditionals inside
      if(idtn.eq.1) call Ddtnmod
      if(myrank.eq. master .and. irscale.ge.0) call DSPEBC
      call Dturbsa          ! always allocated???
!      if((flmpr.ne.0) .or. (flmpl.ne.0)) 
      call Dgmass
      call DpvsQbi ! No conditional??
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

!
!
!

!      subroutine GlobalDestruction_adapt
!!
!      use globalArrays
!      use pointer_data
!  use phcommonvars  
!  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!      !deallocate (y)
!      !deallocate (ac)
!      deallocate (yold)
!      deallocate (acold)
!      !deallocate (ibc)
!      !deallocate (BC)
!      deallocate (solinc)
!      deallocate (rowp)
!      deallocate (colm)
!      !deallocate (x)
!      !deallocate (iper)
!      deallocate (res)
!
!      !deallocate (ilwork)
!      deallocate (ifath)
!      deallocate (nsons)
!      deallocate (velbar)
!!
!!
!      call DShapeTable  ! destroy the shape functions
!!
!! conditional allocates
!!
!      ! the next two are currently passed but should be moved to a
!      ! module so that they can be conditionally allocated to save
!      ! memory when not used
!      deallocate (u)
!      deallocate (uold)
!
!      deallocate (rerr)            ! this one was pulled out because it
!                                   ! is passed down.  Soon, the error
!                                   ! fields should be put in a module and
!                                   ! "used" only where needed to save
!                                   ! this substantial memory not being
!                                   ! allocated when not needed
!      if(ierrcalc.ne.0) then       ! these only live at itrdrv level and thus
!                                   ! safe to conditional allocate
!         deallocate (ybar)
!         deallocate (dummyVar)
!      endif
!      if(ihessian.ne.0) then
!         deallocate (uhess)
!         deallocate (gradu)
!      endif
!!
!!   Now destroy all of the other allocatables
!!
!      if(ilset.gt.0) call DspatVarVol
!      if (iLES/10 .eq. 2) call Daveall
!      if (iLES/10 .eq. 1) then
!         call DlhsGkeep
!         call Drlssave
!      endif
!      call DspecialBC           ! conditionals inside
!      call DconvolImpFlow       ! conditionals inside
!      if(idtn.eq.1) call Ddtnmod
!      if(myrank.eq. master .and. irscale.ge.0) call DSPEBC
!      call Dturbsa          ! always allocated???
!!      if((flmpr.ne.0) .or. (flmpl.ne.0))
!      call Dgmass
!      call DpvsQbi ! No conditional??
!!
!!.... determine how many scalar equations we are going to need to solve
!!
!      nsolflow=mod(impl(1),100)/10  ! 1 if solving flow
!      nsolt=mod(impl(1),2)      ! 1 if solving temperature
!      nsclrsol=nsolt+nsclr
!
!      if (nsolflow.eq.1) then
!         deallocate (aperm)
!         deallocate (atemp)
!         deallocate (lhsP)
!         deallocate (lhsK)
!      endif
!      if(nsclrsol.gt.0) then
!         deallocate (apermS)
!         deallocate (atempS)
!         deallocate (lhsS)
!      endif
!
!!  I got tired of looking at SPEBC  and Andres stuff so it probably
!!               needs a more careful look when we try to get this stuff
!!               going but hopefully we can do some testing without that
!!               stuff first
!!
!!.... loop over the element-blocks
!!
!      do iblk = 1, nelblk  ! I am not sure if this is the correct way to
!                           ! deallocate memory at a pointer
!         deallocate(mien(iblk)%p)
!         deallocate(mmat(iblk)%p)
!         deallocate(mxmudmi(iblk)%p)
!      enddo
!
!!
!!.... loop over the boundary-blocks
!!
!      do iblk = 1, nelblb
!         deallocate(mienb(iblk)%p)
!         deallocate(miBCB(iblk)%p)
!         deallocate(mBCB(iblk)%p)
!         deallocate(mmatb(iblk)%p)
!      enddo
!
!
!      return
!
!      end

