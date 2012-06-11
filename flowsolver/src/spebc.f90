      module spebc
      
      integer  nfin,  nelint, npin, npint, nfint
      real*8   sang
      real*8   xnrml,  ynrml,  znrml,  aR, aI
      
      

      real*8, allocatable :: xyn(:), xynin(:)
      real*8, allocatable :: xcyl(:,:),xintl(:,:,:),xsinfin(:,:,:)


      integer, allocatable :: ien2D(:,:), nen1(:), elcnfin(:,:)
      integer, allocatable :: nrint(:), imax(:) 
  
      end module
      
!-----------------------------------------------------------------------
! deallocate the spebc arrays
!-----------------------------------------------------------------------
      subroutine DSPEBC()
      
      use spebc
      if (allocated(xyn)) deallocate (xyn)
      if (allocated(xynin)) deallocate (xynin)
      if (allocated(xcyl)) deallocate (xcyl)
      if (allocated(xintl)) deallocate (xintl)
      if (allocated(xsinfin)) deallocate (xsinfin)

      if (allocated(ien2D)) deallocate (ien2D)
      if (allocated(nen1)) deallocate (nen1)
      if (allocated(elcnfin)) deallocate (elcnfin)
      if (allocated(nrint)) deallocate (nrint)
      if (allocated(imax)) deallocate (imax)
      

      
      return
      end      

!-----------------------------------------------------------------------
! allocate the spebc arrays
!-----------------------------------------------------------------------
      subroutine setSPEBC(numnp,nsd)
      
      use spebc

      if (.not. allocated(xyn)) allocate (xyn(numnp))
      if (.not. allocated(xynin)) allocate (xynin(numnp))
      if (.not. allocated(xcyl)) allocate (xcyl(numnp,nsd))
      if (.not. allocated(nen1)) allocate (nen1(numnp))

!      allocate (elcnpin(numnp))
!      allocate (xsi(numnp,nsd))
      

      
      return
      end

