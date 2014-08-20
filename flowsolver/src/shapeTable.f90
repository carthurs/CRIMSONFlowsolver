!    module globalArrays contains the arrays that formerly were passed
!    into (or initialized in) itrdrv but now must be reformed on every 
!    adaptivity step by a call
!      to process (note some will appear to be missing because they are
!      already part of other modules that will also be included in itrdrv


      module shapeTable
      use iso_c_binding
      
      real (c_double), target, allocatable :: shp(:,:,:)
      real (c_double), target, allocatable :: shgl(:,:,:,:)
      real (c_double), target, allocatable :: shpb(:,:,:)
      real (c_double), target, allocatable :: shglb(:,:,:,:)
      
      end module

!-----------------------------------------------------------------------
!
!     Initialize:
!
!-----------------------------------------------------------------------
      subroutine initShapeTable
!
      use shapeTable
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      if (.not. allocated(shp)) allocate (shp(MAXTOP,maxsh,MAXQPT))
      if (.not. allocated(shgl)) allocate (shgl(MAXTOP,nsd,maxsh,MAXQPT))
      if (.not. allocated(shpb)) allocate (shpb(MAXTOP,maxsh,MAXQPT))
      if (.not. allocated(shglb)) allocate (shglb(MAXTOP,nsd,maxsh,MAXQPT))

      end

