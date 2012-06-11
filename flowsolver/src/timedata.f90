      module timedata

      integer ntspts, freq, iterat, varcod
      real*8  tolpt
      logical exts

      real*8,  allocatable :: ptts(:,:)
      real*8,  allocatable :: varts(:)

      end module


!-----------------------------------------------------------------------
! allocate the arrays
!-----------------------------------------------------------------------


      subroutine sTD 

      use timedata

      allocate (ptts(ntspts,3))
      allocate (varts(ntspts))

      return
      end
!-----------------------------------------------------------------------
! delete the arrays
!-----------------------------------------------------------------------

      
      subroutine dTD 

      use timedata

      deallocate (ptts)
      deallocate (varts)

      return
      end
