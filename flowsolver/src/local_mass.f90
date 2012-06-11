      module local_mass
      
      parameter (MAXBLK2 = 500000)
      
      integer :: iblock = 0
      integer :: have_local_mass = 0
      
      type r2d
      real*8, pointer :: p(:,:,:)
      end type
      
      type (r2d), dimension(MAXBLK2) :: lmassinv
      
      end module
 
