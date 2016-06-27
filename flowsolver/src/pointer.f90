       module pointer_data
!
!.... maximum number of blocks
!
         parameter ( MAXBLK2 = 500000 ) ! Note compiler was complaining
!                                       because MAXBLK in common.h be careful
!                                       to change both places
!
!.... data type definitions
!

         type r0d
          real*8, pointer :: p
         end type

         type r1d
           real*8, pointer :: p(:)
         end type
!
         type r2d
           real*8, pointer :: p(:,:)
         end type
!
         type r3d
           real*8, pointer :: p(:,:,:)
         end type
!         
         type r4d
           real*8, pointer :: p(:,:,:,:)
         end type
!
         type i1d
           integer, pointer :: p(:)
         end type
!
         type i2d
           integer, pointer :: p(:,:)
         end type
!
         type i3d
           integer, pointer :: p(:,:,:)
         end type
!
!.... pointer declarations
!
         type (i1d), dimension(MAXBLK2) ::  mmat,  mmatb
         ! added target attribute to this type for gfortran
         type (i2d), dimension(MAXBLK2), target ::  mien
         type (i2d), dimension(MAXBLK2) ::  mienb,  miBCB
         type (r2d), dimension(MAXBLK2) ::  mxmudmi
         type (r3d), dimension(MAXBLK2) ::  mBCB

         real*8, allocatable :: gmass(:)
       end module
