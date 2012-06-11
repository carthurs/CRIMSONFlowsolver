
!*******************************************
!      ramg Interface
!      Interface for libLES
!******************************************* 
      subroutine ramg_interface( &
                 colm,rowp,lhsK,lhsP,flowDiag, &
                 mcgR,mcgZ, &
                 ilwork,BC,iBC,iper &
                 )
!$$$      use ramg_data
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      integer,intent(in),dimension(nshg+1)             :: colm
      integer,intent(in),dimension(nnz_tot)            :: rowp
      real(kind=8),intent(in),dimension(9,nnz_tot)     :: lhsK
      real(kind=8),intent(in),dimension(4,nnz_tot)     :: lhsP
      real(kind=8),intent(in),dimension(nshg)          :: mcgR
      real(kind=8),intent(inout),dimension(nshg)       :: mcgZ
      integer, intent(in), dimension(nlwork)           :: ilwork
      integer, intent(in),dimension(nshg)              :: iBC,iper
      real(kind=8),intent(in),dimension(nshg,ndofBC)   :: BC
 
          mcgZ = mcgR
          write (*,*) myrank,' preconditioning ppe'
      end subroutine ramg_interface
    
