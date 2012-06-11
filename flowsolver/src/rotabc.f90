      subroutine rotabc (global, iBC, code)
!---------------------------------------------------------------------
! 
! This subroutine is responsible for rotating 
! the residual and solution vectors for axisymmetric BC's.
!
! input:   
!     global(nshg,n): global vector to be rotated.
!     code:            = 'in' for rotating with the residual
!                      = 'out' for rotating the solution 
!
!  note that the cos and sin of the rotation angles are preprocessed and
!  stored in acs(1 and 2) respectively.
!
!---------------------------------------------------------------------
!
      use specialBC  ! gives us acs, contains (:,1)=cos(theta) (:,2)=sin(theta)
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
 
      dimension global(nshg,2),             iBC(nshg), &
                tmp(nshg)
 
      character*3 code

      if (code .eq. 'in ') then
         where( btest(iBC,10))
            tmp         =  global(:,1)*acs(:,1) - global(:,2)*acs(:,2)
            global(:,2) =  global(:,1)*acs(:,2) + global(:,2)*acs(:,1)
            global(:,1) = tmp
         endwhere
      else  if (code .eq. 'out') then
         where( btest(iBC,10))
            tmp         =  global(:,1)*acs(:,1) + global(:,2)*acs(:,2)
            global(:,2) = -global(:,1)*acs(:,2) + global(:,2)*acs(:,1)
            global(:,1) = tmp
         endwhere
      else 
         call error ('rotabc  ','code    ',0)
      endif

      return
      end
