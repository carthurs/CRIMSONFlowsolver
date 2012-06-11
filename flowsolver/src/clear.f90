        subroutine clear (clr, n)
!
!----------------------------------------------------------------------
!
!  This routine clears a floating point array.
!
! input:
!  n            : number of floating points to be zeroed
!
! output:
!  clr (n)      : the array to be zeroed
!
! Farzin Shakib, Summer 1985.
!----------------------------------------------------------------------
!
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension clr(n)
!
        do i = 1, n
          clr(i) = zero
        enddo
!
        return
        end
