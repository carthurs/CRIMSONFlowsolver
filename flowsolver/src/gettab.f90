        subroutine gettab  (mut,  rhot,  xst)
!
!-----------------------------------------------------------------------
!
!  This subroutine reads the three tables for equilibrium chemistry.
!
!
! output:
!
!    mut  (71,451)      : specific chemical potential function of (p,T)
!    rhot (71,451)      : density function of (p,T)
!    xst  (5,71,451)    : mole fractions functions of (p,T)
!
! Note: These three arrays are always in double precision.
! 
! Frederic Chalot and Zdenek Johan, Fall 1990.
!-----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        real*8  mut(71,451),  rhot(71,451),  xst(5,71,451)
!
!.... open table file
!
        open (unit=itable, file=ftable, form='unformatted', &
                                        status='unknown')
!
!.... read tables
!
        read (itable) mut
!
        read (itable) rhot
!
        read (itable) xst
!
!.... close table file
!
        close(unit=itable)
!
!.... end
!
        return
        end
