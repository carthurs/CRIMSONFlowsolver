        subroutine iclear (iclr, n)
!
!----------------------------------------------------------------------
!
! This routine clears an integer array.
!
! input:
!  n            : number of integers to be zeroed
!
! output:
!  iclr (n)     : the array to be zeroed
!
!
! Farzin Shakib, Summer 1985.
!----------------------------------------------------------------------
!
        dimension iclr(n)
!
        do i = 1, n
          iclr(i) = 0
        enddo
!
        return
        end
