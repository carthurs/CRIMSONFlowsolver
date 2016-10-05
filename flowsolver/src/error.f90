        subroutine error (routin, variab, num)
!
!----------------------------------------------------------------------
!
! This utility routine prints out the error and stops the program.
!
! input:
!  routin       : name of the routine where the error occurred
!  variab       : an 8-character error message
!  num          : any integer number associated with the error
!
! Farzin Shakib, Summer 1985.
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
!
        character(len=*) routin, variab
!
        data ierchk /0/
!
!.... check for redundant error
!
        if (ierchk .eq. 1) stop
        ierchk = 1
!
!.... open file
!
        open (unit=ierror, file=ferror, status='unknown')
!
!.... print the error
!
        write (*,1000) title, routin, variab, num
        if (num .ne. 0) write (ierror,1000) title, routin, variab, num
        if (num .eq. 0) write (ierror,1000) title, routin, variab
!
!.... halt the process
!
        close (ierror)


        WRITE(6,'(A,G14.6)') 'Life: ',death - birth
        if (numpe > 1) then
           call MPI_ABORT(INEWCOMM)
        endif
        
 
1000    format(' ',a80,//, &
               ' ****** Error occurred in routine <',a8,'>',/, &
                '  Error code :',a8,:,' : ',i8,//)
        end
