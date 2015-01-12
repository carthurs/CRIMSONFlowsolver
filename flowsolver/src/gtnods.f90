      subroutine gtnods

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      dimension irecvcount(numpe), numvec(numpe)
      write(*,*) "gtnods location A"
      if(numpe > 1) then
         irecvcount = 1
         numvec = nshg0
         write(*,*) "gtnods location B"
         call MPI_REDUCE_SCATTER (numvec, nshgt, irecvcount, &
                       MPI_INTEGER, MPI_SUM, INEWCOMM, ierr)
         write(*,*) "gtnods location DA", numvec
         write(*,*) "gtnods location DB", irecvcount
         write(*,*) "gtnods location DC",  ierr
      endif
     write(*,*) "gtnods location C"
      if (myrank .eq. master) then
         write(*,*) "gtnods location D", nshgt
         write(6,*) 'Total number of modes = ',nshgt
         write(*,*) "gtnods location E"
      endif
      
      return
      end
