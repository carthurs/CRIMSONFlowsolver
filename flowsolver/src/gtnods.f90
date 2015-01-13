      subroutine gtnods

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      dimension irecvcount(numpe), numvec(numpe)
      if(numpe > 1) then
         irecvcount = 1
         numvec = nshg0
         call MPI_REDUCE_SCATTER (numvec, nshgt, irecvcount, &
                       MPI_INTEGER, MPI_SUM, INEWCOMM, ierr)
      endif
      if (myrank .eq. master) then
         write(6,*) 'Total number of modes = ',nshgt
      endif
      
      return
      end
