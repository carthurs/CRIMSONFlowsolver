      subroutine mpiset

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      logical     reorder

      call MPI_COMM_RANK (INEWCOMM, myrank)  ! we ditched the ierr that fortran 
	                                             ! normally have here to pacify Digital

      return
      end

