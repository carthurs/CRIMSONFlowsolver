!
!--------------
!     drvAllreduce
!--------------
!     
      subroutine drvAllreduce ( eachproc, result, m )
!     
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!     
      dimension eachproc(m), result(m)
!     
      if (numpe > 1) then
         call MPI_ALLREDUCE ( eachproc, result, m,  &
              MPI_DOUBLE_PRECISION, MPI_SUM, INEWCOMM, ierr )
      else
         result = eachproc
      endif
!     
      return
      end
!     
!------------------
!     drvAllreducesclr
!------------------
!     
      subroutine drvAllreducesclr ( eachproc, result ) 
!     
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!     
      if (numpe > 1) then
         call MPI_ALLREDUCE ( eachproc, result, 1, &
              MPI_DOUBLE_PRECISION, MPI_SUM, INEWCOMM, ierr )
      else
         result = eachproc
      endif
!     
      return
      end

!------------------------------------------------------------------------
!
!   sum real*8 array over all processors
!
!------------------------------------------------------------------------
      subroutine sumgat (u, n, summed)

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      dimension u(nshg,n), ilwork(nlwork) 
      dimension sumvec(numpe), irecvcount(numpe)

      summed = sum(u)

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount,  &
              MPI_DOUBLE_PRECISION, MPI_SUM, INEWCOMM, ierr)

      endif

      return
      end

!------------------------------------------------------------------------
!
!   sum real*8 array of length nnp over all processors
!
!------------------------------------------------------------------------
      subroutine sumgatN (u, n, summed, nnp)

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      dimension u(nnp,n), ilwork(nlwork) 
      dimension sumvec(numpe), irecvcount(numpe)

! protect against underflow
!     summed = sum(u)
      summed = sum(u) + 1.e-20

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount,  &
              MPI_DOUBLE_PRECISION, MPI_SUM, INEWCOMM, ierr)

      endif

      return
      end

!------------------------------------------------------------------------
!
!   sum integer array over all processors
!
!------------------------------------------------------------------------
      subroutine sumgatInt (u, n, summed )

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      integer u(n), summed
      integer sumvec(numpe), irecvcount(numpe)

!$$$      ttim(62) = ttim(62) - tmr()

      summed = sum(u)

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount,  &
              MPI_INTEGER, MPI_SUM, INEWCOMM, ierr)

      endif
!$$$      ttim(62) = ttim(62) + tmr()

      return
      end

