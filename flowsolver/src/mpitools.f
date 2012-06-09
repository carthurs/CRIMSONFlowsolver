c
c  Copyright (c) 2000-2007, Stanford University, 
c     Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
c     Charles A. Taylor (see SimVascular Acknowledgements file 
c     for additional contributors to the source code).
c
c  All rights reserved.
c
c  Redistribution and use in source and binary forms, with or without 
c  modification, are permitted provided that the following conditions 
c  are met:
c
c  Redistributions of source code must retain the above copyright notice,
c  this list of conditions and the following disclaimer. 
c  Redistributions in binary form must reproduce the above copyright 
c  notice, this list of conditions and the following disclaimer in the 
c  documentation and/or other materials provided with the distribution. 
c  Neither the name of the Stanford University or Rensselaer Polytechnic
c  Institute nor the names of its contributors may be used to endorse or
c  promote products derived from this software without specific prior 
c  written permission.
c
c  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
c  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
c  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
c  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
c  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
c  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
c  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
c  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
c  DAMAGE.
c
c
c     
c--------------
c     drvAllreduce
c--------------
c     
      subroutine drvAllreduce ( eachproc, result, m )
c     
      include "common.h"
      include "mpif.h"
c     
      dimension eachproc(m), result(m)
c     
      if (numpe > 1) then
         call MPI_ALLREDUCE ( eachproc, result, m, 
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      else
         result = eachproc
      endif
c     
      return
      end
c     
c------------------
c     drvAllreducesclr
c------------------
c     
      subroutine drvAllreducesclr ( eachproc, result ) 
c     
      include "common.h"
      include "mpif.h"
c     
      if (numpe > 1) then
         call MPI_ALLREDUCE ( eachproc, result, 1,
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
      else
         result = eachproc
      endif
c     
      return
      end

c------------------------------------------------------------------------
c
c   sum real*8 array over all processors
c
c------------------------------------------------------------------------
      subroutine sumgat (u, n, summed)

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension u(nshg,n), ilwork(nlwork) 
      dimension sumvec(numpe), irecvcount(numpe)

      summed = sum(u)

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount, 
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      endif

      return
      end

c------------------------------------------------------------------------
c
c   sum real*8 array of length nnp over all processors
c
c------------------------------------------------------------------------
      subroutine sumgatN (u, n, summed, nnp)

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      dimension u(nnp,n), ilwork(nlwork) 
      dimension sumvec(numpe), irecvcount(numpe)

c protect against underflow
c     summed = sum(u)
      summed = sum(u) + 1.e-20

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount, 
     &        MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      endif

      return
      end

c------------------------------------------------------------------------
c
c   sum integer array over all processors
c
c------------------------------------------------------------------------
      subroutine sumgatInt (u, n, summed )

      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

      integer u(n), summed
      integer sumvec(numpe), irecvcount(numpe)

c$$$      ttim(62) = ttim(62) - tmr()

      summed = sum(u)

      if (numpe > 1) then
         irecvcount = 1
         sumvec = summed
         call MPI_REDUCE_SCATTER (sumvec, summed, irecvcount, 
     &        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      endif
c$$$      ttim(62) = ttim(62) + tmr()

      return
      end

