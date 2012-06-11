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
      subroutine genadj (colm, rowp, icnt )
c     
      use pointer_data
c     
      include "common.h"
c     
      integer rowp(nshg*nnz),         colm(nshg+1)
      integer adjcnt(nshg),    row_fill_list(nshg,6*nnz), mloc(1)
c                                          change ^ if overflow
c                                   also change overflow check in asadj TWICE
      integer tmprdim(1)
      real*8, allocatable, dimension(:) :: tmpr

      adjcnt=0
      
      do iblk = 1, nelblk
c     
c.... set up the parameters
c     
         iel    = lcblk(1,iblk)
         lelCat = lcblk(2,iblk)
         lcsyst = lcblk(3,iblk)
         iorder = lcblk(4,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         npro   = lcblk(1,iblk+1) - iel 
         
c     
c.... compute sparse matrix data structures
c     
         call Asadj (row_fill_list,                       
     &               mien(iblk)%p,  adjcnt )
         
      enddo
      
      call sumgatInt ( adjcnt, nshg, nnonzero)
      if ( myrank .eq. master) then
         write (*,*) 'Number of global nonzeros ',nnonzero
      endif

c     
c     build the colm array
c     
      colm(1)=1
      do i=1,nshg
         colm(i+1)=colm(i)+adjcnt(i)
      enddo
c     
c     sort the rowp into increasing order
c     
      ibig=10*nshg
      icnt=0
      tmprdim=maxval(adjcnt)
      allocate (tmpr(tmprdim(1)))
      do i=1,nshg
         ncol=adjcnt(i)
         tmpr(1:ncol)=row_fill_list(i,1:ncol)
         do j=1,ncol
            icnt=icnt+1
            imin=minval(tmpr(1:ncol))
            mloc=minloc(tmpr(1:ncol))
            rowp(icnt)=imin
            tmpr(mloc(1))=ibig
         enddo
      enddo

      maxfill=tmprdim(1)
      write(*,*) 'maxfill=',maxfill
      nnza=icnt/nshg +1
      if(icnt.gt.nnz*nshg) then
         write(*,*) 'increase nnz in genmat to',nnza
         stop
      else
         write(*,*) 'nnz ok  nnz=',nnz,' actually needed',nnza   
         write(*,*) myrank,' is my rank and my nnz_tot is: ',icnt   
      endif
      return
      end









