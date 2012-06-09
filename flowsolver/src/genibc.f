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
        subroutine geniBC (iBC)
c
c----------------------------------------------------------------------
c This routine reads the boundary condition codes.
c
c output: 
c  iBC   (nshg)        : Boundary Condition code
c
c         = 1 * iBC_1 + 2 * iBC_2 + 4 * iBC_3
c              density   temperature   pressure
c
c    if nsd = 3:
c
c        +  8 * iBC_4 +  16 * iBC_5 +  32 * iBC_6
c           x1-velocity   x2-velocity   x3-velocity
c
c        + 64 * iBC_7 + 128 * iBC_8 + 256 * iBC_9 + 512 * iBC_10
c          sclr1         sclr2        sclr3         sclr4
c
c        + 1024 * iBC_11  + 2048* iBC_12 
c          perioidicity     spebc          
c
c  nBC   (nshg)        : Boundary Condition mapping array
c
c----------------------------------------------------------------------
c
c
        use readarrays          ! used to access iBCtmp
        include "common.h"
c
c Arrays in the following 1 line are now dimensioned in readnblk
c        dimension iBCtmp(numpbc)
c
        dimension iBC(nshg)
        dimension itemp(6)
c
c.... set the iBC array
c
        iBC = 0
c
        if(numpbc.eq.0) return  ! sometimes there are no BC's on a partition
        where (nBC(:) .ne. 0) iBC(:) = iBCtmp(nBC(:))
c
c.... echo the input iBC array only if other than zero
c
        if (necho .lt. 3) then
          nn = 0
          do n = 1, nshg
            if (nBC(n) .ne. 0) then
              nb = nBC(n)
              nn = nn + 1
              if (mod(nn,50).eq.1) write(iecho,1000)ititle,(j,j=1,ndof)
              itemp(   1) = mod(iBCtmp(nb)   ,2) - mod(iBCtmp(nb)/ 4,2)
              itemp(   2) = mod(iBCtmp(nb)/ 8,2)
              itemp(   3) = mod(iBCtmp(nb)/16,2)
              itemp(   4) = mod(iBCtmp(nb)/32,2)
              itemp(ndof) = mod(iBCtmp(nb)/ 2,2)
              write(iecho,1100) n,(itemp(i),i=1,ndof)
            endif
          enddo
        endif
        deallocate(iBCtmp)
c
c.... return
c
        return
c
c.... end of file error handling
c
999     call error ('geniBC  ','end file',ibndc)
c
1000    format(a80,//,
     &  ' N o d a l   B o u n d a r y   C o n d i t i o n   C o d e',//,
     &  '    Node   ',13x,6('dof',i1,:,6x))
1100    format(2x,i5,10x,5i10)
c
        end
