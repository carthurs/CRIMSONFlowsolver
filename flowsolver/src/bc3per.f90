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
        subroutine bc3per (iBC,  res, iper, ilwork,nQs)
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of the periodic nodes after Ap product
c
c input:
c  iBC   (nshg)        : Boundary Condition Code
c  iper  (nshg)        : partners of periodic nodes
c  res   (nshg,nQs)   : residual before BC is applied
c
c output:
c  res   (nshg,nQs)   : residual after satisfaction of BC
c
c----------------------------------------------------------------------
c
        use periodicity  ! this gives you rcount(1:nshg) (real*8)
        include "common.h"
c
        dimension iBC(nshg),
     &            res(nshg,nQs),           ilwork(nlwork),
     &            iper(nshg)
c
c.... local periodic (and axisymmetric) boundary conditions (no communications)
c
           do j = 1,nshg
              if ((btest(iBC(j),10)) .or. (btest(iBC(j),12))) then
                 i = iper(j)
                 res(i,:) = res(i,:) + res(j,:)
                 res(j,:) = zero
              endif
           enddo


        if(numpe.gt.1) then
c
c.... nodes treated on another processor are eliminated
c     
           numtask = ilwork(1)
           itkbeg = 1
           
           do itask = 1, numtask
              
              iacc   = ilwork (itkbeg + 2)
              numseg = ilwork (itkbeg + 4)
              
              if (iacc .eq. 0) then
                 do is = 1,numseg
                    isgbeg = ilwork (itkbeg + 3 + 2*is)
                    lenseg = ilwork (itkbeg + 4 + 2*is)
                    isgend = isgbeg + lenseg - 1
                    res(isgbeg:isgend,:) = zero
                 enddo
              endif
              
              itkbeg = itkbeg + 4 + 2*numseg
              
           enddo
        endif
c
c.... return
c
        return
        end





