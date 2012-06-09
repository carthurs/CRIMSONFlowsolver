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
      subroutine getvel (y,ilwork, iBC,
     &                   nsons, ifath, velbar)


      include "common.h"
      include "mpif.h"
      include "auxmpi.h"

c
      dimension nsons(nfath),        rinvsons(nfath),
     &          velo(numnp,nflow),    velf(nfath,nflow),
     &          velft(nfath,nflow),   
     &          velbar(nfath,nflow),
     &          y(numnp,ndof), 
     &          ifath(numnp),
     &          ilwork(nlwork),        iBC(numnp)

c
c  for now keep the compressible numbering in velbar
c
      velo(:,1)=y(:,4)
      velo(:,2)=y(:,1)      
      velo(:,3)=y(:,2)
      velo(:,4)=y(:,3)
      if(nflow.eq.5) velo(:,5)=y(:,5)
      nsonmax=maxval(nsons)
      if ((nsonmax.eq.1)) then  ! we are doing local clipping -no homog dir
         velft=velo
      else
c     
c     zero on processor periodic nodes so that they will not be added twice
c     
         where(btest(iBC,10).or.btest(iBC,12))
            velo(:,1)=zero
            velo(:,2)=zero
            velo(:,3)=zero
            velo(:,4)=zero
         endwhere      
         if(nflow.eq.5) then
            where(btest(iBC,10).or.btest(iBC,12))
                velo(:,5)=zero
	    endwhere
	 endif         
         if (numpe.gt.1) then
            
            numtask = ilwork(1)
            itkbeg = 1
            
c     zero the nodes that are "solved" on the other processors  
            do itask = 1, numtask
               
               iacc   = ilwork (itkbeg + 2)
               numseg = ilwork (itkbeg + 4)
               
               if (iacc .eq. 0) then
                  do is = 1,numseg
                     isgbeg = ilwork (itkbeg + 3 + 2*is)
                     lenseg = ilwork (itkbeg + 4 + 2*is)
                     isgend = isgbeg + lenseg - 1
                     velo(isgbeg:isgend,:) = zero
                     velo(isgbeg:isgend,:) = zero
                  enddo
               endif
               
               itkbeg = itkbeg + 4 + 2*numseg
               
            enddo !itask
            
         endif  ! numpe.gt.1

         
         velf = zero
c     
c     accumulate sum of sons to the fathers
c     
         do i = 1,numnp
            ifathi=ifath(i)
	    velf(ifathi,1:nflow) = velf(ifathi,1:nflow) 
     &                           + velo(i,1:nflow)            
         enddo
         
c     
c     Now  the true fathers and serrogates combine results and update
c     each other.
c     
         if(numpe .gt. 1) then
            call drvAllreduce(velf, velft,nfath*nflow)
         else
            velft=velf
         endif
c     
c     xvelf is the sum of the sons for each father on this processor
c
c     xvelft is the sum of the sons for each father on all processor combined
c     (the same as if we had not partitioned the mesh for each processor)
c     
c     divide by # of sons to get average father for this step
c     
         rinvsons = one/nsons   ! division is expensive
         velft(:,1) = velft(:,1) * rinvsons(:) !  / nsons(:)
         velft(:,2) = velft(:,2) * rinvsons(:) !  / nsons(:)
         velft(:,3) = velft(:,3) * rinvsons(:) !  / nsons(:)
         velft(:,4) = velft(:,4) * rinvsons(:) !  / nsons(:)
	 if(nflow.eq.5) velft(:,5) = velft(:,5) * rinvsons(:)
      endif  ! end of homog direction averaging
      denom=max(one*(lstep),one)
      if(wtavei.lt.0) then
         tavef=one/denom
      else
         tavef=wtavei
      endif
               velbar(:,1)=velft(:,1)
         velbar(:,2)=velft(:,2)
         velbar(:,3)=velft(:,3)
         velbar(:,4)=velft(:,4)

      if(istep.eq.0) then
         velbar(:,:)=velft(:,:)
      else            
         velbar(:,:)=tavef*velft(:,:)+(one-tavef)*velbar(:,:)        
      endif
      
      return
      end



      


