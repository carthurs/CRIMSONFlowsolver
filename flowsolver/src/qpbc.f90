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
      subroutine qpbc( rmass, qres, iBC, iper, ilwork )
c---------------------------------------------------------------------
c
c This routine satisfies the periodic boundary conditions
c on the diffusive flux residual and mass matrix
c
c input:
c     rmass   (nshg)              : mass matrix
c     qres    (nshg,(nflow-1)*nsd) : diffusive flux vector
c 
c output: modified qres and rmass 
c---------------------------------------------------------------------
      include "common.h"
      
      dimension rmass(nshg), qres(nshg,idflx),
     &          iBC(nshg), iper(nshg),uv(nshg,2),
     &          tmpvec(nshg,4), tmp(nshg)  
c
      if(iabc==1) then   !are there any axisym bc's
      do i=1,idflx/nsd
         do j=1,2
            istrt=j+(i-1)*(nflow-1)
            uv(:,j)=qres(:,istrt)
         enddo
         call rotabc(uv, iBC, 'in ')
         do j=1,2
            istrt=j+(i-1)*(nflow-1)
            qres(:,istrt)=uv(:,j)
         enddo
      enddo
      endif
c
c
c.... compute qi for node A, i.e., qres <-- qres/rmass
c
       if (numpe > 1) then
          call commu (qres  , ilwork, idflx  , 'in ')
          call commu (rmass , ilwork,  1            , 'in ')
       endif
c
c  take care of periodic boundary conditions
c  but not on surface tension terms in qres(:,10-12)
c  that are used to compute normal vector
c
        idflow = (nflow-1)*nsd
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
c            qres(i,:) = qres(i,:) + qres(j,:)
            qres(i,1:idflow) = qres(i,1:idflow) + qres(j,1:idflow)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(j) = rmass(i)
c            qres(j,:) = qres(i,:)
            qres(j,1:idflow) = qres(i,1:idflow)
          endif
        enddo
c
c.... invert the diagonal mass matrix and find q
c
        rmass = one/rmass
       
       do i=1, idflx
          qres(:,i) = rmass*qres(:,i)
       enddo
       if (isurf .eq. 1) then
          idflow=(nflow-1)*nsd
c
c.... calculation of the unit normal vector
c
           tmp =  sqrt(qres(:,idflow+1)**2
     &               + qres(:,idflow+2)**2
     &               + qres(:,idflow+3)**2)
           do i = 1, nshg
             if (tmp(i) .lt. 0.0001) tmp(i) = 0.0001
           end do
           tmp = one/tmp

          do i=1, nsd
             qres(:,idflow+i) = tmp*qres(:,idflow+i)
          enddo 
       endif

       if(numpe > 1) then
          call commu (qres, ilwork, idflx, 'out')    
       endif

       if(iabc==1) then         !are there any axisym bc's
c
c       slave has masters value, for abc we need to rotate it
c
          do i=1,idflx/nsd
             do j=1,2
                istrt=j+(i-1)*(nflow-1)
                uv(:,j)=qres(:,istrt)
             enddo
             call rotabc(uv, iBC, 'out')
             do j=1,2
                istrt=j+(i-1)*(nflow-1)
                qres(:,istrt)=uv(:,j)
             enddo
          enddo
       endif
       
c
c.... return
c    
        return
        end


      subroutine qpbcSclr( rmass, qres, iBC, iper, ilwork )
c---------------------------------------------------------------------
c
c This routine satisfies the periodic boundary conditions
c on the diffusive flux residual and mass matrix
c
c input:
c     rmass   (nshg)              : mass matrix
c     qres    (nshg, nsd)         : diffusive flux vector
c 
c output: modified qres and rmass 
c---------------------------------------------------------------------
      include "common.h"
      
      dimension rmass(nshg), qres(nshg,nsd),
     &          iBC(nshg), iper(nshg)

      if(iabc==1) !are there any axisym bc's
     &         call rotabc(qres, iBC,  'in ')

c
c.... compute qi for node A, i.e., qres <-- qres/rmass
c
       if (numpe > 1) then
          call commu (qres  , ilwork, nsd  , 'in ')
          call commu (rmass , ilwork,  1   , 'in ')
       endif

c
c  take care of periodic boundary conditions
c
        do j= 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            qres(i,:) = qres(i,:) + qres(j,:)
          endif
        enddo

        do j= 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(j) = rmass(i)
            qres(j,:) = qres(i,:)
          endif
        enddo
c
c.... invert the diagonal mass matrix and find q
c
        rmass = one/rmass
       
       do i=1, nsd
          qres(:,i) = rmass*qres(:,i)
       enddo

       if(numpe > 1) then
          call commu (qres, ilwork, nsd, 'out')    
       endif

      if(iabc==1) !are there any axisym bc's
     &         call rotabc(qres, iBC, 'out')

c
c.... return
c    
        return
        end

