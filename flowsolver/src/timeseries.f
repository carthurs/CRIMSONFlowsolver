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

      subroutine timeseries(ycl, xl, ien, sgn)

      use timedata
      include "common.h"

      dimension shape(npro,nshl), ycl(npro,nshl,ndofl),
     &     ien(npro,nshl), xl(npro,nenl,nsd),         
     &     sgn(npro,nshl)
      real*8 al(npro,nenl,nsd), 
     &     zi0(npro,nsd), detaij(npro), dzi0(npro,nsd),
     &     m11(npro), m12(npro), m13(npro), m21(npro), m22(npro),
     &     m23(npro), m31(npro), m32(npro), m33(npro), 
     &     r1(npro), r2(npro), r3(npro), shgradl(npro,nshl,nsd)
     
      real*8 xts1, xts2, xts3, soln
      integer e

      do jj = 1, ntspts
            
         xts1 = ptts(jj,1)
         xts2 = ptts(jj,2)
         xts3 = ptts(jj,3)

         if(lcsyst.eq.2) then   ! hex
               
            call get_a_not_hex(xl,al) ! get mapping poly. coeff.
         
c...  get initial guess for Newton Iteration procedure
         
            detaij(:) = -al(:,2,1)*al(:,3,2)*al(:,4,3) + 
     &           al(:,2,1)*al(:,4,2)*al(:,3,3) + al(:,2,2)*
     &           al(:,3,1)*al(:,4,3) - al(:,2,2)*al(:,4,1)*
     &           al(:,3,3) - al(:,2,3)*al(:,3,1)*al(:,4,2)+
     &           al(:,2,3)*al(:,4,1)*al(:,3,2)
            
            detaij = 1./detaij
            
            zi0(:,1) = detaij(:)*((al(:,4,2)*al(:,3,3)
     &           - al(:,3,2)*al(:,4,3))*(xts1-al(:,1,1)) +
     &           (al(:,3,1)*al(:,4,3)
     &           - al(:,4,1)*al(:,4,3))*(xts2-al(:,1,2)) +
     &           (al(:,4,1)*al(:,3,2)
     &           - al(:,3,1)*al(:,4,2))*(xts3-al(:,1,3)))
            
            
            zi0(:,2) = detaij(:)*((al(:,2,2)*al(:,4,3)
     &           - al(:,4,2)*al(:,2,3))*(xts1-al(:,1,1)) +
     &           (al(:,4,1)*al(:,2,3)
     &           - al(:,2,1)*al(:,4,3))*(xts2-al(:,1,2)) +
     &           (al(:,2,1)*al(:,4,2)
     &           - al(:,4,1)*al(:,2,2))*(xts3-al(:,1,3)))
            
            zi0(:,3) = detaij(:)*((al(:,3,2)*al(:,2,3)
     &           - al(:,2,2)*al(:,3,3))*(xts1-al(:,1,1)) +
     &           (al(:,2,1)*al(:,3,3)
     &           - al(:,3,1)*al(:,2,3))*(xts2-al(:,1,2)) +
     &           (al(:,3,1)*al(:,2,2)
     &           - al(:,2,1)*al(:,3,2))*(xts3-al(:,1,3)))
            
            
c...  iterate to convergence
            
            do it = 1, iterat
               
c...  build matrix
               
               m11(:)=al(:,2,1)+al(:,5,1)*zi0(:,2)+al(:,7,1)*zi0(:,3)
     &              +al(:,8,1)*zi0(:,2)*zi0(:,3)
               m12(:)=al(:,3,1)+al(:,5,1)*zi0(:,1)+al(:,6,1)*zi0(:,3)
     &              +al(:,8,1)*zi0(:,1)*zi0(:,3)
               m13(:)=al(:,4,1)+al(:,6,1)*zi0(:,2)+al(:,7,1)*zi0(:,1)
     &              +al(:,8,1)*zi0(:,1)*zi0(:,2)
               
               m21(:)=al(:,2,2)+al(:,5,2)*zi0(:,2)+al(:,7,2)*zi0(:,3)
     &              +al(:,8,2)*zi0(:,2)*zi0(:,3)
               m22(:)=al(:,3,2)+al(:,5,2)*zi0(:,1)+al(:,6,2)*zi0(:,3)
     &              +al(:,8,2)*zi0(:,1)*zi0(:,3)
               m23(:)=al(:,4,2)+al(:,6,2)*zi0(:,2)+al(:,7,2)*zi0(:,1)
     &              +al(:,8,2)*zi0(:,1)*zi0(:,2)
               
               m31(:)=al(:,2,3)+al(:,5,3)*zi0(:,2)+al(:,7,3)*zi0(:,3)
     &              +al(:,8,3)*zi0(:,2)*zi0(:,3)
               m32(:)=al(:,3,3)+al(:,5,3)*zi0(:,1)+al(:,6,3)*zi0(:,3)
     &              +al(:,8,3)*zi0(:,1)*zi0(:,3)
               m33(:)=al(:,4,3)+al(:,6,3)*zi0(:,2)+al(:,7,3)*zi0(:,1)
     &              +al(:,8,3)*zi0(:,1)*zi0(:,2)
               
               
c...  build rhs
               
               r1(:)=al(:,1,1)+al(:,2,1)*zi0(:,1)+al(:,3,1)*zi0(:,2)+
     &              al(:,4,1)*zi0(:,3)+al(:,5,1)*zi0(:,1)*zi0(:,2)+
     &              al(:,6,1)*zi0(:,2)*zi0(:,3)+al(:,7,1)*
     &              zi0(:,1)*zi0(:,3)+al(:,8,1)*zi0(:,1)*
     &              zi0(:,2)*zi0(:,3) - xts1
               
               r2(:)=al(:,1,2)+al(:,2,2)*zi0(:,1)+al(:,3,2)*zi0(:,2)+
     &              al(:,4,2)*zi0(:,3)+al(:,5,2)*zi0(:,1)*zi0(:,2)+
     &              al(:,6,2)*zi0(:,2)*zi0(:,3)+al(:,7,2)*
     &              zi0(:,1)*zi0(:,3)+al(:,8,2)*zi0(:,1)*
     &              zi0(:,2)*zi0(:,3) - xts2
               
               r3(:)=al(:,1,3)+al(:,2,3)*zi0(:,1)+al(:,3,3)*zi0(:,2)+
     &              al(:,4,3)*zi0(:,3)+al(:,5,3)*zi0(:,1)*zi0(:,2)+
     &              al(:,6,3)*zi0(:,2)*zi0(:,3)+al(:,7,3)*
     &              zi0(:,1)*zi0(:,3)+al(:,8,3)*zi0(:,1)*
     &              zi0(:,2)*zi0(:,3) - xts3
               
               
c...  get solution
               
               detaij = m11*m22*m33-m11*m23*m32-m21*m12*m33+
     &              m21*m13*m32+m31*m12*m23-m31*m13*m22  
               
               detaij = 1./detaij
               
               dzi0(:,1) = -detaij(:)*((m22(:)*m33(:)-m23(:)*m32(:))*
     &              r1(:) + (m13(:)*m32(:)-m12(:)*m33(:))*r2(:) + 
     &              (m12(:)*m23(:)-m13(:)*m22(:))*r3(:))
               dzi0(:,2) = -detaij(:)*((m23(:)*m31(:)-m21(:)*m32(:))*
     &              r1(:) + (m11(:)*m33(:)-m13(:)*m31(:))*r2(:) + 
     &              (m13(:)*m21(:)-m11(:)*m23(:))*r3(:))
               dzi0(:,3) = -detaij(:)*((m21(:)*m32(:)-m22(:)*m31(:))*
     &              r1(:) + (m12(:)*m31(:)-m11(:)*m32(:))*r2(:) + 
     &              (m11(:)*m22(:)-m12(:)*m21(:))*r3(:))
               
               zi0(:,:) = zi0(:,:) + dzi0(:,:)
               
            enddo
            
         
            
            do e = 1, npro
               
               if ((abs(zi0(e,1)).lt.(one+tolpt)).and.(abs(zi0(e,2)).lt.
     &              (one+tolpt)).and.(abs(zi0(e,3)).lt.(one+tolpt))) 
     &              then        ! got the element
               
                  call shphex (ipord, zi0(e,:),shape(e,:),
     &                 shgradl(e,:,:))
                  
                  soln = zero
                  do i = 1, nenl
                     soln = soln+ycl(e,i,varcod)*shape(e,i)
                  enddo
                  do i = nenl+1,nshl
                     soln = soln+ycl(e,i,varcod)*shape(e,i)*sgn(e,i)
                  enddo
                  varts(jj) = soln
               
               endif
            
            enddo

         
         elseif (lcsyst.eq.1) then !tet

            call get_a_not_tet(xl,al)
            
            detaij(:) = -al(:,2,1)*al(:,3,2)*al(:,4,3) + 
     &           al(:,2,1)*al(:,4,2)*al(:,3,3) + al(:,2,2)*
     &           al(:,3,1)*al(:,4,3) - al(:,2,2)*al(:,4,1)*
     &           al(:,3,3) - al(:,2,3)*al(:,3,1)*al(:,4,2)+
     &           al(:,2,3)*al(:,4,1)*al(:,3,2)
            
            detaij = 1./detaij
           
c
c solve for r, s, t  for all the elements
c 
            zi0(:,1) = detaij(:)*((al(:,4,2)*al(:,3,3)
     &           - al(:,3,2)*al(:,4,3))*(xts1-al(:,1,1)) +
     &           (al(:,3,1)*al(:,4,3)
     &           - al(:,4,1)*al(:,4,3))*(xts2-al(:,1,2)) +
     &           (al(:,4,1)*al(:,3,2)
     &           - al(:,3,1)*al(:,4,2))*(xts3-al(:,1,3)))
            
            
            zi0(:,2) = detaij(:)*((al(:,2,2)*al(:,4,3)
     &           - al(:,4,2)*al(:,2,3))*(xts1-al(:,1,1)) +
     &           (al(:,4,1)*al(:,2,3)
     &           - al(:,2,1)*al(:,4,3))*(xts2-al(:,1,2)) +
     &           (al(:,2,1)*al(:,4,2)
     &           - al(:,4,1)*al(:,2,2))*(xts3-al(:,1,3)))
            
            zi0(:,3) = detaij(:)*((al(:,3,2)*al(:,2,3)
     &           - al(:,2,2)*al(:,3,3))*(xts1-al(:,1,1)) +
     &           (al(:,2,1)*al(:,3,3)
     &           - al(:,3,1)*al(:,2,3))*(xts2-al(:,1,2)) +
     &           (al(:,3,1)*al(:,2,2)
     &           - al(:,2,1)*al(:,3,2))*(xts3-al(:,1,3)))
            
            
            do e = 1, npro
               
               if (zi0(e,1).lt.(one+tolpt).and.zi0(e,1).gt.(zero-tolpt)
     &              .and.zi0(e,2).lt.(one+tolpt).and.zi0(e,2).gt.(zero
     &              -tolpt).and.zi0(e,3).lt.(one+tolpt).and.zi0(e,3)
     &              .gt.(zero-tolpt)) then
               
                  call shptet (ipord, zi0(e,:),shape(e,:),
     &                 shgradl(e,:,:))
                  
                  soln = zero
                  
                  if (varcod.eq.1) then !pres 
                     do i = 1, nshl
                        soln = soln+ycl(e,i,1)*shape(e,i)
                     enddo
                     varts(jj) = soln
                  elseif (varcod.eq.2) then !u 
                     do i = 1, nshl
                        soln = soln+ycl(e,i,2)*shape(e,i)
                     enddo
                     varts(jj) = soln
                  elseif (varcod.eq.3) then !v 
                     do i = 1, nshl
                        soln = soln+ycl(e,i,3)*shape(e,i)
                     enddo
                     varts(jj) = soln  
                  elseif (varcod.eq.4) then !w 
                     do i = 1, nshl
                        soln = soln+ycl(e,i,4)*shape(e,i)
                     enddo
                     varts(jj) = soln
                  elseif (varcod.eq.5) then !t 
                     do i = 1, nshl
                        soln = soln+ycl(e,i,5)*shape(e,i)
                     enddo
                     varts(jj) = soln
                  endif
                  
                  
               endif
            
            enddo
       
         endif
         
      enddo

      return 
      end
      
