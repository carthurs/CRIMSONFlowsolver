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
        subroutine e3qvar (yl,          shgl,    
     &                     xl,          g1yi,
     &                     g2yi,        g3yi,        shg,
     &                     dxidx,       WdetJ )
c
c----------------------------------------------------------------------
c
c  This routine computes the variables at integration point
c  necessary for the computation of the diffusive flux vector.
c
c input:
c  yl     (npro,nshl,ndof)      : primitive variables
c  shgl   (npro,nsd,nshl)     : element local-grad-shape-functions
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step
c
c output:
c  g1yi   (npro,ndof)           : grad-y in direction 1
c  g2yi   (npro,ndof)           : grad-y in direction 2
c  g3yi   (npro,ndof)           : grad-y in direction 3
c  shg    (npro,nshl,nsd)       : element global grad-shape-functions
c  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient
c  WdetJ  (npro)                : weighted Jacobian
c  u1     (npro)                : x1-velocity component
c  u2     (npro)                : x2-velocity component
c  u3     (npro)                : x3-velocity component
c
c----------------------------------------------------------------------
c
        include "common.h"
c
c  passed arrays
c
        dimension yl(npro,nshl,ndof), 
     &            shgl(npro,nsd,nshl), xl(npro,nenl,nsd),
     &            g1yi(npro,nflow),       g2yi(npro,nflow),
     &            g3yi(npro,nflow),       shg(npro,nshl,nsd), 
     &            dxidx(npro,nsd,nsd),   WdetJ(npro)
c
c  local arrays
c
        dimension tmp(npro),           dxdxi(npro,nsd,nsd)

c
c.... compute the deformation gradient
c
        dxdxi = zero
c
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(:,1,n)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(:,2,n)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(:,3,n)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(:,1,n)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(:,2,n)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(:,3,n)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(:,1,n)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(:,2,n)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(:,3,n)
          enddo
c
c.... compute the inverse of deformation gradient
c
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3) 
     &                 - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3) 
     &                 - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3) 
     &                 - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1) 
     &                       + dxidx(:,1,2) * dxdxi(:,2,1)  
     &                       + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1) 
     &                - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3) 
     &                - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3) 
     &                - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2) 
     &                - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2) 
     &                - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2) 
     &                - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
c
        WdetJ = Qwt(lcsyst,intp)/ tmp

c
c.... --------------------->  Global Gradients  <-----------------------
c
        g1yi = zero
        g2yi = zero
        g3yi = zero
c
c
        do n = 1, nshl
c
c.... compute the global gradient of shape-function
c
c            ! N_{a,x_i}= N_{a,xi_i} xi_{i,x_j}
c
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) + 
     &                 shgl(:,2,n) * dxidx(:,2,1) +
     &                 shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) + 
     &                 shgl(:,2,n) * dxidx(:,2,2) +
     &                 shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) + 
     &                 shgl(:,2,n) * dxidx(:,2,3) +
     &                 shgl(:,3,n) * dxidx(:,3,3) 
c
c.... compute the global gradient of Y-variables
c
c
c  Y_{,x_i}=SUM_{a=1}^nenl (N_{a,x_i}(int) Ya)
c
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
c
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
c
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)

       enddo

c
c.... return
c

       return
       end

c-----------------------------------------------------------------------
c
c     compute the variables for the local scalar diffusion
c
c-----------------------------------------------------------------------
      subroutine e3qvarSclr  (yl,       shgl,         xl, 
     &                        gradT,    dxidx,        WdetJ )

      include "common.h"
c
c  passed arrays
c
      real*8   yl(npro,nshl,ndof),    shp(npro,nshl),
     &         shgl(npro,nsd,nshl),   xl(npro,nenl,nsd),
     &         dxidx(npro,nsd,nsd),   WdetJ(npro),
     &         gradT(npro,nsd)
c
c  local arrays
c
      real*8   shg(npro,nshl,nsd)


      call e3metric( xl,         shgl,       dxidx,  
     &               shg,        WdetJ)

      gradT = zero
      id=5+isclr
c
c  later, when there are more models than SA we will need a 
c  more general function to calculate evisc at a quadrature point
c
      do n = 1, nshl
         gradT(:,1) = gradT(:,1) + shg(:,n,1) * yl(:,n,id)
         gradT(:,2) = gradT(:,2) + shg(:,n,2) * yl(:,n,id)
         gradT(:,3) = gradT(:,3) + shg(:,n,3) * yl(:,n,id)
      enddo
c
c.... return
c

       return
       end
       
            


