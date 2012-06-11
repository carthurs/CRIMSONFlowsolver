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
      subroutine AsBMFG (u,          y,          ac,          
     &                   x,          xdist,      xdnv,
     &                   shpb,       shglb,
     &                   ienb,       materb,     iBCB,        BCB,
     &                   res,        xKebe,      SWB,         TWB,
     &                   EWB,        PS_global,  Kwall_xKebe)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c----------------------------------------------------------------------
c
        use turbSA                ! access to d2wall
        use deformableWall
        use LagrangeMultipliers 
c
        include "common.h"
c
        dimension y(nshg,ndofl),           
     &            ac(nshg,ndofl),          u(nshg,nsd),
     &            shpb(nshl,ngaussb),
     &            x(numnp,nsd),
     &            xdist(numnp),
     &            xdnv(numnp,nsd),
     &            shglb(nsd,nshl,ngaussb),         
     &            ienb(npro,nshl),         materb(npro),
     &            iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB),
     &            res(nshg,nflow),         dwl(npro,nenl),
     &            SWB(npro,nProps),        TWB(npro,2),
     &            EWB(npro,1)     
c
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd),
     &            rl(npro,nshl,nflow),     sgn(npro,nshl),
     &            ul(npro,nshl,nsd),       acl(npro,nshl,ndofl),
     &            xdistl(npro,nshl),       xdnvl(npro,nshl,nsd)
c
        dimension xKebe(npro,9,nshl,nshl) 
        
        real*8    Kwall_global(npro,9,9),  PS_global(npro,9),        
     &            Kwall_xKebe(npro,9,nshl,nshl)
     
        integer iblk
     
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localy(y,      yl,     ienb,   ndofl,  'gather  ')
        call localy(ac,     acl,    ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
        call localx(u,      ul,     ienb,   nsd,    'gather  ')
        
        call localx(xdist,     xdistl,    ienb,   1,      'gather  ')
        call localx(xdnv,      xdnvl,     ienb,   nsd,    'gather  ')
        
        if (iRANS.eq.-2) then
           call local(d2wall, dwl, ienb, 1, 'gather  ')
        endif

c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xKebe = zero
        endif   
        if(Lagrange.gt.zero) then
           loclhsLag = zero
        endif 
c
c.... get the boundary element residuals
c
        rl  = zero
        
        filrhsl = zero
c
c.... 3D
c
        call e3b  (ul,      yl,      acl,     
     &             iBCB,    BCB,     
     &             shpb,    shglb,
     &             xlb,     xdistl,  xdnvl,     
     &             rl,      sgn,     dwl,     xKebe,
     &             SWB,     TWB,     EWB,
     &             PS_global,        
     &             Kwall_xKebe)
c
c.... assemble the residual and the modified residual
c
        call local (res,    rl,     ienb,   nflow,  'scatter ')

c     
c.... end
c
        return
      end


c
c----------------------------------------------------------------------
c
c     This routine computes and assembles the data corresponding to the
c     boundary elements for the temperature equation
c
c----------------------------------------------------------------------
c
      subroutine AsBSclr (y,       x,       shpb,    shglb,
     &                   ienb,    materb,  iBCB,    BCB,
     &                   res)
        use turbSA ! access to d2wall
        include "common.h"
c
        dimension y(nshg,ndofl),           x(numnp,nsd),
     &            shpb(nshl,*),
     &            shglb(nsd,nshl,*),         
     &            ienb(npro,nshl),         materb(npro),
     &            iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB),
     &            res(nshg)         
c
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd),
     &            rl(npro,nshl),     sgn(npro,nshl)
        real*8 dwl(npro,nshl)
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localy(y,      yl,     ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
        if(iRANS.eq.-2) then
           call local(d2wall, dwl, ienb, 1, 'gather  ')
        endif
c
c.... get the boundary element residuals
c
        rl  = zero

        call e3bSclr  (yl,      iBCB,    BCB,     shpb,    shglb,
     &                 xlb,     rl,      sgn,     dwl)
c
c.... assemble the residual and the modified residual
c
        call local (res,    rl,     ienb,   1,  'scatter ')
c     
c.... end
c
        return
        end


