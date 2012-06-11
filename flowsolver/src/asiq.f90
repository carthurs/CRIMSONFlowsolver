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
        subroutine AsIq (y,       x,       shp,
     &                   shgl,    ien,     xmudmi,
     &                   qres,    rmass    )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c input:
c     y     (numnp,ndof)        : Y variables
c     x     (numnp,nsd)         : nodal coordinates
c     shp   (nen,nintg)         : element shape-functions
c     shgl  (nsd,nen,nintg)     : element local shape-function gradients
c     ien   (npro)              : nodal connectivity array
c
c output:
c     qres  (numnp,nsd,nsd)  : residual vector for diffusive flux
c     rmass  (numnp)            : lumped mass matrix
c
c----------------------------------------------------------------------
c
        use turbsa      ! access to d2wall
        include "common.h"
c
        dimension y(nshg,ndof),               x(numnp,nsd),            
     &            shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),      dwl(npro,nenl),
     &            qres(nshg,idflx),    rmass(nshg)
c
        dimension yl(npro,nshl,ndof),          xl(npro,nenl,nsd),         
     &            ql(npro,nshl,idflx),  rmassl(npro,nshl),
     &            xmudmi(npro,ngauss)
c
        dimension sgn(npro,nshl)
c
c.... create the matrix of mode signs for the hierarchic basis 
c     functions. 
c
        do i=1,nshl
           where ( ien(:,i) < 0 )
              sgn(:,i) = -one
           elsewhere
              sgn(:,i) = one
           endwhere
        enddo

c
c.... gather the variables
c

        call localy(y,      yl,     ien,    ndof,   'gather  ')
        call localx (x,      xl,     ien,    nsd,    'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        call e3q  (yl,         dwl,      shp,      shgl,    
     &             xl,         ql,       rmassl,
     &             xmudmi,     sgn  )

c
c.... assemble the diffusive flux residual 
c
        call local (qres,   ql,  ien,  idflx,  'scatter ')
        call local (rmass,  rmassl, ien,  1,          'scatter ')
c
c.... end
c
        return
        end


c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c interior elements for the global reconstruction of the diffusive
c flux vector.
c
c----------------------------------------------------------------------
        subroutine AsIqSclr (y,       x,       shp,
     &                       shgl,    ien,     qres,    
     &                       rmass    )
c
        use turbsa      ! access to d2wall
        include "common.h"
c
        dimension y(nshg,ndof),             x(numnp,nsd),            
     &            shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),      dwl(npro,nenl),
     &            qres(nshg,nsd),           rmass(nshg)
c
        dimension yl(npro,nshl,ndof),       xl(npro,nenl,nsd),         
     &            ql(npro,nshl,nsd),        rmassl(npro,nshl)
c
        dimension sgn(npro,nshl)

        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
c
c.... gather the variables
c
        call localy(y,      yl,     ien,    ndof,   'gather  ')
        call localx (x,      xl,     ien,    nsd,    'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
c
c.... get the element residuals 
c
        ql     = zero
        rmassl = zero

        call e3qSclr  (yl,      dwl,    shp,    shgl,    
     &                 xl,      ql,     rmassl, 
     &                 sgn             )

c
c.... assemble the temperature diffusive flux residual 
c
        call local (qres,   ql,  ien,  nsd,  'scatter ')
        call local (rmass,  rmassl, ien,  1, 'scatter ')
c
c.... end
c
        return
        end

