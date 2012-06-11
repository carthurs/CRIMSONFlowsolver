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
        subroutine AsIGMR (y,       ac,      x,       xmudmi,
     &                     shp,     shgl,    ien,     
     &                     res,     qres,
     &                     xKebe,   xGoC,    rerr, CFLworst)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c----------------------------------------------------------------------
c
      use stats
      use rlssave  ! Use the resolved Leonard stresses at the nodes.
      use timedata    ! time series
      use turbsa                ! access to d2wall
      use LagrangeMultipliers 


      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,ngauss),            shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            res(nshg,nflow),
     &            qres(nshg,idflx)

c
        dimension yl(npro,nshl,ndofl),         acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),           dwl(npro,nenl),      
     &            rl(npro,nshl,nflow), 
     &            ql(npro,nshl,idflx)
c        
        dimension xKebe(npro,9,nshl,nshl), 
     &            xGoC(npro,4,nshl,nshl)
c
        dimension rlsl(npro,nshl,6) 

c
        real*8    lStsVec(npro,nshl,nResDims)
        
        dimension xmudmi(npro,ngauss)
        dimension sgn(npro,nshl)
        dimension CFLworst(npro)
c
        real*8 rerrl(npro,nshl,6), rerr(nshg,10)
c
c.... gather the variables
c
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call local (qres,   ql,     ien,    idflx,  'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
 
        if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina 
           call local (rls, rlsl,     ien,       6, 'gather  ')  
        else
           rlsl = zero
        endif      

c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xKebe = zero
           xGoC  = zero
        endif   
        if(Lagrange.gt.zero) then
           loclhsLag = zero
        endif 
c
c.... get the element residuals, LHS matrix, and preconditioner
c
        rl     = zero

        if(ierrcalc.eq.1) rerrl = zero

        call e3  (yl,      acl,     dwl,     shp,
     &            shgl,    xl,      rl,      
     &            ql,      xKebe,   xGoC,    xmudmi, 
     &            sgn,     rerrl,  rlsl,     CFLworst)
c
c.... assemble the statistics residual
c
        if ( stsResFlg .eq. 1 ) then
           call e3StsRes ( xl, rl, lStsVec )
           call local( stsVec, lStsVec, ien, nResDims, 'scatter ')
        else
c
c.... assemble the residual
c
           call local (res,    rl,     ien,    nflow,  'scatter ')
           
           if ( ierrcalc .eq. 1 ) then
              call local (rerr, rerrl,  ien, 6, 'scatter ')
           endif
        endif
c
c.... end
c
        if (exts) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(yl,xl,ien,sgn)
           endif
        endif
        
        return
        end



C^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c-----------------------------------------------------------------------
c=======================================================================


        subroutine AsIGMRSclr(y,       ac,      x,       
     &                     shp,     shgl,    ien,     
     &                     res,     qres,    xSebe, xmudmi )
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  interior elements.
c
c----------------------------------------------------------------------
c
      use     turbSA  
      include "common.h"
c
        dimension y(nshg,ndofl),              ac(nshg,ndofl),
     &            x(numnp,nsd),              
     &            shp(nshl,ngauss),            shgl(nsd,nshl,ngauss),
     &            ien(npro,nshl),
     &            res(nshg),                  qres(nshg,nsd)

c
        real*8    yl(npro,nshl,ndofl),        acl(npro,nshl,ndofl),
     &            xl(npro,nenl,nsd),         
     &            rl(npro,nshl),              ql(npro,nshl,nsd),
     &            dwl(npro,nenl)            
c        
        real*8    xSebe(npro,nshl,nshl),      xmudmi(npro,ngauss) 
c
c.... gather the variables
c
        real*8 sgn(npro,nshl)
c
c.... get the matrix of mode signs for the hierarchic basis functions. 
c
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        if(iRANS.lt. 0) 
     &  call localx(d2wall, dwl,    ien,    1,      'gather  ')
        call local (qres,   ql,     ien,    nsd,    'gather  ')
c
c.... zero the matrices if they are being recalculated
c
        if (lhs. eq. 1)  then
           xSebe = zero
        endif   
c
c.... get the element residuals, LHS matrix, and preconditioner
c
      rl = zero
      call e3Sclr  (yl,      acl,     shp,
     &              shgl,    xl,      dwl,
     &              rl,      ql,      xSebe,   
     &              sgn, xmudmi)
c
c.... assemble the residual
c
        call local (res,    rl,     ien,    1,  'scatter ')
c
c.... end
c
        return
        end
