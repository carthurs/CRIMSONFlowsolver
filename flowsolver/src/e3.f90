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
        subroutine e3 (yl,      acl,     dwl,     shp,
     &                 shgl,    xl,      rl,      ql,
     &                 xKebe,   xGoC,    xmudmi,  sgn, 
     &                 rerrl, rlsl)
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     UBar formulation of the incompressible Navier Stokes equations.
c
c
c input:    e    a   1..5   when we think of U^e_a  and U is 5 variables
c  yl     (npro,nshl,ndof)      : Y variables (not U)
c  acl    (npro,nshl,ndof)      : Y acceleration (Y_{,t})
c  shp    (nen,ngauss)           : element shape-functions  N_a
c  shgl   (nsd,nen,ngauss)       : element local-shape-functions N_{a,xi}
c  wght   (ngauss)               : element weight (for quadrature)
c  xl     (npro,nenl,nsd)       : nodal coordinates at current step (x^e_a)
c  ql     (npro,nshl,nsd*nsd) : diffusive flux vector (don't worry)
c  rlsl   (npro,nshl,6)       : resolved Leonard stresses
c
c output:
c  rl     (npro,nshl,nflow)      : element RHS residual    (G^e_a)
c  rml    (npro,nshl,nflow)      : element modified residual  (G^e_a tilde)
c  xKebe  (npro,9,nshl,nshl)  : element LHS tangent mass matrix
c  xGoC   (npro,4,nshl,nshl)    : element LHS tangent mass matrix
c
c Note: This routine will calculate the element matrices for the
c        Hulbert's generalized alpha method integrator
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof),
     &            acl(npro,nshl,ndof),       
     &            shp(nshl,ngauss),       shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),      dwl(npro,nenl),
     &            rl(npro,nshl,nflow),     ql(npro,nshl,idflx)
c      

      	dimension xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
c
c.... local declarations
c
        dimension g1yi(npro,ndof),        g2yi(npro,ndof),
     &            g3yi(npro,ndof),        shg(npro,nshl,nsd),
     &            aci(npro,3),            dxidx(npro,nsd,nsd),       
     &            WdetJ(npro),            rho(npro),
     &            pres(npro),             u1(npro),
     &            u2(npro),               u3(npro),
     &            rLui(npro,nsd),         uBar(npro,nsd),
     &            xmudmi(npro,ngauss),     sgn(npro,nshl), 
     &            shpfun(npro,nshl),      shdrv(npro,nsd,nshl),
     &            rmu(npro),              tauC(npro),
     &            tauM(npro),             tauBar(npro),
     &            src(npro,3)

        dimension rlsl(npro,nshl,6),      rlsli(npro,6)

        real*8    rerrl(npro,nshl,6)
        integer   aa

c
c     
c.... local reconstruction of diffusive flux vector for quadratics
c     or greater but NOT for bflux since local mass was not mapped
c
        if ( idiff==2 .and. ires .eq. 1 ) then
           call e3ql (yl,        dwl,       shp,       shgl, 
     &                xl,        ql,        xmudmi, 
     &                sgn)
        endif
c
c.... loop through the integration points
c

        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shp,          shgl,      sgn, 
     &              shpfun,       shdrv)
c
c.... get necessary fluid properties (including eddy viscosity)
c
        call getdiff(dwl,  yl,     shpfun,     xmudmi, xl,   rmu, rho)
c
c.... calculate the integration variables
c
        call e3ivar (yl,          acl,       shpfun,
     &               shdrv,       xl,
     &               aci,         g1yi,      g2yi,    
     &               g3yi,        shg,       dxidx,   
     &               WdetJ,       rho,       pres, 
     &               u1,          u2,        u3,              
     &               ql,          rLui,      src,
     &               rerrl,       rlsl,      rlsli,
     &               dwl) 
c
c.... compute the stabilization terms
c
        call e3stab (rho,          u1,       u2,
     &               u3,           dxidx,    rLui,   
     &               rmu,          tauC,     tauM,   
     &               tauBar,       uBar )  
c
c.... compute the residual contribution at this integration point
c
        call e3Res ( u1,        u2,         u3,
     &               uBar,      aci,        WdetJ,
     &               g1yi,      g2yi,       g3yi,
     &               rLui,      rmu,        rho,
     &               tauC,      tauM,       tauBar,
     &               shpfun,    shg,        src,
     &               rl,        pres,       acl,
     &               rlsli)
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHS ( u1,        u2,         u3,
     &                  uBar,      WdetJ,      rho,
     &                  rLui,      rmu,
     &                  tauC,      tauM,       tauBar,
     &                  shpfun,    shg,        xKebe,
     &                  xGoC )
        endif

c
c.... end of integration loop
c
      enddo

c
c.... symmetrize C
c
      if (lhs .eq. 1) then
         do ib = 1, nshl
            do iaa = 1, ib-1
               xGoC(:,4,iaa,ib) = xGoC(:,4,ib,iaa)
            enddo
         enddo
      endif
c
c.... return
c
      return
      end


c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c###################################################################


      subroutine e3Sclr (yl,      acl,     shp,
     &                     shgl,    xl,      dwl,
     &                     rl,      ql,      xSebe,   
     &                     sgn,     xmudmi)
c                                                                      
c----------------------------------------------------------------------
c
c     This routine calculates the residual and tangent matrix for the 
c     advection - diffusion equation for scalar.
c
c----------------------------------------------------------------------
c
      include "common.h"
c
      real*8    yl(npro,nshl,ndof),     acl(npro,nshl,ndof),       
     &            shp(nshl,ngauss),       shgl(nsd,nshl,ngauss),
     &            xl(npro,nenl,nsd),      rl(npro,nshl),          
     &            ql(npro,nshl,nsd),      xSebe(npro,nshl,nshl),
     &            dwl(npro,nshl)
c
c.... local declarations
c
      real*8    gradS(npro,nsd),        shg(npro,nshl,nsd),
     &            Sdot(npro),             Sclr(npro),
     &            dxidx(npro,nsd,nsd),    WdetJ(npro),      
     &            u1(npro),     u2(npro), u3(npro),
     &            sgn(npro,nshl),         shpfun(npro,nshl),       
     &            shdrv(npro,nsd,nshl),   rLS(npro),
     &            tauS(npro),             diffus(npro),
     &            srcL(npro),             srcR(npro),
     &            gGradS(npro,nsd),       dcFct(npro),
     &            giju(npro,6)
c
c.... Source terms sometimes take the form (beta_i)*(phi,_i).  Since
c     the convective term has (u_i)*(phi,_i), it is useful to treat
c     beta_i as a "correction" to the velocity.  In calculating the
c     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
c     then used in place of the pure velocity for stabilization terms,
c     and the source term sneaks into the RHS and LHS.
      real*8    uMod(npro,nsd), srcRat(npro), xmudmi(npro,ngauss)
c
      integer   aa, b
c     
c.... local reconstruction of diffusive flux vector
c
        if ( idiff==2 ) then
           call e3qlSclr (yl, dwl, shp, shgl, xl, ql, sgn)
        endif
c
c.... loop through the integration points
c
        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shp,          shgl,      sgn, 
     &              shpfun,        shdrv)
c
c.... get necessary fluid properties
c
        call getdiffsclr(shpfun,dwl,yl,diffus)
c
c.... calculate the integration variables
c
        call e3ivarSclr(yl,          acl,       shpfun,
     &                  shdrv,       xl,        xmudmi,
     &                  Sclr,        Sdot,      gradS,
     &                  shg,         dxidx,     WdetJ,       
     &                  u1,          u2,        u3,              
     &                  ql,          rLS,       SrcR,
     &                  SrcL,        uMod,      dwl,
     &                  diffus,      srcRat)


c
c.... compute the stabilization terms
c
        call e3StabSclr (uMod,    dxidx,   tauS, 
     &                   diffus,  srcR,    giju,
     &                   srcRat)
c
c... computing the DC factor for the discontinuity capturing
c
        if (idcsclr(1) .ne. 0) then
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. 
     &          (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
c
              call e3dcSclr ( gradS,    giju,     gGradS,
     &                        rLS,      tauS,     srcR,
     &                        dcFct)
           endif
        endif                   !end of idcsclr
c
c.... compute the residual contribution at this integration point
c
        call e3ResSclr ( uMod,      gGradS,
     &                   Sclr,      Sdot,       gradS,  
     &                   WdetJ,     rLS,        tauS,
     &                   shpfun,    shg,        srcR,
     &                   diffus, 
     &                   rl )
c
c.... compute the tangent matrix contribution
c
        if (lhs .eq. 1) then
           call e3LHSSclr ( uMod,      giju,       dcFct,
     &                      Sclr,      Sdot,       gradS,  
     &                      WdetJ,     rLS,        tauS,
     &                      shpfun,    shg,        srcL,
     &                      diffus,
     &                      xSebe )

        endif

c
c.... end of integration loop
c
      enddo

c
c.... return
c
      return
      end
