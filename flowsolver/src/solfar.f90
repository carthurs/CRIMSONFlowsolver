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
      subroutine SolFlow(y,          ac,         u,
     &                   yold,       acold,      uold,
     &                   x,          xdist,      xdnv,
     &                   iBC,        BC,         res,             
     &                   nPermDims,  nTmpDims,   aperm,
     &                   atemp,      iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       lhsK,       lhsP, 
     &                   solinc,     rerr  )
c
c----------------------------------------------------------------------
c
c This is the 2nd interface routine to the linear equation 
c solver library that uses the CGP and GMRES methods.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  acold   (nshg,ndof)          : Primvar. accel. at beginning of step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  res    (nshg,nflow)           : preconditioned residual
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c The followings are preliminary steps required to use the
c solver library.  New way of writing has to be used such as
c
c          |  K     G | | du |    | Rmom  |
c          |          | |    | =  |       |
c          | G^t    C | | dp |    | Rcon  |
c
c          |     E    | | dT | =  | Rtemp |
c
c     where
c
c      xKebe : K_ab = dRmom_a/du_b    xTe : E_ab = dRtemp_a/dT_b 
c
c              G_ab = dRmom_a/dp_b
c      xGoC  :
c              C_ab = dRcon_a/dp_b       
c
c              resf = Rmon Rcon       rest = Rtemp
c
c----------------------------------------------------------------------
c
      use pointer_data
      use LagrangeMultipliers 
c        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             
     &          xdist(numnp),
     &          xdnv(numnp,nsd),
     &          BC(nshg,ndofBC),
     &          res(nshg,nflow),
     &          flowDiag(nshg,4),
     &          aperm(nshg,nPermDims),    atemp(nshg,nTmpDims),
     &          sclrDiag(nshg,1),         
     &          lhsK(9,nnz_tot),	  lhsP(4,nnz_tot)          
c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg) 
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),         
     &          lesP(nshg,4),             lesQ(nshg,4),
     &          solinc(nshg,ndof)
      
      real*8    rerr(nshg,10),            rtmp(nshg,4)
c
c.... *******************>> Element Data Formation <<******************
c
c
c.... set the parameters for flux and surface tension calculations
c
c
      idflx = 0 
      if(idiff >= 1 )  idflx= (nflow-1) * nsd
      if (isurf == 1) idflx=nflow*nsd
c        
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,    yold,    acold,       
     &                u,       y,       ac,            
     &                uAlpha,  yAlpha,  acAlpha)

c
c.... form the LHS matrices, the residual vector (at alpha)
c
      call ElmGMR ( uAlpha,    yAlpha,     acAlpha,    
     &              x,         xdist,      xdnv,
     &              shp,       shgl,       iBC,       
     &              BC,        shpb,       shglb,
     &              res,       iper,       ilwork,   
     &              rowp,      colm,       lhsK,      
     &              lhsP,      rerr   )

c
c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1)
      eqnType = 1
c
c.... setup the linear algebra solver
c
      rtmp = res(:,1:4)
      call usrNew ( usr,        eqnType,          aperm,
     &              atemp,      rtmp,             solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDims,  
     &              nTmpDims,   rowp,             colm,     
     &              lhsK,       lhsP,             rdtmp,      
     &              nnz_tot )
c
c.... solve linear system
c
      call myfLesSolve ( lesId, usr )
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, nflow, 'out')
      endif

      if(Lagrange .gt. zero) then
         call CalcNANBLagrange(colm, rowp, solinc(:,1:3))
         call LagMultiplyMatrix(solinc, 0, nsrflistLagrange,
     &      numLagrangeSrfs)  
         Lagincr(:,1:3) = (- resL(:,1:3) - AddLag(:,1:3) )
     &      /ScaleFactor(1,1)/alfi/gami/two
      endif
      
      call rstatic (res, y, solinc) ! output flow stats
c     
c.... end
c     
      return
      end

      subroutine SolSclr(y,          ac,         u,
     &                   yold,       acold,      uold,
     &                   x,          iBC,
     &                   BC,         nPermDimsS,  nTmpDimsS,  
     &                   apermS,     atempS,     iper,       
     &                   ilwork,     shp,        shgl, 
     &                   shpb,       shglb,      rowp,     
     &                   colm,       lhsS,       solinc)
c
c----------------------------------------------------------------------
c
c This is the 2nd interface routine to the linear equation 
c solver library.
c
c input:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c  yold   (nshg,ndof)           : Y-variables at beginning of step
c  x      (numnp,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodic nodal information
c
c output:
c  y      (nshg,ndof)           : Y-variables at n+alpha_f
c  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
c
c
c The followings are preliminary steps required to use LesLib
c solver library.  New way of writing has to be used such as
c
c          |     E    | | dS | =  | RScal |
c
c----------------------------------------------------------------------
c
      use pointer_data
        
      include "common.h"
      include "mpif.h"
      include "auxmpi.h"
c     
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          yold(nshg,ndof),          acold(nshg,ndof),
     &          u(nshg,nsd),              uold(nshg,nsd),
     &          x(numnp,nsd),             BC(nshg,ndofBC),
     &          res(nshg,1),
     &          flowDiag(nshg,4),
     &          sclrDiag(nshg,1),           lhsS(nnz_tot),
     &          apermS(nshg,nPermDimsS),  atempS(nshg,nTmpDimsS)

c
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
      integer   usr(100),                 eqnType,
     &          rowp(nshg*nnz),           colm(nshg+1),
     &          iBC(nshg),                ilwork(nlwork),
     &          iper(nshg)
c
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof),
     &          uAlpha(nshg,nsd),
     &          lesP(nshg,1),               lesQ(nshg,1),
     &          solinc(nshg,1)
      
c     
c.... *******************>> Element Data Formation <<******************
c
c.... compute solution at n+alpha
c
      call itrYAlpha( uold,    yold,    acold, 
     &                u,       y,       ac,  
     &                uAlpha,  yAlpha,  acAlpha)
c
c.... form the LHS matrices, the residual vector (at alpha)
c
      call ElmGMRSclr (yAlpha,    acAlpha,    x,
     &                 shp,       shgl,       iBC,       
     &                 BC,        shpb,       shglb,
     &                 res,       iper,       ilwork,   
     &                 rowp,      colm,       lhsS   )

c
c.... lesSolve : main matrix solver
c
      lesId   = numeqns(1+nsolt+isclr)
      eqnType = 2
c
c.... setup the linear algebra solver
c
      call usrNew ( usr,        eqnType,          apermS,
     &              atempS,     res,              solinc,          
     &              flowDiag,   sclrDiag,         lesP,   
     &              lesQ,       iBC,              BC,
     &              iper,       ilwork,           numpe,
     &              nshg,       nshl,             nPermDimsS,  
     &              nTmpDimsS,  rowp,             colm,     
     &              rlhsK,      rlhsP,            lhsS,      
     &              nnz_tot )
c
c.... solve linear system
c
      call myfLesSolve ( lesId, usr )
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, 1, 'out')
      endif
      
      nsolsc=5+isclr
      call rstaticSclr (res, y, solinc, nsolsc) ! output scalar stats
c     
c.... end
c     
      return
      end





