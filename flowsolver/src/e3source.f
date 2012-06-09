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
      subroutine e3source(xx, src)
c-----------------------------------------------------------------------
c
c  this routine computes the body force term.
c
c  currently this computes a swirl body with the axis alligned with 
c  the z-coordinate
c
c-----------------------------------------------------------------------

      include "common.h"
      
      real*8   xx(npro,nsd), src(npro,nsd)

      real*8   nu

      real*8   r, Stheta, dpdz, rP5
      if(datmat(2,5,1).eq. 0) then  ! calc swirl 
c
c  This is the body force which will drive a swirl in a pipe flow
c     
         bigR    = 0.5d0
         dpdz    = datmat(1,5,1)
         do iel = 1, npro
            
            r   = sqrt( xx(iel,1)**2 + xx(iel,2)**2)
            rP5 = (r/bigR)**5
            
            Stheta = dpdz * sin(0.5*pi*rP5)
            
            src(iel,1) = -xx(iel,2)/r * Stheta
            src(iel,2) =  xx(iel,1)/r * Stheta
            src(iel,3) =  dpdz
         enddo
      else  ! contrived test problem

c$$$         nu=datmat(1,2,1)
c$$$         omag=datmat(3,5,1) ! frame rotation rate
c$$$         e1 = one/sqrt(3.0d0) ! axis of rotation
c$$$         e2 = e1
c$$$         e3 = e1
c$$$         om1=omag*e1
c$$$         om2=omag*e2
c$$$         om3=omag*e3
c$$$         w=0
c$$$      
c$$$         do iel = 1, npro
c$$$         
c$$$            x = xx(iel,1)
c$$$            y = xx(iel,2)
c$$$            z = xx(iel,3)
c$$$
c$$$c
c$$$c  The following are MAPLE generated forcing functions for
c$$$c  a contrived flow in a rotating reference frame.
c$$$c
c$$$
c$$$      t1 = x**2
c$$$      t2 = t1**2
c$$$      t5 = dexp(14.D0*x)
c$$$      t6 = t2*x*t5
c$$$      t7 = y**2
c$$$      t8 = t7**2
c$$$      t9 = t8*t7
c$$$      t12 = t2*t5
c$$$      t15 = nu*t1
c$$$      t17 = dexp(7.D0*x)
c$$$      t18 = t7*y
c$$$      t19 = t17*t18
c$$$      t22 = t1*x
c$$$      t23 = nu*t22
c$$$      t24 = t17*y
c$$$      t29 = t17*t7
c$$$      t34 = nu*x
c$$$      t43 = nu*t2
c$$$      t46 = -40.D0*t6*t9-6.D0*t12*t7+92.D0*t15*t19+132.D0*t23*t24+80.D0*
c$$$     #t6*t18-252.D0*t23*t29+168.D0*t23*t19+96.D0*t34*t29-64.D0*t34*t19+2
c$$$     #2.D0*t15*t24-138.D0*t15*t29+294.D0*t43*t29
c$$$      t50 = t2*t22*t5
c$$$      t57 = nu*t17
c$$$      t60 = t8*y
c$$$      t65 = t2**2
c$$$      t66 = t65*t5
c$$$      t73 = t22*t5
c$$$      t77 = t2*t1*t5
c$$$      t80 = -196.D0*t43*t19-96.D0*t50*t9-122.D0*t43*t24-32.D0*t34*t24-4.
c$$$     #D0*t57*y+36.D0*t12*t60+192.D0*t50*t18+98.D0*t66*t8+24.D0*t12*t18+2
c$$$     #8.D0*t66*t9-24.D0*t73*t60-336.D0*t77*t60
c$$$      t106 = -336.D0*t50*t8-12.D0*t12*t9+8.D0*t73*t9-140.D0*t6*t8+28.D0*
c$$$     #t73*t8+12.D0*t15*t17+12.D0*t43*t17+288.D0*t50*t60+112.D0*t77*t9-22
c$$$     #4.D0*t77*t18-56.D0*t66*t18+120.D0*t6*t60
c$$$      t131 = -8.D0*t57*t18-20.D0*t6*t7+56.D0*t77*t7-42.D0*t12*t8-24.D0*t
c$$$     #23*t17-48.D0*t50*t7+4.D0*t73*t7-16.D0*t73*t18+392.D0*t77*t8+14.D0*
c$$$     #t66*t7-84.D0*t66*t60+12.D0*t57*t7
c$$$      fx = t46+t80+t106+t131
c$$$
c$$$
c$$$      t1 = x**2
c$$$      t2 = t1**2
c$$$      t3 = nu*t2
c$$$      t5 = dexp(7.D0*x)
c$$$      t6 = t5*y
c$$$      t9 = y**2
c$$$      t10 = t9**2
c$$$      t11 = nu*t10
c$$$      t18 = t1*x
c$$$      t22 = nu*x
c$$$      t25 = nu*t18
c$$$      t32 = dexp(14.D0*x)
c$$$      t33 = t2*t32
c$$$      t39 = t2*t1*t32
c$$$      t40 = t10*t9
c$$$      t43 = t1*t32
c$$$      t46 = t9*y
c$$$      t47 = t10*t46
c$$$      t50 = -84.D0*t3*t6+343.D0*t11*t5*t2+66.D0*t11*t5*x-98.D0*t11*t5*t1
c$$$     #8-24.D0*t22*t6+120.D0*t25*t6-287.D0*t11*t5*t1-140.D0*t33*t10+14.D0
c$$$     #*t3*t5-56.D0*t39*t40-28.D0*t43*t40+8.D0*t43*t47
c$$$      t52 = t2*x*t32
c$$$      t55 = t18*t32
c$$$      t62 = t10*y
c$$$      t69 = nu*t5
c$$$      t80 = 120.D0*t52*t10-32.D0*t55*t47+56.D0*t33*t47+30.D0*t11*t5-216.
c$$$     #D0*t52*t62-144.D0*t55*t62+72.D0*t39*t62-60.D0*t69*t46+30.D0*t69*t9
c$$$     #-20.D0*t25*t5-20.D0*t43*t10+36.D0*t43*t62
c$$$      t86 = nu*t1
c$$$      t89 = t5*t9
c$$$      t92 = t5*t46
c$$$      t109 = 112.D0*t55*t40+80.D0*t55*t10-12.D0*t86*t6-275.D0*t86*t89+57 
c$$$     #4.D0*t86*t92-218.D0*t25*t89+196.D0*t25*t92+90.D0*t22*t89-132.D0*t2
c$$$     #2*t92+427.D0*t3*t89+8.D0*t39*t46+4.D0*t43*t46
c$$$      t134 = 16.D0*t39*t47+28.D0*t33*t46-48.D0*t52*t47+252.D0*t33*t62+16
c$$$     #8.D0*t52*t40-24.D0*t52*t46-686.D0*t3*t92+2.D0*t86*t5-40.D0*t39*t10
c$$$     #-196.D0*t33*t40-16.D0*t55*t46+4.D0*t22*t5
c$$$      fy = t50+t80+t109+t134
c$$$
c$$$
c$$$      t3 = om3*x  
c$$$      t5 = dexp(7.D0*x)
c$$$      t7 = x**2
c$$$      t12 = y**2
c$$$      t15 = (-1.D0+y)**2
c$$$      fxa = 2.D0*om2*w+2.D0*t3*t5*(2.D0+x-10.D0*t7+7.D0*t7*x)*t12*t15+om
c$$$     #2*(om1*y-om2*x)-om3*(t3-om1*z)
c$$$
c$$$      t1 = x**2
c$$$      t4 = (-1.D0+x)**2
c$$$      t7 = dexp(7.D0*x)
c$$$      t10 = y**2
c$$$      fya = 4.D0*om3*t1*t4*t7*y*(1.D0-3.D0*y+2.D0*t10)-2.D0*om1*w+om3*(o
c$$$     #m2*z-om3*y)-om1*(om1*y-om2*x)
c$$$
c$$$
c$$$      t3 = dexp(7.D0*x)
c$$$      t5 = x**2
c$$$      t10 = y**2
c$$$      t13 = (-1.D0+y)**2
c$$$      t19 = (-1.D0+x)**2
c$$$      fza = -2.D0*om1*x*t3*(2.D0+x-10.D0*t5+7.D0*t5*x)*t10*t13-4.D0*om2*
c$$$     #t5*t19*t3*y*(1.D0-3.D0*y+2.D0*t10)+om1*(om3*x-om1*z)-om2*(om2*z-om
c$$$     #3*y)
c$$$
c$$$      src(iel,1) =  fx + fxa
c$$$      src(iel,2) =  fy + fya
c$$$      src(iel,3) =       fza
c$$$      enddo
!Analytic lid case
   
      do iel = 1, npro
            x = xx(iel,1)
            y = xx(iel,2)
            z = xx(iel,3)
            Re = 1.0/datmat(1,2,1)
c
c  The following are MAPLE generated forcing functions for
c  a Lid Driven cavity flow with an analytic solution
c
            t2 = x**2
            t3 = t2**2
            t4 = t3*x
            t7 = t2*x
            t8 = 8*t7
            t13 = y**2
            t20 = t13**2
            t21 = t20-t13
            t26 = t3**2
            t30 = t3*t2
            t37 = t13*y
            t54 = -8/Re*(24.E0/5.E0*t4-12*t3+t8+2*(4*t7-6*t2+2*x)*(12
     &           *t13-2)+(24*x-12)*t21)-64*(t26/2-2*t3*t7+3*t30-2*t4+t3
     &           /2)*(-24*t20*y+8*t37-4*y)+64*t21*(4*t37-2*y)*(-4*t30+12
     &           *t4-14*t3+t8-2*t2)

            src(iel,1) = 0.0
            src(iel,2) = t54
            src(iel,3) = 0.0
         enddo
   
      endif
      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!******************************************************************
!-------------------------------------------------------------------



      subroutine e3sourceSclr ( Sclr,      Sdot,      gradS,
     &                          dwl,       shape_funct,     shg,   
     &                          yl,        dxidx,     rmu,
     &                          u1,        u2,        u3,
     &                          srcR,      srcL,      uMod,
     &                          srcRat) 


c-----------------------------------------------------------------------
c
c     calculate the coefficients for the Spalart-Allmaras 
c     turbulence model and its jacobian
c
c     input: 
c             Sclr  :   The scalar that is being transported in this pde.
c             gradS :   spatial gradient of Sclr
c             rmu   :   kinematic viscosity
c
c     output:
c             rmu   :   diffusion for eddy viscosity equation
c             srcR  :   residual terms for  
c                         (srcR * turb)
c             srcL :   jacobian of srcR
c
c-----------------------------------------------------------------------
      use     turbSA
      use turbKE
      include "common.h"
c coming in      
      real*8  Sclr(npro),          Sdot(npro),
     &        gradS(npro,nsd),     dwl(npro,nenl),
     &        shape_funct(npro,nshl),    shg(npro,nshl,nsd),
     &        yl(npro,nenl,ndof),  dxidx(npro,nsd,nsd),
     &        rmu(npro),           u1(npro),
     &        u2(npro),            u3(npro)
c going out
      real*8  srcR(npro),          srcL(npro),
     &        uMod(npro,nsd)
c used locally

      real*8  gradv(npro,nsd,nsd), absVort(npro),
     &        dwall(npro)
      real*8  chi,    chiP3,   fv1,   fv2,     st,    r,
     &        g,      fw,      s,     viscInv, k2d2Inv,
     &        dP2Inv, sixth,   tmp,   tmp1,    p,     dp,
     &        fv1p,   fv2p,    stp,         gp,      fwp,   rp,
     &        chiP2,  mytmp(npro),         sign_levelset(npro),
     &        sclr_ls(npro)
c    Kay-Epsilon
      real*8 Ry(npro), Rt(npro), RtP2(npro), f2(npro), f1(npro),
     &       kay(npro), epsilon(npro), fmu(npro), fmui(npro),
     &       srcjac(npro,4), srcRat(npro)
      integer e
c----------------------------------------------------------------------
c
c           --             --              --           --
c           |      cb2      |           1  |             |
c   phi,  + |u  - -----phi, | phi,  - -----|(nu+phi)phi, | 
c       t   | i   sigma    i|     i   sigma|            i|
c           --             --              --           --,
c                                                          i
c
c                  ~              phi^2
c            - cb1*S*phi + cw1*fw*-----
c                                  d^2
c
c----------------------------------------------------------------------
c
      if(iRANS.eq.-1) then    ! spalart almaras model
c
c.... compute the global gradient of u
c
         gradV = zero
         do n = 1, nshl
c
c          du_i/dx_j
c
c           i j   indices match array where V is the velocity (u in our notes)
            gradV(:,1,1) = gradV(:,1,1) + shg(:,n,1) * yl(:,n,2)
            gradV(:,2,1) = gradV(:,2,1) + shg(:,n,1) * yl(:,n,3)
            gradV(:,3,1) = gradV(:,3,1) + shg(:,n,1) * yl(:,n,4)
c
            gradV(:,1,2) = gradV(:,1,2) + shg(:,n,2) * yl(:,n,2)
            gradV(:,2,2) = gradV(:,2,2) + shg(:,n,2) * yl(:,n,3)
            gradV(:,3,2) = gradV(:,3,2) + shg(:,n,2) * yl(:,n,4)
c
            gradV(:,1,3) = gradV(:,1,3) + shg(:,n,3) * yl(:,n,2)
            gradV(:,2,3) = gradV(:,2,3) + shg(:,n,3) * yl(:,n,3)
            gradV(:,3,3) = gradV(:,3,3) + shg(:,n,3) * yl(:,n,4)
c                                             a j     u   a i
c from our notes where we had N_{a,j} = dN_a/dx_j  note that i is off by one because p was first in yl vector
c
         enddo
c
c.... magnitude of vorticity
c
         absVort  = sqrt( (gradV(:,2,3) - gradV(:,3,2)) ** 2
     &                  + (gradV(:,3,1) - gradV(:,1,3)) ** 2
     &		        + (gradV(:,1,2) - gradV(:,2,1)) ** 2 )

         dwall = zero
         do n = 1, nenl
            dwall     = dwall + shape_funct(:,n) * dwl(:,n)
         enddo

         sixth   = 1.0/6.0
c
c.... compute source and its jacobian
c
         do e = 1, npro
         
            if (dwall(e) .ne. 0) dP2Inv = 1.0 / dwall(e)**2
         
            k2d2Inv = dP2Inv * saKappaP2Inv
         
            viscInv = 1.0/datmat(1,2,1)
            chi     = abs( Sclr(e) * viscInv )
            chiP2   = chi * chi 
            chiP3   = chi * chiP2
            
            tmp     = 1.0 / ( chiP3 + saCv1P3 )
            fv1     = chiP3 * tmp
            fv1p    = 3.0 * viscInv * saCv1P3 * chiP2 * tmp**2
            
            tmp     = 1.0 / (1.0 + chi * fv1)
            fv2     = 1.0 - chi * tmp
            fv2p    = (chiP2 * fv1p - viscInv) * tmp**2
            
            s       = absVort(e)
            tmp     = Sclr(e) * k2d2Inv
            st      = s + fv2 * tmp
            stp     = k2d2Inv * (fv2 + fv2p * Sclr(e))
            
            r       = zero
            if (st .ne. 0 ) r = tmp / st
            r       = min( max(r, -8.0d0), 8.0d0)
            rp      = zero
            if (st .ne. 0 ) rp = k2d2Inv / st**2 * (st - Sclr(e)*stp)
            
            rP5     = r * (r * r) ** 2
            tmp     = 1.0 + saCw2 * (rP5 - 1.0)
            g       = r * tmp
            gp      = rp * (tmp + 5.0 * saCw2 * rP5)
            
            gP6     = (g * g * g) ** 2
            tmp     = 1.0 / (gP6 + saCw3P6)
            tmp1    = ( (1.0 + saCw3P6) * tmp ) ** sixth
            fw      = g * tmp1
            fwp     = gp * tmp1 * saCw3P6 * tmp
            
            tmp     = saCw1 * dP2Inv
            
            p       = saCb1 * st - fw * Sclr(e) * tmp 
            dp      = saCb1 * (st + stp * Sclr(e))
     &              - tmp * Sclr(e) * (2.0 * fw + fwp * Sclr(e))
            tmp     = zero

            if ( dp .le. 0 .and. dp - p .le. 0 ) then
               tmp = -dp
            else if ( p .lt. 0 ) then
               tmp = -p
            endif
            p       = p * Sclr(e) 
c
            srcL(e)  = tmp
            srcR(e)   = p
         enddo
c
c.... One source term has the form (beta_i)*(phi,_i).  Since
c     the convective term has (u_i)*(phi,_i), it is useful to treat
c     beta_i as a "correction" to the velocity.  In calculating the
c     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
c     then used in place of the pure velocity for stabilization terms,
c     and the source term sneaks into the RHS and LHS.  Here, the term
c     is given by beta_i=(cb_2)*(phi,_i)/(sigma)
c

         tmp = saCb2 * saSigmaInv
         uMod(:,1) = u1(:) - tmp * gradS(:,1)
         uMod(:,2) = u2(:) - tmp * gradS(:,2)
         uMod(:,3) = u3(:) - tmp * gradS(:,3)

      elseif(iRANS.eq.-2) then    ! K-epsilon
c.... compute the global gradient of u
c
         gradV = zero
         do n = 1, nshl
c
c          du_i/dx_j
c
c           i j   indices match array where V is the velocity (u in our notes)
            gradV(:,1,1) = gradV(:,1,1) + shg(:,n,1) * yl(:,n,2)
            gradV(:,2,1) = gradV(:,2,1) + shg(:,n,1) * yl(:,n,3)
            gradV(:,3,1) = gradV(:,3,1) + shg(:,n,1) * yl(:,n,4)
c
            gradV(:,1,2) = gradV(:,1,2) + shg(:,n,2) * yl(:,n,2)
            gradV(:,2,2) = gradV(:,2,2) + shg(:,n,2) * yl(:,n,3)
            gradV(:,3,2) = gradV(:,3,2) + shg(:,n,2) * yl(:,n,4)
c
            gradV(:,1,3) = gradV(:,1,3) + shg(:,n,3) * yl(:,n,2)
            gradV(:,2,3) = gradV(:,2,3) + shg(:,n,3) * yl(:,n,3)
            gradV(:,3,3) = gradV(:,3,3) + shg(:,n,3) * yl(:,n,4)
c                                             a j     u   a i
c from our notes where we had N_{a,j} = dN_a/dx_j  note that i is off by one because p was first in yl vector
c
         enddo

         dwall = zero
         do n = 1, nenl
            dwall     = dwall + shape_funct(:,n) * dwl(:,n)
         enddo

         kay(:)=zero
         epsilon(:)=zero
         do ii=1,npro
            do jj = 1, nshl
               kay(ii) =  kay(ii) + shape_funct(ii,jj) * yl(ii,jj,6)
               epsilon(ii) =  epsilon(ii) 
     &              + shape_funct(ii,jj) * yl(ii,jj,7)
            enddo
         enddo
c        no source term in the LHS yet
         srcL(:)=zero

         call elm3keps   (kay,     epsilon,
     &                    dwall,   gradV,
     &                    srcRat,  srcR,     srcJac)
         if(isclr.eq.1) srcL = srcJac(:,1)
         if(isclr.eq.2) srcL = srcJac(:,4)
         iadvdiff=0 ! scalar advection-diffusion flag
         if(iadvdiff.eq.1)then
            srcL(:)=zero
            srcR(:)=zero
            srcRat(:)=zero
         endif
c
c.... No source terms with the form (beta_i)*(phi,_i) for K or E
c
         uMod(:,1) = u1(:) - zero
         uMod(:,2) = u2(:) - zero
         uMod(:,3) = u3(:) - zero

      
      elseif (iLSet.ne.0) then
         if (isclr.eq.1)  then
            
            srcR = zero
            srcL = zero

         elseif (isclr.eq.2) then !we are redistancing level-sets

CAD   If Sclr(:,1).gt.zero, result of sign_term function 1
CAD   If Sclr(:,1).eq.zero, result of sign_term function 0
CAD   If Sclr(:,1).lt.zero, result of sign_term function -1

            sclr_ls = zero      !zero out temp variable
        
            do ii=1,npro

               do jj = 1, nshl  ! first find the value of levelset at point ii
                  
                  sclr_ls(ii) =  sclr_ls(ii) 
     &                        + shape_funct(ii,jj) * yl(ii,jj,6)

               enddo
               
               if (sclr_ls(ii) .gt. epsilon_ls)then
            
                  sign_levelset(ii) = one
                  
               elseif  (abs(sclr_ls(ii)) .le. epsilon_ls)then
                  
                  sign_levelset(ii) = sclr_ls(ii)/epsilon_ls 
     &                  + sin(pi*sclr_ls(ii)/epsilon_ls)/pi
                  
               elseif (sclr_ls(ii) .lt. epsilon_ls) then
                  
                  sign_levelset(ii) = - one
                  
               endif
               
               srcR(ii) = sign_levelset(ii)
               
            enddo
            
            srcL =  zero   

cad   The redistancing equation can be written in the following form
cad
cad   d_{,t} + sign(phi)*( d_{,i}/|d_{,i}| )* d_{,i} = sign(phi)
cad
cad   This is rewritten in the form
cad
cad   d_{,t} + u * d_{,i} = sign(phi)
cad

CAD   For the redistancing equation the "velocity" term is calculated as 
CAD   follows



            mytmp = sign_levelset / sqrt( gradS(:,1)*gradS(:,1) 
     &                            + gradS(:,2)*gradS(:,2) 
     &                            + gradS(:,3)*gradS(:,3)) 

            uMod(:,1) = mytmp * gradS(:,1) ! These are for the LHS
            uMod(:,2) = mytmp * gradS(:,2)
            uMod(:,3) = mytmp * gradS(:,3)

            u1 = umod(:,1)      ! These are for the RHS
            u2 = umod(:,2)
            u3 = umod(:,3)


         endif    ! close of scalar 2 of level set

      else        ! NOT turbulence and NOT level set so this is a simple
                  ! scalar. If your scalar equation has a source term
                  ! then add your own like the above but leave an unforced case
                  ! as an option like you see here

         srcR = zero
         srcL = zero
      endif

      return
      end



