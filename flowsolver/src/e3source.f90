      subroutine e3source(xx, src)
!-----------------------------------------------------------------------
!
!  this routine computes the body force term.
!
!  currently this computes a swirl body with the axis alligned with 
!  the z-coordinate
!
!-----------------------------------------------------------------------

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8   xx(npro,nsd), src(npro,nsd)

      real*8   nu

      real*8   r, Stheta, dpdz, rP5
      if(datmat(2,5,1).eq. 0) then  ! calc swirl 
!
!  This is the body force which will drive a swirl in a pipe flow
!     
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

!$$$         nu=datmat(1,2,1)
!$$$         omag=datmat(3,5,1) ! frame rotation rate
!$$$         e1 = one/sqrt(3.0d0) ! axis of rotation
!$$$         e2 = e1
!$$$         e3 = e1
!$$$         om1=omag*e1
!$$$         om2=omag*e2
!$$$         om3=omag*e3
!$$$         w=0
!$$$      
!$$$         do iel = 1, npro
!$$$         
!$$$            x = xx(iel,1)
!$$$            y = xx(iel,2)
!$$$            z = xx(iel,3)
!$$$
!$$$c
!$$$c  The following are MAPLE generated forcing functions for
!$$$c  a contrived flow in a rotating reference frame.
!$$$c
!$$$
!$$$      t1 = x**2
!$$$      t2 = t1**2
!$$$      t5 = dexp(14.D0*x)
!$$$      t6 = t2*x*t5
!$$$      t7 = y**2
!$$$      t8 = t7**2
!$$$      t9 = t8*t7
!$$$      t12 = t2*t5
!$$$      t15 = nu*t1
!$$$      t17 = dexp(7.D0*x)
!$$$      t18 = t7*y
!$$$      t19 = t17*t18
!$$$      t22 = t1*x
!$$$      t23 = nu*t22
!$$$      t24 = t17*y
!$$$      t29 = t17*t7
!$$$      t34 = nu*x
!$$$      t43 = nu*t2
!$$$      t46 = -40.D0*t6*t9-6.D0*t12*t7+92.D0*t15*t19+132.D0*t23*t24+80.D0*
!$$$     #t6*t18-252.D0*t23*t29+168.D0*t23*t19+96.D0*t34*t29-64.D0*t34*t19+2
!$$$     #2.D0*t15*t24-138.D0*t15*t29+294.D0*t43*t29
!$$$      t50 = t2*t22*t5
!$$$      t57 = nu*t17
!$$$      t60 = t8*y
!$$$      t65 = t2**2
!$$$      t66 = t65*t5
!$$$      t73 = t22*t5
!$$$      t77 = t2*t1*t5
!$$$      t80 = -196.D0*t43*t19-96.D0*t50*t9-122.D0*t43*t24-32.D0*t34*t24-4.
!$$$     #D0*t57*y+36.D0*t12*t60+192.D0*t50*t18+98.D0*t66*t8+24.D0*t12*t18+2
!$$$     #8.D0*t66*t9-24.D0*t73*t60-336.D0*t77*t60
!$$$      t106 = -336.D0*t50*t8-12.D0*t12*t9+8.D0*t73*t9-140.D0*t6*t8+28.D0*
!$$$     #t73*t8+12.D0*t15*t17+12.D0*t43*t17+288.D0*t50*t60+112.D0*t77*t9-22
!$$$     #4.D0*t77*t18-56.D0*t66*t18+120.D0*t6*t60
!$$$      t131 = -8.D0*t57*t18-20.D0*t6*t7+56.D0*t77*t7-42.D0*t12*t8-24.D0*t
!$$$     #23*t17-48.D0*t50*t7+4.D0*t73*t7-16.D0*t73*t18+392.D0*t77*t8+14.D0*
!$$$     #t66*t7-84.D0*t66*t60+12.D0*t57*t7
!$$$      fx = t46+t80+t106+t131
!$$$
!$$$
!$$$      t1 = x**2
!$$$      t2 = t1**2
!$$$      t3 = nu*t2
!$$$      t5 = dexp(7.D0*x)
!$$$      t6 = t5*y
!$$$      t9 = y**2
!$$$      t10 = t9**2
!$$$      t11 = nu*t10
!$$$      t18 = t1*x
!$$$      t22 = nu*x
!$$$      t25 = nu*t18
!$$$      t32 = dexp(14.D0*x)
!$$$      t33 = t2*t32
!$$$      t39 = t2*t1*t32
!$$$      t40 = t10*t9
!$$$      t43 = t1*t32
!$$$      t46 = t9*y
!$$$      t47 = t10*t46
!$$$      t50 = -84.D0*t3*t6+343.D0*t11*t5*t2+66.D0*t11*t5*x-98.D0*t11*t5*t1
!$$$     #8-24.D0*t22*t6+120.D0*t25*t6-287.D0*t11*t5*t1-140.D0*t33*t10+14.D0
!$$$     #*t3*t5-56.D0*t39*t40-28.D0*t43*t40+8.D0*t43*t47
!$$$      t52 = t2*x*t32
!$$$      t55 = t18*t32
!$$$      t62 = t10*y
!$$$      t69 = nu*t5
!$$$      t80 = 120.D0*t52*t10-32.D0*t55*t47+56.D0*t33*t47+30.D0*t11*t5-216.
!$$$     #D0*t52*t62-144.D0*t55*t62+72.D0*t39*t62-60.D0*t69*t46+30.D0*t69*t9
!$$$     #-20.D0*t25*t5-20.D0*t43*t10+36.D0*t43*t62
!$$$      t86 = nu*t1
!$$$      t89 = t5*t9
!$$$      t92 = t5*t46
!$$$      t109 = 112.D0*t55*t40+80.D0*t55*t10-12.D0*t86*t6-275.D0*t86*t89+57 
!$$$     #4.D0*t86*t92-218.D0*t25*t89+196.D0*t25*t92+90.D0*t22*t89-132.D0*t2
!$$$     #2*t92+427.D0*t3*t89+8.D0*t39*t46+4.D0*t43*t46
!$$$      t134 = 16.D0*t39*t47+28.D0*t33*t46-48.D0*t52*t47+252.D0*t33*t62+16
!$$$     #8.D0*t52*t40-24.D0*t52*t46-686.D0*t3*t92+2.D0*t86*t5-40.D0*t39*t10
!$$$     #-196.D0*t33*t40-16.D0*t55*t46+4.D0*t22*t5
!$$$      fy = t50+t80+t109+t134
!$$$
!$$$
!$$$      t3 = om3*x  
!$$$      t5 = dexp(7.D0*x)
!$$$      t7 = x**2
!$$$      t12 = y**2
!$$$      t15 = (-1.D0+y)**2
!$$$      fxa = 2.D0*om2*w+2.D0*t3*t5*(2.D0+x-10.D0*t7+7.D0*t7*x)*t12*t15+om
!$$$     #2*(om1*y-om2*x)-om3*(t3-om1*z)
!$$$
!$$$      t1 = x**2
!$$$      t4 = (-1.D0+x)**2
!$$$      t7 = dexp(7.D0*x)
!$$$      t10 = y**2
!$$$      fya = 4.D0*om3*t1*t4*t7*y*(1.D0-3.D0*y+2.D0*t10)-2.D0*om1*w+om3*(o
!$$$     #m2*z-om3*y)-om1*(om1*y-om2*x)
!$$$
!$$$
!$$$      t3 = dexp(7.D0*x)
!$$$      t5 = x**2
!$$$      t10 = y**2
!$$$      t13 = (-1.D0+y)**2
!$$$      t19 = (-1.D0+x)**2
!$$$      fza = -2.D0*om1*x*t3*(2.D0+x-10.D0*t5+7.D0*t5*x)*t10*t13-4.D0*om2*
!$$$     #t5*t19*t3*y*(1.D0-3.D0*y+2.D0*t10)+om1*(om3*x-om1*z)-om2*(om2*z-om
!$$$     #3*y)
!$$$
!$$$      src(iel,1) =  fx + fxa
!$$$      src(iel,2) =  fy + fya
!$$$      src(iel,3) =       fza
!$$$      enddo
!Analytic lid case
   
      do iel = 1, npro
            x = xx(iel,1)
            y = xx(iel,2)
            z = xx(iel,3)
            Re = 1.0/datmat(1,2,1)
!
!  The following are MAPLE generated forcing functions for
!  a Lid Driven cavity flow with an analytic solution
!
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
            t54 = -8/Re*(24.E0/5.E0*t4-12*t3+t8+2*(4*t7-6*t2+2*x)*(12 &
                 *t13-2)+(24*x-12)*t21)-64*(t26/2-2*t3*t7+3*t30-2*t4+t3 &
                 /2)*(-24*t20*y+8*t37-4*y)+64*t21*(4*t37-2*y)*(-4*t30+12 &
                 *t4-14*t3+t8-2*t2)

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



      subroutine e3sourceSclr ( Sclr,      Sdot,      gradS, &
                                dwl,       shape_funct,     shg,    &
                                yl,        dxidx,     rmu, &
                                u1,        u2,        u3, &
                                srcR,      srcL,      uMod, &
                                srcRat) 


!-----------------------------------------------------------------------
!
!     calculate the coefficients for the Spalart-Allmaras 
!     turbulence model and its jacobian
!
!     input: 
!             Sclr  :   The scalar that is being transported in this pde.
!             gradS :   spatial gradient of Sclr
!             rmu   :   kinematic viscosity
!
!     output:
!             rmu   :   diffusion for eddy viscosity equation
!             srcR  :   residual terms for  
!                         (srcR * turb)
!             srcL :   jacobian of srcR
!
!-----------------------------------------------------------------------
      use turbSA
      use turbKE
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
! coming in      
      real*8  Sclr(npro),          Sdot(npro), &
              gradS(npro,nsd),     dwl(npro,nenl), &
              shape_funct(npro,nshl),    shg(npro,nshl,nsd), &
              yl(npro,nenl,ndof),  dxidx(npro,nsd,nsd), &
              rmu(npro),           u1(npro), &
              u2(npro),            u3(npro)
! going out
      real*8  srcR(npro),          srcL(npro), &
              uMod(npro,nsd)
! used locally

      real*8  gradv(npro,nsd,nsd), absVort(npro), &
              dwall(npro)
      real*8  chi,    chiP3,   fv1,   fv2,     st,    r, &
              g,      fw,      s,     viscInv, k2d2Inv, &
              dP2Inv, sixth,   tmp,   tmp1,    p,     dp, &
              fv1p,   fv2p,    stp,         gp,      fwp,   rp, &
              chiP2,  mytmp(npro),         sign_levelset(npro), &
              sclr_ls(npro)
!    Kay-Epsilon
      real*8 Ry(npro), Rt(npro), RtP2(npro), f2(npro), f1(npro), &
             kay(npro), epsilon(npro), fmu(npro), fmui(npro), &
             srcjac(npro,4), srcRat(npro)
      integer e
!----------------------------------------------------------------------
!
!           --             --              --           --
!           |      cb2      |           1  |             |
!   phi,  + |u  - -----phi, | phi,  - -----|(nu+phi)phi, | 
!       t   | i   sigma    i|     i   sigma|            i|
!           --             --              --           --,
!                                                          i
!
!                  ~              phi^2
!            - cb1*S*phi + cw1*fw*-----
!                                  d^2
!
!----------------------------------------------------------------------
!
      if(iRANS.eq.-1) then    ! spalart almaras model
!
!.... compute the global gradient of u
!
         gradV = zero
         do n = 1, nshl
!
!          du_i/dx_j
!
!           i j   indices match array where V is the velocity (u in our notes)
            gradV(:,1,1) = gradV(:,1,1) + shg(:,n,1) * yl(:,n,2)
            gradV(:,2,1) = gradV(:,2,1) + shg(:,n,1) * yl(:,n,3)
            gradV(:,3,1) = gradV(:,3,1) + shg(:,n,1) * yl(:,n,4)
!
            gradV(:,1,2) = gradV(:,1,2) + shg(:,n,2) * yl(:,n,2)
            gradV(:,2,2) = gradV(:,2,2) + shg(:,n,2) * yl(:,n,3)
            gradV(:,3,2) = gradV(:,3,2) + shg(:,n,2) * yl(:,n,4)
!
            gradV(:,1,3) = gradV(:,1,3) + shg(:,n,3) * yl(:,n,2)
            gradV(:,2,3) = gradV(:,2,3) + shg(:,n,3) * yl(:,n,3)
            gradV(:,3,3) = gradV(:,3,3) + shg(:,n,3) * yl(:,n,4)
!                                             a j     u   a i
! from our notes where we had N_{a,j} = dN_a/dx_j  note that i is off by one because p was first in yl vector
!
         enddo
!
!.... magnitude of vorticity
!
         absVort  = sqrt( (gradV(:,2,3) - gradV(:,3,2)) ** 2 &
                        + (gradV(:,3,1) - gradV(:,1,3)) ** 2 &
      		        + (gradV(:,1,2) - gradV(:,2,1)) ** 2 )

         dwall = zero
         do n = 1, nenl
            dwall     = dwall + shape_funct(:,n) * dwl(:,n)
         enddo

         sixth   = 1.0/6.0
!
!.... compute source and its jacobian
!
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
            dp      = saCb1 * (st + stp * Sclr(e)) &
                    - tmp * Sclr(e) * (2.0 * fw + fwp * Sclr(e))
            tmp     = zero

            if ( dp .le. 0 .and. dp - p .le. 0 ) then
               tmp = -dp
            else if ( p .lt. 0 ) then
               tmp = -p
            endif
            p       = p * Sclr(e) 
!
            srcL(e)  = tmp
            srcR(e)   = p
         enddo
!
!.... One source term has the form (beta_i)*(phi,_i).  Since
!     the convective term has (u_i)*(phi,_i), it is useful to treat
!     beta_i as a "correction" to the velocity.  In calculating the
!     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
!     then used in place of the pure velocity for stabilization terms,
!     and the source term sneaks into the RHS and LHS.  Here, the term
!     is given by beta_i=(cb_2)*(phi,_i)/(sigma)
!

         tmp = saCb2 * saSigmaInv
         uMod(:,1) = u1(:) - tmp * gradS(:,1)
         uMod(:,2) = u2(:) - tmp * gradS(:,2)
         uMod(:,3) = u3(:) - tmp * gradS(:,3)

      elseif(iRANS.eq.-2) then    ! K-epsilon
!.... compute the global gradient of u
!
         gradV = zero
         do n = 1, nshl
!
!          du_i/dx_j
!
!           i j   indices match array where V is the velocity (u in our notes)
            gradV(:,1,1) = gradV(:,1,1) + shg(:,n,1) * yl(:,n,2)
            gradV(:,2,1) = gradV(:,2,1) + shg(:,n,1) * yl(:,n,3)
            gradV(:,3,1) = gradV(:,3,1) + shg(:,n,1) * yl(:,n,4)
!
            gradV(:,1,2) = gradV(:,1,2) + shg(:,n,2) * yl(:,n,2)
            gradV(:,2,2) = gradV(:,2,2) + shg(:,n,2) * yl(:,n,3)
            gradV(:,3,2) = gradV(:,3,2) + shg(:,n,2) * yl(:,n,4)
!
            gradV(:,1,3) = gradV(:,1,3) + shg(:,n,3) * yl(:,n,2)
            gradV(:,2,3) = gradV(:,2,3) + shg(:,n,3) * yl(:,n,3)
            gradV(:,3,3) = gradV(:,3,3) + shg(:,n,3) * yl(:,n,4)
!                                             a j     u   a i
! from our notes where we had N_{a,j} = dN_a/dx_j  note that i is off by one because p was first in yl vector
!
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
               epsilon(ii) =  epsilon(ii)  &
                    + shape_funct(ii,jj) * yl(ii,jj,7)
            enddo
         enddo
!        no source term in the LHS yet
         srcL(:)=zero

         call elm3keps   (kay,     epsilon, &
                          dwall,   gradV, &
                          srcRat,  srcR,     srcJac)
         if(isclr.eq.1) srcL = srcJac(:,1)
         if(isclr.eq.2) srcL = srcJac(:,4)
         iadvdiff=0 ! scalar advection-diffusion flag
         if(iadvdiff.eq.1)then
            srcL(:)=zero
            srcR(:)=zero
            srcRat(:)=zero
         endif
!
!.... No source terms with the form (beta_i)*(phi,_i) for K or E
!
         uMod(:,1) = u1(:) - zero
         uMod(:,2) = u2(:) - zero
         uMod(:,3) = u3(:) - zero

      
      elseif (iLSet.ne.0) then
         if (isclr.eq.1)  then
            
            srcR = zero
            srcL = zero

         elseif (isclr.eq.2) then !we are redistancing level-sets

!AD   If Sclr(:,1).gt.zero, result of sign_term function 1
!AD   If Sclr(:,1).eq.zero, result of sign_term function 0
!AD   If Sclr(:,1).lt.zero, result of sign_term function -1

            sclr_ls = zero      !zero out temp variable
        
            do ii=1,npro

               do jj = 1, nshl  ! first find the value of levelset at point ii
                  
                  sclr_ls(ii) =  sclr_ls(ii)  &
                              + shape_funct(ii,jj) * yl(ii,jj,6)

               enddo
               
               if (sclr_ls(ii) .gt. epsilon_ls)then
            
                  sign_levelset(ii) = one
                  
               elseif  (abs(sclr_ls(ii)) .le. epsilon_ls)then
                  
                  sign_levelset(ii) = sclr_ls(ii)/epsilon_ls  &
                        + sin(pi*sclr_ls(ii)/epsilon_ls)/pi
                  
               elseif (sclr_ls(ii) .lt. epsilon_ls) then
                  
                  sign_levelset(ii) = - one
                  
               endif
               
               srcR(ii) = sign_levelset(ii)
               
            enddo
            
            srcL =  zero   

!ad   The redistancing equation can be written in the following form
!ad
!ad   d_{,t} + sign(phi)*( d_{,i}/|d_{,i}| )* d_{,i} = sign(phi)
!ad
!ad   This is rewritten in the form
!ad
!ad   d_{,t} + u * d_{,i} = sign(phi)
!ad

!AD   For the redistancing equation the "velocity" term is calculated as 
!AD   follows



            mytmp = sign_levelset / sqrt( gradS(:,1)*gradS(:,1)  &
                                  + gradS(:,2)*gradS(:,2)  &
                                  + gradS(:,3)*gradS(:,3)) 

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



