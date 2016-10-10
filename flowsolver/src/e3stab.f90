      subroutine e3stab (rho,          u1,       u2, &
                         u3,           &
                         uMesh1, uMesh2, uMesh3, & !ALE variables added MAF 08/10/2016
                         dxidx,    rLui,    &
                         rmu,          tauC,     tauM,    &
                         tauBar,       uBar )  
!
!----------------------------------------------------------------------
!
! This routine computes the diagonal Tau for least-squares operator.  
! Diagonal tau proposed by Shakib.
!
! input:
!  u1     (npro)           : x1-velocity component
!  u2     (npro)           : x2-velocity component
!  u3     (npro)           : x3-velocity component
!  dxidx  (npro,nsd,nsd)   : inverse of deformation gradient
!  rLui   (npro,nsd)      : least-squares residual vector
!
! output:
!  tauC    (npro)          : continuity tau
!  tauM    (npro)          : momentum tau
!  tauBar  (npro)          : additional tau
!  uBar    (npro,nsd)      : modified velocity
!  cfl_loc(npro) 	   : CFL of the element
!
! Zdenek Johan, Summer 1990.  (Modified from e2tau.f)
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        use phcommonvars
        use ale
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension rho(npro),                 u1(npro), &
                  u2(npro),                  u3(npro), &
                  dxidx(npro,nsd,nsd),  &
                  rLui(npro,nsd), &
                  tauC(npro),    tauM(npro), tauBar(npro), &
                  rmu(npro),     uBar(npro,3), unorm(npro)

!
        dimension gijd(npro,6),       fact(npro), rnu(npro), &
             rhoinv(npro)

        !..... ALE variables

      dimension uMesh1(npro), uMesh2(npro), uMesh3(npro)
      integer   i
      logical :: exist

      
      ! if (rigidOn.eq.1) then
      ! write(*,*) "rigidOn"
      ! ! uMesh1(:) = 0.0d0 !no relative stabilization for rigid body motion
      ! ! uMesh2(:) = 0.0d0
      ! ! uMesh3(:) = 0.0d0
      ! uMesh1(:) = -1.0d0*globalRigidVelocity(1)
      ! uMesh2(:) = -1.0d0*globalRigidVelocity(2)
      ! uMesh3(:) = -1.0d0*globalRigidVelocity(3)
      ! else

      !     get mesh velocity KDL, MAF
      ! if (aleRigid.eq.1) then
      !   uMesh1(:) = globalRigidVelocity(1)
      !   uMesh2(:) = globalRigidVelocity(2)
      !   uMesh3(:) = globalRigidVelocity(3)
      ! else
      !   uMesh1(:) = real(0.0,8) 
      !   uMesh2(:) = real(0.0,8)
      !   uMesh3(:) = real(0.0,8)        
      ! endif
      ! endif
      
      ! MAYBE PRECALCULATE OUTSIDE ?
      !if (iALE .eq. 1) then
      !  call calculateMeshVelocity(uMesh1,uMesh2,uMesh3)
      !end if     
!
!
!.... get the metric tensor
!      
      call e3gijd( dxidx, gijd )
!
!... higher order element diffusive correction
!
      if (ipord == 1) then
         fff = 36.0d0
      else if (ipord == 2) then
         fff = 60.0d0
!     fff = 36.0d0
      else if (ipord == 3) then
         fff = 128.0d0
!     fff = 144.0d0
      endif

      omegasq=zero
      if(matflg(6,1).eq.1) omegasq = datmat(1,6,1)**2 &
                                    +datmat(2,6,1)**2 &
                                    +datmat(3,6,1)**2
      rhoinv=one/rho
      rnu=rmu*rhoinv

      if(itau.eq.0)  then  ! original tau
!
!...  momentum tau
! 
         dts=  Dtgl*dtsfct	! Dtgl = (time step)^-1, dtsfct = c1
         tauM = ( (two*dts)**2 &
      		      + ( (u1-uMesh1) * ( gijd(:,1) * (u1-uMesh1) &
      			                      + gijd(:,4) * (u2-uMesh2) &
      			                      + gijd(:,6) * (u3-uMesh3) ) &
      		        + (u2-uMesh2) * ( gijd(:,4) * (u1-uMesh1) &
      			                      + gijd(:,2) * (u2-uMesh2) &
      			                      + gijd(:,5) * (u3-uMesh3) ) &
      		        + (u3-uMesh3) * ( gijd(:,6) * (u1-uMesh1) &
      			                      + gijd(:,5) * (u2-uMesh2) &
      			                      + gijd(:,3) * (u3-uMesh3) ) ) ) &
      		    + fff * rnu** 2 &
      		    * ( gijd(:,1) ** 2 &
      		      + gijd(:,2) ** 2 &
      		      + gijd(:,3) ** 2 &
      		      + 2. &
      		      * ( gijd(:,4) ** 2 &
      		        + gijd(:,5) ** 2 &
      		        + gijd(:,6) ** 2 ) &
                    +omegasq)
        
         fact = sqrt(tauM)
         dtsi=one/dts
         ff=taucfct/dtsfct
         tauC =rho* pt125*fact/(gijd(:,1)+gijd(:,2)+gijd(:,3))*ff
         tauM = one/fact

#if DEBUG_ALE == 1      

      write(*,*) 'printing inside e3stab'    
      inquire(file="taum.dat", exist=exist)
      if (exist) then
        open(794, file="taum.dat", status="old", position="append", action="write")
      else
        open(794, file="taum.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(794,'(1(e40.20))') (tauM(i))      
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      end do 
      close(794)
#endif

      else if(itau.eq.1)  then  ! new tau



!
!  determinant of gijd
!
         fact = gijd(:,1) * gijd(:,2) * gijd(:,3) &
              - gijd(:,2) * gijd(:,6) * gijd(:,6) &
              - gijd(:,1) * gijd(:,5) * gijd(:,5) &
              - gijd(:,3) * gijd(:,4) * gijd(:,4) &
              + gijd(:,6) * gijd(:,4) * gijd(:,5) * two
        
!
! put 1/2u*h1 = sqrt(u_i g^{ij} u_j) into tau_M  note inverse is calculated
! on the fly here from cofactors over the determinent dotted from left and 
! right with u
!
         
         tauM = &
             u1 * ( (gijd(:,2)*gijd(:,3)-gijd(:,5)*gijd(:,5))  * u1 &
           +  two * (gijd(:,5)*gijd(:,6)-gijd(:,4)*gijd(:,3))  * u2 &
           +  two * (gijd(:,4)*gijd(:,5)-gijd(:,6)*gijd(:,2))  * u3) &
           + u2 * ( (gijd(:,1)*gijd(:,3)-gijd(:,6)*gijd(:,6))  * u2 &
           +  two * (gijd(:,4)*gijd(:,6)-gijd(:,1)*gijd(:,5))  * u3) &
           + u3 * ( (gijd(:,1)*gijd(:,2)-gijd(:,4)*gijd(:,4))  * u3)
         tauM=fact/taum  ! here we have (u_i g^{ij} u^j)^{-1} approx 4/u^2h^2
!
!  we can calculate tauC more efficiently now
!
         tauC=tauM*(one+tauM*rmu*rmu)
         tauC=one/tauC
         tauC=taucfct*sqrt(tauC)
!
!
!...  momentum tau
!
!
!     this tau needs a u/h instead of a u*h so we contract with g_{ij} as
!     follows  (i.e. u_i g_{ij} u_j approx u^2/(h^2)/4) 
!
!     the definition of gijd appears to be as follows ... KDL, MAF
!     | g11 g12 g13 |    | gijd(1)  gijd(4) gijd(6) |
!     | ... g22 g23 |  = | .......  gijd(2) gijd(5) |
!     | ... ... g33 |    | .......  ....... gijd(3) |
!

         fact = &
                u1 * ( gijd(:,1) * u1 &      
                     + gijd(:,4) * u2 &
                     + gijd(:,6) * u3 ) &
              + u2 * ( gijd(:,4) * u1 &
                     + gijd(:,2) * u2 &
                     + gijd(:,5) * u3 ) &
              + u3 * ( gijd(:,6) * u1 &
                     + gijd(:,5) * u2 &
                     + gijd(:,3) * u3 )  
! 
! first limit dt effect on tau from causing trouble if user drops CFL below
! .05 (this could cause loss of spatial stability)
!
         velsq=vel*vel
         unorm = (u1*u1+u2*u2+u3*u3)/velsq
         dtsfsq=dtsfct*dtsfct
         dt=one/Dtgl
         taubar=  dtsfsq/( dt*dt + .01*unorm/fact)  ! never gets above (C_1 20*u_inf/h)^2
!
!  this means tau will never get below h/(20*C_1*u) no matter what time step 
!  you choose.  The 0.01 constant comes from minCFL=.05=> .05*.05*4 (where the 
!  4 comes from the bi-unit mapping). If you want to limit sooner the formula
!  would be  ".01-factor"=minCFL^2*4
!

         tauM = rho ** 2 &
      		    * ( four*taubar + fact &
      		    + fff * rmu** 2 &
      		    * ( gijd(:,1) ** 2 &
      		      + gijd(:,2) ** 2 &
      		      + gijd(:,3) ** 2 &
      		      + 2. &
      		      * ( gijd(:,4) ** 2 &
      		        + gijd(:,5) ** 2 &
      		        + gijd(:,6) ** 2 ) )  &
                    +omegasq)
         fact=sqrt(tauM)
!debugcheck         tauBar = pt125*fact/(gijd(:,1)+gijd(:,2)+gijd(:,3)) !*dtsi
      
        tauM=one/fact           ! turn it right side up.
      else if(itau.eq.2)  then  ! new tau different continuity h

         unorm = (u1*u1+u2*u2+u3*u3)
         
         tauM=(gijd(:,1)+gijd(:,2)+gijd(:,3))/unorm ! here we have  4/u^2h^2
!
!  we can calculate tauC more efficiently now
!
         tauC=tauM*(one+tauM*rmu*rmu)
         tauC=one/tauC
         tauC=sqrt(tauC)*taucfct
!
!
!...  momentum tau
!
!
!     this tau needs a u/h instead of a u*h so we contract with g_{ij} as
!     follows  (i.e. u_i g_{ij} u_j approx u^2/(h^2)/4) 
!
         fact = &
                u1 * ( gijd(:,1) * u1 &
                     + gijd(:,4) * u2 &
                     + gijd(:,6) * u3 ) &
              + u2 * ( gijd(:,4) * u1 &
                     + gijd(:,2) * u2 &
                     + gijd(:,5) * u3 ) &
              + u3 * ( gijd(:,6) * u1 &
                     + gijd(:,5) * u2 &
                     + gijd(:,3) * u3 ) 
! 
! first limit dt effect on tau from causing trouble if user drops CFL below
! .05 (this could cause loss of spatial stability)
!
         velsq=vel*vel
         dtsfsq=dtsfct*dtsfct
         dt=one/Dtgl
         unorm=unorm/velsq
         taubar=  dtsfsq/( dt*dt + .01*unorm/fact)  ! never gets above (C_1 20*u_inf/h)^2
!
!  this means tau will never get below h/(20*C_1*u) no matter what time step 
!  you choose.  The 0.01 constant comes from minCFL=.05=> .05*.05*4 (where the 
!  4 comes from the bi-unit mapping). If you want to limit sooner the formula
!  would be  ".01-factor"=minCFL^2*4
!

         tauM = rho ** 2 &
      		    * ( four*taubar + fact &
      		    + fff * rmu** 2 &
      		    * ( gijd(:,1) ** 2 &
      		      + gijd(:,2) ** 2 &
      		      + gijd(:,3) ** 2 &
      		      + 2. &
      		      * ( gijd(:,4) ** 2 &
      		        + gijd(:,5) ** 2 &
      		        + gijd(:,6) ** 2 ) ) &
                    +omegasq)
         fact=sqrt(tauM)
!         tauBar = pt125*fact/(gijd(:,1)+gijd(:,2)+gijd(:,3)) !*dtsi
      
        tauM=one/fact           ! turn it right side up.
      else if(itau.eq.3)  then  ! compressible tau

!
!  determinant of gijd
!
         fact = gijd(:,1) * gijd(:,2) * gijd(:,3) &
              - gijd(:,2) * gijd(:,6) * gijd(:,6) &
              - gijd(:,1) * gijd(:,5) * gijd(:,5) &
              - gijd(:,3) * gijd(:,4) * gijd(:,4) &
              + gijd(:,6) * gijd(:,4) * gijd(:,5) * two
        
!
! put 1/2u*h1 = sqrt(u_i g^{ij} u_j) into tau_M  note inverse is calculated
! on the fly here from cofactors over the determinent dotted from left and 
! right with u
!
         
         tauM = &
             u1 * ( (gijd(:,2)*gijd(:,3)-gijd(:,5)*gijd(:,5))  * u1 &
           +  two * (gijd(:,5)*gijd(:,6)-gijd(:,4)*gijd(:,3))  * u2 &
           +  two * (gijd(:,4)*gijd(:,5)-gijd(:,6)*gijd(:,2))  * u3) &
           + u2 * ( (gijd(:,1)*gijd(:,3)-gijd(:,6)*gijd(:,6))  * u2 &
           +  two * (gijd(:,4)*gijd(:,6)-gijd(:,1)*gijd(:,5))  * u3) &
           + u3 * ( (gijd(:,1)*gijd(:,2)-gijd(:,4)*gijd(:,4))  * u3) 
!
!  we can calculate tauC more efficiently now
!
         tauM=sqrt(tauM/fact)*two
         tauC=pt5*tauM*min(one,pt5*tauM/rmu)*taucfct
!
!
!...  momentum tau
!
!
!     this tau needs a u/h instead of a u*h so we contract with g_{ij} as
!     follows  (i.e. u_i g_{ij} u_j approx u^2/(h^2)/4) 
!
         fact = &
                u1 * ( gijd(:,1) * u1 &
                     + gijd(:,4) * u2 &
                     + gijd(:,6) * u3 ) &
              + u2 * ( gijd(:,4) * u1 &
                     + gijd(:,2) * u2 &
                     + gijd(:,5) * u3 ) &
              + u3 * ( gijd(:,6) * u1 &
                     + gijd(:,5) * u2 &
                     + gijd(:,3) * u3 )  
         fact=one/sqrt(fact)

         unorm = (u1*u1+u2*u2+u3*u3)

         dts= one/( Dtgl*dtsfct)
         tauM =min(dts,min(fact,fact*fact*unorm*pt33/rmu))
      endif
!
!.... calculate tauBar
!
      tauBar = rLui(:,1) * ( gijd(:,1) * rLui(:,1) &
                             + gijd(:,4) * rLui(:,2) &
                             + gijd(:,6) * rLui(:,3) ) &
               + rLui(:,2) * ( gijd(:,4) * rLui(:,1) &
                             + gijd(:,2) * rLui(:,2) &
                             + gijd(:,5) * rLui(:,3) )  &
               + rLui(:,3) * ( gijd(:,6) * rLui(:,1) &
                             + gijd(:,5) * rLui(:,2) &
                             + gijd(:,3) * rLui(:,3) )
      where ( tauBar .ne. 0.0 ) 
         tauBar = tauM / sqrt(tauBar)
      endwhere

!
!.... compute the modified velocity, uBar
!
        uBar(:,1) = u1 - tauM * rLui(:,1)*rhoinv
        uBar(:,2) = u2 - tauM * rLui(:,2)*rhoinv
        uBar(:,3) = u3 - tauM * rLui(:,3)*rhoinv
!     
!.... return
!
        return
        end

!-----------------------------------------------------------------------
!
!  Momentum tau
!
!-----------------------------------------------------------------------
      subroutine e3uBar (rho,          ui,         dxidx,      &
                         rLui,         rmu,        uBar )         

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      real*8     rho(npro),            ui(npro,nsd), &
                 dxidx(npro,nsd,nsd),  rLui(npro,nsd), &
                 rmu(npro),            uBar(npro,nsd)

      real*8     gijd(npro,6),         tauM(npro)

!
!.... get the metric tensor
!      
      call e3gijd( dxidx, gijd )
!
!.... higher order element diffusive correction
!
      if (ipord == 1) then
         fff = 36.0d0
      else if (ipord == 2) then
         fff = 60.0d0
      else if (ipord == 3) then
         fff = 128.0d0
      endif

      dts  =  (Dtgl*dtsfct)
      tauM = rho ** 2 &
      		    * ( (two*dts)**2 &
      		      + ( ui(:,1) * ( gijd(:,1) * ui(:,1) &
      			            + gijd(:,4) * ui(:,2) &
      			            + gijd(:,6) * ui(:,3) ) &
      		        + ui(:,2) * ( gijd(:,4) * ui(:,1) &
      			            + gijd(:,2) * ui(:,2) &
     			            + gijd(:,5) * ui(:,3) ) &
      		        + ui(:,3) * ( gijd(:,6) * ui(:,1) &
      			            + gijd(:,5) * ui(:,2) &
      			            + gijd(:,3) * ui(:,3) ) ) ) &
      		    + fff * rmu** 2 &
      		    * ( gijd(:,1) ** 2 &
      		      + gijd(:,2) ** 2 &
      		      + gijd(:,3) ** 2 &
      		      + 2. &
      		      * ( gijd(:,4) ** 2 &
      		        + gijd(:,5) ** 2 &
      		        + gijd(:,6) ** 2 ) )
        
      tauM = one/sqrt(tauM)
!
!.... compute the modified velocity, uBar
!
      uBar(:,1) = ui(:,1) - tauM * rLui(:,1)
      uBar(:,2) = ui(:,2) - tauM * rLui(:,2)
      uBar(:,3) = ui(:,3) - tauM * rLui(:,3)

      return
      end

!-----------------------------------------------------------------------
! get the metric tensor g_{ij}=xi_{k,i} xi_{k,j}.  
!-----------------------------------------------------------------------
      subroutine e3gijd( dxidx,  gijd )
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8  dxidx(npro,nsd,nsd),  gijd(npro,6), &
              tmp1(npro),           tmp2(npro), &
              tmp3(npro)
!
!  form metric tensor g_{ij}=xi_{k,i} xi_{k,j}.  It is a symmetric
!  tensor so we only form 6 components and use symmetric matrix numbering.
!
      if (lcsyst .ge. 2) then  ! note this makes wedges like hexs..should
!                                be corrected later

         gijd(:,1) = dxidx(:,1,1) * dxidx(:,1,1) &
                   + dxidx(:,2,1) * dxidx(:,2,1) &
                   + dxidx(:,3,1) * dxidx(:,3,1)
!
         gijd(:,4) = dxidx(:,1,1) * dxidx(:,1,2) &
                   + dxidx(:,2,1) * dxidx(:,2,2) &
                   + dxidx(:,3,1) * dxidx(:,3,2)
!
         gijd(:,2) = dxidx(:,1,2) * dxidx(:,1,2) &
                   + dxidx(:,2,2) * dxidx(:,2,2) &
                   + dxidx(:,3,2) * dxidx(:,3,2)
!
         gijd(:,5) = dxidx(:,1,2) * dxidx(:,1,3) &
                   + dxidx(:,2,2) * dxidx(:,2,3) &
                   + dxidx(:,3,2) * dxidx(:,3,3)
!
         gijd(:,6) = dxidx(:,1,1) * dxidx(:,1,3) &
                   + dxidx(:,2,1) * dxidx(:,2,3) &
                   + dxidx(:,3,1) * dxidx(:,3,3)
!
         gijd(:,3) = dxidx(:,1,3) * dxidx(:,1,3) &
                   + dxidx(:,2,3) * dxidx(:,2,3) &
                   + dxidx(:,3,3) * dxidx(:,3,3)
!
      else   if (lcsyst .eq. 1) then
!
!  There is an invariance problem with tets 
!  It is fixed by the following modifications to gijd 
!

         c1 = 1.259921049894873D+00
         c2 = 6.299605249474365D-01
!
         tmp1(:) = c1 * dxidx(:,1,1)+c2 *(dxidx(:,2,1)+dxidx(:,3,1))
         tmp2(:) = c1 * dxidx(:,2,1)+c2 *(dxidx(:,1,1)+dxidx(:,3,1))
         tmp3(:) = c1 * dxidx(:,3,1)+c2 *(dxidx(:,1,1)+dxidx(:,2,1))
         gijd(:,1) = dxidx(:,1,1) * tmp1 &
                    + dxidx(:,2,1) * tmp2 &
                    + dxidx(:,3,1) * tmp3
!
         tmp1(:) = c1 * dxidx(:,1,2)+c2 *(dxidx(:,2,2)+dxidx(:,3,2))
         tmp2(:) = c1 * dxidx(:,2,2)+c2 *(dxidx(:,1,2)+dxidx(:,3,2))
         tmp3(:) = c1 * dxidx(:,3,2)+c2 *(dxidx(:,1,2)+dxidx(:,2,2))
         gijd(:,2) = dxidx(:,1,2) * tmp1 &
                   + dxidx(:,2,2) * tmp2 &
                   + dxidx(:,3,2) * tmp3
!
         gijd(:,4) = dxidx(:,1,1) * tmp1 &
                   + dxidx(:,2,1) * tmp2 &
                   + dxidx(:,3,1) * tmp3
!
         tmp1(:) = c1 * dxidx(:,1,3)+c2 *(dxidx(:,2,3)+dxidx(:,3,3))
         tmp2(:) = c1 * dxidx(:,2,3)+c2 *(dxidx(:,1,3)+dxidx(:,3,3))
         tmp3(:) = c1 * dxidx(:,3,3)+c2 *(dxidx(:,1,3)+dxidx(:,2,3))
         gijd(:,3) = dxidx(:,1,3) * tmp1 &
                   + dxidx(:,2,3) * tmp2 &
                   + dxidx(:,3,3) * tmp3
!
         gijd(:,5) = dxidx(:,1,2) * tmp1 &
                   + dxidx(:,2,2) * tmp2 &
                   + dxidx(:,3,2) * tmp3
!
         gijd(:,6) = dxidx(:,1,1) * tmp1 &
                   + dxidx(:,2,1) * tmp2 &
                   + dxidx(:,3,1) * tmp3
!
      else
         write(*,*) 'lcsyst eq',lcsyst,'not supported'
         stop
      endif

      return
      end

!------------------------------------------------------------------------
!
!     calculate the stabilization for the advection-diffusion equation
!
!------------------------------------------------------------------------
      subroutine e3StabSclr (uMod,  dxidx,  tauT, diffus, srcP, giju, &
                             srcRat )
!
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        real*8    rho(npro),                 uMod(npro,nsd), &
                  dxidx(npro,nsd,nsd),       diffus(npro), &
                  tauT(npro),                srcP(npro)

!
        real*8    gijd(npro,6),       giju(npro,6),    &
                  tmp1(npro),         tmp2(npro), &
                  tmp3(npro),         fact(npro), &
                  srcRat(npro)

        real*8     fff
        if(ivart.eq.1) then
           tauT=zero
           return
        endif
!
!.... get the metric tensor
!      
      call e3gijd( dxidx, gijd )
!
!...  momentum tau
! 
!
!... higher order element diffusive correction
!
        if (ipord == 1) then
           fff = 9.0d0
        else if (ipord == 2) then
           fff = 36.0d0
        else if (ipord == 3) then
           fff = 64.0d0
        endif

        dts=  (Dtgl*dtsfct)
        !if(iRANS.ne.-2) srcRat=srcP
        tauT = &
      	       (two*dts)**2 &
             + srcRat ** 2 &
      	     + uMod(:,1) * ( gijd(:,1) * uMod(:,1) &
      	                   + gijd(:,4) * uMod(:,2) &
      	                   + gijd(:,6) * uMod(:,3) ) &
      	     + uMod(:,2) * ( gijd(:,4) * uMod(:,1) &
      	                   + gijd(:,2) * uMod(:,2) &
      	                   + gijd(:,5) * uMod(:,3) ) &
      	     + uMod(:,3) * ( gijd(:,6) * uMod(:,1) &
      	                   + gijd(:,5) * uMod(:,2) &
      	                   + gijd(:,3) * uMod(:,3) ) &
      	     + fff * diffus(:)** 2 &
      	           * ( gijd(:,1) ** 2 &
      		     + gijd(:,2) ** 2 &
      		     + gijd(:,3) ** 2 &
      		     + 2. &
      		      * ( gijd(:,4) ** 2 &
      		        + gijd(:,5) ** 2 &
      		        + gijd(:,6) ** 2 ) )
        
        tauT = one/sqrt(tauT)
!
        if(idcsclr(1) .ne. 0) then 
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or.  &
                (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
!     
!     determinant of gijd
!     
              fact = one/(gijd(:,1) * gijd(:,2) * gijd(:,3) &
                   - gijd(:,2) * gijd(:,6) * gijd(:,6) &
                   - gijd(:,1) * gijd(:,5) * gijd(:,5) &
                   - gijd(:,3) * gijd(:,4) * gijd(:,4) &
                   + gijd(:,6) * gijd(:,4) * gijd(:,5) * two)
!
! ... note between compressible and incompressible 5 and 6 of giju 
!     are switched        
!
              giju(:,1) = fact * (gijd(:,2)*gijd(:,3)  &
                        - gijd(:,5)**2)
              giju(:,2) = fact * (gijd(:,1)*gijd(:,3)  &
                        - gijd(:,6)**2)
              giju(:,3) = fact * (gijd(:,1)*gijd(:,2) &
                        - gijd(:,4)**2)
              giju(:,4) = fact * (gijd(:,5)*gijd(:,6) &
                        - gijd(:,4)*gijd(:,3) )
              giju(:,5) = fact * (gijd(:,4)*gijd(:,6) &
                        - gijd(:,1)*gijd(:,5) )
              giju(:,6) = fact * (gijd(:,4)*gijd(:,5) &
                        - gijd(:,6)*gijd(:,2) )

!
           endif
        endif                   ! end of idcsclr.ne.0
!     
!.... return
!
        return
        end

