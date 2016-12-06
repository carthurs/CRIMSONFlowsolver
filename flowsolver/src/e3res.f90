      subroutine e3Res ( u1,        u2,         u3, &
                         uMesh1, uMesh2, uMesh3, & !ALE variables added MAF 08/10/2016
                         uBar,      aci,        WdetJ, &
                         g1yi,      g2yi,       g3yi, &
                         rLui,      rmu,        rho, &
                         tauC,      tauM,       tauBar, &
                         shpfun,    shg,        src, &
                         rl,        pres,       acl, &
                         rlsli)
!------------------------------------------------------------------------
! 
!  This routine computes the residual vector at an
!  integration point.
!
!  input:
!     u1(npro)                  : x1-velocity
!     u2(npro)                  : x2-velocity
!     u3(npro)                  : x3-velocity
!     uBar(npro,3)              : u - tauM * Li
!     aci(npro,3)               : acceleration
!     rlsli(npro,6)             : resolved Leonard stresses
!     WdetJ(npro)               : weighted jacobian determinant
!     g1yi(npro,ndof)              : x1-gradient of variables
!     g2yi(npro,ndof)              : x2-gradient of variables
!     g3yi(npro,ndof)              : x3-gradient of variables
!     rLui(npro,3)              : total residual of NS equations
!     rmu(npro)                 : fluid viscosity
!     rho(npro)                 : density
!     tauC(npro)                : continuity tau
!     tauM(npro)                : momentum tau
!     tauBar(npro)              : additional tau
!     shpfun(npro,nshl)         : element shape functions
!     shg(npro,nshl,nsd)        : global grad of element shape functions
!     src(npro,nsd)             : body force term
!
!  output:
!     rl(npro,nshl,nflow)
!
!------------------------------------------------------------------------
      use phcommonvars
      use ale
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension u1(npro),         u2(npro),         u3(npro), &
                uBar(npro,nsd),   aci(npro,nsd),    WdetJ(npro), &
                g1yi(npro,nflow), g2yi(npro,nflow), g3yi(npro,nflow), &
                rLui(npro,nsd),   rmu(npro),        rho(npro), &
                tauC(npro),       tauM(npro),       tauBar(npro), &
                shpfun(npro,nshl),shg(npro,nshl,nsd), src(npro,nsd), &
                pres(npro)
      
      dimension rl(npro,nshl,nflow), &
                acl(npro,nshl,ndof), &
                rlsli(npro,6)
!
!.... local declarations
!
      real*8    tmp1(npro),   tmp2(npro),          tmp3(npro),  &
                tmp(npro),    rGNa(npro,nsd,nsd),  rNa(npro,nsd), &
                locmass(npro,nshl),omega(3)

      integer aa

!
!..... ALE variables
!
      dimension uMesh1(npro), uMesh2(npro), uMesh3(npro)

      !     get mesh velocity KDL, MAF
      ! if (rigidOn.eq.1) then
      ! uMesh1(:) = -1.0d0*globalRigidVelocity(1)
      ! uMesh2(:) = -1.0d0*globalRigidVelocity(2)
      ! uMesh3(:) = -1.0d0*globalRigidVelocity(3)
      ! else
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
!.... initialize multipliers for Na and Na_{,i}
!
      rNa  = zero
      rGNa = zero
!
!.... compute the Na multiplier
!
      tmps     = one-flmpr  ! consistant mass factor

!
! no density yet...it comes later
!
      rNa(:,1) = aci(:,1)  * tmps &
               - src(:,1)
      rNa(:,2) = aci(:,2)  * tmps &
               - src(:,2)
      rNa(:,3) = aci(:,3)  * tmps &
               - src(:,3)

!
!.... rotating frame terms if needed
!

      if(matflg(6,1).eq.1) then ! rotation

         omega(1)=datmat(1,6,1)
         omega(2)=datmat(2,6,1)
         omega(3)=datmat(3,6,1)
!
! no density yet...it comes later
!
         
         rNa(:,1) = rNa(:,1) + (omega(2)-omega(3)) * tauM * rLui(:,1)
         rNa(:,2) = rNa(:,2) + (omega(3)-omega(1)) * tauM * rLui(:,2)
         rNa(:,3) = rNa(:,3) + (omega(1)-omega(2)) * tauM * rLui(:,3)
      endif


!
!.... compute the Na,i multiplier
!
      tmp  = -pres + tauC * (g1yi(:,2) + g2yi(:,3) + g3yi(:,4)) !think this is div(v)
      tmp1 =  rmu * ( g2yi(:,2) + g1yi(:,3) )  ! component 12 of viscous stress, SL, MAF 
      tmp2 =  rmu * ( g3yi(:,3) + g2yi(:,4) )  ! component 23 of viscous stress, SL, MAF 
      tmp3 =  rmu * ( g1yi(:,4) + g3yi(:,2) )  ! component 13 of viscous stress, SL, MAF


      if(iconvflow.eq.2) then  ! advective form (NO IBP either) !IBP = integration by parts, MAF

!
! no density yet...it comes later
!
         rNa(:,1) = rNa(:,1)  &
                  + (ubar(:,1)-uMesh1(:)) * g1yi(:,2) &
                  + (ubar(:,2)-uMesh2(:)) * g2yi(:,2) &
                  + (ubar(:,3)-uMesh3(:)) * g3yi(:,2)
         rNa(:,2) = rNa(:,2) &
                  + (ubar(:,1)-uMesh1(:)) * g1yi(:,3) &
                  + (ubar(:,2)-uMesh2(:)) * g2yi(:,3) &
                  + (ubar(:,3)-uMesh3(:)) * g3yi(:,3)
         rNa(:,3) = rNa(:,3) &
                  + (ubar(:,1)-uMesh1(:)) * g1yi(:,4) &
                  + (ubar(:,2)-uMesh2(:)) * g2yi(:,4) &
                  + (ubar(:,3)-uMesh3(:)) * g3yi(:,4)

         rGNa(:,1,1) = two * rmu * g1yi(:,2) + tmp ! stress term 
         rGNa(:,1,2) = tmp1
         rGNa(:,1,3) = tmp3
         rGNa(:,2,1) = tmp1
         rGNa(:,2,2) = two * rmu * g2yi(:,3) + tmp
         rGNa(:,2,3) = tmp2
         rGNa(:,3,1) = tmp3
         rGNa(:,3,2) = tmp2
         rGNa(:,3,3) = two * rmu * g3yi(:,4) + tmp
      else   ! conservative form (with IBP)

!                                            IBP conservative convection
!                                                      |||||
!                                                      vvvvv
         rGNa(:,1,1) = two * rmu * g1yi(:,2) + tmp - u1(:)*u1(:)*rho(:)
         rGNa(:,1,2) = tmp1                        - u1(:)*u2(:)*rho(:)
         rGNa(:,1,3) = tmp3                        - u1(:)*u3(:)*rho(:)
         rGNa(:,2,1) = tmp1                        - u1(:)*u2(:)*rho(:)
         rGNa(:,2,2) = two * rmu * g2yi(:,3) + tmp - u2(:)*u2(:)*rho(:)
         rGNa(:,2,3) = tmp2                        - u3(:)*u2(:)*rho(:)
         rGNa(:,3,1) = tmp3                        - u1(:)*u3(:)*rho(:)
         rGNa(:,3,2) = tmp2                        - u3(:)*u2(:)*rho(:)
         rGNa(:,3,3) = two * rmu * g3yi(:,4) + tmp - u3(:)*u3(:)*rho(:)
      endif
!      if((iLES.gt.10).and.(iLES.lt.20)) then    ! bard
!         rGNa(:,1,1) = rGNa(:,1,1) - rlsli(:,1)*rho(:)
!         rGNa(:,1,2) = rGNa(:,1,2) - rlsli(:,4)*rho(:)
!         rGNa(:,1,3) = rGNa(:,1,3) - rlsli(:,5)*rho(:)
!         rGNa(:,2,1) = rGNa(:,2,1) - rlsli(:,4)*rho(:)
!         rGNa(:,2,2) = rGNa(:,2,2) - rlsli(:,2)*rho(:)
!         rGNa(:,2,3) = rGNa(:,2,3) - rlsli(:,6)*rho(:)
!         rGNa(:,3,1) = rGNa(:,3,1) - rlsli(:,5)*rho(:)
!         rGNa(:,3,2) = rGNa(:,3,2) - rlsli(:,6)*rho(:)
!         rGNa(:,3,3) = rGNa(:,3,3) - rlsli(:,3)*rho(:)
!      endif
   
      tmp1        = tauM * rLui(:,1) 
      tmp2        = tauM * rLui(:,2) 
      tmp3        = tauM * rLui(:,3)
      
      rGNa(:,1,1) = rGNa(:,1,1) + tmp1 * (u1 - uMesh1)
      rGNa(:,1,2) = rGNa(:,1,2) + tmp1 * (u2 - uMesh2)
      rGNa(:,1,3) = rGNa(:,1,3) + tmp1 * (u3 - uMesh3)
      rGNa(:,2,1) = rGNa(:,2,1) + tmp2 * (u1 - uMesh1)
      rGNa(:,2,2) = rGNa(:,2,2) + tmp2 * (u2 - uMesh2)
      rGNa(:,2,3) = rGNa(:,2,3) + tmp2 * (u3 - uMesh3)
      rGNa(:,3,1) = rGNa(:,3,1) + tmp3 * (u1 - uMesh1)
      rGNa(:,3,2) = rGNa(:,3,2) + tmp3 * (u2 - uMesh2)
      rGNa(:,3,3) = rGNa(:,3,3) + tmp3 * (u3 - uMesh3)

      if(iconvflow.eq.1) then  
!
!... get the u_j w_{i,i} term in there to match A_j^T w_{i,j} tau L_i
!    to match the SUPG of incompressible limit
! 
         rGNa(:,1,1) = rGNa(:,1,1) + tmp1 * u1
         rGNa(:,1,2) = rGNa(:,1,2) + tmp2 * u1
         rGNa(:,1,3) = rGNa(:,1,3) + tmp3 * u1
         rGNa(:,2,1) = rGNa(:,2,1) + tmp1 * u2
         rGNa(:,2,2) = rGNa(:,2,2) + tmp2 * u2
         rGNa(:,2,3) = rGNa(:,2,3) + tmp3 * u2
         rGNa(:,3,1) = rGNa(:,3,1) + tmp1 * u3
         rGNa(:,3,2) = rGNa(:,3,2) + tmp2 * u3
         rGNa(:,3,3) = rGNa(:,3,3) + tmp3 * u3
      endif

      if(iconvflow.eq.2) then  ! advective form has a taubar term to restore con
         tmp1 = tauBar &
           * ( rLui(:,1) * g1yi(:,2) &
             + rLui(:,2) * g2yi(:,2) &
             + rLui(:,3) * g3yi(:,2) )
         tmp2 = tauBar &
           * ( rLui(:,1) * g1yi(:,3) &
             + rLui(:,2) * g2yi(:,3) &
             + rLui(:,3) * g3yi(:,3) )
         tmp3 = tauBar &
           * ( rLui(:,1) * g1yi(:,4) &
             + rLui(:,2) * g2yi(:,4) &
             + rLui(:,3) * g3yi(:,4) )

         rGNa(:,1,1) = rGNa(:,1,1) + tmp1 * rLui(:,1)
         rGNa(:,1,2) = rGNa(:,1,2) + tmp1 * rLui(:,2)
         rGNa(:,1,3) = rGNa(:,1,3) + tmp1 * rLui(:,3)
         rGNa(:,2,1) = rGNa(:,2,1) + tmp2 * rLui(:,1)
         rGNa(:,2,2) = rGNa(:,2,2) + tmp2 * rLui(:,2)
         rGNa(:,2,3) = rGNa(:,2,3) + tmp2 * rLui(:,3)
         rGNa(:,3,1) = rGNa(:,3,1) + tmp3 * rLui(:,1)
         rGNa(:,3,2) = rGNa(:,3,2) + tmp3 * rLui(:,2)
         rGNa(:,3,3) = rGNa(:,3,3) + tmp3 * rLui(:,3)
      endif   ! end of advective form
!
!.... everything that gets multiplied by rNa was supposed
!     to have density multiplying it.  Do it now.

      rNa(:,1) = rNa(:,1) * rho
      rNa(:,2) = rNa(:,2) * rho
      rNa(:,3) = rNa(:,3) * rho

!
!
!.... multiply the residual pieces by the weight space
!
      do aa = 1,nshl
!
!.... continuity
!
         rl(:,aa,4) = rl(:,aa,4) + WdetJ &
                    * ( shg(:,aa,1) * uBar(:,1)  &
                      + shg(:,aa,2) * uBar(:,2)  &
                      + shg(:,aa,3) * uBar(:,3) )
!
!.... momentum
!
         rl(:,aa,1) = rl(:,aa,1) - WdetJ &
                    * ( shpfun(:,aa) * rNa(:,1) &
                      + shg(:,aa,1) * rGNa(:,1,1) &
                      + shg(:,aa,2) * rGNa(:,1,2) &
                      + shg(:,aa,3) * rGNa(:,1,3) )
         rl(:,aa,2) = rl(:,aa,2) - WdetJ &
                    * ( shpfun(:,aa) * rNa(:,2) &
                      + shg(:,aa,1) * rGNa(:,2,1) &
                      + shg(:,aa,2) * rGNa(:,2,2) &
                      + shg(:,aa,3) * rGNa(:,2,3) )
         rl(:,aa,3) = rl(:,aa,3) - WdetJ &
                    * ( shpfun(:,aa) * rNa(:,3) &
                      + shg(:,aa,1) * rGNa(:,3,1) &
                      + shg(:,aa,2) * rGNa(:,3,2) &
                      + shg(:,aa,3) * rGNa(:,3,3) )
      
      enddo                 
!
!.... return
!
      return
      end



!------------------------------------------------------------------------
!
!     calculate the residual for the advection-diffusion equation
!
!------------------------------------------------------------------------
      subroutine e3ResSclr ( uMod,              gGradS, &
                             Sclr,		Sdot,	gradS,   &
                             WdetJ,		rLS,	tauS, &
                             shpfun,            shg,    src, &
                             diffus, &
                             rl )
!
       use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      real*8    uMod(npro,nsd),   gGradS(npro, nsd), &
                Sclr(npro),       Sdot(npro),	gradS(npro,nsd), &
                WdetJ(npro),      rLS(npro),	rho(npro), &
                tauS(npro),       shpfun(npro,nshl), src(npro),  &
                shg(npro,nshl,3), rl(npro,nshl)
      
      real*8    diffus(npro)
!
!.... local declarations
!
      real*8    rGNa(npro,nsd),   rNa(npro),  rcp(npro), tmp(npro)

      integer   aa
!     
!.... initialize multipliers for Na and Na_{,i}
!
      rNa  = zero
      rGNa = zero
!
!.... Na multiplier
!
      tmps     = one-flmpr  ! consistant mass factor
      rcp = one ! rho * cp
      

         rNa = rcp*(tmps*Sdot + uMod(:,1) * gradS(:,1) &
                               + uMod(:,2) * gradS(:,2) &
                               + uMod(:,3) * gradS(:,3) )  &
              - src



      tmp = rcp * tauS * (rLS -src)
!
!.... Na,i multiplier
!
      rGNa(:,1) = diffus * gradS(:,1) + uMod(:,1) * tmp
      rGNa(:,2) = diffus * gradS(:,2) + uMod(:,2) * tmp
      rGNa(:,3) = diffus * gradS(:,3) + uMod(:,3) * tmp
!
      if (idcsclr(1) .ne. 0) then
         if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or.  &
              (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
!
!.... add the contribution of DC to residual
!
            rGNa(:,1) = rGNa(:,1) + gGradS(:,1) ! gGradS is 
            rGNa(:,2) = rGNa(:,2) + gGradS(:,2) ! g^{ij}*Y_{j}*dcFct
            rGNa(:,3) = rGNa(:,3) + gGradS(:,3) ! calculated in e3dc.f
!
         endif
      endif                     ! end of idcsclr
!
!.... multiply the residual pieces by the weight space
!
      do aa = 1,nshl
!
         rl(:,aa) = rl(:,aa)	- WdetJ &
                              * ( shpfun(:,aa) * rNa(:) &
                              + shg(:,aa,1) * rGNa(:,1) &
                              + shg(:,aa,2) * rGNa(:,2) &
                              + shg(:,aa,3) * rGNa(:,3) )

      enddo
!
!.... return
!
      return
      end



!----------------------------------------------------------------------
!     Calculate the strong PDE residual.
!----------------------------------------------------------------------
      subroutine e3resStrongPDE( &
           aci,  u1,   u2,   u3,   Temp, rho,  xx, &
                 g1yi, g2yi, g3yi, &
           rLui, src, divqi)
      use phcommonvars
      use ale
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     INPUTS
      double precision, intent(in), dimension(npro,nsd) ::  &
           aci, xx
      double precision, intent(in), dimension(npro,nflow) ::  &
           g1yi, g2yi, g3yi
      double precision, intent(in), dimension(npro) :: &
           u1, u2, u3, Temp, rho
!     OUTPUTS
      double precision, intent(out), dimension(npro,nsd) :: &
           rLui, src
!     LOCALS
      double precision, dimension(npro) :: &
           divu
      double precision, dimension(npro,nsd) :: &
           divqi
      double precision, dimension(nsd) :: &
           omega
!     MESH VELOCITY TERMS
      double precision, dimension(npro) :: uMesh1, uMesh2, uMesh3

      integer   i
      logical :: exist
!.... compute source term
      src = zero
      if(matflg(5,1) .ge. 1) then
!        body force contribution to src
         bfx      = datmat(1,5,1) ! Boussinesq, g*alfap
         bfy      = datmat(2,5,1)
         bfz      = datmat(3,5,1)
         select case ( matflg(5,1) )
         case ( 1 )             ! standard linear body force
            src(:,1) = bfx
            src(:,2) = bfy
            src(:,3) = bfz
         case ( 2 )             ! boussinesq body force
            Tref = datmat(2,2,1)
            src(:,1) = bfx * (Temp(:)-Tref)
            src(:,2) = bfy * (Temp(:)-Tref)
            src(:,3) = bfz * (Temp(:)-Tref)
         case ( 3 )             ! user specified f(x,y,z)
            !call e3source(xx, src)
         end select
      endif
!     
      if(matflg(6,1).eq.1) then
!        coriolis force contribution to src
         omega(1)=datmat(1,6,1)
         omega(2)=datmat(2,6,1)
         omega(3)=datmat(3,6,1)
!  note that we calculate f as if it contains the usual source
!  plus the Coriolis and the centrifugal forces taken to the rhs (sign change)
!  as long as we are doing SUPG with no accounting for these terms in the
!  LHS this is the only change (which will find its way to the RHS momentum
!  equation (both Galerkin and SUPG parts)).
!
!  uncomment later if you want rotation always about z axis
!                 orig_src - om x om x r       - two om x u
!
!$$$          src(:,1)=src(:,1)+omega(3)*omega(3)*xx(:,1)+two*omega(3)*u2
!$$$          src(:,2)=src(:,2)+omega(3)*omega(3)*xx(:,2)-two*omega(3)*u1
!
! more general for testing
!
         src(:,1)=src(:,1) &
              -omega(2)*(omega(1)*xx(:,2)-omega(2)*xx(:,1)) &
              -omega(3)*(omega(1)*xx(:,3)-omega(3)*xx(:,1)) &
              -two*(omega(2)*u3-omega(3)*u2)
         src(:,2)=src(:,2) &
              -omega(1)*(omega(2)*xx(:,1)-omega(1)*xx(:,2)) &
              -omega(3)*(omega(2)*xx(:,3)-omega(3)*xx(:,2)) &
              -two*(omega(3)*u1-omega(1)*u3)
         src(:,3)=src(:,3) &
              -omega(1)*(omega(3)*xx(:,1)-omega(1)*xx(:,3)) &
              -omega(2)*(omega(3)*xx(:,2)-omega(2)*xx(:,3)) &
              -two*(omega(1)*u2-omega(2)*u1)
      endif
!     get mesh velocity
      ! uMesh1(:) = globalMeshVelocity(1)
      ! uMesh2(:) = globalMeshVelocity(2)
      ! uMesh3(:) = globalMeshVelocity(3)


      ! if (rigidOn.eq.1) then
      ! uMesh1(:) = -1.0d0*globalRigidVelocity(1)
      ! uMesh2(:) = -1.0d0*globalRigidVelocity(2)
      ! uMesh3(:) = -1.0d0*globalRigidVelocity(3)
      ! else
      if (aleRigid.eq.1) then
        uMesh1(:) = globalRigidVelocity(1)
        uMesh2(:) = globalRigidVelocity(2)
        uMesh3(:) = globalRigidVelocity(3)
      else 
        uMesh1(:) = real(0.0,8)
        uMesh2(:) = real(0.0,8)
        uMesh3(:) = real(0.0,8)
      endif
      ! endif

      ! MAYBE PRECALCULATE OUTSIDE ?
      !if (iALE .eq. 1) then
      !  call calculateMeshVelocity(uMesh1,uMesh2,uMesh3)
      !end if


!     calculate momentum residual
      rLui(:,1) =(aci(:,1) + (u1 - uMesh1) * g1yi(:,2) &
                           + (u2 - uMesh2) * g2yi(:,2) &
                           + (u3 - uMesh3) * g3yi(:,2) - src(:,1) ) * rho &
                           + g1yi(:,1) &
                           - divqi(:,1)
      rLui(:,2) =(aci(:,2) + (u1 - uMesh1) * g1yi(:,3) &
                           + (u2 - uMesh2) * g2yi(:,3) &
                           + (u3 - uMesh3) * g3yi(:,3) - src(:,2) ) * rho &
                           + g2yi(:,1) &
                           - divqi(:,2)
      rLui(:,3) =(aci(:,3) + (u1 - uMesh1) * g1yi(:,4) &
                           + (u2 - uMesh2) * g2yi(:,4) &
                           + (u3 - uMesh3) * g3yi(:,4) - src(:,3) ) * rho &
                           + g3yi(:,1) &
                           - divqi(:,3)

      if(iconvflow.eq.1) then

         divu(:)  = (g1yi(:,2) + g2yi(:,3) + g3yi(:,4))*rho
         rLui(:,1)=rlui(:,1)+u1(:)*divu(:)
         rLui(:,2)=rlui(:,2)+u2(:)*divu(:)
         rLui(:,3)=rlui(:,3)+u3(:)*divu(:)
      endif

#if DEBUG_ALE == 1      

      write(*,*) 'printing inside e3resStrongPDE'    
      inquire(file="divu.dat", exist=exist)
      if (exist) then
        open(794, file="divu.dat", status="old", position="append", action="write")
      else
        open(794, file="divu.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(794,'(1(e40.20))') (g1yi(i,2) + g2yi(i,3) + g3yi(i,4))      
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      end do 
      close(794)


      inquire(file="rho_vector.dat", exist=exist)
      if (exist) then
        open(794, file="rho_vector.dat", status="old", position="append", action="write")
      else
        open(794, file="rho_vector.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(794,'(1(e40.20))') rho(i)     
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      end do 
      close(794)

      ! write(*,*) 'printing rlui'    
      inquire(file="rlui.dat", exist=exist)
      if (exist) then
        open(793, file="rlui.dat", status="old", position="append", action="write")
      else
        open(793, file="rlui.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(793,'(3(e40.20))') rlui(i,1), rlui(i,2), rlui(i,3)                                     
      end do 
      close(793)

      inquire(file="aci.dat", exist=exist)
      if (exist) then
        open(795, file="aci.dat", status="old", position="append", action="write")
      else
        open(795, file="aci.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(795,'(3(e40.20))') aci(i,1), aci(i,2), aci(i,3)                                     
      end do 
      close(795)


      inquire(file="vrel.dat", exist=exist)
      if (exist) then
        open(796, file="vrel.dat", status="old", position="append", action="write")
      else
        open(796, file="vrel.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(796,'(3(e40.20))') u1(i) - uMesh1(i), &
                                      u2(i) - uMesh2(i), &
                                      u3(i) - uMesh3(i)           
      end do 
      close(796)


      inquire(file="vfluid.dat", exist=exist)
      if (exist) then
        open(797, file="vfluid.dat", status="old", position="append", action="write")
      else
        open(797, file="vfluid.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(797,'(3(e40.20))') u1(i), &
                                      u2(i), &
                                      u3(i)           
      end do 
      close(797)

      inquire(file="vmesh.dat", exist=exist)
      if (exist) then
        open(798, file="vmesh.dat", status="old", position="append", action="write")
      else
        open(798, file="vmesh.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(798,'(3(e40.20))') uMesh1(i), &
                                      uMesh2(i), &
                                      uMesh3(i)           
      end do 
      close(798)


      !!! here, print divq
      inquire(file="divq.dat", exist=exist)
      if (exist) then
        open(799, file="divq.dat", status="old", position="append", action="write")
      else
        open(799, file="divq.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(799,'(3(e40.20))') divqi(i,1), divqi(i,2), divqi(i,3)                                    
      end do 
      close(799)

      ! write(*,*) 'printing rlui'    g1yi(:,4)
      inquire(file="g1yi.dat", exist=exist)
      if (exist) then
        open(800, file="g1yi.dat", status="old", position="append", action="write")
      else
        open(800, file="g1yi.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(800,'(3(e40.20))') g1yi(i,2), g1yi(i,3), g1yi(i,4)                                     
      end do 
      close(800)

      inquire(file="g2yi.dat", exist=exist)
      if (exist) then
        open(801, file="g2yi.dat", status="old", position="append", action="write")
      else
        open(801, file="g2yi.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(801,'(3(e40.20))') g2yi(i,2), g2yi(i,3), g2yi(i,4)                                     
      end do 
      close(801)

      inquire(file="g3yi.dat", exist=exist)
      if (exist) then
        open(802, file="g3yi.dat", status="old", position="append", action="write")
      else
        open(802, file="g3yi.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(802,'(3(e40.20))') g3yi(i,2), g3yi(i,3), g3yi(i,4)                                     
      end do 
      close(802)


      inquire(file="srcterm.dat", exist=exist)
      if (exist) then
        open(803, file="srcterm.dat", status="old", position="append", action="write")
      else
        open(803, file="srcterm.dat", status="new", action="write")
      end if
      do i = 1, npro
             write(803,'(3(e40.20))') src(i,1), src(i,2), src(i,3)                                     
      end do
      close(803)



#endif
!
      return
      end subroutine e3resStrongPDE




!----------------------------------------------------------------------
!     Calculate the strong PDE residual.
!----------------------------------------------------------------------
      subroutine e3resStrongPDE_ALE( &
           aci,  u1,   u2,   u3,&
              uMesh1,uMesh2,uMesh3,&
              Temp, rho,  xx, &
                 g1yi, g2yi, g3yi, &
           rLui, src, divqi)
      use phcommonvars
      use ale
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     INPUTS
      double precision, intent(in), dimension(npro,nsd) ::  &
           aci, xx
      double precision, intent(in), dimension(npro,nflow) ::  &
           g1yi, g2yi, g3yi
      double precision, intent(in), dimension(npro) :: &
           u1, u2, u3, Temp, rho
!     OUTPUTS
      double precision, intent(out), dimension(npro,nsd) :: &
           rLui, src
!     LOCALS
      double precision, dimension(npro) :: &
           divu
      double precision, dimension(npro,nsd) :: &
           divqi
      double precision, dimension(nsd) :: &
           omega
!     MESH VELOCITY TERMS
      ! double precision, dimension(npro) :: uMesh1, uMesh2, uMesh3
      dimension uMesh1(npro), uMesh2(npro), uMesh3(npro)

      integer   i
      logical :: exist
!.... compute source term
      src = zero
      if(matflg(5,1) .ge. 1) then
!        body force contribution to src
         bfx      = datmat(1,5,1) ! Boussinesq, g*alfap
         bfy      = datmat(2,5,1)
         bfz      = datmat(3,5,1)
         select case ( matflg(5,1) )
         case ( 1 )             ! standard linear body force
            src(:,1) = bfx
            src(:,2) = bfy
            src(:,3) = bfz
         case ( 2 )             ! boussinesq body force
            Tref = datmat(2,2,1)
            src(:,1) = bfx * (Temp(:)-Tref)
            src(:,2) = bfy * (Temp(:)-Tref)
            src(:,3) = bfz * (Temp(:)-Tref)
         case ( 3 )             ! user specified f(x,y,z)
            !call e3source(xx, src)
         end select
      endif
!     
      if(matflg(6,1).eq.1) then
!        coriolis force contribution to src
         omega(1)=datmat(1,6,1)
         omega(2)=datmat(2,6,1)
         omega(3)=datmat(3,6,1)
!  note that we calculate f as if it contains the usual source
!  plus the Coriolis and the centrifugal forces taken to the rhs (sign change)
!  as long as we are doing SUPG with no accounting for these terms in the
!  LHS this is the only change (which will find its way to the RHS momentum
!  equation (both Galerkin and SUPG parts)).
!
!  uncomment later if you want rotation always about z axis
!                 orig_src - om x om x r       - two om x u
!
!$$$          src(:,1)=src(:,1)+omega(3)*omega(3)*xx(:,1)+two*omega(3)*u2
!$$$          src(:,2)=src(:,2)+omega(3)*omega(3)*xx(:,2)-two*omega(3)*u1
!
! more general for testing
!
         src(:,1)=src(:,1) &
              -omega(2)*(omega(1)*xx(:,2)-omega(2)*xx(:,1)) &
              -omega(3)*(omega(1)*xx(:,3)-omega(3)*xx(:,1)) &
              -two*(omega(2)*u3-omega(3)*u2)
         src(:,2)=src(:,2) &
              -omega(1)*(omega(2)*xx(:,1)-omega(1)*xx(:,2)) &
              -omega(3)*(omega(2)*xx(:,3)-omega(3)*xx(:,2)) &
              -two*(omega(3)*u1-omega(1)*u3)
         src(:,3)=src(:,3) &
              -omega(1)*(omega(3)*xx(:,1)-omega(1)*xx(:,3)) &
              -omega(2)*(omega(3)*xx(:,2)-omega(2)*xx(:,3)) &
              -two*(omega(1)*u2-omega(2)*u1)
      endif



!     calculate momentum residual
      rLui(:,1) =(aci(:,1) + (u1 - uMesh1) * g1yi(:,2) &
                           + (u2 - uMesh2) * g2yi(:,2) &
                           + (u3 - uMesh3) * g3yi(:,2) - src(:,1) ) * rho &
                           + g1yi(:,1) &
                           - divqi(:,1)
      rLui(:,2) =(aci(:,2) + (u1 - uMesh1) * g1yi(:,3) &
                           + (u2 - uMesh2) * g2yi(:,3) &
                           + (u3 - uMesh3) * g3yi(:,3) - src(:,2) ) * rho &
                           + g2yi(:,1) &
                           - divqi(:,2)
      rLui(:,3) =(aci(:,3) + (u1 - uMesh1) * g1yi(:,4) &
                           + (u2 - uMesh2) * g2yi(:,4) &
                           + (u3 - uMesh3) * g3yi(:,4) - src(:,3) ) * rho &
                           + g3yi(:,1) &
                           - divqi(:,3)

      if(iconvflow.eq.1) then

         divu(:)  = (g1yi(:,2) + g2yi(:,3) + g3yi(:,4))*rho
         rLui(:,1)=rlui(:,1)+u1(:)*divu(:)
         rLui(:,2)=rlui(:,2)+u2(:)*divu(:)
         rLui(:,3)=rlui(:,3)+u3(:)*divu(:)
      endif

!
      return
      end subroutine e3resStrongPDE_ALE