      subroutine e3ivar (yl,          acl,       shpfun, &
                         shgl,        xl,        &
                         aci,         g1yi,      g2yi,     &
                         g3yi,        shg,       dxidx,    &
                         WdetJ,       rho,       pres,  &
                         u1,          u2,        u3,               &
                         ql,          rLui,      src, &
                         rerrl,       rlsl,      rlsli, &
                         dwl)
!
!----------------------------------------------------------------------
!
!  This routine computes the variables at integration point.
!
! input:
!  yl     (npro,nshl,ndof)      : primitive variables
!  acl    (npro,nshl,ndof)      : prim.var. accel. 
!  shp    (nen)                 : element shape-functions
!  shgl   (nsd,nen)             : element local-grad-shape-functions
!  xl     (npro,nenl,nsd)       : nodal coordinates at current step
!  ql     (npro,nshl,nsd*nsd) : diffusive flux vector
!  rlsl   (npro,nshl,6)       : resolved Leonard stresses
!
! output:
!  aci    (npro,3)              : primvar accel. variables 
!  g1yi   (npro,ndof)           : grad-y in direction 1
!  g2yi   (npro,ndof)           : grad-y in direction 2
!  g3yi   (npro,ndof)           : grad-y in direction 3
!  shg    (npro,nshl,nsd)       : element global grad-shape-functions
!  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient
!  WdetJ  (npro)                : weighted Jacobian
!  rho    (npro)                : density
!  pres   (npro)                : pressure
!  u1     (npro)                : x1-velocity component
!  u2     (npro)                : x2-velocity component
!  u3     (npro)                : x3-velocity component
!  rLui   (npro,nsd)            : xi-momentum residual
!  src    (npro,nsd)            : body force term (not density weighted)
!  rlsli  (npro,6)              : resolved Leonard stresses at quad pt
!
! locally calculated and used
!  divqi  (npro,nsd+isurf)      : divergence of reconstructed quantity
!
! Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
! Zdenek Johan, Winter 1991. (Fortran 90)
! Kenneth Jansen, Winter 1997. Primitive Variables
! Christian Whiting, Winter 1999. (uBar formulation)
!
!----------------------------------------------------------------------
!
      use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
!  passed arrays
!
      dimension yl(npro,nshl,ndof),        dwl(npro,nenl),        &
                  acl(npro,nshl,ndof),       shpfun(npro,nshl), &
                  shgl(npro,nsd,nshl),       xl(npro,nenl,nsd), &
                  aci(npro,nsd),             g1yi(npro,ndof), &
                  g2yi(npro,ndof),           g3yi(npro,ndof), &
                  shg(npro,nshl,nsd),        dxidx(npro,nsd,nsd), &
                  WdetJ(npro),                &
                  rho(npro),                 pres(npro), &
                  u1(npro),                  u2(npro), &
                  u3(npro),                  divqi(npro,nflow-1+isurf), &
                  ql(npro,nshl,idflx),       rLui(npro,nsd), &
                  src(npro,nsd), Temp(npro),xx(npro,nsd)
!
        dimension tmp(npro), dkei(npro),     dist2w(npro)
!
        dimension rlsl(npro,nshl,6),         rlsli(npro,6)
!
        real*8    rerrl(npro,nshl,6), omega(3), divu(npro)
        dimension gyti(npro,nsd),            gradh(npro,nsd), &
                  sforce(npro,3),            weber(npro), &
                  Sclr(npro)
!
!.... ------------->  Primitive variables at int. point  <--------------
!
!.... compute primitive variables
!
       pres = zero
       u1   = zero
       u2   = zero
       u3   = zero
!
       do n = 1, nshl 
          pres = pres + shpfun(:,n) * yl(:,n,1)
          u1   = u1   + shpfun(:,n) * yl(:,n,2)
          u2   = u2   + shpfun(:,n) * yl(:,n,3)
          u3   = u3   + shpfun(:,n) * yl(:,n,4)
       enddo
       if(matflg(5,1).eq.2) then ! boussinesq body force
          Temp = zero
          do n = 1, nshl
             Temp = Temp + shpfun(:,n) * yl(:,n,5)
          enddo
       endif
       if(matflg(5,1).eq.3.or.matflg(6,1).eq.1) then
!         user-specified body force or coriolis force specified
                 xx = zero
          do n  = 1,nenl
             xx(:,1) = xx(:,1)  + shpfun(:,n) * xl(:,n,1)
             xx(:,2) = xx(:,2)  + shpfun(:,n) * xl(:,n,2)
             xx(:,3) = xx(:,3)  + shpfun(:,n) * xl(:,n,3)
          enddo
       endif
!
       if(iRANS.eq.-2) then ! kay-epsilon
          dist2w = zero
          do n = 1, nenl
             dist2w = dist2w + shpfun(:,n) * dwl(:,n)
          enddo
       endif
!
 
       if( (iLES.gt.10).and.(iLES.lt.20))  then  ! bardina
       rlsli = zero
       do n = 1, nshl 

          rlsli(:,1) = rlsli(:,1) + shpfun(:,n) * rlsl(:,n,1)
          rlsli(:,2) = rlsli(:,2) + shpfun(:,n) * rlsl(:,n,2)
          rlsli(:,3) = rlsli(:,3) + shpfun(:,n) * rlsl(:,n,3)
          rlsli(:,4) = rlsli(:,4) + shpfun(:,n) * rlsl(:,n,4)
          rlsli(:,5) = rlsli(:,5) + shpfun(:,n) * rlsl(:,n,5)
          rlsli(:,6) = rlsli(:,6) + shpfun(:,n) * rlsl(:,n,6)

       enddo
       else
          rlsli = zero
       endif
!
!.... ----------------------->  accel. at int. point  <----------------------
!
       aci = zero
       do n = 1, nshl
          aci(:,1) = aci(:,1) + shpfun(:,n) * acl(:,n,2)
          aci(:,2) = aci(:,2) + shpfun(:,n) * acl(:,n,3)
          aci(:,3) = aci(:,3) + shpfun(:,n) * acl(:,n,4)
       enddo
!
!.... --------------------->  Element Metrics  <-----------------------
!
       call e3metric( xl,         shgl,       dxidx,   &
                      shg,        WdetJ)
!
!.... compute the global gradient of u and P
!
!
       g1yi = zero
       g2yi = zero
       g3yi = zero
       do n = 1, nshl
          g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * yl(:,n,1)
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
!
          g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * yl(:,n,1)
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
!
          g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * yl(:,n,1)
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)
       enddo

       divqi = zero
       idflow = 3
       if ( idiff >= 1 .or. isurf==1 ) then
!     
!.... compute divergence of diffusive flux vector, qi,i
!     
          if(idiff >= 1) then
             do n=1, nshl
                divqi(:,1) = divqi(:,1) + shg(:,n,1)*ql(:,n,1 )  &
                                        + shg(:,n,2)*ql(:,n,4 ) &
                                        + shg(:,n,3)*ql(:,n,7 )

                divqi(:,2) = divqi(:,2) + shg(:,n,1)*ql(:,n,2 )  &
                                        + shg(:,n,2)*ql(:,n,5 ) &
                                        + shg(:,n,3)*ql(:,n,8)

                divqi(:,3) = divqi(:,3) + shg(:,n,1)*ql(:,n,3 )  &
                                        + shg(:,n,2)*ql(:,n,6 ) &
                                        + shg(:,n,3)*ql(:,n,9 )

          enddo

          endif                 !end of idiff
!     
          if (isurf .eq. 1) then   
!     .... divergence of normal calculation (curvature)
             do n=1, nshl
                divqi(:,idflow+1) = divqi(:,idflow+1)  &
                     + shg(:,n,1)*ql(:,n,idflx-2) &
                     + shg(:,n,2)*ql(:,n,idflx-1) &
                     + shg(:,n,3)*ql(:,n,idflx)
             enddo 
!     .... initialization of some variables
             Sclr = zero
             gradh= zero
             gyti = zero
             sforce=zero
             do i = 1, npro
                do n = 1, nshl      
                   Sclr(i) = Sclr(i) + shpfun(i,n) * yl(i,n,6) !scalar
!     
!     .... compute the global gradient of Scalar variable
!     
                   gyti(i,1) = gyti(i,1) + shg(i,n,1) * yl(i,n,6) 
                   gyti(i,2) = gyti(i,2) + shg(i,n,2) * yl(i,n,6)
                   gyti(i,3) = gyti(i,3) + shg(i,n,3) * yl(i,n,6)
!     
                enddo

                if (abs (sclr(i)) .le. epsilon_ls) then
                   gradh(i,1) = 0.5/epsilon_ls * (1.0  &
                        + cos(pi*Sclr(i)/epsilon_ls)) * gyti(i,1)
                   gradh(i,2) = 0.5/epsilon_ls * (1.0  &
                        + cos(pi*Sclr(i)/epsilon_ls)) * gyti(i,2) 
                   gradh(i,3) = 0.5/epsilon_ls * (1.0  &
                        + cos(pi*Sclr(i)/epsilon_ls)) * gyti(i,3)
                endif
             enddo              !end of the loop over npro
!     
! .. surface tension force calculation
! .. divide by density now as it gets multiplied in e3res.f, as surface
!    tension force is already in the form of force per unti volume
!     
             weber(:) = Bo
             sforce(:,1) = -(1.0/weber(:)) * divqi(:,idflow+1) & !x-direction 
                  *gradh(:,1) /rho(:)
             sforce(:,2) = -(1.0/weber(:)) * divqi(:,idflow+1) & !y-direction 
                  *gradh(:,2) /rho(:)
             sforce(:,3) = -(1.0/weber(:)) * divqi(:,idflow+1) & !z-direction 
                  *gradh(:,3) /rho(:)          
!
          endif        ! end of the surface tension force calculation
       endif           ! diffusive flux computation
!
! Calculate strong form of pde as well as the source term
!      
       call e3resStrongPDE( &
            aci,  u1,   u2,   u3,   Temp, rho,  xx, &
                  g1yi, g2yi, g3yi, &
            rLui, src, divqi)
!
!.... take care of the surface tension force term here
!
       if (isurf .eq. 1) then  ! note multiplied by density in e3res.f 
          src(:,1) = src(:,1) + sforce(:,1)
          src(:,2) = src(:,2) + sforce(:,2)
          src(:,3) = src(:,3) + sforce(:,3)
       endif       
!
!.... -------------------> error calculation  <-----------------
!     
       if((ierrcalc.eq.1).and.(nitr.eq.iter)) then
          do ia=1,nshl
             tmp=shpfun(:,ia)*WdetJ(:)
             rerrl(:,ia,1) = rerrl(:,ia,1) + &
                             tmp(:)*rLui(:,1)*rLui(:,1)
             rerrl(:,ia,2) = rerrl(:,ia,2) + &
                             tmp(:)*rLui(:,2)*rLui(:,2)
             rerrl(:,ia,3) = rerrl(:,ia,3) + &
                             tmp(:)*rLui(:,3)*rLui(:,3)

             rerrl(:,ia,4) = rerrl(:,ia,4) + &
                             tmp(:)*divqi(:,1)*divqi(:,1)
             rerrl(:,ia,5) = rerrl(:,ia,5) + &
                             tmp(:)*divqi(:,2)*divqi(:,2)
             rerrl(:,ia,6) = rerrl(:,ia,6) + &
                             tmp(:)*divqi(:,3)*divqi(:,3)
          enddo
       endif
!
!.... ----------------------->  dist. kin energy at int. point  <--------------
!
       
       if (ires .ne. 2 .and. iter.eq.1)  then  !only do at beginning of step
!
! calc exact velocity for a channel at quadrature points.
!
       dkei=0.0
!
       do n = 1, nenl 
          dkei = dkei + shpfun(:,n) * (1.0-xl(:,n,2)**2) !u_ex^~ (in FEM space)
       enddo
          dkei = (u1(:)-dkei)**2 +u2(:)**2  ! u'^2+v'^2
          dkei = dkei*WdetJ  ! mult function*W*det of jacobian to
!                              get this quadrature point contribution
          dke  = dke+sum(dkei) ! we move the sum over elements inside of the
!                              sum over quadrature to save memory (we want
!                              a scalar only)
       endif
!     
!.... return
!
       return
       end

!-----------------------------------------------------------------------
! 
!     Calculate the variables for the scalar advection-diffusion
!     equation.
!
!-----------------------------------------------------------------------
      subroutine e3ivarSclr (yl,          acl,       shpfun, &
                            shgl,        xl,        xmudmi, &
                            Sclr,        Sdot,      gradS,   &
                            shg,         dxidx,     WdetJ, &
                            u1,          u2,        u3,               &
                            ql,          rLS ,       SrcR, &
                            SrcL,        uMod,      dwl, &
                            diffus,      srcRat)
!
      use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
!  passed arrays
!
      dimension yl(npro,nshl,ndof),        acl(npro,nshl,ndof),  &
                Sclr(npro),                Sdot(npro), &
                gradS(npro,nsd),           shpfun(npro,nshl), &
                shgl(npro,nsd,nshl),       xl(npro,nenl,nsd), &
                shg(npro,nshl,nsd),        dxidx(npro,nsd,nsd), &
                WdetJ(npro),               &
                u1(npro),                  u2(npro), &
                u3(npro),                  divS(npro), &
                ql(npro,nshl,nsd),         rLS(npro), &
                SrcR(npro),                 SrcL(npro), &
                dwl(npro,nshl),            diffus(npro), &
                umod(npro,nsd), Temp(npro),xx(npro,nsd), &
                divqi(npro)   
!
      dimension tmp(npro), srcRat(npro)
      real*8 rLui(npro,nsd),     aci(npro,nsd), &
             g1yi(npro,nflow),   g2yi(npro,nflow), &
             g3yi(npro,nflow), &
             src(npro,nsd),      rho(npro), &
             rmu(npro)
      real*8 uBar(npro,nsd), xmudmi(npro,ngauss)

!
!.... ------------->  Primitive variables at int. point  <--------------
!
!.... compute primitive variables
!
      u1   = zero
      u2   = zero
      u3   = zero
      Sclr = zero
!
      id=isclr+5
      do n = 1, nshl 
         u1   = u1   + shpfun(:,n) * yl(:,n,2)
         u2   = u2   + shpfun(:,n) * yl(:,n,3)
         u3   = u3   + shpfun(:,n) * yl(:,n,4)
         Sclr = Sclr + shpfun(:,n) * yl(:,n,id)
      enddo
!
!
!.... ----------------------->  dS/dt at int. point  <----------------------
!
      Sdot = zero
      do n = 1, nshl
         Sdot = Sdot + shpfun(:,n) * acl(:,n,id)
      enddo
!
!.... --------------------->  Element Metrics  <-----------------------
!

      call e3metric( xl,         shgl,        dxidx,   &
                     shg,        WdetJ)

!
!.... compute the global gradient of u and P
!
!
       gradS = zero
       do n = 1, nshl
          gradS(:,1) = gradS(:,1) + shg(:,n,1) * yl(:,n,id)
          gradS(:,2) = gradS(:,2) + shg(:,n,2) * yl(:,n,id)
          gradS(:,3) = gradS(:,3) + shg(:,n,3) * yl(:,n,id)
       enddo

       divS = zero
       if ( idiff >= 1 ) then
!
!.... compute divergence of diffusive flux vector, qi,i
!
          do n=1, nshl
             divS(:) = divS(:) + shg(:,n,1)*ql(:,n,1 )  &
                               + shg(:,n,2)*ql(:,n,2 )  &
                               + shg(:,n,3)*ql(:,n,3 ) 
          enddo
       endif                    ! diffusive flux computation

       if(consrv_sclr_conv_vel .eq. 1) then
!         Calculate uBar = u - TauM*L, where TauM is the momentum
!         stabilization factor and L is the momentum residual

          if(matflg(5,1).eq.2) then ! boussinesq body force
             Temp = zero
             do n = 1, nshl
                Temp = Temp + shpfun(:,n) * yl(:,n,5)
             enddo
          endif
          if(matflg(5,1).eq.3.or.matflg(6,1).eq.1) then
!     user-specified body force or coriolis force specified
             xx = zero
             do n  = 1,nenl
                xx(:,1) = xx(:,1)  + shpfun(:,n) * xl(:,n,1)
                xx(:,2) = xx(:,2)  + shpfun(:,n) * xl(:,n,2)
                xx(:,3) = xx(:,3)  + shpfun(:,n) * xl(:,n,3)
             enddo
          endif
          aci = zero
          do n = 1, nshl
             aci(:,1) = aci(:,1) + shpfun(:,n) * acl(:,n,2)
             aci(:,2) = aci(:,2) + shpfun(:,n) * acl(:,n,3)
             aci(:,3) = aci(:,3) + shpfun(:,n) * acl(:,n,4)
          enddo
          g1yi = zero
          g2yi = zero
          g3yi = zero
          do n = 1, nshl
             g1yi(:,1) = g1yi(:,1) + shg(:,n,1) * yl(:,n,1)
             g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
             g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
             g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
!     
             g2yi(:,1) = g2yi(:,1) + shg(:,n,2) * yl(:,n,1)
             g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
             g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
             g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
!     
             g3yi(:,1) = g3yi(:,1) + shg(:,n,3) * yl(:,n,1)
             g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
             g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
             g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)
          enddo
!          
          if (iLSet .eq. 0)then
             rho  = datmat(1,1,1)
             rmu = datmat(1,2,1)
          else
             write(*,*) 'Not sure if we can handle level set with K-E'
             write(*,*) '(different uMods? correct value of rho?)'
          endif
          divqi=zero  ! until we reconstruct q_flow for scalar solve
          call e3resStrongPDE( &
               aci,  u1,   u2,   u3,   Temp, rho,  x, &
                     g1yi, g2yi, g3yi, &
               rLui, src, divqi)
          src(:,1)=u1           !
          src(:,2)=u2           ! store u in src memory
          src(:,3)=u3           !
!         e3uBar calculates Tau_M and assembles uBar
          call getdiff(dwl, yl, shpfun, xmudmi, xl, rmu, rho)
          call e3uBar(rho, src, dxidx, rLui, rmu, uBar)
          u1=ubar(:,1)          ! the entire scalar residual
          u2=ubar(:,2)          ! is based on the modified
          u3=ubar(:,3)          ! velocity for conservation
       endif
!
!.... Initialize uMod, the modified velocity uMod
!      We initialize it to u_i and then calculate
!      the correction in e3sourcesclr
!

       umod(:,1) = u1
       umod(:,2) = u2
       umod(:,3) = u3
!     
!.... compute  source terms
!
!ad
!ad    if we are solving the redistancing equation, the umod(:,:) are 
!AD    modified in e3sourceSclr.  
!AD
!AD  if we are redistancing levelset variable we want to use a use the  
!AD  convective term from the equation.  


       if(nosource.ne.1) then
        call e3sourceSclr ( Sclr,         Sdot,      gradS,  dwl, &
                            shpfun,    shg,       yl,     dxidx, &
                            diffus,       u1,        u2,     u3, &
                            srcR,         srcL,      uMod,    &
                            srcRat)
       else
        srcRat = zero
        srcR   = zero
        srcL   = zero
       endif
!
!.... -------------------> Scalar residual  <-----------------
!

         rLS(:) = ( Sdot(:) +  (u1*gradS(:,1) +  &
                                    u2*gradS(:,2) + &
                                    u3*gradS(:,3)) ) &
              - divS(:)           

!
!.... return
!
       return
       end

