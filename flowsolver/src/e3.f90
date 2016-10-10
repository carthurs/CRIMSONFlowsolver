        subroutine e3 (yl,      acl,     dwl,     shp, &
                       uMeshl, &   !ALE variables added MAF 06/10/2016))
                       shgl,    xl,      rl,      ql, &
                       xKebe,   xGoC,    xmudmi,  sgn, &
                       rerrl, rlsl)
!                                                                      
!----------------------------------------------------------------------
!
!     This routine calculates the residual and tangent matrix for the 
!     UBar formulation of the incompressible Navier Stokes equations.
!
!
! input:    e    a   1..5   when we think of U^e_a  and U is 5 variables
!  yl     (npro,nshl,ndof)      : Y variables (not U)
!  acl    (npro,nshl,ndof)      : Y acceleration (Y_{,t})
!  shp    (nen,ngauss)           : element shape-functions  N_a
!  shgl   (nsd,nen,ngauss)       : element local-shape-functions N_{a,xi}
!  wght   (ngauss)               : element weight (for quadrature)
!  xl     (npro,nenl,nsd)       : nodal coordinates at current step (x^e_a)
!  ql     (npro,nshl,nsd*nsd) : diffusive flux vector (don't worry)
!  rlsl   (npro,nshl,6)       : resolved Leonard stresses
!
! output:
!  rl     (npro,nshl,nflow)      : element RHS residual    (G^e_a)
!  rml    (npro,nshl,nflow)      : element modified residual  (G^e_a tilde)
!  xKebe  (npro,9,nshl,nshl)  : element LHS tangent mass matrix
!  xGoC   (npro,4,nshl,nshl)    : element LHS tangent mass matrix
!  cfl_loc(npro) 		: CFL of the element
!
! Note: This routine will calculate the element matrices for the
!        Hulbert's generalized alpha method integrator
!
! Mathematics done by:  Michel Mallet, Farzin Shakib (version 1)
!                       Farzin Shakib                (version 2)
!
! K. E. Jansen,   Winter 1999.   (advective form formulation)
! C. H. Whiting,  Winter 1999.   (advective form formulation)
!----------------------------------------------------------------------
!
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension yl(npro,nshl,ndof), &
                  acl(npro,nshl,ndof), &      
                  shp(nshl,ngauss),       shgl(nsd,nshl,ngauss), &
                  xl(npro,nenl,nsd),      dwl(npro,nenl), &
                  rl(npro,nshl,nflow),     ql(npro,nshl,idflx)
!      

      	dimension xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
!
!.... local declarations
!
        dimension g1yi(npro,ndof),        g2yi(npro,ndof), &
                  g3yi(npro,ndof),        shg(npro,nshl,nsd), &
                  aci(npro,3),            dxidx(npro,nsd,nsd), &      
                  WdetJ(npro),            rho(npro), &
                  pres(npro),             u1(npro), &
                  u2(npro),               u3(npro), &
                  rLui(npro,nsd),         uBar(npro,nsd), &
                  xmudmi(npro,ngauss),     sgn(npro,nshl), &
                  shpfun(npro,nshl),      shdrv(npro,nsd,nshl), &
                  rmu(npro),              tauC(npro), &
                  tauM(npro),             tauBar(npro), &
                  src(npro,3)

        dimension rlsl(npro,nshl,6),      rlsli(npro,6)

        real*8    rerrl(npro,nshl,6)
        integer   aa, i
        logical :: exist

        dimension uMeshl(npro,nshl,3) !MAF 07/10/2016
        real*8 uMesh1(npro), uMesh2(npro), uMesh3(npro) !MAF 07/10/2016
 
!
!     
!.... local reconstruction of diffusive flux vector for quadratics
!     or greater but NOT for bflux since local mass was not mapped
!
        if ( idiff==2 .and. ires .eq. 1 ) then
           call e3ql (yl,        dwl,       shp,       shgl, &
                      xl,        ql,        xmudmi, &
                      sgn)
        endif
!
!.... loop through the integration points
!

        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
!
!.... get the hierarchic shape functions at this int point
!
        call getshp(shp,          shgl,      sgn, &
                    shpfun,       shdrv)
!
!.... get necessary fluid properties (including eddy viscosity)
!
        call getdiff(dwl,  yl,     shpfun,     xmudmi, xl,   rmu, rho)
!
!.... calculate the integration variables
!
        call e3ivar (yl,          acl,       shpfun, &
                     uMeshl, & !ALE variables added MAF 07/10/2016
                     uMesh1, uMesh2, uMesh3, & !ALE variables added MAF 07/10/2016
                     shdrv,       xl, &
                     aci,         g1yi,      g2yi, &
                     g3yi,        shg,       dxidx, &  
                     WdetJ,       rho,       pres, &
                     u1,          u2,        u3, &          
                     ql,          rLui,      src, &
                     rerrl,       rlsl,      rlsli, &
                     dwl) 

#if DEBUG_ALE == 1

      ! write(*,*) 'printing rlui'
      ! open(793,file='rlui.dat',status='unknown')
      ! do i = 1, npro
      !    write(793,'(3(e20.10))') rlui(i,1), rlui(i,2), rlui(i,3)                                     
      ! end do 
      ! close(793)
      ! stop


  ! write(*,*) 'printing rlui'    
  ! inquire(file="rlui.dat", exist=exist)
  ! if (exist) then
  !   open(793, file="rlui.dat", status="old", position="append", action="write")
  ! else
  !   open(793, file="rlui.dat", status="new", action="write")
  ! end if
  ! do i = 1, npro
  !        write(793,'(3(e20.10))') rlui(i,1), rlui(i,2), rlui(i,3)                                     
  ! end do 
  ! close(793)


#endif 
!
!.... compute the stabilization terms
!
        call e3stab (rho,          u1,       u2, &
                     u3,           &
                     uMesh1, uMesh2, uMesh3, & !ALE variables added MAF 08/10/2016
                     dxidx,    rLui, &   
                     rmu,          tauC,     tauM,  & 
                     tauBar,       uBar )  
!
!.... compute the residual contribution at this integration point
!
        call e3Res ( u1,        u2,         u3, &
                     uMesh1, uMesh2, uMesh3, & !ALE variables added MAF 08/10/2016
                     uBar,      aci,        WdetJ, &
                     g1yi,      g2yi,       g3yi, &
                     rLui,      rmu,        rho, &
                     tauC,      tauM,       tauBar, &
                     shpfun,    shg,        src, &
                     rl,        pres,       acl, &
                     rlsli)
!
!.... compute the tangent matrix contribution
!
        if (lhs .eq. 1) then
           call e3LHS ( u1,        u2,         u3, &
                        uMesh1, uMesh2, uMesh3, & !ALE variables added MAF 09/10/2016
                        uBar,      WdetJ,      rho, &
                        rLui,      rmu, &
                        tauC,      tauM,       tauBar, &
                        shpfun,    shg,        xKebe, &
                        xGoC )
        endif

!
!.... end of integration loop
!
      enddo

!
!.... symmetrize C
!
      if (lhs .eq. 1) then
         do ib = 1, nshl
            do iaa = 1, ib-1
               xGoC(:,4,iaa,ib) = xGoC(:,4,ib,iaa)
            enddo
         enddo
      endif
!
!.... return
!
      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!###################################################################


      subroutine e3Sclr (yl,      acl,     shp, &
                           shgl,    xl,      dwl, &
                           rl,      ql,      xSebe, &   
                           sgn,     xmudmi)
!                                                                      
!----------------------------------------------------------------------
!
!     This routine calculates the residual and tangent matrix for the 
!     advection - diffusion equation for scalar.
!
!----------------------------------------------------------------------
!
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      real*8    yl(npro,nshl,ndof),     acl(npro,nshl,ndof), &
                  shp(nshl,ngauss),       shgl(nsd,nshl,ngauss), &
                  xl(npro,nenl,nsd),      rl(npro,nshl), &
                  ql(npro,nshl,nsd),      xSebe(npro,nshl,nshl), &
                  dwl(npro,nshl)
!
!.... local declarations
!
      real*8    gradS(npro,nsd),        shg(npro,nshl,nsd), &
                  Sdot(npro),             Sclr(npro), &
                  dxidx(npro,nsd,nsd),    WdetJ(npro), &     
                  u1(npro),     u2(npro), u3(npro), &
                  sgn(npro,nshl),         shpfun(npro,nshl), &
                  shdrv(npro,nsd,nshl),   rLS(npro), &
                  tauS(npro),             diffus(npro), &
                  srcL(npro),             srcR(npro), &
                  gGradS(npro,nsd),       dcFct(npro), &
                  giju(npro,6)
!
!.... Source terms sometimes take the form (beta_i)*(phi,_i).  Since
!     the convective term has (u_i)*(phi,_i), it is useful to treat
!     beta_i as a "correction" to the velocity.  In calculating the
!     stabilization terms, the new "modified" velocity (u_i-beta_i) is 
!     then used in place of the pure velocity for stabilization terms,
!     and the source term sneaks into the RHS and LHS.
      real*8    uMod(npro,nsd), srcRat(npro), xmudmi(npro,ngauss)
!
      integer   aa, b
!     
!.... local reconstruction of diffusive flux vector
!
        if ( idiff==2 ) then
           call e3qlSclr (yl, dwl, shp, shgl, xl, ql, sgn)
        endif
!
!.... loop through the integration points
!
        do intp = 1, ngauss

        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
!
!.... get the hierarchic shape functions at this int point
!
        call getshp(shp,          shgl,      sgn, &
                    shpfun,        shdrv)
!
!.... get necessary fluid properties
!
        call getdiffsclr(shpfun,dwl,yl,diffus)
!
!.... calculate the integration variables
!
        call e3ivarSclr(yl,          acl,       shpfun, &
                        shdrv,       xl,        xmudmi, &
                        Sclr,        Sdot,      gradS, &
                        shg,         dxidx,     WdetJ, &     
                        u1,          u2,        u3, &              
                        ql,          rLS,       SrcR, &
                        SrcL,        uMod,      dwl, &
                        diffus,      srcRat)


!
!.... compute the stabilization terms
!
        call e3StabSclr (uMod,    dxidx,   tauS, &
                         diffus,  srcR,    giju, &
                         srcRat)
!
!... computing the DC factor for the discontinuity capturing
!
        if (idcsclr(1) .ne. 0) then
           if ((idcsclr(2).eq.1 .and. isclr.eq.1) .or. &
                (idcsclr(2).eq.2 .and. isclr.eq.2)) then ! scalar with dc
!
              call e3dcSclr ( gradS,    giju,     gGradS, &
                              rLS,      tauS,     srcR, &
                              dcFct)
           endif
        endif                   !end of idcsclr
!
!.... compute the residual contribution at this integration point
!
        call e3ResSclr ( uMod,      gGradS, &
                         Sclr,      Sdot,       gradS, &
                         WdetJ,     rLS,        tauS, &
                         shpfun,    shg,        srcR, &
                         diffus, &
                         rl )
!
!.... compute the tangent matrix contribution
!
        if (lhs .eq. 1) then
           call e3LHSSclr ( uMod,      giju,       dcFct, &
                            Sclr,      Sdot,       gradS, &
                            WdetJ,     rLS,        tauS, &
                            shpfun,    shg,        srcL, &
                            diffus, &
                            xSebe )

        endif

!
!.... end of integration loop
!
      enddo

!
!.... return
!
      return
      end
