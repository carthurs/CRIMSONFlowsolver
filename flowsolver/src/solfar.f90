      subroutine SolFlow(y,          ac,         u, &
                         yold,       acold,      uold, &
                         x,          xdist,      xdnv, &
                         iBC,        BC,         res,              &
                         nPermDims,  nTmpDims,   aperm, &
                         atemp,      iper,        &
                         ilwork,     shp,        shgl,  &
                         shpb,       shglb,      rowp,      &
                         colm,       lhsK,       lhsP,  &
                         solinc,     rerr,              &
                         memLS_lhs,  memLS_ls,   memLS_nFaces  )
!
!----------------------------------------------------------------------
!
! This is the 2nd interface routine to the Farzin's linear equation
! solver library that uses the CGP and GMRES methods.
!
! input:
!  y      (nshg,ndof)           : Y-variables at n+alpha_f
!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
!  yold   (nshg,ndof)           : Y-variables at beginning of step
!  acold   (nshg,ndof)          : Primvar. accel. at beginning of step
!  x      (numnp,nsd)            : node coordinates
!  iBC    (nshg)                : BC codes
!  BC     (nshg,ndofBC)         : BC constraint parameters
!  iper   (nshg)                : periodic nodal information
!
! output:
!  res    (nshg,nflow)           : preconditioned residual
!  y      (nshg,ndof)           : Y-variables at n+alpha_f
!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
!
!
! The followings are preliminary steps required to use Farzin's
! solver library.  New way of writing has to be used such as
!
!          |  K     G | | du |    | Rmom  |
!          |          | |    | =  |       |
!          | G^t    C | | dp |    | Rcon  |
!
!          |     E    | | dT | =  | Rtemp |
!
!     where
!
!      xKebe : K_ab = dRmom_a/du_b    xTe : E_ab = dRtemp_a/dT_b
!
!              G_ab = dRmom_a/dp_b
!      xGoC  :
!              C_ab = dRcon_a/dp_b
!
!              resf = Rmon Rcon       rest = Rtemp
!
!
! Zdenek Johan,  Winter 1991.  (Fortran 90)
! Juin Kim, Summer 1998. (Incompressible flow solver interface)
! Alberto Figueroa.  CMM-FSI
!----------------------------------------------------------------------
!
      use pointer_data
      use LagrangeMultipliers
!
      use phcommonvars
      use memLS

      use multidomain, only: hrt

      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      !include "mpif.h"
      !include "auxmpi.h"

      TYPE(memLS_lhsType) memLS_lhs
      TYPE(memLS_lsType) memLS_ls
!
      real*8    y(nshg,ndof),             ac(nshg,ndof), &
                yold(nshg,ndof),          acold(nshg,ndof), &
                u(nshg,nsd),              uold(nshg,nsd), &
                x(numnp,nsd),              &
                xdist(nshg), &
                xdnv(nshg,nsd), &
                BC(nshg,ndofBC), &
                res(nshg,nflow), &
                flowDiag(nshg,4), &
                aperm(nshg,nPermDims),    atemp(nshg,nTmpDims), &
                sclrDiag(nshg,1),          &
                lhsK(9,nnz_tot),	  lhsP(4,nnz_tot)
!
      real*8    shp(MAXTOP,maxsh,MAXQPT),   &
                shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                shpb(MAXTOP,maxsh,MAXQPT), &
                shglb(MAXTOP,nsd,maxsh,MAXQPT)
!
      integer   usr(100),                 eqnType, &
                rowp(nshg*nnz),           colm(nshg+1), &
                iBC(nshg),                ilwork(nlwork), &
                iper(nshg)
!
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof), &
                uAlpha(nshg,nsd),          &
                lesP(nshg,4),             lesQ(nshg,4), &
                solinc(nshg,ndof)

      real*8    rerr(nshg,10),            rtmp(nshg,4)

      INTEGER dof, memLS_nFaces, i, j, k, l
      INTEGER, ALLOCATABLE :: incL(:)
      REAL*8, ALLOCATABLE :: faceRes(:), Res4(:,:), Val4(:,:)
!
!.... *******************>> Element Data Formation <<******************
!
!
!.... set the parameters for flux and surface tension calculations
!
!
      idflx = 0
      if(idiff >= 1 )  idflx= (nflow-1) * nsd
      if (isurf == 1) idflx=nflow*nsd
!
!.... compute solution at n+alpha
!
      call itrYAlpha( uold,    yold,    acold,        &
                      u,       y,       ac,             &
                      uAlpha,  yAlpha,  acAlpha)

!
!.... form the LHS matrices, the residual vector (at alpha)
!
      call ElmGMR ( uAlpha,    yAlpha,     acAlpha,     &
                    x,         xdist,      xdnv, &
                    shp,       shgl,       iBC,        &
                    BC,        shpb,       shglb, &
                    res,       iper,       ilwork,    &
                    rowp,      colm,       lhsK,       &
                    lhsP,      rerr   )


      IF (memLSFlag .EQ. 1) THEN
!####################################################################
!     Here calling memLS

      ! *** these allocates are now contained in the section below *** !

      ! ALLOCATE(faceRes(memLS_nFaces), incL(memLS_nFaces))
      ! CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces)

      ! ************** !
      ! if heart model ! 
      ! ************** !

      IF (iheart .GT. int(0)) THEN 
        IF (hrt%isavopen()) THEN
          memLS_nFaces_s = memLS_nFaces + int(1)
          ALLOCATE(faceRes(memLS_nFaces_s), incL(memLS_nFaces_s))
          CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces_s)
        ELSE
          ALLOCATE(faceRes(memLS_nFaces), incL(memLS_nFaces))
          CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces)
        END IF
      ! else flow 
      ELSE
        ALLOCATE(faceRes(memLS_nFaces), incL(memLS_nFaces)) ! 
        CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces)
      END IF

      ! ************** !
      ! ************** !

      incL = 1
      dof = 4
      IF (.NOT.ALLOCATED(Res4)) THEN
         ALLOCATE (Res4(dof,nshg), Val4(dof*dof,nnz_tot))
      END IF

      DO i=1, nshg
         Res4(1:dof,i) = res(i,1:dof)
      END DO

      DO i=1, nnz_tot
         Val4(1:3,i)   = lhsK(1:3,i)
         Val4(5:7,i)   = lhsK(4:6,i)
         Val4(9:11,i)  = lhsK(7:9,i)
         Val4(13:15,i) = lhsP(1:3,i)
         Val4(16,i)    = lhsP(4,i)
      END DO

      !Val4(4:12:4,:) = -lhsP(1:3,:)^t
      DO i=1, nshg
         Do j=colm(i), colm(i+1) - 1
            k = rowp(j)
            DO l=colm(k), colm(k+1) - 1
               IF (rowp(l) .EQ. i) THEN
                  Val4(4:12:4,l) = -lhsP(1:3,j)
                  EXIT
               END IF
            END DO
         END DO
      END DO
      CALL memLS_SOLVE(memLS_lhs, memLS_ls, dof, Res4, Val4, incL, faceRes)

      DO i=1, nshg
         solinc(i,1:dof) = Res4(1:dof,i)
      END DO

!####################################################################
      ELSE
!
!.... lesSolve : main matrix solver
!
      lesId   = numeqns(1)
      eqnType = 1
!
!.... setup the linear algebra solver
!
      rtmp = res(:,1:4)

      !aperm = zero
      call usrNew ( usr,        eqnType,          aperm, &
                    atemp,      rtmp,             solinc,           &
                    flowDiag,   sclrDiag,         lesP,    &
                    lesQ,       iBC,              BC, &
                    iper,       ilwork,           numpe, &
                    nshg,       nshl,             nPermDims,   &
                    nTmpDims,   rowp,             colm,      &
                    lhsK,       lhsP,             rdtmp,       &
                    nnz_tot )
!
!.... solve linear system
!
      call myfLesSolve ( lesId, usr )
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, nflow, 'out')
      endif

      if(Lagrange .gt. zero) then
         call CalcNANBLagrange(colm, rowp, solinc(:,1:3))
         call LagMultiplyMatrix(solinc, 0, nsrflistLagrange, &
            numLagrangeSrfs)
         Lagincr(:,1:3) = (- resL(:,1:3) - AddLag(:,1:3) ) &
            /ScaleFactor(1,1)/alfi/gami/two
      endif

      ENDIF

      call rstatic (res, y, solinc) ! output flow stats
!
!.... end
!
      return
      end

      subroutine SolSclr(y,          ac,         u, &
                         yold,       acold,      uold, &
                         x,          iBC, &
                         BC,         nPermDimsS,  nTmpDimsS,   &
                         apermS,     atempS,     iper,        &
                         ilwork,     shp,        shgl,  &
                         shpb,       shglb,      rowp,      &
                         colm,       lhsS,       solinc)
!
!----------------------------------------------------------------------
!
! This is the 2nd interface routine to the linear equation
! solver library.
!
! input:
!  y      (nshg,ndof)           : Y-variables at n+alpha_f
!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
!  yold   (nshg,ndof)           : Y-variables at beginning of step
!  x      (numnp,nsd)            : node coordinates
!  iBC    (nshg)                : BC codes
!  BC     (nshg,ndofBC)         : BC constraint parameters
!  iper   (nshg)                : periodic nodal information
!
! output:
!  y      (nshg,ndof)           : Y-variables at n+alpha_f
!  ac     (nshg,ndof)           : Primvar. accel. variable n+alpha_m
!
!
! The followings are preliminary steps required to use LesLib
! solver library.  New way of writing has to be used such as
!
!          |     E    | | dS | =  | RScal |
!
!----------------------------------------------------------------------
!
      use pointer_data

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"
!
      real*8    y(nshg,ndof),             ac(nshg,ndof), &
                yold(nshg,ndof),          acold(nshg,ndof), &
                u(nshg,nsd),              uold(nshg,nsd), &
                x(numnp,nsd),             BC(nshg,ndofBC), &
                res(nshg,1), &
                flowDiag(nshg,4), &
                sclrDiag(nshg,1),           lhsS(nnz_tot), &
                apermS(nshg,nPermDimsS),  atempS(nshg,nTmpDimsS)

!
      real*8    shp(MAXTOP,maxsh,MAXQPT),   &
                shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                shpb(MAXTOP,maxsh,MAXQPT), &
                shglb(MAXTOP,nsd,maxsh,MAXQPT)
!
      integer   usr(100),                 eqnType, &
                rowp(nshg*nnz),           colm(nshg+1), &
                iBC(nshg),                ilwork(nlwork), &
                iper(nshg)
!
      real*8    yAlpha(nshg,ndof),        acAlpha(nshg,ndof), &
                uAlpha(nshg,nsd), &
                lesP(nshg,1),               lesQ(nshg,1), &
                solinc(nshg,1)

!
!.... *******************>> Element Data Formation <<******************
!
!.... compute solution at n+alpha
!
      call itrYAlpha( uold,    yold,    acold,  &
                      u,       y,       ac,   &
                      uAlpha,  yAlpha,  acAlpha)
!
!.... form the LHS matrices, the residual vector (at alpha)
!
      call ElmGMRSclr (yAlpha,    acAlpha,    x, &
                       shp,       shgl,       iBC,        &
                       BC,        shpb,       shglb, &
                       res,       iper,       ilwork,    &
                       rowp,      colm,       lhsS   )

!
!.... lesSolve : main matrix solver
!
      lesId   = numeqns(1+nsolt+isclr)
      eqnType = 2
!
!.... setup the linear algebra solver
!
      call usrNew ( usr,        eqnType,          apermS, &
                    atempS,     res,              solinc,           &
                    flowDiag,   sclrDiag,         lesP,    &
                    lesQ,       iBC,              BC, &
                    iper,       ilwork,           numpe, &
                    nshg,       nshl,             nPermDimsS,   &
                    nTmpDimsS,  rowp,             colm,      &
                    rlhsK,      rlhsP,            lhsS,       &
                    nnz_tot )
!
!.... solve linear system
!
      call myfLesSolve ( lesId, usr )
      call getSol ( usr, solinc )

      if (numpe > 1) then
         call commu ( solinc, ilwork, 1, 'out')
      endif

      nsolsc=5+isclr
      call rstaticSclr (res, y, solinc, nsolsc) ! output scalar stats
!
!.... end
!
      return
      end





