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
                         memLS_lhs,  memLS_ls,   memLS_nFaces, &
                         dispMesh, dispMeshold, uMesh, uMeshold, &
                         xMeshold)                 !uMesh added MAF 06/10/2016
                                   !rest of ALE variables added MAF 03/11/2016
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
      use cpp_interface
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

      real*8    uMesh(nshg,3) !MAF 06/10/2016
      real*8    uMeshold(nshg,3) !MAF 03/11/2016
      real*8    dispMesh(numnp,nsd),  dispMeshold(numnp,nsd) !MAF 03/11/2016
      real*8    uMeshalpha(nshg,3),  dispMeshalpha(numnp,nsd) !MAF 03/11/2016
      real*8    xMeshold(numnp,nsd) !MAF 03/11/2016

      INTEGER dof, memLS_nFaces, i, j, k, l
      INTEGER, ALLOCATABLE :: incL(:)
      REAL*8, ALLOCATABLE :: faceRes(:), Res4(:,:), Val4(:,:)

      integer numBCsWhichAllowFlow
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
                      uAlpha,  yAlpha,  acAlpha, & 
                      uMeshold, dispMeshold, &
                      uMesh, dispMesh, &
                      uMeshalpha, dispMeshalpha)

      if (aleType.ge.3) then ! update fluid mesh coordinates at time step n+alpha MAF 03/11/2016
           x = xMeshold + (dispMeshalpha-dispMeshold)
      endif

!
!.... form the LHS matrices, the residual vector (at alpha)
!
      call ElmGMR ( uAlpha,    yAlpha,     acAlpha,     &
                    x,         xdist,      xdnv, &
                    shp,       shgl,       iBC,        &
                    BC,        shpb,       shglb, &
                    res,       iper,       ilwork,    &
                    rowp,      colm,       lhsK,       &
                    lhsP,      rerr, & 
                    uMesh   ) !ALE variables added MAF 06/10/2016

#if DEBUG_ALE == 1
      write(*,*) 'printing res after elmgmr'
      open(793,file='resafterelmgmr.dat',status='new')
      do i = 1, nshg
         write(793,'(4(e40.20))') res(i,1), res(i,2), res(i,3),&
                                  res(i,4) 
                                          
      end do 
      close(793)
#endif


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
          IF (hrt%isavopen() .eq. 1) THEN
            memLS_nFaces_s = memLS_nFaces + int(1)
            ALLOCATE(faceRes(memLS_nFaces_s), incL(memLS_nFaces_s))
            CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces_s)
          ELSE
            ALLOCATE(faceRes(memLS_nFaces), incL(memLS_nFaces))
            CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces)
          END IF
        ! else flow 
        ELSE
          ! Count any Netlist boundary which is currently in a state which stops flow
          ! across the boundary, due to closed diodes
          ! numBCsWhichAllowFlow = int(0)
          ! call callCPPGetNumberOfNetlistsWhichCurrentlyAllowFlow(numBCsWhichAllowFlow)

          memLS_nFaces_s = memLS_nFaces! + numBCsWhichAllowFlow
          ALLOCATE(faceRes(memLS_nFaces_s), incL(memLS_nFaces_s)) ! 
          CALL AddElmpvsQFormemLS(faceRes, memLS_nFaces_s)
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

#if DEBUG_ALE == 1

      write(*,*) 'printing res4 before inside MEMLS if'
      open(793,file='solinc.res4before.dat',status='new')
      do i = 1, nshg
         write(793,'(4(e40.20))') res4(1,i), res4(2,i), res4(3,i),&
                                  res4(4,i) 
                                          
      end do 
      close(793)

      write(*,*) 'printing Val4 before inside MEMLS if'
      open(793,file='solinc.Val4.dat',status='new')
      do i = 1, nnz_tot
         write(793,'(16(e40.20))') Val4(1,i), Val4(2,i), Val4(3,i),&
                                   Val4(4,i), Val4(5,i), Val4(6,i),&
                                   Val4(7,i), Val4(8,i), Val4(9,i),&
                                   Val4(10,i), Val4(11,i), Val4(12,i),&
                                   Val4(13,i), Val4(14,i), Val4(15,i),&
                                   Val4(16,i)
                                          
      end do 
      close(793)
      
#endif


      CALL memLS_SOLVE(memLS_lhs, memLS_ls, dof, Res4, Val4, incL, faceRes)

#if DEBUG_ALE == 1

      write(*,*) 'printing res4 after inside MEMLS if'
      open(793,file='solinc.res4after.dat',status='new')
      do i = 1, nshg
         write(793,'(4(e40.20))') res4(1,i), res4(2,i), res4(3,i),&
                                  res4(4,i)
                                          
      end do 
      close(793)
      
#endif      

        DO i=1, nshg
           solinc(i,1:dof) = Res4(1:dof,i)
        END DO

#if DEBUG_ALE == 1

      write(*,*) 'printing solinc inside MEMLS if'
      open(793,file='solinc.memLS.dat',status='new')
      do i = 1, nshg
         write(793,'(4(e40.20))') solinc(i,1), solinc(i,2), solinc(i,3),&
                                  solinc(i,4)
                                          
      end do 
      close(793)
      stop

#endif


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

#if DEBUG_ALE == 1
      open(99,file='rtmp.dat',status='new')
      write(*,*) 'printing rtmp inside solfar OUTSIDE MEMLS if'
      do i=1,nshg
          write(99,'(4(e25.15))') rtmp(i,1:4)
      enddo
      close(99)

      open(91,file='flowDiag.dat',status='new')
      write(*,*) 'printing flowDiag inside solfar OUTSIDE MEMLS if'
      do i=1,nshg
          write(91,'(4(e25.15))') flowDiag(i,1:4)
      enddo
      close(91)          
#endif      

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


#if DEBUG_ALE == 1
      open(99,file='first_residual.dat',status='new')
      write(*,*) 'printing rhs inside solfar'
      do i=1,nshg
          write(99,'(4(e25.15))') res(i,1:4)
      enddo
      close(99)

      write(*,*) 'printing  inside MEMLS if'
      open(100,file='lhsK.dat',status='new')
      do i = 1, nnz_tot
         write(100,'(9(e20.10))') lhsK(1,i),lhsK(2,i),lhsK(3,i), &
                                          lhsK(4,i),lhsK(5,i),lhsK(6,i), &
                                          lhsK(7,i),lhsK(8,i),lhsK(9,i)  
      end do 
      close(100)

      write(*,*) 'printing lhsP inside MEMLS if'
      open(101,file='lhsP.dat',status='new')
      do i = 1, nnz_tot
         write(101,'(4(e20.10))') lhsP(1,i),lhsP(2,i),lhsP(3,i), &
                                          lhsP(4,i)
      end do 
      close(101)


      write(*,*) 'printing solinc inside MEMLS if'
      open(793,file='solinc.dat',status='new')
      do i = 1, nshg
         write(793,'(4(e20.10))') solinc(i,1), solinc(i,2), solinc(i,3),&
                                  solinc(i,4)
                                          
      end do 
      close(793)

      stop
#endif


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
                         colm,       lhsS,       solinc, &
                         dispMesh, dispMeshold, uMesh, uMeshold)
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

      real*8    uMesh(nshg,3) !MAF 03/11/2016
      real*8    uMeshold(nshg,3) !MAF 03/11/2016
      real*8    dispMesh(numnp,nsd),  dispMeshold(numnp,nsd) !MAF 03/11/2016
      real*8    uMeshalpha(nshg,3),  dispMeshalpha(numnp,nsd) !MAF 03/11/2016          

                

!
!.... *******************>> Element Data Formation <<******************
!
!.... compute solution at n+alpha
!
      call itrYAlpha( uold,    yold,    acold,  &
                      u,       y,       ac,   &
                      uAlpha,  yAlpha,  acAlpha, &
                      uMeshold, dispMeshold, &   !ALE variables added MAF 03/11/2016
                      uMesh, dispMesh, &
                      uMeshalpha, dispMeshalpha)

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





