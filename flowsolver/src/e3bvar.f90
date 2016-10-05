  subroutine e3bvar(yl,      acl,     ul, &
                     iBCB,    BCB, &
                     shpb,    shglb, &
                     xlb,     xdistl,  xdnvl, &
                     lnode,   SWB,            &
                     WdetJb,  bnorm,   pres,  &
                     u1,      u2,      u3,     rmu, &
                     unm,     tau1n,   tau2n,  tau3n, &
                     vdot,    usup1,   usup2,  velsup, &
                     rlKwall, xKebe)
!
!----------------------------------------------------------------------
!
!   This routine computes the variables at integration points for 
! the boundary element routine.
!
! input:
!  yl     (npro,nshl,ndof)      : primitive variables (local)
!          ndof: 5[p,v1,v2,v3,T]+number of scalars solved 
!  acl    (npro,nshl,ndof)      : acceleration (local)
!  ul     (npro,nshlb,nsd)       : displacement (local)
!  iBCB   (npro,ndiBCB)         : boundary condition code
!  BCB    (npro,nshlb,ndBCB)    : boundary Condition values
!  shpb   (nen)                 : boundary element shape-functions
!  shglb  (nsd,nen)             : boundary element grad-shape-functions
!  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
!  lnode  (nenb)                : local nodes on the boundary
!
! output:
!  g1yi   (npro,ndof)           : grad-v in direction 1
!  g2yi   (npro,ndof)           : grad-v in direction 2
!  g3yi   (npro,ndof)           : grad-v in direction 3
!  WdetJb (npro)                : weighted Jacobian
!  bnorm  (npro,nsd)            : outward normal
!  pres   (npro)                : pressure
!  u1     (npro)                : x1-velocity component
!  u2     (npro)                : x2-velocity component
!  u3     (npro)                : x3-velocity component
!  unm    (npro)                : BC u dot n
!  p      (npro)                : BC pressure
!  tau1n  (npro)                : BC viscous flux 1
!  tau2n  (npro)                : BC viscous flux 2
!  tau3n  (npro)                : BC viscous flux 3
!  vdot   (npro,nsd)            : acceleration at quadrature points
!  rlKwall(npro,nshlb,nsd)      : wall stiffness contribution to the local residual
!
! Zdenek Johan, Summer 1990.  (Modified from e2bvar.f)
! Zdenek Johan, Winter 1991.  (Fortran 90)
! Alberto Figueroa, Winter 2004.  CMM-FSI
!----------------------------------------------------------------------
!
      !use turbsa
      use pointer_data
      use phcommonvars

      use deformableWall
      use ale

      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension   yl(npro,nshl,ndof),        rmu(npro), &
                  iBCB(npro,ndiBCB), &
                  BCB(npro,nshlb,ndBCB), &
                  shpb(npro,nshl),           shglb(npro,nsd,nshl), &
                  xlb(npro,nenl,nsd), &
                  lnode(27),                 g1yi(npro,ndof), &
                  g2yi(npro,ndof),           g3yi(npro,ndof), &
                  WdetJb(npro),              bnorm(npro,nsd), &
                  pres(npro), &
                  u1(npro),                  u2(npro), &
                  u3(npro), &
                  unm(npro), &                
                  tau1n(npro),               tau2n(npro), &
                  tau3n(npro), &
                  acl(npro,nshl,ndof),       ul(npro,nshl,nsd), &
                  xdistl(npro,nshl),         xdnvl(npro,nshl,nsd), &
                  vdot(npro,nsd), &
                  usup1(npro,nsd),           usup2(npro,nsd), &
                  velsup(npro,nsd), &
                  rlKwall(npro,nshlb,nsd), &
                  SWB(npro,nProps)
!
      dimension   gl1yi(npro,ndof),          gl2yi(npro,ndof), &
                  gl3yi(npro,ndof),          dxdxib(npro,nsd,nsd), &
                  dxidxb(npro,nsd,nsd),      temp(npro), &
                  temp1(npro),               temp2(npro), &
                  temp3(npro), &
                  v1(npro,nsd),              v2(npro,nsd), &
                  v3(npro,nsd)
     
      dimension   xKebe(npro,nsd**2,nshl,nshl)

      real*8      rotel(npro,nsd,nsd)
      real*8      uterm(npro,nsd)
      real*8      ulrot(npro,nshl,nsd)
      real*8      xlbrot(npro,nenl,nsd)
      real*8      dxdxib_el(npro,nsd,nsd)
      real*8      dxidxb_el(npro,nsd,nsd)
      real*8      shgb_el(npro,nsd,nshl)

      real*8      Dmatrix(npro,5,5)

      real*8      strainterm(npro,5,nshlb)
      real*8      Dtimesstrain(npro,5,nshlb)
      real*8      straindisp_g(npro,5,nshlb,nsd)
      real*8      DtimesB(npro,5,nshlb,nsd)

      real*8      LHSwall(npro,nsd**2,nshl,nshl)

      real*8      Wall_LHSfactor(npro), &
                  Wall_LHSfactorSupp(npro), &
                  Wall_LHSfactorDamp(npro), &
                  tempSuppStiff(npro),       tempSuppDamp(npro)

      integer    numrgndslcl
      integer, allocatable :: rgndslcl(:)
      real*8      xmidpoint(nsd), ringsegmid(nsd)
      real*8      tmpvc1(nsd), tmpvc2(nsd), tmpdirec(nsd)
      real*8      tmplen, tmphfact

      dimension  tmp1(npro)
!     
      real*8      Turb(npro),                xki, &
                  xki3,                      fv1
!        
      integer   e, i, j
!      
      integer   aa, b, iblk

      integer logicPassed

      real*8 uMesh1(npro), uMesh2(npro), uMesh3(npro)
      

      !     get mesh velocity KDL, MA
      ! if (rigidOn.eq.1) then
      ! write(*,*) "rigidOn"
      ! ! uMesh1(:) = 0.0d0 !no relative stabilization for rigid body motion
      ! ! uMesh2(:) = 0.0d0
      ! ! uMesh3(:) = 0.0d0
      ! uMesh1(:) = -1.0d0*globalRigidVelocity(1)
      ! uMesh2(:) = -1.0d0*globalRigidVelocity(2)
      ! uMesh3(:) = -1.0d0*globalRigidVelocity(3)
      ! else
      uMesh1(:) = globalMeshVelocity(1)
      uMesh2(:) = globalMeshVelocity(2)
      uMesh3(:) = globalMeshVelocity(3)
      ! endif


!
!.... ------------------->  integration variables  <--------------------
!
!.... compute the primitive variables at the integration point
!
      pres = zero
      u1   = zero
      u2   = zero
      u3   = zero
!     
        
      do n = 1, nshlb
         nodlcl = lnode(n)
!     
         pres = pres + shpb(:,nodlcl) * yl(:,nodlcl,1)
         u1   = u1   + shpb(:,nodlcl) * yl(:,nodlcl,2)
         u2   = u2   + shpb(:,nodlcl) * yl(:,nodlcl,3)
         u3   = u3   + shpb(:,nodlcl) * yl(:,nodlcl,4)

      enddo
!
!.... ---------------------->  Element Metrics  <-----------------------
!
!.... compute the deformation gradient
!
      dxdxib = zero
!
      do n = 1, nenl
         dxdxib(:,1,1) = dxdxib(:,1,1) + xlb(:,n,1) * shglb(:,1,n)
         dxdxib(:,1,2) = dxdxib(:,1,2) + xlb(:,n,1) * shglb(:,2,n)
         dxdxib(:,1,3) = dxdxib(:,1,3) + xlb(:,n,1) * shglb(:,3,n)
         dxdxib(:,2,1) = dxdxib(:,2,1) + xlb(:,n,2) * shglb(:,1,n)
         dxdxib(:,2,2) = dxdxib(:,2,2) + xlb(:,n,2) * shglb(:,2,n)
         dxdxib(:,2,3) = dxdxib(:,2,3) + xlb(:,n,2) * shglb(:,3,n)
         dxdxib(:,3,1) = dxdxib(:,3,1) + xlb(:,n,3) * shglb(:,1,n)
         dxdxib(:,3,2) = dxdxib(:,3,2) + xlb(:,n,3) * shglb(:,2,n)
         dxdxib(:,3,3) = dxdxib(:,3,3) + xlb(:,n,3) * shglb(:,3,n)
      enddo

      ! #if DEBUG_ALE == 1 
      ! write(*,*) 'inside e3bvar then stop'

      ! open(793,file='shglb.dat',status='new')
      ! do i = 1, npro
      !     do j = 1, nsd
      !     write(793,'(2(i10),4(e20.10))') i,j,shglb(i,j,1), shglb(i,j,2), shglb(i,j,3),&
      !                             shglb(i,j,4) ! nenl = 4
      !     end do                          
      ! end do 
      ! close(793)

      ! open(794,file='shpb.dat',status='new')
      ! do i = 1, npro

      !     write(794,'(1(i10),4(e20.10))') i,shpb(i,1), shpb(i,2), shpb(i,3),&
      !                             shpb(i,4)                      
      ! end do 
      ! close(794)

      ! ! stop

      ! #endif

!
!.... compute the normal to the boundary
!
!$$$      if (lcsyst .eq. 4) then   ! wedge-quad
!$$$         temp1 =  dxdxib(:,2,1) * dxdxib(:,3,3) -
!$$$     &            dxdxib(:,2,3) * dxdxib(:,3,1)
!$$$         temp2 =  dxdxib(:,3,1) * dxdxib(:,1,3) -
!$$$     &            dxdxib(:,3,3) * dxdxib(:,1,1)
!$$$         temp3 =  dxdxib(:,1,1) * dxdxib(:,2,3) -
!$$$     &            dxdxib(:,1,3) * dxdxib(:,2,1)
!$$$      elseif( lcyst .eq. 6) then  ! pyr-tri face
!$$$         temp1 =  dxdxib(:,2,1) * dxdxib(:,3,3) -
!$$$     &            dxdxib(:,2,3) * dxdxib(:,3,1)
!$$$         temp2 =  dxdxib(:,3,1) * dxdxib(:,1,3) -
!$$$     &            dxdxib(:,3,3) * dxdxib(:,1,1)
!$$$         temp3 =  dxdxib(:,1,1) * dxdxib(:,2,3) -
!$$$     &            dxdxib(:,1,3) * dxdxib(:,2,1)
!$$$      elseif( lcyst .eq. 1) then !usual wrong way tets
!$$$         temp1 = -dxdxib(:,2,2) * dxdxib(:,3,1) +
!$$$     &            dxdxib(:,2,1) * dxdxib(:,3,2)
!$$$         temp2 = -dxdxib(:,3,2) * dxdxib(:,1,1) +
!$$$     &            dxdxib(:,3,1) * dxdxib(:,1,2)
!$$$         temp3 = -dxdxib(:,1,2) * dxdxib(:,2,1) +
!$$$     &            dxdxib(:,1,1) * dxdxib(:,2,2)
!$$$      else
!$$$         temp1 =  dxdxib(:,2,2) * dxdxib(:,3,1) -
!$$$     &            dxdxib(:,2,1) * dxdxib(:,3,2)
!$$$         temp2 =  dxdxib(:,3,2) * dxdxib(:,1,1) -
!$$$     &            dxdxib(:,3,1) * dxdxib(:,1,2)
!$$$         temp3 =  dxdxib(:,1,2) * dxdxib(:,2,1) -
!$$$     &            dxdxib(:,1,1) * dxdxib(:,2,2)
!$$$      endif
!$$$c
!$$$      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
!$$$      bnorm(:,1) = temp1 * temp
!$$$      bnorm(:,2) = temp2 * temp
!$$$      bnorm(:,3) = temp3 * temp
!$$$c     
!$$$      WdetJb     = Qwtb(lcsyst,intp) / temp
!$$$      if(lcsyst .eq. 3) WdetJb=WdetJb*two
!
!.... compute the normal to the boundary. This is achieved by taking
!     the cross product of two vectors in the plane of the 2-d 
!     boundary face.
!
      if(lcsyst.eq.1) then      ! set to curl into element all others out
         ipt2=2
         ipt3=3
      elseif(lcsyst.eq.2) then
         ipt2=4
         ipt3=2
      elseif(lcsyst.eq.3) then
         ipt2=3
         ipt3=2
      elseif(lcsyst.eq.4) then
         ipt2=2
         ipt3=4
      elseif(lcsyst.eq.5) then
         ipt2=4
         ipt3=2
      elseif(lcsyst.eq.6) then
         ipt2=2
         ipt3=5
      endif
      v1 = xlb(:,ipt2,:) - xlb(:,1,:)
      v2 = xlb(:,ipt3,:) - xlb(:,1,:)
!
! compute cross product
!
      temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
      temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
      temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
!     
! mag is area for quads, twice area for tris
! 
      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      bnorm(:,1) = temp1 * temp
      bnorm(:,2) = temp2 * temp
      bnorm(:,3) = temp3 * temp
!
      if (lcsyst .eq. 1) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 2) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 3) then
         WdetJb     = Qwtb(lcsyst,intp) / (two*temp)
      elseif (lcsyst .eq. 4) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 5) then
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      elseif (lcsyst .eq. 6) then
         WdetJb     = Qwtb(lcsyst,intp) / (two*temp)
      endif
!
!.... -------------------------->  Grad-V  <----------------------------
!
!.... compute grad-v for Navier-Stokes terms
!
      if (Navier .eq. 1) then
!
!.... compute the inverse of deformation gradient
!
         dxidxb(:,1,1) =   dxdxib(:,2,2) * dxdxib(:,3,3) &
              - dxdxib(:,3,2) * dxdxib(:,2,3)
         dxidxb(:,1,2) =   dxdxib(:,3,2) * dxdxib(:,1,3) &
              - dxdxib(:,1,2) * dxdxib(:,3,3)
         dxidxb(:,1,3) =   dxdxib(:,1,2) * dxdxib(:,2,3) &
              - dxdxib(:,1,3) * dxdxib(:,2,2)
         temp          = one / ( dxidxb(:,1,1) * dxdxib(:,1,1) &
              + dxidxb(:,1,2) * dxdxib(:,2,1) &
              + dxidxb(:,1,3) * dxdxib(:,3,1) )
         dxidxb(:,1,1) =  dxidxb(:,1,1) * temp
         dxidxb(:,1,2) =  dxidxb(:,1,2) * temp
         dxidxb(:,1,3) =  dxidxb(:,1,3) * temp
         dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1) &
              - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
         dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3) &
              - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
         dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3) &
              - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
         dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2) &
              - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
         dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2) &
              - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
         dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2) &
              - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp
!
!.... compute local-grad-Y
!
         gl1yi = zero
         gl2yi = zero
         gl3yi = zero
!     
         do n = 1, nshl
            gl1yi(:,1) = gl1yi(:,1) + shglb(:,1,n) * yl(:,n,1)
            gl1yi(:,2) = gl1yi(:,2) + shglb(:,1,n) * yl(:,n,2)
            gl1yi(:,3) = gl1yi(:,3) + shglb(:,1,n) * yl(:,n,3)
            gl1yi(:,4) = gl1yi(:,4) + shglb(:,1,n) * yl(:,n,4)
!     
            gl2yi(:,1) = gl2yi(:,1) + shglb(:,2,n) * yl(:,n,1)
            gl2yi(:,2) = gl2yi(:,2) + shglb(:,2,n) * yl(:,n,2)
            gl2yi(:,3) = gl2yi(:,3) + shglb(:,2,n) * yl(:,n,3)
            gl2yi(:,4) = gl2yi(:,4) + shglb(:,2,n) * yl(:,n,4)
!     
            gl3yi(:,1) = gl3yi(:,1) + shglb(:,3,n) * yl(:,n,1)
            gl3yi(:,2) = gl3yi(:,2) + shglb(:,3,n) * yl(:,n,2)
            gl3yi(:,3) = gl3yi(:,3) + shglb(:,3,n) * yl(:,n,3)
            gl3yi(:,4) = gl3yi(:,4) + shglb(:,3,n) * yl(:,n,4)
         enddo
!     
!.... convert local-grads to global-grads
!     
         g1yi(:,2) = dxidxb(:,1,1) * gl1yi(:,2) + &
              dxidxb(:,2,1) * gl2yi(:,2) + &
              dxidxb(:,3,1) * gl3yi(:,2)
         g2yi(:,2) = dxidxb(:,1,2) * gl1yi(:,2) + &
              dxidxb(:,2,2) * gl2yi(:,2) + &
              dxidxb(:,3,2) * gl3yi(:,2)
         g3yi(:,2) = dxidxb(:,1,3) * gl1yi(:,2) + &
              dxidxb(:,2,3) * gl2yi(:,2) + &
              dxidxb(:,3,3) * gl3yi(:,2)
!     
         g1yi(:,3) = dxidxb(:,1,1) * gl1yi(:,3) + &
              dxidxb(:,2,1) * gl2yi(:,3) + &
              dxidxb(:,3,1) * gl3yi(:,3)
         g2yi(:,3) = dxidxb(:,1,2) * gl1yi(:,3) + &
              dxidxb(:,2,2) * gl2yi(:,3) + &
              dxidxb(:,3,2) * gl3yi(:,3)
         g3yi(:,3) = dxidxb(:,1,3) * gl1yi(:,3) + &
              dxidxb(:,2,3) * gl2yi(:,3) + &
              dxidxb(:,3,3) * gl3yi(:,3)
!     
         g1yi(:,4) = dxidxb(:,1,1) * gl1yi(:,4) + &
              dxidxb(:,2,1) * gl2yi(:,4) + &
              dxidxb(:,3,1) * gl3yi(:,4)
         g2yi(:,4) = dxidxb(:,1,2) * gl1yi(:,4) + &
              dxidxb(:,2,2) * gl2yi(:,4) + &
              dxidxb(:,3,2) * gl3yi(:,4)
         g3yi(:,4) = dxidxb(:,1,3) * gl1yi(:,4) + &
              dxidxb(:,2,3) * gl2yi(:,4) + &
              dxidxb(:,3,3) * gl3yi(:,4)
!     
!.... end grad-v
!     
      endif
!
!.... mass flux
!
      unm = bnorm(:,1) * u1 +bnorm(:,2) * u2  +bnorm(:,3) * u3
      ! unm = bnorm(:,1) * (u1-uMesh1) + &
      !       bnorm(:,2) * (u2-uMesh2) + &
      !       bnorm(:,3) * (u3-uMesh3)



! no rho in continuity eq.


!
!.... viscous flux
!
      tau1n = bnorm(:,1) * two * rmu *  g1yi(:,2)  &
           + bnorm(:,2) *      (rmu * (g2yi(:,2) + g1yi(:,3))) &
           + bnorm(:,3) *      (rmu * (g3yi(:,2) + g1yi(:,4)))
      tau2n = bnorm(:,1) *      (rmu * (g2yi(:,2) + g1yi(:,3))) &
           + bnorm(:,2) * two * rmu *  g2yi(:,3) &
           + bnorm(:,3) *      (rmu * (g3yi(:,3) + g2yi(:,4)))
      tau3n = bnorm(:,1) *      (rmu * (g3yi(:,2) + g1yi(:,4))) &
           + bnorm(:,2) *      (rmu * (g3yi(:,3) + g2yi(:,4))) &
           + bnorm(:,3) * two * rmu *  g3yi(:,4)
!
      temp1 = bnorm(:,1) * tau1n &
           + bnorm(:,2) * tau2n &
           + bnorm(:,3) * tau3n

      pres  = pres - temp1

      tau1n = tau1n - bnorm(:,1) * temp1
      tau2n = tau2n - bnorm(:,2) * temp1
      tau3n = tau3n - bnorm(:,3) * temp1


      !
      !.... -------------------------->  Deformable Wall  <----------------------------
      !

      if (ideformwall.eq.1) then

          ! compute rotation matrix from global to element local flat coordinates
          v3(:,1) = bnorm(:,1)
          v3(:,2) = bnorm(:,2)
          v3(:,3) = bnorm(:,3)

          temp = one / sqrt ( v1(:,1)**2 + v1(:,2)**2 + v1(:,3)**2 )
          v1(:,1) = v1(:,1) * temp
          v1(:,2) = v1(:,2) * temp
          v1(:,3) = v1(:,3) * temp

          ! cross product again for new v2
          temp1 = v3(:,2) * v1(:,3) - v1(:,2) * v3(:,3)
          temp2 = v1(:,1) * v3(:,3) - v3(:,1) * v1(:,3)
          temp3 = v3(:,1) * v1(:,2) - v1(:,1) * v3(:,2)

          temp     = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
          v2(:,1) = temp1 * temp
          v2(:,2) = temp2 * temp
          v2(:,3) = temp3 * temp

          rotel = zero

          rotel(:,1,1:3) = v1(:,:)
          rotel(:,2,1:3) = v2(:,:)
          rotel(:,3,1:3) = v3(:,:)

          ! rotate the element nodal coordinates to the element local coordinate frame
          xlbrot = zero
          do n = 1, nshlb
              nodlcl = lnode(n)
              do i = 1, 3
                  do j = 1, 3
                      xlbrot(:,nodlcl,i) = xlbrot(:,nodlcl,i)+rotel(:,i,j)*xlb(:,nodlcl,j)
                  enddo
              enddo
          enddo

          ! compute the jacobian and inverse
          ! for transformation between element coordinates and element local coordinates
          dxdxib_el = zero
          do n = 1, nenl
              dxdxib_el(:,1,1) = dxdxib_el(:,1,1) + xlbrot(:,n,1) * shglb(:,1,n)
              dxdxib_el(:,1,2) = dxdxib_el(:,1,2) + xlbrot(:,n,1) * shglb(:,2,n)
              dxdxib_el(:,1,3) = dxdxib_el(:,1,3) + xlbrot(:,n,1) * shglb(:,3,n)
              dxdxib_el(:,2,1) = dxdxib_el(:,2,1) + xlbrot(:,n,2) * shglb(:,1,n)
              dxdxib_el(:,2,2) = dxdxib_el(:,2,2) + xlbrot(:,n,2) * shglb(:,2,n)
              dxdxib_el(:,2,3) = dxdxib_el(:,2,3) + xlbrot(:,n,2) * shglb(:,3,n)
              dxdxib_el(:,3,1) = dxdxib_el(:,3,1) + xlbrot(:,n,3) * shglb(:,1,n)
              dxdxib_el(:,3,2) = dxdxib_el(:,3,2) + xlbrot(:,n,3) * shglb(:,2,n)
              dxdxib_el(:,3,3) = dxdxib_el(:,3,3) + xlbrot(:,n,3) * shglb(:,3,n)
          enddo

          ! compute inverse of jacobian
          dxidxb_el(:,1,1) =   dxdxib_el(:,2,2) * dxdxib_el(:,3,3) - dxdxib_el(:,3,2) * dxdxib_el(:,2,3)
          dxidxb_el(:,1,2) =   dxdxib_el(:,3,2) * dxdxib_el(:,1,3) - dxdxib_el(:,1,2) * dxdxib_el(:,3,3)
          dxidxb_el(:,1,3) =   dxdxib_el(:,1,2) * dxdxib_el(:,2,3) - dxdxib_el(:,1,3) * dxdxib_el(:,2,2)
          temp          = one / ( dxidxb_el(:,1,1) * dxdxib_el(:,1,1) + dxidxb_el(:,1,2) * dxdxib_el(:,2,1) + dxidxb_el(:,1,3) * dxdxib_el(:,3,1) )
          dxidxb_el(:,1,1) =  dxidxb_el(:,1,1) * temp
          dxidxb_el(:,1,2) =  dxidxb_el(:,1,2) * temp
          dxidxb_el(:,1,3) =  dxidxb_el(:,1,3) * temp
          dxidxb_el(:,2,1) = (dxdxib_el(:,2,3) * dxdxib_el(:,3,1) - dxdxib_el(:,2,1) * dxdxib_el(:,3,3)) * temp
          dxidxb_el(:,2,2) = (dxdxib_el(:,1,1) * dxdxib_el(:,3,3) - dxdxib_el(:,3,1) * dxdxib_el(:,1,3)) * temp
          dxidxb_el(:,2,3) = (dxdxib_el(:,2,1) * dxdxib_el(:,1,3) - dxdxib_el(:,1,1) * dxdxib_el(:,2,3)) * temp
          dxidxb_el(:,3,1) = (dxdxib_el(:,2,1) * dxdxib_el(:,3,2) - dxdxib_el(:,2,2) * dxdxib_el(:,3,1)) * temp
          dxidxb_el(:,3,2) = (dxdxib_el(:,3,1) * dxdxib_el(:,1,2) - dxdxib_el(:,1,1) * dxdxib_el(:,3,2)) * temp
          dxidxb_el(:,3,3) = (dxdxib_el(:,1,1) * dxdxib_el(:,2,2) - dxdxib_el(:,1,2) * dxdxib_el(:,2,1)) * temp

          ! boundary shape function gradients with respect to element local flat x,y,z
          ! multiply local shape function gradient with transpose of jacobian inverse
          shgb_el = zero
          do n = 1, nenl
              do i = 1 , 3
                  do j = 1 , 3
                      shgb_el(:,i,n) = shgb_el(:,i,n) + dxidxb_el(:,j,i) * shglb(:,j,n)
                  enddo
              enddo
          enddo

          ! Deformable Wall Residual Terms

          ! stiffness term
          !
          ! define the material properties for the enhanced membrane
          ! the membrane thickness is always SWB(:,1)
          Dmatrix = zero

          if (iUseSWB.eq.0 .or. (iUseSWB.gt.0 .and. iUseSWBthickonly.gt.0) ) then

              ! when the legacy SWB field is not used
              ! only the isotropic case
              ! is implemented so far

              tempcoeff = one/(one-rnuvw**2);

              do iel = 1, npro

                  ! the FORTRAN standard does not guarantee short-circuit evaluation
                  ! of logical statements, so we must guard mBET in this manner to ensure
                  ! we don't read from uninitialised memory
                  logicPassed = int(0)
                  if ((iUseBET .eq. 1) .and. (numWallRegions .gt. 0)) then
                    if (mBET(icurrentblk)%p(iel,WallETagID).gt.0) then
                      logicPassed = int(1)
                    endif
                  endif
                  ! regional values
                  ! supersedes default values
                  if (logicPassed .eq. int(1)) then

                      if (iUseSWBthickonly .eq. 0) then
                          SWB(iel,1) = regionWallProps(mBET(icurrentblk)%p(iel,WallhTagID) ,1)
                      end if
                      SWB(iel,7) = regionWallProps(mBET(icurrentblk)%p(iel,WallETagID), 2) * tempcoeff * one
                      SWB(iel,8) = regionWallProps(mBET(icurrentblk)%p(iel,WallETagID), 2) * tempcoeff * rnuvw
                      SWB(iel,9) = regionWallProps(mBET(icurrentblk)%p(iel,WallETagID), 2) * tempcoeff * pt5*(1-rnuvw)
                      SWB(iel,10) = regionWallProps(mBET(icurrentblk)%p(iel,WallETagID), 2) * tempcoeff * pt5*(1-rnuvw)*rshearconstantvw

                      !write(*,*) icurrentblk,intp,iel
                      !if (mBET(icurrentblk)%p(iel,WallETagID) .gt. 0) then
                      !write(*,*) WallETagID, mBET(icurrentblk)%p(iel,WallETagID), ValueListWallE( mBET(icurrentblk)%p(iel,WallETagID) )
                      !endif

                  else
                      ! default values
                      if (iUseSWBthickonly .eq. 0) then
                          SWB(iel,1) = thicknessvw
                      end if
                      SWB(iel,7) = evw * tempcoeff * one
                      SWB(iel,8) = evw * tempcoeff * rnuvw
                      SWB(iel,9) = evw * tempcoeff * pt5*(1-rnuvw)
                      SWB(iel,10) = evw * tempcoeff * pt5*(1-rnuvw)*rshearconstantvw

                  end if

              enddo

!              Dmatrix(:,1,1) = one
!              Dmatrix(:,2,2) = Dmatrix(:,1,1)
!              Dmatrix(:,1,2) = rnuvw
!              Dmatrix(:,2,1) = Dmatrix(:,1,2)
!              Dmatrix(:,3,3) = pt5*(1-rnuvw)
!              Dmatrix(:,4,4) = pt5*(1-rnuvw)*rshearconstantvw
!              Dmatrix(:,5,5) = pt5*(1-rnuvw)*rshearconstantvw

              !Dmatrix(:,:,:) = evw * one/(one-rnuvw**2)*Dmatrix(:,:,:);

              !do iel = 1, npro

              !    Dmatrix(iel,:,:) = ValueListWallE( iBCB(iel,2) ) * &
              !                       one/(one-rnuvw**2) * &
              !                       Dmatrix(iel,:,:)

              !enddo

          endif

          if (nProps.eq.10) then
              ! This is an Isotropic Material
              Dmatrix(:,1,1) = SWB(:,7)
              Dmatrix(:,2,2) = SWB(:,7)
              Dmatrix(:,1,2) = SWB(:,8)
              Dmatrix(:,2,1) = SWB(:,8)
              Dmatrix(:,3,3) = SWB(:,9)
              Dmatrix(:,4,4) = SWB(:,10)
              Dmatrix(:,5,5) = SWB(:,10)
          elseif(nProps.eq.21) then

              ! This is an Orthotropic Material
              Dmatrix(:,1,1) = SWB(:,7)

              Dmatrix(:,2,1) = SWB(:,8)
              Dmatrix(:,1,2) = SWB(:,8)
              Dmatrix(:,2,2) = SWB(:,9)

              Dmatrix(:,3,1) = SWB(:,10)
              Dmatrix(:,1,3) = SWB(:,10)
              Dmatrix(:,3,2) = SWB(:,11)
              Dmatrix(:,2,3) = SWB(:,11)
              Dmatrix(:,3,3) = SWB(:,12)

              Dmatrix(:,4,1) = SWB(:,13)
              Dmatrix(:,1,4) = SWB(:,13)
              Dmatrix(:,4,2) = SWB(:,14)
              Dmatrix(:,2,4) = SWB(:,14)
              Dmatrix(:,4,3) = SWB(:,15)
              Dmatrix(:,3,4) = SWB(:,15)
              Dmatrix(:,4,4) = SWB(:,16)

              Dmatrix(:,5,1) = SWB(:,17)
              Dmatrix(:,1,5) = SWB(:,17)
              Dmatrix(:,5,2) = SWB(:,18)
              Dmatrix(:,2,5) = SWB(:,18)
              Dmatrix(:,5,3) = SWB(:,19)
              Dmatrix(:,3,5) = SWB(:,19)
              Dmatrix(:,5,4) = SWB(:,20)
              Dmatrix(:,4,5) = SWB(:,20)
              Dmatrix(:,5,5) = SWB(:,21)

          else

              write(*,*) 'Number of wall properties not set correctly!'
              stop

          end if

          ! strain term on the right
          strainterm = zero

          ulrot = zero
          Dtimesstrain = zero

          straindisp_g = zero
          DtimesB = zero
          do n = 1, nshlb

              nodlcl = lnode(n)

              ! strain displacement matrix that computes local strain from global nodal displacements
              do i = 1, 3
                  straindisp_g(:,1,nodlcl,i) = shgb_el(:,1,nodlcl) * rotel(:,1,i)
                  straindisp_g(:,2,nodlcl,i) = shgb_el(:,2,nodlcl) * rotel(:,2,i)

                  straindisp_g(:,3,nodlcl,i) = shgb_el(:,2,nodlcl) * rotel(:,1,i) + &
                                               shgb_el(:,1,nodlcl) * rotel(:,2,i)

                  straindisp_g(:,4,nodlcl,i) = shgb_el(:,1,nodlcl) * rotel(:,3,i)
                  straindisp_g(:,5,nodlcl,i) = shgb_el(:,2,nodlcl) * rotel(:,3,i)
              enddo

              ! compute local strain from nodal displacements
              do j = 1, 5
                  do i = 1, 3
                      ! the "prestress" if not added from the SWB field
                      ! is added here via a reference displacement
                      strainterm(:,j,nodlcl) = strainterm(:,j,nodlcl) + &
                                               straindisp_g(:,j,nodlcl,i) * &
                                               (ul(:,nodlcl,i) + mDisp_refl(icurrentblk)%p(:,nodlcl,i))
                  enddo
              enddo

              ! the strain term multiplied by the D matrix
              do i = 1, 5
                  do j = 1, 5
                      Dtimesstrain(:,i,nodlcl) = Dtimesstrain(:,i,nodlcl) + &
                                                 Dmatrix(:,i,j) * strainterm(:,j,nodlcl)
                  enddo
              enddo

              ! the D matrix multiplied by the strain displacement matrix
              do k = 1, 3
                  do i = 1, 5
                      do j = 1, 5
                          DtimesB(:,i,nodlcl,k) = DtimesB(:,i,nodlcl,k) + &
                                                  Dmatrix(:,i,j) * straindisp_g(:,j,nodlcl,k)
                      enddo
                  enddo
              enddo


          enddo   ! end do n=1,nshlb MA



          rlKwall = zero
          LHSwall = zero

          do n = 1, nshlb

              nodlcln = lnode(n)

              do m = 1, nshlb

                  nodlclm = lnode(m)

                  ! multiply by transpose of global strain displacement matrix
                  ! to get the residual contribution
                  do j = 1, 5
                      do i = 1, 3
                          rlKwall(:,nodlcln,i) = &
                          rlKwall(:,nodlcln,i) + &
                          straindisp_g(:,j,nodlcln,i) * Dtimesstrain(:,j,nodlclm)
                      enddo
                  enddo

                  ! now compute the tangent LHS contribution
                  do k = 1, 3
                      do j = 1, 5
                          do i = 1, 3

                              LHSwall(:,3*(i-1)+k,nodlcln,nodlclm) = &
                              LHSwall(:,3*(i-1)+k,nodlcln,nodlclm) + &
                              straindisp_g(:,j,nodlcln,i) * DtimesB(:,j,nodlclm,k) * &
                              WdetJb * SWB(:,1) * iwallstiffactor*betai*Delt(itseq)*Delt(itseq)*alfi

                          enddo
                      enddo
                  enddo

              enddo ! end do m = 1, nshlb

              ! now add the legacy prestress contribution
              if(iUseSWB.gt.0) then

                  ! multiply element local prestress (elements 2 through 6 of SWB)
                  ! by transpose of global strain displacement matrix

                  do j = 1, 5
                      do i = 1, 3

                          rlKwall(:,nodlcln,i) = rlKwall(:,nodlcln,i) + &
                              straindisp_g(:,j,nodlcln,i) * SWB(:,1+j)

                      enddo
                  enddo
              endif

              ! complete the numerical integration
              ! multiply by weighted jacobian and thickness

              do i = 1, 3
                  rlKwall(:,nodlcln,i) = rlKwall(:,nodlcln,i) * WdetJb * SWB(:,1)
              enddo

          enddo ! end do n=1,nshlb, MA 

          ! mass term

          vdot = zero
      
          do n = 1, nshlb
         
              nodlcl = lnode(n)

              vdot(:,1) = vdot(:,1) + shpb(:,nodlcl) * acl(:,nodlcl,2)
              vdot(:,2) = vdot(:,2) + shpb(:,nodlcl) * acl(:,nodlcl,3)
              vdot(:,3) = vdot(:,3) + shpb(:,nodlcl) * acl(:,nodlcl,4)
            
          enddo
         
          vdot(:,1) = vdot(:,1) * SWB(:,1) * rhovw
          vdot(:,2) = vdot(:,2) * SWB(:,1) * rhovw
          vdot(:,3) = vdot(:,3) * SWB(:,1) * rhovw
         
          Wall_LHSfactor = almi*(one-flmpl)*rhovw*SWB(:,1)
         
          ! tissue support stiffness term
          if (iwallsupp.gt.0) then
         
              usup1 = zero
         
              do n = 1, nshlb
         
                  nodlcl = lnode(n)
               
                  usup1(:,1) = usup1(:,1) + shpb(:,nodlcl) * ul(:,nodlcl,1)
                  usup1(:,2) = usup1(:,2) + shpb(:,nodlcl) * ul(:,nodlcl,2)
                  usup1(:,3) = usup1(:,3) + shpb(:,nodlcl) * ul(:,nodlcl,3)
            
              enddo
         
              usup1(:,1) = usup1(:,1) * tissSuppStiffCoeff
              usup1(:,2) = usup1(:,2) * tissSuppStiffCoeff
              usup1(:,3) = usup1(:,3) * tissSuppStiffCoeff
         
              Wall_LHSfactor = Wall_LHSfactor + &
              Delt(itseq) * alfi * betai * Delt(itseq) * tissSuppStiffCoeff

          endif
         
          if (idistancenudge.gt.0) then
         
              usup2 = zero
         
              do n = 1, nshlb
         
                  nodlcl = lnode(n)

                  usup2(:,1) = usup2(:,1) + shpb(:,nodlcl) * xdistl(:,nodlcl) * xdnvl(:,nodlcl,1)
                  usup2(:,2) = usup2(:,2) + shpb(:,nodlcl) * xdistl(:,nodlcl) * xdnvl(:,nodlcl,2)
                  usup2(:,3) = usup2(:,3) + shpb(:,nodlcl) * xdistl(:,nodlcl) * xdnvl(:,nodlcl,3)
     
              enddo

              usup2(:,1) = usup2(:,1) * stateFilterCoeff
              usup2(:,2) = usup2(:,2) * stateFilterCoeff
              usup2(:,3) = usup2(:,3) * stateFilterCoeff
            
          endif
            
          ! tissue support damping term
          if (iwalldamp.gt.0) then
         
              velsup = zero
         
              do n = 1, nshlb
         
                  nodlcl = lnode(n)
                  velsup(:,1) = velsup(:,1) + shpb(:,nodlcl) * yl(:,nodlcl,2)
                  velsup(:,2) = velsup(:,2) + shpb(:,nodlcl) * yl(:,nodlcl,3)
                  velsup(:,3) = velsup(:,3) + shpb(:,nodlcl) * yl(:,nodlcl,4)
            
              enddo

              velsup(:,1) = velsup(:,1) * tissSuppDampCoeff
              velsup(:,2) = velsup(:,2) * tissSuppDampCoeff
              velsup(:,3) = velsup(:,3) * tissSuppDampCoeff

              Wall_LHSfactor = Wall_LHSfactor + Delt(itseq) * alfi * gami * tissSuppDampCoeff
            
          end if

          ! additional contribution at the rings
          ! from wall support

          ! since we aren't using the built-in quadrature rule
          ! we only add to the residual once through the integration loop
          if (intp .eq. 1 .and. (iringsupp.gt.0 .or. iringdamp.gt.0)) then

              allocate(rgndslcl(nenl))

              do iel = 1, npro

                  ! check if element is on deformable wall
                  if (btest(iBCB(iel,1),4)) then

                      numrgndslcl = 0

                      xmidpoint = zero

                      ringsegmid = zero

                      do ielnode = 1, nenbl

                          ! check if node tag is 'ring'
                          if ( btest(mNodeTagl(icurrentblk)%p(iel,ielnode),0) ) then

                              numrgndslcl = numrgndslcl + 1

                              rgndslcl(numrgndslcl) = ielnode

                              ringsegmid = ringsegmid + xlb(iel,ielnode,:)

                          end if

                          xmidpoint = xmidpoint + xlb(iel,ielnode,:)

                      end do

                      xmidpoint = xmidpoint / nenbl

                      ringsegmid = ringsegmid / numrgndslcl

                      if (numrgndslcl .gt. 1) then

                          tmpdirec = ringsegmid - xmidpoint

                          tmpvc2 = xlb(iel,rgndslcl(2),:) - ringsegmid

                          tmpvc2 = tmpvc2 / sqrt(tmpvc2(1)**2+tmpvc2(2)**2+tmpvc2(3)**2)

                          ! compute longitudinal tethering direction
                          tmpdirec = tmpdirec - tmpvc2*(tmpdirec(1)*tmpvc2(1)+tmpdirec(2)*tmpvc2(2)+tmpdirec(3)*tmpvc2(3))

                          tmpdirec = tmpdirec / sqrt(tmpdirec(1)**2+tmpdirec(2)**2+tmpdirec(3)**2)

                          tmpvc2 = xlb(iel,rgndslcl(2),:)-xlb(iel,rgndslcl(1),:)

                          tmplen = sqrt(tmpvc2(1)**2+tmpvc2(2)**2+tmpvc2(3)**2) / real(6.0,8)



                          ! compute the traction
                          tmpvc1 = zero
                          tmpvc2 = zero

                          if (iringsupp .gt. 0) then
                              tmpvc1 = tmpvc1 + tissSuppRingStiffCoeff * ul(iel,rgndslcl(1),:)
                              tmpvc2 = tmpvc2 + tissSuppRingStiffCoeff * ul(iel,rgndslcl(2),:)
                          end if

                          if (iringdamp .gt. 0) then
                              tmpvc1 = tmpvc1 + tissSuppRingDampCoeff * yl(iel,rgndslcl(1),2:4)
                              tmpvc2 = tmpvc2 + tissSuppRingDampCoeff * yl(iel,rgndslcl(1),2:4)
                          end if

                          tmpvc1 = tmpdirec*(tmpdirec(1)*tmpvc1(1)+tmpdirec(2)*tmpvc1(2)+tmpdirec(3)*tmpvc1(3))
                          tmpvc2 = tmpdirec*(tmpdirec(1)*tmpvc2(1)+tmpdirec(2)*tmpvc2(2)+tmpdirec(3)*tmpvc2(3))

                          rlKwall(iel,rgndslcl(1),:) = rlKwall(iel,rgndslcl(1),:) + &
                                                             tmplen * (two * tmpvc1 + tmpvc2)
                          rlKwall(iel,rgndslcl(2),:) = rlKwall(iel,rgndslcl(2),:) + &
                                                             tmplen * (tmpvc1 + two * tmpvc2)

                          ! calculate tangent contribution
                          tmphfact = zero

                          if (iringsupp .gt. 0) then
                              tmphfact = tmphfact + tmplen * tissSuppRingStiffCoeff * betai*Delt(itseq)*Delt(itseq)*alfi
                          end if

                          if (iringdamp .gt. 0) then
                              tmphfact = tmphfact + tmplen * tissSuppRingDampCoeff  * Delt(itseq)*alfi*gami
                          end if


                          do i = 1, 3
                              do k = 1, 3
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(1),rgndslcl(1)) = &
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(1),rgndslcl(1)) + tmpdirec(i) * tmpdirec(k) * tmphfact * two
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(1),rgndslcl(2)) = &
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(1),rgndslcl(2)) + tmpdirec(i) * tmpdirec(k) * tmphfact
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(2),rgndslcl(1)) = &
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(2),rgndslcl(1)) + tmpdirec(i) * tmpdirec(k) * tmphfact
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(2),rgndslcl(2)) = &
                                  LHSwall(iel,3*(i-1)+k,rgndslcl(2),rgndslcl(2)) + tmpdirec(i) * tmpdirec(k) * tmphfact * two
                              end do
                          end do

                          !write(*,*) ringsegmid, tmpvc1 + ringsegmid

                      end if

                  end if

              end do

              if (allocated(rgndslcl)) then
                deallocate(rgndslcl)
              endif

          end if

                   
          ! Wall LHS tangent contributions
         
          Wall_LHSfactor = iwallmassfactor*Wall_LHSfactor*WdetJb

          !.... NOTE:  the wall mass contributions should only have 3 nodal components
          !.... since the fourth node is an interior node... therefore, the loops should
          !.... be done from 1 to nshlb=3...
    
          do b = 1, nshlb
              do aa = 1, nshlb
            
                  tmp1 = Wall_LHSfactor * shpb(:,aa) * shpb(:,b)

                  xKebe(:,1,aa,b) = xKebe(:,1,aa,b) + tmp1
                  xKebe(:,5,aa,b) = xKebe(:,5,aa,b) + tmp1
                  xKebe(:,9,aa,b) = xKebe(:,9,aa,b) + tmp1
               
              enddo
          enddo

          xKebe = xKebe + LHSwall

123       continue

       endif
  !
  !.... return
  !
  return
  end
      
  !---------------------------------------------------------------------
  !
  !     variables for boundary elements
  !
  !---------------------------------------------------------------------
  subroutine e3bvarSclr (yl,        shdrv,    xlb, &
  shapeVar,     WdetJb,   bnorm, &
  flux,      dwl )

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      !
      dimension yl(npro,nshl,ndof),        shdrv(npro,nsd,nshl), &
      xlb(npro,nenl,nsd),        shapeVar(npro,nshl), &
      WdetJb(npro),              bnorm(npro,nsd), &
      flux(npro)
      !
      dimension dxdxib(npro,nsd,nsd), &
      dxidxb(npro,nsd,nsd),      temp(npro), &
      temp1(npro),               temp2(npro), &
      temp3(npro), &
      v1(npro,nsd),              v2(npro,nsd), &
      gradSl(npro,nsd),          gradS(npro,nsd)

      real*8    diffus(npro),              dwl(npro,nshl)
        
      call getdiffsclr(shapeVar,dwl,yl,diffus)
      !
      !.... ---------------------->  Element Metrics  <-----------------------
      !
      !.... compute the deformation gradient
      !
      dxdxib = zero
      !
      do n = 1, nenl
          dxdxib(:,1,1) = dxdxib(:,1,1) + xlb(:,n,1) * shdrv(:,1,n)
          dxdxib(:,1,2) = dxdxib(:,1,2) + xlb(:,n,1) * shdrv(:,2,n)
          dxdxib(:,1,3) = dxdxib(:,1,3) + xlb(:,n,1) * shdrv(:,3,n)
          dxdxib(:,2,1) = dxdxib(:,2,1) + xlb(:,n,2) * shdrv(:,1,n)
          dxdxib(:,2,2) = dxdxib(:,2,2) + xlb(:,n,2) * shdrv(:,2,n)
          dxdxib(:,2,3) = dxdxib(:,2,3) + xlb(:,n,2) * shdrv(:,3,n)
          dxdxib(:,3,1) = dxdxib(:,3,1) + xlb(:,n,3) * shdrv(:,1,n)
          dxdxib(:,3,2) = dxdxib(:,3,2) + xlb(:,n,3) * shdrv(:,2,n)
          dxdxib(:,3,3) = dxdxib(:,3,3) + xlb(:,n,3) * shdrv(:,3,n)
      enddo
      !
      !.... compute the normal to the boundary. This is achieved by taking
      !     the cross product of two vectors in the plane of the 2-d
      !     boundary face.
      !
      v1 = xlb(:,2,:) - xlb(:,1,:)
      v2 = xlb(:,3,:) - xlb(:,1,:)
        
      !
      !.....The following are done in order to correct temp1..3
      !     based on the results from compressible code.  This is done only
      !     for wedges, depending on the bounary face.(tri or quad)
      !
      if (lcsyst .eq. 4) then
          temp1 = dxdxib(:,2,1) * dxdxib(:,3,3) - &
          dxdxib(:,2,3) * dxdxib(:,3,1)
          temp2 = dxdxib(:,3,1) * dxdxib(:,1,3) - &
          dxdxib(:,3,3) * dxdxib(:,1,1)
          temp3 = dxdxib(:,1,1) * dxdxib(:,2,3) - &
          dxdxib(:,1,3) * dxdxib(:,2,1)
             
      elseif (lcsyst .eq. 1) then
          temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
          temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
          temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
      else
          temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
          temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
          temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
      endif
      !
      temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      bnorm(:,1) = temp1 * temp
      bnorm(:,2) = temp2 * temp
      bnorm(:,3) = temp3 * temp
      !
     
      if (lcsyst .eq. 3) then
          WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
      elseif (lcsyst .eq. 4) then
          WdetJb     = Qwtb(lcsyst,intp) / temp
      else
          WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
      endif
      !
      !.... -------------------------->  Grad-V  <----------------------------
      !
      !.... compute grad-v for Navier-Stokes terms
      !
      if (Navier .eq. 1) then
          !
          !.... compute the inverse of deformation gradient
          !
          dxidxb(:,1,1) =   dxdxib(:,2,2) * dxdxib(:,3,3) &
          - dxdxib(:,3,2) * dxdxib(:,2,3)
          dxidxb(:,1,2) =   dxdxib(:,3,2) * dxdxib(:,1,3) &
          - dxdxib(:,1,2) * dxdxib(:,3,3)
          dxidxb(:,1,3) =   dxdxib(:,1,2) * dxdxib(:,2,3) &
          - dxdxib(:,1,3) * dxdxib(:,2,2)
          temp          = one / ( dxidxb(:,1,1) * dxdxib(:,1,1) &
          + dxidxb(:,1,2) * dxdxib(:,2,1) &
          + dxidxb(:,1,3) * dxdxib(:,3,1) )
          dxidxb(:,1,1) =  dxidxb(:,1,1) * temp
          dxidxb(:,1,2) =  dxidxb(:,1,2) * temp
          dxidxb(:,1,3) =  dxidxb(:,1,3) * temp
          dxidxb(:,2,1) = (dxdxib(:,2,3) * dxdxib(:,3,1) &
          - dxdxib(:,2,1) * dxdxib(:,3,3)) * temp
          dxidxb(:,2,2) = (dxdxib(:,1,1) * dxdxib(:,3,3) &
          - dxdxib(:,3,1) * dxdxib(:,1,3)) * temp
          dxidxb(:,2,3) = (dxdxib(:,2,1) * dxdxib(:,1,3) &
          - dxdxib(:,1,1) * dxdxib(:,2,3)) * temp
          dxidxb(:,3,1) = (dxdxib(:,2,1) * dxdxib(:,3,2) &
          - dxdxib(:,2,2) * dxdxib(:,3,1)) * temp
          dxidxb(:,3,2) = (dxdxib(:,3,1) * dxdxib(:,1,2) &
          - dxdxib(:,1,1) * dxdxib(:,3,2)) * temp
          dxidxb(:,3,3) = (dxdxib(:,1,1) * dxdxib(:,2,2) &
          - dxdxib(:,1,2) * dxdxib(:,2,1)) * temp
          !
          !.... compute local-grad-Y
          !
          !
          gradSl = zero
          isc=5+isclr
          do n = 1, nshl
              gradSl(:,1) = gradSl(:,1) + shdrv(:,1,n) * yl(:,n,isc)
              gradSl(:,2) = gradSl(:,2) + shdrv(:,2,n) * yl(:,n,isc)
              gradSl(:,3) = gradSl(:,3) + shdrv(:,3,n) * yl(:,n,isc)
          enddo
          !
          !.... convert local-grads to global-grads
          !
          gradS(:,1) = dxidxb(:,1,1) * gradSl(:,1) + &
          dxidxb(:,2,1) * gradSl(:,2) + &
          dxidxb(:,3,1) * gradSl(:,3)

          !
          gradS(:,2) = dxidxb(:,1,2) * gradSl(:,1) + &
          dxidxb(:,2,2) * gradSl(:,2) + &
          dxidxb(:,3,2) * gradSl(:,3)

          gradS(:,3) = dxidxb(:,1,3) * gradSl(:,1) + &
          dxidxb(:,2,3) * gradSl(:,2) + &
          dxidxb(:,3,3) * gradSl(:,3)
      !
      !.... end grad-T
      !
      endif

      flux = diffus * ( gradS(:,1) * bnorm(:,1) &
      + gradS(:,2) * bnorm(:,2) &
      + gradS(:,3) * bnorm(:,3) )
      !
      !.... return
      !
      return
  end
