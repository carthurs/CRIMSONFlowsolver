  subroutine e3bvar (yl,      acl,     ul, &
                     shpb,    shglb, &
                     xlb,     xdistl,  xdnvl, &
                     lnode,   SWB,     TWB,  &
                     WdetJb,  bnorm,   pres, &
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
      use turbsa
      use pointer_data
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension   yl(npro,nshl,ndof),        rmu(npro), &
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
                  SWB(npro,nProps),          TWB(npro,2)
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
      real*8      strainterm(npro,5,nshlb)
      real*8      uterm(npro,nsd)
      real*8      ulrot(npro,nshl,nsd)
      real*8      xlbrot(npro,nenl,nsd)
      real*8      Dmatrix(npro,5,5)
      real*8      Dtimesstrain(npro,5,nshlb)
      real*8      dxdxib_el(npro,nsd,nsd)
      real*8      dxidxb_el(npro,nsd,nsd)
      real*8      shgb_el(npro,nsd,nshl)
      real*8      LHSwall(npro,nsd**2,nshl,nshl)
      real*8      straindisp_g(npro,5,nshlb,nsd)
      real*8      DtimesB(npro,5,nshlb,nsd)

      real*8      Wall_LHSfactor(npro), &
                  Wall_LHSfactorSupp(npro), &
                  Wall_LHSfactorDamp(npro), &
                  tempSuppStiff(npro),       tempSuppDamp(npro)

      dimension   tmp1(npro)
!     
      real*8      Turb(npro),                xki, &
                  xki3,                      fv1
!        
      integer   e, i, j
!      
      integer   aa, b, iblk


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

          ! ---------------------------------------------------------------------
          !.... -----> Deformable Wall Residual Terms  <-----------
          !

          ! define the material properties for the enhanced membrane
          ! the membrane thickness is always SWB(:,1)
          Dmatrix = zero
          if (iUseSWB.gt.0) then
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

          else

              ! when the legacy SWB field is not used
              ! only the isotropic case
              ! is implemented so far
              Dmatrix(:,1,1) = one
              Dmatrix(:,2,2) = Dmatrix(:,1,1)
              Dmatrix(:,1,2) = rnuvw
              Dmatrix(:,2,1) = Dmatrix(:,1,2)
              Dmatrix(:,3,3) = pt5*(1-rnuvw)
              Dmatrix(:,4,4) = pt5*(1-rnuvw)*rshearconstantvw
              Dmatrix(:,5,5) = pt5*(1-rnuvw)*rshearconstantvw

              Dmatrix(:,:,:) = evw * one/(one-rnuvw**2)*Dmatrix(:,:,:);
          endif


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
                      strainterm(:,j,nodlcl) = strainterm(:,j,nodlcl) + &
                                               straindisp_g(:,j,nodlcl,i) * ul(:,nodlcl,i)
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


          enddo

          rlKwall = zero
          LHSwall = zero

          do n = 1, nshlb

              nodlcln = lnode(n)

              do m = 1, nshlb

                  nodlclm = lnode(m)

                  ! multiply by transpose of global strain displacement matrix
                  do j = 1, 5
                      do i = 1, 3
                          rlKwall(:,nodlcln,i) = &
                          rlKwall(:,nodlcln,i) + &
                          straindisp_g(:,j,nodlcln,i) * Dtimesstrain(:,j,nodlclm)
                      enddo
                  enddo

              enddo

              ! ---------------------------------------------------------------------
              ! now add the prestress contribution
              !
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

              ! ---------------------------------------------------------------------
              ! complete the numerical integration
              ! multiply by weighted jacobian and thickness
              !
              do i = 1, 3
                  rlKwall(:,nodlcln,i) = rlKwall(:,nodlcln,i) * WdetJb * SWB(:,1)
              enddo


              ! ---------------------------------------------------------------------
              ! now compute the tangent LHS contribution
              !
              do m = 1, nshlb

                  nodlclm = lnode(m)

                  do k = 1, 3
                      do j = 1, 5
                          do i = 1, 3

                              LHSwall(:,3*(i-1)+k,nodlcln,nodlclm) = &
                              LHSwall(:,3*(i-1)+k,nodlcln,nodlclm) + &
                              straindisp_g(:,j,nodlcln,i) * DtimesB(:,j,nodlclm,k) * &
                              WdetJb * SWB(:,1)

                          enddo
                      enddo
                  enddo
              enddo


          enddo

          ! ---------------------------------------------------------------------
          ! mass term
          !
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
          if (iwallsupp.eq.1) then
         
              usup1 = zero
         
              do n = 1, nshlb
         
                  nodlcl = lnode(n)
               
                  usup1(:,1) = usup1(:,1) + shpb(:,nodlcl) * ul(:,nodlcl,1)
                  usup1(:,2) = usup1(:,2) + shpb(:,nodlcl) * ul(:,nodlcl,2)
                  usup1(:,3) = usup1(:,3) + shpb(:,nodlcl) * ul(:,nodlcl,3)
            
              enddo
         
              usup1(:,1) = usup1(:,1) * TWB(:,1)
              usup1(:,2) = usup1(:,2) * TWB(:,1)
              usup1(:,3) = usup1(:,3) * TWB(:,1)
         
              Wall_LHSfactor = Wall_LHSfactor + &
              Delt(itseq) * alfi * betai * Delt(itseq) * TWB(:,1)
            
          endif
         
          if (idistancenudge.eq.1) then
         
              usup2 = zero
         
              do n = 1, nshlb
         
                  nodlcl = lnode(n)

                  usup2(:,1) = usup2(:,1) + &
                  shpb(:,nodlcl) * xdistl(:,nodlcl) * xdnvl(:,nodlcl,1)

                  usup2(:,2) = usup2(:,2) + &
                  shpb(:,nodlcl) * xdistl(:,nodlcl) * xdnvl(:,nodlcl,2)

                  usup2(:,3) = usup2(:,3) + &
                  shpb(:,nodlcl) * xdistl(:,nodlcl) * xdnvl(:,nodlcl,3)
     
              enddo
            
          endif
            
          ! tissue support damping term
          if (iwalldamp.eq.1) then
         
              velsup = zero
         
              do n = 1, nshlb
         
                  nodlcl = lnode(n)
                  velsup(:,1) = velsup(:,1) + shpb(:,nodlcl) * yl(:,nodlcl,2)
                  velsup(:,2) = velsup(:,2) + shpb(:,nodlcl) * yl(:,nodlcl,3)
                  velsup(:,3) = velsup(:,3) + shpb(:,nodlcl) * yl(:,nodlcl,4)
            
              enddo

              velsup(:,1) = velsup(:,1) * TWB(:,2)
              velsup(:,2) = velsup(:,2) * TWB(:,2)
              velsup(:,3) = velsup(:,3) * TWB(:,2)

              Wall_LHSfactor = Wall_LHSfactor + Delt(itseq) * alfi * gami * TWB(:,2)
            
          end if
                   
          ! ---------------------------------------------------------------------
          !.... -----> Wall Mass matrix for implicit LHS  <-----------
          !
          !           alpha_m*(1-lmp)*WdetJ*N^aN^b*rho*thickness         mass term
          !          +alpha_fv*gamma_v*DeltaT*WdetJ*N^a*N^b*supportvisc  viscous damping
          !          +alpha_fv*beta*DeltaT^2*WdetJ*N^a*N^b*supportstiff  elastic support
          !
         
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

          ! ---------------------------------------------------------------------
          ! .... -----> Wall Stiffness matrix for implicit LHS
          !

          xKebe = xKebe + LHSwall*iwallstiffactor*betai*Delt(itseq)*Delt(itseq)*alfi

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
