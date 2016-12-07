      subroutine e3ql (yl,      dwl,     shp,     shgl, &
                       xl,      ql,      xmudmi, &
                       sgn )
!                                                                      
!----------------------------------------------------------------------
!
! This routine computes the local diffusive flux vector using a 
! local projection algorithm
!
! input: 
!  yl     (npro,nshl,ndof)       : Y variables
!  shp    (nen,ngauss)           : element shapeVar-functions
!  shgl   (nsd,nen,ngauss)       : element local-grad-shapeVar-functions
!  xl     (npro,nshape,nsd)      : nodal coordinates at current step
!  sgn    (npro,nshl)            : signs for reversed shapeVar functions
!  
! output:
!  ql     (npro,nshl,nsd*nsd) : element RHS diffusion residual 
!
!----------------------------------------------------------------------
!
      use local_mass
      use phcommonvars
      use nnw, only : get_shear_rate, get_mu

      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension yl(npro,nshl,ndof),        dwl(npro,nshl), &
                shp(nshl,ngauss),          shgl(nsd,nshl,ngauss), &
                xl(npro,nenl,nsd),         sgn(npro,nshl), &
                ql(npro,nshl,idflx), xmudmi(npro,ngauss)
!
! local arrays
!
      dimension g1yi(npro,ndof),           g2yi(npro,ndof), &
                g3yi(npro,ndof),           shg(npro,nshl,nsd), &
                dxidx(npro,nsd,nsd),       WdetJ(npro), &
                rmu(npro), &
                rminv(npro,nshl,nshl), &
                qrl(npro,nshl,nsd*nsd)
!
      dimension qdi(npro,nsd*nsd),    shapeVar(npro,nshl), &
                shdrv(npro,nsd,nshl),      indx(nshl), &
                rmass(npro,nshl,nshl)


        real*8 tmp(npro)
        real*8 :: gamma_shear(npro)
!
!.... loop through the integration points
!
      rminv = zero
      rmass = zero
      qrl   = zero
        
      do intp = 1, ngauss

         call getshp(shp, shgl, sgn, shapeVar, shdrv)

         qdi = zero
!
!.... calculate the integration variables 
!    
!
         call e3qvar   (yl,           shdrv,    &
                        xl,           g1yi, &
                        g2yi,         g3yi,         shg, &
                        dxidx,        WdetJ )

         call getdiff(dwl,  yl, shapeVar, xmudmi, xl,rmu, tmp)

         if (nnwType.eq.3) then  ! Recalculating rmu for non-Newtonian flow (We needed the gradients that are calcualted here)
            call get_shear_rate(gamma_shear,npro,g1yi(:,2:4),g2yi(:,2:4),g3yi(:,2:4))
            call get_mu(rmu,gamma_shear,npro)
         endif
!
!.... diffusive flux in x1-direction
!
         qdi(:,1) =  two * rmu *  g1yi(:,2)
         qdi(:,4) =        rmu * (g1yi(:,3) + g2yi(:,2))
         qdi(:,7) =        rmu * (g1yi(:,4) + g3yi(:,2))
!
!.... diffusive flux in x2-direction
!
         qdi(:,2) =        rmu * (g1yi(:,3) + g2yi(:,2))
         qdi(:,5) =  two * rmu *  g2yi(:,3)
         qdi(:,8) =        rmu * (g2yi(:,4) + g3yi(:,3))
!     
!.... diffusive flux in x3-direction
!
         qdi(:,3) =        rmu * (g1yi(:,4) + g3yi(:,2))
         qdi(:,6)=        rmu * (g2yi(:,4) + g3yi(:,3))
         qdi(:,9)=  two * rmu *  g3yi(:,4)
!
!
!.... assemble contribution of qdi to qrl,i.e., contribution to 
!     each element shapeVar function
!
         tmp = Qwt(lcsyst,intp)
         if (lcsyst .eq. 1) then 
            tmp = tmp*(three/four)
         endif
!
! reconsider this when hierarchic wedges come into code WDGCHECK
!
        
         do i=1,nshl
            qrl(:,i,1 ) = qrl(:,i,1 )+ shapeVar(:,i)*tmp*qdi(:,1 )
            qrl(:,i,2 ) = qrl(:,i,2 )+ shapeVar(:,i)*tmp*qdi(:,2 )
            qrl(:,i,3 ) = qrl(:,i,3 )+ shapeVar(:,i)*tmp*qdi(:,3 )
            
            qrl(:,i,4 ) = qrl(:,i,4 )+ shapeVar(:,i)*tmp*qdi(:,4 )
            qrl(:,i,5 ) = qrl(:,i,5 )+ shapeVar(:,i)*tmp*qdi(:,5 )
            qrl(:,i,6 ) = qrl(:,i,6 )+ shapeVar(:,i)*tmp*qdi(:,6 )

            qrl(:,i,7 ) = qrl(:,i,7 )+ shapeVar(:,i)*tmp*qdi(:,7 )
            qrl(:,i,8 ) = qrl(:,i,8 )+ shapeVar(:,i)*tmp*qdi(:,8 )
            qrl(:,i,9 ) = qrl(:,i,9 )+ shapeVar(:,i)*tmp*qdi(:,9 )
         enddo
!
!.... add contribution to local mass matrix
!

         if (have_local_mass .eq. 0) then
            do i=1,nshl
               do j=1,nshl
                  rmass(:,i,j) = rmass(:,i,j)+shapeVar(:,i)*shapeVar(:,j)*tmp
              enddo
           enddo
        endif
!
!.... end of the loop over integration points
!
      enddo

!
!.... find the inverse of the local mass matrix for each element


         if (have_local_mass .eq. 0) then
            allocate (lmassinv(iblock)%p(npro,nshl,nshl))

            do iel=1,npro
               do i=1,nshl      ! form the identy matrix
                  do j=1,nshl
                     lmassinv(iblock)%p(iel,i,j) = 0.0
                  enddo
                  lmassinv(iblock)%p(iel,i,i)=1.0
               enddo
!     
!.... LU factor the mass matrix
!
               call ludcmp(rmass(iel,:,:),nshl,nshl,indx,d)
!     
!.... back substitute with the identy matrix to find the
!     matrix inverse
!          
               do j=1,nshl
                  call lubksb(rmass(iel,:,:),nshl,nshl,indx, &
                              lmassinv(iblock)%p(iel,:,j))
               enddo
            enddo
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         else
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         endif
!
!.... find the modal coefficients of ql by multiplying by the inverse of
!     the local mass matrix
!
      do iel=1,npro
        do j=1,9
!         do j=1, 3*nsd
            ql(iel,:,j) = matmul( rminv(iel,:,:),qrl(iel,:,j) )
         enddo
      enddo
!
!.... return
!
      return
      end




      subroutine e3qlSclr (yl,      dwl,     shp,     shgl, &
                           xl,      ql,      sgn )
!                                                                      
!----------------------------------------------------------------------
!
! This routine computes the local diffusive flux vector using a 
! local projection algorithm: 
!     diffus * phi,i
!
!----------------------------------------------------------------------
!
      use local_mass
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension yl(npro,nshl,ndof),        dwl(npro,nshl), &
                shp(nshl,ngauss),          shgl(nsd,nshl,ngauss), &
                xl(npro,nenl,nsd),         sgn(npro,nshl), &
                ql(npro,nshl,nsd)
!
! local arrays
!
      dimension dxidx(npro,nsd,nsd),       WdetJ(npro), &
                diffus(npro), &
                rminv(npro,nshl,nshl), &
                qrl(npro,nshl,nsd)
!
      dimension qdi(npro,nsd),    shapeVar(npro,nshl), &
                shdrv(npro,nsd,nshl),      indx(nshl), &
                rmass(npro,nshl,nshl),     gradT(npro,nsd), &
                eviscv(npro)

!
!.... loop through the integration points
!
      rminv = zero
      rmass = zero
      qrl   = zero
        
      do intp = 1, ngauss

         call getshp(shp, shgl, sgn, shapeVar, shdrv)

         qdi = zero
!
!.... calculate the integration variables 
!    
!
         call e3qvarSclr  (yl,           shdrv,        xl,            &
                           gradT,        dxidx,        WdetJ )
!
!....  call function to sort out diffusivity (at end of this file)
!
         call getdiffsclr(dwl,shapeVar,yl, diffus)
!
!.... diffusive flux in x1-direction
!
         qdi(:,1) =  diffus * gradT(:,1)
         qdi(:,2) =  diffus * gradT(:,2)
         qdi(:,3) =  diffus * gradT(:,3)

!
!.... assemble contribution of qdi to qrl,i.e., contribution to 
!     each element shapeVar function
!
         tmp = Qwt(lcsyst,intp)
         if (lcsyst .eq. 1) then 
            tmp = tmp*(three/four)
         endif
        
         do i=1,nshl
            qrl(:,i,1 ) = qrl(:,i,1 )+ shapeVar(:,i)*tmp*qdi(:,1 )
            qrl(:,i,2 ) = qrl(:,i,2 )+ shapeVar(:,i)*tmp*qdi(:,2 )
            qrl(:,i,3 ) = qrl(:,i,3 )+ shapeVar(:,i)*tmp*qdi(:,3 )
         enddo
!
!.... add contribution to local mass matrix
!
         if (have_local_mass .eq. 0) then
            do i=1,nshl
               do j=1,nshl
                  rmass(:,i,j)=rmass(:,i,j)+shapeVar(:,i)*shapeVar(:,j)*tmp
               enddo
            enddo
         endif

!.... end of the loop over integration points
!
      enddo

!
!.... find the inverse of the local mass matrix for each element
!
       qrl   = qrl/6.d0
!
!.... Assuming that lmassinv was already computed for flow equations
!     
       rmass = rmass/6.0
!
!.... for cubics, it cannot be precomputed, so compute and
!     save it the first time it is needed
!
         if (have_local_mass .eq. 0) then
            allocate (lmassinv(iblock)%p(npro,nshl,nshl))

            do iel=1,npro
               do i=1,nshl      ! form the identy matrix
                  do j=1,nshl
                     lmassinv(iblock)%p(iel,i,j) = 0.0
                  enddo
                  lmassinv(iblock)%p(iel,i,i)=1.0
               enddo
!     
!.... LU factor the mass matrix
!
               call ludcmp(rmass(iel,:,:),nshl,nshl,indx,d)
!     
!.... back substitute with the identy matrix to find the
!     matrix inverse
!          
               do j=1,nshl
                  call lubksb(rmass(iel,:,:),nshl,nshl,indx, &
                              lmassinv(iblock)%p(iel,:,j))
               enddo
            enddo
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         else
            rminv(:,:,:) = lmassinv(iblock)%p(:,:,:)
         endif
!
!.... find the modal coefficients of ql by multiplying by the inverse of
!     the local mass matrix
!
      do iel=1,npro
         do j=1,nsd
            ql(iel,:,j) = matmul( rminv(iel,:,:),qrl(iel,:,j) )
         enddo
      enddo
!
!.... return
!
      return
      end
