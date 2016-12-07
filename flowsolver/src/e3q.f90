        subroutine e3q (yl,      dwl,     shp,     shgl, &
                        xl,      ql,      rmassl,  &
                        xmudmi,  sgn )
!                                                                      
!----------------------------------------------------------------------
!
! This routine computes the element contribution to the 
! diffusive flux vector and the lumped mass matrix.
!
! input: 
!  yl     (npro,nshl,ndof)       : Y variables
!  shp    (nen,ngauss)            : element shape-functions
!  shgl   (nsd,nen,ngauss)        : element local-grad-shape-functions
!  xl     (npro,nshl,nsd)        : nodal coordinates at current step
!  
! output:
!  ql     (npro,nshl,idflx) : element RHS diffusion residual 
!  rmassl     (npro,nshl)        : element lumped mass matrix
!
!----------------------------------------------------------------------
!
        use phcommonvars
        use nnw, only : get_shear_rate, get_mu

        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension yl(npro,nshl,ndof),     dwl(npro,nenl), &
                  shp(nshl,ngauss),      shgl(nsd,nshl,ngauss), &
                  xl(npro,nenl,nsd), &
                  ql(npro,nshl,idflx),  rmassl(npro,nshl), &
                  xmudmi(npro,ngauss)
!
! local arrays
!
        dimension g1yi(npro,nflow),           g2yi(npro,nflow), &
                  g3yi(npro,nflow),           shg(npro,nshl,nsd), &
                  dxidx(npro,nsd,nsd),       WdetJ(npro), &
                  rmu(npro) 
!
        dimension qdi(npro,idflx),alph1(npro),alph2(npro)
!
        dimension sgn(npro,nshl),          shapeVar(npro,nshl), &
                  shdrv(npro,nsd,nshl),    shpsum(npro)

        real*8 tmp(npro)
        real*8 :: gamma_shear(npro)
!
!.... for surface tension
!     
        dimension g1yti(npro),          g2yti(npro), &
                  g3yti(npro)
        integer idflow
!
!.... loop through the integration points
!
        
        
        alph1 = 0.d0
        alph2 = 0.d0
        
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
!     
        call getshp(shp,          shgl,      sgn,  &
                    shapeVar,        shdrv)
        
!
!.... initialize
!
        qdi = zero
!
!
!.... calculate the integration variables necessary for the
!     formation of q
!

        call e3qvar   (yl,        shdrv,    &
                       xl,           g1yi, &
                       g2yi,      g3yi,         shg, &
                       dxidx,     WdetJ )      
!  
        idflow = 9   ! we ALWAYS save space for tau_{ij} in q_i 
                     ! even if idiff is not greater than 1

        if(idiff >= 1) then   !so taking care of all the idiff=1,3
!
!.... compute diffusive fluxes 
!
!.... compute the viscosity
!
        call getdiff(dwl, yl, shapeVar, xmudmi, xl, rmu, tmp)

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
        qdi(:,6)=         rmu * (g2yi(:,4) + g3yi(:,3))
        qdi(:,9)=  two * rmu *  g3yi(:,4)
!
!
!.... assemble contribution of qdi to ql,i.e., contribution to 
!     each element node
!     
        do i=1,nshl
           ql(:,i,1 ) = ql(:,i,1 )+ shapeVar(:,i)*WdetJ*qdi(:,1 )
           ql(:,i,2 ) = ql(:,i,2 )+ shapeVar(:,i)*WdetJ*qdi(:,2 )
           ql(:,i,3 ) = ql(:,i,3 )+ shapeVar(:,i)*WdetJ*qdi(:,3 )

           ql(:,i,4 ) = ql(:,i,4 )+ shapeVar(:,i)*WdetJ*qdi(:,4 )
           ql(:,i,5 ) = ql(:,i,5 )+ shapeVar(:,i)*WdetJ*qdi(:,5 )
           ql(:,i,6 ) = ql(:,i,6 )+ shapeVar(:,i)*WdetJ*qdi(:,6 )

           ql(:,i,7 ) = ql(:,i,7 )+ shapeVar(:,i)*WdetJ*qdi(:,7 )
           ql(:,i,8 ) = ql(:,i,8 )+ shapeVar(:,i)*WdetJ*qdi(:,8 )
           ql(:,i,9 ) = ql(:,i,9 )+ shapeVar(:,i)*WdetJ*qdi(:,9 )

        enddo
!
!.... compute and assemble the element contribution to the lumped
!     mass matrix
!
!
!.... row sum technique
!
        if ( idiff == 1 ) then
           do i=1,nshl
              rmassl(:,i) = rmassl(:,i) + shapeVar(:,i)*WdetJ
           enddo
        endif
!
!.... "special lumping technique" (Hughes p. 445)
!
        if ( idiff == 3 ) then
           shpsum = zero
           do i=1,nshl
              shpsum = shpsum + shapeVar(:,i)*shapeVar(:,i)
              rmassl(:,i)=rmassl(:,i)+shapeVar(:,i)*shapeVar(:,i)*WdetJ
           enddo
           alph1 = alph1+WdetJ
           alph2 = alph2+shpsum*WdetJ
        endif
      endif                     ! end of idiff=1 .or. 3 
!
      if(isurf .eq. 1) then
!
!.... initialize
!
        g1yti   = zero
        g2yti   = zero
        g3yti   = zero
!
!.... calculate the integration variables necessary for the
!     formation of q
!
!.... compute the global gradient of Yt-variables, assuming 6th entry as 
!.... the phase indicator function 
!
!  Yt_{,x_i}=SUM_{a=1}^nshl (N_{a,x_i}(int) Yta)
!
        do n = 1, nshl
          g1yti(:)  = g1yti(:)  + shg(:,n,1) * yl(:,n,6)
          g2yti(:)  = g2yti(:)  + shg(:,n,2) * yl(:,n,6)
          g3yti(:)  = g3yti(:)  + shg(:,n,3) * yl(:,n,6)
        enddo
!
!    computing N_{b}*N_{a,x_i)*yta*WdetJ
!
        do i=1,nshl
           ql(:,i,idflow+1)  = ql(:,i,idflow+1)   &
                             + shapeVar(:,i)*WdetJ*g1yti
           ql(:,i,idflow+2)  = ql(:,i,idflow+2)   &
                             + shapeVar(:,i)*WdetJ*g2yti
           ql(:,i,idflow+3)  = ql(:,i,idflow+3)   &
                             + shapeVar(:,i)*WdetJ*g3yti
           rmassl(:,i) = rmassl(:,i) + shapeVar(:,i)*WdetJ
        enddo
      endif  !end of the isurf  
!
!.... end of the loop over integration points
!
      enddo
!
!.... normalize the mass matrix for idiff == 3
!
      if ( idiff == 3 ) then
         do i=1,nshl
            rmassl(:,i) = rmassl(:,i)*alph1/alph2
         enddo
      endif
      

!
!.... return
!
       return
       end


        subroutine e3qSclr (yl,      dwl,     shp,     shgl, &
                            xl,      ql,      rmassl,  &
                            sgn )
!                                                                      
!----------------------------------------------------------------------
!
! This routine computes the element contribution to the 
! diffusive flux vector and the lumped mass matrix.
!
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension yl(npro,nshl,ndof),    dwl(npro,nshl), &
                  shp(nshl,ngauss),      shgl(nsd,nshl,ngauss), &
                  xl(npro,nenl,nsd), &
                  ql(npro,nshl,nsd),     rmassl(npro,nshl)
!
! local arrays
!
        dimension gradT(npro,nsd),      &
                  dxidx(npro,nsd,nsd),       WdetJ(npro)
!
        dimension qdi(npro,nsd),alph1(npro),alph2(npro)
!
        dimension sgn(npro,nshl),          shapeVar(npro,nshl), &
                  shdrv(npro,nsd,nshl),    shpsum(npro)

        real*8 diffus(npro)
!
!.... loop through the integration points
!
        
        
        alph1 = 0.d0
        alph2 = 0.d0
        
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
!     
        call getshp(shp,          shgl,      sgn,  &
                    shapeVar,        shdrv)
        
!
!.... initialize
!
        qdi = zero
!
!
!.... calculate the integration variables necessary for the
!     formation of q 
!
        call e3qvarSclr   (yl,        shdrv,    &
                           xl,        gradT, &
                           dxidx,     WdetJ )        

!
!.... compute diffusive flux vector at this integration point
!
        call getdiffsclr(shapeVar, dwl, yl, diffus)

!
!.... diffusive flux 
!
        qdi(:,1) =  diffus * gradT(:,1)
        qdi(:,2) =  diffus * gradT(:,2)
        qdi(:,3) =  diffus * gradT(:,3)
!
!
!.... assemble contribution of qdi to ql,i.e., contribution to 
!     each element node
!     
        do i=1,nshl
           ql(:,i,1 ) = ql(:,i,1 )+ shapeVar(:,i)*WdetJ*qdi(:,1 )
           ql(:,i,2 ) = ql(:,i,2 )+ shapeVar(:,i)*WdetJ*qdi(:,2 )
           ql(:,i,3 ) = ql(:,i,3 )+ shapeVar(:,i)*WdetJ*qdi(:,3 )

        enddo
!
!.... compute and assemble the element contribution to the lumped
!     mass matrix
!
!
!.... row sum technique
!
        if ( idiff == 1 ) then
           do i=1,nshl
              rmassl(:,i) = rmassl(:,i) + shapeVar(:,i)*WdetJ
           enddo
        endif
!
!.... "special lumping technique" (Hughes p. 445)
!
        if ( idiff == 3 ) then
           shpsum = zero
           do i=1,nshl
              shpsum = shpsum + shapeVar(:,i)*shapeVar(:,i)
              rmassl(:,i)=rmassl(:,i)+shapeVar(:,i)*shapeVar(:,i)*WdetJ
           enddo
           alph1 = alph1+WdetJ
           alph2 = alph2+shpsum*WdetJ
        endif
!
!.... end of the loop over integration points
!
      enddo
!
!.... normalize the mass matrix for idiff == 3
!
      if ( idiff == 3 ) then
         do i=1,nshl
            rmassl(:,i) = rmassl(:,i)*alph1/alph2
         enddo
      endif
!
!.... return
!
       return
       end

