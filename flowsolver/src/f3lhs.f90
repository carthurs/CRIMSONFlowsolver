      subroutine f3lhs (shpb,   shglb,  xlb,    flhsl, &
                        fnrml,  sgn )
!
!----------------------------   ------------------------------------------
!
!  This subroutine computes the element LHS matrix and the normal 
! to the boundary for computation of (output) boundary fluxes. 
! 
! input:
!  shpb   (nen,nintg)           : boundary element shape-functions
!  shglb  (nsd,nen,nintg)       : boundary element grad-shape-functions
!  wghtb  (nintg)               : boundary element weight
!  xlb    (npro,nenl,nsd)       : nodal coordinates
!  sgn    (npro,nshl)           : mode signs for hierarchic basis
!
! output:
!  flhsl  (npro,nenl,1)         : element lumped lhs on flux boundary
!  fnrml  (npro,nenl,nsd)       : RHS of LS projection of normal to 
!                                  flux boundary
!
!
! Note: Special lumping technique is used to compute the LHS. 
!       See T.J.R. Hughes, "The Finite Element Method: Linear 
!       Static and Dynamic Finite Element Analysis", page 445.  
!
! Note: Least-squares projection is used to compute the normal to
!       the boundary at the nodes.  This routine provides the element
!       contribution to the RHS of the projection linear system.
!
!
! Zdenek Johan, Summer 1991.
!----------------------------------------------------------------------
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension shpb(nshl,ngaussb),        shglb(nsd,nshl,ngaussb), &
                xlb(npro,nenl,nsd), &
                flhsl(npro,nshl,1),        fnrml(npro,nshl,nsd)
!
      dimension WdetJb(npro), &
                bnorm(npro,nsd),           fmstot(npro), &
                temp(npro),                temp1(npro), &
                temp2(npro),               temp3(npro)

      dimension sgn(npro,nshl),            shapeVar(npro,nshl), &
                shdrv(npro,nsd,nshl), &
                v1(npro,nsd),              v2(npro,nsd)
!
!.... integrate the lumped LHS matrix and normal
!
      fmstot = zero
!
!
!.... compute the normal to the boundary 
!

      v1 = xlb(:,2,:) - xlb(:,1,:)
      v2 = xlb(:,3,:) - xlb(:,1,:)
      
      if (lcsyst .eq. 1) then
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
      
      do intp = 1, ngaussb
!
!.... get the hierarchic shape functions at this int point
!
         call getshp(shpb,        shglb,        sgn,  &
                     shapeVar,       shdrv)
!
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
!
!.... compute the lumped LHS and normal
!          
         do n = 1, nenl ! when changed to nshl solution degraded ipord 10
            flhsl(:,n,1) = flhsl(:,n,1) + WdetJb * shapeVar(:,n)

! for curved geometries the below construct for the normals has to be used
            fnrml(:,n,1) = fnrml(:,n,1) + WdetJb * bnorm(:,1) &
                                                 * shapeVar(:,n)
            fnrml(:,n,2) = fnrml(:,n,2) + WdetJb * bnorm(:,2) &
                                                 * shapeVar(:,n)
            fnrml(:,n,3) = fnrml(:,n,3) + WdetJb * bnorm(:,3) &
                                                 * shapeVar(:,n)
          enddo
!
!  To best represent this case it should be assigned to the vertex 
!  modes and higher entities should get zero as is done below
!
          fmstot = fmstot + WdetJb
!
        enddo
        
!$$$        do i=1,nenl
!$$$           fnrml(:,i,:)=bnorm(:,:)
!$$$        enddo
        if(ipord.gt.1)  fnrml(:,nenl:nshl,:)=zero
!
!.... scale the LHS matrix contribution
!
        temp = zero
        do n = 1, nshl
           temp = temp + flhsl(:,n,1)
        enddo
!
        do n = 1, nshl
           flhsl(:,n,1) = flhsl(:,n,1) * fmstot / temp
        enddo
!
!.... return
!
        return
        end
