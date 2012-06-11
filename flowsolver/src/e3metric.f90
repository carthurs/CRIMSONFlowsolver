!-----------------------------------------------------------------------
!
!  compute the metrics of the mapping from global to local 
!  coordinates and the jacobian of the mapping (weighted by 
!  the quadrature weight
!
!-----------------------------------------------------------------------
      subroutine e3metric(  xl,      shgl,     dxidx, &
                            shg,     WdetJ)

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8     xl(npro,nenl,nsd),    shgl(npro,nsd,nshl), &
                 dxidx(npro,nsd,nsd),  shg(npro,nshl,nsd),  &
                 WdetJ(npro)

      real*8     dxdxi(npro,nsd,nsd),  tmp(npro)

!
!.... compute the deformation gradient
!
      dxdxi = zero
!
       do n = 1, nenl
          dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(:,1,n)
          dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(:,2,n)
          dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(:,3,n)
          dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(:,1,n)
          dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(:,2,n)
          dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(:,3,n)
          dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(:,1,n)
          dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(:,2,n)
          dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(:,3,n)
       enddo
!
!.... compute the inverse of deformation gradient
!
       dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)  &
                      - dxdxi(:,3,2) * dxdxi(:,2,3)
       dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)  &
                      - dxdxi(:,1,2) * dxdxi(:,3,3)
       dxidx(:,1,3) =  dxdxi(:,1,2) * dxdxi(:,2,3)  &
                      - dxdxi(:,1,3) * dxdxi(:,2,2)
       tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)  &
                             + dxidx(:,1,2) * dxdxi(:,2,1)   &
                             + dxidx(:,1,3) * dxdxi(:,3,1) )
       dxidx(:,1,1) = dxidx(:,1,1) * tmp
       dxidx(:,1,2) = dxidx(:,1,2) * tmp
       dxidx(:,1,3) = dxidx(:,1,3) * tmp
       dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)  &
                      - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
       dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)  &
                      - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
       dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)  &
                      - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
       dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)  &
                      - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
       dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)  &
                      - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
       dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)  &
                      - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
!
       WdetJ = Qwt(lcsyst,intp) / tmp
!
!.... compute the global gradient of shape-functions
!
       do n = 1, nshl
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) +  &
                       shgl(:,2,n) * dxidx(:,2,1) + &
                       shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) +  &
                       shgl(:,2,n) * dxidx(:,2,2) + &
                       shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) +  &
                       shgl(:,2,n) * dxidx(:,2,3) + &
                       shgl(:,3,n) * dxidx(:,3,3) 
       enddo

       return
       end



