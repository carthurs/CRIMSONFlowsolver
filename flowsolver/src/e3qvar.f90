        subroutine e3qvar (yl,          shgl,     &
                           xl,          g1yi, &
                           g2yi,        g3yi,        shg, &
                           dxidx,       WdetJ )
!
!----------------------------------------------------------------------
!
!  This routine computes the variables at integration point
!  necessary for the computation of the diffusive flux vector.
!
! input:
!  yl     (npro,nshl,ndof)      : primitive variables
!  shgl   (npro,nsd,nshl)     : element local-grad-shape-functions
!  xl     (npro,nenl,nsd)       : nodal coordinates at current step
!
! output:
!  g1yi   (npro,ndof)           : grad-y in direction 1
!  g2yi   (npro,ndof)           : grad-y in direction 2
!  g3yi   (npro,ndof)           : grad-y in direction 3
!  shg    (npro,nshl,nsd)       : element global grad-shape-functions
!  dxidx  (npro,nsd,nsd)        : inverse of deformation gradient
!  WdetJ  (npro)                : weighted Jacobian
!  u1     (npro)                : x1-velocity component
!  u2     (npro)                : x2-velocity component
!  u3     (npro)                : x3-velocity component
!
!
! Zdenek Johan, Summer 1990. (Modified from e2ivar.f)
! Zdenek Johan, Winter 1991. (Fortran 90)
! Kenneth Jansen, Winter 1997. Primitive Variables
!----------------------------------------------------------------------
!
         use phcommonvars
         IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
!  passed arrays
!
        dimension yl(npro,nshl,ndof),  &
                  shgl(npro,nsd,nshl), xl(npro,nenl,nsd), &
                  g1yi(npro,nflow),       g2yi(npro,nflow), &
                  g3yi(npro,nflow),       shg(npro,nshl,nsd),  &
                  dxidx(npro,nsd,nsd),   WdetJ(npro)
!
!  local arrays
!
        dimension tmp(npro),           dxdxi(npro,nsd,nsd)

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
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)  &
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
        WdetJ = Qwt(lcsyst,intp)/ tmp

!
!.... --------------------->  Global Gradients  <-----------------------
!
        g1yi = zero
        g2yi = zero
        g3yi = zero
!
!
        do n = 1, nshl
!
!.... compute the global gradient of shape-function
!
!            ! N_{a,x_i}= N_{a,xi_i} xi_{i,x_j}
!
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) +  &
                       shgl(:,2,n) * dxidx(:,2,1) + &
                       shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) +  &
                       shgl(:,2,n) * dxidx(:,2,2) + &
                       shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) +  &
                       shgl(:,2,n) * dxidx(:,2,3) + &
                       shgl(:,3,n) * dxidx(:,3,3) 
!
!.... compute the global gradient of Y-variables
!
!
!  Y_{,x_i}=SUM_{a=1}^nenl (N_{a,x_i}(int) Ya)
!
          g1yi(:,2) = g1yi(:,2) + shg(:,n,1) * yl(:,n,2)
          g1yi(:,3) = g1yi(:,3) + shg(:,n,1) * yl(:,n,3)
          g1yi(:,4) = g1yi(:,4) + shg(:,n,1) * yl(:,n,4)
!
          g2yi(:,2) = g2yi(:,2) + shg(:,n,2) * yl(:,n,2)
          g2yi(:,3) = g2yi(:,3) + shg(:,n,2) * yl(:,n,3)
          g2yi(:,4) = g2yi(:,4) + shg(:,n,2) * yl(:,n,4)
!
          g3yi(:,2) = g3yi(:,2) + shg(:,n,3) * yl(:,n,2)
          g3yi(:,3) = g3yi(:,3) + shg(:,n,3) * yl(:,n,3)
          g3yi(:,4) = g3yi(:,4) + shg(:,n,3) * yl(:,n,4)

       enddo

!
!.... return
!

       return
       end

!-----------------------------------------------------------------------
!
!     compute the variables for the local scalar diffusion
!
!-----------------------------------------------------------------------
      subroutine e3qvarSclr  (yl,       shgl,         xl,  &
                              gradT,    dxidx,        WdetJ )

       use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
!  passed arrays
!
      real*8   yl(npro,nshl,ndof),    shp(npro,nshl), &
               shgl(npro,nsd,nshl),   xl(npro,nenl,nsd), &
               dxidx(npro,nsd,nsd),   WdetJ(npro), &
               gradT(npro,nsd)
!
!  local arrays
!
      real*8   shg(npro,nshl,nsd)


      call e3metric( xl,         shgl,       dxidx,   &
                     shg,        WdetJ)

      gradT = zero
      id=5+isclr
!
!  later, when there are more models than SA we will need a 
!  more general function to calculate evisc at a quadrature point
!
      do n = 1, nshl
         gradT(:,1) = gradT(:,1) + shg(:,n,1) * yl(:,n,id)
         gradT(:,2) = gradT(:,2) + shg(:,n,2) * yl(:,n,id)
         gradT(:,3) = gradT(:,3) + shg(:,n,3) * yl(:,n,id)
      enddo
!
!.... return
!

       return
       end
       
            


