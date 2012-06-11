        subroutine shp4T (QuadPt, nQpt,   shp,    shgl,   wght)
!
!----------------------------------------------------------------------
!
!  This subroutine generates shape functions for the 4-node
!   tetrahedra.
!
! input:
!  QuadPt (4,nQpt)              : quadrature points' local coord's
!                                   QuadPt(1,*) : r
!                                   QuadPt(2,*) : s
!                                   QuadPt(3,*) : t
!                                   QuadPt(4,*) : wght
!  nQpt                         : number of quadrature points
!
! output:
!  shp    (nen,nQpt)            : shape functions
!  shgl   (nsd,nen,nQpt)        : local-gradient of shape function 
!  wght   (nQpt)                : quadrature weights
!
!
! shape-functions:
!  N1 = 1 - r - s - t
!  N2 = r
!  N3 = s
!  N4 = t
!
! Note: To be compatible with design of Tau and DC, the local 
!       gradients are divided by 2.  This is equivalent to having
!       r=[-1,1], s=[-1,1] and t=[-1,1], without really changing
!       r, s and t points (i.e., Qpt is for r=[0,1], s=[0,1] and
!       t=[0,1] range)
!
! Zdenek Johan, Summer 1990.
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension QuadPt(4,*),                   shp(nen,*), &
                  shgl(nsd,nen,*),               wght(*)
!
!.... loop through quadrature points
!
        do m = 1, nQpt
!
!.... generate the local-shape-functions
!
          shp(1,m) = one - QuadPt(1,m) - QuadPt(2,m) - QuadPt(3,m)
          shp(2,m) = QuadPt(1,m)
          shp(3,m) = QuadPt(2,m)
          shp(4,m) = QuadPt(3,m)
!
!.... generate the grad-local-shape-functions
!
          shgl(1,1,m) = -pt5
          shgl(2,1,m) = -pt5
          shgl(3,1,m) = -pt5
          shgl(1,2,m) =  pt5
          shgl(2,2,m) =  zero
          shgl(3,2,m) =  zero
          shgl(1,3,m) =  zero
          shgl(2,3,m) =  pt5
          shgl(3,3,m) =  zero
          shgl(1,4,m) =  zero
          shgl(2,4,m) =  zero
          shgl(3,4,m) =  pt5
!
!.... copy the weight
!
          wght(m) = QuadPt(4,m)
!
!.... end of shape-function loop
!
        enddo
!
!.... return
!
        return
        end
