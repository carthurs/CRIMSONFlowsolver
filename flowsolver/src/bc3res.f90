        subroutine bc3Res ( iBC,  BC,  res, iper, ilwork)
!
!----------------------------------------------------------------------
!
! This routine satisfies the BC of the residual vector for 3D elements.
!
! input:
!  iBC   (nshg)        : Boundary Condition Code
!  BC    (nshg,ndofBC) : the boundary condition constraint parameters
!  res   (nshg,nflow)   : residual before BC is applied
!
! output:
!  res   (nshg,nflow)   : residual after satisfaction of BC
!  
!
! Thuc Bui,      Winter 1989.
! Zdenek Johan,  Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension iBC(nshg), & 
                  BC(nshg,ndofBC),    & 
                  res(nshg,nflow),           ilwork(nlwork), & 
                  iper(nshg)
!
!.... local periodic boundary conditions (no communications)
!
        call bc3per(iBC,  res, iper, ilwork, nflow)
! 
!.... pressure 
!
        where (btest(iBC,2))
           res(:,4) = zero
        endwhere
!
!.... velocities
!
! ibits(n1,n2,n3) extracts bits n2+1 through n2+n3 (extending to the left
! as is traditional in binary) of the integer n1
! and returns the base 10 integer. In examples below x y z a b can 
! be 1 or zero without any effect.
!
!.... x1-velocity
!
! if iBC=4   bits of ibc =00000100 => ibits(4,3,3)=0
! if iBC=40  bits of ibc =00101000 => ibits(40,3,3)=5
! if iBC=40  bits of ibc =00101000 => ibits(40,3,2)=1
!
        where (ibits(iBC,3,3) .eq. 1)   ! bits of iBC= xy001zab 
!
!     notice that the extracted 3 bits form the number 1.  below
!     you will see the combinations which make up 2-7, all of the
!     possible velocity combinations
!
          res(:,2) = res(:,2) - BC(:,4) * res(:,1)
          res(:,3) = res(:,3) - BC(:,5) * res(:,1)
          res(:,1) = zero
        endwhere
!
!.... x2-velocity
!
        where (ibits(iBC,3,3) .eq. 2)   ! bits of iBC= xy010zab 
          res(:,1) = res(:,1) - BC(:,4) * res(:,2)
          res(:,3) = res(:,3) - BC(:,5) * res(:,2)
          res(:,2) = zero
        endwhere
!
!.... x1-velocity and x2-velocity
!
        where (ibits(iBC,3,3) .eq. 3)  ! bits of iBC= xy011zab 
          res(:,3) = res(:,3) - BC(:,4) * res(:,1) - BC(:,6) * res(:,2)
          res(:,1) = zero
          res(:,2) = zero
        endwhere
!
!.... x3-velocity
!
        where (ibits(iBC,3,3) .eq. 4)  ! bits of iBC= xy100zab 
          res(:,1) = res(:,1) - BC(:,4) * res(:,3)
          res(:,2) = res(:,2) - BC(:,5) * res(:,3)
          res(:,3) = zero
        endwhere
!
!.... x1-velocity and x3-velocity
!
        where (ibits(iBC,3,3) .eq. 5)  ! bits of iBC= xy101zab 
          res(:,2) = res(:,2) - BC(:,4) * res(:,1) - BC(:,6) * res(:,3)
          res(:,1) = zero
          res(:,3) = zero
        endwhere
!
!.... x2-velocity and x3-velocity
!
        where (ibits(iBC,3,3) .eq. 6)  ! bits of iBC= xy110zab 
          res(:,1) = res(:,1) - BC(:,4) * res(:,2) - BC(:,6) * res(:,3)
          res(:,2) = zero
          res(:,3) = zero
        endwhere
!
!.... x1-velocity, x2-velocity and x3-velocity
!
        where (ibits(iBC,3,3) .eq. 7)  ! bits of iBC= xy111zab 
          res(:,1) = zero
          res(:,2) = zero
          res(:,3) = zero
        endwhere
!
!.... scaled plane extraction boundary condition
!
        if(intpres.eq.1) then  ! interpolating pressure so zero continuity res 
           where (btest(iBC,11))
              res(:,1) = zero
              res(:,2) = zero
              res(:,3) = zero
              res(:,4) = zero
           endwhere
        else  ! leave residual in continuity equation
           where (btest(iBC,11))
              res(:,1) = zero
              res(:,2) = zero
              res(:,3) = zero
           endwhere
        endif
!
!.... return
!
        return
        end


!---------------------------------------------------------------------
!
!     boundary conditions on scalar residual
!
!---------------------------------------------------------------------
        subroutine bc3ResSclr (iBC,  res, iper, ilwork)

  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension iBC(nshg), & 
                  res(nshg),                ilwork(nlwork), & 
                  iper(nshg)


        if(isclr.eq.0) then
!     
!.... temperature
!     
           where (btest(iBC,1)) res(:) = zero
        else
!
!.... turbulence or scalar
!
           is=isclr+5
           where (btest(iBC,is)) res(:) = zero
        endif
!
!.... local periodic boundary conditions (no communications)
!
        call bc3per(iBC,  res, iper, ilwork, 1)
!
!.... return
!
        return
        end

