        subroutine getFld (T,      cp,     rmu,    rlm,    rlm2mu, &
                           con)
!
!----------------------------------------------------------------------
!
! This routine calculates the fluid material properties.
!
! input:
!  T      (npro)        : temperature
!  cp     (npro)        : specific heat at constant pressure
!
! output:
!  rmu    (npro)        : Mu
!  rlm    (npro)        : Lambda
!  rlm2mu (npro)        : Lambda + 2 Mu
!  con    (npro)        : Conductivity
!
! Note: material type flags
!         matflg(2):
!          eq. 0, constant viscosity
!          eq. 1, generalized Sutherland viscosity
!         matflg(3):
!          eq. 0, Stokes approximation
!          eq. 1, shear proportional bulk viscosity
!
!
! Farzin Shakib, Winter 1987.
! Zdenek Johan,  Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension T(npro),                   cp(npro), &
                  rmu(npro),                 rlm(npro), &
                  rlm2mu(npro),              con(npro)
!
!
!.... constant viscosity
!
        if (matflg(2,1) .eq. 0) then
!
          rmu = datmat(1,2,1)
!
        else
!
!.... generalized Sutherland viscosity
!
          rmu = datmat(1,2,1) * (T/datmat(2,2,1))*sqrt(T/datmat(2,2,1)) &
              * ( datmat(2,2,1) + datmat(3,2,1) ) / (T + datmat(3,2,1))
!
        endif
!
!.... calculate the second viscosity coefficient
!
        if (matflg(3,1) .eq. 0) then
!
          rlm = -pt66 * rmu
!
        else
!
          rlm = (datmat(1,3,1) - pt66) * rmu
!
        endif
!
!.... calculate the remaining quantities
!
        cp     = datmat(1,3,1)
        rlm2mu = rlm + two * rmu
        con    = datmat(1,4,1)
!
!.... return
!
        return
        end
