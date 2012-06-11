        subroutine gensvb (ientmp, iBCBtmp, BCBtmp, mattmp, &
                           ienb,   iBCB,    BCB,    materb)
!
!----------------------------------------------------------------------
!
!  This routine saves the boundary element block.
!
! input:
!  ientmp (npro,nshl)           : boundary nodal connectivity
!  iBCBtmp (npro,ndiBCB)         : boundary condition codes
!  BCBtmp (npro,nshlb,ndBCB)    : boundary condition values
!  mattmp (npro)                : material type flag
!
! output:
!  ienb   (npro,nshl)           : boundary nodal connectivity
!  iBCB   (npro,ndiBCB)         : boundary condition codes
!  BCB    (npro,nshlb,ndBCB)    : boundary condition values
!  materb (npro)                : material type flag
!
!
! Zdenek Johan, Winter 1992.
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension   ientmp(npro,nshl), &
                    iBCBtmp(npro,ndiBCB),    BCBtmp(npro,ndBCB)

        dimension   mattmp(npro),           ienb(npro,nshl), &
                    iBCB(npro,ndiBCB),      BCB(npro,nshlb,ndBCB), &
                    materb(npro)
!
!.... generate the boundary element mapping
!
        do i = 1, nshl
          ienb(:,i) = ientmp(:,i)
        enddo
!
!.... save the boundary element data
!
        iBCB   = iBCBtmp
        do i = 1, nenbl ! This is NOT NSHLB as we are just copying the
                        ! piecewise constant data given by NSpre and
                        ! higher order coefficients must be zero
           do j = 1, ndBCB
              BCB(:,i,j)   = BCBtmp(:,j)
           end do
        end do
        do i = nenbl+1, nshlb
           do j = 1, ndBCB
              BCB(:,i,j)   = zero
           end do
        end do

        materb = mattmp
!
!.... return
!
        return
        end
