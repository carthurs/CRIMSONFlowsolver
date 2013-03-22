        subroutine gensvbDef (ientmp, iBCBtmp, BCBtmp, SWBtmp,  TWBtmp, &
                              mattmp, ienb,    iBCB,   BCB,     &
                              SWB,    TWB,  materb)
!
!----------------------------------------------------------------------
!
!  This routine saves the boundary element block.
!
! input:
!  ientmp (npro,nshl)           : boundary nodal connectivity
!  iBCtmp (npro,ndiBCB)         : boundary condition codes
!  BCBtmp (npro,nshlb,ndBCB)    : boundary condition values
!  SWBtmp (npro,nProps)         : Vessel Wall Properties
!  TWBtmp (npro,2)              : Tissue Support Properties
!  mattmp (npro)                : material type flag
!
! output:
!  ienb   (npro,nshl)           : boundary nodal connectivity
!  iBCB   (npro,ndiBCB)         : boundary condition codes
!  BCB    (npro,nshlb,ndBCB)    : boundary condition values
!  SWB    (npro,nProps)         : Vessel Wall Properties
!  TWB    (npro,2)              : Tissue Support Properties
!  materb (npro)                : material type flag
!
!
! Zdenek Johan, Winter 1992.
! Alberto Figueroa, Winter 2007.
!----------------------------------------------------------------------
!
         use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension   ientmp(npro,nshl), &
                    iBCBtmp(npro,ndiBCB),    BCBtmp(npro,ndBCB), &
                    SWBtmp(npro,nProps),     SWB(npro,nProps), &
                    TWBtmp(npro,2),          TWB(npro,2)

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
!.... save the vessel wall properties array
!
        do i = 1, nProps
          SWB(:,i) = SWBtmp(:,i)
        enddo
        
        ! note that due to a legacy compatibility issue,
        ! the expected entries for prestress in SWB are not consistent with
        ! the order of the stress components for the enchanced membrane
        ! as given in e3bvar !!
        ! in particular, the 3 and 4 indices need to be swapped

        SWB(:,3) = SWBtmp(:,4)
        SWB(:,4) = SWBtmp(:,3)

!
!.... save the tissue support properties array
!        
        do i = 1, 2
          TWB(:,i) = TWBtmp(:,i)
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
