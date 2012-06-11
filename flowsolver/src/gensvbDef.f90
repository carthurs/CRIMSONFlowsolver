        subroutine gensvbDef (ientmp, iBCBtmp, BCBtmp, SWBtmp,  TWBtmp, &
                              EWBtmp, mattmp, ienb,    iBCB,   BCB,     &
                              SWB,    TWB,    EWB,     materb)
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
!  EWBtmp (npro,1)              : State Filter Properties
!  mattmp (npro)                : material type flag
!
! output:
!  ienb   (npro,nshl)           : boundary nodal connectivity
!  iBCB   (npro,ndiBCB)         : boundary condition codes
!  BCB    (npro,nshlb,ndBCB)    : boundary condition values
!  SWB    (npro,nProps)         : Vessel Wall Properties
!  TWB    (npro,2)              : Tissue Support Properties
!  EWB    (npro,1)              : State Filter Properties
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
                    TWBtmp(npro,2),          TWB(npro,2), &
                    EWBtmp(npro,1),          EWB(npro,1)

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
        
!
!.... save the tissue support properties array
!        
        do i = 1, 2
          TWB(:,i) = TWBtmp(:,i)
        enddo
        
!
!.... save the state filter properties array
!        
        do i = 1, 1
          EWB(:,i) = EWBtmp(:,i)
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
