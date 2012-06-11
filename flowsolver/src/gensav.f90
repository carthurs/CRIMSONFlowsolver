        subroutine gensav (ientmp, mattmp, ien,    mater)
!
!----------------------------------------------------------------------
!
!  This routine saves the element block data.
!
! input:
!  ientmp (npro,nshl)   : nodal connectivity 
!  mattmp (npro)        : material type flag
!
! output:
!  ien    (npro,nshl)   : nodal connectivity
!  mater  (npro)        : material type flag
!
!
! Zdenek Johan, Winter 1992.
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension   ientmp(npro,nshl), &
                    mattmp(npro),           ien(npro,nshl), &
                    mater(npro)
!
!.... save the element data
!
        do i = 1, nshl
          ien(:,i) = ientmp(:,i)
        enddo
!
        mater = mattmp
!
!.... end
!
        return
        end
