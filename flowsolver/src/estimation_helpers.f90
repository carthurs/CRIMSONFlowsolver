subroutine setstate_comm() bind(C, name="setstate_comm")

    use iso_c_binding
    use shapeTable
    use globalArrays
    use phcommonvars

    implicit none

    if (numpe > 1) then
        call commu ( yold, ilwork, nflow, 'out')
        call commu ( acold, ilwork, nflow, 'out')
        call commu ( uold, ilwork, nsd, 'out')
    endif

end subroutine

subroutine temp_comm() bind(C, name="temp_comm")

    use iso_c_binding
    use shapeTable
    use globalArrays
    use phcommonvars

    implicit none

    if (numpe > 1) then
        call commu ( temporary_array, ilwork, nflow, 'out')
    endif

end subroutine

!
! this code is based off of the timeseries code
!
subroutine locate_soln(y, x, xts, soln)

    use phcommonvars
    use pointer_data  ! brings in the pointers for the blocked arrays
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    dimension y(nshg,ndofl), x(numnp,nsd)

    dimension ycl(npro,nshl,ndofl), xl(npro,nenl,nsd)

    real*8 xts(3), soln(5)

    do iblk = 1, nelblk
        iel    = lcblk(1,iblk)
        lcsyst = lcblk(3,iblk)
        iorder = lcblk(4,iblk)
        nenl   = lcblk(5,iblk) ! no. of vertices per element
        nshl   = lcblk(10,iblk)
        ndofl  = lcblk(8,iblk)
        npro   = lcblk(1,iblk+1) - iel

        call localx(x, xl, mien(iblk)%p, nsd, 'gather  ')
        call localy(y, ycl, mien(iblk)%p, ndofl, 'gather  ')

        call locate_soln_blk(ycl, xl, sgn, xts, soln)

    end do

end subroutine

subroutine locate_soln_blk(ycl, xl,  xts, soln)

    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    dimension shapeVar(npro,nshl), ycl(npro,nshl,ndofl), &
              xl(npro,nenl,nsd)

    real*8 al(npro,nenl,nsd), &
           zi0(npro,nsd), detaij(npro), dzi0(npro,nsd), &
           m11(npro), m12(npro), m13(npro), m21(npro), m22(npro), &
           m23(npro), m31(npro), m32(npro), m33(npro), &
           r1(npro), r2(npro), r3(npro), shgradl(npro,nshl,nsd)

    real*8 xts(3), soln(5)

    real*8 tolval

    integer e

    tolval = 2*2**-53

    if (lcsyst.eq.1) then !tet

        call get_a_not_tet(xl,al)

        detaij(:) = -al(:,2,1)*al(:,3,2)*al(:,4,3) + al(:,2,1)*al(:,4,2)*al(:,3,3) + al(:,2,2)*al(:,3,1)*al(:,4,3) - &
                     al(:,2,2)*al(:,4,1)*al(:,3,3) - al(:,2,3)*al(:,3,1)*al(:,4,2) + al(:,2,3)*al(:,4,1)*al(:,3,2)

        detaij = 1./detaij

        !
        ! solve for r, s, t  for all the elements
        !
        zi0(:,1) = detaij(:)*( ( al(:,4,2)*al(:,3,3) - al(:,3,2)*al(:,4,3) )*( xts(1)-al(:,1,1) ) + &
                               ( al(:,3,1)*al(:,4,3) - al(:,4,1)*al(:,4,3) )*( xts(2)-al(:,1,2) ) + &
                               ( al(:,4,1)*al(:,3,2) - al(:,3,1)*al(:,4,2) )*( xts(3)-al(:,1,3) ) )

        zi0(:,2) = detaij(:)*( ( al(:,2,2)*al(:,4,3) - al(:,4,2)*al(:,2,3) )*( xts(1)-al(:,1,1) ) + &
                               ( al(:,4,1)*al(:,2,3) - al(:,2,1)*al(:,4,3) )*( xts(2)-al(:,1,2) ) + &
                               ( al(:,2,1)*al(:,4,2) - al(:,4,1)*al(:,2,2) )*( xts(3)-al(:,1,3) ) )

        zi0(:,3) = detaij(:)*( ( al(:,3,2)*al(:,2,3) - al(:,2,2)*al(:,3,3) )*( xts(1)-al(:,1,1) ) + &
                               ( al(:,2,1)*al(:,3,3) - al(:,3,1)*al(:,2,3) )*( xts(2)-al(:,1,2) ) + &
                               ( al(:,3,1)*al(:,2,2) - al(:,2,1)*al(:,3,2) )*( xts(3)-al(:,1,3) ) )

        do e = 1, npro

            if (zi0(e,1).lt.(one+tolval).and.&
                zi0(e,1).gt.(zero-tolval).and.&
                zi0(e,2).lt.(one+tolval).and.&
                zi0(e,2).gt.(zero-tolval).and.&
                zi0(e,3).lt.(one+tolval).and.&
                zi0(e,3).gt.(zero-tolval)) then

                call shptet (ipord,zi0(e,:),shapeVar(e,:),shgradl(e,:,:))

                soln = zero

                do i = 1, nshl
                    soln(1) = soln(1)+ycl(e,i,1)*shapeVar(e,i) !pres
                    soln(2) = soln(2)+ycl(e,i,2)*shapeVar(e,i) !u
                    soln(3) = soln(3)+ycl(e,i,3)*shapeVar(e,i) !v
                    soln(4) = soln(4)+ycl(e,i,4)*shapeVar(e,i) !w
                enddo

            endif

        enddo

    endif

    return
end subroutine
