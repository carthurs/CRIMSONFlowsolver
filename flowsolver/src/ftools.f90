!---------------------------------------------------------------------
!
! ftools.f : Bundle of Fortran routines
!
! Each routine is to be called by drvftools.f
!
! Various operations run based on les**.c
!
!---------------------------------------------------------------------
!
!--------------
! flesPrepDiag
!--------------
!
 	subroutine flesPrepDiag ( ien, xKebe, xGoc, flowDiag )
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension xKebe(npro,3*nshl,3*nshl),  &
                  xGoC(npro,4*nshl,nshl)
        dimension Diagl(npro,nshl,nflow),   flowDiag(nshg, 4)
        integer   ien(npro,nshl)
!
!
!.... monentum contribution to diagonal
!
        do i = 1, nshl
          i0 = (nsd) * (i - 1)
          Diagl(:,i,1) = xKebe(1:npro,i0+1,i0+1)
          Diagl(:,i,2) = xKebe(1:npro,i0+2,i0+2)
          Diagl(:,i,3) = xKebe(1:npro,i0+3,i0+3)
        enddo
!
!.... continuity contribution to diagonal
!
        nn2 = nshl * nsd
        do i = 1, nshl
          Diagl(:,i,4) = xGoC(1:npro,nn2+i,i)
        enddo

        call local (flowDiag,  Diagl, ien, nflow, 'scatter ')
!
        return
        end
!
!--------------------------------
! fMtxVdimVecMult 
! row and column index exchanged
!--------------------------------
!
        subroutine fMtxVdimVecMult( a, b, c, na, nb, nc, m, n )
!
!.... Data declaration
!
        implicit none
        integer na,     nb,     nc,     m,      n
        real*8  a(n,na),        b(n,nb),        c(n,nc)
!
        integer i,      j
!
!.... Do the work
!
! WIP: change to F90
!
        if ( m .eq. 1 ) then
!
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
            enddo
!
        else if ( m .eq. 2 ) then
!
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
                c(i,2) = a(i,2) * b(i,2)
            enddo
!
        else if ( m .eq. 3 ) then
!
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
                c(i,2) = a(i,2) * b(i,2)
                c(i,3) = a(i,3) * b(i,3)
            enddo
!
        else if ( m .eq. 4 ) then
!
            do i = 1, n
                c(i,1) = a(i,1) * b(i,1)
                c(i,2) = a(i,2) * b(i,2)
                c(i,3) = a(i,3) * b(i,3)
                c(i,4) = a(i,4) * b(i,4)
            enddo
!
        else
!
            do i = 1, n 
                do j = 1, m 
                    c(i,j) = a(i,j) * b(i,j)
                enddo
            enddo
!
        endif
!
      return
      end
!
!---------- 
! flesZero
!----------
!
	subroutine flesZero ( a, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m)
!
        do i = 1, n
          do j = 1, m
            a(i,j) = 0.0e-0
          enddo
        enddo
!
        return 
        end
!
!--------
! flesCp
!--------
!
	subroutine flesCp ( a, b, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
!
        do i = 1, n
          do j = 1, m
            b(i,j) = a(i,j)
          enddo
        enddo
! 
        return
        end
!
!-----------
! flesScale
!-----------
!
 	subroutine flesScale ( a, s, m, n )
!
        implicit none
        integer  m, n, i, j   
        real*8   a(n,m), s
!
        do i = 1, n
          do j = 1, m
            a(i,j) = a(i,j) * s
          enddo
        enddo
!
        return
        end
!
!-------------
! flesScaleCp
!-------------
!
	subroutine flesScaleCp ( a, b, s, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m), s
!
        do i = 1, n
          do j = 1, m
            b(i,j) = a(i,j) * s
          enddo
        enddo
!
        return 
        end
!
!---------
! flesAdd
!---------
!
	subroutine flesAdd ( a, b, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
!
        do i = 1, n
          do j = 1, m
            b(i,j) = b(i,j) + a(i,j)
          enddo
        enddo
!
        return
        end
!
!---------
! flesSub 
!---------
!
	subroutine flesSub ( a, b, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
!
        do i = 1, n
          do j = 1, m
            b(i,j) = b(i,j) - a(i,j)
          enddo
        enddo
!
        return
        end
!
!----------
! flesDot1
!----------
!
 	real*8 function flesDot1 ( a, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m)
!
	flesDot1 = 0
        do i = 1, n
          do j = 1, m
            flesDot1 = flesDot1 + a(i,j) * a(i,j)
          enddo
        enddo
!
        return 
        end
!
!----------
! flesDot2
!----------
!
	real*8 function flesDot2 ( a, b, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   a(n,m), b(n,m)
!
	flesDot2 = 0
        do i = 1, n
          do j = 1, m
            flesDot2 = flesDot2 + a(i,j) * b(i,j)
          enddo
        enddo
!
        return
        end
!
!-----------
! flesDaxpy
!-----------
!
	subroutine flesDaxpy ( x, y, a, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   x(n,m), y(n,m)
        real*8   a
!
        do i = 1, n
          do j= 1, m
            y(i,j) = y(i,j) + a * x(i,j)
          enddo
        enddo
!
        return
        end
!
!-----------
! flesDxpay
!-----------
!
	subroutine flesDxpay ( x, y, a, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   x(n,m), y(n,m)
        real*8   a
!
        do i = 1, n
          do j = 1, m
            y(i,j) = a * y(i,j) + x(i,j)
          enddo
        enddo
!
        return 
        end
!
!---------
! flesInv
!---------
!
	subroutine flesInv ( x, m, n )
!
        implicit none
        integer  m, n, i, j
        real*8   x(n,m)
!
        do i = 1, n
          do j = 1, m
            if ( x(i,j) .ne. 0 ) x(i,j) = 1. / x(i,j)
          enddo
        enddo
!
        return
        end
!
!--------------------------
! fMtxBlkDot2
! row and column exchanged
!--------------------------
!
        subroutine fMtxBlkDot2( x, y, c, m, n )
!
!.... Data declaration
!
        implicit none
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
!
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
!
!.... Determine the left overs
!
        m1 = mod(m,8) + 1

!
!.... Do the small pieces
!
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
!
1000    continue
        tmp1 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
        enddo
        c(1) = tmp1
        goto 8000
!
2000    continue
        tmp1 = 0
        tmp2 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        goto 8000
!
3000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        goto 8000
!
4000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        goto 8000
!
5000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        tmp5 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
            tmp5 = tmp5 + x(i,5) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        c(5) = tmp5
        goto 8000
!
6000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        tmp5 = 0
        tmp6 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
            tmp5 = tmp5 + x(i,5) * y(i)
            tmp6 = tmp6 + x(i,6) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        c(5) = tmp5
        c(6) = tmp6
        goto 8000
!
7000    continue
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        tmp5 = 0
        tmp6 = 0
        tmp7 = 0
        do i = 1, n
            tmp1 = tmp1 + x(i,1) * y(i)
            tmp2 = tmp2 + x(i,2) * y(i)
            tmp3 = tmp3 + x(i,3) * y(i)
            tmp4 = tmp4 + x(i,4) * y(i)
            tmp5 = tmp5 + x(i,5) * y(i)
            tmp6 = tmp6 + x(i,6) * y(i)
            tmp7 = tmp7 + x(i,7) * y(i)
        enddo
        c(1) = tmp1
        c(2) = tmp2
        c(3) = tmp3
        c(4) = tmp4
        c(5) = tmp5
        c(6) = tmp6
        c(7) = tmp7
        goto 8000
!
!.... Do the remaining part
!
8000    continue
!
        do j = m1, m, 8
            tmp1 = 0
            tmp2 = 0
            tmp3 = 0
            tmp4 = 0
            tmp5 = 0
            tmp6 = 0
            tmp7 = 0
            tmp8 = 0
            do i = 1, n
                tmp1 = tmp1 + x(i,j+0) * y(i)
                tmp2 = tmp2 + x(i,j+1) * y(i)
                tmp3 = tmp3 + x(i,j+2) * y(i)
                tmp4 = tmp4 + x(i,j+3) * y(i)
                tmp5 = tmp5 + x(i,j+4) * y(i)
                tmp6 = tmp6 + x(i,j+5) * y(i)
                tmp7 = tmp7 + x(i,j+6) * y(i)
                tmp8 = tmp8 + x(i,j+7) * y(i)
            enddo
            c(j+0) = tmp1
            c(j+1) = tmp2
            c(j+2) = tmp3
            c(j+3) = tmp4
            c(j+4) = tmp5
            c(j+5) = tmp6
            c(j+6) = tmp7
            c(j+7) = tmp8
        enddo
!
        return
        end
!
!--------------------------
! fMtxBlkDaxpy
! row and column exchanged
!--------------------------
!
        subroutine fMtxBlkDaxpy( x, y, c, m, n )
!
!.... Data declaration
!
        implicit none
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
!
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
!
!.... Determine the left overs
!
        m1 = mod(m,8) + 1
!
!.... Do the small pieces
!
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
!
1000    continue
        tmp1 = c(1)
        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1)
        enddo
        goto 8000
!
2000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1) + tmp2 * x(i,2)
        enddo
        goto 8000
!
3000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)

        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3)
        enddo
        goto 8000
!
4000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4)
        enddo
        goto 8000
!
5000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4) &
                 + tmp5 * x(i,5)
        enddo
        goto 8000
!
6000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4) &
                 + tmp5 * x(i,5) + tmp6 * x(i,6)
        enddo
        goto 8000
!
7000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        tmp7 = c(7)
        do i = 1, n
            y(i) = y(i) &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4) &
                 + tmp5 * x(i,5) + tmp6 * x(i,6) &
                 + tmp7 * x(i,7)
        enddo
        goto 8000
!
!.... Do the remaining part
!
8000    continue
!
        do j = m1, m, 8
            tmp1 = c(j+0)
            tmp2 = c(j+1)
            tmp3 = c(j+2)
            tmp4 = c(j+3)
            tmp5 = c(j+4)
            tmp6 = c(j+5)
            tmp7 = c(j+6)
            tmp8 = c(j+7)
            do i = 1, n
                y(i) = y(i) &
                     + tmp1 * x(i,j+0) + tmp2 * x(i,j+1) &
                     + tmp3 * x(i,j+2) + tmp4 * x(i,j+3) &
                     + tmp5 * x(i,j+4) + tmp6 * x(i,j+5) &
                     + tmp7 * x(i,j+6) + tmp8 * x(i,j+7)
            enddo
        enddo
!
        return 
        end
!
!--------------------------
! fMtxBlkDyeax
! row and column exchanged
!--------------------------
!
        subroutine fMtxBlkDyeax( x, y, c, m, n )
!
!.... Data declaration
!
        implicit none
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
!
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
!
!.... Determine the left overs
!
        m1 = mod(m,8) + 1
!
!.... Do the small pieces
!
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
!
1000    continue
        tmp1 = c(1)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1)
        enddo
        goto 8001
!
2000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1) + tmp2 * x(i,2)
        enddo
        goto 8001
!
3000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3)
        enddo
        goto 8001
!
4000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4)
        enddo
        goto 8001
!
5000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4) &
                 + tmp5 * x(i,5)
        enddo
        goto 8001
!
6000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4) &
                 + tmp5 * x(i,5) + tmp6 * x(i,6)
       enddo
        goto 8001
!
7000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        tmp7 = c(7)
        do i = 1, n
            y(i) = &
                 + tmp1 * x(i,1) + tmp2 * x(i,2) &
                 + tmp3 * x(i,3) + tmp4 * x(i,4) &
                 + tmp5 * x(i,5) + tmp6 * x(i,6) &
                 + tmp7 * x(i,7)
        enddo
        goto 8001
!
8000    continue
        do i = 1, n
            y(i) = 0
        enddo
        goto 8001
!
!.... Do the remaining part
!
8001    continue
!
        do j = m1, m, 8
            tmp1 = c(j+0)
            tmp2 = c(j+1)
            tmp3 = c(j+2)
            tmp4 = c(j+3)
            tmp5 = c(j+4)
            tmp6 = c(j+5)
            tmp7 = c(j+6)
            tmp8 = c(j+7)
            do i = 1, n
                y(i) = y(i) &
                     + tmp1 * x(i,j+0) + tmp2 * x(i,j+1) &
                     + tmp3 * x(i,j+2) + tmp4 * x(i,j+3) &
                     + tmp5 * x(i,j+4) + tmp6 * x(i,j+5) &
                     + tmp7 * x(i,j+6) + tmp8 * x(i,j+7)
            enddo
        enddo
!
        return 
        end
!
!--------------------------
! fMtxBlkDmaxpy
! row and column exchanged 
!--------------------------
!
       subroutine fMtxBlkDmaxpy( x, y, c, m, n )
!
!.... Data declaration
!
        implicit none
        integer m,      n
        real*8  x(n,m), y(n),   c(m)
!
        real*8  tmp1,   tmp2,   tmp3,   tmp4
        real*8  tmp5,   tmp6,   tmp7,   tmp8
        integer i,      j,      m1
!
!.... Determine the left overs
!
        m1 = mod(m,8) + 1
!
!.... Do the small pieces
!
        goto ( 8000, 1000, 2000, 3000, 4000, 5000, 6000, 7000 ) m1
!
1000    continue
        tmp1 = c(1)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1)
        enddo
        goto 8000
!
2000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1) - tmp2 * x(i,2)
        enddo
        goto 8000
!
3000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1) - tmp2 * x(i,2) &
                 - tmp3 * x(i,3)
        enddo
        goto 8000
!
4000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1) - tmp2 * x(i,2) &
                 - tmp3 * x(i,3) - tmp4 * x(i,4) 
        enddo
        goto 8000
!
5000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1) - tmp2 * x(i,2) &
                 - tmp3 * x(i,3) - tmp4 * x(i,4) &
                 - tmp5 * x(i,5)
        enddo
        goto 8000
!
6000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1) - tmp2 * x(i,2) &
                 - tmp3 * x(i,3) - tmp4 * x(i,4) &
                 - tmp5 * x(i,5) - tmp6 * x(i,6)
        enddo
        goto 8000

7000    continue
        tmp1 = c(1)
        tmp2 = c(2)
        tmp3 = c(3)
        tmp4 = c(4)
        tmp5 = c(5)
        tmp6 = c(6)
        tmp7 = c(7)
        do i = 1, n
            y(i) = y(i) &
                 - tmp1 * x(i,1) - tmp2 * x(i,2) &
                 - tmp3 * x(i,3) - tmp4 * x(i,4) &
                 - tmp5 * x(i,5) - tmp6 * x(i,6) &
                 - tmp7 * x(i,7)
        enddo
        goto 8000
!
!.... Do the remaining part
!
8000    continue
!
        do j = m1, m, 8
            tmp1 = c(j+0)
            tmp2 = c(j+1)
            tmp3 = c(j+2)
            tmp4 = c(j+3)
            tmp5 = c(j+4)
            tmp6 = c(j+5)
            tmp7 = c(j+6)
            tmp8 = c(j+7)
            do i = 1, n
                y(i) = y(i) &
                     - tmp1 * x(i,j+0) - tmp2 * x(i,j+1) &
                     - tmp3 * x(i,j+2) - tmp4 * x(i,j+3) &
                     - tmp5 * x(i,j+4) - tmp6 * x(i,j+5) &
                     - tmp7 * x(i,j+6) - tmp8 * x(i,j+7)
            enddo
        enddo
!
        return
        end
!
!--------------------------
! fMtxVdimVecCp
! row and column exchanged 
!--------------------------
!
        subroutine fMtxVdimVecCp( a, b, na, nb, m, n )
!
!.... Data declaration
!
        implicit none
        integer na,     nb,     m,      n
        real*8  a(n,na),        b(n,nb)
!
        integer i,      j
!
!.... Do the work
!
        if ( m .eq. 1 ) then

            do i = 1, n
                b(i,1) = a(i,1)
            enddo

        else if ( m .eq. 2 ) then

            do i = 1, n
                b(i,1) = a(i,1)
                b(i,2) = a(i,2)
            enddo

        else if ( m .eq. 3 ) then

            do i = 1, n
                b(i,1) = a(i,1)
                b(i,2) = a(i,2)
                b(i,3) = a(i,3)
            enddo

        else if ( m .eq. 4 ) then

            do i = 1, n
                b(i,1) = a(i,1)
                b(i,2) = a(i,2)
                b(i,3) = a(i,3)
                b(i,4) = a(i,4)
            enddo

        else

            do i = 1, n
                do j = 1, m 
                    b(i,j) = a(i,j)
                enddo
            enddo

        endif
!
        return
        end
!
!--------------------------
! fMtxVdimVecDot2
! row and column exchanged
!--------------------------
!
        subroutine fMtxVdimVecDot2( a, b, c, na, nb, m, n )
!
!.... Data declaration
!
        implicit none
        integer na,     nb,     m,      n
        real*8  a(n,na),        b(n,nb),        c(m)
!
        integer i,      j
!
!.... Do the work
!
        if ( m .eq. 1 ) then

            c(1) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
            enddo

        else if ( m .eq. 2 ) then

            c(1) = 0
            c(2) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
                c(2) = c(2) + a(i,2) * b(i,2)
            enddo

        else if ( m .eq. 3 ) then

            c(1) = 0
            c(2) = 0
            c(3) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
                c(2) = c(2) + a(i,2) * b(i,2)
                c(3) = c(3) + a(i,3) * b(i,3)
            enddo

        else if ( m .eq. 4 ) then

            c(1) = 0
            c(2) = 0
            c(3) = 0
            c(4) = 0
            do i = 1, n
                c(1) = c(1) + a(i,1) * b(i,1)
                c(2) = c(2) + a(i,2) * b(i,2)
                c(3) = c(3) + a(i,3) * b(i,3)
                c(4) = c(4) + a(i,4) * b(i,4)
            enddo

        else

            do j = 1, m 
                c(j) = 0
                do i = 1, n 
                    c(j) = c(j) + a(i,j) * b(i,j)
                enddo
            enddo

        endif
!
        return
        end
!
!--------------------------
! fMtxVdimVecDaxpy
! row and column exchanged
!--------------------------
!
        subroutine fMtxVdimVecDaxpy( a, b, c, na, nb, m, n )
!
!.... Data declaration
!
        implicit none
        integer na,     nb,     m,      n
        real*8  a(n,na),        b(n,nb),        c(m)
!
        integer i,      j
!
!.... Do the work
!
        if ( m .eq. 1 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
            enddo

        else if ( m .eq. 2 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
                b(i,2) = b(i,2) + c(2) * a(i,2)
            enddo

        else if ( m .eq. 3 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
                b(i,2) = b(i,2) + c(2) * a(i,2)
                b(i,3) = b(i,3) + c(3) * a(i,3)
            enddo

        else if ( m .eq. 4 ) then

            do i = 1, n
                b(i,1) = b(i,1) + c(1) * a(i,1)
                b(i,2) = b(i,2) + c(2) * a(i,2)
                b(i,3) = b(i,3) + c(3) * a(i,3)
                b(i,4) = b(i,4) + c(4) * a(i,4)
            enddo

        else

            do j = 1, m 
                do i = 1, n 
                    b(i,j) = b(i,j) + c(j) * a(i,j)
                enddo
            enddo

        endif
!
        return
        end
!
!---------
! flesApG
!---------
!
	subroutine flesApG ( ien, xGoC, lesP, lesQ, nPs, nQs )
!   
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension xGoC(npro,4*nshl,nshl)
        real*8 lesP(nshg,nPs), lesQ(nshg,nQs) 
        dimension ien(npro,nshl)
        dimension Ptemp(npro,nshl,nPs), Qtemp(npro,nshl,nQs)
!
!.... zero Qtemp
!
	Qtemp = zero
!
!.... localize the lesP for the EBE product
!
        call local ( lesP, Ptemp, ien, nPs, 'gather  ' )
!
!.... Now, product operation
!
    	do i = 1, nshl
           i0 = (nsd) * (i - 1)  
           do j = 1, nshl
!
             Qtemp(:,i,1) = Qtemp(:,i,1)  &
                          + xGoC(1:npro,i0+1,j) * Ptemp(:,j,nPs)
!
             Qtemp(:,i,2) = Qtemp(:,i,2) &
                          + xGoC(1:npro,i0+2,j) * Ptemp(:,j,nPs)
!
             Qtemp(:,i,3) = Qtemp(:,i,3) &
                          + xGoC(1:npro,i0+3,j) * Ptemp(:,j,nPs) 
!
           enddo
        enddo
!
!... assemble the result of the product
!
        call local ( lesQ, Qtemp, ien, nQs, 'scatter ' )
!
        return 
        end
!
!----------
! flesApKG
!----------
!
 	subroutine flesApKG ( ien, xKebe, xGoC, lesP, lesQ, nPs, nQs )
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension xKebe(npro,3*nshl,3*nshl),  &
             xGoC(npro,4*nshl,nshl)
        dimension ien(npro,nshl)
        real*8 lesP(nshg,nPs), lesQ(nshg,nQs)
        dimension Ptemp(npro,nshl,nPs), Qtemp(npro,nshl,nQs)       
!
!.... zero Qtemp
!
	Qtemp = zero
!
!.... localize the lesP for the EBE product
!
        call local ( lesP, Ptemp, ien, nPs, 'gather  ' )
!
!.... Now, product operation
!
!.... K contribution 
!
        do i = 1, nshl
           i0 = (nsd) * (i - 1)
          do j = 1, nshl
             j0 = (nsd) * (j - 1) 
!
            Qtemp(:,i,1) = Qtemp(:,i,1) &
                         + xKebe(1:npro,i0+1,j0+1) * Ptemp(:,j,1) &
                         + xKebe(1:npro,i0+1,j0+2) * Ptemp(:,j,2)     &
                         + xKebe(1:npro,i0+1,j0+3) * Ptemp(:,j,3)
!
            Qtemp(:,i,2) = Qtemp(:,i,2) &
                         + xKebe(1:npro,i0+2,j0+1) * Ptemp(:,j,1) &
                         + xKebe(1:npro,i0+2,j0+2) * Ptemp(:,j,2) &
                         + xKebe(1:npro,i0+2,j0+3) * Ptemp(:,j,3)
            Qtemp(:,i,3) = Qtemp(:,i,3) &
                         + xKebe(1:npro,i0+3,j0+1) * Ptemp(:,j,1) &
                         + xKebe(1:npro,i0+3,j0+2) * Ptemp(:,j,2) &
                         + xKebe(1:npro,i0+3,j0+3) * Ptemp(:,j,3)
!
          enddo
     	enddo
!
!.... G contribution
!
        do i = 1, nshl 
           i0 = (nsd) * (i - 1) 
          do j = 1, nshl
!
            Qtemp(:,i,1) = Qtemp(:,i,1) &
                         + xGoC(1:npro,i0+1,j) * Ptemp(:,j,nPs)
            Qtemp(:,i,2) = Qtemp(:,i,2) &
                         + xGoC(1:npro,i0+2,j) * Ptemp(:,j,nPs)
            Qtemp(:,i,3) = Qtemp(1:,i,3) &
                         + xGoC(1:npro,i0+3,j) * Ptemp(:,j,nPs)
!
          enddo
        enddo
!
!.... assemble the result of the product
!
        call local ( lesQ, Qtemp, ien, nQs, 'scatter ' )
!
        return
        end
!
!-----------
! flesApNGt
!-----------
!
	subroutine flesApNGt ( ien, xGoC, lesP, lesQ, nPs, nQs )
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension ien(npro,nshl), xGoC(npro,4*nshl,nshl)
        real*8 lesP(nshg,nPs), lesQ(nshg,nQs)
        dimension Ptemp(npro,nshl,nPs), Qtemp(npro,nshl,nQs)
!
!.... zero Qtemp
!
	Qtemp = zero 
!
!.... localize the lesP for the EBE product
!
        call local ( lesP, Ptemp, ien, nPs, 'gather  ' )
!
!.... Now, product operation
!
!.... Negative G^t contribution ( not explicitly formed )
!
        do i = 1, nshl
           do j = 1, nshl
              i0 = (nsd) * (j - 1)
!
             Qtemp(:,i,nQs) = Qtemp(:,i,nQs) &
                            - xGoC(1:npro,i0+1,i) * Ptemp(:,j,1) &
                            - xGoC(1:npro,i0+2,i) * Ptemp(:,j,2) &
                            - xGoC(1:npro,i0+3,i) * Ptemp(:,j,3)
!
           enddo
        enddo
!
!... assemble the result of the product
!
        call local ( lesQ, Qtemp, ien, nQs, 'scatter  ' )
!
        return
        end
!
!------------
! flesApNGtC 
!------------
!
	subroutine flesApNGtC ( ien, xGoC, lesP, lesQ, nPs, nQs )
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension ien(npro,nshl), xGoC(npro,4*nshl,nshl)
        real*8 lesP(nshg,nPs), lesQ(nshg,nQs)
        dimension Ptemp(npro,nshl,nPs), Qtemp(npro,nshl,nQs)
!
!.... zero Qtemp
!
	Qtemp = zero
!
!.... localize the lesP for the EBE product
!
	call local ( lesP, Ptemp, ien, nPs, 'gather  ')
!
!.... Now, product operation
!
!.... Negative G^t contribution ( not explicitly formed )
!
        do i = 1, nshl
           do j = 1, nshl
           i0 = (nsd) * (j - 1)
!
             Qtemp(:,i,nQs) = Qtemp(:,i,nQs) &
                            - xGoC(1:npro,i0+1,i) * Ptemp(:,j,1) &
                            - xGoC(1:npro,i0+2,i) * Ptemp(:,j,2) &
                            - xGoC(1:npro,i0+3,i) * Ptemp(:,j,3)
!
           enddo
        enddo
!
!.... C contribution
!
        nnm2 = nshl * (nsd)
!
        do i = 1, nshl
           i0 = nnm2 + i
          do j = 1, nshl
!
             Qtemp(:,i,nQs) = Qtemp(:,i,nQs) &
                            + xGoC(1:npro,i0,j) * Ptemp(:,j,nPs)
!
          enddo
        enddo
!
!... assemble the result of the product
!
        call local ( lesQ, Qtemp, ien, nQs, 'scatter  ' )
!
        return
        end
!
!------------
! flesApFull
!------------
!
	subroutine flesApFull ( ien, xKebe, xGoC, lesP, lesQ, nPs, nQs )
!   
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension ien(npro,nshl)
        dimension xKebe(npro,3*nshl,3*nshl),  &
             xGoC(npro,4*nshl,nshl)
        real*8 lesP(nshg,nPs), lesQ(nshg,nQs)
        dimension Ptemp(npro,nshl,nPs), Qtemp(npro,nshl,nQs)
!
!.... zero Qtemp
!
	Qtemp = zero
!
!.... localize the lesP for the EBE product
!
 	call local ( lesP, Ptemp, ien, nPs, 'gather  ' )
!
!.... Now, product operation
!
!.... K * Du contribution
!
        do i = 1, nshl
           i0 = (nsd) * (i - 1)
          do j = 1, nshl
             j0 = (nsd) * (j - 1)
!
            Qtemp(:,i,1) = Qtemp(:,i,1) &
                         + xKebe(1:npro,i0+1,j0+1) * Ptemp(:,j,1) &
                         + xKebe(1:npro,i0+1,j0+2) * Ptemp(:,j,2) &
                         + xKebe(1:npro,i0+1,j0+3) * Ptemp(:,j,3)
!
            Qtemp(:,i,2) = Qtemp(:,i,2) &
                         + xKebe(1:npro,i0+2,j0+1) * Ptemp(:,j,1) &
                         + xKebe(1:npro,i0+2,j0+2) * Ptemp(:,j,2) &
                         + xKebe(1:npro,i0+2,j0+3) * Ptemp(:,j,3)
            Qtemp(:,i,3) = Qtemp(:,i,3) &
                         + xKebe(1:npro,i0+3,j0+1) * Ptemp(:,j,1) &
                         + xKebe(1:npro,i0+3,j0+2) * Ptemp(:,j,2) &
                         + xKebe(1:npro,i0+3,j0+3) * Ptemp(:,j,3)
!
          enddo
        enddo
!
!.... G * Dp contribution
!
       do i = 1, nshl
           i0 = (nsd) * (i - 1)
          do j = 1, nshl
!
            Qtemp(:,i,1) = Qtemp(:,i,1) &
                         + xGoC(1:npro,i0+1,j) * Ptemp(:,j,nPs)
            Qtemp(:,i,2) = Qtemp(:,i,2) &
                         + xGoC(1:npro,i0+2,j) * Ptemp(:,j,nPs)
            Qtemp(:,i,3) = Qtemp(:,i,3) &
                         + xGoC(1:npro,i0+3,j) * Ptemp(:,j,nPs)
!
          enddo
        enddo
!
!.... -G^t * Du contribution
!
       do i = 1, nshl
           do j = 1, nshl
              i0 = (nsd) * (j - 1)
!
             Qtemp(:,i,nQs) = Qtemp(:,i,nQs) &
                            - xGoC(1:npro,i0+1,i) * Ptemp(:,j,1) &
                            - xGoC(1:npro,i0+2,i) * Ptemp(:,j,2) &
                            - xGoC(1:npro,i0+3,i) * Ptemp(:,j,3)
!
           enddo
        enddo
!
!.... C * Dp contribution
!
        nnm2 = nshl * (nsd)
!
        do i = 1, nshl
           i0 = nnm2 + i
          do j = 1, nshl
!
             Qtemp(:,i,nQs) = Qtemp(:,i,nQs) &
                            + xGoC(1:npro,i0,j) * Ptemp(:,j,nPs)
!
          enddo
        enddo



!
!... assemble the result of the product
!
        call local ( lesQ, Qtemp, ien, nQs, 'scatter ' )
!
        return
        end
!
!-----------
! fsclrDiag
!-----------
!
	subroutine fsclrDiag ( ien, xTe, sclrDiag )
! 
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension xTe(npro,nshl,nshl)
        dimension sclrDiag(nshg,1), Diagl(npro,nshl,1)
        dimension ien(npro,nshl)
!
        do i = 1, nshl
           Diagl(:,i,1) = xTe(1:npro,i,i)
        enddo
!
        call local (sclrDiag, Diagl, ien, 1, 'scatter ')
!  
        return
        end
!
!------------
! flesApSclr
!------------
!
	subroutine flesApSclr ( ien, xTe, lesP, lesQ, nPs, nQs )
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension xTe(npro,nshl,nshl)
        dimension ien(npro,nshl)
        real*8 lesP(nshg,nPs), lesQ(nshg,nQs)
        dimension Ptemp(npro,nshl,nPs), Qtemp(npro,nshl,nQs)
!
!.... zero Qtemp
!
        Qtemp = zero
!
!.... localize the lesP for the EBE product
!
        call local ( lesP, Ptemp, ien, nPs, 'gather  ')
!
!.... Now, product operation
!
        do i = 1, nshl
          do j = 1, nshl
!
            Qtemp(:,i,nQs) = Qtemp(:,i,nQs)  &
                           + xTe(1:npro,i,j) * Ptemp(:,j,nPs)
!
          enddo
        enddo
!
!.... assemble the result of the product
!
  	call local ( lesQ, Qtemp, ien, nQs, 'scatter ' )
!
        return
        end
