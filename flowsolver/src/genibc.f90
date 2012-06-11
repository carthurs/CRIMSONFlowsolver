        subroutine geniBC (iBC)
!
!----------------------------------------------------------------------
! This routine reads the boundary condition codes.
!
! output: 
!  iBC   (nshg)        : Boundary Condition code
!
!         = 1 * iBC_1 + 2 * iBC_2 + 4 * iBC_3
!              density   temperature   pressure
!
!    if nsd = 3:
!
!        +  8 * iBC_4 +  16 * iBC_5 +  32 * iBC_6
!           x1-velocity   x2-velocity   x3-velocity
!
!        + 64 * iBC_7 + 128 * iBC_8 + 256 * iBC_9 + 512 * iBC_10
!          sclr1         sclr2        sclr3         sclr4
!
!        + 1024 * iBC_11  + 2048* iBC_12 
!          perioidicity     spebc          
!
!  nBC   (nshg)        : Boundary Condition mapping array
!
!
! Farzin Shakib, Winter 1986.
! Zdenek Johan,  Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
!
        use readarrays          ! used to access iBCtmp
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
! Arrays in the following 1 line are now dimensioned in readnblk
!        dimension iBCtmp(numpbc)
!
        dimension iBC(nshg)
        dimension itemp(6)
!
!.... set the iBC array
!
        iBC = 0
!
        if(numpbc.eq.0) return  ! sometimes there are no BC's on a partition
        where (nBC(:) .ne. 0) iBC(:) = iBCtmp(nBC(:))
!
!.... echo the input iBC array only if other than zero
!
        if (necho .lt. 3) then
          nn = 0
          do n = 1, nshg
            if (nBC(n) .ne. 0) then
              nb = nBC(n)
              nn = nn + 1
              if (mod(nn,50).eq.1) write(iecho,1000)ititle,(j,j=1,ndof)
              itemp(   1) = mod(iBCtmp(nb)   ,2) - mod(iBCtmp(nb)/ 4,2)
              itemp(   2) = mod(iBCtmp(nb)/ 8,2)
              itemp(   3) = mod(iBCtmp(nb)/16,2)
              itemp(   4) = mod(iBCtmp(nb)/32,2)
              itemp(ndof) = mod(iBCtmp(nb)/ 2,2)
              write(iecho,1100) n,(itemp(i),i=1,ndof)
            endif
          enddo
        endif
        deallocate(iBCtmp)
!
!.... return
!
        return
!
!.... end of file error handling
!
999     call error ('geniBC  ','end file',ibndc)
!
1000    format(a80,//, &
        ' N o d a l   B o u n d a r y   C o n d i t i o n   C o d e',//, &
        '    Node   ',13x,6('dof',i1,:,6x))
1100    format(2x,i5,10x,5i10)
!
        end
