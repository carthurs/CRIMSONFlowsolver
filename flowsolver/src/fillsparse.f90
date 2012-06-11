      subroutine fillsparseI(iens, xKebe, lhsK, &
                             xGoC, lhsP, &
                             row,  col)
!
!
!
      use LagrangeMultipliers 
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      real*8  xKebe(npro,9,nshl,nshl), xGoC(npro,4,nshl,nshl)
      integer ien(npro,nshl),   col(nshg+1), row(nshg*nnz)
      real*8  lhsK(9,nnz_tot),  lhsP(4,nnz_tot)
!
      integer aa, b, c, e, i, k, n, j, l
!
      integer sparseloc

      integer iens(npro,nshl)
!
! prefer to show explicit absolute value needed for cubic modes and
! higher rather than inline abs on pointer as in past versions
! iens is the signed ien array ien is unsigned
!
      ien=abs(iens)
!       
!.... Accumulate the lhs
!
      do e = 1, npro ! loop over the elements
         do aa = 1, nshl ! loop over the local equation numbers
            i = ien(e,aa) ! finds the global equation number or
                          ! block-row of our matrix
            c = col(i)    ! starting point to look for the matching column
            n = col(i+1) - c  !length of the list of entries in rowp
            do b = 1, nshl ! local variable number tangent respect
                           ! to
! function that searches row until it finds the match that gives the
!		   global equation number

               k = sparseloc(row(c),n,ien(e,b))+c-1
!
!                                             *         *
!                   dimension egmass(npro,ndof,nenl,ndof,nenl)
!
! compressible      lhsT(1:5,1:5,k)=lhsT(1:5,1:5,k)+egmass(e,1:5,aa,1:5,b)
!
               lhsK(1,k) = lhsK(1,k) + xKebe(e,1,aa,b)
               lhsK(2,k) = lhsK(2,k) + xKebe(e,2,aa,b)
               lhsK(3,k) = lhsK(3,k) + xKebe(e,3,aa,b)
               lhsK(4,k) = lhsK(4,k) + xKebe(e,4,aa,b)
               lhsK(5,k) = lhsK(5,k) + xKebe(e,5,aa,b)
               lhsK(6,k) = lhsK(6,k) + xKebe(e,6,aa,b)
               lhsK(7,k) = lhsK(7,k) + xKebe(e,7,aa,b)
               lhsK(8,k) = lhsK(8,k) + xKebe(e,8,aa,b)
               lhsK(9,k) = lhsK(9,k) + xKebe(e,9,aa,b)
!
               lhsP(1,k) = lhsP(1,k) + xGoC(e,1,aa,b)
               lhsP(2,k) = lhsP(2,k) + xGoC(e,2,aa,b)
               lhsP(3,k) = lhsP(3,k) + xGoC(e,3,aa,b)
               lhsP(4,k) = lhsP(4,k) + xGoC(e,4,aa,b)
            enddo
         enddo
      enddo
!
!..... lshLagL is required to calculate the LHS and RHS of Lagrange Multipliers
!..... It should be assembled before calling CalcNANBLagrange
!
      if(Lagrange.gt.zero) then
         do e = 1, npro ! loop over the elements
            do aa = 1, nshlb ! loop over the local equation numbers
               i = ien(e,aa) 
               c = col(i)    
               n = col(i+1) - c  
               do b = 1, nshlb
                  k = sparseloc( row(c), n, ien(e,b) ) + c-1
!
                  lhsLagL(1,k,:)=lhsLagL(1,k,:)+ &
                     loclhsLag(e,1,aa,b,:)
                  lhsLagL(2,k,:)=lhsLagL(2,k,:)+ &
                     loclhsLag(e,2,aa,b,:)
                  lhsLagL(3,k,:)=lhsLagL(3,k,:)+ &
                     loclhsLag(e,3,aa,b,:)
                  lhsLagL(4,k,:)=lhsLagL(4,k,:)+ &
                     loclhsLag(e,4,aa,b,:)
                  lhsLagL(5,k,:)=lhsLagL(5,k,:)+ &
                     loclhsLag(e,5,aa,b,:)
                  lhsLagL(6,k,:)=lhsLagL(6,k,:)+ &
                     loclhsLag(e,6,aa,b,:)
                  lhsLagL(7,k,:)=lhsLagL(7,k,:)+ &
                     loclhsLag(e,7,aa,b,:)
                  lhsLagL(8,k,:)=lhsLagL(8,k,:)+ &
                     loclhsLag(e,8,aa,b,:)
                  lhsLagL(9,k,:)=lhsLagL(9,k,:)+ &
                     loclhsLag(e,9,aa,b,:)
               enddo
            enddo
         enddo
      endif
!
!.... end
!
	return
	end


      subroutine fillsparseC(iens, EGmass, lhsK, &
                             row,  col)
!
!-----------------------------------------------------------
! This routine fills up the spasely stored LHS mass matrix
!
! Nahid Razmra, Spring 2000. (Sparse Matrix)
!-----------------------------------------------------------
!
!

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      real*8 EGmass(npro,nedof,nedof)
      integer ien(npro,nshl), col(nshg+1), row(nnz*nshg)
      real*8 lhsK(nflow*nflow,nnz_tot)

!
      integer aa, b, c, e, i, k, n, n1
      integer f, g, r, s, t
!
      integer sparseloc

      integer iens(npro,nshl)
!
! prefer to show explicit absolute value needed for cubic modes and
! higher rather than inline abs on pointer as in past versions
! iens is the signed ien array ien is unsigned
!
      ien=abs(iens)
!
!.... Accumulate the lhsK
!
      do e = 1, npro
         do aa = 1, nshl
            i = ien(e,aa)
            c = col(i)
            n = col(i+1) - c
            do b = 1, nshl
               k = sparseloc( row(c), n, ien(e,b) ) + c-1
!
               do f = 1, nflow
                  do g = 1, nflow
                     r = (aa-1)*nflow + f
                     s = (b-1)*nflow + g
                     t = (g-1)*nflow + f

                     lhsK(t,k) = lhsK(t,k) + EGmass(e,r,s)
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!.... end
!
	return
	end

      subroutine fillsparseSclr(iens, xSebe, lhsS, &
                                row,  col)
!
!
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      real*8  xSebe(npro,nshl,nshl)
      integer ien(npro,nshl), col(nshg+1), row(nshg*nnz)
      real*8  lhsS(nnz_tot)
!
      integer aa, b, c, e, i, k, n
!
      integer sparseloc

      integer iens(npro,nshl)
!
! prefer to show explicit absolute value needed for cubic modes and
! higher rather than inline abs on pointer as in past versions
! iens is the signed ien array ien is unsigned
!
      ien=abs(iens)
!
!.... Accumulate the lhs
!
      do e = 1, npro
         do aa = 1, nshl
            i = ien(e,aa)
            c = col(i)
            n = col(i+1) - c
            do b = 1, nshl
               k = sparseloc( row(c), n, ien(e,b) ) + c-1
!
               lhsS(k) = lhsS(k) + xSebe(e,aa,b)
            enddo
         enddo
      enddo
!
!.... end
!
      return
	end

      integer function sparseloc( list, n, target )

!-----------------------------------------------------------
! This function finds the location of the non-zero elements
! of the LHS matrix in the sparsely stored matrix 
! lhsK(nflow*nflow,nnz*numnp)
!
! Nahid Razmara, Spring 2000. 	(Sparse Matrix)
!-----------------------------------------------------------

      integer list(n), n, target
      integer rowvl, rowvh, rowv

!
!.... Initialize
!
      rowvl = 1
      rowvh = n + 1
!
!.... do a binary search
!
100   if ( rowvh-rowvl .gt. 1 ) then
         rowv = ( rowvh + rowvl ) / 2
         if ( list(rowv) .gt. target ) then
            rowvh = rowv
         else
            rowvl = rowv
         endif
         goto 100
      endif
!
!.... return
!
      sparseloc = rowvl
!
      return
      end
