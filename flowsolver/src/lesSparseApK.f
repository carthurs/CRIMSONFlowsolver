C============================================================================
C
C "fLesSparseApK":
C
C============================================================================

      subroutine fLesSparseApK(	col, row, kLhs,
     &                          p,   q,   nNodes,
     &                          nnz_tot_bdy )
c
c.... Data declaration
c
      include "common.h"
	  
      integer nNodes, nnz_tot
      integer col(nNodes+1), row(nnz_tot_bdy)
      real*8  kLhs(9,nnz_tot_bdy)
      real*8  p(nNodes,3), q(nNodes,3)
c
      real*8  tmp1, tmp2, tmp3, pisave
      integer i,j,k
c
c.... clear the vector
c
      do i = 1, nNodes
	  
         q(i,1) = 0
         q(i,2) = 0
         q(i,3) = 0
		 
      enddo
c
c.... Do an AP product
c
      do i = 1, nNodes
c
         tmp1 = 0
         tmp2 = 0
         tmp3 = 0
c         pisave   = p(i,4)
cdir$ ivdep
         do k = col(i), col(i+1)-1
		 
            j = row(k) 
            tmp1 = tmp1
     &            + kLhs(1,k) * p(j,1)
     &            + kLhs(4,k) * p(j,2)
     &            + kLhs(7,k) * p(j,3)
            tmp2 = tmp2
     1            + kLhs(2,k) * p(j,1)
     2            + kLhs(5,k) * p(j,2)
     3            + kLhs(8,k) * p(j,3)
            tmp3 = tmp3
     1            + kLhs(3,k) * p(j,1)
     2            + kLhs(6,k) * p(j,2)
     3            + kLhs(9,k) * p(j,3)
c
         enddo
		
         q(i,1) = q(i,1) + tmp1
         q(i,2) = q(i,2) + tmp2
         q(i,3) = q(i,3) + tmp3
		 
      enddo

c.... end
c
      return
      end