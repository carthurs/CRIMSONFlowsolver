      subroutine sparseCG (rhsorig, trx, lhsM, row, col, iper,
     &                     ilwork, iBC, BC)
c
c------------------------------------------------------------------------
c
c This subroutine uses Conjugate Gradient,
c to solve the system of equations.
c
c Farzin Shakib,  Summer 1987.
c------------------------------------------------------------------------
c
      include "common.h"
c
      dimension rhsorig(nshg), trx(nshg)
c
      dimension d(nshg),     p(nshg),
     &          q(nshg),     ilwork(nlwork),
     &		    Dinv(nshg),  rhs(nshg),
     &          pp(nshg),    
     &          iBC(nshg),
     &          BC(nshg)


      integer   row(nnz_tot_bdy),        col(nshg+1)
      integer   iBCdumb(nshg),         iper(nshg)

      real*8    BCdumb(nshg,ndofBC)

      real*8	lhsM(9,nnz_tot_bdy)
c
      data nitercf / 100 /
c
c
      BCdumb  = one
      iBCdumb = 1
c
      rhs(:)=rhsorig(:)
c
c Initialize. We have a system Ax=b We have made A as
c well conditionedand as we're willing to go.  Since
c we used left preconditioning (on the old A and old b),
c we don't need to do any back-preconditioning later.
c
      rr = 0
      do n = 1, nshg
         rr  = rr + rhs(n)*rhs(n)
      enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
      dotTot=zero
      call drvAllreduce(rr,dotTot,1)
      rr=dotTot
      rr0 = rr
c
      trx(:) = zero ! x_{0}=0
c                   ! r_{0}=b
c
c                   ! beta_{1}=0
      p(:) = rhs(:) ! p_{1}=r_{0}=b
c
c.... Iterate
c        
      do iter = 1, nitercf      ! 15  ! nitercf
c
c.... calculate alpha
c
         pp=p   ! make a vector that we can copy masters onto slaves
                ! and thereby make ready for an accurate Ap product

         call commOut(pp, ilwork, 1,iper,iBCdumb,BCdumb)  !slaves= master

         call fLesSparseApSclr(	col,	row,	lhsM,	
     &		            		pp,	    q,	    nshg,
     &                          nnz)

         call commIn(q, ilwork, 1,iper,iBC,BC) ! masters=masters+slaves
					                           ! slaves zeroed

         pap = 0
         do  n = 1, nshg
            pap = pap + p(n) * q(n)
         enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
         dotTot=zero
         call drvAllreduce(pap,dotTot,1)
         pap=dotTot
         alpha = rr / pap 
c
c.... calculate the next guess
c
         trx(:) = trx(:) + alpha * p(:)
c
c.... calculate the new residual
c
c
         rhs(:) = rhs(:) - alpha * q(:)
         tmp = rr
         rr = 0
         do n = 1, nshg
            rr = rr + rhs(n)*rhs(n)
         enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
         dotTot=zero
         call drvAllreduce(rr,dotTot,1)
         rr=dotTot
c
c.... check for convergence
c
         if(rr.lt.100.*epsM**2) goto 6000
c
c.... calculate a new search direction
c
         beta = rr / tmp
         p(:) = rhs(:) + beta * p(:)
c
c.... end of iteration
c
      enddo
c
c.... if converged
c
6000  continue

c need a commu(out) on solution (TRX) to get slaves the correct solution AND
c on processor slaves = on processor masters

      if(numpe>1) call commu (trx, ilwork, 1, 'out ')
	  trx(:) = trx(iper(:))

c
      write(*,9000) iter, rr / rr0
c	write(16,9000) iter, rr / rr0
c
c.... return
c
      return
9000  format(20x,'  number of iterations:', i10,/,
     &       20x,'    residual reduction:', 2x,e10.2)
      end
