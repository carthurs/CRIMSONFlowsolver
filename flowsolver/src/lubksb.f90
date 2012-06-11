      subroutine lubksb(aa,n,np,indx,bb)
!---------------------------------------------------------------------
!
! LU back substitution routine.
!
!---------------------------------------------------------------------
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      dimension indx(n),aa(np,np),bb(n)
      
      ii=0
      do i=1,n
         ll=indx(i)
         sum=bb(ll)
         bb(ll)=bb(i)
         if (ii.ne.0)then
            do j=ii,i-1
               sum=sum-aa(i,j)*bb(j)
            enddo
         else if (sum.ne.0.0) then
            ii=i
         endif
         bb(i)=sum
      enddo
      do i=n,1,-1
         sum=bb(i)
         do j=i+1,n
            sum=sum-aa(i,j)*bb(j)
         enddo
         bb(i)=sum/aa(i,i)
      enddo
      
      return
      end
      
            
    
