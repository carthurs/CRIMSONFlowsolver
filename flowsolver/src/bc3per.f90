        subroutine bc3per (iBC,  res, iper, ilwork,nQs)
!
!----------------------------------------------------------------------
!
! This routine satisfies the BC of the periodic nodes after Ap product
!
! input:
!  iBC   (nshg)        : Boundary Condition Code
!  iper  (nshg)        : partners of periodic nodes
!  res   (nshg,nQs)   : residual before BC is applied
!
! output:
!  res   (nshg,nQs)   : residual after satisfaction of BC
!  
!
! Kenneth Jansen,  Winter 1998. 
!----------------------------------------------------------------------
!
  use periodicity  ! this gives you rcount(1:nshg) (real*8)
  use phcommonvars
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension iBC(nshg), & 
                  res(nshg,nQs),           ilwork(nlwork), & 
                  iper(nshg)
!
!.... local periodic (and axisymmetric) boundary conditions (no communications)
!
           do j = 1,nshg
              if ((btest(iBC(j),10)) .or. (btest(iBC(j),12))) then
                 i = iper(j)
                 res(i,:) = res(i,:) + res(j,:)
                 res(j,:) = zero
              endif
           enddo


        if(numpe.gt.1) then
!
!.... nodes treated on another processor are eliminated
!     
           numtask = ilwork(1)
           itkbeg = 1
           
           do itask = 1, numtask
              
              iacc   = ilwork (itkbeg + 2)
              numseg = ilwork (itkbeg + 4)
              
              if (iacc .eq. 0) then
                 do is = 1,numseg
                    isgbeg = ilwork (itkbeg + 3 + 2*is)
                    lenseg = ilwork (itkbeg + 4 + 2*is)
                    isgend = isgbeg + lenseg - 1
                    res(isgbeg:isgend,:) = zero
                 enddo
              endif
              
              itkbeg = itkbeg + 4 + 2*numseg
              
           enddo
        endif
!
!.... return
!
        return
        end





