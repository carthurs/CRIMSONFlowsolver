        module periodicity
       
        real*8, allocatable :: rcount(:)
        end module
        
       subroutine Dperiodicity
       
       use periodicity
       if (allocated(rcount)) then
         deallocate(rcount)
       endif
       
       return
       end


       subroutine setper (nshg)
       
       use periodicity

       if (.not.allocated(rcount)) then
         allocate (rcount(nshg))
       endif
       
       return
       end

       subroutine perprep (iBC, iper,nshg)
       
       use periodicity

       dimension iBC(nshg), &
                 iper(nshg)       
  
!
!..... calculate the inverse of the number of slaves + 1
!
       one=1.00000000000
       rcount=one
       do j = 1,nshg
          if (btest(iBC(j),10)) then
             i = iper(j)
             rcount(i) = rcount(i) + one
          endif
       enddo
       do k=1,nshg
          if(rcount(k).ne.one) then
             rcount(k)=one/rcount(k)
          endif
       enddo
      
!
!.... return
!
       return
       end
