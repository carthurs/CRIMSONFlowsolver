        subroutine restar (code,  qp, acp)
!
!----------------------------------------------------------------------
!  This routine is the restart option.
!
! input:
!  code                 : restart option on the primitive variables
!                          eq. 'in  ', read from [restar.inp]
!                          eq. 'out ', write to  [restar.out]
!
! input or output:
!  q     (nshg,ndof)   : the variables to be read/written
!
!
! Farzin Shakib, Summer 1985.
!----------------------------------------------------------------------
!
        use GlobalArrays          ! used to access acold
        use readarrays          ! used to access qold
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
!
        character*4 code
        character*8 mach2
        character*20 fname1,  fmt1
        character*5  cname

!
        dimension qp(nshg,ndof),acp(nshg,ndof)
        dimension qw(nshg,ndof),acw(nshg,ndof)
!
!.... -------------------------->  'in (unused now since initialization of y and ac moved to readnblk) '  <---------------------------
!
        !if (code .eq. 'in  ') then
!
! incompressible orders velocity, pressure, temperature unlike compressible
! which is what we have our files set up for
!
        !  qp(:,1:3)=qold(:,2:4)
        !  qp(:,4)=qold(:,1)
        !  if(ndof.gt.4)  qp(:,5:ndof)=qold(:,5:ndof)

        !  acp(:,1:3)=acold(:,2:4)
        !  acp(:,4)=acold(:,1)
        !  if(ndof.gt.4)  acp(:,5:ndof)=acold(:,5:ndof)

        !  deallocate(qold)
        !  return
        !endif
!
!.... -------------------------->  'out '  <---------------------------
!
        if (code .eq. 'out ') then
           qw(:,2:4) = qp(:,1:3)
           qw(:,1)   = qp(:,4)
           if(ndof.gt.4) qw(:,5:ndof)   = qp(:,5:ndof)
! 
           acw(:,2:4) = acp(:,1:3)
           acw(:,1)   = acp(:,4)
           if(ndof.gt.4) acw(:,5:ndof)   = acp(:,5:ndof)
!           
           call write_restart(myrank, currentTimestepIndex, nshg, ndof, & 
                qw, acw)

           if (myrank.eq.master) then 
              open(unit=72,file='numstart.dat',status='old')
              write(72,*) currentTimestepIndex
              close(72)

           endif
!$$$          ttim(75) = ttim(75) + tmr()
           return
        endif
!
!.... ---------------------->  Error Handling  <-----------------------
!
!.... Error handling
!
        call error ('restar  ',code//'    ',0)
!
!.... file error handling
!
995     call error ('restar  ','opening ', irstin)
996     call error ('restar  ','opening ', irstou)
!
!.... end
!
        end
