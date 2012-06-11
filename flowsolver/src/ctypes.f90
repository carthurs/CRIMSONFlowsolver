      subroutine ctypes (ilwork)

  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

      parameter (maxseg = 10000)

      integer  sizeofdouble

      dimension ilwork(nlwork)
      dimension isbegin(maxseg), lenseg(maxseg), ioffset(maxseg)

      CALL MPI_TYPE_EXTENT (MPI_DOUBLE_PRECISION,sizeofdouble,ierr)
      lstride = nshg * sizeofdouble
!
!.... maxfront is a common variable being set in this routine
!
      maxfront = 0  
      numtask = ilwork (1)
      itkbeg  = 1

      if (numtask .gt. maxtask)  & 
        call error('ctypes  ','numtask ',numtask)

      nshg0 = nshg

      do itask = 1,numtask !{
!
!.... iacc = 0 ==> this task is a send
!          = 1 ==> this task is a recieve
!
        iacc   = ilwork (itkbeg + 2)
!
!.... numseg : number of data segments to be communicated
!
        numseg = ilwork (itkbeg + 4)
!
!.... adjust the number of the other processor, since processors
!     are numbered here starting from 0, not 1.
!
        ilwork (itkbeg + 3) = ilwork (itkbeg + 3) - 1 
        if (numseg .gt. maxseg)  & 
          call error('ctypes  ','numseg  ',numseg )
!
!.... lfront = total number of nodes involved in this task
!     
        lfront = 0
        do is = 1,numseg ! {
!
!.... isbegin(is): starting node number for each segment
!
          isbegin (is) = ilwork (itkbeg + 3 + 2*is)
!
!.... lenseg(is): length of each segment (number of nodes)
!
          lenseg  (is) = ilwork (itkbeg + 4 + 2*is)
!
!.... increment the total node counter
!
          lfront       = lfront + lenseg(is)
!
!.... nshg0: number of nodes to be assembled on this processor,
!             i.e. subtract the number of nodes which will be 
!             sent to another processor.
!
	  if (iacc .eq. 0) nshg0 = nshg0 - lenseg(is)
        enddo ! }
!
!.... maxfront: number of nodes which will be communicated, including
!               all segments. Note that after the loop over tasks
!               is complete, maxfront will contain the maximum number
!               of nodes for any of the tasks.
!
        maxfront = MAX(maxfront,lfront)
!
!.... ioffset: array offset from the first node in the first segment
!
        ioffset(1:numseg) = isbegin(1:numseg) - isbegin(1)
!
!.... now set up the MPI data types which will be used in commu.f.
!     These data types represent the indexed sets that will be sent
!     and recieved.
! 
!
!.... the following call to MPI_TYPE_INDEXED will create a new data
!     type which will represent the blocks of data we wish to transfer
!     for this task. A handle to the new type is returned 
!     (sevsegtype(itask,1)). This data type describes the blocks of
!     data to be transferred in terms of segments.
!     Input to this routine:
!          numseg: number of segments in this task
!          lenseg: length of each segment (number of nodes)
!          ioffset: where to begin each block with respect to the
!                   first segment
!          MPI_DOUBLE_PRECISION: type to set for each of the blocks
!
        call MPI_TYPE_INDEXED (numseg, lenseg, ioffset, & 
                        MPI_DOUBLE_PRECISION, sevsegtype(itask,1), ierr)
!
!.... now create a new data type for each of the types of arrays we 
!     may wish to communicate with. For example ndof will be used when
!     communicating the residual vector. Each one of these is derived
!     from the first data type defined above, sevsegtype(itask,1).
!
        call MPI_TYPE_HVECTOR(nsd,    1, lstride, sevsegtype(itask,1), & 
                                              sevsegtype(itask,2), ierr)
!
        call MPI_TYPE_HVECTOR(ndof,   1, lstride, sevsegtype(itask,1), & 
                                              sevsegtype(itask,3), ierr)
!
        call MPI_TYPE_HVECTOR(nflow*nflow,1, lstride, & 
                          sevsegtype(itask,1),sevsegtype(itask,4), ierr)
        call MPI_TYPE_HVECTOR((nflow-1)*nsd,1,lstride, & 
                          sevsegtype(itask,1),sevsegtype(itask,5), ierr)
       call MPI_TYPE_HVECTOR(nflow,1,lstride,sevsegtype(itask,1), & 
                                              sevsegtype(itask,6), ierr)
       call MPI_TYPE_HVECTOR(24,1,lstride,sevsegtype(itask,1), & 
                                              sevsegtype(itask,7), ierr)
       call MPI_TYPE_HVECTOR(9,1,lstride,sevsegtype(itask,1), & 
                                              sevsegtype(itask,8), ierr)
       call MPI_TYPE_HVECTOR(11,1,lstride,sevsegtype(itask,1), & 
                                              sevsegtype(itask,9), ierr)
       call MPI_TYPE_HVECTOR(7,1,lstride,sevsegtype(itask,1), & 
                                         sevsegtype(itask,10), ierr)
       call MPI_TYPE_HVECTOR(33,1,lstride,sevsegtype(itask,1), & 
                                         sevsegtype(itask,11), ierr)
       call MPI_TYPE_HVECTOR(22,1,lstride,sevsegtype(itask,1), & 
                                         sevsegtype(itask,12), ierr)
       call MPI_TYPE_HVECTOR(16,1,lstride,sevsegtype(itask,1), & 
                                         sevsegtype(itask,13), ierr)
       call MPI_TYPE_HVECTOR(10,1,lstride,sevsegtype(itask,1), & 
                                         sevsegtype(itask,14), ierr)
       call MPI_TYPE_HVECTOR(nflow*nsd,1,lstride,sevsegtype(itask,1), & 
                                         sevsegtype(itask,15), ierr)
!
!
!.... now this must be done to make MPI recognize each of the data
!     types that were just defined
!
        do kdof = 1,15
          call MPI_TYPE_COMMIT (sevsegtype(itask,kdof), ierr)
        enddo
!
!.... set the counter to the index in ilwork where the next task
!     begins
!

        itkbeg = itkbeg + 4 + 2*numseg
!
!.... end loop over tasks
!
      enddo ! }

      return
      end

subroutine Dctypes(ilwork) ! {
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
  include "mpif.h"
  !include "auxmpi.h"

  dimension ilwork(nlwork)
  
  numtask = ilwork(1)
 
  do itask = 1,numtask 
    do kdof = 1,15
      call MPI_TYPE_FREE(sevsegtype(itask,kdof), ierr)
    enddo
  enddo 
  
end !}
