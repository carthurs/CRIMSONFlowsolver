      subroutine commu (global, ilwork,  n,  code)
!---------------------------------------------------------------------
! 
! This subroutine is responsible for interprocessor communication of
! the residual and solution vectors.
!
! input:   
!     global(nshg,n): global vector to be communicated. Note that
!                      this vector is local to the processor, (i.e.
!                      not distributed across processors)
!     ilwork(nlwork):  this is the local interprocessor work array.
!                      This array is local to the processor, (i.e.
!                      each processor has a unique ilwork array.
!     n:               second dimension of the array to be communicated
!     code:            = 'in' for communicating with the residual
!                      = 'out' for cummunicating the solution 
!
!---------------------------------------------------------------------
!
! The array ilwork describes the details of the communications. 
! Each communication step (call of this routine) consists of a 
! sequence of "tasks", where a task is defined as a communication 
! between two processors where data is exchanged. This would imply 
! that for a given processor, there will be as many tasks as there
! are processors with which it must communicate. Details of the 
! ilwork array appear below.
!
!---------------------------------------------------------------------
!
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"
      integer status(MPI_STATUS_SIZE)
      integer stat(MPI_STATUS_SIZE, 2*numpe), req(2*numpe)

      dimension global(nshg,n), & 
                rtemp(maxfront*n,maxtask), & 
                ilwork(nlwork)
 
      character*3 code


      if (code .ne. 'in ' .and. code .ne. 'out')  & 
        call error ('commu   ','code    ',0)

      if     (n .eq. 1)      then        ! like a scalar
        kdof = 1
      elseif (n .eq. nsd)    then        ! like the normal vectors
        kdof = 2
      elseif (n .eq. ndof)   then        ! res, y, ac, krylov vectors....
        kdof = 3
      elseif (n .eq. nflow*nflow) then     ! bdiag
        kdof = 4
      elseif (n .eq. (nflow-1)*nsd) then  ! qres
        kdof = 5
      elseif (n .eq. nflow) then
        kdof = 6
      elseif (n .eq. 24 ) then
        kdof = 7
      elseif (n .eq. 9) then
        kdof = 8
      elseif (n .eq. 11 ) then
        kdof = 9
      elseif (n .eq. 7 ) then
        kdof = 10
      elseif (n .eq. 33 ) then
         kdof = 11 
      elseif (n .eq. 22 ) then
         kdof = 12
      elseif (n .eq. 16 ) then
         kdof = 13
      elseif (n .eq. 10 ) then
         kdof = 14
       elseif (n .eq. nflow*nsd ) then   !surface tension + qres
         kdof = 15 
      else
        call error ('commu   ','n       ',n)
      endif

!... Note that when adding another kdof to the above set, we must
!... also make changes in ctypes.f and auxmpi.h

!---------------------------------------------------------------------
!  ilwork(1): number of tasks
!
!  The following information is contained in ilwork for each task:
!     itag: tag of the communication
!     iacc: == 0 if task is a send
!           == 1 if task is a recieve
!     iother: rank of processor with which this communication occurs
!     numseg: number of data "segments" to be sent or recieved. A 
!             segment is defined as a continuous section of the global
!             vector to be communicated, (i.e. a group of nodes (or,
!             rather, "shape function coefficients") which occur 
!             sequentially in the array global(nshg,n)).
!     isbeg:  location of the first segment in the array owned by the
!             current processor.
!
! The two types of communication are 'in', where the residual is being
! communicated, and 'out', where the solution is being communicated.
! Note that when the type is 'out', senders recieve and recievers send.
!
! The following comment pertains to a communication of type 'in':
!
!     If the task is a send, then all of the numseg segments are
!     sent with a single call to MPI_SEND. Where these segments live in 
!     the array is built into the array sevsegtype, which is a common 
!     array constructed in the subroutine "ctypes.f". In other words,
!     sevsegtype is a data type that describes the indices of the blocks
!     to be sent, in terms of there beginning index, and the length of 
!     each segment. Using this, we can make a single send to take care of
!     all the segments for this task. 
!      
!     If the task is a recieve, then once the vector is recieved, the
!     recieved segments must be added to the correct locations in the
!     current array. These locations are described in ilwork as the
!     beginning position, then the length of the segment.
!     
!---------------------------------------------------------------------
      numtask = ilwork(1)
      
      itkbeg = 1
      m = 0
      idl=0

      DO itask = 1, numtask
        m      = m + 1
        itag   = ilwork (itkbeg + 1)
        iacc   = ilwork (itkbeg + 2)
        iother = ilwork (itkbeg + 3)
        numseg = ilwork (itkbeg + 4)
        isgbeg = ilwork (itkbeg + 5)
!
!.... if iacc == 0, then this task is a send.
!     slave
!
        if (iacc .EQ. 0) then  
!
!.... residual communication
!
          if (code .eq. 'in ') & 
            call MPI_ISEND(global(isgbeg, 1), 1, sevsegtype(itask,kdof),  & 
                           iother, itag, INEWCOMM, req(m), ierr)
!
!.... solution communication
!
          if (code .eq. 'out') then
            call MPI_IRECV(global(isgbeg, 1), 1, sevsegtype(itask,kdof),  & 
                           iother, itag, INEWCOMM, req(m), ierr)
!            call MPI_RECV(global(isgbeg,1), 1, sevsegtype(itask,kdof),
!     &                    iother, itag, INEWCOMM, status, ierr)
          endif
!
!.... if iacc == 1, then this task is a recieve.
!     master
!
        else
          if (code .eq. 'in ') then
!
!.... determine the number of total number of nodes involved in this
!     communication (lfront), including all segments
!
            lfront = 0
            do is = 1,numseg
              lenseg = ilwork (itkbeg + 4 + 2*is)
              lfront = lfront + lenseg
            enddo
!
!.... recieve all segments for this task in a single step
!
            idl=idl+1 ! stands for i Do Later, the number to fix later
            call MPI_IRECV(rtemp(1,idl), lfront*n, MPI_DOUBLE_PRECISION,  & 
                           iother, itag, INEWCOMM, req(m), ierr)
          endif
          if (code .eq. 'out') then
            call MPI_ISEND(global(isgbeg, 1), 1, sevsegtype(itask,kdof),  & 
                           iother, itag, INEWCOMM, req(m), ierr)
          endif
        endif

        itkbeg = itkbeg + 4 + 2*numseg

      enddo   !! end tasks loop

      call MPI_WAITALL(m, req, stat, ierr)

!
!     Stuff added below is a delayed assembly of that which was communicated
!     above but due to the switch to non-blocking receivves could not be
!     assembled until after the waitall.  Only necessary for commu "in"
!

      if(code .eq. 'in ') then
         itkbeg=1
         jdl=0
         do j=1,numtask         ! time to do all the segments that needed to be
                                ! assembled into the global vector

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)
            isgbeg = ilwork (itkbeg + 5)
            if(iacc.eq.1) then
               jdl=jdl+1  ! keep track of order of rtemp's
!
!... add the recieved data to the global array on the current processor.
!    Note that this involves splitting up the chunk of recieved data
!    into its correct segment locations for the current processor.
!
               itemp = 1
               do idof = 1,n
                  do is = 1,numseg
                 isgbeg = ilwork (itkbeg + 3 + 2*is)
                 lenseg = ilwork (itkbeg + 4 + 2*is)
                 isgend = isgbeg + lenseg - 1
                 global(isgbeg:isgend,idof) = global(isgbeg:isgend,idof) & 
                                      + rtemp (itemp:itemp+lenseg-1,jdl)
                 itemp = itemp + lenseg
                  enddo
               enddo
            endif ! end of receive (iacc=1)
            itkbeg = itkbeg + 4 + 2*numseg
         enddo
      endif  ! commu "in"
      return
      end



