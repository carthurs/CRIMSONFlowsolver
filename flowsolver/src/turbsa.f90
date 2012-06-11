c
c  Copyright (c) 2000-2007, Stanford University, 
c     Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
c     Charles A. Taylor (see SimVascular Acknowledgements file 
c     for additional contributors to the source code).
c
c  All rights reserved.
c
c  Redistribution and use in source and binary forms, with or without 
c  modification, are permitted provided that the following conditions 
c  are met:
c
c  Redistributions of source code must retain the above copyright notice,
c  this list of conditions and the following disclaimer. 
c  Redistributions in binary form must reproduce the above copyright 
c  notice, this list of conditions and the following disclaimer in the 
c  documentation and/or other materials provided with the distribution. 
c  Neither the name of the Stanford University or Rensselaer Polytechnic
c  Institute nor the names of its contributors may be used to endorse or
c  promote products derived from this software without specific prior 
c  written permission.
c
c  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
c  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
c  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
c  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
c  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
c  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
c  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
c  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
c  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
c  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
c  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
c  DAMAGE.
c
c
c-----------------------------------------------------------------------
c
c     Spalart-Allmaras turbulence model constants 
c
c-----------------------------------------------------------------------
      module turbSA

      real*8, allocatable ::  d2wall(:)
      real*8, allocatable ::  wnrm(:,:)
      integer, allocatable :: otwn(:)
      real*8  saCb1, saCb2, saCw1, saCw2, saCw3, saCv1, saSigma,
     &        saKappa, saKappaP2Inv, saCv2P3, saCw3P6, saSigmaInv
      integer, allocatable :: sidmapg(:) ! list of all surfID's, low to high
      
      parameter ( 
     &     saCb1        = 0.1355d0,
     &     saCb2        = 0.622d0,
     &     saCw1        = 3.239067817d0,
     &     saCw2        = 0.3d0,
     &     saCw3        = 2.0d0,
     &     saKappa      = 0.41d0,
     &     saSigma      = 0.666666666666666667d0,
     &     saCv1        = 7.1d0,
     &     saKappaP2Inv = 5.94883997620464d0,
     &     saCv1P3      = 3.579109999999999d+02,
     &     saCw3P6      = 64.0d0,
     &     saSigmaInv   = 1.50d0
     &     )

      
      end module

c-----------------------------------------------------------------------
c
c     Initialize: compute the distance to the nearest wall for 
c     each vertex in the mesh.
c
c-----------------------------------------------------------------------
      subroutine initTurb( x )
      
      use     pointer_data
      use     turbSA
      include "common.h"
      include "mpif.h"
      
      character*20 fname1,  fmt1
      real*8   x(numnp,nsd)
      integer  nwall(numpe),      idisp(numpe)
      character*5  cname      
      real*8, allocatable :: xwi(:,:,:), xw(:,:,:)
      
      allocate ( d2wall(numnp) )
c
c.... Find out if d2wall is already calculated
c
      open(unit=72,file='distcalc.dat',status='old')
      read(72,*) idistcalc
c      idistcalc=1
      close(72)
      if(idistcalc.eq.1) then   ! calculate the distance
c
c.... Count the welts (wall-elements)
c
         nwalli=0
         do iblk = 1, nelblb    ! loop over boundary elt blocks
            npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
            do j = 1, npro
               if(btest(miBCB(iblk)%p(j,1),4)) nwalli=nwalli+1
            enddo
         enddo
c
c.... Create wallnode-coord list for welts on processor
c
         if (nwalli.eq.0) nwalli=1 !  patch for mpi's lack of imagination
         allocate (xwi(nwalli,nenb+1,nsd))
         xwi = 1.0d32
         xwi(:,nenb+1,:)=zero
         nwalli = 0
         do iblk = 1, nelblb    ! loop over boundary elt blocks
c
            iel    = lcblkb(1,iblk)
            nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
            npro   = lcblkb(1,iblk+1) - iel 
c
            do j = 1, npro      ! loop over belts in this blk
               if(btest(miBCB(iblk)%p(j,1),4)) then
                  nwalli=nwalli+1
c assemble local coord list
                  do node = 1, nenbl
                     xwi(nwalli,node,1:3)=x(mienb(iblk)%p(j,node),:)
                  enddo
c put the centroid coordinates in the last slot
                  do node = 1, nenbl
                     xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)
     &                    +xwi(nwalli,node,:)
                  enddo
                  xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)/nenbl
c
               endif
            enddo               ! loop over belts in this blk
c
         enddo                  ! loop over boundary elt blocks
c
         if (nwalli.eq.0) xwi=1.0e32 ! fix for mpi's lack of imagination
         if (nwalli.eq.0) nwalli=1 !  patch for mpi's lack of imagination
c
c  Pool "number of welts" info from all processors
c
cMPI_ALLGATHER(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm) 
c[ IN sendbuf] starting address of send buffer (choice) 
c[ IN sendcount] number of elements in send buffer (integer) 
c[ IN sendtype] data type of send buffer elements (handle) 
c[ OUT recvbuf] address of receive buffer (choice) 
c[ IN recvcount] number of elements received from any process (integer) 
c[ IN recvtype] data type of receive buffer elements (handle) 
c[ IN comm] communicator (handle)
         if (numpe.gt.1) then   ! multiple processors
c write the number of wall elts on the jth processor to slot j of nwall
            call MPI_ALLGATHER(nwalli,1,MPI_INTEGER,nwall,1,
     &          MPI_INTEGER,MPI_COMM_WORLD,ierr)
c count up the total number of wall elts among all processes
            nwallt=0
            do j=1,numpe
               nwallt=nwallt+nwall(j)
            enddo
         else                   ! single processor
c the local information is the global information for single-processor
            nwall=nwalli
            nwallt=nwalli
         endif                  ! if-else for multiple processors
c
c  Make all-processor wallnode-coord collage
c
         allocate (xw(nwallt,nenb+1,nsd))
         if (numpe.gt.1) then   ! multiple processors
c we will gather coordinates from local on-proc sets to a global set
c we will stack each processor's coordinate list atop that of the
c previous processor.  If the coordinate list for processor i is
c called xwi, then our global coordinate list xw will look like this:
c ---------------------------------------------------------------
c | xw1            | xw2                | xw3        |   ...    |
c ---------------------------------------------------------------
c  <---nwall(1)---> <-----nwall(2)-----> <-nwall(3)->
c  <------------------------nwallt-----------------------...---->
c To accomplish this with MPI, we use MPI_ALLGATHERV, summarized as:
cMPI_ALLGATHERV(sendbuf,sendcount,sendtype,recvbuf,recvcount,disp,recvtype,comm) 
c[ IN sendbuf] starting address of send buffer (choice) 
c[ IN sendcount] number of elements in send buffer (integer) 
c[ IN sendtype] data type of send buffer elements (handle) 
c[ OUT recvbuf] address of receive buffer (choice) 
c[ IN recvcount] number of elements received from any process (integer) 
c[ IN disp] displacement array
c[ IN recvtype] data type of receive buffer elements (handle) 
c[ IN comm] communicator (handle)
c The displacement array disp is an array of integers in which the jth
c entry indicates which slot of xw marks beginning of xwj
c So, first we will build this displacement array
            idisp(:)=0 ! starting with zero, since MPI likes C-numbering
            do j=2,numpe
               idisp(j)=idisp(j-1)+nwall(j-1) ! see diagram above
            enddo
c Now, we gather the data one slice at a time (1:nwalli)
            do j=1,nenb+1
               do k=1,nsd
                  call MPI_ALLGATHERV(xwi(:,j,k),nwalli,
     &                 MPI_DOUBLE_PRECISION,xw(:,j,k),nwall,idisp,
     &                 MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
               enddo
            enddo
         else                   ! single-processor
c global data is local data in single processor case
            xw=xwi
         endif

c
c  For each point, loop over wall nodes and calculate distance; 
c  save the distance in this node's position of d2wall if it's 
c  shorter than currently stored distance
c
         d2wall=1.0e32
         do i=1,numnp
            do j=1, nwallt
               do k=1,nenb+1
                  distance =  ( x(i,1) - xw(j,k,1) )**2
     &                 +( x(i,2) - xw(j,k,2) )**2
     &                 +( x(i,3) - xw(j,k,3) )**2
                  if ( d2wall(i).gt.distance ) d2wall(i) = distance
               enddo
            enddo
         enddo
         d2wall=sqrt(d2wall)
c
         deallocate(xwi)
         deallocate(xw)
c
c.... write d2wall to a file so we don't have to do this again
c
         write (fmt1,"('(''d2wall.'',i',i1,',1x)')") 1
         write (fname1,fmt1) 0
         fname1 = trim(fname1) // cname(myrank+1)
         open (unit=72, file=fname1, status='unknown',
     &                                    form='unformatted', err=996)

         write (72) d2wall
         close (72)
c         write(*,*) "make sure to: echo 0 > distcalc.dat"
c         call MPI_BARRIER (MPI_COMM_WORLD,ierr)
c         call error ('distcalc','complete',0)
      endif
      if (idistcalc.eq.0) then ! d2wall is already done, don't calculate it
        write (fmt1,"('(''d2wall.'',i',i1,',1x)')") 1
        write (fname1,fmt1) 0
        fname1 = trim(fname1) // cname(myrank+1)
        write (*,*) 'Reading dist file : ', fname1

        open (unit=72, file=fname1, status='old',
     &        form='unformatted', err=995)

        read (72) d2wall
        close (72)
      endif

      return
995     call error ('d2wall  ','opening ', 72)
996     call error ('d2wall  ','opening ', 72)

      end
