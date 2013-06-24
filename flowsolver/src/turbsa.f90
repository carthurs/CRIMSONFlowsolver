!-----------------------------------------------------------------------
!
!     Spalart-Allmaras turbulence model constants 
!
!-----------------------------------------------------------------------
      module turbSA

      real*8, allocatable ::  d2wall(:)
      real*8, allocatable ::  wnrm(:,:)
      integer, allocatable :: otwn(:)
      real*8  saCb1, saCb2, saCw1, saCw2, saCw3, saCv1, saSigma, &
              saKappa, saKappaP2Inv, saCv2P3, saCw3P6, saSigmaInv
      integer, allocatable :: sidmapg(:) ! list of all surfID's, low to high
      
      parameter (  &
           saCb1        = 0.1355d0, &
           saCb2        = 0.622d0, &
           saCw1        = 3.239067817d0, &
           saCw2        = 0.3d0, &
           saCw3        = 2.0d0, &
           saKappa      = 0.41d0, &
           saSigma      = 0.666666666666666667d0, &
           saCv1        = 7.1d0, &
           saKappaP2Inv = 5.94883997620464d0, &
           saCv1P3      = 3.579109999999999d+02, &
           saCw3P6      = 64.0d0, &
           saSigmaInv   = 1.50d0 &
           )

      
      end module
      
      subroutine DturbSA
      use turbSA
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      if (allocated(wnrm)) deallocate (wnrm)
      if (allocated(otwn)) deallocate (otwn)
      if(iRANS.lt.0 .and. allocated(d2wall)) deallocate (d2wall)
      if(allocated(sidmapg)) deallocate (sidmapg)
      return
      end

!-----------------------------------------------------------------------
!
!     Initialize: compute the distance to the nearest wall for 
!     each vertex in the mesh.
!
!     Michael Yaworski (fall 1998)
!
!-----------------------------------------------------------------------
      subroutine initTurb( x )
      
      use     pointer_data
      use     turbSA
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      
      character*20 fname1,  fmt1
      real*8   x(numnp,nsd)
      integer  nwall(numpe),      idisp(numpe)
      character*5  cname      
      real*8, allocatable :: xwi(:,:,:), xw(:,:,:)
      
      allocate ( d2wall(numnp) )
!
!.... Find out if d2wall is already calculated
!
      open(unit=72,file='distcalc.dat',status='old')
      read(72,*) idistcalc
!      idistcalc=1
      close(72)
      if(idistcalc.eq.1) then   ! calculate the distance
!
!.... Count the welts (wall-elements)
!
         nwalli=0
         do iblk = 1, nelblb    ! loop over boundary elt blocks
            npro = lcblkb(1,iblk+1) - lcblkb(1,iblk)
            do j = 1, npro
               if(btest(miBCB(iblk)%p(j,1),4)) nwalli=nwalli+1
            enddo
         enddo
!
!.... Create wallnode-coord list for welts on processor
!
         if (nwalli.eq.0) nwalli=1 !  patch for mpi's lack of imagination
         allocate (xwi(nwalli,nenb+1,nsd))
         xwi = 1.0d32
         xwi(:,nenb+1,:)=zero
         nwalli = 0
         do iblk = 1, nelblb    ! loop over boundary elt blocks
!
            iel    = lcblkb(1,iblk)
            nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
            npro   = lcblkb(1,iblk+1) - iel 
!
            do j = 1, npro      ! loop over belts in this blk
               if(btest(miBCB(iblk)%p(j,1),4)) then
                  nwalli=nwalli+1
! assemble local coord list
                  do node = 1, nenbl
                     xwi(nwalli,node,1:3)=x(mienb(iblk)%p(j,node),:)
                  enddo
! put the centroid coordinates in the last slot
                  do node = 1, nenbl
                     xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:) &
                          +xwi(nwalli,node,:)
                  enddo
                  xwi(nwalli,nenb+1,:)=xwi(nwalli,nenb+1,:)/nenbl
!
               endif
            enddo               ! loop over belts in this blk
!
         enddo                  ! loop over boundary elt blocks
!
         if (nwalli.eq.0) xwi=1.0e32 ! fix for mpi's lack of imagination
         if (nwalli.eq.0) nwalli=1 !  patch for mpi's lack of imagination
!
!  Pool "number of welts" info from all processors
!
!MPI_ALLGATHER(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,comm) 
![ IN sendbuf] starting address of send buffer (choice) 
![ IN sendcount] number of elements in send buffer (integer) 
![ IN sendtype] data type of send buffer elements (handle) 
![ OUT recvbuf] address of receive buffer (choice) 
![ IN recvcount] number of elements received from any process (integer) 
![ IN recvtype] data type of receive buffer elements (handle) 
![ IN comm] communicator (handle)
         if (numpe.gt.1) then   ! multiple processors
! write the number of wall elts on the jth processor to slot j of nwall
            call MPI_ALLGATHER(nwalli,1,MPI_INTEGER,nwall,1, &
                MPI_INTEGER,INEWCOMM,ierr)
! count up the total number of wall elts among all processes
            nwallt=0
            do j=1,numpe
               nwallt=nwallt+nwall(j)
            enddo
         else                   ! single processor
! the local information is the global information for single-processor
            nwall=nwalli
            nwallt=nwalli
         endif                  ! if-else for multiple processors
!
!  Make all-processor wallnode-coord collage
!
         allocate (xw(nwallt,nenb+1,nsd))
         if (numpe.gt.1) then   ! multiple processors
! we will gather coordinates from local on-proc sets to a global set
! we will stack each processor's coordinate list atop that of the
! previous processor.  If the coordinate list for processor i is
! called xwi, then our global coordinate list xw will look like this:
! ---------------------------------------------------------------
! | xw1            | xw2                | xw3        |   ...    |
! ---------------------------------------------------------------
!  <---nwall(1)---> <-----nwall(2)-----> <-nwall(3)->
!  <------------------------nwallt-----------------------...---->
! To accomplish this with MPI, we use MPI_ALLGATHERV, summarized as:
!MPI_ALLGATHERV(sendbuf,sendcount,sendtype,recvbuf,recvcount,disp,recvtype,comm) 
![ IN sendbuf] starting address of send buffer (choice) 
![ IN sendcount] number of elements in send buffer (integer) 
![ IN sendtype] data type of send buffer elements (handle) 
![ OUT recvbuf] address of receive buffer (choice) 
![ IN recvcount] number of elements received from any process (integer) 
![ IN disp] displacement array
![ IN recvtype] data type of receive buffer elements (handle) 
![ IN comm] communicator (handle)
! The displacement array disp is an array of integers in which the jth
! entry indicates which slot of xw marks beginning of xwj
! So, first we will build this displacement array
            idisp(:)=0 ! starting with zero, since MPI likes C-numbering
            do j=2,numpe
               idisp(j)=idisp(j-1)+nwall(j-1) ! see diagram above
            enddo
! Now, we gather the data one slice at a time (1:nwalli)
            do j=1,nenb+1
               do k=1,nsd
                  call MPI_ALLGATHERV(xwi(:,j,k),nwalli, &
                       MPI_DOUBLE_PRECISION,xw(:,j,k),nwall,idisp, &
                       MPI_DOUBLE_PRECISION,INEWCOMM,ierr)
               enddo
            enddo
         else                   ! single-processor
! global data is local data in single processor case
            xw=xwi
         endif

!
!  For each point, loop over wall nodes and calculate distance; 
!  save the distance in this node's position of d2wall if it's 
!  shorter than currently stored distance
!
         d2wall=1.0e32
         do i=1,numnp
            do j=1, nwallt
               do k=1,nenb+1
                  distance =  ( x(i,1) - xw(j,k,1) )**2 &
                       +( x(i,2) - xw(j,k,2) )**2 &
                       +( x(i,3) - xw(j,k,3) )**2
                  if ( d2wall(i).gt.distance ) d2wall(i) = distance
               enddo
            enddo
         enddo
         d2wall=sqrt(d2wall)
!
         deallocate(xwi)
         deallocate(xw)
!
!.... write d2wall to a file so we don't have to do this again
!
         write (fmt1,"('(''d2wall.'',i',i1,',1x)')") 1
         write (fname1,fmt1) 0
         fname1 = trim(fname1) // cname(myrank+1)
         open (unit=72, file=fname1, status='unknown', &
                                          form='unformatted', err=996)

         write (72) d2wall
         close (72)
!         write(*,*) "make sure to: echo 0 > distcalc.dat"
!         call MPI_BARRIER (INEWCOMM,ierr)
!         call error ('distcalc','complete',0)
      endif
      if (idistcalc.eq.0) then ! d2wall is already done, don't calculate it
        write (fmt1,"('(''d2wall.'',i',i1,',1x)')") 1
        write (fname1,fmt1) 0
        fname1 = trim(fname1) // cname(myrank+1)
        write (*,*) 'Reading dist file : ', fname1

        open (unit=72, file=fname1, status='old', &
              form='unformatted', err=995)

        read (72) d2wall
        close (72)
      endif

      return
995     call error ('d2wall  ','opening ', 72)
996     call error ('d2wall  ','opening ', 72)

      end
