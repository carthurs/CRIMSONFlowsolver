!
!     boundary data module
!     not implicit, requires data from simV
!
      module boundarymodule      
      
      use pointer_data, only : r3d

      integer                :: maxSurfID       ! global max surface IDs
      integer                :: numSurfID       ! local number of surfaces on this processor
      integer, allocatable   :: surfID(:)       ! local surface IDs on this processor           
      integer, allocatable   :: lstSurfID(:)    ! global surface IDs

      real*8, pointer        :: surfFlow(:) => null() !
      real*8, pointer        :: surfPres(:) => null()


      real*8, allocatable    :: area(:)         ! surface area of each surface
      real*8, allocatable    :: centroid(:,:)   ! centroid of each surface
      integer, allocatable   :: edsurf(:)       ! edge nodes
      integer, allocatable   :: edpair(:,:)     ! local edge pairs on this processor
      type(r3d), allocatable :: edgePairs(:)    ! global coordinates of each edge pair
      real*8, allocatable    :: nodeRadius(:)    ! local radii the nodes on each surface
      real*8, allocatable    :: nodeProfile(:)

      real*8, allocatable    :: qprf(:)


      contains

      !!
      !! master subroutine called from proces.f
      !!

      subroutine setupBoundaryModule(x)

      use pvsQbi, only : ndsurf 

      ! include "common.h"
      use phcommonvars

      real*8  :: x(numnp,nsd) ! xyz coordinates


      ! find max surf IDs
      call findSurfIDs(ndsurf, nshg, MAXSURF)

      ! find edge nodes 
      call findEdgeNodes(ndsurf, nshg, nelblb, lcblkb, MAXBLK)

      ! calculate area and centroid of each surface
      call calculateArea(ndsurf, nshg, x, numnp, nsd)

!       ! find edge pairs
!       call findEdgePairs(ndsurf, nshg, nelblb, lcblkb, MAXBLK)

!       ! share edge pairs
!       call shareEdgePairs(ndsurf, nshg, x, numnp, nsd)

!       ! calculate local radius 
!       call calculateRadii(ndsurf, nshg, x, numnp, nsd)

!       ! calculate velocity profile 
!       call calculateProfile(ndsurf, nshg, nprofile)

! c      open(unit=834,file='radius.dat',status='replace')
! c      do i = 1, nshg
! c        if (nodeRadius(i) .gt. real(0.0,8)) then
! c          if (ndsurf(i) .gt. int(1)) then
! c          write(834,'(i5,5(f12.6))') ndsurf(i), 
! c     &                               x(i,1:3), 
! c     &                               nodeRadius(i),
! c     &                               nodeProfile(i)  
! c          else          
! c          write(834,'(i5,4(f12.6))') edsurf(i), 
! c     &                               x(i,1:3), 
! c     &                               nodeRadius(i),
! c     &                               nodeProfile(i)          
! c          end if          
! c
! c        end if
! c      end do 
! c      close(834)

      ! write details to file 
      call outputBoundaryStatistics(ndsurf, nshg, x, numnp, nsd, iprofile, nprofile)

      return
      end subroutine setupboundarymodule

      !!
      !! subroutine to find the maximum surface id using ndsurf 
      !!
      
      subroutine findSurfIDs(ndsurf, nshg, MAXSURF) 

      use mpi, only : MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD

      integer :: nshg
      integer :: ndsurf(nshg) 
      integer :: MAXSURF

      integer :: lmaxval, curSurfID
      integer :: mpierr, mpirank
      integer :: i, j 
      integer :: locHasID, glbHasID

      ! find the local maximum surface ID 
      lmaxval = MAXVAL(ndsurf)
      
      ! find maximum surface ID over all processors
      call MPI_ALLREDUCE(lmaxval, maxSurfID, int(1), MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpierr) 

      ! count number of surfaces on this processor
      numSurfID = 0
      do i = 2, maxSurfID 
        do j = 1, nshg
          if (ndsurf(j) .eq. i) then
            numSurfID = numSurfID + 1
            exit
          end if 
        end do
      end do

      ! if number of surfaces on this processor > 0, then save the surface ids
      if (numSurfID .gt. int(0)) then
        
        allocate(surfID(numSurfID))
        curSurfID = 1
        
        ! maxSurfID loop
        do i = 2, maxSurfID
          do j = 1, nshg
            if (ndsurf(j) .eq. i) then
              surfID(curSurfID) = ndsurf(j)
              curSurfID = curSurfID + 1
              exit 
            end if   
          end do 
        end do
        
      end if

      ! set up global surface list and global surface count
      allocate(lstSurfID(0:MAXSURF))
      lstSurfID(:) = int(0)

      do i = 1, maxSurfID      
        lstSurfID(i) = i
      end do

      ! allocate pointer arrays
      allocate(surfFlow(maxSurfID))
      surfFlow(:) = real(0.0,8)

      allocate(surfPres(maxSurfID))
      surfPres(:) = real(0.0,8)


      return
      end subroutine findSurfIDs
  
      !!
      !!
      !!

      subroutine findEdgeNodes(ndsurf, nshg, nelblb, lcblkb, MAXBLK)
      
      use pointer_data, only : miBCB    ! pointer to iBCB array
      use pointer_data, only : mienb    ! pointer to ienb array
      integer :: ndsurf(nshg)           ! global node number array 
      integer :: nelblb                 ! number of boundary element block
      integer :: lcblkb(10,MAXBLK+1)    ! blocking data for the boundary elements

      integer :: iblk, iel, npro, iedge(3), ipair(2), e, n
      integer :: elemcount, edgecount
      integer, allocatable :: edgetemp(:,:)
     
      ! allocate an array the size of the total number of nodes 
      allocate(edsurf(nshg))
      edsurf(:) = int(0)

      ! loop over blocks
      do iblk = 1, nelblb

        ! calculate number of boundary elements in this block 
        iel = lcblkb(1,iblk)        
        npro = lcblkb(1,iblk+1) - iel

        ! loop over elements
        do e = 1, npro
               
          ! check if element is tagged with surface number > 1
          if(miBCB(iblk)%p(e,2) .gt. int(1)) then

              ! check if any of it's nodes are edge nodes
              iedge(1:3) = isEdge(mienb(iblk)%p(e,1:3), ndsurf, nshg)

              ! if they are, set edsurf(nodeid) = surf id
              do n = 1, 3
                
                ! set element surface id to this edge id
                if (iedge(n) .gt. int(0)) then                
                  edsurf(mienb(iblk)%p(e,n)) = miBCB(iblk)%p(e,2)
                end if

              end do 

          end if                  

        ! end of element loop
        end do

      ! end of block loop
      end do     

      end subroutine findEdgeNodes


      !!
      !! function to check if the boundary nodes passed in are edge nodes
      !! this assumes the wall is always set to 1 - KDL NAN - 21/08/14
      !! 

      function isEdge(nodeid, ndsurf, nshg)
      
      integer :: isedge(3), nodeid(3), ndsurf(nshg) 

      ! set all nodes to false 
      isedge(:) = int(0)   
      
      ! loop over the 3 node, if the ndsurf is 1 then set true
      do i = 1, 3
        if (ndsurf(nodeid(i)) .eq. int(1)) then
          isedge(i) = int(1)
        end if
      end do
   
      return
      end function isEdge    

      !!
      !! calculate area
      !! 

      subroutine calculateArea(ndsurf, nshg, x, numnp, nsd)

      use pvsQbi, only : NASC   ! enables access to NASC arrays
      use mpi, only : MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD


      integer :: nshg, numnp, nsd                 ! 
      integer :: ndsurf(nshg)                     ! global node number array 
      real*8  :: x(numnp,nsd)                     ! xyz coordinates
      real*8  :: area_integral(maxSurfID), &      !
                 surf_integral(maxSurfID,3), &    !
                 surf_integral_sum(maxSurfID,3)   !
      integer :: i, indx, mpierr

      area_integral(:) = 0.0d0
      surf_integral(:,:) = 0.0d0

      ! loop over all the nodes on this processor
      do i = 1, nshg

         ! find nodes within the boundary using ndsurf > 1
         if (ndsurf(i) .gt. int(1)) then
            
            ! surface number of this node
            indx = ndsurf(i)
            
            ! area and centroid integral
            area_integral(indx) = area_integral(indx) + NASC(i)*real(1.0d0,8)    
            surf_integral(indx,1:3) = surf_integral(indx,1:3) + NASC(i)*x(i,1:3)
         
         ! find nodes on the boundary using edsurf > 1
         else if (edsurf(i) .gt. int(1)) then
            
            ! surface number of this node
            indx = edsurf(i)
            
            ! area and centroid integral
            area_integral(indx) = area_integral(indx) + NASC(i)*real(1.0d0,8)    
            surf_integral(indx,1:3) = surf_integral(indx,1:3) + NASC(i)*x(i,1:3)

         end if 
      end do

      ! here calculate area

      allocate(area(maxSurfID))
      area(:) = 0.0d0

      call MPI_ALLREDUCE(area_integral, area, maxSurfID, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) 

      ! here calculate centriod 
      surf_integral_sum(:,:) = 0.0d0

      call MPI_ALLREDUCE(surf_integral(:,1), surf_integral_sum(:,1), maxSurfID, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) 

      call MPI_ALLREDUCE(surf_integral(:,2), surf_integral_sum(:,2), maxSurfID, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) 

      call MPI_ALLREDUCE(surf_integral(:,3), surf_integral_sum(:,3), maxSurfID, & 
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)       


      allocate(centroid(maxSurfID,3))
      centroid(:,:) = 0.0d0

      centroid(2:maxSurfID,1) = surf_integral_sum(2:maxSurfID,1) &
                              / area(2:maxSurfID)
      
      centroid(2:maxSurfID,2) = surf_integral_sum(2:maxSurfID,2) &
                              / area(2:maxSurfID)
      
      centroid(2:maxSurfID,3) = surf_integral_sum(2:maxSurfID,3) &
                              / area(2:maxSurfID)
      
      return
      end subroutine calculateArea

      !!
      !! find boundary pairs
      !!

      subroutine findEdgePairs(ndsurf, nshg, nelblb, lcblkb, MAXBLK)

      use pointer_data, only : miBCB, & ! pointer to iBCB array
                               mienb    ! pointer to ienb array
      integer :: ndsurf(nshg)           ! global node number array 
      integer :: nelblb                 ! number of boundary element block
      integer :: lcblkb(10,MAXBLK+1)    ! blocking data for the boundary elements

      integer :: iblk, iel, npro, ipair(2), e, n
      integer :: elemcount, edgecount
      integer, allocatable :: edgetemp(:,:)
     
      ! allocate temporary edge array with the number of boundary elements
      elemcount = int(0)
      do iblk = 1, nelblb
        elemcount = elemcount + lcblkb(1,iblk+1) - lcblkb(1,iblk)  
      end do
      allocate(edgetemp(elemcount,3))
      edgecount = int(0)

      ! loop over blocks
      do iblk = 1, nelblb

        ! calculate number of boundary elements in this block 
        iel = lcblkb(1,iblk)        
        npro = lcblkb(1,iblk+1) - iel

        ! loop over elements
        do e = 1, npro
               
          ! check if element is tagged with surface number > 1
          if(miBCB(iblk)%p(e,2) .gt. int(1)) then

              ipair(:) = int(0)
              ! check if this element's nodes are an edge pair using ndsurf
              ipair(1:2) = isEdgePair(mienb(iblk)%p(e,1:3), ndsurf, nshg)
              !write(*,*) ipair
              
              ! store surf ID and node IDs in edgetemp
              if (ipair(1) .ne. int(0)) then
                if (ipair(2) .ne. int(0)) then
                  edgecount = edgecount + 1
                  edgetemp(edgecount,1) = miBCB(iblk)%p(e,2)   
                  edgetemp(edgecount,2:3) = ipair(1:2)     
                end if
              end if 

          end if                  

        ! end of element loop
        end do

      ! end of block loop
      end do     

      ! allocate edpair array and store temp array
      allocate(edpair(edgecount,3))
      do e = 1, edgecount
        edpair(e,1:3) = edgetemp(e,1:3) 
        ! write(*,*) edpair(e,:) 
      end do

      ! deallocate edgetemp array
      deallocate(edgetemp)
      
      return

      end subroutine findEdgePairs


      !!
      !! function to check if the boundary nodes passed in are an edge pair
      !! 

      function isEdgePair(nodeid, ndsurf, nshg)
      
      integer :: isEdgePair(2), nodeid(3), ndsurf(nshg), count, n

      ! isEdgePair is an array containing the node IDs for this pair 
      isEdgePair(:) = int(0)   
      
      ! loop over the 3 node ids and count edge nodes
      count = int(0)    
      do i = 1, 3      
        if (ndsurf(nodeid(i)) .eq. int(1)) then
          count = count + int(1)
        end if
      end do


      ! if the number of edge nodes on this element = 2, store node IDs
      if (count .eq. int(2)) then
        n = 1
        do i = 1, 3      
          if (ndsurf(nodeid(i)) .eq. int(1)) then
            isEdgePair(n) = nodeid(i)

! c            write(*,*) 'n ',n,' nID ', isEdgePair(n)
            n = n + 1            
          end if
        end do
      end if 
   
      return

      end function isEdgePair   


      !!
      !! share edge pairs 
      !!

      subroutine shareEdgePairs(ndsurf, nshg, x, numnp, nsd)

      use mpi, only : MPI_COMM_RANK, MPI_COMM_WORLD, MPI_COMM_SIZE, &
                      MPI_INTEGER, MPI_DOUBLE_PRECISION 

      integer :: nshg
      integer :: numnp
      integer :: nsd
      integer :: ndsurf(nshg)      
      real*8  :: x(numnp, nsd)

      integer :: edgesPerSurf(maxSurfID) ! number of pairs
      
      integer :: nProcs 
      integer :: i, j, k
      integer :: surf, pair, count
      integer, allocatable :: procRanks(:)
      integer :: recvindx, sendindx, status
      integer :: sendsize, offset
      integer, allocatable :: edgesPerProc(:) !               
      real*8, allocatable :: recvbuff(:), sendbuff(:)
      
      integer :: idim, ipair, itag
      integer :: mpirank, mpierr

      ! get processor rank and size
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)        
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpisize, mpierr)

      ! count the number of pairs for each surface on each processor
      edgesPerSurf(:) = int(0)
      if (SIZE(edpair,1) .gt. int(0)) then
        do i = 2, maxSurfID
          do j = 1, SIZE(edpair,1)
            if (i .eq. edpair(j,1)) then
              edgesPerSurf(i) = edgesPerSurf(i) + 1
            end if 
          end do
        end do
      end if

      ! allocate r3d pointer for the edge pair coordinates
      allocate(edgePairs(maxSurfID))

      ! loop over all the surfaces
      do surf = 2, maxSurfID

        ! allocate parallel count
        if (allocated(edgesPerProc)) then
          deallocate(edgesPerProc)
        end if 
        allocate(edgesPerProc(mpisize))
        edgesPerProc(:) = int(0)
        
        ! all gather to share the number of pairs on each rank 
        call MPI_ALLGATHER(edgesPerSurf(surf), &
                           int(1), &
                           MPI_INTEGER, &
                           edgesPerProc, &
                           int(1), &
                           MPI_INTEGER, &
                           MPI_COMM_WORLD, &
                           mpierr)


        ! count how many processors have this surface
        nProcs = 0
        do i = 1, mpisize
          if (edgesPerProc(i) .gt. int(0)) then
            nProcs = nProcs + 1
          end if 
        end do 


        ! allocate number of edge pairs if on this processor
        if (edgesPerSurf(surf) .gt. int(0)) then
          allocate(edgePairs(surf)%p(SUM(edgesPerProc), 2, 3))
        end if 


        ! PARALLEL/SERIAL  

        ! here the edges are shared with the relevant processors
        if (nProcs .gt. int(1)) then

          ! find out which processor ranks share this surface
          if (allocated(procRanks)) then
            deallocate(procRanks)
          end if 
          allocate(procRanks(nProcs))
          procRanks(:) = int(0)

          count = 0
          do i = 1, mpisize
            if (edgesPerProc(i) .gt. int(0)) then
              count = count + 1
              procRanks(count) = i - 1 ! rank 0 -> mpisize
            end if 
          end do 

          ! loop through the sending and receiving processor pairs
          do sendindx = 1, nProcs
            do recvindx = 1, nProcs

              ! avoid sending data to oneself
              if (sendindx .ne. recvindx) then

                ! allocate send and recv buffers before MPI operations
                sendsize = edgesPerProc(procRanks(sendindx) + 1)                                              

                if (allocated(sendbuff))then
                  deallocate(sendbuff)
                end if 
                allocate(sendbuff(sendsize))
                if (allocated(recvbuff)) then
                  deallocate(recvbuff)
                end if                 
                allocate(recvbuff(sendsize))

                ! check if sending 
                if (mpirank .eq. procRanks(sendindx)) then

                  ! loop over each pair and all 3 dimensions
                  do ipair = 1, 2
                    do idim = 1, 3

                      ! copy coordinates to send buffer
                      j = 0          
                      do i = 1, SIZE(edpair,1)
                        if (edpair(i,1) .eq. surf) then              
                          j = j + 1
                          sendbuff(j) = x(edpair(i,ipair+1),idim)
                        end if 
                      end do 
                            
                      ! set tag for MPI, itag = 1 -> 6                            
                      itag = 3*(ipair-1)+idim 

                      ! MPI send
                      call MPI_SEND(sendbuff, &             ! send buffer
                                    sendsize, &             ! size of buffer
                                    MPI_DOUBLE_PRECISION, & ! MPI data type  
                                    procRanks(recvindx), &  ! receiving rank 
                                    itag, &                 ! tag/key number
                                    MPI_COMM_WORLD, &       ! communicator
                                    mpierr)                 ! error handle

                      ! save sent data  
                      if (sendindx .eq. int(1)) then
                        do i = 1, sendsize
                          edgePairs(surf)%p(i, ipair, idim) = sendbuff(i)
                        end do 
                      else                         
                        offset = int(0)
                        do i = 1, sendindx-1
                          offset = offset + edgesPerProc(procRanks(i)+1)
                        end do                         
                        do i = 1, sendsize
                          edgePairs(surf)%p(offset+i, ipair, idim) = sendbuff(i)
                        end do 
                      end if 

                    ! end of ipair and idim loops  
                    end do
                  end do 

                ! check if recieving 
                else if (mpirank .eq. procRanks(recvindx)) then

                  ! loop over each pair and all 3 dimensions
                  do ipair = 1, 2
                    do idim = 1, 3

                      ! zero recieving buffer
                      recvbuff(:) = real(0.0,8)

                      ! set tag for MPI, itag = 1 -> 6                            
                      itag = 3*(ipair-1)+idim 

                      ! MPI receive
                      call MPI_RECV(recvbuff, &             ! recv buffer
                                    sendsize, &             ! size of buffer
                                    MPI_DOUBLE_PRECISION, & ! MPI datatype
                                    procRanks(sendindx), &  ! sending rank
                                    itag, &                 ! tag/key number
                                    MPI_COMM_WORLD, &       ! communicator
                                    status, &               ! MPI status
                                    mpierr)                 ! error handel     

                      ! save received data  
                      if (sendindx .eq. int(1)) then
                        do i = 1, sendsize
                          edgePairs(surf)%p(i, ipair, idim) = recvbuff(i)
                        end do 
                      else 
                        offset = int(0)
                        do i = 1, sendindx-1
                          offset = offset + edgesPerProc(procRanks(i)+1)
                        end do                         
                        do i = 1, sendsize
                          edgePairs(surf)%p(offset+i, ipair, idim) = recvbuff(i)
                        end do 
                      end if   

                    ! end of ipair and idim loops
                    end do  
                  end do

                ! end of send/recv if   
                end if                     
              end if 

            ! end sendindx and recvindx loops  
            end do 
          end do 

        ! else if edges are only on one processor
        else 
          
          ! check if on this processor and store data
          if (edgesPerSurf(surf) .gt. int(0)) then
            count = 0          
            do pair = 1, SIZE(edpair,1) 
              if (edpair(pair,1) .eq. surf) then              
                count = count + 1
                edgePairs(surf)%p(count, 1, 1:3) = x(edpair(pair,2),1:3)
                edgePairs(surf)%p(count, 2, 1:3) = x(edpair(pair,3),1:3) 
                end if 
            end do 
          end if 

        ! end of 
        end if 

! c        ! store coordinates if on this processor
! c        if (edgesPerSurf(surf) .gt. int(0)) then
! c          do j = 1, SUM(edgesPerProc)
! c            write(*,'(i5,i5,6(f10.5))') mpirank, j, 
! c     &                                  edgePairs(surf)%p(j, 1, 1:3), 
! c     &                                  edgePairs(surf)%p(j, 2, 1:3)     
! c          end do
! c        end if 
    
      ! end surface loop  
      end do

      ! deallocate edpair array
      deallocate(edpair)

      return

      end subroutine shareEdgePairs

      !!
      !!
      !!


      subroutine calculateRadii(ndsurf, nshg, x, numnp, nsd)

      integer :: nshg
      integer :: numnp
      integer :: nsd
      integer :: ndsurf(nshg)      
      real*8  :: x(numnp, nsd)

      integer :: i, j
      real*8  :: dist2edge(4), minDist, minXYZ(3)
      real*8 :: edge2centroid(3), length, node2centriod(3)
      real*8 :: localRadius, nodeDistance
      integer :: rnum = 764

      ! allocate nodeRadius array
      allocate(nodeRadius(nshg))
      nodeRadius(:) = real(0.0d0,8)

! cc      open(unit=rnum,file='radius.dat',status='replace')

      ! loop over all the nodes on this processor
      do i = 1, nshg

        ! check if each node is surface node
        if (ndsurf(i) .gt. int(1)) then
              
          ! find the distance from this node to the first edge pair
          dist2edge = distanceToEdge(x(i,1:3), &
                                     edgePairs(ndsurf(i))%p(1,1,1:3), &
                                     edgePairs(ndsurf(i))%p(1,2,1:3)) 

          minDist = dist2edge(1)
          minXYZ(1:3) = dist2edge(2:4)

          ! loop over remaining edge pairs and find the minimum distance
          do j = 2, SIZE(edgePairs(ndsurf(i))%p, 1)

            dist2edge = distanceToEdge(x(i,1:3), &
                                       edgePairs(ndsurf(i))%p(j,1,1:3), &
                                       edgePairs(ndsurf(i))%p(j,2,1:3)) 
            

            if (dist2edge(1) .lt. minDist) then
              minDist = dist2edge(1)
              minXYZ(1:3) = dist2edge(2:4)
            end if 
            
          end do 


          ! calculate the local radius 
          edge2centroid = minXYZ(1:3) - centroid(ndsurf(i),1:3)
          localRadius = real(0d0,8)
          localRadius = DOT_PRODUCT(edge2centroid,edge2centroid)
          localRadius = SQRT(localRadius)

          ! calculate the node location relative to the centroid 
          node2centriod = x(i,1:3) - centroid(ndsurf(i),1:3)
          nodeDistance = real(0d0,8)
          nodeDistance = DOT_PRODUCT(node2centriod,node2centriod)
          nodeDistance = SQRT(nodeDistance)

          ! store the normalised radial distance
          nodeRadius(i) =  nodeDistance/localRadius

! cc          write(rnum,'(i4,4(e16.6))') ndsurf(i), x(i,1:3),
! cc     &                             real(1,8)- nodeDistance/localRadius

        else if (edsurf(i) .gt. int(1)) then

          ! if edge node set normalised radius to 1
          nodeRadius(i) = real(1.0,8)

        end if
      end do

! cc      close(rnum)




      end subroutine calculateRadii
      
      !!
      !! calculate the distance between point p and the segment made by x1 and x2.
      !! 

      function  distanceToEdge(p, x1, x2) 

      real*8 :: p(3), x1(3), x2(3), t
      real*8 :: distance(3), point(3)
      real*8 :: distanceToEdge(4) 
      integer :: i

      distance = x2 - x1
! c      write(*,'(6(f12.3))') x2, x1

      ! calculate the t that minimizes the distance.
      t = real(0.0d0,8)
      do i = 1, 3
        t = t + ( p(i) - x1(i) ) * distance(i)
      end do 
      t = t / ( distance(1)*distance(1) + &
                distance(2)*distance(2) + &
                distance(3)*distance(3) )  
      
      ! see if this represents one of the segment's
      ! end points or a point in the middle.
      if (t .lt. real(0,8)) then        
        distance = p - x1
        distanceToEdge(2:4) = x1(1:3)
      else if (t .gt. real(1.0,8)) then        
        distance = p - x2
        distanceToEdge(2:4)  = x2(1:3)
      else        
        distanceToEdge(2:4)  = x1 + t * distance
        distance = p - distanceToEdge(2:4)
      end if 

      distanceToEdge(1) = distance(1)*distance(1) &
                        + distance(2)*distance(2) &
                        + distance(3)*distance(3)  
      distanceToEdge(1) = sqrt(distanceToEdge(1))

      end function distanceToEdge

      !!
      subroutine calculateProfile(ndsurf, nshg ,order)

      use pvsQbi, only : NASC   ! enables access to NASC arrays
      use mpi, only : MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD

      integer :: nshg, order
      integer :: ndsurf(nshg)
      integer :: i, indx
      real*8  :: coeff 
      real*8  :: nodeIntegral(maxSurfID), profileSum(maxSurfID)

      coeff = real(order,8) + real(2,8)
      coeff = coeff/real(order,8)

      ! allocate normalised profile array
      allocate(nodeProfile(nshg))
      nodeProfile(:) = real(0.0,8)
      
      ! also allocate dimensionalised array, although not used here
      allocate(qprf(nshg))
      qprf(:) = real(0.0,8)


      do i = 1, nshg
        if (ndsurf(i) .gt. int(1)) then
           nodeProfile(i) = real(1.0,8) - (nodeRadius(i)**order)
           nodeProfile(i) = coeff*nodeProfile(i) 
        else if (edsurf(i) .gt. int(1)) then
           nodeProfile(i) = real(1.0,8) - (nodeRadius(i)**order)
           nodeProfile(i) = coeff*nodeProfile(i)         
        end if
      end do 

      nodeIntegral(:) = 0.0d0
      do i = 1, nshg
         if (ndsurf(i) .gt. int(1)) then
            indx = ndsurf(i)
            nodeIntegral(indx) = nodeIntegral(indx) + NASC(i)*nodeProfile(i)
         else if (edsurf(i) .gt. int(1)) then
            indx = edsurf(i)
            nodeIntegral(indx) = nodeIntegral(indx) + NASC(i)*nodeProfile(i)
         end if 
      end do

      ! integrate node profile over its surface
      profileSum(:) = 0.0d0
      
      call MPI_ALLREDUCE(nodeIntegral, profileSum, maxSurfID, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) 

      profileSum(2:maxSurfID) = profileSum(2:maxSurfID) / area(2:maxSurfID)

! cc      write(*,*) 'PROFILESUM: ',profileSum(2:maxSurfID)

      ! normalise velocity profile to 1
      do i = 1, nshg
         if (ndsurf(i) .gt. int(1)) then
            indx = ndsurf(i)
            nodeProfile(i) = nodeProfile(i) / profileSum(indx)
         end if 
      end do

      end subroutine calculateProfile

      !!
      !! output boundary surface data in parallel
      !!

      subroutine outputBoundaryStatistics(ndsurf, nshg, x, numnp, nsd, iprfl, nprfl)

      use mpi, only : MPI_COMM_WORLD
      integer :: nshg
      integer :: numnp
      integer :: nsd
      integer :: ndsurf(nshg)      
      real*8  :: x(numnp, nsd)
      integer :: iprfl
      integer :: nprfl

      integer :: mpirank, mpierr, i, j, k
      integer :: fnum = 947
      character(len=100) :: chararray1, chararray2

      if (numSurfID .gt. int(0)) then

        call MPI_COMM_RANK(MPI_COMM_WORLD, mpirank, mpierr)    

        ! at Nan's request, here we follow the Fortran indexing
        write (chararray1,'(i10)') mpirank + int(1)
        write (chararray2,'(a,a)') 'boundary.data.', trim(adjustl(chararray1))

        open(unit=fnum, file=chararray2, status='replace')
        write(fnum,'(a)')   '---------------------------'
        write(fnum,'(a)')   'BOUNDARY DATA MODULE OUTPUT'
        write(fnum,'(a,/)') '---------------------------'

        write(fnum,'(a)') 'MPI processor rank'
        write(fnum,'(i12,3x,/)') mpirank + int(1)

! cc        write(fnum,'(a)') 'Number of surfaces in this model'
! cc        write(fnum,'(i12,3x,/)') surfCount (:)     
! cc
! cc        write(fnum,'(a)') 'Surface IDs in this model'
! cc        do i = 1, surfCount
! cc          write(fnum,'(i12,3x,/)') surfList(i)
! cc        end do 

        write(fnum,'(a)') 'Number of surfaces on this rank'
        write(fnum,'(i12,3x,/)') numSurfID

        write(fnum,'(a)') 'Surface IDs on this rank'
        do i = 1, numSurfID
          write(fnum,'(i12,3x)') surfID(i)
        end do

        if (iprfl .gt. int(0)) then
          write(fnum,'(/,a)') 'Velocity profile order'
          write(fnum,'(i12,3x)') nprfl
        end if 

        do i = 1, numSurfID
          write(chararray1,'(i)') surfID(i)
          
          write(fnum,'(/,a)') '---------------------------'
          write(fnum,'(a,a)') 'Surface number ',trim(adjustl(chararray1))
          write(fnum,'(a,/)') '---------------------------'

          write(fnum,'(a)') 'Area'
          write(fnum,'(e12.5,3x,/)') area(surfID(i))

          write(fnum,'(a)') 'Centroid'
          write(fnum,'(3(e12.5,3x),/)') centroid(surfID(i),:)

          ! count number of surface nodes 
          k = 0
          do j = 1, nshg
            if (ndsurf(j) .eq. surfID(i) .or. edsurf(j) .eq. surfID(i)) then
              k = k + 1
            end if
          end do 
          write(fnum,'(a,/,i12,3x,/)') 'Boundary nodes on this rank', k 
    
          ! ! write out the boundary nodes on this rank
          ! write(fnum,'(a11,34x,a12,3x,a12,3x)') 'Coordinates','Radius ratio', 'Flow profile'
          ! write(fnum,'(5(a1,14x))') 'X','Y','Z','R','Q'

          ! do j = 1, nshg
          !   if (ndsurf(j) .eq. surfID(i) .or. edsurf(j) .eq. surfID(i)) then    
          !     write(fnum,'(5(e12.5,3x))') x(j,1:3), nodeRadius(j), nodeProfile(j)
          !   end if
          ! end do  

! c          write(fnum,'(/,a,/,i12,3x,/)') 'Boundary edges [global]', 
! c     &                                    SIZE(edgePairs(surfID(i))%p,1) 

! c          write(fnum,'(a)') 'Coordinates'
! c          write(fnum,'(6(a2,13x))') 'X1','Y1','Z1','X2','Y2','Z2'  
! c          do j = 1, SIZE(edgePairs(surfID(i))%p,1)
! c            write(fnum,'(6(e12.5,3x))') edgePairs(surfID(i))%p(j,1,1:3),
! c     &                                  edgePairs(surfID(i))%p(j,2,1:3)      
! c          end do  

        end do 

        ! close file
        close(fnum)

      end if 

      return
      end subroutine outputBoundaryStatistics


      !
      subroutine updateFlow(y, nshg, MAXSURF)

      real*8  :: y(nshg,1:3)
      integer :: nshg
      integer :: MAXSURF

      real*8  :: flowQ(0:MAXSURF)
      integer :: i

      call GetFlowQ(flowQ, y, lstSurfID, maxSurfID)

      do i = 2, maxSurfID
        surfFlow(i) = flowQ(i)
! cc        write(*,*) 'SURFID: ',i, ' FLOW ', surfFlow(i)        
      end do 

      call setFlowProfile(nshg)

      return
      end subroutine


      !
      subroutine updatePressure()
      return
      end subroutine

      !
      ! set flow profile for surface
      !

      subroutine setFlowProfile(nshg)

      use pvsQbi, only : ndsurf 

      integer :: surf

      ! loop over the nodes on this processor      
      do i = 1, nshg

        ! check if ndsurf or edsurf is > 1
        if (ndsurf(i) .gt. int(1)) then

          surf = ndsurf(i)          

          ! check if flow is entering domain          
          if (surfFlow(surf) .lt. real(0.0,8)) then
            qprf(i) = surfFlow(surf) * nodeProfile(i)
          else
            qprf(i) = real(0.0,8)
          end if 

        else if (edsurf(i) .gt. int(1)) then 
          
          surf = edsurf(i)          

          ! check if flow is entering domain          
          if (surfFlow(surf) .lt. real(0.0,8)) then
            qprf(i) = surfFlow(surf) * nodeProfile(i)
          else 
            qprf(i) = real(0.0,8)
          end if 
          
        else 

          qprf(i) = real(0.0,8)

        end if 

      end do 
  
      return
      end subroutine setFlowProfile

!
! *** New function to integrate scalars over a surface 
!     KDL and NAN 21/08/14
!
      subroutine integrScalar(scalInt, scal, srfIdList, numSrfs)

      use pvsQbi ! ndsurf, NASC
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      
      integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled, i, k
      real*8    scal(nshg), scalInt(0:MAXSURF), scalIntProc(0:MAXSURF)
      
      scalIntProc = zero

      do i = 1,nshg

        if (numSrfs.gt.zero) then
          
          do k = 1,numSrfs
          
            irankCoupled = 0
          
            if (srfIdList(k) .eq. ndsurf(i)) then

              irankCoupled=k
              scalIntProc(irankCoupled) = scalIntProc(irankCoupled) + NASC(i)*scal(i)

            else if (srfIdList(k) .eq. edsurf(i)) then

              irankCoupled=k
              scalIntProc(irankCoupled) = scalIntProc(irankCoupled) + NASC(i)*scal(i)


            end if      
      
          end do       
      
        end if
      
      end do
!      
!     at this point, each scalint has its "nodes" contributions to the scalar
!     accumulated into scalIntProc. Note, because NASC is on processor this
!     will NOT be the scalar for the surface yet
!
!.... reduce integrated scalar for each surface, push on scalInt
!
      npars=MAXSURF+1
      call MPI_ALLREDUCE (scalIntProc, scalInt(:), npars, &
                          MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)  
   
      return
      end subroutine

!
! *** New function to integrate velocity in the normal direction 
!     KDL and NAN 21/08/14
!
      subroutine GetFlowQ (qsurf, y, srfIdList, numSrfs)

      use pvsQbi  ! brings in NABI
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF), qsurfProc(0:MAXSURF)
      integer numSrfs, srfIdList(0:MAXSURF)
!
! note we only need the first three entries (u) from y

      if (numSrfs .eq. 0) return

      qsurfProc=zero
      do i = 1, nshg          
        do k = 1, numSrfs            
            if (srfIdList(k) .eq. ndsurf(i)) then               
               do j = 1, 3
                  qsurfProc(k) = qsurfProc(k) + NABI(i,j)*y(i,j)
               enddo
            else if (srfIdList(k) .eq. edsurf(i)) then               
               do j = 1, 3
                  qsurfProc(k) = qsurfProc(k) + NABI(i,j)*y(i,j)
               enddo
            endif
          enddo
      enddo
!
!     at this point, each qsurf has its "nodes" contributions to Q
!     accumulated into qsurf. Note, because NABI is on processor this
!     will NOT be Q for the surface yet
!
!.... reduce integrated Q for each surface, push on qsurf
!
      npars=MAXSURF+1
      call MPI_ALLREDUCE (qsurfProc, qsurf, npars, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
!
!.... return
!
      return
      end subroutine

      end module boundarymodule