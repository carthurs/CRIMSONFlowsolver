module measureWallDistance
      
    implicit none
     
    !     AABB tree child
    type AABB
        !       bounding box xyz coordinates
        real*8, dimension(3) :: bmin,bmax
        
        !       maximum radius in this AABB
        real*8 maxr
        
        !       number of triangles in this AABB
        integer nTris
        
        !       number of children in subtree of this AABB
        integer nChild

        !       array of triangle indices
        integer, dimension(:), allocatable :: tris
        
        !       array of child indices
        integer, dimension(:), allocatable :: child
    end type AABB
      
    !     one level of the octree
    type AABBlevel
        type(AABB), dimension(:), allocatable :: bb
    end type AABBlevel

    !     define structure for a single data frame
    type dataframe
      
        integer :: nPoints,nTris,nEdges
        
        !       xyz coordinates of triangle vertices
        real*8, dimension(:,:), allocatable :: coords
        
        !       pseudonormals for points (see Baerentzen et al)
        real*8, dimension(:,:), allocatable :: pntNrmls
        
        !       face normals
        real*8, dimension(:,:), allocatable :: faceNrmls
        
        !       pseudonormals for edges
        real*8, dimension(:,:,:), allocatable :: edgeNrmls
        
        !       array of triangle connectivity
        integer, dimension(:,:), allocatable :: tris
        
        !       triangle bounding spheres center and radius
        real*8, dimension(:,:), allocatable :: sphc
        real*8, dimension(:), allocatable :: sphr
        
        !       arrays to represent queue of AABB or triangles
        !       key array
        real*8, dimension(:), allocatable :: qA
        !       index, array
        integer, dimension(:), allocatable :: qI
         
        !       array of AABBs to check
        integer, dimension(:), allocatable :: toChBB
        
        !       array of triangles to check
        integer, dimension(:), allocatable :: toChTri
                
        !       AABBs
        type(AABBlevel), dimension(:), allocatable :: octree
                
        !       number of AABBs per level of the tree
        integer, dimension(:), allocatable :: nPerLev
        integer, dimension(:), allocatable :: nPerLevG
        
        !       how many levels in octree
        integer nLev       
        
        !       size of queues
        integer qSize         
                
    end type dataframe
      
    !     data structure for multiple data frames
    type(dataframe), dimension(:), allocatable, target :: dtfrs
      
    !     the total number of data frames
    integer numDataFrames
      
    !     length of the cardiac cycle
    real*8 cycleLength
      
contains
      
    ! /////////////////////////////////////////////////////////////////////////
    !  initialization function
    ! /////////////////////////////////////////////////////////////////////////

    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine dm_initialize
    ! load data and make octrees
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine dm_initialize()
      
        implicit none

        integer kk
      
        do kk=1,numDataFrames

            call sphbmesh(kk)
            call makebb(kk,5)

        end do
      
    end subroutine

      
    ! /////////////////////////////////////////////////////////////////////////
    !  functions dealing with one vector in R-3 and a triangular mesh
    ! /////////////////////////////////////////////////////////////////////////
      
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine dm_cpmeshp3
    ! Find a closest point on the triangular mesh
    ! to the specified closestpnt.
    ! Compute the signed distance from
    ! the specified point to a closest point on mesh.
    ! Only called after mesh is loaded.
    ! Uses octree with bounding spheres
    ! Only returns one of possibly many closest points
    ! since we are only interested in distance
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine dm_cpmeshp3(currFr,p,closestpnt,distanceSq,nv)

        implicit none

        integer currFr
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv 
        real*8 distanceSq   
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      currFr : current frame number
        !      p : the reference closestpnt

        ! OUTPUTS
        !      closestpnt : a closest point
        !      distanceSq: distance^2 to a closest point on mesh to p
        !      nv : normal to the surface at the closest point
        ! -------------------------------------------

        integer ii,jj,kk,mm,nn ! counters
        
        integer qBack           ! back of queue
        integer nCheck,nCheckTri
        
        real*8 ub,lb ! upper bound and lower bound
        real*8 tlb,tub ! temporary lower, upper bound
        
        integer, pointer :: nLev
        
        type(AABB), pointer :: AABBcurr ! pointer to an AABB
        
        integer, allocatable :: tempQ(:);

        integer duplicateFlag

        nLev=>dtfrs(currFr)%nLev ! pointers for easier syntax
        
        !       initialize the checking array of AABB indices
        nCheck = 8
        
        do kk=1,nCheck
            dtfrs(currFr)%toChBB(kk)=kk
        end do
        
        !       start at the uppermost level of the tree
        !       compute upper and lower bounds based on distances to the AABB
        do kk=2,nLev-1

            !         pointer to the first BB
            AABBcurr=>dtfrs(currFr)%octree(1)%bb(1)
        
            !         initialize the bounds
            call fpaabb(p,AABBcurr,ub)
            lb = ub+AABBcurr%maxr
            ub = ub+100
          
            !         reset priority queue
            qBack = 0
           
            !         insert the first BB in the priority queue
            call pqInsert(dtfrs(currFr)%qI, &
                dtfrs(currFr)%qA, &
                qBack, &
                1,lb)
          
            !         loop through all the BBs to check on this level
            do jj=1,nCheck
          
                !           get index of next BB to check
                mm = dtfrs(currFr)%toChBB(jj)
            
                !           pointer to the BB
                AABBcurr=>dtfrs(currFr)%octree(kk)%bb(mm)

                !           obtain the nearest and farthest possible distance to any point
                !           in the AABB
                call cpaabb(p,AABBcurr,tlb)
                call fpaabb(p,AABBcurr,tub)
            
                !           add the maximum radius inside the AABB
                tlb = tlb-AABBcurr%maxr
                tub = tub+AABBcurr%maxr
                        
                !           prune the AABBs that are too far away
                do while(dtfrs(currFr)%qA(1) .ge. ub)
                    call pqRemoveMax(dtfrs(currFr)%qI, &
                        dtfrs(currFr)%qA, &
                        qBack)
                end do

                if (tlb .lt. ub) then

                    !             add AABB to the queue
                    call pqInsert(dtfrs(currFr)%qI, &
                        dtfrs(currFr)%qA, &
                        qBack, &
                        mm,tlb)

                    !             update the bounds if they become smaller
                    if (tlb .lt. lb) then
                        lb = tlb
                    end if
                    if (tub .lt. ub) then
                        ub = tub
                    end if
              
                end if

            end do
                    
            !         copy the BB children indices in the queue to the
            !         array of BBs to check
            nCheck = 0

            allocate(tempQ(qBack))
            tempQ = 0

            do jj=1,qBack
                AABBcurr=> &
                    dtfrs(currFr)%octree(kk)%bb(dtfrs(currFr)%qI(jj))

                tempQ(jj) = dtfrs(currFr)%qI(jj)

                duplicateFlag = 0
                do nn=1,jj-1
                    if (dtfrs(currFr)%qI(jj).eq.tempQ(nn)) then
                        duplicateFlag = 1
                        exit
                    end if
                end do
     
                if (duplicateFlag.ne.1) then
                    do ii=1,AABBcurr%nChild
                        nCheck = nCheck + 1

                        !if (nCheck.gt.dtfrs(currFr)%nPerLev( dtfrs(currFr)%nLev)) then
                        !  nCheck = dtfrs(currFr)%nPerLev( dtfrs(currFr)%nLev)
                        !end if

                        dtfrs(currFr)%toChBB(nCheck)=AABBcurr%child(ii)
                    end do
                end if

            end do
          
            if (allocated(tempQ)) then
              deallocate(tempQ)
            endif

        end do
                          
        !       copy the triangle indices contained in the leftover BBs
        nCheckTri = 0

        do kk=1,nCheck

            AABBcurr=> &
                dtfrs(currFr)%octree(nLev)%bb(dtfrs(currFr)%toChBB(kk))
     
            do jj=1,AABBcurr%nTris
                nCheckTri = nCheckTri + 1
                dtfrs(currFr)%toChTri(nCheckTri) = AABBcurr%tris(jj)
            end do

        end do
        
        !        write(*,*) dtfrs(currFr)%toCheck(1:nCheckTri)
      
        call cpmeshpSph(currFr, &
            dtfrs(currFr)%toChTri(1:nCheckTri),nCheckTri, &
            p,closestpnt,distanceSq,nv)

    end subroutine

    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine cpmeshpSph
    ! Find a closest point on the triangular mesh
    ! to the specified closestpnt.
    ! Compute the signed distance from
    ! the specified point to a closest point on mesh.
    ! Only called after mesh is loaded.
    ! Uses brute force with bounding spheres
    ! Only returns one of possibly many closest points
    ! since we are only interested in distance
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine dm_cpmeshp2(currFr,p,closestpnt,distanceSq,nv)
        
        implicit none

        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        real*8 distanceSq   

        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      currFr : current frame number
        !      p : the reference closestpnt

        ! OUTPUTS
        !      closestpnt : a closest point
        !      distanceSq: distance^2 to a closest point on mesh to p
        ! -------------------------------------------

        integer kk
       
        do kk=1,dtfrs(currFr)%nTris
            dtfrs(currFr)%toChTri(kk)=kk
        end do

        call cpmeshpSph(currFr, &
            dtfrs(currFr)%toChTri, &
            dtfrs(currFr)%nTris, &
            p,closestpnt,distanceSq,nv)
     
    end subroutine

    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine cpmeshpSph
    ! Find a closest point on the triangular mesh
    ! to the specified closestpnt.
    ! Compute the signed distance from
    ! the specified point to a closest point on mesh.
    ! Only called after mesh is loaded.
    ! Uses brute force with bounding spheres
    ! Only returns one of possibly many closest points
    ! since we are only interested in distance
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine cpmeshpSph(currFr,toCheck,nCheck, &
        p,closestpnt,distanceSq,nv)

        implicit none
        
        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        integer, dimension(:) :: toCheck
        integer nCheck  
        real*8 distanceSq   

        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      currFr : current frame number
        !      p : the reference closestpnt
        !      toCheck : array of triangles to check
        !      nCheck : number of triangles to check

        ! OUTPUTS
        !      closestpnt : a closest point
        !      distanceSq: distance^2 to a closest point on mesh to p
        !      nv : normal to the surface at the closest point
        ! -------------------------------------------

        real*8 currDist ! current distance
        real*8, dimension(3) :: currCP,rayToRef
        integer loc
          
        integer ii,jj,kk ! counters
        
        integer qBack ! back of the queue
        
        real*8 ub,lb ! upper bound and lower bound
        real*8 tlb,tub ! temporary lower, upper bound
        
        qBack = 0

        !       initially the index of triangle to be checked carefully is
        !       going to be the first triangle in the sequence
        jj=dtfrs(currFr)%toChTri(1)
            
        !       get the initial smallest lower and upper bounds
        ub = sqrt(dot_product(p-dtfrs(currFr)%sphc(:,jj), &
            p-dtfrs(currFr)%sphc(:,jj))) &
            +dtfrs(currFr)%sphr(jj)
        lb = ub-2*dtfrs(currFr)%sphr(jj)
        
        !       store the lower bound for this triangle
        !       in the priority queue
        call pqInsert(dtfrs(currFr)%qI, &
            dtfrs(currFr)%qA, &
            qBack, &
            jj,lb)
        
        !       loop through the triangles
        do jj=2,nCheck
        
            kk = dtfrs(currFr)%toChTri(jj)

            !         compute lower and upper bounds for current triangle
            tub =sqrt(dot_product(p-dtfrs(currFr)%sphc(:,kk), &
                p-dtfrs(currFr)%sphc(:,kk))) &
                +dtfrs(currFr)%sphr(kk)
          
            tlb = tub-2*dtfrs(currFr)%sphr(kk)
          
            !         prune triangles from queue that we know are too far away
            do while(dtfrs(currFr)%qA(1) .ge. ub)
                call pqRemoveMax(dtfrs(currFr)%qI, &
                    dtfrs(currFr)%qA, &
                    qBack)
            end do
          
            !         if current lower bound is lower than established upper bound
            if (tlb .lt. ub) then
  
                !           add triangle to the queue
                call pqInsert(dtfrs(currFr)%qI, &
                    dtfrs(currFr)%qA,qBack, &
                    kk,tlb)
         
                !           update bounds
                if (tlb .lt. lb) then
                    lb = tlb ! update smallest lower bound
                end if
              
                if (tub .lt. ub) then
                    ub = tub ! update smallest upper bound
                end if

            end if

        end do        
        
        !        write(*,*) dtfrs(currFr)%qI(1:qBack)
        !        write(*,*) dtfrs(currFr)%qA(1:qBack)
        
        !      loop through the triangles to check
        call cpmeshp(currFr,dtfrs(currFr)%qI(1:qBack),qBack, &
            p,closestpnt,distanceSq,nv)
        
    end subroutine
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine cpmeshp1
    ! Find a closest point on the triangular mesh
    ! to the specified closestpnt.
    ! Compute the signed distance from
    ! the specified point to a closest point on mesh.
    ! Only called after mesh is loaded.
    ! Uses naive brute-force search.
    ! Only returns one of possibly many closest points
    ! since we are only interested in distance
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine dm_cpmeshp1(currFr,p,closestpnt,distanceSq,nv)
        
        implicit none

        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        real*8 distanceSq 
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      currFr : current frame number
        !      p : the reference closestpnt

        ! OUTPUTS
        !      closestpnt : a closest point
        !      distanceSq: distance^2 to a closest point on mesh to p
        !      nv : normal to the surface at the closest point
        ! -------------------------------------------

        integer kk
       
        do kk=1,dtfrs(currFr)%nTris
            dtfrs(currFr)%toChTri(kk)=kk
        end do

        call cpmeshp(currFr, &
            dtfrs(currFr)%toChTri, &
            dtfrs(currFr)%nTris, &
            p,closestpnt,distanceSq,nv)
          
    end subroutine
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine cpmeshp
    ! Find a closest point on the triangular mesh
    ! to the specified closestpnt.
    ! Compute the signed distance from
    ! the specified point to a closest point on mesh.
    ! Only called after mesh is loaded.
    ! Uses naive brute-force search.
    ! Specify the triangles to check.
    ! Only returns one of possibly many closest points
    ! since we are only interested in distance
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine cpmeshp(currFr,toCheck,nCheck, &
        p,closestpnt,distanceSq,nv)

        implicit none
        
        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        integer, dimension(:) :: toCheck
        integer nCheck 
        real*8 distanceSq 
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      currFr : current frame number
        !      p : the reference closestpnt
        !      toCheck : array of triangles to check
        !      nCheck : number of triangles to check

        ! OUTPUTS
        !      closestpnt : a closest point
        !      distanceSq: distance^2 to a closest point on mesh to p
        !      nv : normal to the surface at the closest point
        ! -------------------------------------------

        real*8 currDist ! current distance
        real*8, dimension(3) :: currCP,rayToRef
        integer loc
          
        integer jj,kk ! counters
          
        !       loop through the triangles
        do jj=1,nCheck        
           
            kk = toCheck(jj)
        
            call cptrip(p, &
                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(1,kk)), &
                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(2,kk)), &
                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(3,kk)), &
                currCP,loc)
     
            rayToRef = p-currCP;
    
            !         compute distance
            currDist = sum(rayToRef**2,1)
            
            !         determine which normal vector to use
            select case(loc)
                case(1)
                    nv = dtfrs(currFr)%pntNrmls(:,dtfrs(currFr)%tris(1,kk))
                case(2)
                    nv = dtfrs(currFr)%pntNrmls(:,dtfrs(currFr)%tris(2,kk))
                case(3)
                    nv = dtfrs(currFr)%pntNrmls(:,dtfrs(currFr)%tris(3,kk))
                case(4)
                    nv = dtfrs(currFr)%edgeNrmls(:,1,kk)
                case(5)
                    nv = dtfrs(currFr)%edgeNrmls(:,2,kk)
                case(6)
                    nv = dtfrs(currFr)%edgeNrmls(:,3,kk)
                case(7)
                    nv = dtfrs(currFr)%faceNrmls(:,kk)
                case default
            end select
            
            !         determine the sign of the distance
            currDist = sign(currDist,dot_product(rayToRef,nv))
     
            !         update distance and closest point
            if (abs(currDist) .lt. abs(distanceSq) .or. jj .eq. 1) then
                distanceSq = currDist
                closestpnt = currCP
            !            write(*,*) kk
            end if
        end do             
        
    end subroutine
        
    ! /////////////////////////////////////////////////////////////////////////
    !  functions dealing with building spatial heirarchy
    ! /////////////////////////////////////////////////////////////////////////
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine makebb
    ! Construct the AABB octree for the mesh
    ! Assumes that triangle bounding volumes
    ! have been constructed already
    ! Assume at least one level of refinement
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine makebb(currFr,refine)
        
        implicit none

        integer currFr,refine
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      currFr : current frame number
        !      refine : desired level of refinement
        ! -------------------------------------------
        
        integer ii,jj,kk
        integer parent
        integer nChildG
        
        integer, dimension(:), allocatable :: temp ! temporary array
        
        type(AABB), dimension(:), allocatable :: tempLev
        
        integer, dimension(3) :: xyzInd
        
        real*8, dimension(3) :: pmin,pmax,tmin,tmax,pmid,tcent
        
        real*8, dimension(3) :: boundDim,boxDim
                
        real*8 maxr
        
        integer, dimension(:,:,:), allocatable :: boxHere
               
        integer, pointer :: nPerLev, nPerLevG
        
        type(AABB), pointer :: AABBcurr,AABBnew
                 
        !       allocate enough space
        allocate(dtfrs(currFr)%nPerLev(1+refine))
        allocate(dtfrs(currFr)%nPerLevG(1+refine))
        allocate(dtfrs(currFr)%octree(1+refine))
        
        !       set the box presence grid to all zeros
        allocate(boxHere(2**refine,2**refine,2**refine))
        boxHere = 0
        
        !       set initial array sizes at each level of the tree
        do kk=1,refine+1
        
            dtfrs(currFr)%nPerLev(kk) = 0
            dtfrs(currFr)%nPerLevG(kk) = 8
          
        end do
        
        !       compute initial bounding box
        !       look for the maximum extent of x,y,z coordinates of the
        !       triangle bounding spheres
        pmin = dtfrs(currFr)%sphc(:,1)
        pmax = pmin
                
        maxr = dtfrs(currFr)%sphr(1)

        do kk = 2,dtfrs(currFr)%nTris
        
            if (dtfrs(currFr)%sphc(1,kk) .lt. pmin(1)) then
                pmin(1) = dtfrs(currFr)%sphc(1,kk)
            end if
          
            if (dtfrs(currFr)%sphc(1,kk) .gt. pmax(1)) then
                pmax(1) = dtfrs(currFr)%sphc(1,kk)
            end if
          
            if (dtfrs(currFr)%sphc(2,kk) .lt. pmin(2)) then
                pmin(2) = dtfrs(currFr)%sphc(2,kk)
            end if
          
            if (dtfrs(currFr)%sphc(2,kk) .gt. pmax(2)) then
                pmax(2) = dtfrs(currFr)%sphc(2,kk)
            end if
          
            if (dtfrs(currFr)%sphc(3,kk) .lt. pmin(3)) then
                pmin(3) = dtfrs(currFr)%sphc(3,kk)
            end if
          
            if (dtfrs(currFr)%sphc(3,kk) .gt. pmax(3)) then
                pmax(3) = dtfrs(currFr)%sphc(3,kk)
            end if
          
            if (dtfrs(currFr)%sphr(kk) .gt. maxr) then
                maxr = dtfrs(currFr)%sphr(kk)
            end if
          
        end do

        !        write(*,*) "Bounds"
        !        write(*,*) pmin
        !        write(*,*) pmax
        !        write(*,*) "Max radius"
        !        write(*,*) maxr

        !       these are the dimensions of the whole bounding volume
        boundDim = pmax-pmin
                
        !       these are the dimensions of the smallest AABB
        boxDim = boundDim / 2**refine
        
        !        write(*,*) "Box dim"
        !        write(*,*) boxDim

        !       pointers for easier access
        nPerLev => dtfrs(currFr)%nPerLev(1+refine)
        nPerLevG => dtfrs(currFr)%nPerLevG(1+refine)
        
        !       allocate enough space for each level
        allocate(tempLev(nPerLevG))
        allocate(dtfrs(currFr)%octree(1+refine)%bb(nPerLevG))
        
        !       define the number of levels in the tree to be one greater
        !       than the refine parameter
        dtfrs(currFr)%nLev = refine+1        
        
        !       we loop through each triangle and assign them to the correct
        !       AABB at the bottommost level
        do kk=1,dtfrs(currFr)%nTris
        
            !         get the "xyz index" of the sphere
            !         these are just integer values which will
            !         help us define the bounds for the AABB
            xyzInd = ceiling((dtfrs(currFr)%sphc(:,kk)-pmin)/boxdim)
          
            do ii=1,3
          
                if (xyzInd(ii) .eq. 0) then
                    xyzInd(ii) = 1
                end if
            
                if (xyzInd(ii) .eq. 2**refine+1) then
                    xyzInd(ii) = 2**refine
                end if
            
            end do
          
            !         if no previous AABB has been made in this spatial location
            !         then we make one
            if (boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) .eq. 0) then
          
                !           compute the center of the box
                tcent = xyzInd*boxDim-0.5*boxdim+pmin
                !           compute the bounds of the box
                tmin = tcent-0.5*boxdim
                tmax = tcent+0.5*boxdim
                    
                !           increment the number of AABBs at his level
                nPerLev = nPerLev + 1
            
                !           we store the number of the new AABB
                !           by doing this, we make sure that no other AABBs can
                !           be created at this spatial location
                boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) = nPerLev
              
                !            write(*,*) xyzInd
              
                !           resize the AABB array at this level
                !           if necessary
                if (nPerLev .gt. nPerLevG) then
              
                    nPerLevG = nPerLevG + 8!**(refine-1)
              
                    !              write(*,*) nPerLevG
              
                    if (allocated(tempLev)) then
                        deallocate(tempLev)
                    end if
              
                    allocate(tempLev(nPerLevG))
              
                    do ii=1,nPerLev-1
                        tempLev(ii)=dtfrs(currFr)%octree(1+refine)%bb(ii)
                    end do
              
                    if (allocated(dtfrs(currFr)%octree(1+refine)%bb)) then
                        deallocate(dtfrs(currFr)%octree(1+refine)%bb)
                    end if
          
                    allocate(dtfrs(currFr)%octree(1+refine)%bb(nPerLevG))
          
                    do ii=1,nPerLev-1
                        dtfrs(currFr)%octree(1+refine)%bb(ii)=tempLev(ii)
                    end do
            
                end if
                        
                !           pointer to the new box for easier access
                AABBcurr=>dtfrs(currFr)%octree(1+refine)%bb(nPerLev)
            
                !           assign bounds
                AABBcurr%bmin = tmin
                AABBcurr%bmax = tmax
            
                !           this AABB has a triangle count of 1
                AABBcurr%nTris = 1

                !           since we are at the bottommost level,
                !           this AABB has no children
                AABBcurr%nChild = 0
            
                !           allocate space for the triangle array
                allocate(AABBcurr%tris(1))
            
                !           assign the triangle index
                AABBcurr%tris(1) = kk
            
                !           the maximum radius in this AABB
                AABBcurr%maxr = dtfrs(currFr)%sphr(kk)
 
            else
            
                !           there is already an AABB at this spatial location
                !           we merely have to find its index
                AABBcurr=>dtfrs(currFr)%octree(1+refine)%bb( &
                    boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) )
     
                !           increment triangle count
                AABBcurr%nTris = AABBcurr%nTris+1
     
                !           resize the triangle array
                allocate(temp(AABBcurr%nTris-1))
            
                do ii=1,AABBcurr%nTris-1
            
                    temp(ii)=AABBcurr%tris(ii)
              
                end do
            
                if (allocated(AABBcurr%tris)) then
                    deallocate(AABBcurr%tris)
                end if
            
                allocate(AABBcurr%tris(AABBcurr%nTris))
            
                do ii=1,AABBcurr%nTris-1
            
                    AABBcurr%tris(ii)=temp(ii)
              
                end do
            
                if (allocated(temp)) then
                    deallocate(temp)
                end if
            
                !           assign the triangle index
                AABBcurr%tris(AABBcurr%nTris) = kk
            
                !           if the sphere bounding the triangle
                !           has a larger radius than the current
                !           maximum radius in the AABB, then it becomes the
                !           new max radius in the AABB
                if (dtfrs(currFr)%sphr(kk) .gt. AABBcurr%maxr) then
                    AABBcurr%maxr = dtfrs(currFr)%sphr(kk)
                end if

            end if
                   
        end do
        
        !do ii=1,dtfrs(currFr)%nPerLev(3)
        !  write(*,*) dtfrs(currFr)%octree(3)%bb(ii)%nTris
        !end do
        
        !       deallocate arrays
        if (allocated(tempLev)) then
            deallocate(tempLev)
        end if
        
        if (allocated(boxHere)) then
            deallocate(boxHere)
        end if
        
        !       we have now created the AABBs that bound the triangle volumes
        !       at the bottommost level of the tree
        !       now we go up the levels
        do kk=refine,2,-1
        
            !       reallocate the box presence grid
            allocate(boxHere(2**(kk-1),2**(kk-1),2**(kk-1)))
        
            !       initialize the box presence grid to zero
            boxHere = 0
        
            !       compute the dimensions of the AABBs in this level
            boxDim = boundDim / 2**(kk-1)
        
            !       pointers for easier access
            nPerLev=>dtfrs(currFr)%nPerLev(kk)
            nPerLevG=>dtfrs(currFr)%nPerLevG(kk)
        
            !       allocate enough space for the level
            allocate(tempLev(nPerLevG))
            allocate(dtfrs(currFr)%octree(kk)%bb(nPerLevG))
        
            !       loop through the AABBs at this level
            do jj=1,dtfrs(currFr)%nPerLev(1+kk)

                !         get the "xyz index" of the AABB
                !         just like we did for the triangle volumes
                AABBcurr=>dtfrs(currFr)%octree(1+kk)%bb(jj)
                tcent = 0.5*(AABBcurr%bmin+AABBcurr%bmax)
                xyzInd = ceiling((tcent-pmin)/boxdim)
                do ii=1,3
          
                    if (xyzInd(ii) .eq. 0) then
                        xyzInd(ii) = 1
                    end if
            
                    if (xyzInd(ii) .eq. 2**(kk-1)+1) then
                        xyzInd(ii) = 2**(kk-1)
                    end if
            
                end do
          
                !         if no AABB exists here, then create one
                if (boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) .eq. 0) then
            
                    !           compute center of AABB
                    tcent = xyzInd*boxDim-0.5*boxdim+pmin
            
                    !           compute bounds of the AABB
                    tmin = tcent-0.5*boxdim
                    tmax = tcent+0.5*boxdim
                    
                    !           increment the AABB array at this level by one
                    nPerLev = nPerLev + 1
            
                    !           store AABB index
                    boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) = nPerLev
            
                    !           resize AABB array at this level
                    if (nPerLev .gt. nPerLevG) then
            
                        nPerLevG = nPerLevG + 8**(kk-2)
            
                        if (allocated(tempLev)) then
                            deallocate(tempLev)
                        end if
              
                        allocate(tempLev(nPerLevG))
              
                        do ii=1,nPerLev-1
                            tempLev(ii)=dtfrs(currFr)%octree(kk)%bb(ii)
                        end do
              
                        if (allocated(dtfrs(currFr)%octree(kk)%bb)) then
                            deallocate(dtfrs(currFr)%octree(kk)%bb)
                        end if
              
                        allocate(dtfrs(currFr)%octree(kk)%bb(nPerLevG))
              
                        do ii=1,nPerLev-1
                            dtfrs(currFr)%octree(kk)%bb(ii)=tempLev(ii)
                        end do
            
                    end if
            
                    !           pointer to the new AABB
                    AABBnew=>dtfrs(currFr)%octree(kk)%bb(nPerLev)
            
                    !           assign bounds
                    AABBnew%bmin = tmin
                    AABBnew%bmax = tmax
            
                    !           assign number of triangles
                    !           inherit this value from the AABB encompassed
                    !           in the lower level
                    AABBnew%nTris = AABBcurr%nTris

                    !           same for the max radius
                    AABBnew%maxr = AABBcurr%maxr
            
                    !           set number of children to 1
                    AABBnew%nChild = 1
            
                    !           allocate space for the child array
                    allocate(AABBnew%child(1))
            
                    !           assign index of the AABB in the lower level
                    AABBnew%child(1) = jj
            
                else
          
                    !           get the pointer to the box already here
                    AABBnew=>dtfrs(currFr)%octree(kk)%bb( &
                        boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) )
            
                    !           increase triangle count
                    !           add the triangle count of the AABB in the lower level
                    !           to the triangle count of the current AABB
                    AABBnew%nTris = AABBnew%nTris+AABBcurr%nTris
            
                    !           update maximum radius
                    if (AABBcurr%maxr .gt. AABBnew%maxr) then
                        AABBnew%maxr = AABBcurr%maxr
                    end if
            
                    !           increment number of children
                    AABBnew%nChild = AABBnew%nChild + 1
            
                    !           resize child array
                    allocate(temp(AABBnew%nChild-1))
            
                    do ii=1,AABBnew%nChild-1
                        temp(ii)=AABBnew%child(ii)
                    end do
            
                    if (allocated(AABBnew%child)) then
                        deallocate(AABBnew%child)
                    end if
                        
                    allocate(AABBnew%child(AABBnew%nChild))
            
                    do ii=1,AABBnew%nChild-1
                        AABBnew%child(ii)=temp(ii)
                    end do
            
                    if (allocated(temp)) then
                        deallocate(temp)
                    end if
            
                    !           assign index of child
                    AABBnew%child(AABBnew%nChild) = jj
            
                end if
          
            end do
        
            !       deallocate arrays
            if (allocated(tempLev)) then
                deallocate(tempLev)
            end if
        
            if (allocated(boxHere)) then
                deallocate(boxHere)
            end if
        
        
        end do
        
        !       write the information
        !       for the level 1 AABB which just encompasses everything
        !       this is for mostly debugging
        !       since we can check if the max radius and number of triangles
        !       was correctly propagated up the tree
        allocate(dtfrs(currFr)%octree(1)%bb(1))        
        dtfrs(currFr)%nPerLev(1) = 1
             
        AABBnew=>dtfrs(currFr)%octree(1)%bb(1)

        AABBnew%bmin = tmin
        AABBnew%bmax = tmax
        AABBnew%nTris = 0
        AABBnew%maxr = 0
        AABBnew%nChild = dtfrs(currFr)%nPerLev(2)
        
        allocate(AABBnew%child(dtfrs(currFr)%nPerLev(2)))
        
        do jj=1,dtfrs(currFr)%nPerLev(2)
        
            AABBcurr=>dtfrs(currFr)%octree(2)%bb(jj)
            AABBnew%nTris = AABBnew%nTris+AABBcurr%nTris
          
            if (AABBcurr%maxr .gt. AABBnew%maxr) then
                AABBnew%maxr = AABBcurr%maxr
            end if
          
            AABBnew%child(jj) = jj
          
        end do
        
        !       temporary the toCheck array AABBs
        allocate(dtfrs(currFr)%toChBB( dtfrs(currFr)%nPerLev( dtfrs(currFr)%nLev) ) )
     
        !       find the proper size for allocation of the queues
        dtfrs(currFr)%qSize = dtfrs(currFr)%nTris
        
        if ( dtfrs(currFr)%nPerLev(dtfrs(currFr)%nLev) .gt. &
            dtfrs(currFr)%qSize ) then
     
            dtfrs(currFr)%qSize = &
                dtfrs(currFr)%nPerLev(dtfrs(currFr)%nLev)
     
        end if
        
        !       allocate space for priority queue
        allocate(dtfrs(currFr)%qA(dtfrs(currFr)%qSize))
        allocate(dtfrs(currFr)%qI(dtfrs(currFr)%qSize))
                
    !       for debugging purposes
        
    !        write(*,*) AABBnew%nTris
    !        write(*,*) AABBnew%maxr
        
    !        do kk=1,refine+1
    !        write(*,*) "Level",kk
    !        do jj=1,dtfrs(currFr)%nPerLev(kk)
    !          AABBcurr=>dtfrs(currFr)%octree(kk)%bb(jj)
    !          write(*,*) jj,AABBcurr%nChild,AABBcurr%nTris
    !        end do
    !        end do

    end subroutine
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine sphbmesh
    ! Compute bounding spheres for triangular mesh
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine sphbmesh(currFr)
        
        implicit none

        integer currFr
        integer kk
        
        allocate( dtfrs(currFr)%sphc &
            (3,dtfrs(currFr)%nTris) )
        allocate( dtfrs(currFr)%sphr &
            (dtfrs(currFr)%nTris) )
     
        !       loop through the triangles
        do kk=1,dtfrs(currFr)%nTris
            call sphbt1(dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(1,kk)), &
                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(2,kk)), &
                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(3,kk)), &
                dtfrs(currFr)%sphc(:,kk),dtfrs(currFr)%sphr(kk) )
        end do
        
    end subroutine
      	  
    ! /////////////////////////////////////////////////////////////////////////
    !  functions dealing with AABB
    ! /////////////////////////////////////////////////////////////////////////
	  
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine cpaabb
    ! Compute distance to AABB
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine cpaabb(p,box,distance)
        
        implicit none

        real*8, dimension(3) :: p
        type(AABB) box
        real*8 distance,v
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      p : point p
        !      box : AABB

        ! OUTPUTS
        !      distanceSq  : distance squared
        ! -------------------------------------------
        
        integer ii
        
        distance = 0.0;
        
        do ii=1,3
            v = p(ii)
          
            if (v .lt. box%bmin(ii)) then
                distance = distance + (box%bmin(ii)-v)*(box%bmin(ii)-v)
            end if
            if (v .gt. box%bmax(ii)) then
                distance = distance + (v-box%bmax(ii))*(v-box%bmax(ii))
            end if
          
        end do
        
        distance=sqrt(distance)
        
    end subroutine
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine fpaabb
    ! Compute farthest distance to AABB
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine fpaabb(p,box,distance)
        
        implicit none

        real*8, dimension(3) :: p,fp,cent
        type(AABB) box
        real*8 distance
        integer, dimension(3) :: side
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      p : point p
        !      box : AABB

        ! OUTPUTS
        !      distanceSq  : distance squared
        ! -------------------------------------------

        distance = 0.0
        side = 1
        
        cent = 0.5*(box%bmin+box%bmax)

        if (p(1) .lt. cent(1) ) then
            side(1) = 0
        end if
        
        if (p(2) .lt. cent(2) ) then
            side(2) = 0
        end if
        
        if (p(3) .lt. cent(3) ) then
            side(3) = 0
        end if
        
        side = abs(side-1)
        
        fp = box%bmin+side*(box%bmax-box%bmin)
        
        distance = sqrt(sum((p-fp)*(p-fp)))
        
    end subroutine
        
    ! /////////////////////////////////////////////////////////////////////////
    !  functions dealing with triangles
    ! /////////////////////////////////////////////////////////////////////////

    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine cptrip
    ! Find the closest point on a triangle to
    ! specified point
    ! -------------------------------------------
    ! -------------------------------------------
    subroutine cptrip(p,a,b,c,closestpnt,loc) bind(C, name="cptrip")

        implicit none
       
        real*8, dimension(3) :: a 
        real*8, dimension(3) :: b 
        real*8, dimension(3) :: c 
        real*8, dimension(3) :: p 
        real*8, dimension(3) :: closestpnt 
        integer loc 

        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE

        !      p : the reference closestpnt
        !      a : point a of triangle
        !      b : point b of triangle
        !      c : point c of triangle

        ! OUTPUTS
        !      point: closest point on triangle to p
        !      loc: 1 denotes closest point is on triangle vertex a
        !            2 denotes closest point is on triangle vertex b
        !            3 denotes closest point is on triangle vertex c
        !            4 denotes closest point is on loc ab
        !            5 denotes closest point is on loc ac
        !            6 denotes closest point is on loc bc
        !            7 denotes closest point is on the face
        ! -------------------------------------------
         
        real*8, dimension(3) :: ab,ac,ap,bp,cp
         
        real*8 d1,d2,d3,d4,d5,d6,v,w,va,vb,vc,denom
         
        !       Check if p in vertex region outside A
        ab = b - a
        ac = c - a
        ap = p - a
        d1 = dot_product(ab,ap)
        d2 = dot_product(ac,ap)
        if (d1 .le. 0.0 .and. d2 .le. 0) then
            !         Barycentric coordinates (1,0,0)
            closestpnt = a
            loc = 1
            return
        end if
	   
        !       Check if p in vertex region outside B
        bp = p - b
        d3 = dot_product(ab,bp)   
        d4 = dot_product(ac,bp)
        if (d3 .ge. 0.0 .and. d4 .le. d3) then
            !         Barycentric coordinates (0,1,0)
            closestpnt = b
            loc = 2
            return
        end if
         
        !       Check if p in loc region of AB, if so return projection of P onto AB
        vc = d1*d4 - d3*D2
        if (vc .le. 0.0 .and. d1 .ge. 0.0 .and. d3 .le. 0.0) then
            v = d1 / (d1 - d3)
            closestpnt = a + v*ab
            loc = 4
            return
        end if
         
        !       Check if p in vertex region outside C
        cp = p - c
        d5 = dot_product(ab,cp)
        d6 = dot_product(ac,cp)
        if (d6 .ge. 0.0 .and. d5 .le. d6) then
            !          Barycentric coordinates (0,0,1)
            closestpnt = c
            loc = 3
            return
        end if
         
        !       Check if p in loc region of AC, if so return projection of P onto AC
        vb = d5*d2 - d1*d6
        if (vb .le. 0.0 .and. d2 .ge. 0.0 .and. d6 .le. 0.0) then
            w = d2 / (d2 - d6)
            closestpnt = a + w*ac
            loc = 5
            return
        end if
         
        !       Check if p in loc region of BC, if so return projection of P onto B!
        va = d3*d6 - d5*d4
        if (va .le. 0.0 .and. (d4 - d3) .ge. 0.0 .and. &
            (d5 - d6) .ge. 0.0) then
            w = (d4 - d3) / ( (d4 - d3) + (d5 - d6) )
            closestpnt = b + w*(c - b)
            loc = 6
            return
        end if
         
        !       p inside face region. compute closestpnt using barycentric coordinates
        denom = 1.0 / (va+vb+vc)
        v = vb*denom
        w = vc*denom
        closestpnt = a + ab*v + ac*w
        loc = 7
         
    end subroutine
                 
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine sphbt2
    ! Compute bounding sphere for triangle
    ! Version 2
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine sphbt2(a,b,c,gcent,radius)
        
        implicit none

        real*8, dimension(3) :: a,b,c,gcent
        real*8 radius,tr
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      a : point a
        !      b : point b
        !      c : point c
        ! OUTPUTS
        !      gcent  : center of the bounding sphere
        !      radius : radius of the bounding sphere
	  
        ! -------------------------------------------

        !      compute the centroid of the triangle
        gcent = (a+b+c)/3
       
        radius = 0.0
       
        tr = sqrt(sum((a-gcent)*(a-gcent)))
       
        if (tr .gt. radius) then
            radius = tr
        end if
       
        tr = sqrt(sum((b-gcent)*(b-gcent)))
       
        if (tr .gt. radius) then
            radius = tr
        end if
       
        tr = sqrt(sum((c-gcent)*(c-gcent)))
       
        if (tr .gt. radius) then
            radius = tr
        end if

    end subroutine
       
    ! -------------------------------------------
    ! subroutine sphbt1
    ! Compute bounding sphere for triangle
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine sphbt1(a,b,c,gcent,radius)
        
        implicit none

        real*8, dimension(3) :: a,b,c,gcent
        real*8 radius
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      a : point a
        !      b : point b
        !      c : point c
        ! OUTPUTS
        !      gcent  : center of the bounding sphere
        !      radius : radius of the bounding sphere
        ! -------------------------------------------
        
        real*8, dimension(3) :: al,bl,cl,cent,tr
        real*8, dimension(3,3) :: T
        real*8 D,la,lb,lc
        
        la = sqrt(dot_product(b-a,b-a))
        lb = sqrt(dot_product(c-a,c-a))
        lc = sqrt(dot_product(b-c,b-c))
        
        radius = la*lb*lc / &
            sqrt((la+lb+lc)*(-la+lb+lc)*(la-lb+lc)*(la+lb-lc))
        
        call rottriloc(a,b,c,al,bl,cl,T)
        
        !       ignore the z coordinate
        !       translate local point a to the origin
        tr = al
        al(1:2) = 0
        bl(1:2) = bl(1:2) - tr(1:2)
        cl(1:2) = cl(1:2) - tr(1:2)
        
        !       compute center of circumscribed circle
        D = 2*(bl(1)*cl(2)-bl(2)*cl(1))
        cent(1) = ( cl(2)*(bl(1)**2+bl(2)**2)- &
            bl(2)*(cl(1)**2+cl(2)**2) ) / D
        cent(2) = ( bl(1)*(cl(1)**2+cl(2)**2)- &
            cl(1)*(bl(1)**2+bl(2)**2) ) / D
        cent(3) = al(3)
        
        !       translate back
        cent(1:2) = cent(1:2) + tr(1:2)
              
        !       rotate back to global coordinate system
        gcent(1) = sum(T(:,1)*cent,1)
        gcent(2) = sum(T(:,2)*cent,1)
        gcent(3) = sum(T(:,3)*cent,1)    
        
    end subroutine

    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine rottriloc
    ! Rotate triangle to local 2D cartesian system
    ! -------------------------------------------
    ! -------------------------------------------
      
    subroutine rottriloc(a,b,c,al,bl,cl,T)
        
        implicit none

        real*8, dimension(3) :: a,b,c,al,bl,cl
        real*8, dimension(3,3) :: T
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !      a : point a
        !      b : point b
        !      c : point c
        ! OUTPUTS
        !      al : local point al
        !      bl : local point bl
        !      cl : local point cl
        !      T  : 3x3 transformation matrix
        ! -------------------------------------------

        real*8, dimension(3) :: v1,v2,v3

        !       compute local basis vectors
        v1 = b - a
        v1 = v1/sqrt(dot_product(v1,v1))
        v2 = c - a
        v3 = cross_product(v1, v2)
        v3 = v3/sqrt(dot_product(v3,v3))
        v2 = cross_product(v3, v1)
        
        !       compute transformation matrix
        T(1,:) = v1
        T(2,:) = v2
        T(3,:) = v3
        
        !       rotate to local cartesian coordinates
        al(1) = sum(T(1,:)*a,1)
        al(2) = sum(T(2,:)*a,1)
        al(3) = sum(T(3,:)*a,1)
        
        bl(1) = sum(T(1,:)*b,1)
        bl(2) = sum(T(2,:)*b,1)
        bl(3) = sum(T(3,:)*b,1)
        
        cl(1) = sum(T(1,:)*c,1)
        cl(2) = sum(T(2,:)*c,1)
        cl(3) = sum(T(3,:)*c,1)
                
    end subroutine
	  
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine dataread
    ! Read in the mesh data
    ! Call this first before cpmeshp
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine dm_dataread(filename)

        implicit none

        character(*) :: filename
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE

        !      filename : the name of the file to be read

        ! OUTPUTS
        !      no outputs
        ! -------------------------------------------

        character dummy
        integer dummyint
        
        integer currFr ! current data frame
        
        integer :: ii,jj,kk ! counters
        
        open(unit=2,file=filename)

        !       Skip the first identifier line and
        !       read the number of data frames
        read(unit=2,fmt='(A1)') dummy
        read(2,*) numDataFrames
        
        !       Skip the first identifier line and
        !       read the length of the cycle
        read(unit=2,fmt='(A1)') dummy
        read(2,*) cycleLength
        
        !       Allocate space for the data frames
        allocate(dtfrs(numDataFrames))     
        
        do kk=1,numDataFrames
            !         Skip the next identifier line and
            !         read the data frame number
            read(unit=2,fmt='(A1)') dummy
            read(2,*) currFr
        
            !         Skip the next identifier line and
            !         read the number of points
            read(unit=2,fmt='(A1)') dummy
            read(2,*) dtfrs(currFr)%nPoints

            !         Skip the next identifier line and
            !         read the number of triangles
            read(unit=2,fmt='(A1)') dummy
            read(2,*) dtfrs(currFr)%nTris

            !         Skip the next identifier line and
            !         read the number of edges
            read(unit=2,fmt='(A1)') dummy
            read(2,*) dtfrs(currFr)%nEdges

            !         Now allocate memory
            allocate( dtfrs(currFr)%coords &
                (3,dtfrs(currFr)%nPoints) )
     
            allocate( dtfrs(currFr)%tris &
                (3,dtfrs(currFr)%nTris) )
     
            allocate( dtfrs(currFr)%pntNrmls &
                (3,dtfrs(currFr)%nPoints) )
     
            allocate( dtfrs(currFr)%edgeNrmls &
                (3,3,dtfrs(currFr)%nTris) ) ! index by dof, edge, then tri
        
            allocate( dtfrs(currFr)%faceNrmls &
                (3,dtfrs(currFr)%nTris) )
             
            !         Skip the next identifier line and
            !         read the coordinates
            read(unit=2,fmt='(A1)') dummy
            do ii=1,dtfrs(currFr)%nPoints
                read(2,*) dtfrs(currFr)%coords(1,ii), &
                    dtfrs(currFr)%coords(2,ii), &
                    dtfrs(currFr)%coords(3,ii)
            end do
        
            !         Skip the next identifier line and
            !         read the surface triangle connectivity
            read(unit=2,fmt='(A1)') dummy
            do ii=1,dtfrs(currFr)%nTris
                read(2,*) dtfrs(currFr)%tris(1,ii), &
                    dtfrs(currFr)%tris(2,ii), &
                    dtfrs(currFr)%tris(3,ii)
            end do
        
            !         Skip the next identifier line and
            !         read the face normals
            read(unit=2,fmt='(A1)') dummy
            do ii=1,dtfrs(currFr)%nTris
                read(2,*) dtfrs(currFr)%faceNrmls(1,ii), &
                    dtfrs(currFr)%faceNrmls(2,ii), &
                    dtfrs(currFr)%faceNrmls(3,ii)
            end do
        
            !         Skip the next identifier line and
            !         read the edge normals
            read(unit=2,fmt='(A1)') dummy
            do ii=1,dtfrs(currFr)%nTris
                do jj=1,3
                    read(2,*) dummyint, &
                        dtfrs(currFr)%edgeNrmls(1,jj,ii), &
                        dtfrs(currFr)%edgeNrmls(2,jj,ii), &
                        dtfrs(currFr)%edgeNrmls(3,jj,ii)
                end do
            end do
        
            !         Skip the next identifier line and
            !         read the closestpnt normals
            read(unit=2,fmt='(A1)') dummy
            do ii=1,dtfrs(currFr)%nPoints
                read(2,*) dtfrs(currFr)%pntNrmls(1,ii), &
                    dtfrs(currFr)%pntNrmls(2,ii), &
                    dtfrs(currFr)%pntNrmls(3,ii)
            end do
          
            !         allocate enough space for toCheck
            !         for safety, we will allocate the with the number of triangles
          
            !         array of triangles to check
            allocate(dtfrs(currFr)%toChTri(dtfrs(currFr)%nTris))
          
        end do

        close(2)

    !        write(*,*) dummy
    !        write(*,*) numDataFrames
    !        write(*,*) dtfrs(2)%nPoints
    !        write(*,*) dtfrs(2)%nTris
    !        write(*,*) dtfrs(2)%nEdges
        
    !        do ii=1,dtfrs(currFr)%nTris
    !          write(*,*) dtfrs(currFr)%tris(1,ii),
    !     &               dtfrs(currFr)%tris(2,ii),
    !     &               dtfrs(currFr)%tris(3,ii)
    !        end do
      
    end subroutine
        
    ! /////////////////////////////////////////////////////////////////////////
    !  functions dealing with vectors in R-3
    ! /////////////////////////////////////////////////////////////////////////

    ! -------------------------------------------
    ! -------------------------------------------
    ! function cross_product
    ! computes cross product of two vectors
    ! -------------------------------------------
    ! -------------------------------------------

    function cross_product(a, b)
        
        implicit none

        real*8, dimension(3) :: a, b, cross_product

        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE

        !      a : 3x1
        !      b : 3x1

        ! OUTPUTS
        !      cross_product
        ! -------------------------------------------

        cross_product(1) = a(2) * b(3) - a(3) * b(2)
        cross_product(2) = a(3) * b(1) - a(1) * b(3)
        cross_product(3) = a(1) * b(2) - a(2) * b(1)
            
    end function cross_product
        
    ! /////////////////////////////////////////////////////////////////////////
    !  Functions that deal with priority queue insertion, removal
    ! /////////////////////////////////////////////////////////////////////////
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine pqInsert
    ! insert a new item in the priority queue
    ! -------------------------------------------
    ! -------------------------------------------

    subroutine pqInsert(qI,qA,qBack,ival,pval)
        
        implicit none

        integer, dimension(:) :: qI
        real*8, dimension(:) :: qA
        integer qBack,ival
        real*8 pval
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        ! NONE OF THESE ARE CHANGED BY THE SUBROUTINE
        !
        !      ival : index value
        !      pval : key value
        !
        ! THESE ARE CHANGED
        !
        !      qI : queue of indices
        !      qA : queue of keys
        !      qBack : index of the back of the queue
        !
        ! -------------------------------------------

        integer index,parInd    
        real*8 tmpSwpR
        integer tmpSwpI    
        
        !       add to last level

        !       add key
        qBack = qBack + 1
        qA(qBack) = pval
        
        !       add index value
        qI(qBack) = ival        
        
        index = qBack

        do while (index .gt. 1)
        
            !         compare with parent
            !         if necessary swap
            parInd = index/2
          
            if (qA(index) .le. qA(parInd)) then
                exit
            end if
          
            !         swap values of keys
            tmpSwpR = qA(parInd)
            qA(parInd) = qA(index)
            qA(index) = tmpSwpR
          
            !         swap values of indices
            tmpSwpI = qI(parInd)
            qI(parInd) = qI(index)
            qI(index) = tmpSwpI
          
            index = parInd
          
        end do
                 
    end subroutine
        
    ! -------------------------------------------
    ! -------------------------------------------
    ! subroutine pqRemoveMax
    ! removes the index in the priority queue associated
    ! with the largest key
    ! -------------------------------------------
    ! -------------------------------------------
      
    subroutine pqRemoveMax(qI,qA,qBack)

        implicit none
          
        integer, dimension(:) :: qI  
        real*8, dimension(:) :: qA
        integer qBack
        
        ! -------------------------------------------
        ! INPUT PARAMETERS
        !
        ! THESE ARE CHANGED
        !
        !      qI : queue of indices
        !      qA : queue of keys
        !      qBack : index of the back of the queue
        !
        ! -------------------------------------------

        integer c1Ind,c2Ind,lcInd,index      
        real*8 tmpSwpR
        integer tmpSwpI
        
        !       do not remove max if only one element left
        if (qBack .gt. 1) then
            !         assign the value of the last child to the root
            qA(1) = qA(qBack)
            qI(1) = qI(qBack)
            !         decrease size of queue
            qBack = qBack - 1
        end if
        
        index = 1
        c1Ind = index*2
        c2Ind = index*2 + 1
        
        do while(c1Ind .le. qBack)
            !         find the largest child
            lcInd = c1Ind
            if (qA(c2Ind) .gt. qA(lcInd)) then
                lcInd = c2Ind
            end if
          
            !         compare with child
            !         if necessary swap
            if ( qA(index) .ge. qA(lcInd) ) then
                exit
            end if
          
            !         swap keys
            tmpSwpR = qA(lcInd)
            qA(lcInd) = qA(index)
            qA(index) = tmpSwpR
          
            !         swap indices
            tmpSwpI = qI(lcInd)
            qI(lcInd) = qI(index)
            qI(index) = tmpSwpI
          
            index = lcInd
          
            !         compute child indices
            c1Ind = index*2
            c2Ind = index*2 + 1
          
        end do
          
    end subroutine
      
end module measureWallDistance
