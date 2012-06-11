      module measureWallDistance
      
      implicit none
     
c     AABB tree child   
      type AABB
c       bounding box xyz coordinates      
        real*8, dimension(3) :: bmin,bmax
        
c       maximum radius in this AABB
        real*8 maxr
        
c       number of triangles in this AABB
        integer nTris
        
c       number of children in subtree of this AABB
        integer nChild

c       array of triangle indices        
        integer, dimension(:), allocatable :: tris
        
c       array of child indices        
        integer, dimension(:), allocatable :: child
      end type AABB
      
c     one level of the octree      
      type AABBlevel
        type(AABB), dimension(:), allocatable :: bb
      end type AABBlevel

c     define structure for a single data frame      
      type dataframe
      
        integer :: nPoints,nTris,nEdges
        
c       xyz coordinates of triangle vertices        
        real*8, dimension(:,:), allocatable :: coords
        
c       pseudonormals for points (see Baerentzen et al)        
        real*8, dimension(:,:), allocatable :: pntNrmls
        
c       face normals        
        real*8, dimension(:,:), allocatable :: faceNrmls
        
c       pseudonormals for edges        
        real*8, dimension(:,:,:), allocatable :: edgeNrmls
        
c       array of triangle connectivity        
        integer, dimension(:,:), allocatable :: tris
        
c       triangle bounding spheres center and radius        
        real*8, dimension(:,:), allocatable :: sphc
        real*8, dimension(:), allocatable :: sphr
        
c       arrays to represent queue of AABB or triangles
c       key array
        real*8, dimension(:), allocatable :: qA
c       index, array       
        integer, dimension(:), allocatable :: qI
         
c       array of AABBs to check
        integer, dimension(:), allocatable :: toChBB
        
c       array of triangles to check        
        integer, dimension(:), allocatable :: toChTri
                
c       AABBs        
        type(AABBlevel), dimension(:), allocatable :: octree
                
c       number of AABBs per level of the tree                
        integer, dimension(:), allocatable :: nPerLev
        integer, dimension(:), allocatable :: nPerLevG
        
c       how many levels in octree
        integer nLev       
        
c       size of queues
        integer qSize         
                
      end type dataframe
      
c     data structure for multiple data frames      
      type(dataframe), dimension(:), allocatable, target :: dtfrs
      
c     the total number of data frames      
      integer numDataFrames
      
c     length of the cardiac cycle      
      real*8 cycleLength
      
      contains
      
C /////////////////////////////////////////////////////////////////////////
C  initialization function
C /////////////////////////////////////////////////////////////////////////             

C -------------------------------------------      
C -------------------------------------------
C subroutine dm_initialize
C load data and make octrees 
C -------------------------------------------
C -------------------------------------------     

      subroutine dm_initialize()
      
      integer kk
      
      do kk=1,numDataFrames

        call sphbmesh(kk)
        call makebb(kk,6)

      end do
      
      end subroutine

      
C /////////////////////////////////////////////////////////////////////////
C  functions dealing with one vector in R-3 and a triangular mesh
C /////////////////////////////////////////////////////////////////////////       
      
C -------------------------------------------      
C -------------------------------------------
C subroutine dm_cpmeshp3
C Find a closest point on the triangular mesh
C to the specified closestpnt.  
C Compute the signed distance from 
C the specified point to a closest point on mesh.
C Only called after mesh is loaded.
C Uses octree with bounding spheres
C Only returns one of possibly many closest points
C since we are only interested in distance
C -------------------------------------------
C -------------------------------------------     

        subroutine dm_cpmeshp3(currFr,p,closestpnt,distanceSq,nv)    

        integer currFr
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv 
        real*8 distanceSq   
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      currFr : current frame number
C      p : the reference closestpnt

C OUTPUTS
C      closestpnt : a closest point
C      distanceSq: distance^2 to a closest point on mesh to p
C      nv : normal to the surface at the closest point
C -------------------------------------------           

        integer ii,jj,kk,mm ! counters 
        
        integer qBack           ! back of queue
        integer nCheck,nCheckTri
        
        real*8 ub,lb ! upper bound and lower bound
        real*8 tlb,tub ! temporary lower, upper bound
        
        integer, pointer :: nLev
        
        type(AABB), pointer :: AABBcurr ! pointer to an AABB
        
        nLev=>dtfrs(currFr)%nLev ! pointers for easier syntax
        
c       initialize the checking array of AABB indices
        nCheck = 8
        
        do kk=1,nCheck
          dtfrs(currFr)%toChBB(kk)=kk
        end do
        
c       start at the uppermost level of the tree
c       compute upper and lower bounds based on distances to the AABB
        do kk=2,nLev-1

c         pointer to the first BB          
          AABBcurr=>dtfrs(currFr)%octree(1)%bb(1)
        
c         initialize the bounds
          call fpaabb(p,AABBcurr,ub)
          lb = ub+AABBcurr%maxr
          ub = ub+100
          
c         reset priority queue
          qBack = 0          
           
c         insert the first BB in the priority queue           
          call pqInsert(dtfrs(currFr)%qI,
     &                  dtfrs(currFr)%qA,
     &                  qBack,
     &                  1,lb)          
          
c         loop through all the BBs to check on this level                   
          do jj=1,nCheck
          
c           get index of next BB to check          
            mm = dtfrs(currFr)%toChBB(jj)
            
c           pointer to the BB            
            AABBcurr=>dtfrs(currFr)%octree(kk)%bb(mm)

c           obtain the nearest and farthest possible distance to any point
c           in the AABB
            call cpaabb(p,AABBcurr,tlb)
            call fpaabb(p,AABBcurr,tub)
            
c           add the maximum radius inside the AABB          
            tlb = tlb-AABBcurr%maxr
            tub = tub+AABBcurr%maxr
                        
c           prune the AABBs that are too far away
            do while(dtfrs(currFr)%qA(1) .ge. ub)
              call pqRemoveMax(dtfrs(currFr)%qI,
     &                         dtfrs(currFr)%qA,
     &                         qBack)
            end do

            if (tlb .lt. ub) then   

c             add AABB to the queue  
              call pqInsert(dtfrs(currFr)%qI,
     &                      dtfrs(currFr)%qA,
     &                      qBack,
     &                      mm,tlb)   

c             update the bounds if they become smaller
              if (tlb .lt. lb) then              
                lb = tlb 
              end if
              if (tub .lt. ub) then
                ub = tub 
              end if
              
            end if

          end do
                    
c         copy the BB children indices in the queue to the 
c         array of BBs to check
          nCheck = 0

          do jj=1,qBack
            AABBcurr=>
     &        dtfrs(currFr)%octree(kk)%bb(dtfrs(currFr)%qI(jj))
     
            do ii=1,AABBcurr%nChild
              nCheck = nCheck + 1
              dtfrs(currFr)%toChBB(nCheck)=AABBcurr%child(ii)
            end do
            
          end do
          
        end do
                          
c       copy the triangle indices contained in the leftover BBs
        nCheckTri = 0
        
        do kk=1,nCheck
        
          AABBcurr=>
     &      dtfrs(currFr)%octree(nLev)%bb(dtfrs(currFr)%toChBB(kk))
     
          do jj=1,AABBcurr%nTris
            nCheckTri = nCheckTri + 1
            dtfrs(currFr)%toChTri(nCheckTri) = AABBcurr%tris(jj)
          end do

        end do
        
c        write(*,*) dtfrs(currFr)%toCheck(1:nCheckTri)
      
        call cpmeshpSph(currFr,
     &                  dtfrs(currFr)%toChTri(1:nCheckTri),nCheckTri,
     &                  p,closestpnt,distanceSq,nv)

        end subroutine

C -------------------------------------------      
C -------------------------------------------
C subroutine cpmeshpSph
C Find a closest point on the triangular mesh
C to the specified closestpnt.  
C Compute the signed distance from 
C the specified point to a closest point on mesh.
C Only called after mesh is loaded.
C Uses brute force with bounding spheres
C Only returns one of possibly many closest points
C since we are only interested in distance
C -------------------------------------------
C -------------------------------------------     

        subroutine dm_cpmeshp2(currFr,p,closestpnt,distanceSq,nv) 
        
        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        real*8 distanceSq   

C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      currFr : current frame number
C      p : the reference closestpnt

C OUTPUTS
C      closestpnt : a closest point
C      distanceSq: distance^2 to a closest point on mesh to p
C -------------------------------------------

        integer kk
       
        do kk=1,dtfrs(currFr)%nTris
         dtfrs(currFr)%toChTri(kk)=kk
        end do

        call cpmeshpSph(currFr,
     &                  dtfrs(currFr)%toChTri,
     &                  dtfrs(currFr)%nTris,
     &                  p,closestpnt,distanceSq,nv)
     
        end subroutine  

C -------------------------------------------      
C -------------------------------------------
C subroutine cpmeshpSph
C Find a closest point on the triangular mesh
C to the specified closestpnt.  
C Compute the signed distance from 
C the specified point to a closest point on mesh.
C Only called after mesh is loaded.
C Uses brute force with bounding spheres
C Only returns one of possibly many closest points
C since we are only interested in distance
C -------------------------------------------
C -------------------------------------------   

        subroutine cpmeshpSph(currFr,toCheck,nCheck,
     &                        p,closestpnt,distanceSq,nv) 
        
        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        integer, dimension(:) :: toCheck
        integer nCheck  
        real*8 distanceSq   

C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      currFr : current frame number
C      p : the reference closestpnt
C      toCheck : array of triangles to check
C      nCheck : number of triangles to check

C OUTPUTS
C      closestpnt : a closest point
C      distanceSq: distance^2 to a closest point on mesh to p
C      nv : normal to the surface at the closest point
C -------------------------------------------

        real*8 currDist ! current distance
        real*8, dimension(3) :: currCP,rayToRef
        integer loc
          
        integer ii,jj,kk ! counters
        
        integer qBack ! back of the queue
        
        real*8 ub,lb ! upper bound and lower bound
        real*8 tlb,tub ! temporary lower, upper bound
        
        qBack = 0

c       initially the index of triangle to be checked carefully is
c       going to be the first triangle in the sequence                
        jj=dtfrs(currFr)%toChTri(1)
            
c       get the initial smallest lower and upper bounds  
        ub = sqrt(dot_product(p-dtfrs(currFr)%sphc(:,jj),
     &                        p-dtfrs(currFr)%sphc(:,jj)))      
     &       +dtfrs(currFr)%sphr(jj)
        lb = ub-2*dtfrs(currFr)%sphr(jj)
        
c       store the lower bound for this triangle
c       in the priority queue
        call pqInsert(dtfrs(currFr)%qI,
     &                dtfrs(currFr)%qA,
     &                qBack,
     &                jj,lb)
        
c       loop through the triangles
        do jj=2,nCheck
        
          kk = dtfrs(currFr)%toChTri(jj)

c         compute lower and upper bounds for current triangle
          tub =sqrt(dot_product(p-dtfrs(currFr)%sphc(:,kk),
     &                          p-dtfrs(currFr)%sphc(:,kk)))
     &         +dtfrs(currFr)%sphr(kk)
          
          tlb = tub-2*dtfrs(currFr)%sphr(kk)
          
c         prune triangles from queue that we know are too far away    
          do while(dtfrs(currFr)%qA(1) .ge. ub)
            call pqRemoveMax(dtfrs(currFr)%qI,
     &                       dtfrs(currFr)%qA,
     &                       qBack)
          end do
          
c         if current lower bound is lower than established upper bound         
          if (tlb .lt. ub) then
  
c           add triangle to the queue  
            call pqInsert(dtfrs(currFr)%qI,
     &                    dtfrs(currFr)%qA,qBack,
     &                    kk,tlb)
         
c           update bounds     
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
        
c       loop through the triangles to check
        call cpmeshp(currFr,dtfrs(currFr)%qI(1:qBack),qBack,
     &               p,closestpnt,distanceSq,nv)
        
        end subroutine 
        
C -------------------------------------------        
C -------------------------------------------
C subroutine cpmeshp1
C Find a closest point on the triangular mesh
C to the specified closestpnt.  
C Compute the signed distance from 
C the specified point to a closest point on mesh.
C Only called after mesh is loaded.
C Uses naive brute-force search.
C Only returns one of possibly many closest points
C since we are only interested in distance
C -------------------------------------------
C -------------------------------------------     

        subroutine dm_cpmeshp1(currFr,p,closestpnt,distanceSq,nv) 
        
        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        real*8 distanceSq 
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      currFr : current frame number
C      p : the reference closestpnt

C OUTPUTS
C      closestpnt : a closest point
C      distanceSq: distance^2 to a closest point on mesh to p
C      nv : normal to the surface at the closest point
C -------------------------------------------   

       integer kk
       
       do kk=1,dtfrs(currFr)%nTris
         dtfrs(currFr)%toChTri(kk)=kk
       end do

       call cpmeshp(currFr,
     &              dtfrs(currFr)%toChTri,
     &              dtfrs(currFr)%nTris,
     &              p,closestpnt,distanceSq,nv)
          
       end subroutine
        
C -------------------------------------------        
C -------------------------------------------
C subroutine cpmeshp
C Find a closest point on the triangular mesh
C to the specified closestpnt.  
C Compute the signed distance from 
C the specified point to a closest point on mesh.
C Only called after mesh is loaded.
C Uses naive brute-force search.
C Specify the triangles to check.
C Only returns one of possibly many closest points
C since we are only interested in distance
C -------------------------------------------
C -------------------------------------------     

        subroutine cpmeshp(currFr,toCheck,nCheck,
     &                     p,closestpnt,distanceSq,nv) 
        
        integer currFr 
        real*8, dimension(3) :: p   
        real*8, dimension(3) :: closestpnt,nv
        integer, dimension(:) :: toCheck
        integer nCheck 
        real*8 distanceSq 
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      currFr : current frame number
C      p : the reference closestpnt
C      toCheck : array of triangles to check
C      nCheck : number of triangles to check

C OUTPUTS
C      closestpnt : a closest point
C      distanceSq: distance^2 to a closest point on mesh to p
C      nv : normal to the surface at the closest point
C -------------------------------------------

        real*8 currDist ! current distance
        real*8, dimension(3) :: currCP,rayToRef
        integer loc
          
        integer jj,kk ! counters
          
c       loop through the triangles
        do jj=1,nCheck        
           
          kk = toCheck(jj)
        
          call cptrip(p,
     &           dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(1,kk)),
     &           dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(2,kk)),
     &           dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(3,kk)),
     &           currCP,loc)
     
          rayToRef = p-currCP;
    
c         compute distance
          currDist = sum(rayToRef**2,1)
            
c         determine which normal vector to use
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
            
c         determine the sign of the distance
          currDist = sign(currDist,dot_product(rayToRef,nv))
     
c         update distance and closest point     
          if (abs(currDist) .lt. abs(distanceSq) .or. jj .eq. 1) then
            distanceSq = currDist
            closestpnt = currCP
c            write(*,*) kk
          end if
        end do             
        
        end subroutine
        
C /////////////////////////////////////////////////////////////////////////
C  functions dealing with building spatial heirarchy
C /////////////////////////////////////////////////////////////////////////          
        
C -------------------------------------------
C -------------------------------------------
C subroutine makebb
C Construct the AABB octree for the mesh
C Assumes that triangle bounding volumes
C have been constructed already
C Assume at least one level of refinement
C -------------------------------------------
C -------------------------------------------

        subroutine makebb(currFr,refine)
        
        integer currFr,refine
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      currFr : current frame number
C      refine : desired level of refinement
C -------------------------------------------        
        
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
                 
c       allocate enough space                 
        allocate(dtfrs(currFr)%nPerLev(1+refine))
        allocate(dtfrs(currFr)%nPerLevG(1+refine))
        allocate(dtfrs(currFr)%octree(1+refine))
        
c       set the box presence grid to all zeros        
        allocate(boxHere(2**refine,2**refine,2**refine))
        boxHere = 0
        
c       set initial array sizes at each level of the tree
        do kk=1,refine+1
        
          dtfrs(currFr)%nPerLev(kk) = 0
          dtfrs(currFr)%nPerLevG(kk) = 8
          
        end do
        
c       compute initial bounding box
c       look for the maximum extent of x,y,z coordinates of the
c       triangle bounding spheres
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

c        write(*,*) "Bounds"        
c        write(*,*) pmin
c        write(*,*) pmax
c        write(*,*) "Max radius"
c        write(*,*) maxr

c       these are the dimensions of the whole bounding volume        
        boundDim = pmax-pmin
                
c       these are the dimensions of the smallest AABB        
        boxDim = boundDim / 2**refine
        
c        write(*,*) "Box dim"
c        write(*,*) boxDim

c       pointers for easier access
        nPerLev => dtfrs(currFr)%nPerLev(1+refine)
        nPerLevG => dtfrs(currFr)%nPerLevG(1+refine)
        
c       allocate enough space for each level        
        allocate(tempLev(nPerLevG))
        allocate(dtfrs(currFr)%octree(1+refine)%bb(nPerLevG))
        
c       define the number of levels in the tree to be one greater
c       than the refine parameter         
        dtfrs(currFr)%nLev = refine+1        
        
c       we loop through each triangle and assign them to the correct
c       AABB at the bottommost level        
        do kk=1,dtfrs(currFr)%nTris
        
c         get the "xyz index" of the sphere
c         these are just integer values which will
c         help us define the bounds for the AABB        
          xyzInd = ceiling((dtfrs(currFr)%sphc(:,kk)-pmin)/boxdim)
          
          do ii=1,3  
          
            if (xyzInd(ii) .eq. 0) then
              xyzInd(ii) = 1
            end if
            
            if (xyzInd(ii) .eq. 2**refine+1) then
              xyzInd(ii) = 2**refine
            end if
            
          end do
          
c         if no previous AABB has been made in this spatial location
c         then we make one
          if (boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) .eq. 0) then
          
c           compute the center of the box          
            tcent = xyzInd*boxDim-0.5*boxdim+pmin
c           compute the bounds of the box            
            tmin = tcent-0.5*boxdim
            tmax = tcent+0.5*boxdim
                    
c           increment the number of AABBs at his level
            nPerLev = nPerLev + 1
            
c           we store the number of the new AABB            
c           by doing this, we make sure that no other AABBs can
c           be created at this spatial location            
            boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) = nPerLev
              
c            write(*,*) xyzInd  
              
c           resize the AABB array at this level
c           if necessary              
            if (nPerLev .gt. nPerLevG) then       
              
              nPerLevG = nPerLevG + 8**(refine-1)
              
c              write(*,*) nPerLevG
              
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
                        
c           pointer to the new box for easier access
            AABBcurr=>dtfrs(currFr)%octree(1+refine)%bb(nPerLev)
            
c           assign bounds            
            AABBcurr%bmin = tmin
            AABBcurr%bmax = tmax
            
c           this AABB has a triangle count of 1            
            AABBcurr%nTris = 1

c           since we are at the bottommost level,
c           this AABB has no children  
            AABBcurr%nChild = 0
            
c           allocate space for the triangle array            
            allocate(AABBcurr%tris(1))
            
c           assign the triangle index            
            AABBcurr%tris(1) = kk
            
c           the maximum radius in this AABB            
            AABBcurr%maxr = dtfrs(currFr)%sphr(kk)
 
          else
            
c           there is already an AABB at this spatial location
c           we merely have to find its index            
            AABBcurr=>dtfrs(currFr)%octree(1+refine)%bb(
     &        boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) )
     
c           increment triangle count     
            AABBcurr%nTris = AABBcurr%nTris+1
     
c           resize the triangle array     
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
            
c           assign the triangle index            
            AABBcurr%tris(AABBcurr%nTris) = kk
            
c           if the sphere bounding the triangle 
c           has a larger radius than the current
c           maximum radius in the AABB, then it becomes the
c           new max radius in the AABB            
            if (dtfrs(currFr)%sphr(kk) .gt. AABBcurr%maxr) then
              AABBcurr%maxr = dtfrs(currFr)%sphr(kk)
            end if
      
          end if
                   
        end do
        
        !do ii=1,dtfrs(currFr)%nPerLev(3)
        !  write(*,*) dtfrs(currFr)%octree(3)%bb(ii)%nTris
        !end do
        
c       deallocate arrays        
        if (allocated(tempLev)) then
          deallocate(tempLev)
        end if
        
        if (allocated(boxHere)) then
          deallocate(boxHere)
        end if
        
c       we have now created the AABBs that bound the triangle volumes
c       at the bottommost level of the tree   
c       now we go up the levels
        do kk=refine,2,-1
        
c       reallocate the box presence grid       
        allocate(boxHere(2**(kk-1),2**(kk-1),2**(kk-1)))
        
c       initialize the box presence grid to zero        
        boxHere = 0
        
c       compute the dimensions of the AABBs in this level        
        boxDim = boundDim / 2**(kk-1)
        
c       pointers for easier access        
        nPerLev=>dtfrs(currFr)%nPerLev(kk)
        nPerLevG=>dtfrs(currFr)%nPerLevG(kk)
        
c       allocate enough space for the level        
        allocate(tempLev(nPerLevG))
        allocate(dtfrs(currFr)%octree(kk)%bb(nPerLevG))
        
c       loop through the AABBs at this level      
        do jj=1,dtfrs(currFr)%nPerLev(1+kk)

c         get the "xyz index" of the AABB
c         just like we did for the triangle volumes
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
          
c         if no AABB exists here, then create one          
          if (boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) .eq. 0) then
            
c           compute center of AABB            
            tcent = xyzInd*boxDim-0.5*boxdim+pmin
            
c           compute bounds of the AABB            
            tmin = tcent-0.5*boxdim
            tmax = tcent+0.5*boxdim
                    
c           increment the AABB array at this level by one          
            nPerLev = nPerLev + 1
            
c           store AABB index            
            boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) = nPerLev
            
c           resize AABB array at this level            
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
            
c           pointer to the new AABB
            AABBnew=>dtfrs(currFr)%octree(kk)%bb(nPerLev)
            
c           assign bounds            
            AABBnew%bmin = tmin
            AABBnew%bmax = tmax
            
c           assign number of triangles
c           inherit this value from the AABB encompassed 
c           in the lower level            
            AABBnew%nTris = AABBcurr%nTris            

c           same for the max radius         
            AABBnew%maxr = AABBcurr%maxr    
            
c           set number of children to 1            
            AABBnew%nChild = 1
            
c           allocate space for the child array            
            allocate(AABBnew%child(1))
            
c           assign index of the AABB in the lower level            
            AABBnew%child(1) = jj
            
          else
          
c           get the pointer to the box already here          
            AABBnew=>dtfrs(currFr)%octree(kk)%bb(
     &        boxHere(xyzInd(1),xyzInd(2),xyzInd(3)) ) 
            
c           increase triangle count 
c           add the triangle count of the AABB in the lower level            
c           to the triangle count of the current AABB
            AABBnew%nTris = AABBnew%nTris+AABBcurr%nTris
            
c           update maximum radius            
            if (AABBcurr%maxr .gt. AABBnew%maxr) then
              AABBnew%maxr = AABBcurr%maxr
            end if
            
c           increment number of children            
            AABBnew%nChild = AABBnew%nChild + 1
            
c           resize child array            
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
            
c           assign index of child            
            AABBnew%child(AABBnew%nChild) = jj
            
          end if    
          
        end do
        
c       deallocate arrays             
        if (allocated(tempLev)) then
          deallocate(tempLev)
        end if
        
        if (allocated(boxHere)) then
          deallocate(boxHere)
        end if
        
        
        end do
        
c       write the information 
c       for the level 1 AABB which just encompasses everything   
c       this is for mostly debugging
c       since we can check if the max radius and number of triangles
c       was correctly propagated up the tree
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
        
c       temporary the toCheck array AABBs
        allocate(dtfrs(currFr)%toChBB(
     &           dtfrs(currFr)%nPerLev(
     &           dtfrs(currFr)%nLev)))
     
c       find the proper size for allocation of the queues
        dtfrs(currFr)%qSize = dtfrs(currFr)%nTris
        
        if ( dtfrs(currFr)%nPerLev(dtfrs(currFr)%nLev) .gt. 
     &    dtfrs(currFr)%qSize ) then
     
          dtfrs(currFr)%qSize = 
     &      dtfrs(currFr)%nPerLev(dtfrs(currFr)%nLev)
     
        end if
        
c       allocate space for priority queue
        allocate(dtfrs(currFr)%qA(dtfrs(currFr)%qSize))
        allocate(dtfrs(currFr)%qI(dtfrs(currFr)%qSize))
                
c       for debugging purposes                
        
c        write(*,*) AABBnew%nTris
c        write(*,*) AABBnew%maxr
        
c        do kk=1,refine+1
c        write(*,*) "Level",kk
c        do jj=1,dtfrs(currFr)%nPerLev(kk)
c          AABBcurr=>dtfrs(currFr)%octree(kk)%bb(jj)  
c          write(*,*) jj,AABBcurr%nChild,AABBcurr%nTris
c        end do
c        end do

        end subroutine
        
C -------------------------------------------        
C -------------------------------------------
C subroutine sphbmesh
C Compute bounding spheres for triangular mesh
C -------------------------------------------
C -------------------------------------------      

        subroutine sphbmesh(currFr)
        
        integer currFr
        integer kk
        
        allocate( dtfrs(currFr)%sphc
     &            (3,dtfrs(currFr)%nTris) )
        allocate( dtfrs(currFr)%sphr
     &            (dtfrs(currFr)%nTris) )
     
c       loop through the triangles
        do kk=1,dtfrs(currFr)%nTris
          call sphbt1(dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(1,kk)),
     &                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(2,kk)),
     &                dtfrs(currFr)%coords(:,dtfrs(currFr)%tris(3,kk)),
     &                dtfrs(currFr)%sphc(:,kk),dtfrs(currFr)%sphr(kk) )
        end do
        
        end subroutine
      	  
C /////////////////////////////////////////////////////////////////////////
C  functions dealing with AABB
C /////////////////////////////////////////////////////////////////////////  	  
	  
C -------------------------------------------
C -------------------------------------------
C subroutine cpaabb
C Compute distance to AABB
C -------------------------------------------
C -------------------------------------------      	  

        subroutine cpaabb(p,box,distance)
        
        real*8, dimension(3) :: p
        type(AABB) box
        real*8 distance,v
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      p : point p
C      box : AABB

C OUTPUTS
C      distanceSq  : distance squared
C -------------------------------------------          
        
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
        
C -------------------------------------------
C -------------------------------------------
C subroutine fpaabb
C Compute farthest distance to AABB
C -------------------------------------------
C -------------------------------------------      	  

        subroutine fpaabb(p,box,distance)
        
        real*8, dimension(3) :: p,fp,cent
        type(AABB) box
        real*8 distance
        integer, dimension(3) :: side
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      p : point p
C      box : AABB

C OUTPUTS
C      distanceSq  : distance squared
C ------------------------------------------- 

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
        
C /////////////////////////////////////////////////////////////////////////
C  functions dealing with triangles
C /////////////////////////////////////////////////////////////////////////  

C -------------------------------------------      
C -------------------------------------------
C subroutine cptrip
C Find the closest point on a triangle to
C specified point   
C -------------------------------------------
C -------------------------------------------
        subroutine cptrip(p,a,b,c,closestpnt,loc)
       
        real*8, dimension(3) :: a 
        real*8, dimension(3) :: b 
        real*8, dimension(3) :: c 
        real*8, dimension(3) :: p 
        real*8, dimension(3) :: closestpnt 
        integer loc 

C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE

C      p : the reference closestpnt
C      a : point a of triangle
C      b : point b of triangle
C      c : point c of triangle

C OUTPUTS
C      point: closest point on triangle to p
C      loc: 1 denotes closest point is on triangle vertex a
C            2 denotes closest point is on triangle vertex b
C            3 denotes closest point is on triangle vertex c
C            4 denotes closest point is on loc ab
C            5 denotes closest point is on loc ac
C            6 denotes closest point is on loc bc
C            7 denotes closest point is on the face
C -------------------------------------------
         
        real*8, dimension(3) :: ab,ac,ap,bp,cp
         
        real*8 d1,d2,d3,d4,d5,d6,v,w,va,vb,vc,denom
         
C       Check if p in vertex region outside A	   
	  ab = b - a
	  ac = c - a
	  ap = p - a
	  d1 = dot_product(ab,ap)
	  d2 = dot_product(ac,ap)	   
	  if (d1 .le. 0.0 .and. d2 .le. 0) then
C         Barycentric coordinates (1,0,0)	   
	    closestpnt = a 
	    loc = 1
	    return
	  end if
	   
C       Check if p in vertex region outside B
        bp = p - b
        d3 = dot_product(ab,bp)   
        d4 = dot_product(ac,bp)
        if (d3 .ge. 0.0 .and. d4 .le. d3) then
C         Barycentric coordinates (0,1,0)	            
          closestpnt = b
          loc = 2
          return
        end if
         
C       Check if p in loc region of AB, if so return projection of P onto AB
        vc = d1*d4 - d3*D2
        if (vc .le. 0.0 .and. d1 .ge. 0.0 .and. d3 .le. 0.0) then
          v = d1 / (d1 - d3)
          closestpnt = a + v*ab
          loc = 4
          return
        end if
         
C       Check if p in vertex region outside C
        cp = p - c
        d5 = dot_product(ab,cp)
        d6 = dot_product(ac,cp)
        if (d6 .ge. 0.0 .and. d5 .le. d6) then
C          Barycentric coordinates (0,0,1)           
          closestpnt = c
          loc = 3
          return
        end if
         
C       Check if p in loc region of AC, if so return projection of P onto AC
        vb = d5*d2 - d1*d6
        if (vb .le. 0.0 .and. d2 .ge. 0.0 .and. d6 .le. 0.0) then
          w = d2 / (d2 - d6)
          closestpnt = a + w*ac
          loc = 5
          return
        end if
         
C       Check if p in loc region of BC, if so return projection of P onto BC         
        va = d3*d6 - d5*d4
        if (va .le. 0.0 .and. (d4 - d3) .ge. 0.0 .and. 
     &      (d5 - d6) .ge. 0.0) then
          w = (d4 - d3) / ( (d4 - d3) + (d5 - d6) )
          closestpnt = b + w*(c - b)
          loc = 6
          return
        end if
         
C       p inside face region. compute closestpnt using barycentric coordinates
        denom = 1.0 / (va+vb+vc)
        v = vb*denom
        w = vc*denom
        closestpnt = a + ab*v + ac*w
        loc = 7
         
	  end subroutine  
                 
C -------------------------------------------
C -------------------------------------------
C subroutine sphbt2
C Compute bounding sphere for triangle
C Version 2
C -------------------------------------------
C -------------------------------------------      

        subroutine sphbt2(a,b,c,gcent,radius)
        
        real*8, dimension(3) :: a,b,c,gcent
        real*8 radius,tr
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      a : point a
C      b : point b
C      c : point c
C OUTPUTS
C      gcent  : center of the bounding sphere
C      radius : radius of the bounding sphere        
	  
C -------------------------------------------

c      compute the centroid of the triangle
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
       
C -------------------------------------------
C subroutine sphbt1
C Compute bounding sphere for triangle
C -------------------------------------------
C -------------------------------------------      

        subroutine sphbt1(a,b,c,gcent,radius)
        
        real*8, dimension(3) :: a,b,c,gcent
        real*8 radius
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      a : point a
C      b : point b
C      c : point c
C OUTPUTS
C      gcent  : center of the bounding sphere
C      radius : radius of the bounding sphere
C -------------------------------------------           
        
        real*8, dimension(3) :: al,bl,cl,cent,tr
        real*8, dimension(3,3) :: T
        real*8 D,la,lb,lc
        
        la = sqrt(dot_product(b-a,b-a))
        lb = sqrt(dot_product(c-a,c-a))
        lc = sqrt(dot_product(b-c,b-c))
        
        radius = la*lb*lc / 
     &           sqrt((la+lb+lc)*(-la+lb+lc)*(la-lb+lc)*(la+lb-lc))
        
        call rottriloc(a,b,c,al,bl,cl,T)
        
c       ignore the z coordinate
c       translate local point a to the origin
        tr = al
        al(1:2) = 0
        bl(1:2) = bl(1:2) - tr(1:2)
        cl(1:2) = cl(1:2) - tr(1:2)
        
c       compute center of circumscribed circle
        D = 2*(bl(1)*cl(2)-bl(2)*cl(1))
        cent(1) = ( cl(2)*(bl(1)**2+bl(2)**2)-
     &              bl(2)*(cl(1)**2+cl(2)**2) ) / D
        cent(2) = ( bl(1)*(cl(1)**2+cl(2)**2)-
     &              cl(1)*(bl(1)**2+bl(2)**2) ) / D
        cent(3) = al(3)
        
c       translate back
        cent(1:2) = cent(1:2) + tr(1:2)
              
c       rotate back to global coordinate system
        gcent(1) = sum(T(:,1)*cent,1)
        gcent(2) = sum(T(:,2)*cent,1)
        gcent(3) = sum(T(:,3)*cent,1)    
        
        end subroutine

C -------------------------------------------
C -------------------------------------------
C subroutine rottriloc
C Rotate triangle to local 2D cartesian system
C -------------------------------------------
C -------------------------------------------           
      
        subroutine rottriloc(a,b,c,al,bl,cl,T)
        
        real*8, dimension(3) :: a,b,c,al,bl,cl
        real*8, dimension(3,3) :: T
        
C -------------------------------------------	   
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C      a : point a
C      b : point b
C      c : point c
C OUTPUTS
C      al : local point al
C      bl : local point bl
C      cl : local point cl
C      T  : 3x3 transformation matrix
C -------------------------------------------

        real*8, dimension(3) :: v1,v2,v3

c       compute local basis vectors
        v1 = b - a
        v1 = v1/sqrt(dot_product(v1,v1))
        v2 = c - a
        v3 = cross_product(v1, v2)
        v3 = v3/sqrt(dot_product(v3,v3))
        v2 = cross_product(v3, v1)
        
c       compute transformation matrix       
        T(1,:) = v1
        T(2,:) = v2
        T(3,:) = v3
        
c       rotate to local cartesian coordinates        
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
	  
C -------------------------------------------
C -------------------------------------------
C subroutine dataread
C Read in the mesh data 
C Call this first before cpmeshp
C -------------------------------------------
C -------------------------------------------

        subroutine dm_dataread(filename)

        character(*) :: filename
        
C -------------------------------------------
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE

C      filename : the name of the file to be read

C OUTPUTS
C      no outputs
C -------------------------------------------        

        character dummy
        integer dummyint
        
        integer currFr ! current data frame
        
        integer :: ii,jj,kk ! counters
        
        open(unit=2,file=filename)

c       Skip the first identifier line and 
c       read the number of data frames
        read(unit=2,fmt='(A1)') dummy
        read(2,*) numDataFrames
        
c       Skip the first identifier line and 
c       read the length of the cycle
        read(unit=2,fmt='(A1)') dummy
        read(2,*) cycleLength
        
c       Allocate space for the data frames
        allocate(dtfrs(numDataFrames))     
        
        do kk=1,numDataFrames
c         Skip the next identifier line and
c         read the data frame number        
          read(unit=2,fmt='(A1)') dummy
          read(2,*) currFr
        
c         Skip the next identifier line and
c         read the number of points        
          read(unit=2,fmt='(A1)') dummy
          read(2,*) dtfrs(currFr)%nPoints

c         Skip the next identifier line and
c         read the number of triangles 
          read(unit=2,fmt='(A1)') dummy
          read(2,*) dtfrs(currFr)%nTris

c         Skip the next identifier line and
c         read the number of edges
          read(unit=2,fmt='(A1)') dummy
          read(2,*) dtfrs(currFr)%nEdges

c         Now allocate memory 
          allocate( dtfrs(currFr)%coords
     &            (3,dtfrs(currFr)%nPoints) )
     
          allocate( dtfrs(currFr)%tris
     &              (3,dtfrs(currFr)%nTris) )
     
          allocate( dtfrs(currFr)%pntNrmls
     &              (3,dtfrs(currFr)%nPoints) )
     
          allocate( dtfrs(currFr)%edgeNrmls
     &              (3,3,dtfrs(currFr)%nTris) ) ! index by dof, edge, then tri
        
          allocate( dtfrs(currFr)%faceNrmls
     &              (3,dtfrs(currFr)%nTris) )
             
c         Skip the next identifier line and
c         read the coordinates        
          read(unit=2,fmt='(A1)') dummy
          do ii=1,dtfrs(currFr)%nPoints
            read(2,*) dtfrs(currFr)%coords(1,ii),
     &                dtfrs(currFr)%coords(2,ii),
     &                dtfrs(currFr)%coords(3,ii)
          end do
        
c         Skip the next identifier line and
c         read the surface triangle connectivity
          read(unit=2,fmt='(A1)') dummy
          do ii=1,dtfrs(currFr)%nTris
            read(2,*) dtfrs(currFr)%tris(1,ii),
     &                dtfrs(currFr)%tris(2,ii),
     &                dtfrs(currFr)%tris(3,ii)
          end do        
        
c         Skip the next identifier line and
c         read the face normals      
          read(unit=2,fmt='(A1)') dummy
          do ii=1,dtfrs(currFr)%nTris
            read(2,*) dtfrs(currFr)%faceNrmls(1,ii),
     &                dtfrs(currFr)%faceNrmls(2,ii),
     &                dtfrs(currFr)%faceNrmls(3,ii)
          end do        
        
c         Skip the next identifier line and
c         read the edge normals    
          read(unit=2,fmt='(A1)') dummy
          do ii=1,dtfrs(currFr)%nTris
            do jj=1,3
              read(2,*) dummyint,
     &                  dtfrs(currFr)%edgeNrmls(1,jj,ii),
     &                  dtfrs(currFr)%edgeNrmls(2,jj,ii),
     &                  dtfrs(currFr)%edgeNrmls(3,jj,ii)
            end do
          end do               
        
c         Skip the next identifier line and
c         read the closestpnt normals   
          read(unit=2,fmt='(A1)') dummy
          do ii=1,dtfrs(currFr)%nPoints
            read(2,*) dtfrs(currFr)%pntNrmls(1,ii),
     &                dtfrs(currFr)%pntNrmls(2,ii),
     &                dtfrs(currFr)%pntNrmls(3,ii)
          end do               
          
c         allocate enough space for toCheck
c         for safety, we will allocate the with the number of triangles  
          
c         array of triangles to check          
          allocate(dtfrs(currFr)%toChTri(dtfrs(currFr)%nTris))
          
        end do

c        write(*,*) dummy
c        write(*,*) numDataFrames
c        write(*,*) dtfrs(2)%nPoints
c        write(*,*) dtfrs(2)%nTris
c        write(*,*) dtfrs(2)%nEdges
        
c        do ii=1,dtfrs(currFr)%nTris
c          write(*,*) dtfrs(currFr)%tris(1,ii),
c     &               dtfrs(currFr)%tris(2,ii),
c     &               dtfrs(currFr)%tris(3,ii)
c        end do
      
        end subroutine
        
C /////////////////////////////////////////////////////////////////////////
C  functions dealing with vectors in R-3
C /////////////////////////////////////////////////////////////////////////        

C -------------------------------------------        
C -------------------------------------------
C function cross_product
C computes cross product of two vectors
C -------------------------------------------        
C -------------------------------------------

        function cross_product(a, b) 
        
        real*8, dimension(3) :: a, b, cross_product

C -------------------------------------------
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE

C      a : 3x1
C      b : 3x1

C OUTPUTS
C      cross_product
C -------------------------------------------          

        cross_product(1) = a(2) * b(3) - a(3) * b(2)
        cross_product(2) = a(3) * b(1) - a(1) * b(3)
        cross_product(3) = a(1) * b(2) - a(2) * b(1)
            
        end function cross_product
        
C /////////////////////////////////////////////////////////////////////////
C  Functions that deal with priority queue insertion, removal
C /////////////////////////////////////////////////////////////////////////
        
C -------------------------------------------        
C -------------------------------------------
C subroutine pqInsert
C insert a new item in the priority queue
C -------------------------------------------        
C -------------------------------------------                

        subroutine pqInsert(qI,qA,qBack,ival,pval)
        
        integer, dimension(:) :: qI
        real*8, dimension(:) :: qA
        integer qBack,ival
        real*8 pval
        
C -------------------------------------------
C INPUT PARAMETERS
C NONE OF THESE ARE CHANGED BY THE SUBROUTINE
C
C      ival : index value
C      pval : key value
C
C THESE ARE CHANGED
C      
C      qI : queue of indices
C      qA : queue of keys
C      qBack : index of the back of the queue
C
C -------------------------------------------  

        integer index,parInd    
        real*8 tmpSwpR
        integer tmpSwpI    
        
c       add to last level        

c       add key
        qBack = qBack + 1
        qA(qBack) = pval
        
c       add index value
        qI(qBack) = ival        
        
        index = qBack

        do while (index .gt. 1)
        
c         compare with parent
c         if necessary swap       
          parInd = index/2   
          
          if (qA(index) .le. qA(parInd)) then  
            exit
          end if
          
c         swap values of keys          
          tmpSwpR = qA(parInd)
          qA(parInd) = qA(index)
          qA(index) = tmpSwpR
          
c         swap values of indices          
          tmpSwpI = qI(parInd)
          qI(parInd) = qI(index)
          qI(index) = tmpSwpI
          
          index = parInd
          
        end do
                 
        end subroutine
        
C -------------------------------------------        
C -------------------------------------------
C subroutine pqRemoveMax
C removes the index in the priority queue associated
C with the largest key
C -------------------------------------------        
C -------------------------------------------            
      
        subroutine pqRemoveMax(qI,qA,qBack)
          
        integer, dimension(:) :: qI  
        real*8, dimension(:) :: qA
        integer qBack
        
C -------------------------------------------
C INPUT PARAMETERS
C
C THESE ARE CHANGED
C      
C      qI : queue of indices
C      qA : queue of keys
C      qBack : index of the back of the queue
C
C -------------------------------------------    

        integer c1Ind,c2Ind,lcInd,index      
        real*8 tmpSwpR
        integer tmpSwpI
        
c       do not remove max if only one element left        
        if (qBack .gt. 1) then
c         assign the value of the last child to the root        
          qA(1) = qA(qBack)
          qI(1) = qI(qBack)
c         decrease size of queue    
          qBack = qBack - 1
        end if
        
        index = 1
        c1Ind = index*2
        c2Ind = index*2 + 1
        
        do while(c1Ind .le. qBack)
c         find the largest child        
          lcInd = c1Ind
          if (qA(c2Ind) .gt. qA(lcInd)) then
            lcInd = c2Ind
          end if
          
c         compare with child
c         if necessary swap          
          if ( qA(index) .ge. qA(lcInd) ) then
            exit
          end if
          
c         swap keys          
          tmpSwpR = qA(lcInd)
          qA(lcInd) = qA(index)
          qA(index) = tmpSwpR
          
c         swap indices
          tmpSwpI = qI(lcInd)
          qI(lcInd) = qI(index)
          qI(index) = tmpSwpI                    
          
          index = lcInd
          
c         compute child indices
          c1Ind = index*2
          c2Ind = index*2 + 1
          
        end do
          
        end subroutine
      
      end module measureWallDistance