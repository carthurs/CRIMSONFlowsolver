      subroutine genBC (iBC,  BC,   x,   ilwork, iper)
!
!----------------------------------------------------------------------
!  This routine generates the essential prescribed boundary conditions.
!
! input:
!  iBC   (nshg)        : boundary condition code
!  nBC   (nshg)        : boundary condition mapping array
!
! output:
!  BC    (nshg,ndofBC) : The constraint data for prescribed BC 
!
!
! Note: genBC1 reduces the input data for the velocity. In the
!       case of varying velocity direction in the generation, the 
!       results may not be correct. (since a linearity assumption is 
!       made in the generation).
!
!
! Farzin Shakib, Spring 1986.
! Zdenek Johan,  Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
      use readarrays            ! used to access BCinp, nBC
      use specialBC ! filling acs here
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension iBC(nshg),                nsurf(nshg), &
                BC(nshg,ndofBC), &
                x(numnp,nsd),             ilwork(nlwork), &
                iper(nshg)
!
! BCinp for each point has:
!   D T P c11 c12 c13 M1 c21 c22 c23 M2 theta S1 S2 S3...
!   1 2 3 4   5   6   7  8   9   10  11 12    13 14 15...
! Remember, ndof=nsd+2+nsclr
!
! Arrays in the following 1 line are now dimensioned in readnblk
!        dimension BCinp(numpbc,ndof+7)
!  
      dimension BCtmp(nshg,ndof+7)
!
! ndof+7= 3(thermos) + (nsd-1)*(nsd+1) + nscalars + 1 (theta)
!                       #vect *(vec dir +mag)
!
!.... --------------------------->  Input  <---------------------------
!
!.... convert boundary condition data
!
      BCtmp = zero
!
      if(numpbc.ne.0) then  
         do i = 1, ndof+7
            where (nBC(:) .ne. 0) BCtmp(:,i) = BCinp(nBC(:),i)
         enddo
         deallocate(BCinp)
      endif
            
!
      if(any(BCtmp(:,12).ne.0)) then
         iabc=1
         allocate (acs(nshg,2))
         where (btest(iBC,10))
            acs(:,1) = cos(BCtmp(:,12)) 
            acs(:,2) = sin(BCtmp(:,12)) 
         endwhere
      endif
           
!
!.... ----------------------> Wall Normals  <--------------------------
! (calculate the normal and adjust BCinp to the true normal as needed)
!
!
      call genwnm (iBC,  BCtmp,   x,   ilwork, iper, nsurf)
!
!  determine the first point off the wall for each wall node
!
      call genotwn (x,BCtmp, iBC, nsurf)
!
!.... ------------------------>  Conversion  <-------------------------
!
!.... convert the input boundary conditions to condensed version
!
      BC = zero
!
      call genBC1 (BCtmp,  iBC,  BC)
!
!.... --------------------------->  Echo  <----------------------------
!
!.... echo the input data
!
      if (necho .lt. 3) then
         nn = 0
         do n = 1, nshg
            if (nBC(n) .ne. 0) then
               nn = nn + 1
               if(mod(nn,50).eq.1)  &
                    write(iecho,1000)ititle,(j,j=1,ndofBC)
               write (iecho,1100) n, (BC(n,i),i=1,ndofBC)
            endif
         enddo
      endif
!     
!.... return
!
      return
!
 1000 format(a80,//, &
      ' P r e s c r i b e d   B o u n d a r y   C o n d i t i o n s',//, &
        '    Node  ',/, &
        '   Number ',5x,6('BC',i1,:,10x))
 1100 format(1p,2x,i5,3x,6(e12.5,1x))
!
      end










      subroutine genwnm (iBC, BCtmp,    x,   ilwork, iper, nsurf)
!----------------------------------------------------------------------
!  This routine generates the normal to a wall
!
! input:
!  iBC   (nshg)        : boundary condition code
!  nBC   (nshg)        : boundary condition mapping array
!
! output:
!  BCtmp    (nshg,ndof+6) : The constraint data for prescribed BC 
!
!
!----------------------------------------------------------------------
!
      use turbSA
      use pointer_data          ! used for mienb, mibcb
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      character*20 fname1,  fmt1
      character*5  cname         
      dimension iBC(nshg),                  iper(nshg), &
                  x(numnp,nsd),             ilwork(nlwork)
!
      dimension  BCtmpSAV(nshg,ndof+7)
      dimension  BCtmp(nshg,ndof+7),      fBC(nshg,ndofBC), &
                  e1(3),                    e2(3), &
                  elnrm(3),                 asum(numnp)
!
      integer sid, nsidg
      integer nsurf(nshg), ivec(nsd)
      logical :: firstvisit(nshg)
      real*8 BCvecs(2,nsd)
      integer, allocatable :: ienb(:)
      dimension wnmdb(nshg,nsd)
!
!  wnrm is dimensioned nshg but the math is only done for straight sided
!  elements at this point so wnrm will not be calculated for the hierarchic 
!  modes.  Note that the wall model creates a p.w. linear representation
!  only at this time.
!
      allocate ( wnrm(nshg,3) )
!
!.... ----------------------> Wall Normals  <--------------------------
! (calculate the normal and adjust BCinp to the true normal as needed)
!
!
      asum = zero
      wnrm = zero

!
!....  Save a copy of BCtmp so that after we calculate the normals we
!      can recover the comp3 information.
!
      BCtmpSAV=BCtmp
!
! Count out the number of surface ID's on-processor, and map them
      call gensidcount(nsidg)

! 
      if(nsidg.gt.0) then       ! if there are any surfID's
         nsurf(:) = 0
         do k = 1, nsidg        ! loop over Surface ID's
            sid = sidmapg(k)
            firstvisit(:)=.true.
            wnrm(:,1:3)=zero
            do iblk=1, nelblb   ! loop over boundary element blocks
               npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
               nenbl = lcblkb(6,iblk)
               nshl = lcblkb(9,iblk)
               allocate( ienb(nshl) )
               do i = 1, npro   ! loop over boundary elements
                  iBCB1=miBCB(iblk)%p(i,1)
                  iBCB2=miBCB(iblk)%p(i,2)
                  ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
! don't bother with elements that aren't on modeled surfaces
!              if ( not          (wall set  and   traction set)    )
                  if (.not.(btest(iBCB1,2).and.btest(iBCB1,4))) &
                       cycle
! don't bother with elements that don't lie on the current surface
                  if (iBCB2.ne.sid) cycle
!     
!.... calculate this element's area-weighted normal vector
!     
                  e1 = x(ienb(2),:)-x(ienb(1),:)
                  e2 = x(ienb(3),:)-x(ienb(1),:)
                  elnrm(1) = e1(2)*e2(3)-e1(3)*e2(2)
                  elnrm(2) = e1(3)*e2(1)-e1(1)*e2(3)
                  elnrm(3) = e1(1)*e2(2)-e1(2)*e2(1)
! Tetrahedral elements have negative volumes in phastaI, so
! the normals calculated from the boundary faces must be inverted
! to point into the fluid
                  if(nenbl.eq.3) elnrm(:)=-elnrm(:)
!     
!.... add area-weighted normals to the nodal tallies for this surface
!     
                  do j = 1, nenbl ! loop over elt boundary nodes
                     nn=ienb(j) ! global node number
                     if(firstvisit(nn)) then
                        firstvisit(nn)=.false.
                        nsurf(nn)=nsurf(nn)+1
                        if(nsurf(nn).eq.1) BCtmp(nn,4:6)=zero
                        if(nsurf(nn).eq.2) BCtmp(nn,8:10)=zero
                     endif
                     wnrm(nn,:)=wnrm(nn,:)+elnrm
                  enddo         ! loop over elt boundary nodes
               enddo            ! end loop over boundary elements in block
               deallocate(ienb)
            enddo               ! end loop over boundary element blocks
! Now we have all of this surface's contributions to wall normals
! for all nodes, along with an indication of how many surfaces
! each node has encountered so far.  Contributions from other processors
! should now be accumulated for this surface
!     
! axisymm BC's need BC (and we have not built it yet) so we need to create
! the entries it needs.
!
            if(iabc==1) then
               where (btest(iBC,10))
                  fBC(:,1) = cos(BCtmp(:,12)) 
                  fBC(:,2) = sin(BCtmp(:,12))
               endwhere

! before the commu we need to rotate the residual vector for axisymmetric
! boundary conditions (so that off processor periodicity is a dof add instead
! of a dof combination).  Take care of all nodes now so periodicity, like
! commu is a simple dof add.
               call rotabc(wnrm, iBC, 'in ')
            endif

            if (numpe.gt.1) then
!.... add areas and norms contributed from other processors     
               call commu (wnrm, ilwork, 3, 'in ')
            endif
!.... account for periodicity and axisymmetry
            call bc3per(iBC,wnrm, iper,ilwork, 3)
! at this point the master has all the information (slaves are zeroed and
! waiting to be told what the master has...lets do that).
!
!.... local periodic (and axisymmetric) boundary conditions (no communications)
            wnmdb=wnrm
            do i = 1,numnp      ! only use the vertices in the normal calc
               wnrm(i,:) = wnrm(iper(i),:)
            enddo
            wnmdb=wnrm
            if (numpe.gt.1) then
!.... tell other processors the resulting (and now complete) sums
               call commu (wnrm, ilwork, 3, 'out')
            endif
! Axisymmetric slaves need to be rotated back
            if(iabc==1) then    !are there any axisym bc's
               call rotabc(wnrm, iBC, 'out')
            endif
! Now all nodes have all the normal contributions for this surface,
! including those from off-processor and periodic nodes.
! We distribute these summed vectors to the proper slots in BCtmp,
! according to how many surfaces each node has seen so far
            do nn = 1, nshg
               if(nsurf(nn).eq.1)  &
                    BCtmp(nn,4:6)=BCtmp(nn,4:6)+wnrm(nn,:)
               if(nsurf(nn).eq.2) &
                    BCtmp(nn,8:10)=BCtmp(nn,8:10)+wnrm(nn,:)
               if(nsurf(nn).gt.2) &
                    BCtmp(nn,4:6)=BCtmp(nn,4:6)+wnrm(nn,:)
            enddo
! That's all for this surface; move on to the next
         enddo                  ! end loop over surface ID's

!     
!.... complete the averaging, adjust iBC, adjust BCtmp
!     
         wnrm(:,1:3)=zero
         do nn = 1, numnp       ! loop over all nodes
! leave nodes without wall-modeled surfaces unchanged
            if(nsurf(nn).eq.0) cycle
!
! mark this as a node with comp3
!
            ikp=0
            if(ibits(iBC(nn),3,3).eq.7 ) ikp=1
!         if(ibits(ibc(nn),3,3).eq.7 .and. BCtmp(nn,7).eq.zero) cycle
! find out which velocity BC's were set by the user, and cancel them
            ixset=ibits(iBC(nn),3,1)
            iyset=ibits(iBC(nn),4,1)
            izset=ibits(iBC(nn),5,1)
            iBC(nn)=iBC(nn)-ixset*8-iyset*16-izset*32
!
! save this stripped iBC for later un-fixing
!
            iBCSAV=iBC(nn)
!     
            if(abs(itwmod).eq.1) then ! slip velocity wall-model
! If we're using the slip-velocity wall-model, then the velocity
! boundary condition for this node will depend on how many modeled
! surfaces share it...
               if(nsurf(nn).eq.1) then ! 1 modeled wall
!   If this node is part of only one modeled wall, then the component
!   of the velocity normal to the surface is set to zero.  In this case
!   only the first vector of BCtmp received normal contributions
                  iBC(nn)=iBC(nn)+8 ! assume normal is x-dominated first
                  wnrm(nn,:)=BCtmp(nn,4:6)
                  if(abs(wnrm(nn,3)).gt.abs(wnrm(nn,2)))then ! z beats y
                     if(abs(wnrm(nn,3)).gt.abs(wnrm(nn,1)))then ! z beats x
                        iBC(nn)=iBC(nn)+24
                     endif      ! z beats x
                  endif         ! z beats y
                  if(abs(wnrm(nn,2)).ge.abs(wnrm(nn,3)))then ! y beats z
                     if(abs(wnrm(nn,2)).gt.abs(wnrm(nn,1)))then ! y beats x
                        iBC(nn)=iBC(nn)+8
                     endif      ! y beats x
                  endif         ! y beats z           !(adjusted iBC)
                  BCtmp(nn,7)=zero
               else if(nsurf(nn).eq.2) then
!   If this node is shared by two modeled walls, then each wall
!   provides a normal vector along which the velocity must be zero.
!   This leaves only one "free" direction, along which the flow is nonzero.
!   The two normal vectors define a plane.  Any pair of non-parallel
!   vectors in this plane can also be specified, leaving the same free 
!   direction.  Area-weighted average normal vectors for the two surfaces
!   sharing this node have been stored in BCtmp(4:6) and BCtmp(8:10)
!   We normalize the first of these, and then orthogonalize the second
!   against the first.  After then normalizing the second, we choose which
!   cartesian direction dominates each of the vectors, and adjust iBC 
!   and BCtmp accordingly
                  e1=BCtmp(nn,4:6)
                  wmag=e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3)
                  wmag=1./sqrt(wmag)
                  e1=e1*wmag    ! first vector is normalized
!
                  e2=BCtmp(nn,8:10)
                  wmag=e2(1)*e1(1)+e2(2)*e1(2)+e2(3)*e1(3) ! wmag=e2.n1
                  e2=e2-wmag*e1 ! second vector is orthogonalized
!     
                  wmag=e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3)
                  wmag=1./sqrt(wmag)
                  e2=e2*wmag    ! second vector is normalized
!
               ! idom tells us which component is currently dominant
               ! ivec(n) tells us which vector is dominated by comp n
                  ivec(:)=0     ! initialize
                  idom=1        ! assume x-comp dominates e1
                  bigcomp=abs(e1(1))
                  ivec(idom)=1
                  do i = 2, nsd
                     if(abs(e1(i)).gt.bigcomp) then
                        ivec(idom)=0
                        bigcomp=abs(e1(i))
                        idom=i
                        ivec(i)=1
                     endif
                  enddo
                  if(idom.ne.1) then
                     idom=1     ! assume x-comp dominates e2
                  else
                     idom=3     ! assume z-comp dominates e2
                  endif
                  bigcomp=abs(e2(idom))
                  
                  ivec(idom)=2
                  do i = 1, nsd
                     if(abs(e2(i)).gt.bigcomp) then
                        if(ivec(i).eq.1) cycle
                        ivec(idom)=0
                        bigcomp=abs(e2(i))
                        idom=i
                        ivec(i)=2
                     endif
                  enddo
               ! now we know which components dominate each vector
                  ixset=0       !
                  iyset=0       ! initialize
                  izset=0       !
                  if(ivec(1).ne.0) ixset=1
                  if(ivec(2).ne.0) iyset=1
                  if(ivec(3).ne.0) izset=1
                  ncomp=ixset+iyset+izset
                  if(ncomp.ne.2) write(*,*) 'WARNING: TROUBLE IN GENBC'
                  BCvecs(1,1:3)=e1(:)
                  BCvecs(2,1:3)=e2(:)
                  if((ixset.eq.1).and.(iyset.eq.1)) then
                     BCtmp(nn,4:6)=BCvecs(ivec(1),1:3)
                     BCtmp(nn,8:10)=BCvecs(ivec(2),1:3)
                  else if((ixset.eq.1).and.(izset.eq.1)) then
                     BCtmp(nn,4:6)=BCvecs(ivec(1),1:3)
                     BCtmp(nn,8:10)=BCvecs(ivec(3),1:3)
                  else if((iyset.eq.1).and.(izset.eq.1)) then
                     BCtmp(nn,4:6)=BCvecs(ivec(2),1:3)
                     BCtmp(nn,8:10)=BCvecs(ivec(3),1:3)
                  endif
                  iBC(nn)=iBC(nn)+ixset*8+iyset*16+izset*32
                  BCtmp(nn,7)=zero
                  BCtmp(nn,11)=zero
                  wnrm(nn,:)=BCtmp(nn,4:6)+BCtmp(nn,8:10)
               else if(nsurf(nn).gt.2) then
!     If this node is shared by more than two modeled walls, then
!     it is a corner node and it will be treated as having no slip
!     The normals for all surfaces beyond the first two were 
!     contributed to the first vector of BCtmp
                  wnrm(nn,:)=BCtmp(nn,4:6)+BCtmp(nn,8:10)
                  iBC(nn)=iBC(nn)+56
                  BCtmp(nn,7)=zero
               endif            ! curved surface
            else if(abs(itwmod).eq.2) then ! effective viscosity wall-model
! Otherwise, we're using the effective-viscosity wall-model and the 
! nodes on modeled surfaces are treated as no-slip
               iBC(nn)=iBC(nn)+56 ! set 3-comp
               if(itwmod.eq.-2) BCtmp(nn,7)=zero ! no slip
               if(nsurf(nn).eq.1) &
                    wnrm(nn,:)=BCtmp(nn,4:6)
               if(nsurf(nn).ge.2) &
                    wnrm(nn,:)=BCtmp(nn,4:6)+BCtmp(nn,8:10)
            endif
! Now normalize the wall normal for this node
            if(itwmod.ne.0) then
               wmag=sqrt(wnrm(nn,1)*wnrm(nn,1) &
                    +wnrm(nn,2)*wnrm(nn,2)+wnrm(nn,3)*wnrm(nn,3))
               wnrm(nn,:)=wnrm(nn,:)/wmag
            endif
!
!.... put back the comp3 info for bctmp
!
            if(ikp.eq.1) then
               iBC(nn)=iBCSAV+56 ! put it back to a comp3
               BCtmp(nn,:)=BCtmpSAV(nn,:) ! ditto
            endif
         enddo                  ! loop over all nodes
      endif                     ! end "if there are any surfID's"
!
!  If you are using the turbulence wall with axisymmetry you
!  need to modify the axisymmetry angle to account for the discrete
!  normals at the wall being different than the exact normals
!
! find the my normal, my partners normal and correct the angle
!$$$
!$$$        do i=1,numnp
!$$$           wmag = wnrm(i,1) * wnrm(i,1)
!$$$     &          + wnrm(i,2) * wnrm(i,2)
!$$$     &          + wnrm(i,3) * wnrm(i,3)
!$$$c
!$$$c  only "fix" wall nodes...other nodes still have a zero in wnrm
!$$$c
!$$$           if((btest(iBC(i),12)).and.(wmag.ne.zero)) then  
!$$$              BCtmp(i,1)=acos( wnrm(i,1) * wnrm(iper(i),1)
!$$$     &                       + wnrm(i,2) * wnrm(iper(i),2)
!$$$     &                       + wnrm(i,3) * wnrm(iper(i),3) )
!$$$           endif
!$$$        enddo
!
!.... return
!
      return
!
      end


      subroutine gensidcount(nsidg)
!---------------------------------------------------------------------
!
! This routine counts up the total number of surface ID's across
! all processors and makes a list of them
!
! Inputs:
!     iBCB        natural boundary condition switches and surfIDs
!
! Outputs:
!     nsidg       number of surface ID's globally (including all procs)
!     sidmapg     global list of surface ID's, lowest to highest
!
!---------------------------------------------------------------------
      use pointer_data ! access to miBCB
      use turbSA ! access to sidmapg
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      integer newflag, i, j
      integer nsidl             ! number of surface IDs on-proc
      integer nsidt             ! sum of all nsidl's
      integer nsidg             ! number of surface IDs globally
      integer nsid(numpe)       ! array of nsidl's
      integer idisp(numpe)      ! needed by mpi: see note below
      type llnod                ! structure for one node in a linked list
        integer :: value
        type (llnod), pointer :: next
      end type llnod
      type (llnod), pointer :: sidlist ! points to first elt of linked list
      type (llnod), pointer :: sidelt  ! points to generic elt of linked list
      integer, allocatable :: sidmapl(:) ! list of surfID's on-proc
      integer, allocatable :: temp(:)    ! temp space
!      integer iBCB(numelb,ndiBCB) ! column 1: naturalBC switches
                                   ! column 2: surface ID's
! Since we don't know a priori how many surface ID's there are, 
! on-processor or globally, we will store the ID's as an open-ended
! link list while we determine the total number of distinct ID's
      allocate (sidlist) ! allocate the first element of the sid 
                         ! linked list and point sidlist to it
      nsidl=0            ! initialize # of sids on this processor
      nsidg=0
      do iblk=1, nelblb  ! loop over boundary element blocks
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         do i = 1, npro         ! loop over boundary elements (belts)
            iBCB2=miBCB(iblk)%p(i,2)
            if(iBCB2.ne.zero) then ! if a sid is given for this belt
               if(nsidl.eq.0) then     !   if this is the first sid we've seen
                  nsidl=1              !       increment our count and
                  sidlist % value = iBCB2    ! record its value
               else                    !   if this isn't the first sid
                  newflag = 1          !     assume this is a new sid
                  sidelt => sidlist    !     starting with the first sid
                  do j = 1, nsidl      !     check the assumption
                     if((sidelt % value).eq.iBCB2) then
                        newflag=0      !        ...
                     endif
                     if(j.ne.nsidl) then !      ...
                        sidelt => sidelt % next
                     endif!                     ...
                  enddo
                  if(newflag.eq.1) then!     if it really is new to us
                     nsidl = nsidl + 1 !         increment our count
                     allocate (sidelt % next)!   tack a new element to the end
                     sidelt => sidelt % next!    point to the new element
                     sidelt % value = iBCB2    ! record the new sid
                  endif ! (really is a new sid)
               endif ! (first sid)
            endif ! (belt has a sid)
         enddo ! (loop over belts)
      enddo ! (loop over boundary element blocks)
! Copy the data from the linked list to a more easily-used array form
      if(nsidl.gt.0) then
         allocate( sidmapl(nsidl) )
         sidelt => sidlist      !     starting with the first sid
         do j = 1, nsidl
            sidmapl(j)=sidelt%value
            if(j.ne.nsidl) sidelt => sidelt%next
         enddo
      endif
! Gather the number of surface ID's on each processor
      if (numpe.gt.1) then      ! multiple processors
! write the number of surfID's on the jth processor to slot j of nsid
         call MPI_ALLGATHER(nsidl,1,MPI_INTEGER,nsid,1, &
                MPI_INTEGER,INEWCOMM,ierr)
! count up the total number of surfID's among all processes
         nsidt=0
         do j=1,numpe
            nsidt=nsidt+nsid(j)
         enddo
      else                      ! single processor
! the local information is the global information for single-processor
         nsid=nsidl
         nsidt=nsidl
      endif                     ! if-else for multiple processors
      if(nsidt.gt.0) then
!
!  Make all-processor surfID collage
!
! there will be some duplicate surface ID's when we gather, so we
! will use a temporary array
         allocate (temp(nsidt))
         if (numpe.gt.1) then   ! multiple processors
! we will gather surfID's from local on-proc sets to a global set
! we will stack each processor's surfID list atop that of the previous 
! processor.  If the list for processor i is called sidmapi, then our
! global coordinate list sidmap will look like this:
! ---------------------------------------------------------------
! | sidmap1       | sidmap2           | sidmap3   |   ...       |
! ---------------------------------------------------------------
!  <---nsid(1)---> <-----nsid(2)-----> <-nsid(3)->
!  <------------------------nsidt-----------------------...---->
! To accomplish this with MPI, we use MPI_ALLGATHERV, summarized as:
!MPI_ALLGATHERV(sendbuf,sendcount,sendtype,recvbuf,recvcount,disp,recvtype,comm) 
![ IN sendbuf] starting address of send buffer (choice) 
![ IN sendcount] number of elements in send buffer (integer) 
![ IN sendtype] data type of send buffer elements (handle) 
![ OUT recvbuf] address of receive buffer (choice) 
![ IN recvcount] number of elements received from each process (int array) 
![ IN disp] displacement array
![ IN recvtype] data type of receive buffer elements (handle) 
![ IN comm] communicator (handle)
! The displacement array disp is an array of integers in which the jth
! entry indicates which slot of sidmap marks beginning of sidmapj
! So, first we will build this displacement array
            idisp(:)=0      ! starting with zero, since MPI likes C-numbering
            do j=2,numpe
               idisp(j)=idisp(j-1)+nsid(j-1) ! see diagram above
            enddo
! Now, we gather the data
            call MPI_ALLGATHERV(sidmapl(:),nsidl, &
                 MPI_INTEGER,temp(:),nsid,idisp, &
                 MPI_INTEGER,INEWCOMM,ierr)
! sort surfID's, lowest to highest
            isorted = 0
            do while (isorted.eq.0) ! while the list isn't sorted
               isorted = 1      ! assume the list is sorted this time
               do j = 2, nsidt  ! loop over ID's
                  if(temp(j).lt.temp(j-1)) then ! ID exceeds predecessor
                     itmp=temp(j-1) 
                     temp(j-1)=temp(j)
                     temp(j)=itmp
                     isorted=0  ! not sorted yet
                  endif
               enddo            !loop over ID's
            enddo               ! while not sorted
! determine the total number of distinct surfID's, globally
            nsidg=nsidt         ! assume there are no duplicate SurfID's
            do j = 2, nsidt
               if(temp(j).eq.temp(j-1)) nsidg = nsidg - 1 ! correction
            enddo
! create duplicate-free surfID list
            allocate( sidmapg(nsidg) )
            sidmapg(1)=temp(1)
            nsidg = 1
            do j = 2, nsidt
               if(temp(j).ne.temp(j-1)) then
                  nsidg = nsidg + 1
                  sidmapg(nsidg)=temp(j)
               endif
            enddo
            deallocate( temp )
         else                   ! single-processor
! global data is local data in single processor case
            nsidg=nsidl
            allocate( sidmapg(nsidg) )
            sidmapg=sidmapl
         endif                  ! if-else multiple processor
!     
      endif                     ! if at least one surfid
!
      return
!
      end



      subroutine genotwn(x, BCtmp, iBC, nsurf)
!
!----------------------------------------------------------------------
!  This routine determines which node to use as the first node off the
!  wall in near-wall modeling traction calculations for each wall node.  
!  Each wall node has a normal, calculated from the wall elements to
!  which that node belongs.  We seek the node within the boundary
!  element that lies closest to the line defined by the normal vector.
!  We create normalized vectors pointing from the wall node in 
!  question to each of the nodes in the boundary element. The vector
!  that has the largest projection onto the wall node's normal vector
!  points to the node we want.  Nodes that are not on a wall point to
!  themselves as their own off-the-wall node.
!
! input:
!  x     (nshg,3)     :        : nodal position vectors
!  wnrm  (nshg,3)     : (mod)  : normal vectors for each node
!  iBCB5 (nshg)       : (file) : wall flag for belt
!  ienb  (numelb,nen) : (file) : connectivity for belts
!
! output:
!  otwn  (nshg)       : (mod)  : off-the-wall nodes for each node
!
!----------------------------------------------------------------------
!
      use pointer_data          ! used for mienb, miBCB
      use turbSA
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      integer iel, nod, can
      real*8 vec(3), leng, dp, bigdp, lil
      real*8 x(numnp,nsd),BCtmp(nshg,ndof+7)
      integer iBC(nshg), nsurf(nshg)
      integer gbits
      integer, allocatable :: ienb(:)

      allocate( otwn(nshg) )
!
! initialize otwn to oneself
!
      do nod = 1, nshg
         otwn(nod)=nod
      enddo
!
! determine otwn for each wall node
!
      do iblk=1, nelblb         ! loop over boundary element blocks
         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
         nenl  = lcblkb(5,iblk)
         nenbl = lcblkb(6,iblk)
         nshl = lcblkb(9,iblk)
         allocate( ienb(nshl) )
         do i = 1, npro         ! loop over boundary elements
            iBCB1=miBCB(iblk)%p(i,1)
            iBCB2=miBCB(iblk)%p(i,2)
            ienb(1:nshl)=mienb(iblk)%p(i,1:nshl)
            if (btest(iBCB1,4)) then ! (on the wall)
               do nod = 1, nenbl !   for each wall node in this belt
                  bigdp = zero  !     initialize largest dot product
                  do can=nenbl+1,nenl !  loop over off-the-wall candidates
                     nn=ienb(can)
                               !       don't bother with wall nodes
                     if(nsurf(nn).ne.0) cycle
                               !       don't bother with no-slip nodes
                     if(ibits(iBC(nn),3,3).eq.7 .and.  &
                          BCtmp(nn,7).eq.zero) cycle
!                              !       calculate candidate vector
                     vec(:)=x(ienb(can),:)-x(ienb(nod),:)
                     leng=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
                     vec(:)=vec(:)/leng
!                              !       calculate dot product with wnrm
                     vec(:)=vec(:)*wnrm(ienb(nod),:)
                     dp=vec(1)+vec(2)+vec(3)
!                              !       set candidate as otwn if necessary
!                              !       (wnrm points into fluid, so
!                              !        we want the most positive dp)
                     if (dp.gt.bigdp) then
                        otwn(ienb(nod)) = ienb(can)
                        bigdp=dp
                     endif
                  enddo         !(loop over off-the-wall candidates)
               enddo            !(loop over wall nodes in current belt)
            endif
         enddo                  !(loop over elts in block)
         deallocate(ienb)
      enddo                     !(loop over boundary element blocks)
      do nn = 1, nshg
         if((otwn(nn).eq.nn).and.(nsurf(nn).ne.0)) then ! if a node on a
                                                    ! modeled surface
                                                    ! didn't find a node
                                                    ! off the wall
            lil = 1.0e32 ! distance to current closest prospect
            do can = 1, nshg ! loop over nodes
               if(nsurf(can).eq.0) then ! if this candidate is off the
                                        ! wall
                  if(ibits(iBC(can),3,3).eq.7 .and.  &
                     BCtmp(can,7).eq.zero) then  ! no slip nodes not allowed
                  else          ! not a forbidden node
                     vec(:)=x(nn,:)-x(can,:)
                     leng=vec(1)**2+vec(2)**2+vec(3)**2
                     if(leng.lt.lil) then ! this candidate is closest so far
                        lil=leng
                        otwn(nn)=can
                     endif      ! closest so far
                  endif  ! end of no slip nodes
               endif ! off-the-wall check
            enddo ! loop over nodes
         endif ! lonely wall-node check
      enddo ! loop over nodes
!
      return
!
      end

