        subroutine gendat (y,       ac,       x,      iBC,     BC, &
                           iper,    ilwork, &
                           shp,     shgl,    shpb,    shglb, &
                           ifath,   velbar,   nsons ) 
!
!----------------------------------------------------------------------
!
! This routine inputs the geometry and the boundary conditions.
!
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
      
        use dtnmod
        use pointer_data
        use deformableWall
        use phcommonvars
        use multidomain     
        use ale    
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
! arrays in the following line are now dimensioned in readnblk
!
        dimension y(nshg,ndof),      ac(nshg,ndof), &
                  x(numnp,nsd),      iBC(nshg), &
                  BC(nshg,ndofBC), &
                  nodflx(numflx),    ilwork(nlwork), &
                  iper(nshg)
!
!.... shape function declarations
!     
        dimension shp(MAXTOP,maxsh,MAXQPT),   &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                  shpb(MAXTOP,maxsh,MAXQPT), &
                  shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!
!  stuff for dynamic model s.w.avg and wall model
!
        dimension ifath(numnp),    velbar(nfath,nflow), nsons(nfath)
        
        dimension lnode(27) 
        
        integer wnodestp(numel*nshlb)
        integer wnodesgtlmap(nshg)
        integer wnodefnd
        
!
!.... start the timer
!
        
!AD        call timer ('PrProces')
!
!.... ---------------------------->  Nodes  <--------------------------
!
!.... compute length scales
!
        call xyzbound(x)
!
!.... echo the coordinates
!
        if ((necho .lt. 2).and.(myrank.eq.master)) then
          do n = 1, numnp
            if (mod(n,50) .eq. 1) write (iecho,1000) ititle,(i,i=1,nsd)
            write (iecho,1100) n, (x(n,i),i=1,nsd)
          enddo
        endif
!
!.... prepare periodic boundary conditions
!
        do i = 1,nshg
          if (iper(i) .ne. 0) then
            nshg0 = nshg0 - 1
          else
            iper(i) = i
          endif
        enddo
!
!.... ---------------------->  Interior Elements  <--------------------
!
        ibound = 0
!
!.... generate the interior nodal mapping
!
        call genshp ( shp, shgl, nshape, nelblk)
!
!.... --------------------->  Boundary Conditions  <-------------------
!
!.... read and generate the boundary condition codes (iBC array)
!
        call geniBC (iBC)
!
!.... read and generate the essential boundary conditions (BC array)
!
        call genBC  (iBC,   BC,   x, &
                     ilwork, iper)
!
!.... ---------------------->  Boundary Elements  <--------------------
!
        ibound = 1
        call gtnods
!
!  We now take care of Direchlet to Neumann BC's.  It had to move here
!  so that the IBC array was of size nshg and ready to be marked.
!

        if(nsclr.gt.0) then 
           call initDtN         ! Dirichlet to Neumann module: 
                                     ! initialize this once only
           do iblk = 1, nelblb  ! number of blocks
              iel    = lcblkb(1,iblk)
              npro   = lcblkb(1,iblk+1) - iel
!
!  for the DtN BC we need to mark all of the nodes that are involved.
!
              do i=1,npro
!
! if this element has the BCB AND it has not been found yet then mark it
!
                 if(miBCB(iblk)%p(i,2).lt.0) then  
                    idtn = 1    !set the flag for dtn bc's
                    do j=1,nshapeb
                       do isclr=1,nsclr
                          ignd=mienb(iblk)%p(i,j)
                             ifeature(ignd) = abs(miBCB(iblk)%p(i,2))       
                             iBC(ignd)=ior(iBC(ignd),2**13)
                                ! must mark this as a Neumann BC now
                             miBCB(iblk)%p(i,1)= &
                             ior(miBCB(iblk)%p(i,1),2**(4+isclr))
                       end do
                    end do
                 endif
              end do
           end do
        endif
!
!.... generate the boundary element shape functions
!
        call genshpb ( shpb, shglb, nshapeb, nelblb)
!.... Evaluate the shape funcs. and their gradients at the desired quadrature
!.... for filtering. Save these evaluations using a module
!
! KEJ moved them to this point because cdelsq now passed with module
!     and it is read in with velb.<stepnum>.<proc#> now
!
!        if (iLES .gt. 0) then
!
!           call setfilt         ! For setting quad. rule to use for integrating
!           call filtprep        ! the hat filter.
!           if(iLES/10 .eq. 2) then
!              call setave       ! For averaging cdelsq computed at quad pts
!              call aveprep(shp,x)
!           endif
!        endif
!
! User sets request pzero in solver.inp now
!
!        call genpzero(iBC,iper)
!
!      if((myrank.eq.master).and.(irscale.ge.0)) then
!         call setSPEBC(numnp,nsd)
!	 call eqn_plane(x, iBC)
!      endif
      
!
! Here we find the nodes on the deformable wall
!      
      nwnp = 0
      wnodestp = 0

      call getbnodes(lnode)
      do iblk = 1, nelblb
         iel = lcblkb(1,iblk)
         npro = lcblkb(1,iblk+1) - iel
         do i=1,npro
!            write(*,*) btest(miBCB(iblk)%p(i,1),4)
            if (btest(miBCB(iblk)%p(i,1),4)) then ! check element deformable
               do j=1,nshlb
               
                  n = lnode(j)
               
                  wnodefnd = 0
                  do k=1,nwnp
                     if (wnodestp(k).eq.mienb(iblk)%p(i,n)) then
                        wnodefnd = 1  
                        exit
                     end if
                  end do
               
                  if (wnodefnd.eq.0) then
                     nwnp = nwnp+1
                     wnodestp(nwnp) = mienb(iblk)%p(i,n)
!
! this mapping takes a global node number and returns a wall node number
! from 1 to the nwnp (the number of nodes on the wall)
!                     
                     wnodesgtlmap(mienb(iblk)%p(i,n)) = nwnp

                  end if
               
               end do                  
            end if
         end do
      end do     
      
!    
! Copy to the global array
!          
      if (nwnp .gt. 0) then
        allocate(mWNodes%p(nwnp))
      endif
      allocate(mWNodes_gtlmap%p(nshg))
      
      do i=1,nwnp
         mWNodes%p(i) = wnodestp(i)
      end do
      
      do i=1,nshg
         mWNodes_gtlmap%p(i) = wnodesgtlmap(i)
      end do
!
! *** start of intialise multidomain models
!
      call startmultidomain()     
      call startmultidomain(iheart,'heart')       
      !call startmultidomain(isystemic,'systemic')
      if (numControlledCoronarySrfs .gt. int(0)) then
        call startmultidomain(numControlledCoronarySrfs,'coronary')
      endif
      if (numNetlistLPNSrfs .gt. int(0)) then
        call startmultidomain(numNetlistLPNSrfs,'netlist')
      endif
      
      ! write out status 
      call multidomainstatus()                              
!
! *** initialise the multidomain container
!
      if (multidomainactive) then
         multidom = multidomconstructor()   
      end if
!
! *** end of intialise multidomain models   
!
!
! *** initialise ALE 
!
      call readGlobalMeshVelocity() !
      write (*,*) "globalMeshVelocity = ",globalMeshVelocity
      ! write (*,*) "rigidOn =",rigidOn
      ! write (*,*) "globalRigidVelocity = ",globalRigidVelocity

      call addGlobalMeshVelocityToSolution(y,nshg,ndof)
!
! *** end of initialise ALE 
!
!
!.... --------------------->  Initial Conditions  <--------------------
!
!.... generate the initial conditions and initialize time varying BC
!
        call genini (iBC,      BC,         y,  &
                     ac,       iper,  &
                     ilwork,   ifath,      velbar,   &
                     nsons,    x, &
                     shp,     shgl,    shpb,    shglb) 
!
!.... close the geometry, boundary condition and material files
!
        close (igeom)
        close (ibndc)
        if (mexist.eq.1) close (imat)
!
!.... return
!
!AD        call timer ('Back    ')
        return
!
!.... end of file error handling
!
999     call error ('gendat  ','end file',igeom)
!
1000    format(a80,//, &
        ' N o d a l   C o o r d i n a t e s                  ',//, &
        '    Node     ',12x,3('x',i1,:,17x))
1100    format(1p,2x,i5,13x,3(1e12.5,7x))
2000    format(a80,//, &
        ' B o u n d a r y   F l u x   N o d e s              '//, &
        '   index          Node          ')
2100    format(1x,i5,5x,i10)
!
        end subroutine gendat


        subroutine xyzbound(x)

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
        !include "auxmpi.h"

        dimension x(numnp,3)

        real*8   Forout(3), Forin(3)

        xlngth=maxval(x(:,1))
        ylngth=maxval(x(:,2))
        zlngth=maxval(x(:,3))
        if(numpe .gt. 1) then
           Forin=(/xlngth,ylngth,zlngth/)
           call MPI_ALLREDUCE (Forin, Forout, 3, &
             MPI_DOUBLE_PRECISION,MPI_MAX, INEWCOMM,ierr)
           xmax = Forout(1)
           ymax = Forout(2)
           zmax = Forout(3)
        else
           xmax = xlngth
           ymax = ylngth
           zmax = zlngth
        endif
        xlngth=minval(x(:,1))
        ylngth=minval(x(:,2))
        zlngth=minval(x(:,3))
        if(numpe .gt. 1) then
           Forin=(/xlngth,ylngth,zlngth/)
           call MPI_ALLREDUCE (Forin, Forout, 3, &
             MPI_DOUBLE_PRECISION,MPI_MIN, INEWCOMM,ierr)
        else
           Forout(1) = xlngth
           Forout(2) = ylngth
           Forout(3) = zlngth
        endif

        xlngth = xmax-Forout(1)
        ylngth = ymax-Forout(2)
        zlngth = zmax-Forout(3)

        if(myrank.eq.master) then
           print 108,  xlngth,ylngth,zlngth
        endif
 108    format(' Domain size (x,y,z):',2x,3f15.10)
        return
        end
