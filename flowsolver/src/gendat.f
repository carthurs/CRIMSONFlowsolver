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
        subroutine gendat (y,       ac,       x,      iBC,     BC,
     &                     iper,    ilwork,
     &                     shp,     shgl,    shpb,    shglb,
     &                     ifath,   velbar,   nsons ) 
c
c----------------------------------------------------------------------
c
c This routine inputs the geometry and the boundary conditions.
c
c----------------------------------------------------------------------
c
      
        use readarrays          ! used to acess nBC
        use dtnmod
        use pointer_data
        use deformableWall
        include "common.h"
c
c arrays in the following line are now dimensioned in readnblk
c        dimension nBC(nshg)
c
        dimension y(nshg,ndof),      ac(nshg,ndof),
     &            x(numnp,nsd),      iBC(nshg),
     &            BC(nshg,ndofBC),
     &            nodflx(numflx),    ilwork(nlwork),
     &            iper(nshg)
c
c.... shape function declarations
c     
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,nflow), nsons(nfath)
        
        dimension lnode(27) 
        
        integer wnodestp(numel*nshlb)
        integer wnodesgtlmap(nshg)
        integer wnodefnd
        
c
c.... start the timer
c
        
CAD        call timer ('PrProces')
c
c.... ---------------------------->  Nodes  <--------------------------
c
c.... compute length scales
c
        call xyzbound(x)
c
c.... echo the coordinates
c
        if ((necho .lt. 2).and.(myrank.eq.master)) then
          do n = 1, numnp
            if (mod(n,50) .eq. 1) write (iecho,1000) ititle,(i,i=1,nsd)
            write (iecho,1100) n, (x(n,i),i=1,nsd)
          enddo
        endif
c
c.... prepare periodic boundary conditions
c
        do i = 1,nshg
          if (iper(i) .ne. 0) then
            nshg0 = nshg0 - 1
          else
            iper(i) = i
          endif
        enddo
c
c.... ---------------------->  Interior Elements  <--------------------
c
        ibound = 0
c
c.... generate the interior nodal mapping
c
        call genshp ( shp, shgl, nshape, nelblk)
c
c.... --------------------->  Boundary Conditions  <-------------------
c
c.... read and generate the boundary condition codes (iBC array)
c
        call geniBC (iBC)
c
c.... read and generate the essential boundary conditions (BC array)
c
        call genBC  (iBC,   BC,   point2x,
     &               point2ilwork, point2iper)
        deallocate(nBC)
c
c.... ---------------------->  Boundary Elements  <--------------------
c
        ibound = 1
        call gtnods
c
c  We now take care of Direchlet to Neumann BC's.  It had to move here
c  so that the IBC array was of size nshg and ready to be marked.
c

        if(nsclr.gt.0) then 
           call initDtN         ! Dirichlet to Neumann module: 
                                     ! initialize this once only
           do iblk = 1, nelblb  ! number of blocks
              iel    = lcblkb(1,iblk)
              npro   = lcblkb(1,iblk+1) - iel
c
c  for the DtN BC we need to mark all of the nodes that are involved.
c
              do i=1,npro
c
c if this element has the BCB AND it has not been found yet then mark it
c
                 if(miBCB(iblk)%p(i,2).lt.0) then  
                    idtn = 1    !set the flag for dtn bc's
                    do j=1,nshapeb
                       do isclr=1,nsclr
                          ignd=mienb(iblk)%p(i,j)
                             ifeature(ignd) = abs(miBCB(iblk)%p(i,2))       
                             iBC(ignd)=ior(iBC(ignd),2**13)
                                ! must mark this as a Neumann BC now
                             miBCB(iblk)%p(i,1)=
     &                       ior(miBCB(iblk)%p(i,1),2**(4+isclr))
                       end do
                    end do
                 endif
              end do
           end do
        endif
c
c.... generate the boundary element shape functions
c
        call genshpb ( shpb, shglb, nshapeb, nelblb)
c.... Evaluate the shape funcs. and their gradients at the desired quadrature
c.... for filtering. Save these evaluations using a module
c
c KEJ moved them to this point because cdelsq now passed with module
c     and it is read in with velb.<stepnum>.<proc#> now
c
        if (iLES .gt. 0) then

           call setfilt         ! For setting quad. rule to use for integrating
           call filtprep        ! the hat filter.
           if(iLES/10 .eq. 2) then
              call setave       ! For averaging cdelsq computed at quad pts
              call aveprep(shp,x)
           endif
        endif
c
c User sets request pzero in solver.inp now
c
c        call genpzero(iBC,iper)
c
      if((myrank.eq.master).and.(irscale.ge.0)) then
         call setSPEBC(numnp,nsd) 	 
	 call eqn_plane(point2x, iBC)
      endif
      
c
c Here we find the nodes on the deformable wall
c      
      nwnp = 0
      wnodestp = 0
      
      call getbnodes(lnode)

      do iblk = 1, nelblb
         iel = lcblkb(1,iblk)
         npro = lcblkb(1,iblk+1) - iel
         do i=1,npro
c            write(*,*) btest(miBCB(iblk)%p(i,1),4)
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
c
c this mapping takes a global node number and returns a wall node number
c from 1 to the nwnp (the number of nodes on the wall)
c                     
                     wnodesgtlmap(mienb(iblk)%p(i,n)) = nwnp
                  end if
               
               end do                  
            end if
         end do
      end do     
      
c    
c Copy to the global array
c          
      allocate(mWNodes%p(nwnp))
      allocate(mWNodes_gtlmap%p(nshg))
      
      do i=1,nwnp
         mWNodes%p(i) = wnodestp(i)
      end do
      
      do i=1,nshg
         mWNodes_gtlmap%p(i) = wnodesgtlmap(i)
      end do
      
c
c.... --------------------->  Initial Conditions  <--------------------
c
c.... generate the initial conditions and initialize time varying BC
c
        call genini (iBC,      BC,         y, 
     &               ac,       iper, 
     &               ilwork,   ifath,      velbar,  
     &               nsons,    x,
     &               shp,     shgl,    shpb,    shglb) 
c
c.... close the geometry, boundary condition and material files
c
        close (igeom)
        close (ibndc)
        if (mexist) close (imat)
c
c.... return
c
CAD        call timer ('Back    ')
        return
c
c.... end of file error handling
c
999     call error ('gendat  ','end file',igeom)
c
1000    format(a80,//,
     &  ' N o d a l   C o o r d i n a t e s                  ',//,
     &  '    Node     ',12x,3('x',i1,:,17x))
1100    format(1p,2x,i5,13x,3(1e12.5,7x))
2000    format(a80,//,
     &  ' B o u n d a r y   F l u x   N o d e s              '//,
     &  '   index          Node          ')
2100    format(1x,i5,5x,i10)
c
        end


        subroutine xyzbound(x)

        include "common.h"
        include "mpif.h"
        include "auxmpi.h"

        dimension x(numnp,3)

        real*8   Forout(3), Forin(3)

        xlngth=maxval(x(:,1))
        ylngth=maxval(x(:,2))
        zlngth=maxval(x(:,3))
        if(numpe. gt. 1) then
           Forin=(/xlngth,ylngth,zlngth/)
           call MPI_ALLREDUCE (Forin, Forout, 3,
     &       MPI_DOUBLE_PRECISION,MPI_MAX, MPI_COMM_WORLD,ierr)
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
           call MPI_ALLREDUCE (Forin, Forout, 3,
     &       MPI_DOUBLE_PRECISION,MPI_MIN, MPI_COMM_WORLD,ierr)
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
