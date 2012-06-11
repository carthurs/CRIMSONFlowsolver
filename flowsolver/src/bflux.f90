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
      subroutine Bflux ( y,          ac,        u,      
     &                   x,          xdist,     xdnv, 
     &                   shp,        shgl,      shpb,   
     &                   shglb,      ilwork,    iBC,
     &                   BC,         iper  )
c
c----------------------------------------------------------------------
c
c This routine :
c   1. computes the boundary fluxes
c   2. prints the results in the file [FLUX.lstep]
c
c output:
c  in file flux.<lstep>.n (similar to restart file):
c     machin  nshg  lstep 
c     normal_1 ... normal_nsd            ! outward normal direction
c     tau_1n   ... tau_nsd n             ! boundary viscous flux
c
c----------------------------------------------------------------------
c
      
      use pointer_data
      use deformableWall
      use LagrangeMultipliers 
      
      include "common.h"
      include "mpif.h"

      character*5  cname
      
      real*8    y(nshg,ndof),             ac(nshg,ndof),
     &          u(nshg,nsd),              
     &          x(numnp,nsd),
     &          xdist(numnp),              
     &          xdnv(numnp,nsd)
      dimension iBC(nshg),           
     &          BC(nshg,ndofBC),  
     &          iper(nshg)
     
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &          shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &          shpb(MAXTOP,maxsh,MAXQPT),
     &          shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c    
c
      real*8    flxres(nshg,nflow),
     &          flxLHS(nshg,1),           flxnrm(nshg,nsd),
     &          temp(nshg),               rtmp(nshg,ndof),
     &          flxTot(nflow),            wallssVec(nshg,ndof)

      real*8    qres(nshg,nsd*nsd)

c
      integer   ilwork(nlwork),
     &          invflx(nshg),             nodflx(nshg)             
c
      character*20 fname1,  fmt1, fmt2, fnamer
      character*25 fname2
      integer irstin, isize, nitems
      integer iarray(50)  ! integers read from headers

      real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC
      integer, allocatable, dimension(:,:)    :: ien2
      integer, allocatable, dimension(:)      :: map
      real*8, allocatable, dimension(:,:)     :: xmu2
c
c....  calculate the flux nodes
c
      numflx  = 0
      invflx  = 0
      nodflx  = 0
      do iblk = 1, nelblb 
         iel    = lcblkb(1,iblk)
         lcsyst = lcblkb(3,iblk)
         nenl   = lcblkb(5,iblk) 
         nenbl  = lcblkb(6,iblk) 
         ndofl  = lcblkb(8,iblk)
         nshl   = lcblkb(9,iblk)
         nshlb  = lcblkb(10,iblk)
         npro   = lcblkb(1,iblk+1) - iel 
         call flxNode (mienb(iblk)%p,   miBCB(iblk)%p,  invflx) 
      enddo 
c
c      do i = 1, nshg
      do i = 1, numnp
         if (invflx(i) .ne. 0) then
            numflx = numflx + 1
            nodflx(numflx) = i
         endif
      enddo
c     
c.... -------------------->   interior elements   <--------------------
c     
c.... initialize the arrays
c     
      flxres = zero
      flxLHS = zero
      flxnrm = zero

      if (numflx .ne. 0)  then !we have flux nodes
         qres   = zero
c     
c.... loop over the element-blocks
c
         lhs    = 0

         ires=2  ! shield e3ql from an unmapped lmassinv
         ierrcalcsave=ierrcalc
         ierrcalc=0
         do iblk = 1, nelblk
c     
c.... set up the parameters
c     
            iel    = lcblk(1,iblk)
            nenl   = lcblk(5,iblk) ! no. of vertices per element
            nshl   = lcblk(10,iblk)
            ndofl  = lcblk(8,iblk)
            lcsyst = lcblk(3,iblk)
            npro   = lcblk(1,iblk+1) - iel 
            ngauss = nint(lcsyst)       
            allocate ( ien2(npro,nshl) )
            allocate ( xmu2(npro,maxsh))
            allocate ( map(npro) )
c
c.... get the elements touching the boundary
c         
            call mapConn( mien(iblk)%p,    ien2,    invflx,
     &                    map,             nshl,    npro,    
     &                    npro2,           nshg )

            nprold = npro
            npro = npro2
         
            if (npro .ne. 0) then

               call mapArray( mxmudmi(iblk)%p, xmu2,    map,
     &                        maxsh,           nprold)
c
c.... allocate the element matrices (though they're not needed)
c
               allocate ( xKebe(npro,9,nshl,nshl) )
               allocate ( xGoC (npro,4,nshl,nshl) )
               if(Lagrange.gt.zero) then
                  allocate(loclhsLag(npro,9,nshlb,nshlb,3))
               endif 
c     
c.... compute and assemble the residuals
c     
               call AsIGMR (y,                    ac,
     &                      x,                    xmu2(1:npro,:),
     &                      shp(lcsyst,1:nshl,:),
     &                      shgl(lcsyst,:,1:nshl,:),
     &                      ien2(1:npro,:),       
     &                      flxres,               qres,
     &                      xKebe,                xGoC,
     &                      rtmp)
c     
               deallocate ( xKebe )
               deallocate ( xGoC  )
               if(Lagrange.gt.zero) then
                  deallocate(loclhsLag)
               endif
            endif
            deallocate ( ien2  )
            deallocate ( xmu2  )
            deallocate ( map   )
c     
         enddo ! iblk = 1, nelblk
         ierrcalc=ierrcalcsave
c     
c.... -------------------->   boundary elements   <--------------------
c     
         do iblk = 1, nelblb
c     
c.... set up the parameters
c
            iel    = lcblkb(1,iblk)
            lcsyst = lcblkb(3,iblk)
            nenl   = lcblkb(5,iblk)
            nshl   = lcblkb(9,iblk)
            nenbl  = lcblkb(6,iblk)
            nshlb  = lcblkb(10,iblk)
            npro   = lcblkb(1,iblk+1) - iel 
 
            if(lcsyst.eq.3) lcsyst=nenbl
c     
            if(lcsyst.eq.3 .or. lcsyst.eq.4) then
               ngaussb = nintb(lcsyst)
            else
               ngaussb = nintb(lcsyst)
            endif
c
c.... allocate the element matrices (though they're not needed)
c
            allocate ( xKebe(npro,9,nshl,nshl) )   
c.... compute and assemble the residuals
c
            call AsBFlx (u,                       y,
     &                   ac,                      
     &                   x,
     &                   xdist,
     &                   xdnv,
     &                   shpb(lcsyst,1:nshl,:),
     &                   shglb(lcsyst,:,1:nshl,:),
     &                   mienb(iblk)%p,
     &                   miBCB(iblk)%p,           mBCB(iblk)%p,
     &                   invflx,                  flxres,
     &                   flxLHS,                  flxnrm,
     &                   xKebe,                   
     &                   mSWB(iblk)%p,            mTWB(iblk)%p,
     &                   mEWB(iblk)%p, 
     &                   mPS_global(iblk)%p,
     &                   mKwall_xKebe(iblk)%p)
c     
            deallocate ( xKebe )
c     
c.... end of boundary element loop
c
         enddo !iblk = 1, nelblb

      else
c         print *, "in Bflux: partition ", myrank, " has no flux nodes!"
      endif  ! make sure the zero numflux processors still commu

c.... Communication needed before we take care of periodicity and
c     division of RHS by LHS ???
cpf: note that the domains that do not have flux nodes,
c    have zero flxres, flxLHS, and flxnrm vectors
c
      if ( numpe > 1 ) then
         call commu (flxres, ilwork, nflow, 'in ')
         call commu (flxLHS, ilwork, 1   , 'in ')
         call commu (flxnrm, ilwork, nsd , 'in ')
      endif
c
c  take care of periodic boundary conditions
cpf: check this!
c
      do j= 1,nshg
         if ((btest(iBC(j),10))) then
            i = iper(j)
            flxLHS(i,1) = flxLHS(i,1) + flxLHS(j,1)
            flxres(i,:) =  flxres(i,:) + flxres(j,:)
         endif
      enddo
c
      do j= 1,nshg
         if ((btest(iBC(j),10))) then
            i = iper(j)
            flxLHS(j,1) = flxLHS(i,1)
            flxres(j,:) = flxres(i,:)
         endif
      enddo
c
c        call bc3per(iBC,  flxres, iper, ilwork, nflow)

c
c.... integrated fluxes (aerodynamic forces update)
c
      flxTot = zero
      do n = 1, numflx
         flxTot = flxTot + flxres(nodflx(n),:)
      enddo
      Force(1) = flxTot(1)
      Force(2) = flxTot(2)
      Force(3) = flxTot(3)
c
c.... only need to commu if we are going to print surface flux since
c     the force calculation just sums flxres (and each on processor node
c     has his "piece" of the sum already).
c
      ntoutv=ntout
      if ( (irs .ge. 1) 
     &     .and. (mod(lstep, ntoutv) .eq. 0) 
     &     .or.  (istep .eq. nstep(itseq)) ) then

c
c  need to zero the slaves to prevent counting twice
c  (actually unnecessary since flxres of boundary nodes will be counted n
c  times while flxlhs will be counted n times-> the ratio is still
c  correct
c      
         wallssVec=rtmp

         if (numflx .eq. 0) then   !no flux nodes
            rtmp=zero
            wallssVec = zero
         else
c     
c.... ---------------------------->  Solve  <---------------------------
c
c.... compute the viscous and heat fluxes
c     
c
c.... ---------------------------->  Print  <---------------------------
c
c.... nodal fluxes
c
            do i = 1, 3
               where ( (invflx .ne. 0) .and. (flxLHS(:,1) .ne. zero) )
     &              flxres(:,i) = flxres(:,i) / flxLHS(:,1)
            enddo
c     
c.... normalize the outward normal
c     
            temp = sqrt( flxnrm(:,1)**2 
     &                 + flxnrm(:,2)**2 
     &                 + flxnrm(:,3)**2 )
            where ( (invflx .ne. 0) .and. (temp .ne. zero) )
               flxnrm(:,1) = flxnrm(:,1) / temp
               flxnrm(:,2) = flxnrm(:,2) / temp
               flxnrm(:,3) = flxnrm(:,3) / temp
            endwhere
         endif !no flux nodes
         
c     
c.... ---------------------------->  Communications <-------------------
c
         if(numpe > 1) then
c           print *, "in Bflux: ", myrank, 
c    &               " is calling commu (3 times)"
            call commu (flxres, ilwork, nflow, 'out')
            call commu (flxLHS, ilwork, 1   , 'out')
            call commu (flxnrm, ilwork, nsd , 'out')
         endif
c
         rtmp = zero
         wallssVec  = zero

         do i=1, numnp
            if (invflx(i) .ne. 0) then
               rtmp(i,2:4) = flxres(i,1:3) !viscous flux
c     calculate the WSS
               tn = flxres(i,1) * flxnrm(i,1)
     &            + flxres(i,2) * flxnrm(i,2)
     &            + flxres(i,3) * flxnrm(i,3)

                wallssVec(i,1) = flxres(i,1) - tn * flxnrm(i,1)
                wallssVec(i,2) = flxres(i,2) - tn * flxnrm(i,2)
                wallssVec(i,3) = flxres(i,3) - tn * flxnrm(i,3)
            endif
         enddo

c         itmp = 1
c         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
c         write (fmt1,"('(''flux.'',i',i1,',1x)')") itmp
c         write (fname1,fmt1) lstep
      
c         fname1 = trim(fname1) // cname(myrank+1)
   
c        print *, "in Bflux: ",myrank," is opening flux file", fname1

c         open (unit=iflux, file=fname1, status='unknown', 
c     &         form='formatted',err=997)

c      write (iflux) machin, nshg, lstep
c      write (iflux) rtmp(:,1:6)
c
c.... output the results
c     
c         do n = 1, numflx
c            k = nodflx(n)
c            write (iflux,2000) k, (x(k,i), i=1,3), 
c     &           (flxnrm(k,i),  i=1,3),
c     &           (flxres(k,i),  i=1,3)
c         enddo
c         close (iflux)

c        print *, "in Bflux: ",myrank," closed flux file", fname1

c... output the results in the new format in restart.step#.proc# file

         itmp = 1
         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
         write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
         write (fname2,fmt2) lstep

         fname2 = trim(fname2) // cname(myrank+1)
c
c.... open input files
c
         call openfile(  fname2,  'append?', irstin )
         
         fnamer = 'boundary flux'        
         isize = nshg*ndof
         nitems = 3;
         iarray(1) = nshg
         iarray(2) = ndof
         iarray(3) = lstep
         call writeheader(irstin, fnamer,iarray, nitems, isize, 
     &        'double', iotype )
    
c         fnamer = 'boundary flux'        
         nitems = nshg*ndof
         call writedatablock(irstin, fnamer,rtmp, nitems, 
     &        'double', iotype)
        
         call closefile( irstin, "append" )
c         call Write_boundaryflux(myrank,lstep,nshg,ndof,rtmp(:,1:ndof))

c     wallss vectors into the restart file(s)
         if( iowflux .eq. 1) then
            call openfile(  fname2,  'append?', irstin )
            
            fnamer = 'wall shear stresses'        
            isize = nshg*ndof
            
            nitems = 3
            iarray(1) = nshg
            iarray(2) = ndof
            iarray(3) = lstep
            call writeheader(irstin, fnamer,iarray, nitems, isize, 
     &           'double', iotype )
         
c     fnamer = 'boundary flux'        
            nitems = nshg*ndof

         
c     wall shear stresses vectors  are in wallssVec
            call writedatablock(irstin, fnamer,wallssVec, nitems, 
     &           'double', iotype)
            
            call closefile( irstin, "append" )         
         endif! iowflux

c     else
c        print *, "in Bflux: ", myrank, " no printing surface flux"
      endif
c     
      return
c
c.... file error handling
c
997     call error ('bflux   ','opening ', iflux)
c
c$$$1000    format(' ',a80,/,1x,i10,1p,3e20.7)
 2000   format(i6,9(2x,E12.5e2))
c$$$2001    format(1p,1x,i6,3e15.7)
c
c.... end
c
        end


      subroutine flxNode(ienb, iBCB, flg)
c---------------------------------------------------------------------
c
c     This routine flags the flux nodes
c
c----------------------------------------------------------------------
      include "common.h"
c
      integer   flg(nshg),        iBCB(npro,ndiBCB),     
     &          ienb(npro, nshl), lnode(27)

c
c.... compute the nodes which lie on the boundary (hierarchic)
c
      call getbnodes(lnode)

      do i=1, npro 
         if (nsrflist(iBCB(i,2)).eq.1) then
            do j=1, nshlb
               nodelcl = lnode(j)
               flg(abs(ienb(i,nodelcl)))=flg(abs(ienb(i,nodelcl)))+1  
            enddo
         endif
      enddo
c
      return
      end

      
      subroutine mapConn( ien,      ien2,    mask,
     &                    map,      nshl,    npro,    
     &                    npro2,    nshg )
c-----------------------------------------------------------------------
c
c  Create a condensed connectivity array based on the nodes in
c  mask.
c
c-----------------------------------------------------------------------
      
      integer ien(npro,nshl),  ien2(npro,nshl), mask(nshg),
     &        map(npro)

      integer nshl, nshg, npro, npro2, i, iel

c
c.... first build the map
c      
      map = 0
      do i = 1, nshl
         do iel = 1, npro
            map(iel) = map(iel) + mask( abs(ien(iel,i)) )
         enddo
      enddo
      
      npro2 = 0
      do iel = 1, npro
         if ( map(iel) .gt. 0 ) then
            npro2 = npro2 + 1
            map(iel) = npro2
         else
            map(iel) = npro
         endif
      enddo
c
c.... create the condensed connectivity array
c
      if ( npro2 .gt. 0 ) then
         do i = 1, nshl
            do iel = 1, npro
               ien2(map(iel),i) = ien(iel,i)
            enddo
         enddo
      endif
      
      return 
      end
         
      
      subroutine mapArray( x,      x2,      map,
     &                     nshl,   nprold)
c-----------------------------------------------------------------------
c
c  Maps array x into array x2 based on the given map
c
c-----------------------------------------------------------------------
      real*8   x(nprold,nshl),    x2(nprold,nshl)
      
      integer  map(nprold)

      integer   nprold, nshl,  i
      
c
c.... map the array
c
      do i = 1, nshl
         x2(map(:),i) = x(:,i)
      enddo

      return 
      end
