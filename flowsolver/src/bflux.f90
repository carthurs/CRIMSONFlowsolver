      subroutine Bflux ( y,          ac,        u,      &
                         x,          xdist,     xdnv,   &
                         shp,        shgl,      shpb,   &
                         shglb,      ilwork,    iBC,    &
                         BC,         iper  )
!
!----------------------------------------------------------------------
!
! This routine :
!   1. computes the boundary fluxes
!   2. prints the results in the file [FLUX.lstep]
!
! output:
!  in file flux.<lstep>.n (similar to restart file):
!     machin  nshg  lstep 
!     normal_1 ... normal_nsd            ! outward normal direction
!     tau_1n   ... tau_nsd n             ! boundary viscous flux
!
! Zdenek Johan, Summer 1991.
!----------------------------------------------------------------------
!
      
      use pointer_data
      use deformableWall
      use LagrangeMultipliers 
      
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"

      character*5  cname
      
      real*8    y(nshg,ndof),             ac(nshg,ndof), &
                u(nshg,nsd),              &
                x(numnp,nsd), &
                xdist(nshg), &
                xdnv(nshg,nsd)
      dimension iBC(nshg),           &
                BC(nshg,ndofBC),     &
                iper(nshg)
     
      real*8    shp(MAXTOP,maxsh,MAXQPT),  &
                shgl(MAXTOP,nsd,maxsh,MAXQPT), &
                shpb(MAXTOP,maxsh,MAXQPT),  &
                shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!    
!
      real*8    flxres(nshg,nflow), &
                flxLHS(nshg,1),           flxnrm(nshg,nsd), &
                temp(nshg),               rtmp(nshg,ndof), &
                flxTot(nflow),            wallssVec(nshg,ndof)

      real*8    qres(nshg,nsd*nsd)

!
      integer   ilwork(nlwork), &
                invflx(nshg),             nodflx(nshg)             
!
      character*20 fname1,  fmt1, fmt2, fnamer
      character*25 fname2
      !integer irstin, isize, nitems
      integer isize, nitems
      integer iarray(50)  ! integers read from headers

      real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC
      integer, allocatable, dimension(:,:)    :: ien2
      integer, allocatable, dimension(:)      :: map
      real*8, allocatable, dimension(:,:)     :: xmu2
!
!....  calculate the flux nodes
!
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
!
!      do i = 1, nshg
      do i = 1, numnp
         if (invflx(i) .ne. 0) then
            numflx = numflx + 1
            nodflx(numflx) = i
         endif
      enddo
!     
!.... -------------------->   interior elements   <--------------------
!     
!.... initialize the arrays
!     
      flxres = zero
      flxLHS = zero
      flxnrm = zero

      if (numflx .ne. 0)  then !we have flux nodes
         qres   = zero
!     
!.... loop over the element-blocks
!
         lhs    = 0

         ires=2  ! shield e3ql from an unmapped lmassinv
         ierrcalcsave=ierrcalc
         ierrcalc=0
         do iblk = 1, nelblk
!     
!.... set up the parameters
!     
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
!
!.... get the elements touching the boundary
!         
            call mapConn( mien(iblk)%p,    ien2,    invflx, &
                          map,             nshl,    npro,   &
                          npro2,           nshg )

            nprold = npro
            npro = npro2
         
            if (npro .ne. 0) then

               call mapArray( mxmudmi(iblk)%p, xmu2,    map, &
                              maxsh,           nprold)
!
!.... allocate the element matrices (though they're not needed)
!
               allocate ( xKebe(npro,9,nshl,nshl) )
               allocate ( xGoC (npro,4,nshl,nshl) )
               if(Lagrange.gt.zero) then
                  allocate(loclhsLag(npro,9,nshlb,nshlb,3))
               endif 
!     
!.... compute and assemble the residuals
!     
               call AsIGMR (y,                    ac, &
                            x,                    xmu2(1:npro,:), &
                            shp(lcsyst,1:nshl,:), &
                            shgl(lcsyst,:,1:nshl,:), &
                            ien2(1:npro,:),          &
                            flxres,               qres, &
                            xKebe,                xGoC, &
                            rtmp)
!     
               deallocate ( xKebe )
               deallocate ( xGoC  )
               if(Lagrange.gt.zero) then
                  deallocate(loclhsLag)
               endif
            endif
            deallocate ( ien2  )
            deallocate ( xmu2  )
            deallocate ( map   )
!     
         enddo ! iblk = 1, nelblk
         ierrcalc=ierrcalcsave
!     
!.... -------------------->   boundary elements   <--------------------
!     
         do iblk = 1, nelblb
!     
!.... set up the parameters
!
            iel    = lcblkb(1,iblk)
            lcsyst = lcblkb(3,iblk)
            nenl   = lcblkb(5,iblk)
            nshl   = lcblkb(9,iblk)
            nenbl  = lcblkb(6,iblk)
            nshlb  = lcblkb(10,iblk)
            npro   = lcblkb(1,iblk+1) - iel 
 
            if(lcsyst.eq.3) lcsyst=nenbl
!     
            if(lcsyst.eq.3 .or. lcsyst.eq.4) then
               ngaussb = nintb(lcsyst)
            else
               ngaussb = nintb(lcsyst)
            endif

            icurrentblk = iblk  ! current block
!
!.... allocate the element matrices (though they're not needed)
!
            allocate ( xKebe(npro,9,nshl,nshl) )   
!.... compute and assemble the residuals
!
            call AsBFlx (u,                       y, &
                         ac,                         &
                         x,                          &
                         xdist,                      &
                         xdnv,                       &
                         shpb(lcsyst,1:nshl,:),      &
                         shglb(lcsyst,:,1:nshl,:),   &
                         mienb(iblk)%p,              &
                         miBCB(iblk)%p,           mBCB(iblk)%p, &
                         invflx,                  flxres,       &
                         flxLHS,                  flxnrm,       &
                         xKebe,                                 &
                         mSWB(iblk)%p )
!     
            deallocate ( xKebe )
!     
!.... end of boundary element loop
!
         enddo !iblk = 1, nelblb

      else
!         print *, "in Bflux: partition ", myrank, " has no flux nodes!"
      endif  ! make sure the zero numflux processors still commu

!.... Communication needed before we take care of periodicity and
!     division of RHS by LHS ???
!pf: note that the domains that do not have flux nodes,
!    have zero flxres, flxLHS, and flxnrm vectors
!
      if ( numpe > 1 ) then
         call commu (flxres, ilwork, nflow, 'in ')
         call commu (flxLHS, ilwork, 1   , 'in ')
         call commu (flxnrm, ilwork, nsd , 'in ')
      endif
!
!  take care of periodic boundary conditions
!pf: check this!
!
      do j= 1,nshg
         if ((btest(iBC(j),10))) then
            i = iper(j)
            flxLHS(i,1) = flxLHS(i,1) + flxLHS(j,1)
            flxres(i,:) =  flxres(i,:) + flxres(j,:)
         endif
      enddo
!
      do j= 1,nshg
         if ((btest(iBC(j),10))) then
            i = iper(j)
            flxLHS(j,1) = flxLHS(i,1)
            flxres(j,:) = flxres(i,:)
         endif
      enddo
!
!        call bc3per(iBC,  flxres, iper, ilwork, nflow)

!
!.... integrated fluxes (aerodynamic forces update)
!
      flxTot = zero
      do n = 1, numflx
         flxTot = flxTot + flxres(nodflx(n),:)
      enddo
      Force(1) = flxTot(1)
      Force(2) = flxTot(2)
      Force(3) = flxTot(3)
!
!.... only need to commu if we are going to print surface flux since
!     the force calculation just sums flxres (and each on processor node
!     has his "piece" of the sum already).
!
      ntoutv=ntout
      if ( (irs .ge. 1) &
           .and. (mod(lstep, ntoutv) .eq. 0) &
           .or.  (istep .eq. nstep(itseq)) ) then

!
!  need to zero the slaves to prevent counting twice
!  (actually unnecessary since flxres of boundary nodes will be counted n
!  times while flxlhs will be counted n times-> the ratio is still
!  correct
!      
         wallssVec=rtmp

         if (numflx .eq. 0) then   !no flux nodes
            rtmp=zero
            wallssVec = zero
         else
!     
!.... ---------------------------->  Solve  <---------------------------
!
!.... compute the viscous and heat fluxes
!     
!
!.... ---------------------------->  Print  <---------------------------
!
!.... nodal fluxes
!
            do i = 1, 3
               where ( (invflx .ne. 0) .and. (flxLHS(:,1) .ne. zero) ) &
                    flxres(:,i) = flxres(:,i) / flxLHS(:,1)
            enddo
!     
!.... normalize the outward normal
!     
            temp = sqrt( flxnrm(:,1)**2 &
                       + flxnrm(:,2)**2 &
                       + flxnrm(:,3)**2 )
            where ( (invflx .ne. 0) .and. (temp .ne. zero) )
               flxnrm(:,1) = flxnrm(:,1) / temp
               flxnrm(:,2) = flxnrm(:,2) / temp
               flxnrm(:,3) = flxnrm(:,3) / temp
            endwhere
         endif !no flux nodes
         
!     
!.... ---------------------------->  Communications <-------------------
!
         if(numpe > 1) then
!           print *, "in Bflux: ", myrank, 
!    &               " is calling commu (3 times)"
            call commu (flxres, ilwork, nflow, 'out')
            call commu (flxLHS, ilwork, 1   , 'out')
            call commu (flxnrm, ilwork, nsd , 'out')
         endif
!
         rtmp = zero
         wallssVec  = zero

         do i=1, numnp
            if (invflx(i) .ne. 0) then
               rtmp(i,2:4) = flxres(i,1:3) !viscous flux
!     calculate the WSS
               tn = flxres(i,1) * flxnrm(i,1) &
                  + flxres(i,2) * flxnrm(i,2) &
                  + flxres(i,3) * flxnrm(i,3)

                wallssVec(i,1) = flxres(i,1) - tn * flxnrm(i,1)
                wallssVec(i,2) = flxres(i,2) - tn * flxnrm(i,2)
                wallssVec(i,3) = flxres(i,3) - tn * flxnrm(i,3)
            endif
         enddo

!         itmp = 1
!         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
!         write (fmt1,"('(''flux.'',i',i1,',1x)')") itmp
!         write (fname1,fmt1) lstep
      
!         fname1 = trim(fname1) // cname(myrank+1)
   
!        print *, "in Bflux: ",myrank," is opening flux file", fname1

!         open (unit=iflux, file=fname1, status='unknown', 
!     &         form='formatted',err=997)

!      write (iflux) machin, nshg, lstep
!      write (iflux) rtmp(:,1:6)
!
!.... output the results
!     
!         do n = 1, numflx
!            k = nodflx(n)
!            write (iflux,2000) k, (x(k,i), i=1,3), 
!     &           (flxnrm(k,i),  i=1,3),
!     &           (flxres(k,i),  i=1,3)
!         enddo
!         close (iflux)

!        print *, "in Bflux: ",myrank," closed flux file", fname1

!... output the results in the new format in restart.step#.proc# file

         itmp = 1
         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
         write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
         write (fname2,fmt2) lstep

         fname2 = trim(fname2) // cname(myrank+1)
!
!.... open input files
!
         call openfile(  fname2//c_null_char,  c_char_"append?"//c_null_char, irstin )
         
         fnamer = 'boundary flux'        
         isize = nshg*ndof
         nitems = 3;
         iarray(1) = nshg
         iarray(2) = ndof
         iarray(3) = lstep
         call writeheader(irstin, fnamer//c_null_char,iarray, nitems, isize, &
              c_char_"double"//c_null_char, iotype )
    
!         fnamer = 'boundary flux'        
         nitems = nshg*ndof
         call writedatablock(irstin, fnamer//c_null_char,rtmp, nitems, &
              c_char_"double"//c_null_char, iotype)
        
         call closefile( irstin, c_char_"append"//c_null_char )
!         call Write_boundaryflux(myrank,lstep,nshg,ndof,rtmp(:,1:ndof))

!     wallss vectors into the restart file(s)
         if( iowflux .eq. 1) then
            call openfile(  fname2//c_null_char,  c_char_"append?"//c_null_char, irstin )
            
            fnamer = 'wall shear stresses'        
            isize = nshg*ndof
            
            nitems = 3
            iarray(1) = nshg
            iarray(2) = ndof
            iarray(3) = lstep
            call writeheader(irstin, fnamer//c_null_char,iarray, nitems, isize, &
                 c_char_"double"//c_null_char, iotype )
         
!     fnamer = 'boundary flux'        
            nitems = nshg*ndof

         
!     wall shear stresses vectors  are in wallssVec
            call writedatablock(irstin, fnamer//c_null_char,wallssVec, nitems, &
                 c_char_"double"//c_null_char, iotype)
            
            call closefile( irstin, c_char_"append"//c_null_char )
         endif! iowflux

!     else
!        print *, "in Bflux: ", myrank, " no printing surface flux"
      endif
!     
      return
!
!.... file error handling
!
997     call error ('bflux   ','opening ', iflux)
!
!$$$1000    format(' ',a80,/,1x,i10,1p,3e20.7)
 2000   format(i6,9(2x,E12.5e2))
!$$$2001    format(1p,1x,i6,3e15.7)
!
!.... end
!
        end


      subroutine flxNode(ienb, iBCB, flg)
!---------------------------------------------------------------------
!
!     This routine flags the flux nodes
!
!----------------------------------------------------------------------
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      integer   flg(nshg),        iBCB(npro,ndiBCB), &
                ienb(npro, nshl), lnode(27)

!
!.... compute the nodes which lie on the boundary (hierarchic)
!
      call getbnodes(lnode)

      do i=1, npro 
         if (nsrflist(iBCB(i,2)).eq.1) then
            do j=1, nshlb
               nodelcl = lnode(j)
               flg(abs(ienb(i,nodelcl)))=flg(abs(ienb(i,nodelcl)))+1  
            enddo
         endif
      enddo
!
      return
      end

      
      subroutine mapConn( ien,      ien2,    mask, &
                          map,      nshl,    npro, &  
                          npro2,    nshg )
!-----------------------------------------------------------------------
!
!  Create a condensed connectivity array based on the nodes in
!  mask.
!
!-----------------------------------------------------------------------
      
      integer ien(npro,nshl),  ien2(npro,nshl), mask(nshg), &
              map(npro)

      integer nshl, nshg, npro, npro2, i, iel

!
!.... first build the map
!      
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
!
!.... create the condensed connectivity array
!
      if ( npro2 .gt. 0 ) then
         do i = 1, nshl
            do iel = 1, npro
               ien2(map(iel),i) = ien(iel,i)
            enddo
         enddo
      endif
      
      return 
      end
         
      
      subroutine mapArray( x,      x2,      map, &
                           nshl,   nprold)
!-----------------------------------------------------------------------
!
!  Maps array x into array x2 based on the given map
!
!-----------------------------------------------------------------------
      real*8   x(nprold,nshl),    x2(nprold,nshl)
      
      integer  map(nprold)

      integer   nprold, nshl,  i
      
!
!.... map the array
!
      do i = 1, nshl
         x2(map(:),i) = x(:,i)
      enddo

      return 
      end
