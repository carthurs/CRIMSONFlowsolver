!  readnblk.f (pronounce "Reed and Block Dot Eff") contains:
!
!    module readarrays ("Red Arrays") -- contains the arrays that
!     are read in from binary files but not immediately blocked 
!     through pointers.
!
!    subroutine readnblk ("Reed and Block") -- allocates space for
!     and reads data to be contained in module readarrays.  Reads
!     all remaining data and blocks them with pointers.
!


      module readarrays
      
      use, intrinsic :: iso_c_binding
      !real*8, allocatable :: x(:,:)
      !real*8, allocatable :: qold(:,:)
      !real*8, allocatable :: uold(:,:)
      !real*8, allocatable :: acold(:,:)
      integer, allocatable :: iBCtmp(:)
      real*8, allocatable, target :: BCinp(:,:)

      !integer, allocatable :: ilwork(:)
      integer, allocatable, target :: nBC(:)
      !integer, allocatable :: iper(:)
      !integer, allocatable :: ifath(:)
      !integer, allocatable :: nsons(:)
      
      end module



      subroutine readnblk
!     
      use, intrinsic :: iso_c_binding
      
      use readarrays
      use globalArrays
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      real (c_double), allocatable :: xread(:,:), qread(:,:), acread(:,:)
      integer, allocatable :: ntagread(:)
      real (c_double), allocatable :: uread(:,:)
      real (c_double), allocatable :: BCinpread(:,:)
      integer, allocatable :: iperread(:), iBCtmpread(:)
      integer, allocatable :: ilworkread(:), nBCread(:)
      character(c_char) :: cname*5
      character(c_char) :: fmt1r1*20
      character(c_char) :: fname1*255,fnameg*255,fnamer*255,fnamelr*255,actname*255
      character*255 warning
      !integer igeom, ibndc, irstin, 
      integer ierr, ixsiz
      integer intfromfile(50) ! integers read from headers
       
!
!
!.... determine the step number to start with
!
      open(unit=72,file='numstart.dat',status='old')
      read(72,*) irstart
      close(72)
!
      fnameg='geombc.dat'
      fnameg= trim(fnameg)  // cname(myrank+1)

      itmp=1
      if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
      write (fmt1r1,"('(''restart.'',i',i1,',1x)')") itmp
      write (fnamer,fmt1r1) irstart
      fnamer = trim(fnamer) // cname(myrank+1)
      
      !fnamelr='restart.latest'
      !fnamelr = trim(fnamelr) // cname(myrank+1)

!
!.... open input files
!
      call openfile( fnameg//c_null_char, c_char_"read?"//c_null_char, igeom );
!
!.... try opening restart.latest.proc before trying restart.stepno.proc
!
      !call openfile(  fnamelr,  'read?', irstin );
      !if ( irstin .eq. 0 ) call openfile( fnamer, 'read?', irstin );
      call openfile( fnamer//c_null_char, c_char_"read?"//c_null_char, irstin );
! either one will work
!
!.... input the geometry parameters
!

      ieleven=11
      ione=1
      !fname1='number of spatial dimensions?'
      !call readheader(igeom,fname1,nsd,ione,c_char_"integer", iotype)
      !if(nsd.ne.2) nsd=3        ! in case it is an old geombc file 
      !                          ! that does not have nsd in it
                                
      fname1='number of nodes?'
      call readheader(igeom,fname1//c_null_char,numnp,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='number of modes?'
      call readheader(igeom,fname1//c_null_char,nshg,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='number of interior elements?'
      call readheader(igeom,fname1//c_null_char,numel,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='number of boundary elements?'
      call readheader(igeom,fname1//c_null_char,numelb,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='maximum number of element nodes?'
      call readheader(igeom,fname1//c_null_char,nen,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='number of interior tpblocks?'
      call readheader(igeom,fname1//c_null_char,nelblk,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='number of boundary tpblocks?'
      call readheader(igeom,fname1//c_null_char,nelblb,ione,c_char_"integer"//c_null_char, iotype)
                              
      fname1='number of nodes with Dirichlet BCs?'
      call readheader(igeom,fname1//c_null_char,numpbc,ione,c_char_"integer"//c_null_char, iotype)
           
      fname1='number of shape functions?'
      call readheader(igeom,fname1//c_null_char,ntopsh,ione,c_char_"integer"//c_null_char, iotype)

      fname1='number of boundary element tag IDs?'
      call readheader(igeom,fname1//c_null_char,numBETFields,ione,c_char_"integer"//c_null_char, iotype)

!
!.... calculate the maximum number of boundary element nodes
!     
      nenb = 0
      do i = 1, melCat
         if (nen .eq. nenCat(i,nsd)) nenb = max(nenCat(i,nsd-1), nenb)
      enddo
!     
      if (myrank == master) then
         if (nenb .eq. 0) call error ('input   ','nen     ',nen)
      endif
!
!.... setup some useful constants
!
      I3nsd  = nsd / 3          ! nsd=3 integer flag
      E3nsd  = float(I3nsd)     ! nsd=3 real    flag
!    
      if(matflg(1,1).lt.0) then
         nflow = nsd + 1
      else
         nflow = nsd + 2
      endif 
      ndof   = nsd + 2
      nsclr=impl(1)/100
      ndof=ndof+nsclr           ! number of sclr transport equations to solve
      
      ndofBC = ndof + I3nsd     ! dimension of BC array
      ndiBCB = 2                ! dimension of iBCB array
      ndBCB  = ndof + 1         ! dimension of BCB array
!     
      nsymdf = (ndof*(ndof + 1)) / 2 ! symm. d.o.f.'s
!
! now that we have all of the constants set, initialize all of the
! arrays
      call initGlobalArrays
      
!
!.... ----------------------> Communication tasks <--------------------
!
      if(numpe > 1) then

         fname1='size of ilwork array?'
         call readheader(igeom,fname1//c_null_char,nlwork,ione,c_char_"integer"//c_null_char, iotype)

         ione=1
         fname1='ilwork?'
         call readheader(igeom,fname1//c_null_char,nlwork,ione,c_char_"integer"//c_null_char, iotype)

         if (.not.ALLOCATED(ilwork)) then
           allocate( ilwork(nlwork) )
         endif
         if (.not.ALLOCATED(ilworkread)) then
           allocate( ilworkread(nlwork) )
         endif
         call readdatablock(igeom,fname1//c_null_char,ilworkread, &
                            nlwork,c_char_"integer"//c_null_char, iotype)
         ilwork = ilworkread
         call ctypes (ilwork)
      else
           nlwork=1
           if (.not.ALLOCATED(ilwork)) then
            allocate( ilwork(1)) 
           endif
      endif

!     
!.... read the node coordinates
!
      itwo=2
      fname1='co-ordinates?'
      call readheader(igeom,fname1//c_null_char,intfromfile,itwo,c_char_"double"//c_null_char, iotype)
      numnp=intfromfile(1)
!      nsd=intfromfile(2)
!      allocate( x(numnp,nsd) )
      if (.not.allocated(xread)) then
        allocate( xread(numnp,nsd) )
      endif
      ixsiz=numnp*nsd
      call readdatablock(igeom,fname1//c_null_char,xread,ixsiz, c_char_"double"//c_null_char,iotype)
      x = xread

!
!.... read the node tags
!
      if (geombcHasNodeTags) then
        itwo=2
        fname1='node tags?'
        call readheader(igeom,fname1//c_null_char,intfromfile,itwo,c_char_"integer"//c_null_char, iotype)
        numnp=intfromfile(1)
        if (.not.allocated(ntagread)) then
          allocate( ntagread(numnp) )
        endif
        ixsiz=numnp
        call readdatablock(igeom,fname1//c_null_char,ntagread,ixsiz, c_char_"integer"//c_null_char,iotype)
        nodetagfield = ntagread
      endif

!
!.... read in and block out the connectivity
!  
      call genblk (IBKSIZ)      
!
!.... read the boundary condition mapping array
!
      ione=1
      fname1='bc mapping array?'
      call readheader(igeom,fname1//c_null_char,nshg, &
           ione,c_char_"integer"//c_null_char, iotype)
      if (.not.ALLOCATED(nBC)) then
       allocate( nBC(nshg) )
      endif

      if (.not.allocated(nBCread)) then
        allocate( nBCread(nshg) )
      endif
      call readdatablock(igeom,fname1//c_null_char,nBCread,nshg,c_char_"integer"//c_null_char,iotype)
      nBC=nBCread


!
!.... read the temporary iBC array
!
      ione = 1
      fname1='bc codes array?'
      call readheader(igeom,fname1//c_null_char,numpbc, &
           ione, c_char_"integer"//c_null_char, iotype)
      if ( numpbc > 0 ) then
         if (.not.allocated(iBCtmp)) then
           allocate( iBCtmp(numpbc) )
         endif

         if (.not.allocated(iBCtmpread)) then
           allocate( iBCtmpread(numpbc) )
         endif
         call readdatablock(igeom,fname1//c_null_char,iBCtmpread,numpbc, &
                            c_char_"integer"//c_null_char,iotype)
         iBCtmp=iBCtmpread

      else  ! sometimes a partition has no BC's
         if (.not.allocated(iBCtmp)) then
           allocate( iBCtmp(1) )
         endif
         iBCtmp=0

      endif
!
!.... read boundary condition data
!
      ione=1
      fname1='boundary condition array?'
      call readheader(igeom,fname1//c_null_char,intfromfile, &
           ione, c_char_"integer"//c_null_char, iotype)
! here intfromfile(1) contains (ndof+7)*numpbc
      if ( numpbc > 0 ) then
         if(intfromfile(1).ne.(ndof+7)*numpbc) then
           warning='WARNING more data in BCinp than needed: keeping 1st'
           write(*,*) warning, ndof+7
         endif
         if (.not.allocated(BCinp)) then
           allocate( BCinp(numpbc,ndof+7) )
         endif
         nsecondrank=intfromfile(1)/numpbc
         if (.not.allocated(BCinpread)) then
           allocate( BCinpread(numpbc,nsecondrank) )
         endif
         iBCinpsiz=intfromfile(1)
         call readdatablock(igeom,fname1//c_null_char,BCinpread,iBCinpsiz, &
                            c_char_"double"//c_null_char,iotype)
         BCinp(:,1:(ndof+7))=BCinpread(:,1:(ndof+7))

      else  ! sometimes a partition has no BC's
         if (.not.allocated(BCinp)) then
           allocate( BCinp(1,ndof+7) )
         endif
         BCinp=0

      endif

!
!.... read periodic boundary conditions
!
      fname1='periodic masters array?'
      call readheader(igeom,fname1//c_null_char,nshg, &
           ione, c_char_"integer"//c_null_char, iotype)
      !allocate( iper(nshg) )
      if (.not.allocated(iperread)) then
        allocate( iperread(nshg) )
      endif
      call readdatablock(igeom,fname1//c_null_char,iperread,nshg, &
                            c_char_"integer"//c_null_char,iotype)
      iper=iperread


      !
      !.... read in the local index of unique nodes
      !
      if(numpe > 1) then
          fname1='local index of unique nodes?'
          call readheader(igeom,fname1//c_null_char,nshguniq, &
          ione, c_char_"integer"//c_null_char, iotype)
          if (.not.allocated(inodesuniq)) then
            allocate( inodesuniq(nshguniq) )
          endif
          call readdatablock(igeom,fname1//c_null_char,inodesuniq,nshguniq, &
          c_char_"integer"//c_null_char,iotype)

      else
          ! if (.not.allocated(inodesuniq)) then
          !   allocate( inodesuniq(nshg) )
          ! endif
          if (allocated(inodesuniq)) then
            deallocate(inodesuniq)
            allocate( inodesuniq(nshg) )
          else
            allocate( inodesuniq(nshg) )
          endif
          nshguniq = nshg
          do ii=1,nshg
              inodesuniq(ii) = ii
          end do

      endif


      !
      !.... read in the simple observation function arrays
      !
      if (geombcHasObservaionFields) then
        itwo=2
        fname1='observation function solution?'
        call readheader(igeom,fname1//c_null_char,intfromfile, &
        itwo,c_char_"integer"//c_null_char, iotype)
        nshg2=intfromfile(1)
        ndof2=intfromfile(2)
        iisiz=nshg2*ndof2
        if (.not.allocated(ilinobsfunc_sol)) then
          allocate( ilinobsfunc_sol(nshg2,ndof2) )
        endif
        call readdatablock(igeom,fname1//c_null_char,ilinobsfunc_sol,iisiz, &
        c_char_"integer"//c_null_char,iotype)

        fname1='observation function time derivative of solution?'
        call readheader(igeom,fname1//c_null_char,intfromfile, &
        itwo,c_char_"integer"//c_null_char, iotype)
        nshg2=intfromfile(1)
        ndof2=intfromfile(2)
        iisiz=nshg2*ndof2
        if (.not.allocated(ilinobsfunc_acc)) then
          allocate( ilinobsfunc_acc(nshg2,ndof2) )
        endif
        call readdatablock(igeom,fname1//c_null_char,ilinobsfunc_acc,iisiz, &
        c_char_"integer"//c_null_char,iotype)

        if (ideformwall.eq.1) then
            fname1='observation function displacement?'
            call readheader(igeom,fname1//c_null_char,intfromfile, &
            itwo,c_char_"integer"//c_null_char, iotype)
            nshg2=intfromfile(1)
            nsd2=intfromfile(2)
            iisiz=nshg2*nsd2
            if (.not.allocated(ilinobsfunc_disp)) then
              allocate( ilinobsfunc_disp(nshg2,nsd2) )
            endif
            call readdatablock(igeom,fname1//c_null_char,ilinobsfunc_disp,iisiz, &
            c_char_"integer"//c_null_char,iotype)

            fname1='observation function distance?'
            call readheader(igeom,fname1//c_null_char,intfromfile, &
            itwo,c_char_"integer"//c_null_char, iotype)
            nshg2=intfromfile(1)
            iisiz=nshg2
            if (.not.allocated(obsfunc_dist)) then
              allocate( obsfunc_dist(nshg2) )
            endif
            call readdatablock(igeom,fname1//c_null_char,obsfunc_dist,iisiz, &
            c_char_"integer"//c_null_char,iotype)

        endif
      endif

!
!.... generate the boundary element blocks
!
      call genbkb (ibksiz)


!
!  Read in the nsons and ifath arrays if needed
!
!  There is a fundamental shift in the meaning of ifath based on whether
!  there exist homogenous directions in the flow.  
!
!  HOMOGENOUS DIRECTIONS EXIST:  Here nfath is the number of inhomogenous
!  points in the TOTAL mesh.  That is to say that each partition keeps a 
!  link to  ALL inhomogenous points.  This link is furthermore not to the
!  sms numbering but to the original structured grid numbering.  These 
!  inhomogenous points are thought of as fathers, with their sons being all
!  the points in the homogenous directions that have this father's 
!  inhomogeneity.  The array ifath takes as an arguement the sms numbering
!  and returns as a result the father.
!
!  In this case nsons is the number of sons that each father has and ifath
!  is an array which tells the 
!
!  NO HOMOGENOUS DIRECTIONS.  In this case the mesh would grow to rapidly
!  if we followed the above strategy since every partition would index its
!  points to the ENTIRE mesh.  Furthermore, there would never be a need
!  to average to a node off processor since there is no spatial averaging.
!  Therefore, to properly account for this case we must recognize it and
!  inerrupt certain actions (i.e. assembly of the average across partitions).
!  This case is easily identified by noting that maxval(nsons) =1 (i.e. no
!  father has any sons).  Reiterating to be clear, in this case ifath does
!  not point to a global numbering but instead just points to itself.
!
      nfath=1  ! some architectures choke on a zero or undeclared
                 ! dimension variable.  This sets it to a safe, small value.
!      if(((iLES .lt. 20) .and. (iLES.gt.0)) &
!                         .or. (itwmod.gt.0)  ) then ! don't forget same
!                                                    ! conditional in proces.f
!
!!           read (igeom) nfath  ! nfath already read in input.f,
!                                     ! needed for alloc
!         ione=1
!!         call creadlist(igeom,ione,nfath)
!!         fname1='keyword sonfath?'
!         if(nohomog.gt.0) then
!            fname1='number of father-nodes?'
!            call readheader(igeom,fname1//c_null_char,nfath,ione,c_char_"integer"//c_null_char, iotype)
!!
!!     fname1='keyword nsons?'
!            fname1='number of son-nodes for each father?'
!            call readheader(igeom,fname1//c_null_char,nfath,ione,c_char_"integer"//c_null_char, iotype)
!            allocate (nsons(nfath))
!            call readdatablock(igeom,fname1//c_null_char,nsons,nfath, &
!                            c_char_"integer"//c_null_char,iotype)
!!
!            fname1='keyword ifath?'
!            call readheader(igeom,fname1//c_null_char,nshg,ione,c_char_"integer"//c_null_char, iotype)
!            allocate (ifath(nshg))
!            call readdatablock(igeom,fname1//c_null_char,ifath,nshg, &
!                            c_char_"integer"//c_null_char,iotype)
!!
!            nsonmax=maxval(nsons)
!!
!         else  ! this is the case where there is no homogeneity
!               ! therefore ever node is a father (too itself).  sonfath
!               ! (a routine in NSpre) will set this up but this gives
!               ! you an option to avoid that.
!            nfath=nshg
!            allocate (nsons(nfath))
!            nsons=1
!            allocate (ifath(nshg))
!            do i=1,nshg
!               ifath(i)=i
!            enddo
!            nsonmax=1
!!
!         endif
!      else
!         allocate (nsons(1))
!         allocate (ifath(1))
!      endif
      if (allocated(velbar)) then
        deallocate (velbar)
      endif
      allocate (velbar(nfath,ndof))


!
!  renumber the master partition for SPEBC
!
!      if((myrank.eq.master).and.(irscale.ge.0)) then
!         call setSPEBC(numnp, nfath, nsonmax)
!         call renum(x,ifath,nsons)
!      endif
!
!.... Read restart files
!
!.... read the header and check it against the run data
!
      ithree=3
!      call creadlist(irstin,ithree,nshg2,ndof2,lstep)
      fname1='solution?'
      call readheader(irstin,fname1//c_null_char,intfromfile, &
           ithree,c_char_"integer"//c_null_char, iotype)
      nshg2=intfromfile(1)
      ndof2=intfromfile(2)
      lstep=intfromfile(3)
      if(ndof2.ne.ndof) then
        warning='WARNING more data in restart than needed: keeping 1st '
        write(*,*) warning , ndof
      endif
!
      if (nshg2 .ne. nshg) call error ('restar  ', 'nshg   ', nshg)
      
!
!.... read the values of primitive variables into q
!
      !allocate( qold(nshg,ndof) )
      if (allocated(qread)) then
        deallocate( qread )
      endif
      allocate( qread(nshg,ndof2) )

      iqsiz=nshg*ndof2
      call readdatablock(irstin,fname1//c_null_char,qread,iqsiz, &
                            c_char_"double"//c_null_char,iotype)
      yold(:,1:3)=qread(:,2:4)
      yold(:,4)=qread(:,1)
      if(ndof.gt.4) then
          yold(:,5:ndof)=qread(:,5:ndof)
      else
          yold(:,5:ndof)=zero
      endif
      y = yold ! initialization of y moved here from genini
! 
      fname1='time derivative of solution?'
      intfromfile=0
      call readheader(irstin,fname1//c_null_char,intfromfile, &
           ithree,c_char_"integer"//c_null_char, iotype)
      !allocate( acold(nshg,ndof) )

      if(intfromfile(1).ne.0) then 
         nshg2=intfromfile(1)
         ndof2=intfromfile(2)
         lstep=intfromfile(3)
         
         if (nshg2 .ne. nshg) call error ('restar  ', 'nshg   ', nshg)

         if (.not.allocated(acread)) then
           allocate( acread(nshg,ndof2) )
         endif
         acread=zero

         iacsiz=nshg*ndof2
         call readdatablock(irstin,fname1//c_null_char,acread,iacsiz, &
                         c_char_"double"//c_null_char,iotype)
         acold(:,1:3)=acread(:,2:4)
         acold(:,4)=acread(:,1)
         ac = acold ! initialization of ac moved here from genini
         if(ndof.gt.4)  acold(:,5:ndof)=acread(:,5:ndof)
         if (allocated(acread)) then
           deallocate(acread)
         endif
      else
         warning='Time derivative of solution is set to zero (SAFE)'
         write(*,*) warning
         acold=zero
         ac=zero
      endif

!
!....  read the displacement
!

      if (ideformwall.eq.1) then
         fname1='displacement?'
         call readheader(irstin,fname1//c_null_char,intfromfile, &
              ithree,c_char_"integer"//c_null_char, iotype)
         nshg2=intfromfile(1)
         ndisp=intfromfile(2)
         lstep=intfromfile(3)
         if(ndisp.ne.nsd) then
            warning='WARNING ndisp not equal nsd'
            write(*,*) warning , ndisp
         endif
!
         if (nshg2 .ne. nshg) call error ('restar  ', 'nshg   ', nshg)

!
!.... read the values of primitive variables into uold
!
         !allocate( uold(nshg,nsd) )
         if (.not.allocated(uread)) then
           allocate( uread(nshg,nsd) )
         endif
         
         iusiz=nshg*nsd
         
         call readdatablock(irstin,fname1//c_null_char,uread,iusiz, &
              c_char_"double"//c_null_char,iotype)
         uold(:,1:nsd)=uread(:,1:nsd)
         u = uold


         !
         !....  read the reference displacement
         !
         if (iinitialprestress.eq.1) then
             fname1='displacement_ref?'
             call readheader(irstin,fname1//c_null_char,intfromfile, &
                 ithree,c_char_"integer"//c_null_char, iotype)
             nshg2=intfromfile(1)
             ndisp=intfromfile(2)
             lstep=intfromfile(3)
             if(ndisp.ne.nsd) then
                 warning='WARNING ndisp not equal nsd'
                 write(*,*) warning , ndisp
             endif
             !
             if (nshg2 .ne. nshg) call error ('restar  ', 'nshg   ', nshg)

             iusiz=nshg*nsd

             call readdatablock(irstin,fname1//c_null_char,uread,iusiz, &
                 c_char_"double"//c_null_char,iotype)
             uref(:,1:nsd)=uread(:,1:nsd)
         else
             uref(:,1:nsd)=zero
         endif

      else
         !allocate( uold(nshg,nsd) )
         uold(:,1:nsd) = zero
         uref(:,1:nsd) = zero
      endif

      temporary_array = zero

      call PhAssignPointerInt(c_loc(inodesuniq), c_char_"local index of unique nodes"//c_null_char)
      if (geombcHasObservaionFields) then
        call PhAssignPointerInt(c_loc(ilinobsfunc_sol), c_char_"observation function solution"//c_null_char)
        call PhAssignPointerInt(c_loc(ilinobsfunc_acc), c_char_"observation function time derivative of solution"//c_null_char)
        call PhAssignPointerInt(c_loc(ilinobsfunc_disp), c_char_"observation function displacement"//c_null_char)
        call PhAssignPointerInt(c_loc(obsfunc_dist), c_char_"observation function distance"//c_null_char)
      endif

      call PhAssignPointerDP(c_loc(yold), c_char_"solution"//c_null_char)
      call PhAssignPointerDP(c_loc(acold), c_char_"time derivative of solution"//c_null_char)
      call PhAssignPointerDP(c_loc(uold), c_char_"displacement"//c_null_char)
      call PhAssignPointerDP(c_loc(x), c_char_"coordinates"//c_null_char)
      call PhAssignPointerDP(c_loc(temporary_array), c_char_"temporary array"//c_null_char)
      call PhAssignPointerDP(c_loc(xdist), c_char_"distance"//c_null_char)
      call PhAssignPointerDP(c_loc(df_fem), c_char_"M*distance"//c_null_char)
      call PhAssignPointerDP(c_loc(xdnv), c_char_"distance normal"//c_null_char)

! 
!
!.... close c-binary files
!
      call closefile( irstin, c_char_"read" )
      call closefile( igeom,  c_char_"read" )
!
      if (allocated(xread)) then
        deallocate(xread)
      endif
      if (allocated(qread)) then
        deallocate(qread)
      endif
      if ( numpbc > 0 )  then
         if (allocated(bcinpread)) then
           deallocate(bcinpread)
         endif
         if (allocated(ibctmpread)) then
           deallocate(ibctmpread)
         endif
      endif
      if (allocated(iperread)) then
        deallocate(iperread)
      endif
      if(numpe.gt.1) then
           if (allocated(ilworkread)) then
             deallocate(ilworkread)
           endif
      endif
      if (allocated(nbcread)) then
        deallocate(nbcread)
      endif


      return
!
 994  call error ('input   ','opening ', igeom)
 995  call error ('input   ','opening ', igeom)
 997  call error ('input   ','end file', igeom)
 998  call error ('input   ','end file', igeom)
!
      end

!!
!! No longer called but kept around in case....
!!
!      subroutine genpzero(iBC)
!
!      use pointer_data
!!
!       use phcommonvars
! IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!      integer iBC(nshg)
!!
!!....  check to see if any of the nodes have a dirichlet pressure
!!
!      pzero=1
!      if (any(btest(iBC,2))) pzero=0
!!
!      do iblk = 1, nelblb
!         npro = lcblkb(1,iblk+1)-lcblkb(1,iblk)
!         do i=1, npro
!            iBCB1=miBCB(iblk)%p(i,1)
!!
!!.... check to see if any of the nodes have a Neumann pressure
!!     but not periodic (note that
!!
!            if(btest(iBCB1,1)) pzero=0
!         enddo
!!
!!.... share results with other processors
!!
!         pzl=pzero
!         if (numpe .gt. 1) &
!              call MPI_ALLREDUCE (pzl, pzero, 1, &
!              MPI_DOUBLE_PRECISION,MPI_MIN, INEWCOMM,ierr)
!
!      enddo
!!
!!.... return
!!
!      return
!!
!      end
