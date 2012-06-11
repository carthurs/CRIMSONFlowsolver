      module dtnmod
      integer, allocatable :: ifeature(:)
      end module


      subroutine Ddtnmod
      use dtnmod
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      deallocate (ifeature)
      end

      subroutine initDtN
      use dtnmod
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      allocate (ifeature(nshg))
      end


      subroutine DtN(iBC,BC,y)
      use dtnmod
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      real*8 BC(nshg,ndofBC),y(nshg,ndof),tmp(nsclr)
      integer iBC(nshg)
      do i=1,nshg
         itype=ifeature(i)
         if(btest(iBC(i),13)) then
            do j=1,nsclr
               tmp(j)=y(i,5+j)
            end do
            call Dirichlet2Neumann(nsclr, itype, tmp)
!
!  put the value in the position a Dirichlet value would be in BC.
!  later we will localize this value to the BCB array.  
!  this is not dangerous because we should NEVER need to set Dirichlet
!  on the same node as a DtN condition
!
            do j=1,nsclr
               BC(i,6+j)=-tmp(j)
            end do
         endif
      end do
      return
      end

      subroutine Dirichlet2Neumann_faux(nscalar, itype, tmp)
!
! This is a fake routine, designed to do nothing but feed back
! fluxes as if there were a fast reaction at the surface and the
! thickness of the stagnant BL were given by "distance".
! It can handle up to 24 different scalars.
!
! If itype is zero, the flux is arbitrarily set to half what it would
! be for any other itype.
!
! The assumption of "units" is that the initial concentrations are in
! moles/cubic-meter and the fluxes are in moles/(sec square-meter)
!
! The listed diffusivities are characteristic of metal-ions in room
! temperature water, in units of square-meter/sec.
!
! 
      integer itype, nscalar
      real*8 tmp(nscalar)

      integer i
      real*8 distance
      real*8 units

!  Completely fake diffusivities
      real*8 D(24)
      data D/ & 
             5.0d-05, 1.0d-5, 8.0d-10, 7.0d-10, 6.0d-10, 5.0d-10, & 
             4.0d-10, 3.0d-10, 2.0d-10, 1.0d-10, 0.9d-10, 0.8d-10,  & 
             1.0d-10, 9.0d-10, 8.0d-10, 7.0d-10, 6.0d-10, 5.0d-10, & 
             4.0d-10, 3.0d-10, 2.0d-10, 1.0d-10, 0.9d-10, 0.8d-10/ 

      do i=1,nscalar 
         tmp(i) = 0.0d0
      enddo
      return
      distance = 10.0d-03
      units = 1.0d-3
      units = 1.0
      if(nscalar.gt.24) then
         write(*,*) 'Problem in Dir2Neu: nscalar larger than 24!'
         stop
      endif
      
      do i=1,nscalar 
         tmp(i) = D(i) * ( tmp(i) - 0.0 ) / distance  * units
         if(itype.eq.2) tmp(i) = tmp(i) / 2.0d+00
      enddo
!      tmp(1)=1.0d-5
      return
      end









      subroutine dtnl(iBC,BC,ienb,iBCB,BCB)
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      integer ienb(npro,nshl), iBC(nshg),iBCB(npro,ndiBCB)
      real*8  BCB(npro,nshlb,ndBCB), tmpBCB(npro,nshlb,nsclr), & 
              BC(nshg,ndofBC),        tmpBC(nshg,nsclr)

      nstart=ndofBC-nsclr
!      tmpBC=zero
!      do i=1,nshg
!         if(btest(iBC(i),13)) then
            do j=1,nsclr
!               tmpBC(i,j)=BC(i,nstart+j)
               tmpBC(:,j)=BC(:,nstart+j)
            enddo
!         endif
!      enddo
      
      call localb(tmpBC,tmpBCB,ienb,nsclr,'gather  ')
      
      do i=1,npro
         do j=1,nsclr
            if(iBCB(i,2).lt.0) then  !this is a face with dtn
               do k=1,nshlb
                  BCB(i,k,6+j)=tmpBCB(i,k,j)
               enddo
            endif
         enddo
      enddo
      return
      end

!
! This routine just calls the appropriate version of D2N for the number 
! of scalars used
!
      subroutine Dirichlet2Neumann(nscalar, itype, tmp)
      integer nscalar, itype
      real*8 tmp(nscalar),foo
      
! Just short circuit the routine for a little bit.
!      tmp(1)=0.0d0
!      return
      if(nscalar .eq. 1) then
!         write(*,*) 'Entering D2N1'
!          foo= rand(0)
         call Dirichlet2Neumann_1(nscalar,itype,tmp)
!         write(*,*) 'Returning from D2N after DTN1'
!         return
      elseif(nscalar.eq.2) then
         call Dirichlet2Neumann_2(nscalar,itype,tmp)
      else
         write(*,*) 'FATAL ERROR: cannont handle ',nscalar,' scalars'
         stop
      endif
            
      return
      end

      subroutine Dirichlet2Neumann_2(nscalar, itype, tmp)
!
! This is an interface routine, designed to call return a value for
! the flux to a point on the wafer due to electrochemical deposition
! to Ken Jansen's PHASTA given a boundary conditions and an index for
! a particular feature.
!
! There is an inherent assumption that we are going to be doing
! electroplating. This routine sets up the filenames and the 
! top-of-the-domain boundary conditions.
!
      implicit none

      integer maxdata,maxtypes
      parameter(maxdata=100,maxtypes=5)

      integer itype, nscalar
      real*8 tmp(nscalar)
! For each table up to maxtypes, we have 4 pieces of data--two independent,
! two dependent--for each point, up to maxdata+1.
      real*8 table(4,0:maxdata,0:maxdata,maxtypes)
      save table

      integer i,j,n
      logical readfile(maxtypes)
      save readfile
      data (readfile(i),i=1,maxtypes) / maxtypes*.false./

      real*8  dx(2,maxtypes)
      integer numdata(2,maxtypes)
      save dx
      save numdata
      
      real*8  x,y, z(3,2)
! We can only deal with two parameter models for now.
      if(nscalar .ne. 2) then
         write(*,*) 'Sorry, Dirichlet2Neumann handles 2 scalars!'
         write(*,*) 'You asked for ', nscalar
         write(*,*) 'STOPPING...'
         stop
      endif

! If we haven't read in our parameters for this featuretype yet...

      if( .not. readfile(itype)) then
         readfile(itype) = .true.
         call readtable_2(itype,table,numdata,dx, & 
              maxdata,maxtypes)
      endif

      x = tmp(1)
      y = tmp(2)
      
      if(.false.) then
         if( x .gt. table(1,0,0,itype) .or.  & 
              x .lt. table(1,numdata(1,itype)-1,0,itype) ) then
            write(*,*) 'Sorry, concentration 1 asked for: ', x
            write(*,*) '  is out of the table bounds.'
            write(*,*)  '#1  [ ',table(1,0,0,itype), ' , ', & 
                 table(1,numdata(1,itype)-1,0,itype), ' ] ', & 
                 numdata(1,itype)-1
            
            write(*,*) '  STOPPING...'
            stop
         endif
         if( y .gt. table(2,0,0,itype) .or.  & 
              y .lt. table(2,0,numdata(2,itype)-1,itype) ) then
            write(*,*) 'Sorry, concentration 2 asked for: ', y
            write(*,*) '  is out of the table bounds.'
            write(*,*)  '#2   [ ',table(2,0,0,itype), ' , ', & 
                 table(2,0,numdata(2,itype)-1,itype), ' ] ', & 
                 numdata(2,itype)-1
            write(*,*) '  STOPPING...'
            stop
         endif
      endif

      i = int ( (x - table(1,0,0,itype) ) / dx(1,itype))
      j = int ( (y - table(2,0,0,itype) ) / dx(2,itype))
!      write(*,*) 'i,j,x,y: ',i,j,x,y
      if(i .lt. 0) then
         i = 0
!         x = table(1,0,0,itype)
!         write(*,*) 'Reseting i low: ',i,j,x,y
      endif
      if(j .lt. 0) then
         j = 0
         y = table(2,0,0,itype)
!         write(*,*) 'Reseting j low: ',i,j,x,y
      endif
      if(i .ge. numdata(1,itype)) then
         i = numdata(1,itype)-2
!         x = table(1,i+1,0,itype)
!         write(*,*) 'Reseting i high: ',i,j,x,y
      endif
      if(j .ge. numdata(2,itype)) then
         j = numdata(2,itype)-2
         y = table(1,0,j+1,itype)
!         write(*,*) 'Reseting j high: ',i,j,x,y
      endif

      do n=3,4
         
         z(1,1) = table(n,i,j,itype)
         z(3,1) = table(n,i+1,j,itype)
         z(1,2) = table(n,i,j+1,itype)
         z(3,2) = table(n,i+1,j+1,itype)

         z(2,1) = (z(3,1) - z(1,1)) / dx(1,itype) & 
              *(x-table(1,i,j,itype)) + z(1,1)
         z(2,2) = (z(3,2) - z(1,2)) / dx(1,itype) & 
              *(x-table(1,i,j,itype)) + z(1,2)
         tmp(n-2) = (z(2,2) - z(2,1))/dx(2,itype) & 
              *(y-table(2,i,j,itype)) + z(2,1)
      
      enddo
!      write(*,*) 'Interpolation from ',x,y,' to:', tmp(1),tmp(2)
      return

      end
!--------------------------------------------------------------------

      subroutine readtable_2(islot,table,numdata,dx,maxdata,maxslots)
!
!    Read a table of ordered quadruplets and place them into the slot in
!    TABLE that is assosciated with ISLOT. Store the number of 
!    data in NUMDATA and the spacing in DX.  The file to be read 
!    is 'TABLE.[islot]' The data but be in a rectangular regular grid.
!
!     AUTHOR: Max Bloomfield, May 2000
!
      implicit none
!
      integer islot
      integer maxslots,numdata(2,maxslots), maxdata
!
      real*8 table(4,0:maxdata,0:maxdata,maxslots), dx(2,maxslots)
      real*8 x1,x2,y1,y2,x2old
!
      character*250 linein,filename
!
      integer i,j,k
      logical verbose

      verbose = .true.

      i=0
      j=0
      write(filename,1066) islot
 1066 format('TABLE.',i1)

      open(file=filename,unit=1066)

      write(*,*) 'Opening ', filename

 1    continue 
      read (unit=1066,fmt='(a)',end=999) linein

      if (linein(1:1).eq.'#') then
         write (*,'(a)') linein
         goto 1
      endif
!
      if (i.gt.maxdata**2+maxdata+1) then
         write(*,*)  & 
              'reached the maximum number of data points allowed'
         write(*,*) 'FATAL ERROR: stopping'
         stop
      endif
!
      read (linein,*) x1,x2,y1,y2
      if(i .gt. 0 .and. x2 .ne. x2old) then
!        increment the outer index in this nested loop. That is, go on
!        to the next "row" (not in fortran speak, but in table speak.)
         j = j + 1
         i=0
      endif
         
      table(1,i,j,islot) = x1*1.d-0
      table(2,i,j,islot) = x2*1.d-0
      table(3,i,j,islot) = y1*1.d-0
      table(4,i,j,islot) = y2*1.d-0
!      
      i=i+1
      x2old = x2

      goto 1
!
 999  continue
!
      numdata(1,islot) = I
      numdata(2,islot) = j+1
!      
      dx(1,islot) = table(1,2,1,islot) - table(1,1,1,islot)
      dx(2,islot) = table(2,1,2,islot) - table(2,1,1,islot)

      if(verbose) then
         write(*,*) 'Table is ',i,' by ',j+1
         
         write(*,*) 'there are ',i*(j+1),' flux data points'
         write(*,*) 'closing unit 1066'
         close(1066)
!     
         write(*,*) 'The abscissa are ', & 
              dx(1,islot),' and ',dx(2,islot),' apart.'
         
         write(*,*) 'the flux data are '
         do i=0,numdata(1,islot)-1
            do j=0,numdata(2,islot)-1
               write(*,*) i,j,(table(k,i,j,islot), k=1,4)
            end do
         end do
         
      endif
      return
      end

      

      subroutine Dirichlet2Neumann_1(nscalar, itype, tmp)
!
! This is an interface routine, designed to call return a value for
! the flux to a point on the wafer due to electrochemical deposition
! to Ken Jansen's PHASTA given a boundary conditions and an index for
! a particular feature.
!
! There is an inherent assumption that we are going to be doing
! electroplating. This routine sets up the filenames and the 
! top-of-the-domain boundary conditions.
!
      implicit none

      integer maxdata,maxtypes
      parameter(maxdata=200,maxtypes=2)

      integer itype, nscalar
      real*8 tmp(nscalar)
      real*8 table(0:maxdata,2,maxtypes)
      save table

      integer i
      logical readfile(maxtypes)
      save readfile
      data (readfile(i),i=1,maxtypes) / maxtypes*.false./

      real*8  dx(maxtypes)
      save dx
      integer numdata(maxtypes)
      save numdata
      
      real*8 dt, conc_BC, flux_BC
! We can only deal with one parameter models for now.

      if(nscalar .ne. 1) then
         write(*,*) 'Sorry, Dirichlet2Neumann can only handle 1 scalar!'
         write(*,*) 'You asked for ', nscalar
         write(*,*) 'STOPPING...'
         stop
      endif

! If we haven't read in our parameters for this featuretype yet...

      if( .not. readfile(itype)) then
         readfile(itype) = .true.
         call readtable_1(itype,table,numdata(itype),dx(itype), & 
              maxdata,maxtypes)
!         write(*,*) 'back from readtable'
         if(dx(itype) .eq. 0.0d0) then
            write(*,*) 'DX for table ',itype,' is zero (I think!)'
            stop
         endif
      endif
!      write(*,*) 'returning from D2N'

      conc_BC = tmp(1)
      
!      if( conc_BC .lt. table(0,1,itype) .or. 
!     &    conc_BC .gt. table(numdata(itype),1,itype) ) then
!         write(*,*) 'Sorry, concentration asked for: ', conc_BC
!         write(*,*) '  is out of the table bounds.'
!         write(*,*) '[',table(0,1,itype),',
!     &        ',table(numdata(itype),1,itype),']'
!         write(*,*) '  STOPPING...'
!         stop
!      endif

      i = int ( (conc_BC - table(0,1,itype) ) / dx(itype))

      if( conc_BC .lt. table(0,1,itype))then
         i = 0
         conc_BC =  table(i,1,itype)
      elseif( conc_BC .gt. table(numdata(itype),1,itype)) then
         i = numdata(itype)
         conc_BC =  table(i,1,itype)
      endif
         

      dt = conc_BC - table(i,1,itype)
      flux_BC = dt * (table(i+1,2,itype) - table(i,2,itype)) + & 
           table(i,2,itype)


      tmp(1) = flux_BC
      

      end
!--------------------------------------------------------------------

      subroutine readtable_1(islot,table,numdata,dx,maxdata,maxslots)
!
!     Read a table of ordered pairs and place them into the slot in
!     TABLE that is assosciated with ISLOT. Store the number of 
!     data in NUMDATA and the spacing in DX.  The file to be read 
!     is 'TABLE.[islot]'
!
!     AUTHOR: Max Bloomfield, May 2000
!
      implicit none
!
      integer islot
      integer numdata, maxdata, maxslots
!
      real*8 table(0:maxdata,2,maxslots),dx
!
      character*80 linein,filename
!
      integer i,j
      logical verbose
      verbose = .true.

      i=-1
      
      write(filename,1066) islot
 1066 format('TABLE.',i1)
      open(file=filename,unit=1066)
      if(verbose) write(*,*) 'Opening ', filename

 1    continue 
      read (unit=1066,fmt='(a)',end=999) linein

      if (linein(1:1).eq.'#') then
         write (*,'(a)') linein
         goto 1
      endif
!
      i=i+1
      if (i.ge.maxdata) then
         write(*,*)  & 
              'reached the maximum number of data points allowed'
         write(*,*) 'FATAL ERROR: stopping'
         stop
      endif
!
      read (linein,*) table(i,1,islot), table(i,2,islot)
      table(i,1,islot)= table(i,1,islot)*1.0d-0
      table(i,2,islot)= table(i,2,islot)*1.0d-0
!      
      goto 1
!
 999  continue
!
      numdata = i
      dx = table(1,1,islot)-table(0,1,islot)
!      
      if(verbose) then
         write(*,*) 'there are ',numdata,' flux data points'
         write(*,*) 'closing unit 1066'
         close(1066)
!     
         write(*,*) 'the flux data are '
         do 101 j=0,i
            write(*,*) j,table(j,1,islot), table(j,2,islot)
 101     continue
      endif
      return
      end
