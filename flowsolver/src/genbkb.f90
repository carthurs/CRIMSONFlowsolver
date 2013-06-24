        subroutine genbkb (ibksz)
!
!----------------------------------------------------------------------
!
!  This routine reads the boundary elements, reorders them and
!  generates traces for the gather/scatter operations.
!
! Zdenek Johan, Fall 1991.
!----------------------------------------------------------------------
!
        use dtnmod
        use deformableWall
        use pointer_data
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!

        integer, allocatable :: ientp(:,:),iBCBtp(:,:)
        real*8, allocatable :: BCBtp(:,:), SWBtp(:,:)
        integer materb(ibksz)
        integer intfromfile(50) ! integers read from headers
        character*255 fname1

        iel=1
        itpblk=nelblb
        nelblb=0
        mattyp=0
        ndofl = ndof
        do iblk = 1, itpblk
           ieight=8
           fname1='connectivity boundary?'
           call readheader(igeom,fname1//c_null_char,intfromfile,ieight, &
                           c_char_"integer"//c_null_char,iotype)
           neltp =intfromfile(1)
           nenl  =intfromfile(2)
           ipordl=intfromfile(3)
           nshl  =intfromfile(4)
           nshlb =intfromfile(5)
           nenbl =intfromfile(6)
           lcsyst=intfromfile(7)
           numnbc=intfromfile(8)
           
           !call phSolverUpdateField( &
           !     c_char_"connectivity boundary linear tetrahedron"//c_null_char, &
           !     c_char_"Surface"//c_null_char, &
           !     c_char_"Integer"//c_null_char, &
           !     neltp, nshl, c_null_ptr, c_null_ptr)
!
           allocate (ientp(neltp,nshl))
           allocate (iBCBtp(neltp,ndiBCB))
           allocate (BCBtp(neltp,ndBCB))
           if (ideformwall.eq.1) then
             allocate (SWBtp(neltp,nProps))
           end if
           iientpsiz=neltp*nshl
           call readdatablock(igeom,fname1//c_null_char,ientp,iientpsiz, &
                           c_char_"integer"//c_null_char,iotype)
!     
!.... Read the boundary flux codes
!     
           fname1='nbc codes?'
           call readheader(igeom,fname1//c_null_char,intfromfile,ieight, &
                           c_char_"integer"//c_null_char,iotype)
           !call phSolverUpdateField( &
           !     c_char_"nbc codes linear tetrahedron"//c_null_char, &
           !     c_char_"Surface"//c_null_char, &
           !     c_char_"Integer"//c_null_char, &
           !     neltp, ndiBCB, c_null_ptr, c_null_ptr)
                           
           iiBCBtpsiz=neltp*ndiBCB
           call readdatablock(igeom,fname1//c_null_char,iBCBtp,iiBCBtpsiz, &
                           c_char_"integer"//c_null_char,iotype)
!     
!.... read the boundary condition data
!     
           fname1='nbc values?'
           call readheader(igeom,fname1//c_null_char,intfromfile,ieight, &
                           c_char_"integer"//c_null_char,iotype)
           !call phSolverUpdateField( &
           !     c_char_"nbc values linear tetrahedron"//c_null_char, &
           !     c_char_"Surface"//c_null_char, &
           !     c_char_"Double"//c_null_char, &
           !     neltp, ndBCB, c_null_ptr, c_null_ptr)
                           
           BCBtp    = zero
           iBCBtpsiz=neltp*ndBCB
           call readdatablock(igeom,fname1//c_null_char,BCBtp,iBCBtpsiz, &
                           c_char_"double"//c_null_char,iotype)
!
! This is a temporary fix until NSpre properly zeros this array where it
! is not set.  DEC has indigestion with these arrays though the
! result is never used (never effects solution).
!

           where(.not.btest(iBCBtp(:,1),0)) BCBtp(:,1)=zero
           where(.not.btest(iBCBtp(:,1),1)) BCBtp(:,2)=zero
           where(.not.btest(iBCBtp(:,1),3)) BCBtp(:,6)=zero
           if(ndBCB.gt.6) then
		      do i=6,ndof
		         where(.not.btest(iBCBtp(:,1),i-1)) BCBtp(:,i+1)=zero
	        enddo
	     endif
           where(.not.btest(iBCBtp(:,1),2)) 
              BCBtp(:,3)=zero
              BCBtp(:,4)=zero
              BCBtp(:,5)=zero
           endwhere

!
!.... Special arrays for deformable wall
!

!
!.... read the temporary SWB array 
!

           SWBtp = zero
           
           if(ideformwall.eq.1.and.iUseSWB.gt.0) then
              itwo=2                      
              fname1='SWB array?'
              call readheader(igeom,fname1//c_null_char,intfromfile,itwo,c_char_"double"//c_null_char, &
                              iotype)
              numelb=intfromfile(1)
              nProps=intfromfile(2)
              
              call readdatablock(igeom,fname1//c_null_char,SWBtp,numelb*nProps, &
                             c_char_"double"//c_null_char,iotype)
           endif

           do n=1,neltp,ibksz 
              nelblb=nelblb+1
              npro= min(IBKSZ, neltp - n + 1)
!
              lcblkb(1,nelblb)  = iel
!              lcblkb(2,nelblb)  = iopen ! available for later use
              lcblkb(3,nelblb)  = lcsyst
              lcblkb(4,nelblb)  = ipordl
              lcblkb(5,nelblb)  = nenl
              lcblkb(6,nelblb)  = nenbl
              lcblkb(7,nelblb)  = mattyp
              lcblkb(8,nelblb)  = ndofl
              lcblkb(9,nelblb)  = nshl 
              lcblkb(10,nelblb) = nshlb ! # of shape functions per elt
!
!.... save the element block
!
              n1=n
              n2=n+npro-1
              materb=1   ! all one material for now
!
!.... allocate memory for arrays
!

              allocate (mienb(nelblb)%p(npro,nshl))
!
              allocate (miBCB(nelblb)%p(npro,ndiBCB))
!
              allocate (mBCB(nelblb)%p(npro,nshlb,ndBCB))
!              
              if (ideformwall.eq.1) then
                allocate (mSWB(nelblb)%p(npro,nProps))
              end if
                           
!
              allocate (mmatb(nelblb)%p(npro))
!
!.... save the boundary element block
!
              if (ideformwall.eq.1) then
                call gensvbDef (ientp(n1:n2,1:nshl), &
                                iBCBtp(n1:n2,:),      BCBtp(n1:n2,:),  &
                                SWBtp(n1:n2,:),                        &
                                materb,               mienb(nelblb)%p, &
                                miBCB(nelblb)%p,      mBCB(nelblb)%p,  &
                                mSWB(nelblb)%p,                        &
                                mmatb(nelblb)%p)
              
              else  
                call gensvb (ientp(n1:n2,1:nshl), &
                             iBCBtp(n1:n2,:),      BCBtp(n1:n2,:), &
                             materb,        mienb(nelblb)%p, &
                             miBCB(nelblb)%p,        mBCB(nelblb)%p, &
                             mmatb(nelblb)%p)
!
              endif
              iel=iel+npro
           enddo
           
           
           deallocate(ientp)
           deallocate(iBCBtp)
           deallocate(BCBtp)
           if (ideformwall.eq.1) then
             deallocate(SWBtp)
           end if
           
           
        enddo
        lcblkb(1,nelblb+1) = iel

!
!.... return
!
        return
!
!.... end of file error handling
!
 911    call error ('genbcb  ','end file',igeom)
!
1000    format(a80,//, &
        ' B o u n d a r y   E l e m e n t   C o n n e c t i v i t y',//, &
        '   Elem   BC codes',/, &
        '  Number  C P V H ',5x,27('Node',i1,:,2x))
1100    format(2x,i5,2x,4i2,3x,27i7)
!$$$2000    format(a80,//,
!$$$     &  ' B o u n d a r y   E l e m e n t   B C   D a t a ',//,
!$$$     &  '   Node   ',3x,'mass',/,
!$$$     &  '  Number  ',3x,'flux',6x,'Pressure',6x,'Heat',6x,
!$$$     &  3('Viscous',i1,:,4x))
2100    format(2x,i5,1p,1x,6e12.4)
!
        end




