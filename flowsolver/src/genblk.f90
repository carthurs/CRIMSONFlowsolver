        subroutine genblk (IBKSZ)
!
!----------------------------------------------------------------------
!
!  This routine reads the interior elements and generates the
!  appropriate blocks.
!
! Zdenek Johan, Fall 1991.
!----------------------------------------------------------------------
!
        use pointer_data
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        integer, allocatable :: ientp(:,:)
        integer mater(ibksz)
        integer intfromfile(50) ! integers read from headers
        character*255 fname1
!
        iel=1
        itpblk=nelblk
        nelblk=0
        mattyp = 0
        ndofl = ndof
        nsymdl = nsymdf
        
        do iblk = 1, itpblk
!
!           read(igeom) neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst
!           call creadlist(igeom,7,
!     &          neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst)
           fname1='connectivity interior?'
           call readheader(igeom,fname1//c_null_char,intfromfile,7, &
                           c_char_"integer"//c_null_char, iotype)
           neltp  =intfromfile(1)
           nenl   =intfromfile(2)
           ipordl =intfromfile(3)
           nshl   =intfromfile(4)
           ijunk  =intfromfile(5)
           ijunk  =intfromfile(6)
           lcsyst =intfromfile(7)
           
           allocate (ientp(neltp,nshl))
!           read(igeom) ientp
           iientpsiz=neltp*nshl
           call readdatablock(igeom,fname1//c_null_char,ientp,iientpsiz, &
                           c_char_"integer"//c_null_char, iotype)

           ! Check if we're about to overrun the preallocated (MAXBLK-sized) data arrays:
           ! (Intel compilers using the -CB flags would notice this automatically; without such
           ! bounds-checking we would be in the realm of undefined behaviour, so guard it here
           ! and provide a more helpful error message)
           if (neltp/ibksz .gt. MAXBLK) then
            write(*,*) "Too much data per CPU (MAXBLK exceeded). Try running with more CPUs."
            write(*,*) "ibksz:", ibksz
            write(*,*) "neltp:", neltp
            write(*,*) "MAXBLK:", MAXBLK
            call exit(1)
           endif
           
           do n=1,neltp,ibksz 
              nelblk=nelblk+1
              npro= min(IBKSZ, neltp - n + 1)
!
              lcblk(1,nelblk)  = iel
!              lcblk(2,nelblk)  = iopen ! available for later use
              lcblk(3,nelblk)  = lcsyst
              lcblk(4,nelblk)  = ipordl
              lcblk(5,nelblk)  = nenl
              lcblk(6,nelblk)  = nfacel
              lcblk(7,nelblk)  = mattyp
              lcblk(8,nelblk)  = ndofl
              lcblk(9,nelblk)  = nsymdl 
              lcblk(10,nelblk) = nshl ! # of shape functions per elt
!
!.... allocate memory for stack arrays
!
              allocate (mmat(nelblk)%p(npro))
!
              allocate (mien(nelblk)%p(npro,nshl))

              call phglobalblockedarrayassignpointer (npro, nshl, c_loc(mien(nelblk)%p))

!              call phSolverUpdateBlockField( &
!              c_char_"connectivity interior linear tetrahedron"//c_null_char, &
!              npro, nshl, c_loc(mien(nelblk)%p) )

              allocate (mxmudmi(nelblk)%p(npro,maxsh))
!
!.... save the element block
!
              n1=n
              n2=n+npro-1
              mater=1   ! all one material for now
              call gensav (ientp(n1:n2,1:nshl), &
                           mater,           mien(nelblk)%p, &
                           mmat(nelblk)%p)
              iel=iel+npro
!
           enddo
           if (allocated(ientp)) then
             deallocate(ientp)
           endif
        enddo
        lcblk(1,nelblk+1) = iel
!
!.... return
!
!AD        call timer ('Back    ')
!
        return
!
1000    format(a80,//, &
        ' N o d a l   C o n n e c t i v i t y',//, &
        '   Elem  ',/, &
        '  Number  ',7x,27('Node',i2,:,2x))
1100    format(2x,i5,6x,27i8)
        end
