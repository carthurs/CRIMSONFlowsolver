        subroutine proces() bind(C, name="proces")
!
!----------------------------------------------------------------------
!
! This subroutine generates the problem data and calls the solution
!  driver.
!
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        !use readarrays          ! used to access x, iper, ilwork
        use iso_c_binding
        use globalArrays
        use shapeTable
        use turbsa          ! used to access d2wall

        use boundarymodule, only : setupBoundaryModule  ! used to setup boundary module - added by kdl 09/04/14
!        use pvsQbi, only : ndsurf                       ! used to setup boundary module - added by kdl 09/04/14        

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
!
! arrays in the following 2 lines are now dimensioned in readnblk
!        dimension x(numnp,nsd)
!        dimension iper(nshg), ilwork(nlwork)
!
!        dimension y(nshg,ndof), &
!                  iBC(nshg), &
!                  BC(nshg,ndofBC), &
!                  ac(nshg,ndof)
!
!.... shape function declarations
!     
!        dimension shp(MAXTOP,maxsh,MAXQPT),   &
!                  shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
!                  shpb(MAXTOP,maxsh,MAXQPT), &
!                  shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!
!  stuff for dynamic model s.w.avg and wall model
!
!        dimension ifath(numnp),    velbar(nfath,nflow),
!     &            nsons(nfath) 

!        dimension velbar(nfath,nflow)
!
! stuff to interpolate profiles at inlet
!
        real*8 bcinterp(100,ndof+1),interp_mask(ndof)
        logical exlog
        
        call initShapeTable

!        if ((irscale .ge. 0) .and. (myrank.eq.master)) then
!           call setSPEBC(numnp, point2nfath, nsonmax)
!        endif


!     
!.... generate the geometry and boundary conditions data
!
        allocate(ifath(numpe)) ! \todo-binary WARNING - remove if you ever want to actually use ifath -this is an empty dummy to avoid memory issues for now
        allocate(nsons(nfath)) ! \todo-binary WARNING - remove if you ever want to actually use nsons -this is an empty dummy to avoid memory issues for now
        call gendat (y,              ac,             x, &
                     iBC,            BC, &
                     iper,     ilwork,   shp, &
                     shgl,           shpb,           shglb, &
                     ifath,    velbar,         nsons )      
        deallocate(ifath) ! \todo-binary WARNING - remove if you ever want to actually use ifath -this is an empty dummy to avoid memory issues for now
        deallocate(nsons) ! \todo-binary WARNING - remove if you ever want to actually use nsons -this is an empty dummy to avoid memory issues for now

        call setper(nshg)
        
        call perprep(iBC,iper,nshg)
        
!        if (iLES/10 .eq. 1) then
!           call keeplhsG ! Allocating space for the mass (Gram) matrix to be
!                         ! used for projection filtering and reconstruction
!                         ! of the strain rate tensor.
!
!           call setrls   ! Allocating space for the resolved Leonard stresses
!                         ! See bardmc.f
!        endif
!
!.... time averaged statistics
!
        if (ioform .eq. 2) then
           call initStats(x, iBC, iper, ilwork) 
        endif
!
!.... RANS turbulence model
!
!        if (iRANS .lt. 0) then
!           call initTurb( x )
!        endif
!
!.... p vs. Q boundary
!
        call initNABI( x, shpb )
!
! *** setup boundary data module 
!
        call setupBoundaryModule(x)
!     
!.... check for execution mode
!
        if (iexec .eq. 0) then
           lstep = 0
           call restar ('out ',  y  ,ac)
           return
        endif
!
!.... initialize AutoSponge
!
        if(matflg(5,1).ge.4) then ! cool case (sponge)
           call initSponge( y,x) 
        endif
!
!
!.... adjust BC's to interpolate from file
!
        inquire(file="inlet.dat",exist=exlog)
        if(exlog) then
           open (unit=654,file="inlet.dat",status="old")
           read(654,*) ninterp,ibcmatch,(interp_mask(j),j=1,ndof)
           do i=1,ninterp
              read(654,*) (bcinterp(i,j),j=1,ndof+1) ! distance to wall+
                        ! ndof but note order of BC's
                        ! p,t,u,v,w,scalars. Also note that must be in
                        ! increasing distance from the wall order.

           enddo
           do i=1,nshg  ! only correct for linears at this time
              if(mod(iBC(i),1024).eq.ibcmatch) then
                 iupper=0
                 do j=2,ninterp
                    if(bcinterp(j,1).gt.d2wall(i)) then !bound found
                       xi=(d2wall(i)-bcinterp(j-1,1))/ &
                          (bcinterp(j,1)-bcinterp(j-1,1))
                       iupper=j
                       exit
                    endif
                 enddo
                 if(iupper.eq.0) then
                    write(*,*) 'failure in finterp, ynint=',d2wall(i)
                    stop
                 endif
                 if(interp_mask(1).ne.zero) then 
                    BC(i,1)=(xi*bcinterp(iupper,2) &
                      +(one-xi)*bcinterp(iupper-1,2))*interp_mask(1)
                 endif
                 if(interp_mask(2).ne.zero) then 
                    BC(i,2)=(xi*bcinterp(iupper,3) &
                      +(one-xi)*bcinterp(iupper-1,3))*interp_mask(2)
                 endif
                 if(interp_mask(3).ne.zero) then 
                    BC(i,3)=(xi*bcinterp(iupper,4) &
                      +(one-xi)*bcinterp(iupper-1,4))*interp_mask(3)
                 endif
                 if(interp_masK(4).ne.zero) then 
                    BC(i,4)=(xi*bcinterp(iupper,5) &
                      +(one-xi)*bcinterp(iupper-1,5))*interp_mask(4)
                 endif
                 if(interp_mask(5).ne.zero) then 
                    BC(i,5)=(xi*bcinterp(iupper,6) &
                      +(one-xi)*bcinterp(iupper-1,6))*interp_mask(5)
                 endif
                 if(interp_mask(6).ne.zero) then 
                    BC(i,7)=(xi*bcinterp(iupper,7) &
                      +(one-xi)*bcinterp(iupper-1,7))*interp_mask(6)
                 endif
                 if(interp_mask(7).ne.zero) then 
                    BC(i,8)=(xi*bcinterp(iupper,8) &
                      +(one-xi)*bcinterp(iupper-1,8))*interp_mask(7)
                 endif
              endif
           enddo
        endif
!$$$$$$$$$$$$$$$$$$$$

!
!
!.... call the semi-discrete predictor multi-corrector iterative driver
!
!        call itrdrv (y,              ac,              &
!                     uold,           x, &
!                     iBC,            BC, &
!                     iper,     ilwork,   shp, &
!                     shgl,           shpb,           shglb, &
!                     ifath,    velbar,         nsons ) 
!        call itrdrv
!
!.... return
!
!
!.... stop CPU-timer
!
!AD        call timer ('End     ')
!
!.... close echo file
!
        close (iecho)
!
!.... end of the program
!
!AD        write(6,*) 'Life: ', second(0) - ttim(100)
!        deallocate(iper)
!        if(numpe.gt.1) deallocate(ilwork)
!        deallocate(x)
!
!        if((irscale.ge.0).or. ((iLES .lt. 20) .and. (iLES.gt.0)) &
!                         .or. (itwmod.gt.0)  ) then ! don't forget same
!                                                    ! conditional in
!                                                    ! readnblk2.f
!           deallocate(nsons)
!           deallocate(ifath)
!        endif
        return
        end


