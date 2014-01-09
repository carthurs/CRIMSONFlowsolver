        subroutine genini (iBC, BC, y, ac, iper, ilwork, &
                     ifath, velbar,  nsons,x, &
                     shp,     shgl,    shpb,    shglb ) 
!
!----------------------------------------------------------------------
! This routine reads the initial values in primitive form (density,
! velocity and temperature), satisfies the boundary conditions and 
! converts them to Y-variables.
!
! input:
!  iBC    (nshg)               : boundary condition code
!  BC     (nshg,ndofBC)        : boundary condition constrain data
!   x     (numnp,nsd)	       : locations of nodes, numnp-> # of node
!				  nsd-> space dimension, 1=x, 2=y, 3=z 
!			
!
! output:
!  y      (nshg,ndof)          : initial values of Y variables
!
!
! Farzin Shakib, Winter 1986.
! Zdenek Johan,  Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        use specialBC   ! gets itvn from here
        use convolImpFlow ! brings in ntimeptpT and other variables
        use convolRCRFlow ! brings RCR variables
        use convolTRCRFlow ! brings time-varying RCR variables
        use convolCORFlow ! brings Coranary variables
        use incpBC        ! brings in INCP variables
        use LagrangeMultipliers ! brings in variables for Lagrange Multipliers
        use calcFlowPressure

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
        !include "auxmpi.h"
!
        dimension iBC(nshg),                iper(nshg), &
                  BC(nshg,ndofBC),          y(nshg,ndof), &
                  ac(nshg,ndof),            x(numnp,nsd), &
                  shp(MAXTOP,maxsh,MAXQPT), &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT), &
                  shpb(MAXTOP,maxsh,MAXQPT), &
                  shglb(MAXTOP,nsd,maxsh,MAXQPT)

!
        dimension ilwork(nlwork)
!
!  stuff for dynamic model s.w.avg and wall model
!
        dimension ifath(numnp),    velbar(nfath,nflow), &
                 nsons(nfath) 
!
!.... -------------------------->  Restart  <---------------------------
!
!.... read q from [RESTAR.INP], reset LSTEP
!
        !call restar ('in  ',  y,  ac)
!
!        if((itwmod.gt.0)  &
!           .or. (nsonmax.eq.1 .and. iLES.gt.0) ) then
!           call rwvelb('in  ',velbar,ifail)
!
! if the read failed calculate velbar
!
!           if(ifail.eq.1) then
!              call getvel (y,     ilwork, iBC, &
!                           nsons, ifath, velbar)
!           endif
!
!        endif   ! for the itwmod or irscale
!
!.... time varying boundary conditions as set from file bct.dat and impt.dat 
!     (see function for format in file bctint.f)
!
        if (itvn .gt. 0 ) then !for inlet velocities
           call initBCt( x, iBC, BC)
           call BCint(lstep*Delt(1),shp,shgl,shpb,shglb,x, BC, iBC)
        endif
        if (impfile .gt. 0 ) then !for impedance BC
           do irank=1, numpe
              call MPI_BARRIER (INEWCOMM,ierr)
              if((irank-1).eq.myrank) then 
                 write(*,*) 'reading Qhistor.dat'
                 open(unit=816, file='Qhistor.dat',status='old')
                 read (816,*) ntimeptpT
                 allocate (QHistImp(ntimeptpT+1,numImpSrfs)) 
                 allocate (QHistTry(ntimeptpT,numImpSrfs)) 
                 allocate (QHistTryF(ntimeptpT,numImpSrfs)) 
                 do j=1,ntimeptpT+1
                    read(816,*) (QHistImp(j,n),n=1,numImpSrfs) !read flow history
                 enddo
                 close(816)              
                 call initImpt() !read impedance data and initialize begin/end values
                 do i=2,ntimeptpT
                    call Impint((ntimeptpT-i+1)*Delt(1),i) !return Imp values in reverse order ZN->Z0
                 enddo

                 allocate (poldImp(0:MAXSURF)) !for pressure part that depends on the history only
              else 
                 continue
              endif
           enddo
        endif
        if (ircrfile .gt. 0 ) then !for RCR BC
           call initRCRt()      !read RCR data
           dtRCR(:) = Delt(1)/(ValueListRCR(2,:)*ValueListRCR(3,:))
           !last one needed only to have array of same size as imp BC
           allocate (poldRCR(0:MAXSURF)) !for pressure part that depends on the history only
           allocate (HopRCR(0:MAXSURF)) !for H operator contribution
        endif
        if (itrcrfile .gt. 0 ) then !for time-varying RCR BC
           call initTRCRt()      !read time-varying RCR data
        endif
        if (icorfile .gt. 0 ) then !for Coronary BC
           call initCORt()
           allocate (CORArea(numCORSrfs))
           allocate (CORic(2,numCORSrfs))
           allocate (poldCOR(0:MAXSURF)) !for pressure part that depends on the history only
           allocate (plvoldCOR(0:MAXSURF))
           allocate (HopCOR(0:MAXSURF)) !for H operator contribution
           CORArea = zero
           poldCOR = zero
           plvoldCOR = zero
           HopCOR = zerol
        endif
        if (incpfile .gt. 0) then !for INCP BC
           call initINCPt()  !read incp.dat, Elastance.dat
           allocate (INCPConvCoef(lstep+nstep(1)+2, numINCPSrfs))
           allocate (INCPCoef(2, numINCPSrfs))
           allocate (InflowArea(numINCPSrfs))
           allocate (poldINCP(0:MAXSURF))
           allocate (QHistINCP(lstep+nstep(1)+1, numINCPSrfs))
           allocate (Eadjust(numINCPSrfs))
           allocate (inactive(0:MAXSURF))
           inactive = 0
           QHistINCP = zero
           poldINCP = zero
        endif
        if (numCalcSrfs .gt. 0) then !for CalcSurfaces
           allocate (CalcArea(numCalcSrfs))
           CalcArea = zero
           call initCalcSrfst()
        endif
        if (iLagfile .gt. 0) then !for Lagrange multipliers
           call initLagrange()
           allocate(QLagrange(numLagrangeSrfs,2))
           allocate(PQLagrange(numLagrangeSrfs,2))
           allocate(IPLagrange(numLagrangeSrfs,6))
           allocate(NANBLagrange(6,nshg,3))
           QLagrange = zero
           PQLagrange = zero
           IPLagrange = zero
           NANBLagrange = zero
        endif
!
!
!.... satisfy the boundary conditions
!

!        call itrBC (y, ac,  iBC, BC, iper, ilwork)

        itempr=mod(impl(1),2)  ! tempr solve if impl odd
        if(itempr.eq.1) then
           isclr=0
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        endif
        do isclr=1,nsclr
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        enddo

!        if((irscale.ge.0) .and. (myrank.eq.master)) then
!           call genscale(y, x, iBC)
!        endif
!
!.... --------------------------->  Echo  <----------------------------
!
!.... echo the initial data
!
        if ((necho .lt. 0).and.(myrank.eq.master)) then
          do n = 1, nshg
            if (mod(n,50) .eq. 1) write(iecho,1000) ititle,(i,i=1,ndof)
            write (iecho,1100) n, (y(n,i),i=1,ndof)
          enddo
        endif
!
!.... return
!
        return
!
1000    format(a80,//, &
        ' I n i t i a l   V a l u e s                        ',//, &
        '    Node ',/, &
        '   Number ',6x,6('dof',i1,:,10x))
1100    format(1p,2x,i5,5x,5(e12.5,2x))
!
        end

