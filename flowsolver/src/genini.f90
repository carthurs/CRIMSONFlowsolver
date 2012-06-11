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
        subroutine genini (iBC, BC, y, ac, iper, ilwork,
     &               ifath, velbar,  nsons,x,
     &               shp,     shgl,    shpb,    shglb ) 
c
c----------------------------------------------------------------------
c This routine reads the initial values in primitive form (density,
c velocity and temperature), satisfies the boundary conditions and 
c converts them to Y-variables.
c
c input:
c  iBC    (nshg)               : boundary condition code
c  BC     (nshg,ndofBC)        : boundary condition constrain data
c   x     (numnp,nsd)	       : locations of nodes, numnp-> # of node
c				  nsd-> space dimension, 1=x, 2=y, 3=z 
C			
c
c output:
c  y      (nshg,ndof)          : initial values of Y variables
c
c----------------------------------------------------------------------
c
        use specialBC   ! gets itvn from here
        use convolImpFlow ! brings in ntimeptpT and other variables
        use convolRCRFlow ! brings RCR variables
        use convolCORFlow ! brings Coranary variables
        use incpBC        ! brings in INCP variables
        use LagrangeMultipliers ! brings in variables for Lagrange Multipliers
        use calcFlowPressure
        
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension iBC(nshg),                iper(nshg),
     &            BC(nshg,ndofBC),          y(nshg,ndof),
     &            ac(nshg,ndof),            x(numnp,nsd),
     &            shp(MAXTOP,maxsh,MAXQPT),
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT),
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT)

c
        dimension ilwork(nlwork)
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,nflow),
     &           nsons(nfath) 
c
c.... -------------------------->  Restart  <---------------------------
c
c.... read q from [RESTAR.INP], reset LSTEP
c
        call restar ('in  ',  y,  ac)
c
        if((itwmod.gt.0) 
     &     .or. (nsonmax.eq.1 .and. iLES.gt.0) ) then 
           call rwvelb('in  ',velbar,ifail)
c
c if the read failed calculate velbar
c
           if(ifail.eq.1) then
              call getvel (y,     ilwork, iBC,
     &                     nsons, ifath, velbar)
           endif
 
        endif   ! for the itwmod or irscale
c
c.... time varying boundary conditions as set from file bct.dat and impt.dat 
c     (see function for format in file bctint.f)
c
        if (itvn .gt. 0 ) then !for inlet velocities
           call initBCt( x, iBC, BC)
           call BCint(lstep*Delt(1),shp,shgl,shpb,shglb,x, BC, iBC)
        endif
        if (impfile .gt. 0 ) then !for impedance BC
           do irank=1, numpe
              call MPI_BARRIER (MPI_COMM_WORLD,ierr)
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
           HopCOR = zero
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
           allocate (FlowHist(lstep+nstep(1)+1,numCalcSrfs)) !for flow history
           allocate (PressHist(lstep+nstep(1)+1,numCalcSrfs)) !for pressure history
           FlowHist = zero
           PressHist = zero
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
c
c
c.... satisfy the boundary conditions
c

c        call itrBC (y, ac,  iBC, BC, iper, ilwork)

        itempr=mod(impl(1),2)  ! tempr solve if impl odd
        if(itempr.eq.1) then
           isclr=0
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        endif
        do isclr=1,nsclr
           call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
        enddo

        if((irscale.ge.0) .and. (myrank.eq.master)) then
           call genscale(y, x, iBC)
        endif
c
c.... --------------------------->  Echo  <----------------------------
c
c.... echo the initial data
c
        if ((necho .lt. 0).and.(myrank.eq.master)) then
          do n = 1, nshg
            if (mod(n,50) .eq. 1) write(iecho,1000) ititle,(i,i=1,ndof)
            write (iecho,1100) n, (y(n,i),i=1,ndof)
          enddo
        endif
c
c.... return
c
        return
c
1000    format(a80,//,
     &  ' I n i t i a l   V a l u e s                        ',//,
     &  '    Node ',/,
     &  '   Number ',6x,6('dof',i1,:,10x))
1100    format(1p,2x,i5,5x,5(e12.5,2x))
c
        end

