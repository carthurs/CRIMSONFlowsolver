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
        subroutine proces
c
c----------------------------------------------------------------------
c
c This subroutine generates the problem data and calls the solution
c  driver.
c
c----------------------------------------------------------------------
c
        use readarrays          ! used to access x, iper, ilwork
        use turbsa          ! used to access d2wall
        use measureWallDistance
        include "common.h"
        include "mpif.h"
c
c arrays in the following 2 lines are now dimensioned in readnblk
c        dimension x(numnp,nsd)
c        dimension iper(nshg), ilwork(nlwork)
c
        dimension y(nshg,ndof),
     &            iBC(nshg),
     &            BC(nshg,ndofBC),
     &            ac(nshg,ndof)
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
c        dimension ifath(numnp),    velbar(nfath,nflow),
c     &            nsons(nfath) 

        dimension velbar(nfath,nflow)
c
c stuff to interpolate profiles at inlet
c
        real*8 bcinterp(100,ndof+1),interp_mask(ndof)
        logical exlog

c        if ((irscale .ge. 0) .and. (myrank.eq.master)) then
c           call setSPEBC(numnp, point2nfath, nsonmax)
c        endif

c
c.... initialize distance evaluation
c      
        if(ideformwall.eq.1 .and. imeasdist.eq.1) then
c         read "observed mesh data"
          write(*,*) "reading observed mesh at rank ",myrank     
          call dm_dataread("observed.dat")

c         build octrees for observed meshes
          write(*,*) "building octrees at rank ",myrank
          call dm_initialize()
        end if

c     
c.... generate the geometry and boundary conditions data
c
        call gendat (y,              ac,             point2x,
     &               iBC,            BC,
     &               point2iper,     point2ilwork,   shp,
     &               shgl,           shpb,           shglb,
     &               point2ifath,    velbar,         point2nsons )
        
        call setper(nshg)
        
        call perprep(iBC,point2iper,nshg)
        
        if (iLES/10 .eq. 1) then
           call keeplhsG ! Allocating space for the mass (Gram) matrix to be
                         ! used for projection filtering and reconstruction
                         ! of the strain rate tensor.

           call setrls   ! Allocating space for the resolved Leonard stresses
                         ! See bardmc.f 
        endif
c
c.... time averaged statistics
c
        if (ioform .eq. 2) then
           call initStats(point2x, iBC, point2iper, point2ilwork) 
        endif
c
c.... RANS turbulence model
c
        if (iRANS .lt. 0) then
           call initTurb( point2x )
        endif
c
c.... p vs. Q boundary
c
           call initNABI( point2x, shpb )
c     
c.... check for execution mode
c
        if (iexec .eq. 0) then
           lstep = 0
           call restar ('out ',  y  ,ac)
           return
        endif
c
c.... initialize AutoSponge
c
        if(matflg(5,1).ge.4) then ! cool case (sponge)
           call initSponge( y,point2x) 
        endif
c
c
c.... adjust BC's to interpolate from file
c
        
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
                       xi=(d2wall(i)-bcinterp(j-1,1))/
     &                    (bcinterp(j,1)-bcinterp(j-1,1))
                       iupper=j
                       exit
                    endif
                 enddo
                 if(iupper.eq.0) then
                    write(*,*) 'failure in finterp, ynint=',d2wall(i)
                    stop
                 endif
                 if(interp_mask(1).ne.zero) then 
                    BC(i,1)=(xi*bcinterp(iupper,2)
     &                +(one-xi)*bcinterp(iupper-1,2))*interp_mask(1)
                 endif
                 if(interp_mask(2).ne.zero) then 
                    BC(i,2)=(xi*bcinterp(iupper,3)
     &                +(one-xi)*bcinterp(iupper-1,3))*interp_mask(2)
                 endif
                 if(interp_mask(3).ne.zero) then 
                    BC(i,3)=(xi*bcinterp(iupper,4)
     &                +(one-xi)*bcinterp(iupper-1,4))*interp_mask(3)
                 endif
                 if(interp_masK(4).ne.zero) then 
                    BC(i,4)=(xi*bcinterp(iupper,5)
     &                +(one-xi)*bcinterp(iupper-1,5))*interp_mask(4)
                 endif
                 if(interp_mask(5).ne.zero) then 
                    BC(i,5)=(xi*bcinterp(iupper,6)
     &                +(one-xi)*bcinterp(iupper-1,6))*interp_mask(5)
                 endif
                 if(interp_mask(6).ne.zero) then 
                    BC(i,7)=(xi*bcinterp(iupper,7)
     &                +(one-xi)*bcinterp(iupper-1,7))*interp_mask(6)
                 endif
                 if(interp_mask(7).ne.zero) then 
                    BC(i,8)=(xi*bcinterp(iupper,8)
     &                +(one-xi)*bcinterp(iupper-1,8))*interp_mask(7)
                 endif
              endif
           enddo
        endif
c$$$$$$$$$$$$$$$$$$$$

c
c
c.... call the semi-discrete predictor multi-corrector iterative driver
c
        call itrdrv (y,              ac,             
     &               uold,           point2x,
     &               iBC,            BC,
     &               point2iper,     point2ilwork,   shp,
     &               shgl,           shpb,           shglb,
     &               point2ifath,    velbar,         point2nsons ) 
c
c.... return
c
c
c.... stop CPU-timer
c
CAD        call timer ('End     ')
c
c.... close echo file
c
        close (iecho)
c
c.... end of the program
c
CAD        write(6,*) 'Life: ', second(0) - ttim(100)
        deallocate(point2iper)
        if(numpe.gt.1) deallocate(point2ilwork)
        deallocate(point2x)

        if((irscale.ge.0).or. ((iLES .lt. 20) .and. (iLES.gt.0))
     &                   .or. (itwmod.gt.0)  ) then ! don't forget same
                                                    ! conditional in
                                                    ! readnblk2.f
           deallocate(point2nsons)
           deallocate(point2ifath)
        endif
        return
        end


