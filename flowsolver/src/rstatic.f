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
        subroutine rstatic (res, y, Dy)
c
c----------------------------------------------------------------------
c
c This subroutine calculates the statistics of the residual.
c
c input:
c  res   (nshg,nflow)   : preconditioned residual
c
c output:
c  The time step, cpu-time and entropy-norm of the residual
c     are printed in the file HISTOR.DAT.
c
c----------------------------------------------------------------------
c
        use ResidualControl 
        
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension res(nshg,nflow),    mproc(1),  rvec(numpe)
        dimension rtmp(nshg),        nrsmax(1)
        dimension irecvcount(numpe), resvec(numpe)

        real*8    y(nshg,ndof),    Dy(nshg,4)
        integer tmrc
c
c$$$	ttim(68) = ttim(68) - tmr()
c
c.... compute max delta y
c
        rdy1 = zero
        rdy2 = zero
        rdy4 = zero
        rdy5 = zero
        call sumgatN( abs(gami*Delt(itseq) 
     &                * Dy(1:numnp,1:3)),3,rdy1, numnp)
        call sumgatN( abs( y(1:numnp,1:3)),3,rdy2,numnp)
        call sumgatN( abs(gami*alfi*Delt(itseq)
     &                * Dy(1:numnp,4)),1,rdy4,numnp)
        call sumgatN( abs( y(1:numnp,4)),  1,rdy5,numnp)
        rmaxdyU = rdy1/rdy2
        rmaxdyP = rdy4/rdy5
	
c
c..... Signal to quit if delta is very small. look in itrdrv.f for the
c      completion of the hack.
c
	if( rmaxdyU .lt. dtol(1) .and. rmaxdyP .lt. dtol(2)) then
           istop = 1000
        endif

        if (numpe == 1) nshgt=nshg   ! global = this processor
c
c
c.... ----------------------->  Convergence  <-------------------------
c
c.... compute the maximum residual and the corresponding node number
c
        rtmp = zero
        if (impl(itseq) .ge. 9) then
          do i = 1, nflow
            rtmp = rtmp + res(:,i)**2    ! only add continuity and momentum
          enddo
        endif

        call sumgat (rtmp, 1, resnrm)
        
        resmaxl = maxval(rtmp)

        irecvcount = 1
        resvec = resmaxl
        if (numpe > 1) then
           call MPI_REDUCE_SCATTER (resvec, resmax, irecvcount,
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
           if (resmax .eq. resvec(1) ) then
              mproc(1) = myrank
              nrsmax   = maxloc(rtmp)
           else
              nrsmax(1) = -1
              mproc(1)  = -1
           endif
           resvec = nrsmax
           call MPI_REDUCE_SCATTER (resvec, rvec, irecvcount,
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
           nrsmax = rvec(1)
           resvec = mproc
           call MPI_REDUCE_SCATTER (resvec, rvec, irecvcount,
     &          MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
           mproc = rvec(1)
      else
          resmax   = resmaxl
          nrsmax   = maxloc(rtmp)
          mproc(1) = 0
      endif
c
c.... correct the residuals
c
        if (loctim(itseq) .eq. 0) then
          resnrm = resnrm 
          resmax = resmax
        else
          resnrm = resnrm
          resmax = resmax
        endif
c
c.... approximate the number of entries
c
        totres = resnrm / float(nshgt)
        totres = sqrt(totres)
        resmax = sqrt(resmax)
        if (resfrt .eq. zero) resfrt = totres
        jtotrs = int  ( 10.d0 * log10 ( totres / resfrt ) )
        jresmx = int  ( 10.d0 * log10 ( resmax / totres ) )
        
        if(rescontrol .gt. 0) then
           controlResidual = totres
           CurrentIter = CurrentIter + 1
        endif
c     
c.... get the CPU-time
c
CAD        cputme = (second(0) - ttim(100))
        intsec=TMRC()
        cputme = (intsec - ttim(100))
c
c.... output the result
c
        if (numpe > 1) call MPI_BARRIER (MPI_COMM_WORLD, ierr)
        
        if (myrank .eq. master) then
c
c.... results of continuity and momentum 
c
           
           print 2000, lstep+1, cputme, totres, jtotrs, rmaxdyU,
     &          rmaxdyP,nrsmax,
     &          mproc(1)+1, jresmx, int(statsflow(4)),
     &          int(statsflow(1))
           write (ihist,2000) lstep+1, cputme, totres, jtotrs, 
     &          rmaxdyU, rmaxdyP, nrsmax,
     &          mproc(1)+1,jresmx,int(statsflow(4)),
     &          int(statsflow(1))
           
           call flush(ihist)
        endif
        if(numpe>1) call MPI_BARRIER (MPI_COMM_WORLD,ierr)

c$$$	ttim(68) = ttim(68) + tmr()

c
c.... return
c
        return
c
 1000   format(1p,i6,5e13.5)
 2000   format(1p,i6,e10.3,e10.3,2x,'(',i4,')',2x,e10.3,2x,e10.3,
     &       2x,'<',i6,'-',i2,'|',
     &       i4,'>', ' [', i4,' -',i4,']')
 3000   format(1p,i6,e10.3,e10.3,3x,'(',i4,')',3x,'<',i6,'-',i2,'|',
     &       i4,'>', ' [', i4,' -',i4,' -',i4,']')

c
        end


        subroutine rstaticSclr (res, y, Dy, icomp)
c
c----------------------------------------------------------------------
c
c This subroutine calculates the statistics of the residual
c
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        dimension res(nshg)
        dimension rtmp(nshg)
        real*8    y(nshg,ndof),    Dy(nshg), nrm
        integer tmrc
c
c.... compute max delta y
c
        rdy1 = zero
        rdy2 = zero
c
c.... normalize turbulence with molecular viscosity
c        
        if ( (icomp .eq. 6).and. (iRANS.eq.-1) ) then
           nrm = datmat(1,2,1)
        else 
           nrm = zero
        endif
        call sumgat( abs(gami*Delt(itseq)*Dy(:)),1,rdy1)
        call sumgat( abs( y(:,icomp)),1,rdy2)
        rmaxdyT = rdy1/(rdy2+nrm)
c
c.... compute the maximum residual and the corresponding node number
c
        rtmp = zero
        rtmp = rtmp + res**2 ! add temperature also
        call sumgat (rtmp, 1, resnrm)

        if (numpe == 1) nshgt=nshg ! global = this processor

        totres = resnrm / float(nshgt)
        totres = sqrt(totres)

c        if (mod(impl(1),100)/10 .eq. 0) then  !not solving flow
           if (myrank .eq. master) then
c     
c.... get the CPU-time
c
              intsec=TMRC()
              cputme = (intsec - ttim(100))

           print 802, lstep+1, cputme, totres, rmaxdyT,
     &                int(statssclr(1))
           write (ihist,802) lstep+1, cputme, totres, 
     &          rmaxdyT,int(statssclr(1))
           
               call flush(ihist)
           endif
c        else 
c           if (myrank .eq. master) then
c              print 803, totres, rmaxdyT, int(statssclr(1))
c              write(ihist,803) totres, rmaxdyT, int(statssclr(1))
c           endif
c        endif

        return
        
 802    format(1p,i6,e10.3,e10.3,10X,e10.3,31X'[',i6,']')
 803    format(1p,16x,e10.3,10x,e10.3,31X,'[',i10,']')
    
        end



