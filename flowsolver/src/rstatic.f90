        subroutine rstatic (res, y, Dy)
!
!----------------------------------------------------------------------
!
! This subroutine calculates the statistics of the residual.
!
! input:
!  res   (nshg,nflow)   : preconditioned residual
!
! output:
!  The time step, cpu-time and entropy-norm of the residual
!     are printed in the file HISTOR.DAT.
!  
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        use ResidualControl 
        
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
        !include "auxmpi.h"
!
        dimension res(nshg,nflow),    mproc(1),  rvec(numpe)
        dimension rtmp(nshg),        nrsmax(1)
        dimension irecvcount(numpe), resvec(numpe)

        real*8    y(nshg,ndof),    Dy(nshg,4)

!
!$$$	ttim(68) = ttim(68) - tmr()
!
!.... compute max delta y
!
        rdy1 = zero
        rdy2 = zero
        rdy4 = zero
        rdy5 = zero
        call sumgatN( abs(gami*Delt(itseq)  &
                      * Dy(1:numnp,1:3)),3,rdy1, numnp)
        call sumgatN( abs( y(1:numnp,1:3)),3,rdy2,numnp)
        call sumgatN( abs(gami*alfi*Delt(itseq) &
                      * Dy(1:numnp,4)),1,rdy4,numnp)
        call sumgatN( abs( y(1:numnp,4)),  1,rdy5,numnp)
        rmaxdyU = rdy1/rdy2
        rmaxdyP = rdy4/rdy5
	
!
!..... Signal to quit if delta is very small. look in itrdrv.f for the
!      completion of the hack.
!
	if( rmaxdyU .lt. dtol(1) .and. rmaxdyP .lt. dtol(2)) then
           istop = 1000
        endif

        if (numpe == 1) nshgt=nshg   ! global = this processor
!
!
!.... ----------------------->  Convergence  <-------------------------
!
!.... compute the maximum residual and the corresponding node number
!
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
           call MPI_REDUCE_SCATTER (resvec, resmax, irecvcount, &
                MPI_DOUBLE_PRECISION, MPI_MAX, INEWCOMM, ierr)
           if (resmax .eq. resvec(1) ) then
              mproc(1) = myrank
              nrsmax   = maxloc(rtmp)
           else
              nrsmax(1) = -1
              mproc(1)  = -1
           endif
           resvec = nrsmax
           call MPI_REDUCE_SCATTER (resvec, rvec, irecvcount, &
                MPI_DOUBLE_PRECISION, MPI_MAX, INEWCOMM, ierr)
           nrsmax = rvec(1)
           resvec = mproc
           call MPI_REDUCE_SCATTER (resvec, rvec, irecvcount, &
                MPI_DOUBLE_PRECISION, MPI_MAX, INEWCOMM, ierr)
           mproc = rvec(1)
      else
          resmax   = resmaxl
          nrsmax   = maxloc(rtmp)
          mproc(1) = 0
      endif
!
!.... correct the residuals
!
        if (loctim(itseq) .eq. 0) then
          resnrm = resnrm 
          resmax = resmax
        else
          resnrm = resnrm
          resmax = resmax
        endif
!
!.... approximate the number of entries
!
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
!     
!.... get the CPU-time
!
!AD        cputme = (second(0) - ttim(100))
        realsec=TMRC()
        cputme = (realsec - ttim(100))
!
!.... output the result
!
        if (numpe > 1) call MPI_BARRIER (INEWCOMM, ierr)
        
        if (myrank .eq. master) then
!
!.... results of continuity and momentum 
!
           
           write (*,2000) lstep+1, cputme, totres, jtotrs, rmaxdyU, &
                rmaxdyP,nrsmax, &
                mproc(1)+1, jresmx, int(statsflow(4)), &
                int(statsflow(1))
           write (ihist,2000) lstep+1, cputme, totres, jtotrs,  &
                rmaxdyU, rmaxdyP, nrsmax, &
                mproc(1)+1,jresmx,int(statsflow(4)), &
                int(statsflow(1))
           
           call flush(ihist)
        endif
        if(numpe>1) call MPI_BARRIER (INEWCOMM,ierr)

!$$$	ttim(68) = ttim(68) + tmr()

!
!.... return
!
        return
!
 1000   format(1p,i6,5e13.5)
 2000   format(1p,i6,e10.3,e10.3,2x,'(',i4,')',2x,e10.3,2x,e10.3, &
             2x,'<',i6,'-',i2,'|', &
             i4,'>', ' [', i4,' -',i4,']')
 3000   format(1p,i6,e10.3,e10.3,3x,'(',i4,')',3x,'<',i6,'-',i2,'|', &
             i4,'>', ' [', i4,' -',i4,' -',i4,']')

!
        end


        subroutine rstaticSclr (res, y, Dy, icomp)
!
!----------------------------------------------------------------------
!
! This subroutine calculates the statistics of the residual
!
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
        !include "auxmpi.h"
!
        dimension res(nshg)
        dimension rtmp(nshg)
        real*8    y(nshg,ndof),    Dy(nshg), nrm
        !integer tmrc
!
!.... compute max delta y
!
        rdy1 = zero
        rdy2 = zero
!
!.... normalize turbulence with molecular viscosity
!        
        if ( (icomp .eq. 6).and. (iRANS.eq.-1) ) then
           nrm = datmat(1,2,1)
        else 
           nrm = zero
        endif
        call sumgat( abs(gami*Delt(itseq)*Dy(:)),1,rdy1)
        call sumgat( abs( y(:,icomp)),1,rdy2)
        rmaxdyT = rdy1/(rdy2+nrm)
!
!.... compute the maximum residual and the corresponding node number
!
        rtmp = zero
        rtmp = rtmp + res**2 ! add temperature also
        call sumgat (rtmp, 1, resnrm)

        if (numpe == 1) nshgt=nshg ! global = this processor

        totres = resnrm / float(nshgt)
        totres = sqrt(totres)

!        if (mod(impl(1),100)/10 .eq. 0) then  !not solving flow
           if (myrank .eq. master) then
!     
!.... get the CPU-time
!
              realsec=TMRC()
              cputme = (realsec - ttim(100))

           print 802, lstep+1, cputme, totres, rmaxdyT, &
                      int(statssclr(1))
           write (ihist,802) lstep+1, cputme, totres,  &
                rmaxdyT,int(statssclr(1))
           
               call flush(ihist)
           endif
!        else 
!           if (myrank .eq. master) then
!              print 803, totres, rmaxdyT, int(statssclr(1))
!              write(ihist,803) totres, rmaxdyT, int(statssclr(1))
!           endif
!        endif

        return
        
 802    format(1p,i6,e10.3,e10.3,10X,e10.3,31X'[',i6,']')
 803    format(1p,16x,e10.3,10x,e10.3,31X,'[',i10,']')
    
        end



