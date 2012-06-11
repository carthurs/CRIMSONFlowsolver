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
      subroutine itrdrv (y,         ac,         
     &                   uold,      x,         
     &                   iBC,       BC,         
     &                   iper,      ilwork,     shp,       
     &                   shgl,      shpb,       shglb,
     &                   ifath,     velbar,     nsons ) 
c
c----------------------------------------------------------------------
c
c This iterative driver is the semi-discrete, predictor multi-corrector 
c algorithm. It contains the Hulbert Generalized Alpha method which
c is 2nd order accurate for Rho_inf from 0 to 1.  The method can be
c made  first-order accurate by setting Rho_inf=-1. It uses CGP and
c GMRES iterative solvers.
c
c working arrays:
c  y      (nshg,ndof)           : Y variables
c  x      (nshg,nsd)            : node coordinates
c  iBC    (nshg)                : BC codes
c  BC     (nshg,ndofBC)         : BC constraint parameters
c  iper   (nshg)                : periodicity table
c
c----------------------------------------------------------------------
c
      use pvsQbi     !gives us splag (the spmass at the end of this run 
      use specialBC !gives us itvn
      use timedata   !allows collection of time series
      use convolImpFlow !for Imp bc
      use convolRCRFlow !for RCR bc
      use convolCORFlow !for Coronary bc
      use incpBC        !for INCP bc
      use calcFlowPressure !to save history of flow and pressure of bc surfaces 
      use LagrangeMultipliers 
      use deformableWall
      use ResidualControl 
c      use measureWallDistance !to measure distances to segmented meshes
c      use readarrays !reads in uold and acold
      
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c

        
        real*8    y(nshg,ndof),              ac(nshg,ndof),           
     &            yold(nshg,ndof),           acold(nshg,ndof),
     &            u(nshg,nsd),               uold(nshg,nsd),
     &            x(numnp,nsd),              solinc(nshg,ndof),
     &            BC(nshg,ndofBC),           tf(nshg,ndof),
     &            xdist(numnp),
     &            xdnv(numnp,nsd)

c
        real*8    res(nshg,ndof)
c     
        real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        integer   rowp(nshg,nnz),        colm(nshg+1),
     &            iBC(nshg),             ilwork(nlwork),  
     &            iper(nshg),            ifuncs(6)

        integer stopjob
        character*10 cname2
        character*5  cname
c
c  stuff for dynamic model s.w.avg and wall model
c
        dimension ifath(numnp),    velbar(nfath,ndof),  nsons(nfath)

        dimension wallubar(2),walltot(2)
c     
c.... For linear solver Library
c
        integer eqnType, prjFlag, presPrjFlag, verbose
c
        real*8, allocatable, dimension(:,:) :: aperm,  atemp, atempS
        real*8, allocatable, dimension(:,:,:) :: apermS

        real*8, allocatable, dimension(:,:) :: lhsP, lhsK, lhsS
        real*8, allocatable, dimension(:,:) :: iBCs
        real*8, allocatable, dimension(:)   :: iBCd
        real*8  almit, alfit, gamit
c
        character*1024    servername
        character*20    fname1,fmt1
        character*20    fname2,fmt2,fnamer2
        integer         iarray(50) ! integers for headers
        character*20        license_f_name

        real*8 rerr(nshg,10), ybar(nshg,5) ,  dummyVar(nshg),
     &         uhess(nshg,27), gradu(nshg,9)
     
     
c
c find the machine name so that we set the license key properly

        license_f_name='license.dat'

        call SolverLicenseServer(servername)
c
c only master should be verbose
c

        if(numpe.gt.0 .and. myrank.ne.master)iverbose=0  
c

        inquire(file='xyzts.dat',exist=exts)
       
        if(exts) then
           
           open(unit=626,file='xyzts.dat',status='old')
           read(626,*) ntspts, freq, tolpt, iterat, varcod
           call sTD             ! sets data structures
           
           do jj=1,ntspts       ! read coordinate data where solution desired
              read(626,*) ptts(jj,1),ptts(jj,2),ptts(jj,3)
           enddo

           varts = zero

        endif
c
c.... open history and aerodynamic forces files
c
        if (myrank .eq. master) then
           open (unit=ihist,  file=fhist,  status='unknown')
           open (unit=iforce, file=fforce, status='unknown')
           open (unit=76, file="fort.76", status='unknown')
        endif
c
c.... initialize
c     
        ifuncs(:)  = 0              ! func. evaluation counter
        istep  = 0
        yold   = y
        acold  = ac

        rerr = zero
        ybar = y
c
c.... ---------------> initialize LesLib Library <---------------
c
c.... assign parameter values
c     
        do i = 1, 100
           numeqns(i) = i
        enddo
        nKvecs       = Kspace
        prjFlag      = iprjFlag
        presPrjFlag  = ipresPrjFlag
        verbose      = iverbose
c
c.... determine how many scalar equations we are going to need to solve
c
      nsolt=mod(impl(1),2)      ! 1 if solving temperature
      nsclrsol=nsolt+nsclr      ! total number of scalars solved At
                                ! some point we probably want to create
                                ! a map, considering stepseq(), to find
                                ! what is actually solved and only
                                ! dimension lhs to the appropriate
                                ! size. (see 1.6.1 and earlier for a
                                ! "failed" attempt at this).


      nsolflow=mod(impl(1),100)/10  ! 1 if solving flow
      
c
c.... Now, call lesNew routine to initialize
c     memory space
c
      call genadj(colm, rowp, icnt )  ! preprocess the adjacency list
      
      nnz_tot=icnt ! this is exactly the number of non-zero blocks on
                   ! this proc

      if (nsolflow.eq.1) then
         lesId   = numeqns(1)
         eqnType = 1
         nDofs   = 4
         call myfLesNew( lesId,   41994,
     &                 eqnType,
     &                 nDofs,          minIters,       maxIters,
     &                 nKvecs,         prjFlag,        nPrjs,
     &                 presPrjFlag,    nPresPrjs,      epstol(1),
     &                 prestol,        verbose,        statsflow,
     &                 nPermDims,      nTmpDims,      servername  )
         
         allocate (aperm(nshg,nPermDims))
         allocate (atemp(nshg,nTmpDims))
         allocate (lhsP(4,nnz_tot))
         allocate (lhsK(9,nnz_tot))

         call readLesRestart( lesId,  aperm, nshg, myrank, lstep,
     &                        nPermDims )

      else
         nPermDims = 0
         nTempDims = 0
      endif


      if(nsclrsol.gt.0) then
       do isolsc=1,nsclrsol
         lesId       = numeqns(isolsc+1)
         eqnType     = 2
         nDofs       = 1
         presPrjflag = 0        
         nPresPrjs   = 0       
         prjFlag     = 1
         indx=isolsc+2-nsolt ! complicated to keep epstol(2) for
                             ! temperature followed by scalars
         call myfLesNew( lesId,            41994,
     &                 eqnType,
     &                 nDofs,          minIters,       maxIters,
     &                 nKvecs,         prjFlag,        nPrjs,
     &                 presPrjFlag,    nPresPrjs,      epstol(indx),
     &                 prestol,        verbose,        statssclr,
     &                 nPermDimsS,     nTmpDimsS,   servername )
       enddo
c
c  Assume all scalars have the same size needs
c
       allocate (apermS(nshg,nPermDimsS,nsclrsol))
       allocate (atempS(nshg,nTmpDimsS))  !they can all share this
       allocate (lhsS(nnz_tot,nsclrsol))
c
c actually they could even share with atemp but leave that for later
c
      else
         nPermDimsS = 0
         nTmpDimsS  = 0
      endif
c
c...  prepare lumped mass if needed
c
      if((flmpr.ne.0).or.(flmpl.ne.0)) call genlmass(x, shp,shgl)
c
c.... -----------------> End of initialization <-----------------
c
c.....open the necessary files to gather time series
c
      lstep0 = lstep+1
c
c.... make two copies of iBC array to use a switch for INCP boundary condition
c    
      if (incp .gt. zero) then
         allocate(iBCd(nshg))
         iBCd = iBC
         if (numINCPSrfs .eq. one) then
            allocate(iBCs(1,nshg))
            iBCs(1,:) = iBC
            do i=1, nshg
               if(ndsurf(i).eq.nsrflistINCP(1)) then
                  iBCs(1,i)=0
               endif
            enddo
         elseif (numINCPSrfs .eq. two) then
            allocate(iBCs(3,nshg))
            do k=1, 3
               iBCs(k,:) = iBC
            enddo

            do i=1, nshg
               if(ndsurf(i).eq.nsrflistINCP(1)) then
                  iBCs(1,i)=0
                  iBCs(2,i)=0
               elseif(ndsurf(i).eq.nsrflistINCP(2)) then
                  iBCs(1,i)=0
                  iBCs(3,i)=0                              
               endif
            enddo
         endif            
      endif
c
c...initialize the initial condition for INCP BC
c   
      if(numINCPSrfs.gt.zero) then
         call calcINCPic(Delt(1), y, nsrflistINCP, numINCPSrfs)
      endif

      if (incp.gt.zero) then  
         if (numINCPSrfs .eq. one) then
            if(Qaorta(lstep+1,1) .lt. zero) then
               iBC = iBCs(1,:)   !systole
               inactive(1)=1
               INCPSwitch = 1
            endif
         elseif (numINCPSrfs .eq. two) then
            if(Qaorta(lstep+1,1) .lt. zero .and. 
     &         Qaorta(lstep+1,2) .lt. zero) then
               iBC = iBCs(1,:)   !systole
               inactive(1:2)=1
               INCPSwitch = 1
            endif
         endif
      endif
      
c
c
c.... satisfy the boundary conditions
c
      call itrBC (y, ac,  iBC, BC, iper, ilwork)

c     initialize distances
      xdist = zero
      xdnv = zero

c
c.... loop through the time sequences
c
      do 3000 itsq = 1, ntseq
         itseq = itsq

CAD         tcorecp1 = second(0)
CAD         tcorewc1 = second(-1)
c
c.... set up the time integration parameters
c         
         nstp   = nstep(itseq)
         nitr   = niter(itseq)
         LCtime = loctim(itseq)
         dtol(:)= deltol(itseq,:)

         call itrSetup ( y, acold )
c
c...initialize the coefficients for the impedance convolution,
c   which are functions of alphaf so need to do it after itrSetup
         if(numImpSrfs.gt.zero) then
            call calcImpConvCoef (numImpSrfs, ntimeptpT)
         endif
c
c...initialize the initial condition P(0)-RQ(0)-Pd(0) for RCR BC
c   need ndsurf so should be after initNABI
         if(numRCRSrfs.gt.zero) then
           call calcRCRic(y,nsrflistRCR,numRCRSrfs)
         endif
c
c...calculate area and initial pressure and flow for CalcSurfaces
c
         if(numCalcSrfs.gt.zero) then
           call calcCalcic(y,nsrflistCalc,numCalcSrfs)
         endif
c
c...initialize the initial condition for Coronary BC
c   need ndsurf so should be after initNABI
         if(numCORSrfs.gt.zero) then
            call calcCORic(y,nsrflistCOR,numCORSrfs)
         endif
c
c.... allocate LHS and RHS arrays required for the constrained surfaces
c   
         if(Lagrange.gt.zero) then
           call calcLagrangeic(nsrflistLagrange, 
     &        numLagrangeSrfs)
         endif   
         
c
c.... precompute the deformable wall stiffness matrix
c   
         if(ideformwall.eq.1) then
           call vlmwStTri(x,iBC,BC)
c           call solveWallProb(rowp,colm,ilwork,iBC,BC,iper)
         end if         
c
c  find the last solve of the flow in the step sequence so that we will
c         know when we are at/near end of step
c
c         ilast=0
         nitr=0  ! count number of flow solves in a step (# of iterations)
         do i=1,seqsize
            if(stepseq(i).eq.0) nitr=nitr+1
         enddo
c
c.... loop through the time steps
c
         istop=0
         rmub=datmat(1,2,1)
         if(rmutarget.gt.0) then
            rmue=rmutarget
         else
            rmue=datmat(1,2,1) ! keep constant
         endif
         do 2000 istp = 1, nstp

          if (incp.gt.zero) then  ! works only when there is one "INCP" srf
            if (numINCPSrfs .eq. one) then
               if(PLV(lstep+1,1) .gt. Paorta(lstep+1,1)) then
                  iBC = iBCs(1,:)   !systole
                  inactive(1)=1
                  INCPSwitch = 1
               elseif (INCPSwitch .gt. 0 .and. 
     &            Qaorta(lstep+1,1) .le. zero) then  
                  iBC = iBCs(1,:)   !systole
                  inactive(1)=1
                  INCPSwitch = 1
               else
                  iBC=iBCd   !diastole
                  inactive(1)=nsrflistINCP(1)
                  INCPSwitch = 0
               endif
            elseif (numINCPSrfs .eq. two) then
               if(PLV(lstep+1,1) .gt. Paorta(lstep+1,1) .and. 
     &            PLV(lstep+1,2) .gt. Paorta(lstep+1,2)) then
                  iBC = iBCs(1,:)   !systole
                  inactive(1:2)=1
                  INCPSwitch = 1
               elseif(PLV(lstep+1,1) .gt. Paorta(lstep+1,1) .and. 
     &            PLV(lstep+1,2) .le. Paorta(lstep+1,2)) then
                  iBC = iBCs(2,:)   !systole
                  inactive(1)=1
                  inactive(2)=nsrflistINCP(2)
                  INCPSwitch = 1
               elseif(PLV(lstep+1,1) .le. Paorta(lstep+1,1) .and. 
     &            PLV(lstep+1,2) .gt. Paorta(lstep+1,2)) then
                  iBC = iBCs(3,:)   !systole
                  inactive(2)=1
                  inactive(1)=nsrflistINCP(1)
                  INCPSwitch = 1
               elseif (INCPSwitch.gt.0 .and. Qaorta(lstep+1,1).le.zero 
     &            .and. Qaorta(lstep+1,2) .le. zero) then  
                  iBC = iBCs(1,:)   !systole
                  inactive(1:2)=1
                  INCPSwitch = 1
               elseif (INCPSwitch.gt.0 .and. Qaorta(lstep+1,1).le.zero 
     &            .and. Qaorta(lstep+1,2) .gt. zero) then  
                  iBC = iBCs(2,:)   !systole
                  inactive(1)=1
                  inactive(2)=nsrflistINCP(2)
                  INCPSwitch = 1
               elseif (INCPSwitch.gt.0 .and. Qaorta(lstep+1,1).gt.zero 
     &            .and. Qaorta(lstep+1,2) .le. zero) then  
                  iBC = iBCs(3,:)   !systole
                  inactive(2)=1
                  inactive(1)=nsrflistINCP(1)
                  INCPSwitch = 1
               else
                  iBC=iBCd   !diastole
                  inactive(1)=nsrflistINCP(1)
                  INCPSwitch = 0
               endif
            endif
          endif
 
          call rerun_check(stopjob)
           if(stopjob.ne.0) goto 2001

            xi=istp*1.0/nstp
            datmat(1,2,1)=rmub*(1.0-xi)+xi*rmue
c            write(*,*) "current mol. visc = ", datmat(1,2,1)
c.... if we have time varying boundary conditions update the values of BC.
c     these will be for time step n+1 so use lstep+1
c     
            if(itvn.gt.0) call BCint((lstep+1)*Delt(1), shp, shgl, 
     &                               shpb, shglb, x, BC, iBC)

c
c ... calc the pressure contribution that depends on the history for the imp BC
c     
            if(numImpSrfs.gt.0) call pHist(poldImp,QHistImp,ImpConvCoef,
     &                                          ntimeptpT,numImpSrfs)
c
c ... calc the pressure contribution that depends on the history for the RCR BC
c     
            if(numRCRSrfs.gt.0) then 
               call CalcHopRCR (Delt(itseq), lstep, numRCRSrfs) 
               call CalcRCRConvCoef(lstep,numRCRSrfs) 
               call pHist(poldRCR,QHistRCR,RCRConvCoef,
     &            nstep+nptsRCR,numRCRSrfs)
            endif
c
c ... calc the pressure contribution that depends on the history for the Coronary BC
c     
            if(numCORSrfs.gt.0) then 
               call CalcCORConvCoef(lstep,numCORSrfs) 
               call pHist(poldCOR, QHistCOR, CORConvCoef,
     &            nstep+nptsCOR,numCORSrfs)
               call CalcHopCOR (Delt(itseq), lstep, nsrflistCOR, 
     &                 numCORSrfs,y)
            endif
c
c.... calculate the coefficients
c
           if (numINCPSrfs .gt. 0) then
              call CalcINCPConvCoef(lstep, numINCPSrfs)
              call pHist(poldINCP, QHistINCP,
     &            INCPConvCoef, nstep+nptsINCP, numINCPSrfs)
              call CalcINCPCoef(Delt(itseq), lstep,
     &                nsrflistINCP, numINCPSrfs, y)
           endif
c
c Decay of scalars
c
           if(nsclr.gt.0 .and. tdecay.ne.1) then
              yold(:,6:ndof)=y(:,6:ndof)*tdecay
              BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*tdecay
           endif

           if(nosource.eq.1) BC(:,7:6+nsclr)= BC(:,7:6+nsclr)*0.8


            if(iLES.gt.0) then  !complicated stuff has moved to
                                        !routine below
               call lesmodels(yold,  acold,     shgl,      shp, 
     &                        iper,  ilwork,    rowp,      colm,
     &                        nsons, ifath,     x,   
     &                        iBC,   BC)

            
            endif

c.... set traction BCs for modeled walls
c
            if (itwmod.ne.0) then
               call asbwmod(yold,   acold,   x,      BC,     iBC,
     &                      iper,   ilwork,  ifath,  velbar)
            endif
c
c.... -----------------------> predictor phase <-----------------------
c
            call itrPredict(yold, y,   acold,  ac ,  uold,  u)
            call itrBC (y,  ac,  iBC,  BC,  iper,ilwork)

            if(nsolt.eq.1) then
               isclr=0
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
            endif
            do isclr=1,nsclr
               call itrBCSclr (y, ac,  iBC, BC, iper, ilwork)
            enddo
            
            iter=0
            ilss=0  ! this is a switch thrown on first solve of LS redistance
            
c           interface to compute distances to observed data      
            if (imeasdist.eq.1) then
               call ElmDist(u,x,xdist,xdnv)
c                  write(*,*) xdist
            end if
            
            do istepc=1,seqsize
               icode=stepseq(istepc)
               if(mod(icode,5).eq.0) then ! this is a solve
                  isolve=icode/10
                  if(icode.eq.0) then ! flow solve (encoded as 0)
c
                     iter   = iter+1
                     ifuncs(1)  = ifuncs(1) + 1
c     
                     Force(1) = zero
                     Force(2) = zero
                     Force(3) = zero
                     HFlux    = zero
                     lhs = 1 - min(1,mod(ifuncs(1)-1,LHSupd(1))) 

                     call SolFlow(y,             ac,        u,
     &                            yold,          acold,     uold,
     &                            x,             xdist,     xdnv,
     &                            iBC,           BC,        res,
     &                            nPermDims,     nTmpDims,  aperm,
     &                            atemp,         iper,          
     &                            ilwork,        shp,       shgl,
     &                            shpb,          shglb,     rowp,     
     &                            colm,          lhsK,      lhsP,
     &                            solinc,        rerr)
                  
                  else          ! scalar type solve
                     if (icode.eq.5) then ! Solve for Temperature
                                ! (encoded as (nsclr+1)*10)
                        isclr=0
                        ifuncs(2)  = ifuncs(2) + 1
                        j=1
                     else       ! solve a scalar  (encoded at isclr*10)
                        isclr=isolve  
                        ifuncs(isclr+2)  = ifuncs(isclr+2) + 1
                        j=isclr+nsolt
                        if((iLSet.eq.2).and.(ilss.eq.0)
     &                       .and.(isclr.eq.2)) then 
                           ilss=1 ! throw switch (once per step)
                           y(:,7)=y(:,6) ! redistance field initialized
                           ac(:,7)   = zero
                           call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                          ilwork)
c     
c....store the flow alpha, gamma parameter values and assigm them the 
c....Backward Euler parameters to solve the second levelset scalar
c     
                           alfit=alfi
                           gamit=gami
                           almit=almi
                           Deltt=Delt(1)
                           Dtglt=Dtgl
                           alfi = 1
                           gami = 1
                           almi = 1
c     Delt(1)= Deltt ! Give a pseudo time step
                           Dtgl = one / Delt(1)
                        endif  ! level set eq. 2
                     endif ! deciding between temperature and scalar

                     lhs = 1 - min(1,mod(ifuncs(isclr+2)-1,
     &                                   LHSupd(isclr+2))) 

                     call SolSclr(y,             ac,        u,
     &                            yold,          acold,     uold,
     &                            x,             iBC,
     &                            BC,            nPermDimsS,nTmpDimsS,  
     &                            apermS(1,1,j), atempS,    iper,
     &                            ilwork,        shp,       shgl,
     &                            shpb,          shglb,     rowp,     
     &                            colm,          lhsS(1,j), 
     &                            solinc(1,isclr+5))
                        
                        
                  endif         ! end of scalar type solve

               else ! this is an update  (mod did not equal zero)
c
c.... -----------------------> corrector phase <-----------------------
c               
                  iupdate=icode/10  ! what to update
                  if(icode.eq.1) then !update flow  
                     call itrCorrect ( y,    ac,    u,   solinc)
                     call itrBC (y,  ac,  iBC,  BC, iper, ilwork)
                  else  ! update scalar
                     isclr=iupdate  !unless
                     if(icode.eq.6) isclr=0
                     if(iRANS.lt.0)then  ! RANS
                        call itrCorrectSclrPos(y,ac,solinc(1,isclr+5))
                     else
                        call itrCorrectSclr (y, ac, solinc(1,isclr+5))
                     endif
                     if (ilset.eq.2 .and. isclr.eq.2)  then
                        if (ivconstraint .eq. 1) then
                           call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                          ilwork)
c                    
c ... applying the volume constraint on second level set scalar
c
                           call solvecon (y,    x,      iBC,  BC, 
     &                          iper, ilwork, shp,  shgl)
c
                        endif   ! end of volume constraint calculations
                     endif      ! end of redistance calculations
c                     
                     call itrBCSclr (  y,  ac,  iBC,  BC, iper,
     &                    ilwork)
                  endif
               endif         !end of switch between solve or update

               if(rescontrol .gt. 0) then
                  if (controlResidual .lt. ResCriteria .and. 
     &               CurrentIter .ge. MinNumIter) then
                     CurrentIter = 0
                     goto 1009
                     exit
                  endif      
               endif                
               
            enddo            ! loop over sequence in step
c     
c
c.... obtain the time average statistics
c
1009        if (ioform .eq. 2) then   

               call stsGetStats( y,      yold,     ac,     acold,
     &                           u,      uold,     
     &                           x,      xdist,    xdnv,
     &                           shp,    shgl,     shpb,   shglb,
     &                           iBC,    BC,       iper,   ilwork,
     &                           rowp,   colm,     lhsK,   lhsP )

            endif

c     
c  Find the solution at the end of the timestep and move it to old
c
c  
c ...First to reassign the parameters for the original time integrator scheme
c
            if((iLSet.eq.2).and.(ilss.eq.1)) then 
               alfi =alfit
               gami =gamit
               almi =almit 
               Delt(1)=Deltt
               Dtgl =Dtglt
            endif          
            call itrUpdate( yold,  acold,   uold,  y,    ac,   u)
            call itrBC (yold, acold,  iBC,  BC,  iper,ilwork)

            istep = istep + 1
            lstep = lstep + 1
c
c ... write out the solution
c
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
               call restar ('out ',  yold  ,ac)
               if(ideformwall.eq.1) then 
                  call write_displ(myrank, lstep, nshg, 3, uold ) 
                  if (imeasdist.eq.1) then
                     call write_distl(myrank, lstep, numnp, 1, xdist ) ! should use nshg or numnp?
                  end if
               end if
            endif
c
c ... update the flow history for the INCP convolution
c
            if(numINCPSrfs.gt.zero) then
               call UpdHistConv(y,nsrflistINCP,numINCPSrfs)
               call UpdHeartModel(Delt(itseq), 
     &            y,nsrflistINCP,numINCPSrfs, lstep)
            endif
c 
c ... update the flow history for the impedance convolution, filter it and write it out
c    
            if(numImpSrfs.gt.zero) then
               call UpdHistConv(y,nsrflistImp,numImpSrfs) !uses Delt(1)
            endif

c 
c ... update the flow history for the RCR convolution
c    
            if(numRCRSrfs.gt.zero) then
               call UpdHistConv(y,nsrflistRCR,numRCRSrfs) !uses lstep
               call UpdRCR(y,nsrflistRCR,numRCRSrfs)
            endif
c 
c ... update the flow history for the Coronary convolution
c    
            if(numCORSrfs.gt.zero) then
               call UpdHistConv(y,nsrflistCOR,numCORSrfs) !uses lstep
               call UpdHistPlvConv(y, Delt(itseq), lstep, 
     &            nsrflistCOR, numCORSrfs) 
            endif
c 
c ... update the flow history for the CalcSurfaces
c    
            if(numCalcSrfs.gt.zero) then
               call Updcalc(y,nsrflistCalc,numCalcSrfs)
            endif
c
c.... calculate the values of constraint functions and write Lagrange Multipliers
c
           if (Lagrange .gt. 0) then
              call UpdateLagrangeCoef(y, colm, rowp, nsrflistLagrange,
     &           numLagrangeSrfs)
           endif               
           
c
c.... compute the consistent boundary flux
c
            if(abs(itwmod).ne.1)
     &         call Bflux ( yold,      acold,      uold,     
     &                      x,         xdist,      xdnv,
     &                      shp,       shgl,       shpb,   
     &                      shglb,     ilwork,     iBC,
     &                      BC,        iper)


c...  dump TIME SERIES
            
            if (exts) then
               
               if (mod(lstep-1,freq).eq.0) then
                  
                  do jj = 1, ntspts
                     
                     if (numpe > 1) then
                        
                        soln = varts(jj)
                        asoln = abs(soln)
                        
c     if(jj.eq.24) then
c     write(*,*) soln
c     write(*,*)"and..."
c     endif
                        
                        if (asoln.ne.zero) then
                           sgn = soln/asoln
                        else
                           sgn = 1
                        endif
                        
                        call MPI_ALLREDUCE ( asoln, asolng, 1, 
     &                       MPI_DOUBLE_PRECISION, MPI_MAX,
     &                       MPI_COMM_WORLD,ierr)
                        varts(jj) = sgn * asolng
                        
                     endif
                     
                     if (myrank.eq.zero) then
                        ifile = 1000+jj
                        write(ifile,555) varts(jj)
                        call flush(ifile)
                     endif
                     
                  enddo
                  
                  
                  varts = zero  ! reset the array for next step
                  
 555              format(e18.11)
                  
               endif
               
            endif


c
c.... update and the aerodynamic forces
c
            call forces ( yold,  ilwork )
            
            if((irscale.ge.0).or.(itwmod.gt.0)) 
     &           call getvel (yold,     ilwork, iBC,
     &                        nsons,    ifath, velbar)

            if((irscale.ge.0).and.(myrank.eq.master)) then
               call genscale(yold,       x,       iper, 
     &                       iBC,     ifath,   velbar,
     &                       nsons)
            endif
c
c....  print out results.
c
            ntoutv=max(ntout,100)   ! velb is not needed so often
            if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
               if( (mod(lstep, ntoutv) .eq. 0) .and.
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.(iLES.gt.0))))
     &              call rwvelb  ('out ',  velbar  ,ifail)
            endif
c
c.... end of the NSTEP and NTSEQ loops
c
c
c.... -------------------> error calculation  <-----------------
c 
            if(ierrcalc.eq.1 .or. ioybar.eq.1) then
c$$$c
c$$$c compute average
c$$$c
c$$$               tfact=one/istep
c$$$               ybar =tfact*yold + (one-tfact)*ybar

c compute average
c ybar(:,1) - ybar(:,3) is average velocity components
c ybar(:,4) is average pressure
c ybar(:,5) is average speed
c averaging procedure justified only for identical time step sizes
c istep is number of time step
c
               tfact=one/istep

c ybar to contain the averaged ((u,v,w),p)-field
c and speed average, i.e sqrt(u^2+v^2+w^2)

               ybar(:,1) = tfact*yold(:,1) + (one-tfact)*ybar(:,1)
               ybar(:,2) = tfact*yold(:,2) + (one-tfact)*ybar(:,2)
               ybar(:,3) = tfact*yold(:,3) + (one-tfact)*ybar(:,3)
               ybar(:,4) = tfact*yold(:,4) + (one-tfact)*ybar(:,4)
c    
               dummyVar  = sqrt(yold(:,1)**2+yold(:,2)**2+yold(:,3)**2)

               if (istep .eq. 1) then
                  ybar(:,5) = dummyVar
               else
                  ybar(:,5) = tfact*dummyVar + (one-tfact)*ybar(:,5)
               endif
c
c compute rms
c
               rerr(:, 7)=rerr(:, 7)+(yold(:,1)-ybar(:,1))**2
               rerr(:, 8)=rerr(:, 8)+(yold(:,2)-ybar(:,2))**2
               rerr(:, 9)=rerr(:, 9)+(yold(:,3)-ybar(:,3))**2
               rerr(:,10)=rerr(:,10)+(yold(:,4)-ybar(:,4))**2
            endif
            
            if(istop.eq.1000) exit ! stop when delta small (see rstatic)
 2000    continue
 2001    continue
        

CAD         tcorecp2 = second(0)
CAD         tcorewc2 = second(-1)
         
CAD         write(6,*) 'T(core) cpu-wallclock = ',tcorecp2-tcorecp1,
CAD     &                                        tcorewc2-tcorewc1

 3000 continue
 
c
c.... ---------------------->  Post Processing  <----------------------
c
c.... print out the last step
c
      if ((irs .ge. 1) .and. ((mod(lstep, ntout) .ne. 0) .or.
     &     (nstp .eq. 0))) then
         if(
     &              ((irscale.ge.0).or.(itwmod.gt.0) .or. 
     &              ((nsonmax.eq.1).and.iLES.gt.0)))
     &              call rwvelb  ('out ',  velbar  ,ifail)
         call restar ('out ',  yold  ,ac)
         if(ideformwall.eq.1) then
            call write_displ(myrank, lstep, nshg, 3, u ) 
            if (imeasdist.eq.1) then
               call ElmDist(u,x,xdist)
               call write_distl(myrank, lstep, numnp,1,xdist)  
            end if
         end if
      endif


         lesId   = numeqns(1)
         call saveLesRestart( lesId,  aperm , nshg, myrank, lstep,
     &                        nPermDims )


      if(ierrcalc.eq.1) then

c
c.....smooth the error indicators
c
        do i=1,ierrsmooth
            call errsmooth( rerr, x, iper, ilwork, shp, shgl, iBC )
        end do
c
c.... open the output file
c
           iqoldsiz=nshg*ndof*2
           call write_error(myrank, lstep, nshg, 10, rerr ) 
                         
                         
      endif

      if(ioybar.eq.1) then

         itmp = 1
         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
         write (fmt2,"('(''restart.'',i',i1,',1x)')") itmp
         write (fname2,fmt2) lstep

         fname2 = trim(fname2) // cname(myrank+1)
c
c.... open  files
c
         call openfile(  fname2,  'append?', irstin )

         fnamer2 = 'ybar'
         isize = nshg*5
         nitems = 3
         iarray(1) = nshg
         iarray(2) = 5
         iarray(3) = lstep
         call writeheader(irstin, fnamer2,iarray, nitems, isize,
     &        'double', iotype )

         nitems = nshg*5
         call writedatablock(irstin, fnamer2,ybar, nitems,
     &        'double', iotype)

         call closefile( irstin, "append" )

      endif


      if ( ( ihessian .eq. 1 ) .and. ( numpe < 2 )  )then

          uhess = zero
          gradu = zero
          tf = zero

          do ku=1,nshg
c           tf(ku,1) = x(ku,1)**2+2*x(ku,1)*x(ku,2)
            tf(ku,1) = x(ku,1)**3
          end do

          call hessian( yold, x,     shp,  shgl,   iBC, 
     &                  shpb, shglb, iper, ilwork, uhess, gradu )

          call write_hessian( uhess, gradu, nshg )

      endif
c
c.... close history and aerodynamic forces files
c
      if (myrank .eq. master) then
         close (ihist)
         close (iforce)
      endif
 5    format(1X,F15.10,3X,F15.10,3X,F15.10,3X,F15.10)
 444  format(6(2x,e14.7))
c
c.... end
c
      if(nsolflow.eq.1) then
         deallocate (lhsK)
         deallocate (lhsP)
         deallocate (aperm)
         deallocate (atemp)
         if (Lagrange .gt. 0) then
            deallocate (lhsLagL)
         endif         
      endif
      if(nsclrsol.gt.0) then
         deallocate (lhsS)
         deallocate (apermS)
         deallocate (atempS)
      endif
      
      if(iabc==1) deallocate(acs)

      return
      end
      
      subroutine lesmodels(y,     ac,        shgl,      shp, 
     &                     iper,  ilwork,    rowp,      colm,    
     &                     nsons, ifath,     x,   
     &                     iBC,   BC)
      
      include "common.h"

      real*8    y(nshg,ndof),              ac(nshg,ndof),           
     &            x(numnp,nsd),
     &            BC(nshg,ndofBC)
      real*8    shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT)

c
      integer   rowp(nshg,nnz),         colm(nshg+1),
     &            iBC(nshg),
     &            ilwork(nlwork),
     &            iper(nshg)
      dimension ifath(numnp),    nsons(nfath)

      real*8, allocatable, dimension(:) :: fwr2,fwr3,fwr4
      real*8, allocatable, dimension(:) :: stabdis,cdelsq1
      real*8, allocatable, dimension(:,:) :: xavegt, xavegt2,xavegt3

      if( (iLES.gt.1) )   then ! Allocate Stuff for advanced LES models
         allocate (fwr2(nshg))
         allocate (fwr3(nshg))
         allocate (fwr4(nshg))
         allocate (xavegt(nfath,12))
         allocate (xavegt2(nfath,12))
         allocate (xavegt3(nfath,12))
         allocate (stabdis(nfath))
      endif

c.... get dynamic model coefficient
c
      ilesmod=iLES/10  
c
c digit bit set filter rule, 10 bit set model
c
      if (ilesmod.eq.0) then    ! 0 < iLES< 10 => dyn. model calculated
                                ! at nodes based on discrete filtering


         if(isubmod.eq.2) then
            call SUPGdis(y,      ac,        shgl,      shp, 
     &                   iper,   ilwork,    
     &                   nsons,  ifath,     x,   
     &                   iBC,    BC, stabdis, xavegt3)
         endif

         if( ((isubmod.eq.0).or.(isubmod.eq.2)))then ! If no
                                                     ! sub-model
                                                     ! or SUPG
                                                     ! model wanted

            if(i2filt.eq.0)then ! If simple filter
              
               if(modlstats .eq. 0) then ! If no model stats wanted
                  call getdmc (y,       shgl,      shp, 
     &                         iper,       ilwork,    nsons,
     &                         ifath,      x)
               else             ! else get model stats 
                  call stdfdmc (y,       shgl,      shp, 
     &                          iper,       ilwork,    nsons,
     &                          ifath,      x)
               endif            ! end of stats if statement  

            else                ! else if twice filtering

               call widefdmc(y,       shgl,      shp, 
     &                       iper,       ilwork,    nsons,
     &                       ifath,      x)

               
            endif               ! end of simple filter if statement

         endif                  ! end of SUPG or no sub-model if statement


         if( (isubmod.eq.1) ) then ! If DFWR sub-model wanted
            call cdelBHsq (y,       shgl,      shp, 
     &                     iper,       ilwork,    nsons,
     &                     ifath,      x,         cdelsq1)
            call FiltRat (y,       shgl,      shp, 
     &                    iper,       ilwork,    nsons,
     &                    ifath,      x,         cdelsq1,
     &                    fwr4,       fwr3)

            
            if (i2filt.eq.0) then ! If simple filter wanted
               call DFWRsfdmc(y,       shgl,      shp, 
     &                        iper,       ilwork,    nsons,
     &                        ifath,      x,         fwr2, fwr3) 
            else                ! else if twice filtering wanted 
               call DFWRwfdmc(y,       shgl,      shp, 
     &                        iper,       ilwork,    nsons,
     &                        ifath,      x,         fwr4, fwr4) 
            endif               ! end of simple filter if statement
             
         endif                  ! end of DFWR sub-model if statement

         if( (isubmod.eq.2) )then ! If SUPG sub-model wanted
            call dmcSUPG (y,           ac,         shgl,      
     &                    shp,         iper,       ilwork,    
     &                    nsons,       ifath,      x,
     &                    iBC,    BC,  rowp,       colm,
     &                    xavegt2,    stabdis)
         endif

         if(idis.eq.1)then      ! If SUPG/Model dissipation wanted
            call ediss (y,        ac,      shgl,      
     &                  shp,      iper,       ilwork,    
     &                  nsons,    ifath,      x,
     &                  iBC,      BC,  xavegt)
         endif

      endif                     ! end of ilesmod
      
      if (ilesmod .eq. 1) then  ! 10 < iLES < 20 => dynamic-mixed
                                ! at nodes based on discrete filtering
         call bardmc (y,       shgl,      shp, 
     &                iper,    ilwork,    
     &                nsons,   ifath,     x) 
      endif
      
      if (ilesmod .eq. 2) then  ! 20 < iLES < 30 => dynamic at quad
                                ! pts based on lumped projection filt. 

         if(isubmod.eq.0)then
            call projdmc (y,       shgl,      shp, 
     &                    iper,       ilwork,    x) 
         else
            call cpjdmcnoi (y,      shgl,      shp, 
     &                      iper,   ilwork,       x,
     &                      rowp,   colm, 
     &                      iBC,    BC)
         endif

      endif

      if( (iLES.gt.1) )   then ! Deallocate Stuff for advanced LES models
         deallocate (fwr2)
         deallocate (fwr3)
         deallocate (fwr4)
         deallocate (xavegt)
         deallocate (xavegt2)
         deallocate (xavegt3)
         deallocate (stabdis)
      endif
      return
      end

c
c...initialize the coefficients for the impedance convolution
c
      subroutine CalcImpConvCoef (numISrfs, numTpoints)

      use convolImpFlow !uses flow history and impedance for convolution
      
      include "common.h" !for alfi
      
      integer numISrfs, numTpoints      

      allocate (ConvCoef(numTpoints+2,3)) !same time discret. for all imp. BC
      do j=1,numTpoints+2
         ConvCoef(j,:)=0.5/numTpoints !dt/2 divided by period T=N*dt
         ConvCoef(j,1)=ConvCoef(j,1)*(1.0-alfi)*(1.0-alfi)
         ConvCoef(j,2)=ConvCoef(j,2)*(1.0+2*alfi*(1.0-alfi))
         ConvCoef(j,3)=ConvCoef(j,3)*alfi*alfi
      enddo
      ConvCoef(1,2)=zero
      ConvCoef(1,3)=zero
      ConvCoef(2,3)=zero
      ConvCoef(numTpoints+1,1)=zero
      ConvCoef(numTpoints+2,2)=zero
      ConvCoef(numTpoints+2,1)=zero  
c
c...calculate the coefficients for the impedance convolution
c 
      allocate (ImpConvCoef(numTpoints+2,numISrfs))

c..coefficients below assume Q linear in time step, Z constant
c            do j=3,numTpoints
c                ImpConvCoef(j,:) = ValueListImp(j-1,:)*ConvCoef(j,3)
c     &                             + ValueListImp(j,:)*ConvCoef(j,2)    
c     &                             + ValueListImp(j+1,:)*ConvCoef(j,1)  
c            enddo
c            ImpConvCoef(1,:) = ValueListImp(2,:)*ConvCoef(1,1)
c            ImpConvCoef(2,:) = ValueListImp(2,:)*ConvCoef(2,2)    
c     &                       + ValueListImp(3,:)*ConvCoef(2,1)
c            ImpConvCoef(numTpoints+1,:) =
c     &           ValueListImp(numTpoints,:)*ConvCoef(numTpoints+1,3)
c     &         + ValueListImp(numTpoints+1,:)*ConvCoef(numTpoints+1,2) 
c            ImpConvCoef(numTpoints+2,:) = 
c     &           ValueListImp(numTpoints+1,:)*ConvCoef(numTpoints+2,3)

c..try easiest convolution Q and Z constant per time step
      do j=3,numTpoints+1
         ImpConvCoef(j,:) = ValueListImp(j-1,:)/numTpoints
      enddo
      ImpConvCoef(1,:) =zero
      ImpConvCoef(2,:) =zero
      ImpConvCoef(numTpoints+2,:) = 
     &           ValueListImp(numTpoints+1,:)/numTpoints
c compensate for yalpha passed not y in Elmgmr()
      ImpConvCoef(numTpoints+1,:)= ImpConvCoef(numTpoints+1,:)
     &                  - ImpConvCoef(numTpoints+2,:)*(1.0-alfi)/alfi 
      ImpConvCoef(numTpoints+2,:)= ImpConvCoef(numTpoints+2,:)/alfi 
      return
      end

c 
c ... update the flow rate history for the impedance convolution, filter it and write it out
c    
      subroutine UpdHistConv(y,nsrfIdList,numSrfs)
      
      use convolImpFlow !brings ntimeptpT, QHistImp, QHistTry, QHistTryF, numImpSrfs
      use convolRCRFlow !brings QHistRCR, numRCRSrfs
      use convolCORFlow 
      use incpBC
c
      include "common.h" !needed?
      include "mpif.h" !needed?
      
      integer   nsrfIdList(0:MAXSURF), numSrfs
      real*8    y(nshg,3) !velocity at time n+1   
      real*8    NewQ(0:MAXSURF)

      call GetFlowQ(NewQ,y,nsrfIdList,numSrfs) !new flow at time n+1
c
c... for imp BC: shift QHist, add new constribution, filter and write out
c      
      if(numImpSrfs.gt.zero .and. nsrfIdList(1).eq.nsrflistImp(1)) then
         do j=1, ntimeptpT
            QHistImp(j,1:numSrfs)=QHistImp(j+1,1:numSrfs)
         enddo
         QHistImp(ntimeptpT+1,1:numSrfs) = NewQ(1:numSrfs)
         QHistImp(1,:)=zero

c
c....filter the flow rate history
c
c         cutfreq = 10           !hardcoded cutting frequency of the filter
c         do j=1, ntimeptpT
c            QHistTry(j,:)=QHistImp(j+1,:)
c         enddo
c         call Filter(QHistTryF,QHistTry,ntimeptpT,Delt(1),cutfreq)
c         QHistImp(1,:)=zero
c         do j=1, ntimeptpT
c            QHistImp(j+1,:)=QHistTryF(j,:)
c         enddo
c
c.... write out the new history of flow rates to Qhistor.dat
c      
         if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &               (myrank .eq. zero)) then
            open(unit=816, file='Qhistor.dat',status='replace')
            write(816,*) ntimeptpT
            do j=1,ntimeptpT+1
               write(816,*) (QHistImp(j,n),n=1, numSrfs)
            enddo
            close(816)
         endif
      endif 

c
c... for RCR bc just add the new contribution
c
      if(numRCRSrfs.gt.zero .and. nsrfIdList(1).eq.nsrflistRCR(1)) then
         QHistRCR(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      endif      
c
c... for Coronary bc just add the new contribution
c
      if(numCORSrfs.gt.zero.and.nsrfIdList(1).eq.nsrflistCOR(1)) then
         QHistCOR(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      endif      
c
c... for INCP bc just add the new contribution
c
      if(numINCPSrfs.gt.zero.and.nsrfIdList(1).eq.nsrflistINCP(1)) then
         QHistINCP(lstep+1,1:numSrfs) = NewQ(1:numSrfs)  
      endif 
      
      return
      end

c
c...calculate the time varying coefficients for the RCR convolution
c
      subroutine CalcRCRConvCoef (stepn, numSrfs)

      use convolRCRFlow !brings in ValueListRCR, dtRCR
      
      include "common.h" !brings alfi
      
      integer numSrfs, stepn    

      RCRConvCoef = zero
      if (stepn .eq. 0) then
        RCRConvCoef(1,:) = ValueListRCR(1,:)*(1.0-alfi) +
     &   ValueListRCR(3,:)*(-alfi + 1.0 + 1/dtRCR(:) 
     &     - exp(-alfi*dtRCR(:))*(1 + 1/dtRCR(:)))
        RCRConvCoef(2,:) = ValueListRCR(1,:)*alfi 
     &     + ValueListRCR(3,:)
     &     *(alfi - 1/dtRCR(:) + exp(-alfi*dtRCR(:))/dtRCR(:))
      endif
      if (stepn .ge. 1) then
        RCRConvCoef(1,:) =-ValueListRCR(3,:)*exp(-dtRCR(:)*(stepn+alfi))
     &        *(1 + (1 - exp(dtRCR(:)))/dtRCR(:))
        RCRConvCoef(stepn+1,:) = ValueListRCR(1,:)*(1-alfi) 
     &     - ValueListRCR(3,:)*(alfi - 1 - 1/dtRCR(:) 
     &     + exp(-alfi*dtRCR(:))/dtRCR(:)*(2 - exp(-dtRCR(:))))
        RCRConvCoef(stepn+2,:) = ValueListRCR(1,:)*alfi 
     &     + ValueListRCR(3,:)
     &     *(alfi - 1/dtRCR(:) + exp(-alfi*dtRCR(:))/dtRCR(:))
      endif
      if (stepn .ge. 2) then
        do j=2,stepn
         RCRConvCoef(j,:) = ValueListRCR(3,:)/dtRCR(:)*
     &        exp(-dtRCR(:)*(stepn + alfi + 2 - j))*
     &        (1 - exp(dtRCR(:)))**2
        enddo
      endif

c compensate for yalpha passed not y in Elmgmr()
      RCRConvCoef(stepn+1,:)= RCRConvCoef(stepn+1,:)
     &                  - RCRConvCoef(stepn+2,:)*(1.0-alfi)/alfi 
      RCRConvCoef(stepn+2,:)= RCRConvCoef(stepn+2,:)/alfi 

      return
      end

c
c...calculate the time dependent H operator for the RCR convolution
c
      subroutine CalcHopRCR (timestepRCR, stepn, numSrfs)

      use convolRCRFlow !brings in HopRCR, dtRCR

      include "common.h"
      include "mpif.h" !needed?

      integer numSrfs, stepn      
      real*8  PdistCur(0:MAXSURF), timestepRCR
      
      HopRCR=zero
      call RCRint(timestepRCR*(stepn + alfi),PdistCur)
      HopRCR(1:numSrfs) = RCRic(1:numSrfs) 
     &     *exp(-dtRCR(1:numSrfs)*(stepn + alfi)) + PdistCur(1:numSrfs)
      return
      end

c
c.... This subroutine writes FlowHist.dat and PressHist.dat files
c
      subroutine UpdRCR(y, srfIDList, numSrfs)

      use convolRCRFlow 

      include "common.h"
c      
      real*8   y(nshg, ndof), NewP(0:MAXSURF)
      integer  srfIDList(0:MAXSURF),  numSrfs
      
      call integrScalar(NewP,y(:,4),srfIdList,numSrfs)
         PHistRCR(lstep+1,1:numSrfs)=NewP(1:numSrfs)/RCRArea(1:numSrfs)
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &   (myrank .eq. zero)) then
         call OutputDataFile(QHistRCR(1:lstep+1,:),lstep+1,numSrfs,
     &      'QHistRCR.dat',870)
         call OutputDataFile(PHistRCR(1:lstep+1,:),lstep+1,numSrfs,
     &      'PHistRCR.dat',871)
      endif 

      return
      end
      
c
c...calculate the time varying coefficients of pressure 
c...for the Coronary convolution
c
      subroutine CalcCORConvCoef (stepn, numSrfs)

      use convolCORFlow !brings in dtCOR, COR, CoefCOR

      include "common.h"
      include "mpif.h"

      integer numSrfs, stepn    

      CORConvCoef = zero
      if (stepn. eq. 0) then
        CORConvCoef(1,:)=real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)/alfi*
     &              ((-one+alfi*dtCOR(1,:))*exp(dtCOR(1,:)*alfi)+one)
     &                  - CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &               ((-one+alfi*dtCOR(2,:))*exp(dtCOR(2,:)*alfi)+one))
        CORConvCoef(2,:)= real(CoefCOR(5,:) +
     &                    CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)/alfi*
     &                  (exp(dtCOR(1,:)*alfi)-(one+alfi*dtCOR(1,:)))-
     &                    CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &                  (exp(dtCOR(2,:)*alfi)-(one+alfi*dtCOR(2,:))))
      endif
      if (stepn. ge. 1) then  
         CORConvCoef(1,:) =real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)*
     &                  (exp(dtCOR(1,:)*(stepn+alfi-one))-
     &                (one-dtCOR(1,:))*exp(dtCOR(1,:)*(stepn+alfi)))
     &                 -CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)*
     &                 (exp(dtCOR(2,:)*(stepn+alfi-one))-
     &                (one-dtCOR(2,:))*exp(dtCOR(2,:)*(stepn+alfi))))
         CORConvCoef(stepn+1,:) = real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)
     &                        /alfi*(alfi*exp(dtCOR(1,:)*(alfi+one))-
     &                        (alfi+one)*exp(dtCOR(1,:)*(alfi))+one)-
     &                        CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &                        (alfi*exp(dtCOR(2,:)*(alfi+1))-
     &                        (alfi+one)*exp(dtCOR(2,:)*(alfi))+one))
         CORConvCoef(stepn+2,:) = real(CoefCOR(5,:)+
     &                        CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)/alfi*
     &                     (exp(dtCOR(1,:)*alfi)-one-alfi*dtCOR(1,:))-
     &                        CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &                      (exp(dtCOR(2,:)*alfi)-one-alfi*dtCOR(2,:)))
      endif
      if (stepn. ge. 2) then
      do j=2,stepn
         CORConvCoef(j,:) = real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)*
     &                      (exp(dtCOR(1,:)*(stepn+alfi-j))-
     &                      two*exp(dtCOR(1,:)*(stepn+alfi-j+one))+
     &                     exp(dtCOR(1,:)*(stepn+alfi-j+two)))-
     &                      CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)*
     &                      (exp(dtCOR(2,:)*(stepn+alfi-j))-
     &                      two*exp(dtCOR(2,:)*(stepn+alfi-j+one))+
     &                      exp(dtCOR(2,:)*(stepn+alfi-j+two))))
      enddo
      endif
      
      return
      end

c
c...calculate the time varying coefficients of left ventricular pressure
c...for the Coronary convolution
c...need to do for t=0 and 1
c

      subroutine CalcCORPlvConvCoef (stepn, numSrfs)

      use convolCORFlow !brings in dtCOR, COR, CoefCOR

      include "common.h"

      integer numSrfs, stepn    
      
      CORPlvConvCoef = zero
      if (stepn. eq. 0) then
         CORPlvConvCoef(1,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)/alfi*
     &             ((-one+alfi*dtCOR(1,:))*exp(dtCOR(1,:)*alfi)+one)
     &                 - CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &             ((-one+alfi*dtCOR(2,:))*exp(dtCOR(2,:)*alfi)+one))
         CORPlvConvCoef(2,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)/alfi*
     &                (exp(dtCOR(1,:)*alfi)-(one+alfi*dtCOR(1,:)))-
     &                   CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &                 (exp(dtCOR(2,:)*alfi)-(one+alfi*dtCOR(2,:))))
      endif
      if (stepn. ge. 1) then
       CORPlvConvCoef(1,:) = real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)*
     &                  (exp(dtCOR(1,:)*(stepn+alfi-one))-
     &             (one-dtCOR(1,:))*exp(dtCOR(1,:)*(stepn+alfi)))
     &                -CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)*
     &                 (exp(dtCOR(2,:)*(stepn+alfi-one))-
     &               (one-dtCOR(2,:))*exp(dtCOR(2,:)*(stepn+alfi))))
      CORPlvConvCoef(stepn+1,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)
     &                       /alfi*(alfi*exp(dtCOR(1,:)*(alfi+one))-
     &                       (alfi+one)*exp(dtCOR(1,:)*(alfi))+one)-
     &                        CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &                        (alfi*exp(dtCOR(2,:)*(alfi+1))-
     &                        (alfi+one)*exp(dtCOR(2,:)*(alfi))+one))
       CORPlvConvCoef(stepn+2,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)
     &                /alfi*(exp(dtCOR(1,:)*alfi)-one-alfi*dtCOR(1,:))-
     &                        CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi*
     &                     (exp(dtCOR(2,:)*alfi)-one-alfi*dtCOR(2,:)))
      endif
      if (stepn. ge. 2) then
      do j=2,stepn
         CORPlvConvCoef(j,:) = real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)*
     &                      (exp(dtCOR(1,:)*(stepn+alfi-j))-
     &                      two*exp(dtCOR(1,:)*(stepn+alfi-j+one))+
     &                     exp(dtCOR(1,:)*(stepn+alfi-j+two)))-
     &                      CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)*
     &                      (exp(dtCOR(2,:)*(stepn+alfi-j))-
     &                      two*exp(dtCOR(2,:)*(stepn+alfi-j+one))+
     &                      exp(dtCOR(2,:)*(stepn+alfi-j+two))))
      enddo
      endif

      return
      end

c
c...calculate the time dependent H operator for the Coronary convolution
c
      subroutine CalcHopCOR (timestepCOR, stepn, srfIdList, numSrfs, y)

      use convolCORFlow !brings in HopCOR, dtCOR, COR, CoefCOR
      
      include "common.h" !needed?
      include "mpif.h" !needed?

      integer   srfIdList(0:MAXSURF), numSrfs, stepn
      real*8    y(nshg,4), timestepCOR, PlvistNext(0:MAXSURF)
      real*8    CoupleArea(0:MAXSURF), POnly(nshg)
      
      HopCOR=zero
      PlvistNext=zero
      
      call CalcCORPlvConvCoef (stepn, numSrfs)
      call pHist(plvoldCOR, PlvHistCOR, CORPlvConvCoef,
     &    nstep+nptsCOR,numSrfs) 
      call CORint(timestepCOR*(stepn + alfi),PlvistNext,1)
      POnly(:)=y(:,4) ! pressure
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get initial pressure integral
      PHistCOR(stepn+1,1:numSrfs) = CoupleArea(1:numSrfs)
     &   /CORArea(1:numSrfs)
      HopCOR(1:numSrfs) =real(plvoldCOR(1:numSrfs)+ 
     &      CORic(1,1:numSrfs)*exp(dtCOR(1,1:numSrfs)*(stepn+alfi))-
     &      CORic(2,1:numSrfs)*exp(dtCOR(2,1:numSrfs)*(stepn+alfi))+ 
     &      CorPlvConvCoef(stepn+2,1:numSrfs)*PlvistNext(1:numSrfs))

      return
      end
c
c...update time history of left ventricular pressure
c
      subroutine UpdHistPlvConv(y,timestepCOR,stepn,srfIdList,numSrfs) 
      
      use convolCORFlow 
      
      include "common.h" 
      include "mpif.h" 

      integer   srfIdList(0:MAXSURF), numSrfs, stepn
      real*8    timestepCOR, PlvistNext(0:MAXSURF)
      real*8    y(nshg,4)
      real*8    CoupleArea(0:MAXSURF), POnly(nshg)
      
      PlvistNext=zero      
      call CORint(timestepCOR*stepn,PlvistNext,zero)
      PlvHistCOR(stepn+1,1:numSrfs)=PlvistNext(1:numSrfs)   
      POnly(:)=y(:,4) ! pressure
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get initial pressure integral
      PHistCOR(stepn+1,1:numSrfs) = CoupleArea(1:numSrfs)
     &   /CORArea(1:numSrfs)

      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &   (myrank .eq. zero)) then
         call OutputDataFile(QHistCOR(1:lstep+1,:),lstep+1,numSrfs,
     &      'QHistCOR.dat',876)
         call OutputDataFile(PHistCOR(1:lstep+1,:),lstep+1,numSrfs,
     &      'PHistCOR.dat',877)
         call OutputDataFile(PlvHistCOR(1:lstep+1,:),lstep+1,numSrfs,
     &      'PlvHistCOR.dat',879)
      endif 

      return
      end
c 
c ... calculate initial conditions for the CalcSurfaces
c      
      subroutine calcCalcic(y,srfIdList,numSrfs)
      
      use calcFlowPressure
c
      include "common.h" !needed?
      include "mpif.h" !needed?
c      
      integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled
      real*8    y(nshg,4)   !need velocity and pressure
      real*8    Qini(0:MAXSURF) !initial flow rate
      real*8    PdistIni(0:MAXSURF)!initial distal pressure
      real*8    Pini(0:MAXSURF),CoupleArea(0:MAXSURF) ! initial pressure
      real*8    VelOnly(nshg,3), POnly(nshg)
c
      POnly(:)= one ! one to get area
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get surf area
      CalcArea(1:numSrfs) = CoupleArea(1:numSrfs)
      VelOnly(:,1:3)=y(:,1:3)
      call GetFlowQ(Qini,VelOnly,srfIdList,numSrfs) !get initial flow
      FlowHist(lstep+1,1:numSrfs)=Qini(1:numSrfs) !initialize QHistRCR
      POnly(:)=y(:,4) ! pressure
      call integrScalar(Pini,POnly,srfIdList,numSrfs) !get initial pressure integral
      Pini(1:numSrfs) = Pini(1:numSrfs)/CalcArea(1:numSrfs)
      PressHist(lstep+1,1:numSrfs)=Pini(1:numSrfs)
     
      return
      end
c 
c ... initialize the influence of the initial conditions for the RCR BC
c    
      subroutine calcRCRic(y,srfIdList,numSrfs)
      
      use convolRCRFlow    !brings RCRic, ValueListRCR, ValuePdist

      include "common.h" !needed?
      include "mpif.h" !needed?
      
      integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled
      real*8    y(nshg,4)   !need velocity and pressure
      real*8    Qini(0:MAXSURF) !initial flow rate
      real*8    PdistIni(0:MAXSURF)!initial distal pressure
      real*8    Pini(0:MAXSURF),CoupleArea(0:MAXSURF) ! initial pressure
      real*8    VelOnly(nshg,3), POnly(nshg)

      allocate (RCRic(0:MAXSURF))
      call RCRint(lstep,PdistIni) !get initial distal P 
      POnly(:)= one ! one to get area
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get surf area
      RCRArea(1:numSrfs) = CoupleArea(1:numSrfs)
      if (lstep .eq. zero) then
         VelOnly(:,1:3)=y(:,1:3)
         call GetFlowQ(Qini,VelOnly,srfIdList,numSrfs) !get initial flow
         QHistRCR(1,1:numSrfs)=Qini(1:numSrfs) !initialize QHistRCR
         POnly(:)=y(:,4) ! pressure
         call integrScalar(Pini,POnly,srfIdList,numSrfs) !get initial pressure integral
         Pini(1:numSrfs) = Pini(1:numSrfs)/RCRArea(1:numSrfs)
         PHistRCR(1,1:numSrfs)=Pini(1:numSrfs)
         RCRic(1:numSrfs) = Pini(1:numSrfs) 
     &          - ValueListRCR(1,:)*Qini(1:numSrfs)-PdistIni(1:numSrfs)
      elseif (lstep .gt. zero) then
          RCRic(1:numSrfs) = PHistRCR(1,1:numSrfs) 
     &     -ValueListRCR(1,1:numSrfs)*QHistRCR(1,1:numSrfs)
     &     -PdistIni(1:numSrfs)
      endif
      
      return
      end
c 
c ... initialize the influence of the initial conditions for the Coronary BC
c    
      subroutine calcCORic(y,srfIdList,numSrfs)
      
      use convolCORFlow    

      include "common.h" !needed?
      include "mpif.h" !needed?
      
      integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled
      real*8    y(nshg,4)   !need velocity and pressure
      real*8    CoupleArea(0:MAXSURF)
      real*8    VelOnly(nshg,3), POnly(nshg)
      real*8    Qini(0:MAXSURF), Pini(0:MAXSURF)
      real*8    PlvistIni(0:MAXSURF)
           
      call calcCOR()
      call calcCoefCOR()  
      POnly(:)= one ! one to get area
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get surf area
      CORArea(1:numSrfs)=CoupleArea(1:numSrfs)
      call CORint(Delt(1)*lstep,PlvistIni,zero)
      CORic = zero
      if (lstep .eq. zero) then
         VelOnly(:,1:3)=y(:,1:3)
         call GetFlowQ(Qini,VelOnly,srfIdList,numSrfs) !get initial flow
         QHistCOR(1,1:numSrfs)=Qini(1:numSrfs)
         POnly(:)=y(:,4) ! pressure
         call integrScalar(Pini,POnly,srfIdList,numSrfs) !get initial pressure integral
         PHistCOR(1,1:numSrfs) = Pini(1:numSrfs)/CORArea(1:numSrfs)
         Pini(1:numSrfs)=Pini(1:numSrfs)/CORArea(1:numSrfs)
         PlvHistCOR(1,1:numSrfs)=PlvistIni(1:numSrfs)
      elseif (lstep .gt. zero) then
         Qini(1:numSrfs) = QHistCOR(1,1:numSrfs)
         Pini(1:numSrfs) = PHistCOR(1,1:numSrfs)
         PlvistIni(1:numSrfs)=PlvHistCOR(1,1:numSrfs)
      endif
      
      do k=1, numSrfs
      CORic(1,k) = (ValueListCOR(6,k)*dPinidT(k)-
     &           ValueListCOR(6,k)*COR(2,k)*Pini(k)-
     &         ValueListCOR(3,k)*dQinidT(k)-(ValueListCOR(3,k)*COR(1,k)+
     &         ValueListCOR(2,k))*Qini(k)-
     &         ValueListCOR(8,k)*PlvistIni(k)*CORScaleFactor(k))
     &         /DetCOR(k)
      CORic(2,k) = (ValueListCOR(6,k)*dPinidT(k)-
     &           ValueListCOR(6,k)*COR(1,k)*Pini(k)-
     &         ValueListCOR(3,k)*dQinidT(k)-(ValueListCOR(3,k)*COR(2,k)+
     &         ValueListCOR(2,k))*Qini(k)-
     &         ValueListCOR(8,k)*PlvistIni(k)*CORScaleFactor(k))
     &         /DetCOR(k)
       enddo       
      
      return
      end
c
c...  calculates the coefficients needed for beta calculation in the Coronary BC
c

      subroutine calcCoefCOR()
      use convolCORFlow

      include "common.h"

      do k=1, numCORSrfs
      CoefCOR(1,k)= (ValueListCOR(3,k)*COR(1,k)*COR(1,k)
     &   +ValueListCOR(2,k)*COR(1,k)+ValueListCOR(1,k))/DetCOR(k)
      CoefCOR(2,k)=(ValueListCOR(3,k)*COR(2,k)*COR(2,k)
     &   +ValueListCOR(2,k)*COR(2,k)+ValueListCOR(1,k))/DetCOR(k)
      CoefCOR(3,k)=(ValueListCOR(8,k)*COR(1,k)+ValueListCOR(7,k))
     &             /DetCOR(k)*CORScaleFactor(k)
      CoefCOR(4,k)=(ValueListCOR(8,k)*COR(2,k)+ValueListCOR(7,k))
     &             /DetCOR(k)*CORScaleFactor(k)
      CoefCOR(5,k)=ValueListCOR(3,k)/ValueListCOR(6,k)
      enddo

      return
      end

c
c...  calculates dtCOR, the exponents if the exponentials for the Coronary BC
c
      subroutine calcCOR()

      use convolCORFlow ! brings ValueListCOR
 
      include "common.h"

      do k=1, numCORSrfs
         DetCOR(k)=sqrt(ValueListCOR(5,k)*ValueListCOR(5,k)
     &          -four*ValueListCOR(4,k)*ValueListCOR(6,k))
         COR(2,k)=-(ValueListCOR(5,k)+DetCOR(k))/two/ValueListCOR(6,k)
         COR(1,k)= ValueListCOR(4,k)/ValueListCOR(6,k)/COR(2,k) 
         dtCOR(1,k)=Delt(1)*COR(1,k)
         dtCOR(2,k)=Delt(1)*COR(2,k)   
      enddo

      return
      end

c
c.... calculate the time integral coefficients for the coupled inflow BC
c

      subroutine CalcINCPConvCoef (stepn, numSrfs)
      
      use incpBC
      
      include "common.h"
      
      integer  stepn, numSrfs, k, j
      real*8   Exponent
      
      INCPConvCoef = zero
c
c.... the following coefficients are used when the inductor is connected in series.
c      
       do k=1, numSrfs
         if (stepn .eq. 0) then
            INCPConvCoef(1,k) = alfi*Delt(1)/two
            INCPConvCoef(2,k) = alfi*Delt(1)/two 
         endif
         if (stepn .ge. 1) then
            INCPConvCoef(1,k) = Delt(1)/two
            INCPConvCoef(stepn+1,k) = (one+alfi)*Delt(1)/two 
            INCPConvCoef(stepn+2,k) = alfi*Delt(1)/two 
         endif
         if (stepn .ge. 2) then
            do j=2, stepn
               INCPConvCoef(j,k) = Delt(1)
            enddo
         endif
      enddo

c
c.... the following coefficients are used when the inductor is connected in parallel.
c      
c       do k=1, numSrfs        
c         if (ValueVv(7,k) .ne. zero) then
c            Exponent=Delt(1)*ValueVv(2,k)/ValueVv(7,k)
c         endif        
c         if (stepn .eq. 0) then
c            INCPConvCoef(1,k)=(one-(alfi*Exponent+one)
c     &         *exp(-alfi*Exponent))/alfi/Exponent
c            INCPConvCoef(2,k)=(alfi*Exponent-one+exp(-alfi*Exponent))
c     &         /alfi/Exponent 
c         endif
c         if (stepn .ge. 1) then
c            INCPConvCoef(1,k)=(exp(-(stepn-one+alfi)*Exponent)-(Exponent
c     &         +one)*exp(-(stepn+alfi)*Exponent))/Exponent
c            INCPConvCoef(stepn+1,k)=((Exponent-1)*exp(-alfi*Exponent)
c     &         +exp(-(one+alfi)*Exponent))/Exponent+(one-(alfi*Exponent
c     &         +one)*exp(-alfi*Exponent))/alfi/Exponent 
c            INCPConvCoef(stepn+2,k)=
c     &         (alfi*Exponent-one+exp(-alfi*Exponent))/alfi/Exponent  
c         endif
c         if (stepn .ge. 2) then
c            do j=2, stepn
c               INCPConvCoef(j,k)=((Exponent-one)*exp(-(stepn-j+one+alfi)
c     &            *Exponent)+exp(-(stepn-j+two+alfi)*Exponent))/Exponent
c     &            +(exp(-(stepn-j+alfi)*Exponent)-(Exponent+one)
c     &            *exp(-(stepn-j+one+alfi)*Exponent))/Exponent
c            enddo
c         endif
c      enddo
                
      return
      end
      
c
c.... calculate the time dependent INCPCoef for the coupled inflow BC
c
      subroutine CalcINCPCoef(timeINCP, stepn, srfIDList, numSrfs, y)   
      
      use incpBC
      
      include "common.h"
      
      integer  srfIDList(0:MAXSURF), stepn, numSrfs, k
      real*8   y(nshg, ndof), timeINCP
      real*8   ENext(0:MAXSURF), PvenousNext(0:MAXSURF)
      
      INCPCoef = zero
      Enext = zero
      PvenousNext = zero
      call INCPint(timeINCP*(stepn+alfi), ENext, PvenousNext)
      Eadjust(1:numSrfs) = ENext(1:numSrfs)
c
c.... the following coefficients are used when the inductor is connected in series.
c      
      do k=1, numSrfs
         INCPCoef(1,k)=ValueVv(2,k)
     &      +INCPConvCoef(stepn+2,k)*Eadjust(k)
     &      +ValueVv(7,k)/timeINCP/alfi
         INCPCoef(2,k)=Eadjust(k)*VLV(stepn+1,k)
     &      +Eadjust(k)*Qaorta(stepn+1,k)*INCPConvCoef(stepn+2,k)
     &      -ValueVv(7,k)/timeINCP/alfi*Qaorta(stepn+1,k)
      enddo
c
c.... the following coefficients are used when the inductor is connectied in parallel.
c
c      do k=1, numSrfs
c         INCPCoef(1,k)=ValueVv(2,k)+alfi*Delt(1)/two *Eadjust(k)
c         INCPCoef(2,k)=Eadjust(k)*VLV(stepn+1,k)
c     &      +Eadjust(k)*Qaorta(stepn+1,k)*alfi*Delt(1)/two 
c         if (ValueVv(7,k) .ne. zero) then
c            INCPCoef(1,k) = INCPCoef(1,k)
c     &         -ValueVv(2,k)*INCPConvCoef(stepn+2,k)
c            INCPCoef(2,k) = INCPCoef(2,k)-ValueVv(2,k)*poldINCP(k)
c     &         -(PLV(1,k)-Paorta(1,k)-ValueVv(2,k)*Qaorta(1,k))
c     &         *exp(-timeINCP*(stepn+alfi)*ValueVv(2,k)/ValueVv(7,k)) 
c         endif
c      enddo      
      
      return 
      end
      
      
      subroutine calcINCPic(timeINCP, y, srfIDList, numSrfs)
      
      use incpBC      
      
      include "common.h" 
      include "mpif.h"   
            
      integer  srfIDList(0:MAXSURF), numSrfs, k
      real*8   y(nshg, 4), timeINCP
      real*8   initPvenous(0:MAXSURF)
      real*8   initElast(0:MAXSURF), Integral(0:MAXSURF)
      real*8   CoupleArea(0:MAXSURF), Area(nshg)
      
      initElast = zero
      initPvenous = zero
      Integral = zero
      InflowArea = zero
      Area = one
      CoupleArea = zero
      INCPResidual = zero
      call integrScalar(CoupleArea, Area, srfIDList, numSrfs) !calculate inlet area
      InflowArea(1:numSrfs) = CoupleArea(1:numSrfs)
      call INCPint(zero, initElast, initPvenous)
      Eadjust(1:numSrfs) = initElast(1:numSrfs)
      if (lstep .eq. zero) then
         call integrScalar(Integral, y(:,4), srfIDList, numSrfs)
         Paorta(1,1:numSrfs)=
     &      Integral(1:numSrfs)/InflowArea(1:numSrfs)
         Integral = zero    
         call GetFlowQ(Integral, y(:,1:3), srfIDList, numSrfs)
         Qaorta(1,1:numSrfs)=Integral(1:numSrfs)
         QHistINCP(1,1:numSrfs)=Qaorta(1,1:numSrfs) 
         do k=1, numSrfs
            PLV(1,k)=Eadjust(k)*VLV(1,k)
            QAV(1,k)=zero
            if (PLV(1,k) .lt. initPvenous(k)) then
               QAV(1,k)=(initPvenous(k)-PLV(1,k))/ValueVv(1,k)
            endif              
         enddo
      elseif (lstep .gt. zero) then
         do k=1, numINCPSrfs
            do j=1, lstep+1
               if (Qaorta(j,k) .ne. zero) then
                  QHistINCP(j,k) = Qaorta(j,k)
               elseif (QAV(j,k) .ne. zero) then
                  QHistINCP(j,k) = QAV(j,k)
               else
                  QHistINCP(j,k) = zero
               endif
            enddo
         enddo
      endif
      
      return
      end

c
c.... calculate the time dependent INCPCoef for the coupled inflow BC
c
      subroutine UpdHeartModel(timeINCP, y, srfIDList,
     &   numSrfs, stepn) 

      use incpBC

      include "common.h"
      
      integer  srfIDList(0:MAXSURF), stepn, numSrfs, k
      real*8   y(nshg, ndof), timeINCP
      real*8   ENext(0:MAXSURF), Integral(0:MAXSURF)
      real*8   PvenousNext(0:MAXSURF)
      
      Enext = zero
      Integral = zero
      PvenousNext = zero
      call integrScalar(Integral, y(:,4), srfIDList, numSrfs)
      Paorta(stepn+1,1:numSrfs)=
     &   Integral(1:numSrfs)/InflowArea(1:numSrfs)
      Qaorta(stepn+1,1:numSrfs)=QHistINCP(stepn+1,1:numSrfs)
      call INCPint(timeINCP*stepn, ENext, PvenousNext)
      Eadjust(1:numSrfs) = ENext(1:numSrfs)

      do k=1, numSrfs
         if (Qaorta(stepn+1,k) .ne. zero) then
            VLV(stepn+1,k)=VLV(stepn,k)+timeINCP/2*Qaorta(stepn,k)
     &             +timeINCP/2*Qaorta(stepn+1,k)
            PLV(stepn+1,k)=Eadjust(k)*VLV(stepn+1,k)
            QAV(stepn+1,k)=zero
         elseif (Qaorta(stepn+1,k) .eq. zero) then
            VLV(stepn+1,k)=VLV(stepn,k)
     &         +timeINCP/2*(Qaorta(stepn,k)+QAV(stepn,k))
            PLV(stepn+1,k)=Eadjust(k)*VLV(stepn,k)
            QAV(stepn+1,k)=zero
            if (PLV(stepn+1,k) .lt. PvenousNext(k)) then
               QAV(stepn+1,k)=(PvenousNext(k)-PLV(stepn+1,k)
     &            +ValueVv(6,k)*QAV(stepn,k)/timeINCP)/
     &            (ValueVv(1,k)+timeINCP/2*Eadjust(k)
     &            +ValueVv(6,k)/timeINCP)
               VLV(stepn+1,k)=VLV(stepn+1,k)+timeINCP/2*QAV(stepn+1,k)
               PLV(stepn+1,k)=Eadjust(k)*VLV(stepn+1,k)
               QHistINCP(stepn+1,k)=QAV(stepn+1,k)
            elseif (QAV(stepn,k) .gt. zero) then
                QAV(stepn+1,k)=(PvenousNext(k)-PLV(stepn+1,k)
     &            +ValueVv(6,k)*QAV(stepn,k)/timeINCP)/
     &            (ValueVv(1,k)+timeINCP/2*Eadjust(k)
     &            +ValueVv(6,k)/timeINCP)
               VLV(stepn+1,k)=VLV(stepn+1,k)+timeINCP/2*QAV(stepn+1,k)
               PLV(stepn+1,k)=Eadjust(k)*VLV(stepn+1,k)
               QHistINCP(stepn+1,k)=QAV(stepn+1,k)              
            endif
         endif
      enddo
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &      (myrank .eq. zero)) then
         call OutputDataFile(Paorta(1:lstep+1,:),stepn+1,numSrfs,
     &      'Paorta.dat',821)
         call OutputDataFile(Qaorta(1:lstep+1,:),stepn+1,numSrfs,
     &      'Qaorta.dat',822)
         call OutputDataFile(PLV(1:lstep+1,:),stepn+1,numSrfs,
     &      'PLV.dat',823)
         call OutputDataFile(VLV(1:lstep+1,:),stepn+1,numSrfs,
     &      'VLV.dat',824)
         call OutputDataFile(QAV(1:lstep+1,:),stepn+1,numSrfs,
     &      'QAV.dat',825)
      endif
        
      return
      end

c
c.... This subroutine writes FlowHist.dat and PressHist.dat files
c
      subroutine Updcalc(y, srfIDList, numSrfs)

      use calcFlowPressure

      include "common.h"
c      
      real*8   y(nshg, ndof), NewP(0:MAXSURF), NewQ(0:MAXSURF)
      integer  srfIDList(0:MAXSURF),  numSrfs

      call GetFlowQ(NewQ, y(:,1:3), srfIDList,numSrfs) !new flow at time n+1
      FlowHist(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      call integrScalar(NewP, y(:,4), srfIDList, numSrfs)
      PressHist(lstep+1,1:numSrfs) = NewP(1:numSrfs)/CalcArea(1:numSrfs)
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &   (myrank .eq. zero)) then
         call OutputDataFile(FlowHist(1:lstep+1,:),lstep+1,numSrfs,
     &      'FlowHist.dat',1004)
         call OutputDataFile(PressHist(1:lstep+1,:),lstep+1,numSrfs,
     &      'PressHist.dat',1005)
      endif 

      return
      end

c
c.... This subroutine writes Lagrange Multipliers and errors in 
c.... LagrangeMultipliers.dat and LagrangeErrors.dat
c
      subroutine UpdateLagrangeCoef(y, col, row, srfIDList, numSrfs)

      use LagrangeMultipliers

      include "common.h"
      
      real*8   y(nshg, ndof)
	integer  col(nshg+1),	          row(nnz_tot)
      integer  srfIDList(0:MAXSURF),  numSrfs, NumOfData
      real*8   Integral(0:MAXSURF),   InnerProduct(0:MAXSURF,3)
 
      Integral = zero     
      InnerProduct = zero   
      call GetFlowQ(Integral, y(:,1:3), srfIDList, numSrfs)  
      QLagrange(1:numSrfs,1)=Integral(1:numSrfs)
      Integral = zero
      call GetProfileFlowQ(Integral, y(:,1:3), srfIDList, numSrfs) 
      PQLagrange(1:numSrfs,1)=Integral(1:numSrfs)
     &   /LagProfileArea(1:numSrfs) 
      Integral = zero
      LagSwitch = 0 
    	call CalcNANBLagrange(col, row, y(:,1:3))
      call GetInnerProduct(InnerProduct, y(:,1:3), srfIDList, numSrfs)
      IPLagrange(1:numSrfs,1:3)=InnerProduct(1:numSrfs,1:3)
      do k=1, numSrfs
         NumOfData = (k-1)*3+1
         LagErrorHist(lstep+1,NumOfData)=abs(IPLagrange(k,1)
     &      -two*QLagrange(k,1)*PQLagrange(k,1)
     &      +QLagrange(k,1)**2*ProfileDelta(k))
         LagErrorHist(lstep+1,NumOfData+1)=abs(IPLagrange(k,2))
         LagErrorHist(lstep+1,NumOfData+2)=abs(IPLagrange(k,3))
            LagErrorHist(lstep+1,NumOfData:NumOfData+2)=
     &      LagErrorHist(lstep+1,NumOfData:NumOfData+2)
     &      *LagMeanFlow(k)
         LagHist(lstep+1,NumOfData:NumOfData+2)=Lag(k,1:3)
      enddo    

      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and.
     &      (myrank .eq. zero)) then
         NumOfData = numLagrangeSrfs*3
         call OutputDataFile(LagHist(1:lstep+1,:),lstep+1, NumOfData,
     &      'LagrangeMultipliers.dat',801)
         call OutputDataFile(LagErrorHist(1:lstep+1,:),lstep+1,
     &      NumOfData,'LagrangeErrors.dat',802)
      endif
c
      return
      end  
c
c.... this function calculates an initial condition of a constrained surface
c
      subroutine calcLagrangeic(srfIDList, numSrfs)
c      
      use LagrangeMultipliers
      
      include "common.h" 

      integer  srfIDList(0:MAXSURF),  numSrfs 
      
      LagSwitch = 0 
      allocate(lhsLagL(9,nnz_tot,3))
      allocate(resL(numSrfs,3))
      allocate(LagAPproduct(nshg,3))
      lhsLagL = zero
      resL = zero   
      LagAPproduct = zero
      call MergeLagrangeParameters(srfIDList, numSrfs)
      ProfileDelta(1:numSrfs)=ProfileDelta(1:numSrfs)
     &   /LagProfileArea(1:numSrfs)/LagProfileArea(1:numSrfs)
      do k=1, numSrfs
         LagMeanFlow(k)=two*LagProfileArea(k)/LagMeanFlow(k)
     &      /LagMeanFlow(k)
         if (lstep .eq. zero) then
            LagHist(1,(k-1)*3+1:(k-1)*3+3)=Lagold(k,1:3)
         elseif (lstep .gt. zero) then
            Lagold(k,1:3)=LagHist(lstep+1,(k-1)*3+1:(k-1)*3+3)
         endif
      enddo
c            
      return 
      end
c
c.... this function calculates an area and plane vectors of a constrained surface
c
      subroutine MergeLagrangeParameters(srfIDList, numSrfs)
c      
      use LagrangeMultipliers
c      
      include "common.h" 
      include "mpif.h"   
c            
      integer  srfIDList(0:MAXSURF),  numSrfs 
      real*8   VectMag(3), Inplane1, Inplane2, Inplane3, InplaneNorm
      real*8, allocatable, dimension (:) :: TotalArea
      real*8, allocatable, dimension (:,:,:) :: InplaneVectors
c      
      allocate(TotalArea(numSrfs))
      allocate(InplaneVectors(3,3,numSrfs))
      TotalArea = zero
      InplaneVectors = zero
      call MPI_ALLREDUCE (LagProfileArea, TotalArea, numSrfs,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)  
      LagProfileArea(1:numSrfs)=TotalArea(1:numSrfs)
      TotalArea = zero
      call MPI_ALLREDUCE (ProfileDelta, TotalArea, numSrfs,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)  
      ProfileDelta(1:numSrfs)=TotalArea(1:numSrfs)
      InplaneVectors = zero
      
      do i=1,3
         do j=1,3
            call MPI_ALLREDUCE(LagInplaneVectors(i,j,:),
     &         InplaneVectors(i,j,:), numSrfs,  
     &         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr)
         enddo
      enddo
      LagInplaneVectors = InplaneVectors
      do k=1, numSrfs
         do i=1,3
            VectMag(i)=sqrt(LagInplaneVectors(1,i,k)**2+
     &         LagInplaneVectors(2,i,k)**2+LagInplaneVectors(3,i,k)**2)
         enddo
         if ( VectMag(1) .gt. zero .and. VectMag(2) .gt. zero) then
            LagInplaneVectors(1:3,1,k) = LagInplaneVectors(1:3,1,k)
     &         /VectMag(1)
            LagInplaneVectors(1:3,2,k)=LagInplaneVectors(1:3,2,k)
     &         /VectMag(2)
            Inplane1=-LagInplaneVectors(2,2,k)*LagInplaneVectors(3,1,k)
     &         +LagInplaneVectors(2,1,k)*LagInplaneVectors(3,2,k)
            Inplane2=-LagInplaneVectors(1,1,k)*LagInplaneVectors(3,2,k)
     &         +LagInplaneVectors(1,2,k)*LagInplaneVectors(3,1,k)
            Inplane3=-LagInplaneVectors(1,2,k)*LagInplaneVectors(2,1,k)
     &         +LagInplaneVectors(1,1,k)*LagInplaneVectors(2,2,k)
            InplaneNorm=one/sqrt(Inplane1**2+Inplane2**2+Inplane3**2)
            LagInplaneVectors(1,3,k)=Inplane1*InplaneNorm
            LagInplaneVectors(2,3,k)=Inplane2*InplaneNorm
            LagInplaneVectors(3,3,k)=Inplane3*InplaneNorm
         endif
      enddo
c            
      return 
      end      
c  
c.........function that integrates a scalar over a boundary
c
      subroutine integrScalar(scalInt,scal,srfIdList,numSrfs)

      use pvsQbi !brings ndsurf, NASC

      include "common.h"
      include "mpif.h"
      
      integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled, i, k
      real*8    scal(nshg), scalInt(0:MAXSURF), scalIntProc(0:MAXSURF)
      
      scalIntProc = zero
      do i = 1,nshg
        if(numSrfs.gt.zero) then
          do k = 1,numSrfs
            irankCoupled = 0
            if (srfIdList(k).eq.ndsurf(i)) then
              irankCoupled=k
              scalIntProc(irankCoupled) = scalIntProc(irankCoupled)
     &                            + NASC(i)*scal(i)
            endif      
          enddo       
        endif
      enddo
c      
c     at this point, each scalint has its "nodes" contributions to the scalar
c     accumulated into scalIntProc. Note, because NASC is on processor this
c     will NOT be the scalar for the surface yet
c
c.... reduce integrated scalar for each surface, push on scalInt
c
        npars=MAXSURF+1
       call MPI_ALLREDUCE (scalIntProc, scalInt(:), npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)  
   
      return
      end

c  
c.........function that outputs an input data array
c
      subroutine OutputDataFile(DataFile, nrows, ncolms, Filename,
     &   UnitNumber)

      include "common.h"
      
      character*40 Filename
      real*8    DataFile(nrows,ncolms)
      integer   nrows, ncolms, UnitNumber
      
      open(unit=UnitNumber, file=Filename,status='replace')
         write(UnitNumber,*) nrows
         do i=1, nrows
            write(UnitNumber,*) (DataFile(i,n),n=1, ncolms)
         enddo
      close(UnitNumber)
   
      return
      end





