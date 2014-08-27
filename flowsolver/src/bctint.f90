!-----------------------------------------------------------------------
!
!  This module conveys temporal BC data.  Below functions read in the data
!  and interpolate it to the current time level. 
!
!-----------------------------------------------------------------------
      module specialBC

      real*8, allocatable ::  BCt(:,:,:), acs(:,:), spamp(:)
      real*8, allocatable ::  ytarget(:,:)
      integer, allocatable :: nBCt(:), numBCt(:)
     

      integer ntv,nptsmax
!$$$      integer itvn
      end module

!-----------------------------------------------------------------------
!
!  This module conveys coupled inflow BC(INCP) data.  Below functions read
! in the data and interpolate it to the current time level. 
!
!-----------------------------------------------------------------------
      module incpBC

      real*8, allocatable  ::  ValueVv(:,:),         Tmax(:) 
      real*8, allocatable  ::  Period(:),            Pvenous(:,:,:) 
      real*8, allocatable  ::  Enormal(:,:),         Emax(:)
      real*8, allocatable  ::  INCPConvCoef(:,:),    INCPCoef(:,:)
      real*8, allocatable  ::  poldINCP(:),          QHistINCP(:,:)
      real*8, allocatable  ::  InflowArea(:)
      real*8, allocatable  ::  Paorta(:,:),          PLV(:,:)
      real*8, allocatable  ::  VLV(:,:),             QAV(:,:)
      real*8, allocatable  ::  Qaorta(:,:),          Eadjust(:)
      integer, allocatable ::  inactive(:)
      
      !(actually allocated in itrdrv.f but destroyed
      !with global destructor)
      ! *** chill Winston, this is now in itrdrv somewhere - kdl & nan 26/08/14 *** 
      ! real*8, allocatable :: iBCs(:,:)
      ! real*8, allocatable :: iBCd(:)
      
      integer nptsElast,   INCPSwitch
      integer nptsPvenous, nptsINCP
      
      end module
      
!-----------------------------------------------------------------------
!
!  This module conveys flow rate history for the different impedance outlets
!  over one period. Below functions read in the data and store it for the
!  current time level. 
!
!-----------------------------------------------------------------------
      module convolImpFlow

      real*8, allocatable ::  QHistImp(:,:), ValueImpt(:,:,:)
      real*8, allocatable ::  ValueListImp(:,:), ConvCoef(:,:)
      real*8, allocatable ::  ImpConvCoef(:,:), poldImp(:)
      integer ntimeptpT,numDataImp
      integer, allocatable :: nImpt(:), numImpt(:)
      integer nptsImpmax
      real*8, allocatable ::  QHistTry(:,:), QHistTryF(:,:) !filter
      integer cutfreq !filter
      end module

!-----------------------------------------------------------------------
!
!  This module conveys the parameters for the different RCR outlets.
!  Below functions read in the inputs (proximal resistance, capacitance, 
!  distal resistance and distal pressure) and store it for the
!  current time level. 
!
!-----------------------------------------------------------------------
      module convolRCRFlow

      real*8, allocatable ::  ValueListRCR(:,:), ValuePdist(:,:,:) !inputs
      real*8, allocatable ::  QHistRCR(:,:), PHistRCR(:,:), HopRCR(:) !calc
      real*8, allocatable ::  RCRConvCoef(:,:), poldRCR(:) !calc
      real*8, allocatable ::  dtRCR(:), RCRArea(:) !scaled timestep: deltat/RdC
      real*8, allocatable ::  RCRic(:) !(P(0)-RQ(0)-Pd(0))
      integer nptsRCRmax,numDataRCR, nptsRCR !to read inputs
      integer, allocatable :: numRCRt(:) !to read inputs
      end module
      
!-----------------------------------------------------------------------
!
!  This module conveys the parameters for time-varying RCR outlets.
!  Below functions read in the inputs (proximal resistance, capacitance,
!  distal resistance and distal pressure) and store it for the
!  current time level.
!
!-----------------------------------------------------------------------
      module convolTRCRFlow

      real*8, allocatable ::  ValueTRCR(:,:,:),  CurrValueTRCR(:,:,:)
      real*8, allocatable ::  QHistTRCR(:,:),    PHistTRCR(:,:)
      real*8, allocatable ::  TRCRConvCoef(:,:), TRCRArea(:)
      real*8, allocatable ::  dtTRCR(:,:),       HopTRCR(:)
      real*8, allocatable ::  TRCRic(:),         poldTRCR(:)
      real*8, allocatable ::  QmeanReg(:,:),     CtrlParamsReg(:,:)
      real*8, allocatable ::  CoeffReg(:,:),     RegIntCoeff(:,:)
      integer, allocatable :: numTRCRt(:)
      integer nptsTRCRmax,numDataTRCR, nptsTRCR
      integer RegStartStep, RegIniPeriod, RegFlowSwitch
      end module
      
!-----------------------------------------------------------------------
!
! This module conveys the parameters for the different Coronary outlets.
! Below functions read in the inputs (coronary resistances, capacitances,
! and left ventricular pressure) and store it for the current time level.
!
!-----------------------------------------------------------------------
      module convolCORFlow

      real*8, allocatable  :: ValuePlvist(:,:,:),    PlvHistCOR(:,:)
      real*8, allocatable  :: QHistCOR(:,:),         HopCOR(:) 
      real*8, allocatable  :: poldCOR(:),            plvoldCOR(:)
      real*8, allocatable  :: dQinidT(:),            dPinidT(:)
      real*8, allocatable  :: CORArea(:),            PHistCOR(:,:)
      complex, allocatable :: COR(:,:),              dtCOR(:,:) 
      complex, allocatable :: CORic(:,:),            ValueListCOR(:,:) 
      complex, allocatable :: CoefCOR(:,:),          DetCOR(:) 
      real*8, allocatable  :: CORConvCoef(:,:),      CORPlvConvCoef(:,:)
      real*8, allocatable  :: CORScaleFactor(:)
      integer, allocatable :: numCORt(:)
      integer nptsCORmax, numDataCOR, nptsCOR
      end module

!-----------------------------------------------------------------------
!
! This module conveys the parameters to save flow and pressure history.
!
!-----------------------------------------------------------------------
      module calcFlowPressure

      real*8, allocatable  :: FlowHist(:,:),         PressHist(:,:)
      real*8, allocatable  :: CalcArea(:)
      end module

!-----------------------------------------------------------------------
!
! This module conveys the parameters to save residuals.
!
!-----------------------------------------------------------------------
      module ResidualControl 

      real*8   controlResidual
      integer  CurrentIter
      end module
      
!-----------------------------------------------------------------------
!
!  This module conveys parameters for Lagrange multipliers. Below 
!  function reads in the inputs (LagCenter, LagRadius and ProfileOrder).
!  Defined variables are used to construct LHS and RHS of the solver. 
!
!-----------------------------------------------------------------------
      module LagrangeMultipliers

      real*8, allocatable  :: QLagrange(:,:),        PQLagrange(:,:)
      real*8, allocatable  :: NANBLagrange(:,:,:),   IPLagrange(:,:)   
      real*8, allocatable  :: Lag(:,:),              Lagold(:,:)
      real*8, allocatable  :: Lagincr(:,:),          Lagalpha(:,:)
      real*8, allocatable  :: LagCenter(:,:),        LagRadius(:)
      real*8, allocatable  :: ProfileDelta(:),       LagProfileArea(:)
      real*8, allocatable  :: loclhsLag(:,:,:,:,:),  lhsLagL(:,:,:)
      real*8, allocatable  :: LagErrorHist(:,:),     LagHist(:,:)
      real*8, allocatable  :: PenaltyCoeff(:,:),     Penalty(:,:)
      real*8, allocatable  :: ScaleFactor(:,:),      AddLag(:,:)
      real*8, allocatable  :: LagAPproduct(:,:),     resL(:,:)
      real*8, allocatable  :: LagInplaneVectors(:,:,:)
      real*8, allocatable  :: LagMeanFlow(:)
      integer, allocatable :: ProfileOrder(:) 
      integer LagSwitch  
      end module

!-----------------------------------------------------------------------
!
!     Initialize:
!
!-----------------------------------------------------------------------
      subroutine initSponge( y,x)
      
      use     specialBC
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8   y(nshg,nflow), x(numnp,3)
      allocate (ytarget(nshg,nflow))  
      
      if(matflg(5,1).eq.5) then
         write(*,*) 'calculating IC sponge'
         ytarget = y
      else
         write(*,*) 'calculating Analytic sponge'

!
! OLD style sponge pushed onto target.  You need to be sure that your
! solver.inp entries for start and stop of sponge match as well as the
! growth rates
!
      vcl=datmat(1,5,1)         ! velocity on centerline
      rslc=datmat(2,5,1)        ! shear layer center radius
      bfz=datmat(3,5,1)
      we=3.0*29./682.
      rsteep=3.0
      zstart=30.0
      radst=10.0
      radsts=radst*radst
      do id=1,numnp
         radsqr=x(id,2)**2+x(id,1)**2
!         if((x(id,3).gt. zstart) .or. (radsqr.gt.radsts))  then
            rad=sqrt(radsqr)
            radc=max(rad,radst)
            zval=max(x(id,3),zstart)
            utarget=(tanh(rsteep*(rslc-rad))+one)/two* &
                          (vcl-we) + we
            Ttarget  = press/(ro*Rgas)
            ptarget= press
            ytarget(id,1) = zero
            ytarget(id,2) = zero
            ytarget(id,3) = utarget
            ytarget(id,4) = ptarget
            ytarget(id,5) = Ttarget            
!         endif
      enddo
      endif
      return
      end


!-----------------------------------------------------------------------
!
!     Initialize:time varying boundary condition
!
!-----------------------------------------------------------------------
      subroutine initBCt( x, iBC, BC )
      
      use     specialBC
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8   x(numnp,nsd), BC(nshg,ndofBC), rj1,rj2,rj3,rj4,distd,epsd
      integer  iBC(numnp)
      character*80 card
      real*8 distds
      real*8 dd
!
!  This one should be used for boundary layer meshes where bct.dat must
!  be given to greater precision than is currently being generated.
!
!      epsd=1.0d-12    ! this is distance SQUARED to save square root

!      epsd=1.0d-8              ! this is distance SQUARED to save square root
      epsd=1.0d-6              ! new distance to avoid problems when using mm

      ic=0                      !count the number on this processor
     
      if(any(ibits(iBC,3,3).eq.7)) then
         write(*,*) 'opening bct.dat'
!         open(unit=567, file='bct.dat',status='old')
         open(unit=567, file='bct.dat',ACTION='READ',STATUS='old')
         read (567,'(a80)') card
           read (card,*) ntv, nptsmax
!        read(567,*) ntv,nptsmax
         allocate (nBCt(numnp))  
         allocate (numBCt(ntv))  
         allocate (BCt(ntv,nptsmax,4))  
         do k=1,ntv
            read(567,*) x1,x2,x3,ntpts
!
! Find the point on the boundary (if it is on this processor)
! that matches this point
!
            do i=1,numnp
               if(ibits(ibc(i),3,3) .eq.7) then
                  dd= distds(x1,x2,x3,x(i,1),x(i,2),x(i,3))
                  if(dd.lt.epsd) then
                     ic=ic+1
                     nBCt(ic)=i ! the pointer to this point
                     numBCt(ic)=ntpts ! the number of time series
                     do j=1,ntpts
!                        read(567,*) BCt(ic,j,4),(BCt(ic,j,n),n=1,3)
                        read(567,*) (BCt(ic,j,n),n=1,4)
                     enddo
                     exit
                  endif
               endif
            enddo
            if(i.eq.numnp+1) then
!
!  if we get here the point was not found.  It must be on another
!  processor so we read past this record and move on
!
               do j=1,ntpts
                  read(567,*) rj1,rj2,rj3,rj4
               enddo
            endif
         enddo                  ! end of the loop over ntv
      BCt(:,:,4)=BCt(:,:,4)*bcttimescale
      endif                     ! any 3 component nodes
      itvn=ic
      close(567)
      write(*,*)'myrank=',myrank,' and I found ',ic,' nodes.'


      return
      end


      subroutine BCint(timel,shp,shgl,shpb,shglb,x,BC,iBC)

      use     specialBC ! brings in itvn,nbct, bct, numbct, nptsmax

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      real*8   BC(nshg,ndofBC), timel,t
      real*8   x(numnp,nsd), &
               shp(MAXTOP,maxsh,MAXQPT), &
               shgl(MAXTOP,nsd,maxsh,MAXQPT), &
               shpb(MAXTOP,maxsh,MAXQPT), &
               shglb(MAXTOP,nsd,maxsh,MAXQPT)

      integer  iBC(numnp),nlast,i,j,nper 

      do i =1,itvn ! itvn is the number of varying nodes on this proc 

         nlast=numBCt(i)     ! number of time series to interpolate from
         nper=timel/BCt(i,nlast,4)! number of periods completed to shift off


         t=timel-nper*BCt(i,nlast,4)  ! now time in periodic domain

         do j=2,nlast   !loop to find the interval that we are in

            if(BCt(i,j,4).gt.t) then  ! this is upper bound, j-1 is lower

               wr=(t-BCt(i,j-1,4))/(BCt(i,j,4)-BCt(i,j-1,4))
               BC(nbct(i),3:5)= BCt(i,j-1,1:3)*(one-wr) &
                              + BCt(i,j,1:3)*wr
               exit

            endif
         enddo
      enddo
      return
      end

      function distds(x1,y1,z1,x2,y2,z2)
      real*8 distds 
      real*8 x1,y1,z1,x2,y2,z2,x,y,z
      x=x1-x2
      y=y1-y2
      z=z1-z2
      distds=x*x+y*y+z*z
      return
      end
      
      subroutine initCalcSrfst()
      
      use calcFlowPressure
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      allocate (FlowHist(lstep+nstep(1)+1,numCalcSrfs)) !for flow history
      allocate (PressHist(lstep+nstep(1)+1,numCalcSrfs)) !for pressure history
      FlowHist = zero
      PressHist = zero
      if (lstep .gt. 0) then
         call ReadDataFile(FlowHist(1:lstep+1,:),lstep+1,numCalcSrfs, &
            'FlowHist.dat',1004)
         call ReadDataFile(PressHist(1:lstep+1,:),lstep+1,numCalcSrfs, &
            'PressHist.dat',1005)
      endif
      
      end
      
!-----------------------------------------------------------------------
!   initialize the impedance boundary condition:
!   read the data in initImpt
!   interpolate the data to match the process time step in Impint
!-----------------------------------------------------------------------
      subroutine initImpt()
      
      use convolImpFlow
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      open(unit=817, file='impt.dat',status='old')
         read (817,*) nptsImpmax
         allocate (numImpt(numImpSrfs))  
         allocate (ValueImpt(nptsImpmax,2,numImpSrfs))
         ValueImpt=0
         do k=1,numImpSrfs
            read (817,*) numDataImp
            numImpt(k) = numDataImp
            do j=1,numDataImp
               read(817,*) (ValueImpt(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
         enddo
      close(817)
      
      allocate (ValueListImp(ntimeptpT+1,numImpSrfs))
      ValueListImp(ntimeptpT+1,:) = ValueImpt(1,2,:) !Z(time=0), last entry
      ValueListImp(1,:) = ValueImpt(1,2,:) !Z(time=0)=Z(time=T)
      return
      end
      
      
      
      subroutine Impint(ctime,jstep)
      
      use convolImpFlow
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8 ctime, ptime
      integer nlast, nper, k, j , jstep
      
         
      do k =1,numImpSrfs
         nlast=numImpt(k)     ! number of time series to interpolate from
         nper=ctime/ValueImpt(nlast,1,k)!number of periods completed to shift off
         ptime = ctime-nper*ValueImpt(nlast,1,k)  ! now time in periodic domain
            
         do j=2,nlast   !loop to find the interval that we are in

            if(ValueImpt(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
               wr=(ptime-ValueImpt(j-1,1,k)) &
                   / ( ValueImpt(j,1,k)-ValueImpt(j-1,1,k) ) 
               ValueListImp(jstep,k)= ValueImpt(j-1,2,k)*(one-wr)  &
                              + ValueImpt(j,2,k)*wr
               exit
            endif

         enddo
      enddo
      return
      end

!-----------------------------------------------------------------------------
!     time filter for a periodic function (sin cardinal + window function)     
!     is used for the impedance and the flow rate history
!-----------------------------------------------------------------------------
      subroutine Filter(Filtered,DataHist,nptf,timestep,cutfreq)
      
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      integer nptf, cutfreq, j, k, m, s, Filtime(nptf)
      real*8  DataHist(nptf,numImpSrfs), Window(nptf)
      real*8  Sinc(nptf), FilterSW(nptf), Filtered(nptf,numImpSrfs)
      real*8  windK, timestep

      windK = cutfreq*2 + 1
      do j=1,nptf
         Filtime(j) = j-1
         Window(j) = 0.42+0.5*cos(2*pi*Filtime(j)/nptf) &
                    +0.08*cos(4*pi*Filtime(j)/nptf)
         Sinc(j) = sin(pi*Filtime(j)*windK/nptf)/sin(pi*Filtime(j)/nptf) 
      enddo          
      Sinc(1) = windK
            
      do j=1,nptf     
         FilterSW(j) = Window(nptf+1-j)*Sinc(nptf+1-j) !filter for convolution
      enddo     
      
      Filtered = zero
      do m=1,nptf
         do j=1,nptf
            s=modulo(m-nptf+j,nptf)
            if(s.eq.zero) then
               s=nptf
            endif
            Filtered(m,:) = Filtered(m,:) &
                    +FilterSW(j)*DataHist(s,:)/nptf !filter convolution
         enddo
      enddo
      
      return
      end
!-----------------------------------------------------------------------
!   initialize the RCR boundary condition:
!   read the data in initRCRt
!   interpolate the data to match the process time step in RCRint
!-----------------------------------------------------------------------
      subroutine initRCRt()
      
      use convolRCRFlow
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      open(unit=818, file='rcrt.dat',status='old')
         read (818,*) nptsRCRmax
         allocate (numRCRt(numRCRSrfs))  
         allocate (RCRArea(numRCRSrfs))  
         allocate (ValuePdist(nptsRCRmax,2,numRCRSrfs))
         allocate (ValueListRCR(3,numRCRSrfs))
         RCRArea = zero
         ValuePdist=0
         ValueListRCR=0
         do k=1,numRCRSrfs
            read (818,*) numDataRCR
            numRCRt(k) = numDataRCR
            do j=1,3
               read(818,*) ValueListRCR(j,k) ! reads Rp,C,Rd
            enddo
            do j=1,numDataRCR
               read(818,*) (ValuePdist(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
         enddo
      close(818)

      allocate (dtRCR(numRCRSrfs))
      
      nptsRCR = lstep            
      allocate (QHistRCR(lstep+nstep(1)+1,numRCRSrfs))
      allocate (RCRConvCoef(lstep+nstep(1)+2,numRCRSrfs)) !for convolution coeff
      allocate (PHistRCR(lstep+nstep(1)+1,numRCRSrfs))
      PHistRCR = zero
      QHistRCR = zero
      if (lstep .gt. 0) then
         call ReadDataFile(QHistRCR(1:lstep+1,:),lstep+1,numRCRSrfs, &
            'QHistRCR.dat',870)
         call ReadDataFile(PHistRCR(1:lstep+1,:),lstep+1,numRCRSrfs, &
            'PHistRCR.dat',871)
      endif
      
      return
      end
           
      
      subroutine RCRint(ctime,Pdist)
      
      use convolRCRFlow ! brings numRCRSrfs, ValuePdist
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8  ctime, ptime
      integer nlast, nper, k, j
      real*8  Pdist(0:MAXSURF)      
         
      do k =1,numRCRSrfs
         nlast=numRCRt(k)     ! number of time series to interpolate from
         nper=ctime/ValuePdist(nlast,1,k)!number of periods completed to shift off
         ptime = ctime-nper*ValuePdist(nlast,1,k)  ! now time in periodic domain
            
         do j=2,nlast   !loop to find the interval that we are in

            if(ValuePdist(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
               wr=(ptime-ValuePdist(j-1,1,k)) &
                   / ( ValuePdist(j,1,k)-ValuePdist(j-1,1,k) )
               Pdist(k)= ValuePdist(j-1,2,k)*(one-wr) &
                              + ValuePdist(j,2,k)*wr
               exit
            endif

         enddo
      enddo
      return
      end
      
!-----------------------------------------------------------------------
!   initialize the time-varying RCR boundary condition:
!   read the data in initTRCRt
!   interpolate the data to match the process time step in TRCRint
!-----------------------------------------------------------------------
      subroutine initTRCRt()

      use convolTRCRFlow
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      open(unit=881, file='timercrt.dat',status='old')
         read (881,*) nptsTRCRmax
         allocate (numTRCRt(numTRCRSrfs))
         allocate (TRCRArea(numTRCRSrfs))
         allocate (ValueTRCR(nptsTRCRmax,5,numTRCRSrfs))
         numTRCRt = zero
         TRCRArea = zero
         ValueTRCR = zero
         do k=1,numTRCRSrfs
            read (881,*) numDataTRCR
            numTRCRt(k) = numDataTRCR
            do j=1,numDataTRCR
               read(881,*) (ValueTRCR(j,n,k), n=1,5) ! n=1 time, n=2 Rp, n=3 C, n=4 Rd, n=5 Pd
            enddo
         enddo
      close(881)

      allocate (poldTRCR(0:MAXSURF))
      allocate (HopTRCR(0:MAXSURF))
      allocate (dtTRCR(lstep+nstep(1)+2,numTRCRSrfs))
      allocate (QHistTRCR(lstep+nstep(1)+1,numTRCRSrfs))
      allocate (PHistTRCR(lstep+nstep(1)+1,numTRCRSrfs))
      allocate (TRCRic(0:MAXSURF))
      allocate (TRCRConvCoef(lstep+nstep(1)+2,numTRCRSrfs)) !for convolution coeff
      allocate (CurrValueTRCR(lstep+nstep(1)+2,numTRCRSrfs,4))
      poldTRCR = zero
      HopTRCR = zero
      dtTRCR = zero
      QHistTRCR = zero
      PHistTRCR = zero
      TRCRic = zero
      TRCRConvCoef = zero
      CurrValueTRCR = zero
      if (lstep .eq. 0) then
         nptsTRCR = 0
      elseif (lstep .gt. 0) then
         nptsTRCR = lstep
         call ReadDataFile(QHistTRCR(1:lstep+1,:),lstep+1,numTRCRSrfs, &
            'QHistTRCR.dat',882)
         call ReadDataFile(PHistTRCR(1:lstep+1,:),lstep+1,numTRCRSrfs, &
            'PHistTRCR.dat',883)
         call ReadDataFile(CurrValueTRCR(1:lstep+1,:,1),lstep+1, &
            numTRCRSrfs,'ValueTRCRRp.dat',984)
         call ReadDataFile(CurrValueTRCR(1:lstep+1,:,2),lstep+1, &
            numTRCRSrfs,'ValueTRCRC.dat',985)
         call ReadDataFile(CurrValueTRCR(1:lstep+1,:,3),lstep+1, &
            numTRCRSrfs,'ValueTRCRRd.dat',986)
         call ReadDataFile(CurrValueTRCR(1:lstep+1,:,4),lstep+1, &
            numTRCRSrfs,'ValueTRCRPd.dat',987)
         dtTRCR(1:lstep+1,:)=Delt(1) &
           /CurrValueTRCR(1:lstep+1,:,2)/CurrValueTRCR(1:lstep+1,:,3)
      endif

!      if (regflow.gt.0) then
!         allocate(QmeanReg(lstep+nstep(1)+1,numRegSrfs))
!         allocate(CtrlParamsReg(lstep+nstep(1)+1,numRegSrfs))
!         allocate(CoeffReg(5,numRegSrfs))
!         allocate(RegIntCoeff(lstep+nstep(1)+1,numRegSrfs))
!         QmeanReg = zero
!         CtrlParamsReg = zero
!         CoeffReg = zero
!         RegFlowSwitch = zero
!         open(unit=981, file='regflow.dat',status='old')
!            read(981,*)
!            read(981,*) RegStartStep
!            read(981,*)
!            read(981,*) RegIniPeriod
!            do k=1,numRegSrfs
!               read(981,*) (CoeffReg(n,k), n=1,5) ! n=1 desired flow, n=2 scaling factor, n=3 minimum resistance
!            enddo
!         close(981)
!         if (lstep .gt. 0) then
!            call ReadDataFile(QmeanReg(1:lstep+1,:),lstep+1,numRegSrfs,
!     &         'qmeanreg.dat',989)
!            call ReadDataFile(CtrlParamsReg(1:lstep+1,:),lstep+1,
!     &         numRegSrfs,'ctrlparamsreg.dat',988)
!         endif
!      endif

      return
      end


      subroutine TRCRint(curstep)

      use convolTRCRFlow
!      use FeedbackSystem
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      real*8  ctime, TotResist
      integer nlast, k, j, i, curstep, n

      ctime = curstep*Delt(1)
      do k =1,numTRCRSrfs
         nlast=numTRCRt(k)
         do j=2,nlast
            if(ValueTRCR(j,1,k) .gt. ctime) then  ! this is upper bound, j-1 is lower
               wr=(ctime-ValueTRCR(j-1,1,k)) &
                  / ( ValueTRCR(j,1,k)-ValueTRCR(j-1,1,k) )
               do i=1,4
                  CurrValueTRCR(curstep+1,k,i)=ValueTRCR(j-1,i+1,k) &
                             *(one-wr)+ ValueTRCR(j,i+1,k)*wr
               enddo
               dtTRCR(curstep+1,k)=Delt(1) &
                 /CurrValueTRCR(curstep+1,k,2) &
                 /CurrValueTRCR(curstep+1,k,3)
!               if (feedback .gt. zero .and. curstep .gt. WaitStep) then
!                  do n=1, numEfferent
!                     if (nsrflistEfferent(n) .eq. nsrflistTRCR(k)) then
!                        CurrValueTRCR(curstep+1,k,3)=
!     &                     CurrValueTRCR(curstep+1,k,3)
!     &                     *ControlParameterValues(curstep,3)
!                     endif
!                  enddo
!               endif

!               if (regflow .gt. 0 .and. curstep .gt. RegStartStep) then
!                  do n=1, numRegSrfs
!                     if (nsrflistReg(n) .eq. nsrflistTRCR(k)) then
!                        CurrValueTRCR(curstep+1,k,3)=
!     &                     CurrValueTRCR(curstep+1,k,3)
!     &                     *CtrlParamsReg(curstep,n)
!                        TotResist=CurrValueTRCR(curstep+1,k,1)
!     &                     +CurrValueTRCR(curstep+1,k,3)
!                        if (TotResist .lt. CoeffReg(3,n)) then
!                           CurrValueTRCR(curstep+1,k,3)=CoeffReg(3,n)
!     &                        -CurrValueTRCR(curstep+1,k,1)
!                        endif
!                     endif
!                  enddo
!               endif

               exit
            endif
         enddo
      enddo

      return
      end

!-----------------------------------------------------------------------
!   initialize the Coronary boundary condition:
!   read the data in initCORt
!   interpolate the data to match the process time step in CORint
!-----------------------------------------------------------------------
      subroutine initCORt()
      
      use convolCORFlow
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      open(unit=815, file='cort.dat',status='old')
         read (815,*)
         read (815,*) nptsCORmax
         allocate (numCORt(numCorSrfs))  
         allocate (ValuePlvist(nptsCORmax,2,numCORSrfs))
         allocate (ValueListCOR(9,numCORSrfs))
         allocate (dQinidT(numCORSrfs))
         allocate (dPinidT(numCORSrfs))
         allocate (CORScaleFactor(numCORSrfs))
         numCORt = zero
         ValuePlvist = zero
         ValueListCOR = zero
         dQinidT = zero
         dPinidT = zero
         do k=1,numCORSrfs
            read (815,*)
            read (815,*) numDataCOR
            numCORt(k) = numDataCOR
            read (815,*)
            do j=1,9
               read (815,*)
               read (815,*) ValueListCOR(j,k) ! reads q0, q1, q2, p0, p1, p2, b0, b1, b2
            enddo
            read (815,*)
            read (815,*) dQinidT(k)
            read (815,*)
            read (815,*) dPinidT(k)
            read (815,*)
            read (815,*) CORScaleFactor(k)
            read (815,*)
            do j=1,numDataCOR
               read(815,*) (ValuePlvist(j,n,k),n=1,2) ! n=1 time, 2 value
            enddo
         enddo
      close(815)
       
      allocate (dtCOR(2,numCORSrfs))
      allocate (COR(2, numCORSrfs))
      allocate (CoefCOR(5,numCORSrfs))
      allocate (DetCOR(numCORSrfs))

      if (lstep .eq. 0) then
         nptsCOR = 0
         allocate (CORConvCoef(nstep(1)+2,numCORSrfs)) !for convolution coeff
         allocate (CORPlvConvCoef(nstep(1)+2,numCORSrfs))
         allocate (QHistCOR(nstep(1)+1,numCORSrfs)) !for flow history
         allocate (PlvHistCOR(nstep(1)+1,numCORSrfs))
         allocate (PHistCOR(nstep(1)+1,numCORSrfs)) !for pressure history
         QHistCOR = zero
         PlvHistCOR = zero
         PHistCOR = zero
      elseif (lstep .gt. 0) then   
         nptsCOR = lstep            
         allocate (CORConvCoef(lstep+nstep(1)+2,numCORSrfs)) !for convolution coeff
         allocate (CORPlvConvCoef(lstep+nstep(1)+2,numCORSrfs))
         allocate (QHistCOR(lstep+nstep(1)+1,numCORSrfs)) !for flow history
         allocate (PlvHistCOR(lstep+nstep(1)+1,numCORSrfs))
         allocate (PHistCOR(lstep+nstep(1)+1,numCORSrfs)) !for pressure history
         QHistCOR = zero
         PlvHistCOR = zero
         PHistCOR = zero
         call ReadDataFile(QHistCOR(1:lstep+1,:),lstep+1,numCORSrfs, &
            'QHistCOR.dat',876)
         call ReadDataFile(PHistCOR(1:lstep+1,:),lstep+1,numCORSrfs, &
            'PHistCOR.dat',877)
         call ReadDataFile(PlvHistCOR(1:lstep+1,:),lstep+1,numCORSrfs, &
            'PlvHistCOR.dat',879)
      endif

      return
      end
           
      
      subroutine CORint(ctime,Plvist,Switch)
      
      use convolCORFlow ! brings numCORSrfs, ValuePlvist
      use incpBC
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8  ctime, ptime
      integer nlast, nper, k, j, curstep, Switch
      real*8  Plvist(0:MAXSURF)      
         
      do k =1,numCORSrfs
         nlast=numCORt(k)     ! number of time series to interpolate from
         nper=ctime/ValuePlvist(nlast,1,k)!number of periods completed to shift off
         ptime = ctime-nper*ValuePlvist(nlast,1,k)  ! now time in periodic domain
         if (incp .eq. zero) then
            do j=2,nlast   !loop to find the interval that we are in
               if(ValuePlvist(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
                  wr=(ptime-ValuePlvist(j-1,1,k)) &
                      / ( ValuePlvist(j,1,k)-ValuePlvist(j-1,1,k) )
                  Plvist(k)= ValuePlvist(j-1,2,k)*(one-wr) &
                              + ValuePlvist(j,2,k)*wr
                  exit
               endif
            enddo
         else 
            curstep = ctime/Delt(1)
            if (Switch .eq. zero) then
               Plvist(k)=PLV(curstep+1,1) 
            elseif(Switch .eq. one .and. INCPSwitch .lt. two) then
               Plvist(k)=PLV(curstep,1)
            else
              Plvist(k)=PLV(curstep,1)*(one-alfi)+ &
                 alfi*PLV(curstep+1,1)
            endif
         endif
      enddo
      
      return
      end
!-----------------------------------------------------------------------
!   initialize the coupled inflow boundary condition:
!   read the data in initINCPt
!   interpolate the data to match the process time step in INCPint
!-----------------------------------------------------------------------
      subroutine initINCPt()
      
      use incpBC
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      allocate(ValueVv(7,numINCPSrfs))
      allocate(Period(numINCPSrfs))
      allocate(Tmax(numINCPSrfs))
      allocate(Emax(numINCPSrfs))
      ValueVv = zero
      Period = zero
      Tmax = zero
      Emax = zero

      open(unit=819, file='incp.dat', status='old')
      read(819,*) nptsPvenous
      allocate(Pvenous(nptsPvenous,2,numINCPSrfs))
      do k=1, numINCPSrfs
         read(819,*) 
         read(819,*) ValueVv(1,k)   !RA-V
         read(819,*)
         read(819,*) ValueVv(2,k)   !RV-art
         read(819,*)
         read(819,*) ValueVv(3,k)   !Vv(0), initial ventricular volume
         read(819,*)
         read(819,*) ValueVv(4,k)   !Vo, correction volume
         read(819,*)
         read(819,*) ValueVv(5,k)   !shifted time
         read(819,*)
         read(819,*) ValueVv(6,k)   !LA-V
         read(819,*)
         read(819,*) ValueVv(7,k)   !LV-art
         read(819,*)
         read(819,*) Period(k)      !period
         read(819,*)
         read(819,*) Tmax(k)        !time from the onset of systole to the end systole
         read(819,*)
         read(819,*) Emax(k)        !maximum elastance
         read(819,*)
         do j=1, nptsPvenous
            read(819,*) (Pvenous(j,n,k), n=1,2)
         enddo
      enddo
      close(819)
      
      open(unit=820, file='Elastance.dat', status='old')
         read (820,*) nptsElast
         allocate (Enormal(nptsElast,2))
         do k=1, nptsElast
            read(820,*) (Enormal(k,n), n=1,2)  ! n=1, time, n=2, value
         enddo
      close(820)
      
      if (lstep .eq. 0) then
         nptsINCP = 0
         allocate (Paorta(nstep(1)+1,numINCPSrfs))
         allocate (Qaorta(nstep(1)+1,numINCPSrfs))
         allocate (PLV(nstep(1)+1,numINCPSrfs))
         allocate (VLV(nstep(1)+1,numINCPSrfs))
         allocate (QAV(nstep(1)+1,numINCPSrfs))
         Paorta = zero
         Qaorta = zero
         PLV = zero
         VLV = zero
         QAV = zero        
      elseif (lstep .gt. 0) then       
         nptsINCP = lstep        
         allocate (Paorta(nstep(1)+nptsINCP+1,numINCPSrfs))
         allocate (Qaorta(nstep(1)+nptsINCP+1,numINCPSrfs))
         allocate (PLV(nstep(1)+nptsINCP+1,numINCPSrfs))
         allocate (VLV(nstep(1)+nptsINCP+1,numINCPSrfs))
         allocate (QAV(nstep(1)+nptsINCP+1,numINCPSrfs))
         call ReadDataFile(Paorta(1:lstep+1,:),lstep+1,numINCPSrfs, &
            'Paorta.dat',821)
         call ReadDataFile(Qaorta(1:lstep+1,:),lstep+1,numINCPSrfs, &
            'Qaorta.dat',822)
         call ReadDataFile(PLV(1:lstep+1,:),lstep+1,numINCPSrfs, &
            'PLV.dat',823)
         call ReadDataFile(VLV(1:lstep+1,:),lstep+1,numINCPSrfs, &
            'VLV.dat',824)
         call ReadDataFile(QAV(1:lstep+1,:),lstep+1,numINCPSrfs, &
            'QAV.dat',825)
      endif
      VLV(1,1:numINCPSrfs) = ValueVv(3,1:numINCPSrfs)- &
              ValueVv(4,1:numINCPSrfs)
      
      return
      end
      
!
!.... This function scales the normalized elastance function according
!.... to the input parameters and interpolates the data to match the 
!.... process time step

      subroutine INCPint(ctime, Elastance, curPvenous)
      
      use incpBC   !Need Period, Tmax, Emax, Enormal
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8   ctime, ptime, shifttime, Elast
      integer  nlast, nper, j, k, l
      real*8   Elastance(0:MAXSURF), curPvenous(0:MAXSURF)
      real*8   Etime(nptsElast), Evalue(nptsElast), wr
      
      Elastance = zero
      nlast=nptsElast   !number of time series to interpolate from
      do k=1, numINCPSrfs
         shifttime=Tmax(k)*ValueVv(5,k)
         nper=(ctime+shifttime)/Period(k)  !number of periods completed to shift off
         ptime = ctime+shifttime-nper*Period(k)   !now time in periodic domain, 0 to Period(k)
         Etime = Enormal(:,1)*Tmax(k) 
         Evalue = Enormal(:,2)*Emax(k)*1333.2237
!         if (Etime(nlast) .gt. Period(k)) then
!            do l=2, nlast
!               if (Etime(l) .gt. Period(k)) then
!                  wr=(Period(k)-Etime(l-1))/(Etime(l)-Etime(l-1))
!                  Elast = Evalue(l-1)*(one-wr)+Evalue(l)*wr
!                  exit
!               endif
!            enddo
!         else 
!            Elast = Evalue(nlast)
!         endif 
         Elast = Evalue(nlast)

!
!.... Here I assume the simulation starts at the onset of systole
!
         do j=2, nlast
            if (Etime(j) .gt. ptime) then 
               wr=(ptime-Etime(j-1))/(Etime(j)-Etime(j-1))
               Elastance(k) = Evalue(j-1)*(one-wr)+Evalue(j)*wr
               if (nper .gt. 0 .and. j .eq. 2 .or. j .eq. nlast) then
                  wr=(ptime-Etime(j-1))/(Etime(j)-Etime(j-1))
                  Elastance(k) = Elast*(one-wr)+Evalue(j)*wr 
               endif 
               exit
            elseif (ptime .gt. Etime(nlast)) then
               Elastance(k) = Elast    
               exit
            endif
         enddo
         nlast=nptsPvenous     ! number of time series to interpolate from 
         do j=2,nlast   !loop to find the interval that we are in
            if(Pvenous(j,1,k).gt.ptime) then  ! this is upper bound, j-1 is lower
               wr=(ptime-Pvenous(j-1,1,k)) &
                   / ( Pvenous(j,1,k)-Pvenous(j-1,1,k) )
               curPvenous(k)= Pvenous(j-1,2,k)*(one-wr) &
                              + Pvenous(j,2,k)*wr
               exit
            endif
         enddo
      enddo
      
      return
      end
!-----------------------------------------------------------------------
!   Read data for Lagrange multipliers: read input data in initLagrange
!   This data is required to generate profile functions
!-----------------------------------------------------------------------
      subroutine initLagrange()
      
      use LagrangeMultipliers 
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      integer NumOfData
      
      allocate(LagCenter(3,numLagrangeSrfs))
      allocate(LagInplaneVectors(3,3,numLagrangeSrfs))
      allocate(LagRadius(numLagrangeSrfs))
      allocate(LagProfileArea(numLagrangeSrfs))
      allocate(Lagold(numLagrangeSrfs,3))
      allocate(Lag(numLagrangeSrfs,3))
      allocate(Lagincr(numLagrangeSrfs,3))
      allocate(Lagalpha(numLagrangeSrfs,3))
      allocate(ProfileOrder(numLagrangeSrfs))
      allocate(ProfileDelta(numLagrangeSrfs))
      allocate(Penalty(numLagrangeSrfs,3))
      allocate(PenaltyCoeff(numLagrangeSrfs,3))
      allocate(ScaleFactor(numLagrangeSrfs,3))
      allocate(AddLag(numLagrangeSrfs,3))
      allocate(LagMeanFlow(numLagrangeSrfs))
      NumOfData = numLagrangeSrfs*3
      allocate(LagHist(lstep+nstep(1)+1, NumOfData))
      allocate(LagErrorHist(lstep+nstep(1)+1, NumOfData))
      LagCenter = zero
      LagInplaneVectors = zero
      LagRadius = zero
      LagProfileArea = zero
      Lagold = zero
      Lag = zero
      Lagincr = zero
      Lagalpha = zero
      ProfileOrder = 0
      LagSwitch = 0
      ProfileDelta = zero
      Penalty = zero
      PenaltyCoeff = zero
      ScaleFactor = zero
      AddLag = zero
      LagMeanFlow = zero
      LagHist = zero
      LagErrorHist = zero
      open(unit=800, file='LagrangeData.dat', status='old')
      do k=1, numLagrangeSrfs
         read(800,*)
         read(800,*) (LagCenter(n,k), n=1,3)  !Center of a constrained surface
         read(800,*)
         read(800,*) LagRadius(k)        !Surface radius
         read(800,*)
         read(800,*) ProfileOrder(k)      !Profile order
         read(800,*)
         read(800,*) LagMeanFlow(k)      !Mean flow
         read(800,*)
         read(800,*) (Lagold(k,n), n=1,3) !Initial Lagrange Multipliers 
         read(800,*)
         read(800,*) (PenaltyCoeff(k,n), n=1,3) !Penalty numbers
         read(800,*)
         read(800,*) (ScaleFactor(k,n), n=1,3) !Scaling factors
      enddo
      close(800)
      if (lstep .gt. zero) then
         call ReadDataFile(LagHist(1:lstep+1,:),lstep+1,NumOfData, &
           'LagrangeMultipliers.dat',801)
         call ReadDataFile(LagErrorHist(1:lstep+1,:),lstep+1,NumOfData, &
           'LagrangeErrors.dat',802)
      endif
      
      return
      end         
!----------------------------------------------------------------------- 
!     returns in pold the history dependent part of the pressure in the
!     impedance/flow rate convolution for the impedance, RCR, Coronary 
!     and INCP BC      
!-----------------------------------------------------------------------      
      subroutine pHist(pressHist,QHist,betas,nTimePoint,nSrfs)

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer  nTimePoint,nSrfs
      real*8   pressHist(0:MAXSURF)
      real*8   QHist(nTimePoint+1,nSrfs),betas(nTimePoint+2,nSrfs)
      !don't need here betas(ntimePoint+2)
      !but pb of array passing if cut at nTimePoint+1
      pressHist=zero
      do k=1,nSrfs
        do j=1,nTimePoint+1
            pressHist(k) = pressHist(k) + QHist(j,k)*betas(j,k)
        enddo
      enddo
      return
      end


!----------------------------------------------------------------------- 
! This subroutine reads a data file and copies to a data array
!----------------------------------------------------------------------- 
      subroutine ReadDataFile(DataFile,nrows,ncolms,Filename,UnitNumber)

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      character(*) Filename
      real*8    DataFile(nrows,ncolms)
      integer   nrows, ncolms, UnitNumber
      
      open(unit=UnitNumber, file=Filename, status='old')
         read(UnitNumber,*) 
         do i=1, nrows
            read(UnitNumber,*) (DataFile(i,j), j=1, ncolms)
         enddo
      close(UnitNumber)
   
      return
      end
