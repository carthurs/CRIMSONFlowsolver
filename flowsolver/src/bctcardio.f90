!
!...initialize the coefficients for the impedance convolution
!
      subroutine CalcImpConvCoef (numISrfs, numTpoints)

      use convolImpFlow !uses flow history and impedance for convolution
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !for alfi
      
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
!
!...calculate the coefficients for the impedance convolution
! 
      allocate (ImpConvCoef(numTpoints+2,numISrfs))

!..coefficients below assume Q linear in time step, Z constant
!            do j=3,numTpoints
!                ImpConvCoef(j,:) = ValueListImp(j-1,:)*ConvCoef(j,3)
!     &                             + ValueListImp(j,:)*ConvCoef(j,2)    
!     &                             + ValueListImp(j+1,:)*ConvCoef(j,1)  
!            enddo
!            ImpConvCoef(1,:) = ValueListImp(2,:)*ConvCoef(1,1)
!            ImpConvCoef(2,:) = ValueListImp(2,:)*ConvCoef(2,2)    
!     &                       + ValueListImp(3,:)*ConvCoef(2,1)
!            ImpConvCoef(numTpoints+1,:) =
!     &           ValueListImp(numTpoints,:)*ConvCoef(numTpoints+1,3)
!     &         + ValueListImp(numTpoints+1,:)*ConvCoef(numTpoints+1,2) 
!            ImpConvCoef(numTpoints+2,:) = 
!     &           ValueListImp(numTpoints+1,:)*ConvCoef(numTpoints+2,3)

!..try easiest convolution Q and Z constant per time step
      do j=3,numTpoints+1
         ImpConvCoef(j,:) = ValueListImp(j-1,:)/numTpoints
      enddo
      ImpConvCoef(1,:) =zero
      ImpConvCoef(2,:) =zero
      ImpConvCoef(numTpoints+2,:) =  &
                 ValueListImp(numTpoints+1,:)/numTpoints
! compensate for yalpha passed not y in Elmgmr()
      ImpConvCoef(numTpoints+1,:)= ImpConvCoef(numTpoints+1,:) &
                        - ImpConvCoef(numTpoints+2,:)*(1.0-alfi)/alfi 
      ImpConvCoef(numTpoints+2,:)= ImpConvCoef(numTpoints+2,:)/alfi 
      return
      end

! 
! ... update the flow rate history for the impedance convolution, filter it and write it out
!    
      subroutine UpdHistConv(y,nsrfIdList,numSrfs)
      
      use convolImpFlow !brings ntimeptpT, QHistImp, QHistTry, QHistTryF, numImpSrfs
      use convolRCRFlow !brings QHistRCR, numRCRSrfs
      use convolTRCRFlow
      use convolCORFlow 
      use incpBC
      use boundarymodule,only: GetFlowQ
!
      use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !needed?
      include "mpif.h" !needed?
      
      integer   nsrfIdList(0:MAXSURF), numSrfs
      real*8    y(nshg,3) !velocity at time n+1   
      real*8    NewQ(0:MAXSURF)

      call GetFlowQ(NewQ,y,nsrfIdList,numSrfs) !new flow at time n+1
!
!... for imp BC: shift QHist, add new constribution, filter and write out
!      
      if(numImpSrfs.gt.zero .and. nsrfIdList(1).eq.nsrflistImp(1)) then
         do j=1, ntimeptpT
            QHistImp(j,1:numSrfs)=QHistImp(j+1,1:numSrfs)
         enddo
         QHistImp(ntimeptpT+1,1:numSrfs) = NewQ(1:numSrfs)
         QHistImp(1,:)=zero

!
!....filter the flow rate history
!
!         cutfreq = 10           !hardcoded cutting frequency of the filter
!         do j=1, ntimeptpT
!            QHistTry(j,:)=QHistImp(j+1,:)
!         enddo
!         call Filter(QHistTryF,QHistTry,ntimeptpT,Delt(1),cutfreq)
!         QHistImp(1,:)=zero
!         do j=1, ntimeptpT
!            QHistImp(j+1,:)=QHistTryF(j,:)
!         enddo
!
!.... write out the new history of flow rates to Qhistor.dat
!      
         if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
                     (myrank .eq. zero)) then
            open(unit=816, file='Qhistor.dat',status='replace')
            write(816,*) ntimeptpT
            do j=1,ntimeptpT+1
               write(816,*) (QHistImp(j,n),n=1, numSrfs)
            enddo
            close(816)
         endif
      endif 

!
!... for RCR bc just add the new contribution
!
      if(numRCRSrfs.gt.zero .and. nsrfIdList(1).eq.nsrflistRCR(1)) then
         QHistRCR(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      endif
!
!... for time-varying RCR bc just add the new contribution
!
      if(numTRCRSrfs.gt.zero.and.nsrfIdList(1).eq.nsrflistTRCR(1)) then
         QHistTRCR(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      endif
!
!... for Coronary bc just add the new contribution
!
      if(numCORSrfs.gt.zero.and.nsrfIdList(1).eq.nsrflistCOR(1)) then
         QHistCOR(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      endif      
!
!... for INCP bc just add the new contribution
!
      if(numINCPSrfs.gt.zero.and.nsrfIdList(1).eq.nsrflistINCP(1)) then
         QHistINCP(lstep+1,1:numSrfs) = NewQ(1:numSrfs)  
      endif 
      
      return
      end

!
!...calculate the time varying coefficients for the RCR convolution
!
      subroutine CalcRCRConvCoef (stepn, numSrfs)

      use convolRCRFlow !brings in ValueListRCR, dtRCR
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !brings alfi
      
      integer numSrfs, stepn    

      RCRConvCoef = zero
      if (stepn .eq. 0) then
        RCRConvCoef(1,:) = ValueListRCR(1,:)*(1.0-alfi) + &
         ValueListRCR(3,:)*(-alfi + 1.0 + 1/dtRCR(:)  &
           - exp(-alfi*dtRCR(:))*(1 + 1/dtRCR(:)))
        RCRConvCoef(2,:) = ValueListRCR(1,:)*alfi  &
           + ValueListRCR(3,:) &
           *(alfi - 1/dtRCR(:) + exp(-alfi*dtRCR(:))/dtRCR(:))
      endif
      if (stepn .ge. 1) then
        RCRConvCoef(1,:) =-ValueListRCR(3,:)*exp(-dtRCR(:)*(stepn+alfi)) &
              *(1 + (1 - exp(dtRCR(:)))/dtRCR(:))
        RCRConvCoef(stepn+1,:) = ValueListRCR(1,:)*(1-alfi)  &
           - ValueListRCR(3,:)*(alfi - 1 - 1/dtRCR(:)  &
           + exp(-alfi*dtRCR(:))/dtRCR(:)*(2 - exp(-dtRCR(:))))
        RCRConvCoef(stepn+2,:) = ValueListRCR(1,:)*alfi  &
           + ValueListRCR(3,:) &
           *(alfi - 1/dtRCR(:) + exp(-alfi*dtRCR(:))/dtRCR(:))
      endif
      if (stepn .ge. 2) then
        do j=2,stepn
         RCRConvCoef(j,:) = ValueListRCR(3,:)/dtRCR(:)* &
              exp(-dtRCR(:)*(stepn + alfi + 2 - j))* &
              (1 - exp(dtRCR(:)))**2
        enddo
      endif

! compensate for yalpha passed not y in Elmgmr()
      RCRConvCoef(stepn+1,:)= RCRConvCoef(stepn+1,:) &
                        - RCRConvCoef(stepn+2,:)*(1.0-alfi)/alfi 
      RCRConvCoef(stepn+2,:)= RCRConvCoef(stepn+2,:)/alfi 

      return
      end

!
!...calculate the time dependent H operator for the RCR convolution
!
      subroutine CalcHopRCR (timestepRCR, stepn, numSrfs)

      use convolRCRFlow !brings in HopRCR, dtRCR

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h" !needed?

      integer numSrfs, stepn      
      real*8  PdistCur(0:MAXSURF), timestepRCR
      
      HopRCR=zero
      call RCRint(timestepRCR*(stepn + alfi),PdistCur)
      HopRCR(1:numSrfs) = RCRic(1:numSrfs)  &
           *exp(-dtRCR(1:numSrfs)*(stepn + alfi)) + PdistCur(1:numSrfs)
      return
      end

!
!...calculate the time varying coefficients for time-varying RCR convolution
!
      subroutine CalcTRCRConvCoef (stepn, numSrfs)

      use convolTRCRFlow !brings in ValueListRCR, dtRCR

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      integer numSrfs, stepn

!      if (regflow.gt.0) then
!         call UpdRegParams(stepn)
!      endif

      call TRCRint(stepn+1)
      TRCRConvCoef = zero
      if (stepn .eq. 0) then
        TRCRConvCoef(1,:) = CurrValueTRCR(1,:,1)*(one-alfi) + &
          CurrValueTRCR(1,:,3)*(-alfi + one + 1/dtTRCR(1,:) &
          - exp(-alfi*dtTRCR(1,:))*(1 + 1/dtTRCR(1,:)))
        TRCRConvCoef(2,:) = CurrValueTRCR(2,:,1)*alfi &
          + CurrValueTRCR(2,:,3) &
          *(alfi - 1/dtTRCR(2,:) + exp(-alfi*dtTRCR(2,:))/dtTRCR(2,:))
      endif
      if (stepn .ge. 1) then
        TRCRConvCoef(1,:)=-CurrValueTRCR(1,:,3) &
          *exp(-dtTRCR(1,:)*(stepn+alfi)) &
          *(1 + (1 - exp(dtTRCR(1,:)))/dtTRCR(1,:))
        TRCRConvCoef(stepn+1,:) = CurrValueTRCR(stepn+1,:,1)*(1-alfi) &
          - CurrValueTRCR(stepn+1,:,3)*(alfi - 1 - 1/dtTRCR(stepn+1,:) &
          + exp(-alfi*dtTRCR(stepn+1,:))/dtTRCR(stepn+1,:) &
          * (2 - exp(-dtTRCR(stepn+1,:))))
        TRCRConvCoef(stepn+2,:) = CurrValueTRCR(stepn+2,:,1)*alfi &
          + CurrValueTRCR(stepn+2,:,3) &
          *(alfi - 1/dtTRCR(stepn+2,:) + exp(-alfi*dtTRCR(stepn+2,:)) &
          /dtTRCR(stepn+2,:))
      endif
      if (stepn .ge. 2) then
        do j=2,stepn
         TRCRConvCoef(j,:) = CurrValueTRCR(j,:,3)/dtTRCR(j,:)* &
             exp(-dtTRCR(j,:)*(stepn + alfi + 2 - j))* &
             (1 - exp(dtTRCR(j,:)))**2
        enddo
      endif

! compensate for yalpha passed not y in Elmgmr()
      TRCRConvCoef(stepn+1,:)= TRCRConvCoef(stepn+1,:) &
                       - TRCRConvCoef(stepn+2,:)*(1.0-alfi)/alfi
      TRCRConvCoef(stepn+2,:)= TRCRConvCoef(stepn+2,:)/alfi

      return
      end

!
!...calculate the time dependent H operator for time-varying RCR convolution
!
      subroutine CalcHopTRCR (timestepTRCR, stepn, numSrfs)

      use convolTRCRFlow !brings in HopRCR, dtRCR

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      include "mpif.h" !needed?

      integer numSrfs, stepn
      real*8  PdistCur(0:MAXSURF), timestepTRCR

      HopTRCR=zero
      HopTRCR(1:numSrfs) = TRCRic(1:numSrfs) &
        *exp(-dtTRCR(stepn+2,1:numSrfs)*(stepn+alfi)) &
        +CurrValueTRCR(stepn+2,:,4)

      return
      end

!
!.... This subroutine writes FlowHist.dat and PressHist.dat files
!
      subroutine UpdRCR(y, srfIDList, numSrfs)

      use convolRCRFlow 
      use boundarymodule, only: integrScalar
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!      
      real*8   y(nshg, ndof), NewP(0:MAXSURF)
      integer  srfIDList(0:MAXSURF),  numSrfs
      
      call integrScalar(NewP,y(:,4),srfIdList,numSrfs)
         PHistRCR(lstep+1,1:numSrfs)=NewP(1:numSrfs)/RCRArea(1:numSrfs)
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
         (myrank .eq. zero)) then
         call OutputDataFile(QHistRCR(1:lstep+1,:),lstep+1,numSrfs, &
            'QHistRCR.dat',870)
         call OutputDataFile(PHistRCR(1:lstep+1,:),lstep+1,numSrfs, &
            'PHistRCR.dat',871)
      endif 

      return
      end
      
!
!.... This subroutine writes FlowHist.dat and PressHist.dat files
!
      subroutine UpdTRCR(y, srfIDList, numSrfs)

      use convolTRCRFlow
      use boundarymodule, only: integrScalar
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      real*8   y(nshg, ndof), NewP(0:MAXSURF)
      integer  srfIDList(0:MAXSURF),  numSrfs

      call integrScalar(NewP,y(:,4),srfIDList,numSrfs)
      PHistTRCR(lstep+1,1:numSrfs)=NewP(1:numSrfs)/TRCRArea(1:numSrfs)
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
        (myrank .eq. zero)) then
         call OutputDataFile(QHistTRCR(1:lstep+1,:),lstep+1,numSrfs, &
           'QHistTRCR.dat',882)
         call OutputDataFile(PHistTRCR(1:lstep+1,:),lstep+1,numSrfs, &
           'PHistTRCR.dat',883)
         call OutputDataFile(CurrValueTRCR(1:lstep+1,:,1),lstep+1, &
           numSrfs,'ValueTRCRRp.dat',984)
         call OutputDataFile(CurrValueTRCR(1:lstep+1,:,2),lstep+1, &
           numSrfs,'ValueTRCRC.dat',985)
         call OutputDataFile(CurrValueTRCR(1:lstep+1,:,3),lstep+1, &
           numSrfs,'ValueTRCRRd.dat',986)
         call OutputDataFile(CurrValueTRCR(1:lstep+1,:,4),lstep+1, &
           numSrfs,'ValueTRCRPd.dat',987)
      endif

      return
      end

!
!...calculate the time varying coefficients of pressure 
!...for the Coronary convolution
!
      subroutine CalcCORConvCoef (stepn, numSrfs)

      use convolCORFlow !brings in dtCOR, COR, CoefCOR

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"

      integer numSrfs, stepn    

      CORConvCoef = zero
      if (stepn .eq. 0) then
        CORConvCoef(1,:)=real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)/alfi* &
                    ((-one+alfi*dtCOR(1,:))*exp(dtCOR(1,:)*alfi)+one) &
                        - CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                     ((-one+alfi*dtCOR(2,:))*exp(dtCOR(2,:)*alfi)+one))
        CORConvCoef(2,:)= real(CoefCOR(5,:) + &
                          CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)/alfi* &
                        (exp(dtCOR(1,:)*alfi)-(one+alfi*dtCOR(1,:)))- &
                          CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                        (exp(dtCOR(2,:)*alfi)-(one+alfi*dtCOR(2,:))))
      endif
      if (stepn .ge. 1) then  
         CORConvCoef(1,:) =real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)* &
                        (exp(dtCOR(1,:)*(stepn+alfi-one))- &
                      (one-dtCOR(1,:))*exp(dtCOR(1,:)*(stepn+alfi))) &
                       -CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)* &
                       (exp(dtCOR(2,:)*(stepn+alfi-one))- &
                      (one-dtCOR(2,:))*exp(dtCOR(2,:)*(stepn+alfi))))
         CORConvCoef(stepn+1,:) = real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:) &
                              /alfi*(alfi*exp(dtCOR(1,:)*(alfi+one))- &
                              (alfi+one)*exp(dtCOR(1,:)*(alfi))+one)- &
                              CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                              (alfi*exp(dtCOR(2,:)*(alfi+1))- &
                              (alfi+one)*exp(dtCOR(2,:)*(alfi))+one))
         CORConvCoef(stepn+2,:) = real(CoefCOR(5,:)+ &
                              CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)/alfi* &
                           (exp(dtCOR(1,:)*alfi)-one-alfi*dtCOR(1,:))- &
                              CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                            (exp(dtCOR(2,:)*alfi)-one-alfi*dtCOR(2,:)))
      endif
      if (stepn .ge. 2) then
      do j=2,stepn
         CORConvCoef(j,:) = real(CoefCOR(1,:)/dtCOR(1,:)/COR(1,:)* &
                            (exp(dtCOR(1,:)*(stepn+alfi-j))- &
                            two*exp(dtCOR(1,:)*(stepn+alfi-j+one))+ &
                           exp(dtCOR(1,:)*(stepn+alfi-j+two)))- &
                            CoefCOR(2,:)/dtCOR(2,:)/COR(2,:)* &
                            (exp(dtCOR(2,:)*(stepn+alfi-j))- &
                            two*exp(dtCOR(2,:)*(stepn+alfi-j+one))+ &
                            exp(dtCOR(2,:)*(stepn+alfi-j+two))))
      enddo
      endif
      
      return
      end

!
!...calculate the time varying coefficients of left ventricular pressure
!...for the Coronary convolution
!...need to do for t=0 and 1
!

      subroutine CalcCORPlvConvCoef (stepn, numSrfs)

      use convolCORFlow !brings in dtCOR, COR, CoefCOR

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      integer numSrfs, stepn    
      
      CORPlvConvCoef = zero
      if (stepn .eq. 0) then
         CORPlvConvCoef(1,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)/alfi* &
                   ((-one+alfi*dtCOR(1,:))*exp(dtCOR(1,:)*alfi)+one) &
                       - CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                   ((-one+alfi*dtCOR(2,:))*exp(dtCOR(2,:)*alfi)+one))
         CORPlvConvCoef(2,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)/alfi* &
                      (exp(dtCOR(1,:)*alfi)-(one+alfi*dtCOR(1,:)))- &
                         CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                       (exp(dtCOR(2,:)*alfi)-(one+alfi*dtCOR(2,:))))
      endif
      if (stepn .ge. 1) then
       CORPlvConvCoef(1,:) = real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)* &
                        (exp(dtCOR(1,:)*(stepn+alfi-one))- &
                   (one-dtCOR(1,:))*exp(dtCOR(1,:)*(stepn+alfi))) &
                      -CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)* &
                       (exp(dtCOR(2,:)*(stepn+alfi-one))- &
                     (one-dtCOR(2,:))*exp(dtCOR(2,:)*(stepn+alfi))))
      CORPlvConvCoef(stepn+1,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:) &
                             /alfi*(alfi*exp(dtCOR(1,:)*(alfi+one))- &
                             (alfi+one)*exp(dtCOR(1,:)*(alfi))+one)- &
                              CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                              (alfi*exp(dtCOR(2,:)*(alfi+1))- &
                              (alfi+one)*exp(dtCOR(2,:)*(alfi))+one))
       CORPlvConvCoef(stepn+2,:)=real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:) &
                      /alfi*(exp(dtCOR(1,:)*alfi)-one-alfi*dtCOR(1,:))- &
                              CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)/alfi* &
                           (exp(dtCOR(2,:)*alfi)-one-alfi*dtCOR(2,:)))
      endif
      if (stepn .ge. 2) then
      do j=2,stepn
         CORPlvConvCoef(j,:) = real(CoefCOR(3,:)/dtCOR(1,:)/COR(1,:)* &
                            (exp(dtCOR(1,:)*(stepn+alfi-j))- &
                            two*exp(dtCOR(1,:)*(stepn+alfi-j+one))+ &
                           exp(dtCOR(1,:)*(stepn+alfi-j+two)))- &
                            CoefCOR(4,:)/dtCOR(2,:)/COR(2,:)* &
                            (exp(dtCOR(2,:)*(stepn+alfi-j))- &
                            two*exp(dtCOR(2,:)*(stepn+alfi-j+one))+ &
                            exp(dtCOR(2,:)*(stepn+alfi-j+two))))
      enddo
      endif

      return
      end

!
!...calculate the time dependent H operator for the Coronary convolution
!
      subroutine CalcHopCOR (timestepCOR, stepn, srfIdList, numSrfs, y)

      use convolCORFlow !brings in HopCOR, dtCOR, COR, CoefCOR
      use boundarymodule, only: integrScalar      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !needed?
      include "mpif.h" !needed?

      integer   srfIdList(0:MAXSURF), numSrfs, stepn
      real*8    y(nshg,4), timestepCOR, PlvistNext(0:MAXSURF)
      real*8    CoupleArea(0:MAXSURF), POnly(nshg)
      
      HopCOR=zero
      PlvistNext=zero
      
      call CalcCORPlvConvCoef (stepn, numSrfs)
      call pHist(plvoldCOR, PlvHistCOR, CORPlvConvCoef, &
          nstep+nptsCOR,numSrfs) 
      call CORint(timestepCOR*(stepn + alfi),PlvistNext,1)
      POnly(:)=y(:,4) ! pressure
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get initial pressure integral
      PHistCOR(stepn+1,1:numSrfs) = CoupleArea(1:numSrfs) &
         /CORArea(1:numSrfs)
      HopCOR(1:numSrfs) =real(plvoldCOR(1:numSrfs)+  &
            CORic(1,1:numSrfs)*exp(dtCOR(1,1:numSrfs)*(stepn+alfi))- &
            CORic(2,1:numSrfs)*exp(dtCOR(2,1:numSrfs)*(stepn+alfi))+  &
            CorPlvConvCoef(stepn+2,1:numSrfs)*PlvistNext(1:numSrfs))

      return
      end
!
!...update time history of left ventricular pressure
!
      subroutine UpdHistPlvConv(y,timestepCOR,stepn,srfIdList,numSrfs) 
      
      use convolCORFlow 
      use boundarymodule, only: integrScalar      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision 
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
      PHistCOR(stepn+1,1:numSrfs) = CoupleArea(1:numSrfs) &
         /CORArea(1:numSrfs)

      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
         (myrank .eq. zero)) then
         call OutputDataFile(QHistCOR(1:lstep+1,:),lstep+1,numSrfs, &
            'QHistCOR.dat',876)
         call OutputDataFile(PHistCOR(1:lstep+1,:),lstep+1,numSrfs, &
            'PHistCOR.dat',877)
         call OutputDataFile(PlvHistCOR(1:lstep+1,:),lstep+1,numSrfs, &
            'PlvHistCOR.dat',879)
      endif 

      return
      end
! 
! ... calculate initial conditions for the CalcSurfaces
!      
      subroutine calcCalcic(y,srfIdList,numSrfs)
      
      use calcFlowPressure
      use boundarymodule, only: area, integrScalar, GetFlowQ      
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !needed?
      include "mpif.h" !needed?
!      
      integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled
      real*8    y(nshg,4)   !need velocity and pressure
      real*8    Qini(0:MAXSURF) !initial flow rate
      real*8    PdistIni(0:MAXSURF)!initial distal pressure
      real*8    Pini(0:MAXSURF),CoupleArea(0:MAXSURF) ! initial pressure
      real*8    VelOnly(nshg,3), POnly(nshg)
!

      POnly(:)= one ! one to get area
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs) !get surf area
      CalcArea(1:numSrfs) = CoupleArea(1:numSrfs)

      do i = 1, numSrfs
        write(*,*) i,' ',srfIdList(i),' ',CalcArea(i)       
        CalcArea(i) = area(srfIdList(i))
        write(*,*) i,' ',srfIdList(i),' ',CalcArea(i) 
      end do



      VelOnly(:,1:3)=y(:,1:3)
      call GetFlowQ(Qini,VelOnly,srfIdList,numSrfs) !get initial flow
      FlowHist(lstep+1,1:numSrfs)=Qini(1:numSrfs) !initialize QHistRCR
      POnly(:)=y(:,4) ! pressure
      call integrScalar(Pini,POnly,srfIdList,numSrfs) !get initial pressure integral
      Pini(1:numSrfs) = Pini(1:numSrfs)/CalcArea(1:numSrfs)
      PressHist(lstep+1,1:numSrfs)=Pini(1:numSrfs)
     
      return
      end
! 
! ... initialize the influence of the initial conditions for the RCR BC
!    
      subroutine calcRCRic(y,srfIdList,numSrfs)
      
      use convolRCRFlow    !brings RCRic, ValueListRCR, ValuePdist
      use boundarymodule, only: integrScalar
      use boundarymodule, only: GetFlowQ      
      use phcommonvars

      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !needed?
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
         RCRic(1:numSrfs) = Pini(1:numSrfs)  &
                - ValueListRCR(1,:)*Qini(1:numSrfs)-PdistIni(1:numSrfs)
      elseif (lstep .gt. zero) then
          RCRic(1:numSrfs) = PHistRCR(1,1:numSrfs)  &
           -ValueListRCR(1,1:numSrfs)*QHistRCR(1,1:numSrfs) &
           -PdistIni(1:numSrfs)
      endif
      
      return
      end

!
! ... compute area and initial flow and pressure for the time-varying RCR BC
!
      subroutine calcTRCRic(y,srfIdList,numSrfs)

      use convolTRCRFlow
      use boundarymodule, only: integrScalar, GetFlowQ      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !needed?

      integer   srfIdList(0:MAXSURF), numSrfs
      real*8    y(nshg,4), VelOnly(nshg,3), POnly(nshg)
      real*8    Qini(0:MAXSURF), Pini(0:MAXSURF)
      real*8    CoupleArea(0:MAXSURF)

      POnly(:) = one
      call integrScalar(CoupleArea,POnly,srfIdList,numSrfs)
      call TRCRint(zero)
      TRCRArea(1:numSrfs) = CoupleArea(1:numSrfs)
      if (lstep .eq. zero) then
         VelOnly(:,1:3)=y(:,1:3)
         call GetFlowQ(Qini,VelOnly,srfIdList,numSrfs)
         QHistTRCR(1,1:numSrfs)=Qini(1:numSrfs)
         POnly(:)=y(:,4)
         call integrScalar(Pini,POnly,srfIdList,numSrfs)
         Pini(1:numSrfs) = Pini(1:numSrfs)/TRCRArea(1:numSrfs)
         PHistTRCR(1,1:numSrfs)=Pini(1:numSrfs)
      endif
      TRCRic(1:numSrfs)=PHistTRCR(1,1:numSrfs) &
        -CurrValueTRCR(1,1:numSrfs,1)*QHistTRCR(1,1:numSrfs) &
        -CurrValueTRCR(1,1:numSrfs,4)

      return
      end

! 
! ... initialize the influence of the initial conditions for the Coronary BC
!    
      subroutine calcCORic(y,srfIdList,numSrfs)
      
      use convolCORFlow    
      use boundarymodule, only: integrScalar, GetFlowQ      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision !needed?
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
      CORic(1,k) = (ValueListCOR(6,k)*dPinidT(k)- &
                 ValueListCOR(6,k)*COR(2,k)*Pini(k)- &
               ValueListCOR(3,k)*dQinidT(k)-(ValueListCOR(3,k)*COR(1,k)+ &
               ValueListCOR(2,k))*Qini(k)- &
               ValueListCOR(8,k)*PlvistIni(k)*CORScaleFactor(k)) &
               /DetCOR(k)
      CORic(2,k) = (ValueListCOR(6,k)*dPinidT(k)- &
                 ValueListCOR(6,k)*COR(1,k)*Pini(k)- &
               ValueListCOR(3,k)*dQinidT(k)-(ValueListCOR(3,k)*COR(2,k)+ &
               ValueListCOR(2,k))*Qini(k)- &
               ValueListCOR(8,k)*PlvistIni(k)*CORScaleFactor(k)) &
               /DetCOR(k)
       enddo       
      
      return
      end
!
!...  calculates the coefficients needed for beta calculation in the Coronary BC
!

      subroutine calcCoefCOR()
      use convolCORFlow

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      do k=1, numCORSrfs
      CoefCOR(1,k)= (ValueListCOR(3,k)*COR(1,k)*COR(1,k) &
         +ValueListCOR(2,k)*COR(1,k)+ValueListCOR(1,k))/DetCOR(k)
      CoefCOR(2,k)=(ValueListCOR(3,k)*COR(2,k)*COR(2,k) &
         +ValueListCOR(2,k)*COR(2,k)+ValueListCOR(1,k))/DetCOR(k)
      CoefCOR(3,k)=(ValueListCOR(8,k)*COR(1,k)+ValueListCOR(7,k)) &
                   /DetCOR(k)*CORScaleFactor(k)
      CoefCOR(4,k)=(ValueListCOR(8,k)*COR(2,k)+ValueListCOR(7,k)) &
                   /DetCOR(k)*CORScaleFactor(k)
      CoefCOR(5,k)=ValueListCOR(3,k)/ValueListCOR(6,k)
      enddo

      return
      end

!
!...  calculates dtCOR, the exponents if the exponentials for the Coronary BC
!
      subroutine calcCOR()

      use convolCORFlow ! brings ValueListCOR
 
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      do k=1, numCORSrfs
         DetCOR(k)=sqrt(ValueListCOR(5,k)*ValueListCOR(5,k) &
                -four*ValueListCOR(4,k)*ValueListCOR(6,k))
         COR(2,k)=-(ValueListCOR(5,k)+DetCOR(k))/two/ValueListCOR(6,k)
         COR(1,k)= ValueListCOR(4,k)/ValueListCOR(6,k)/COR(2,k) 
         dtCOR(1,k)=Delt(1)*COR(1,k)
         dtCOR(2,k)=Delt(1)*COR(2,k)   
      enddo

      return
      end

!
!.... calculate the time integral coefficients for the coupled inflow BC
!

      subroutine CalcINCPConvCoef (stepn, numSrfs)
      
      use incpBC
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer  stepn, numSrfs, k, j
      real*8   Exponent
      
      INCPConvCoef = zero
!
!.... the following coefficients are used when the inductor is connected in series.
!      
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

!
!.... the following coefficients are used when the inductor is connected in parallel.
!      
!       do k=1, numSrfs        
!         if (ValueVv(7,k) .ne. zero) then
!            Exponent=Delt(1)*ValueVv(2,k)/ValueVv(7,k)
!         endif        
!         if (stepn .eq. 0) then
!            INCPConvCoef(1,k)=(one-(alfi*Exponent+one)
!     &         *exp(-alfi*Exponent))/alfi/Exponent
!            INCPConvCoef(2,k)=(alfi*Exponent-one+exp(-alfi*Exponent))
!     &         /alfi/Exponent 
!         endif
!         if (stepn .ge. 1) then
!            INCPConvCoef(1,k)=(exp(-(stepn-one+alfi)*Exponent)-(Exponent
!     &         +one)*exp(-(stepn+alfi)*Exponent))/Exponent
!            INCPConvCoef(stepn+1,k)=((Exponent-1)*exp(-alfi*Exponent)
!     &         +exp(-(one+alfi)*Exponent))/Exponent+(one-(alfi*Exponent
!     &         +one)*exp(-alfi*Exponent))/alfi/Exponent 
!            INCPConvCoef(stepn+2,k)=
!     &         (alfi*Exponent-one+exp(-alfi*Exponent))/alfi/Exponent  
!         endif
!         if (stepn .ge. 2) then
!            do j=2, stepn
!               INCPConvCoef(j,k)=((Exponent-one)*exp(-(stepn-j+one+alfi)
!     &            *Exponent)+exp(-(stepn-j+two+alfi)*Exponent))/Exponent
!     &            +(exp(-(stepn-j+alfi)*Exponent)-(Exponent+one)
!     &            *exp(-(stepn-j+one+alfi)*Exponent))/Exponent
!            enddo
!         endif
!      enddo
                
      return
      end
      
!
!.... calculate the time dependent INCPCoef for the coupled inflow BC
!
      subroutine CalcINCPCoef(timeINCP, stepn, srfIDList, numSrfs, y)   
      
      use incpBC
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer  srfIDList(0:MAXSURF), stepn, numSrfs, k
      real*8   y(nshg, ndof), timeINCP
      real*8   ENext(0:MAXSURF), PvenousNext(0:MAXSURF)
      
      INCPCoef = zero
      Enext = zero
      PvenousNext = zero
      call INCPint(timeINCP*(stepn+alfi), ENext, PvenousNext)
      Eadjust(1:numSrfs) = ENext(1:numSrfs)
!
!.... the following coefficients are used when the inductor is connected in series.
!      
      do k=1, numSrfs
         INCPCoef(1,k)=ValueVv(2,k) &
            +INCPConvCoef(stepn+2,k)*Eadjust(k) &
            +ValueVv(7,k)/timeINCP/alfi
         INCPCoef(2,k)=Eadjust(k)*VLV(stepn+1,k) &
            +Eadjust(k)*Qaorta(stepn+1,k)*INCPConvCoef(stepn+2,k) &
            -ValueVv(7,k)/timeINCP/alfi*Qaorta(stepn+1,k)
      enddo
!
!.... the following coefficients are used when the inductor is connectied in parallel.
!
!      do k=1, numSrfs
!         INCPCoef(1,k)=ValueVv(2,k)+alfi*Delt(1)/two *Eadjust(k)
!         INCPCoef(2,k)=Eadjust(k)*VLV(stepn+1,k)
!     &      +Eadjust(k)*Qaorta(stepn+1,k)*alfi*Delt(1)/two 
!         if (ValueVv(7,k) .ne. zero) then
!            INCPCoef(1,k) = INCPCoef(1,k)
!     &         -ValueVv(2,k)*INCPConvCoef(stepn+2,k)
!            INCPCoef(2,k) = INCPCoef(2,k)-ValueVv(2,k)*poldINCP(k)
!     &         -(PLV(1,k)-Paorta(1,k)-ValueVv(2,k)*Qaorta(1,k))
!     &         *exp(-timeINCP*(stepn+alfi)*ValueVv(2,k)/ValueVv(7,k)) 
!         endif
!      enddo      
      
      return 
      end
      
      
      subroutine calcINCPic(timeINCP, y, srfIDList, numSrfs)
      
      use incpBC      
      use boundarymodule, only: integrScalar, GetFlowQ           
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision 
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
         Paorta(1,1:numSrfs)= &
            Integral(1:numSrfs)/InflowArea(1:numSrfs)
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

!
!.... calculate the time dependent INCPCoef for the coupled inflow BC
!
      subroutine UpdHeartModel(timeINCP, y, srfIDList, &
         numSrfs, stepn) 

      use incpBC
      use boundarymodule, only: integrScalar
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer  srfIDList(0:MAXSURF), stepn, numSrfs, k
      real*8   y(nshg, ndof), timeINCP
      real*8   ENext(0:MAXSURF), Integral(0:MAXSURF)
      real*8   PvenousNext(0:MAXSURF)
      
      Enext = zero
      Integral = zero
      PvenousNext = zero
      call integrScalar(Integral, y(:,4), srfIDList, numSrfs)
      Paorta(stepn+1,1:numSrfs)= &
         Integral(1:numSrfs)/InflowArea(1:numSrfs)
      Qaorta(stepn+1,1:numSrfs)=QHistINCP(stepn+1,1:numSrfs)
      call INCPint(timeINCP*stepn, ENext, PvenousNext)
      Eadjust(1:numSrfs) = ENext(1:numSrfs)

      do k=1, numSrfs
         if (Qaorta(stepn+1,k) .ne. zero) then
            VLV(stepn+1,k)=VLV(stepn,k)+timeINCP/2*Qaorta(stepn,k) &
                   +timeINCP/2*Qaorta(stepn+1,k)
            PLV(stepn+1,k)=Eadjust(k)*VLV(stepn+1,k)
            QAV(stepn+1,k)=zero
         elseif (Qaorta(stepn+1,k) .eq. zero) then
            VLV(stepn+1,k)=VLV(stepn,k) &
               +timeINCP/2*(Qaorta(stepn,k)+QAV(stepn,k))
            PLV(stepn+1,k)=Eadjust(k)*VLV(stepn,k)
            QAV(stepn+1,k)=zero
            if (PLV(stepn+1,k) .lt. PvenousNext(k)) then
               QAV(stepn+1,k)=(PvenousNext(k)-PLV(stepn+1,k) &
                  +ValueVv(6,k)*QAV(stepn,k)/timeINCP)/ &
                  (ValueVv(1,k)+timeINCP/2*Eadjust(k) &
                  +ValueVv(6,k)/timeINCP)
               VLV(stepn+1,k)=VLV(stepn+1,k)+timeINCP/2*QAV(stepn+1,k)
               PLV(stepn+1,k)=Eadjust(k)*VLV(stepn+1,k)
               QHistINCP(stepn+1,k)=QAV(stepn+1,k)
            elseif (QAV(stepn,k) .gt. zero) then
                QAV(stepn+1,k)=(PvenousNext(k)-PLV(stepn+1,k) &
                  +ValueVv(6,k)*QAV(stepn,k)/timeINCP)/ &
                  (ValueVv(1,k)+timeINCP/2*Eadjust(k) &
                  +ValueVv(6,k)/timeINCP)
               VLV(stepn+1,k)=VLV(stepn+1,k)+timeINCP/2*QAV(stepn+1,k)
               PLV(stepn+1,k)=Eadjust(k)*VLV(stepn+1,k)
               QHistINCP(stepn+1,k)=QAV(stepn+1,k)              
            endif
         endif
      enddo
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
            (myrank .eq. zero)) then
         call OutputDataFile(Paorta(1:lstep+1,:),stepn+1,numSrfs, &
            'Paorta.dat',821)
         call OutputDataFile(Qaorta(1:lstep+1,:),stepn+1,numSrfs, &
            'Qaorta.dat',822)
         call OutputDataFile(PLV(1:lstep+1,:),stepn+1,numSrfs, &
            'PLV.dat',823)
         call OutputDataFile(VLV(1:lstep+1,:),stepn+1,numSrfs, &
            'VLV.dat',824)
         call OutputDataFile(QAV(1:lstep+1,:),stepn+1,numSrfs, &
            'QAV.dat',825)
      endif
        
      return
      end

!
!.... This subroutine writes FlowHist.dat and PressHist.dat files
!
      subroutine Updcalc(y, srfIDList, numSrfs)

      use calcFlowPressure
      use boundarymodule, only: integrScalar, GetFlowQ      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!      
      real*8   y(nshg, ndof), NewP(0:MAXSURF), NewQ(0:MAXSURF)
      integer  srfIDList(0:MAXSURF),  numSrfs

      call GetFlowQ(NewQ, y(:,1:3), srfIDList,numSrfs) !new flow at time n+1
      FlowHist(lstep+1,1:numSrfs) = NewQ(1:numSrfs)
      call integrScalar(NewP, y(:,4), srfIDList, numSrfs)
      PressHist(lstep+1,1:numSrfs) = NewP(1:numSrfs)/CalcArea(1:numSrfs)
      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
         (myrank .eq. zero)) then
         call OutputDataFile(FlowHist(1:lstep+1,:),lstep+1,numSrfs, &
            'FlowHist.dat',1004)
         call OutputDataFile(PressHist(1:lstep+1,:),lstep+1,numSrfs, &
            'PressHist.dat',1005)
      endif 

      return
      end

!
!.... This subroutine writes Lagrange Multipliers and errors in 
!.... LagrangeMultipliers.dat and LagrangeErrors.dat
!
      subroutine UpdateLagrangeCoef(y, col, row, srfIDList, numSrfs)

      use LagrangeMultipliers
      use boundarymodule, only: GetFlowQ      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
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
      PQLagrange(1:numSrfs,1)=Integral(1:numSrfs) &
         /LagProfileArea(1:numSrfs) 
      Integral = zero
      LagSwitch = 0 
    	call CalcNANBLagrange(col, row, y(:,1:3))
      call GetInnerProduct(InnerProduct, y(:,1:3), srfIDList, numSrfs)
      IPLagrange(1:numSrfs,1:3)=InnerProduct(1:numSrfs,1:3)
      do k=1, numSrfs
         NumOfData = (k-1)*3+1
         LagErrorHist(lstep+1,NumOfData)=abs(IPLagrange(k,1) &
            -two*QLagrange(k,1)*PQLagrange(k,1) &
            +QLagrange(k,1)**2*ProfileDelta(k))
         LagErrorHist(lstep+1,NumOfData+1)=abs(IPLagrange(k,2))
         LagErrorHist(lstep+1,NumOfData+2)=abs(IPLagrange(k,3))
            LagErrorHist(lstep+1,NumOfData:NumOfData+2)= &
            LagErrorHist(lstep+1,NumOfData:NumOfData+2) &
            *LagMeanFlow(k)
         LagHist(lstep+1,NumOfData:NumOfData+2)=Lag(k,1:3)
      enddo    

      if (((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)).and. &
            (myrank .eq. zero)) then
         NumOfData = numLagrangeSrfs*3
         call OutputDataFile(LagHist(1:lstep+1,:),lstep+1, NumOfData, &
            'LagrangeMultipliers.dat',801)
         call OutputDataFile(LagErrorHist(1:lstep+1,:),lstep+1, &
            NumOfData,'LagrangeErrors.dat',802)
      endif
!
      return
      end  
!
!.... this function calculates an initial condition of a constrained surface
!
      subroutine calcLagrangeic(srfIDList, numSrfs)
!      
      use LagrangeMultipliers
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision 

      integer  srfIDList(0:MAXSURF),  numSrfs 
      
      LagSwitch = 0 
      allocate(lhsLagL(9,nnz_tot,3))
      allocate(resL(numSrfs,3))
      allocate(LagAPproduct(nshg,3))
      lhsLagL = zero
      resL = zero   
      LagAPproduct = zero
      call MergeLagrangeParameters(srfIDList, numSrfs)
      ProfileDelta(1:numSrfs)=ProfileDelta(1:numSrfs) &
         /LagProfileArea(1:numSrfs)/LagProfileArea(1:numSrfs)
      do k=1, numSrfs
         LagMeanFlow(k)=two*LagProfileArea(k)/LagMeanFlow(k) &
            /LagMeanFlow(k)
         if (lstep .eq. zero) then
            LagHist(1,(k-1)*3+1:(k-1)*3+3)=Lagold(k,1:3)
         elseif (lstep .gt. zero) then
            Lagold(k,1:3)=LagHist(lstep+1,(k-1)*3+1:(k-1)*3+3)
         endif
      enddo
!            
      return 
      end
!
!.... this function calculates an area and plane vectors of a constrained surface
!
      subroutine MergeLagrangeParameters(srfIDList, numSrfs)
!      
      use LagrangeMultipliers
!      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision 
      include "mpif.h"   
!            
      integer  srfIDList(0:MAXSURF),  numSrfs 
      real*8   VectMag(3), Inplane1, Inplane2, Inplane3, InplaneNorm
      real*8, allocatable, dimension (:) :: TotalArea
      real*8, allocatable, dimension (:,:,:) :: InplaneVectors
!      
      allocate(TotalArea(numSrfs))
      allocate(InplaneVectors(3,3,numSrfs))
      TotalArea = zero
      InplaneVectors = zero
      call MPI_ALLREDUCE (LagProfileArea, TotalArea, numSrfs, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)  
      LagProfileArea(1:numSrfs)=TotalArea(1:numSrfs)
      TotalArea = zero
      call MPI_ALLREDUCE (ProfileDelta, TotalArea, numSrfs, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)  
      ProfileDelta(1:numSrfs)=TotalArea(1:numSrfs)
      InplaneVectors = zero
      
      do i=1,3
         do j=1,3
            call MPI_ALLREDUCE(LagInplaneVectors(i,j,:), &
               InplaneVectors(i,j,:), numSrfs,   &
               MPI_DOUBLE_PRECISION, MPI_SUM, INEWCOMM,ierr)
         enddo
      enddo
      LagInplaneVectors = InplaneVectors
      do k=1, numSrfs
         do i=1,3
            VectMag(i)=sqrt(LagInplaneVectors(1,i,k)**2+ &
               LagInplaneVectors(2,i,k)**2+LagInplaneVectors(3,i,k)**2)
         enddo
         if ( VectMag(1) .gt. zero .and. VectMag(2) .gt. zero) then
            LagInplaneVectors(1:3,1,k) = LagInplaneVectors(1:3,1,k) &
               /VectMag(1)
            LagInplaneVectors(1:3,2,k)=LagInplaneVectors(1:3,2,k) &
               /VectMag(2)
            Inplane1=-LagInplaneVectors(2,2,k)*LagInplaneVectors(3,1,k) &
               +LagInplaneVectors(2,1,k)*LagInplaneVectors(3,2,k)
            Inplane2=-LagInplaneVectors(1,1,k)*LagInplaneVectors(3,2,k) &
               +LagInplaneVectors(1,2,k)*LagInplaneVectors(3,1,k)
            Inplane3=-LagInplaneVectors(1,2,k)*LagInplaneVectors(2,1,k) &
               +LagInplaneVectors(1,1,k)*LagInplaneVectors(2,2,k)
            InplaneNorm=one/sqrt(Inplane1**2+Inplane2**2+Inplane3**2)
            LagInplaneVectors(1,3,k)=Inplane1*InplaneNorm
            LagInplaneVectors(2,3,k)=Inplane2*InplaneNorm
            LagInplaneVectors(3,3,k)=Inplane3*InplaneNorm
         endif
      enddo
!            
      return 
      end      
!  
! !.........function that integrates a scalar over a boundary
! !
!       subroutine integrScalar(scalInt,scal,srfIdList,numSrfs)

!       use pvsQbi !brings ndsurf, NASC

!       use phcommonvars
!       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!       include "mpif.h"
      
!       integer   srfIdList(0:MAXSURF), numSrfs, irankCoupled, i, k
!       real*8    scal(nshg), scalInt(0:MAXSURF), scalIntProc(0:MAXSURF)
      
!       scalIntProc = zero
!       do i = 1,nshg
!         if(numSrfs.gt.zero) then
!           do k = 1,numSrfs
!             irankCoupled = 0
!             if (srfIdList(k).eq.ndsurf(i)) then
!               irankCoupled=k
!               scalIntProc(irankCoupled) = scalIntProc(irankCoupled) &
!                                   + NASC(i)*scal(i)
!             endif      
!           enddo       
!         endif
!       enddo
! !      
! !     at this point, each scalint has its "nodes" contributions to the scalar
! !     accumulated into scalIntProc. Note, because NASC is on processor this
! !     will NOT be the scalar for the surface yet
! !
! !.... reduce integrated scalar for each surface, push on scalInt
! !
!         npars=MAXSURF+1
!       call MPI_ALLREDUCE (scalIntProc, scalInt(:), npars, &
!               MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)  
   
!       return
!       end

!  
!.........function that outputs an input data array
!
      subroutine OutputDataFile(DataFile, nrows, ncolms, Filename, &
         UnitNumber)

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      character(*) Filename
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
