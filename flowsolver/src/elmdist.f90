      subroutine ElmDist(u, x, xdist, normvect)
      
         use pointer_data  ! brings in the pointers for the blocked arrays
         use measureWallDistance       
         use deformableWall
         
         use phcommonvars
         IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
         
         dimension x(numnp,nsd),u(nshg,nsd)
         dimension xdist(nshg)
         dimension lnode(27)
         
         real*8 tempnv1(nsd)
         real*8 tempnv2(nsd)
         real*8 normvect(numnp,nsd)
         
         real*8 closestPt1(nsd),closestPt2(nsd)
         
         real*8 tempPt(nsd)
         real*8 tempDistSq1,tempDistSq2,cycleTime,obsInterval,intTime
         real*8 alphaObs
         
         integer ii,jj,nn,nPer,nObsInt,obsFr1,obsFr2
         
!         
!.....compute the number of periods that has passed
!     according to the cycle length parameter that was set
!     in observed.dat         
!
         nPer = (lstep)*Delt(1)/cycleLength
         
!
!.... compute the physical time
!
                  
         cycleTime = (lstep)*Delt(1)-nPer*cycleLength      
         
!
!.... compute the interval that we are in right now
!         

         obsInterval = cycleLength / numDataFrames ! interval between data        
         nObsInt = cycleTime / obsInterval
         obsFr1 = nObsInt + 1
         obsFr2 = nObsInt + 2
         intTime = cycleTime - nObsInt*obsInterval
         
         if (obsFr2 .gt. numDataFrames) then
            obsFr2 = 1
         end if
         
         alphaObs = 1-intTime/obsInterval
         
         if (myrank .eq. 0) then 
            write(*,*) "dtime",obsFr1,obsFr2,(1-alphaObs)
         end if    

!         
!        loop over the deformable wall nodes
!

         do ii = 1, nwnp_EWB

            ! TODO: assumption is that numnp = nshg
            tempPt = x(mWNodes_EWB%p(ii),:)+ &
                     u(mWNodes_EWB%p(ii),:)
     
            call dm_cpmeshp3(obsFr1,tempPt, &
                             closestpt1,tempDistSq1, &
                             tempnv1)
            call dm_cpmeshp3(obsFr2,tempPt, &
                             closestpt2,tempDistSq2, &
                             tempnv2)
                      
            tempDistSq1 = sign(sqrt(abs(tempDistSq1)),tempDistSq1) 
     
            tempDistSq2 = sign(sqrt(abs(tempDistSq2)),tempDistSq2) 
     
            xdist(mWNodes_EWB%p(ii)) =  &
               alphaObs*tempDistSq1+(1-alphaObs)*tempDistSq2
     
            normvect(mWNodes_EWB%p(ii),:) =  &
               alphaObs*tempnv1+(1-alphaObs)*tempnv2
   
         end do
                  
      end subroutine
