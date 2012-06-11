      subroutine getdmc (y,      shgl,      shp,  &
                         iper,   ilwork,     &
                         nsons,  ifath,     x)

      use pointer_data

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
!                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
!                    Shpf and shglf are the shape funciotns and their 
!                    gradient evaluated using the quadrature rule desired 
!                    for computing the dmod. Qwt contains the weights of the 
!                    quad. points.  



      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

!
      dimension fres(nshg,24),         fwr(nshg), &
                strnrm(nshg),         cdelsq(nshg), &
                xnum(nshg),           xden(nshg), &
                xmij(nshg,6),         xlij(nshg,6), &
                xnude(nfath,2),        xnuder(nfath,2), &
                nsons(nfath), &
                strl(numel,maxnint),            &
                y(nshg,5),  &
                ifath(nshg),          iper(nshg), &
                ilwork(nlwork), & !        xmudmi(numel,ngauss), &
                x(numnp,3), &
                shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT)    
!$$$     &          ,xnutf(nfath)  must be uncommmented for diags at bottom
!
!
!   setup the weights for time averaging of cdelsq (now in quadfilt module)
!
      denom=max(1.0d0*(lstep),one)
      if(dtavei.lt.0) then
         wcur=one/denom
      else
         wcur=dtavei
      endif  
      whist=1.0-wcur
!
!  hack in an interesting velocity field (uncomment to test dmod)
!
!      y(:,5) = 1.0  
!      y(:,2) = 2.0*x(:,1) - 3*x(:,2) 
!      y(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
!      y(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3) 
!      y(:,1) = Rgas * y(:,5) ! Necessary to make model suitable suitable
                             ! for the incompressible case.
!

     
      fres = zero
!
!       y(:,5) = 1.0  
!       y(:,1) = 2.0*x(:,1) - 3*x(:,2) 
!       y(:,2) = 3.0*x(:,1) + 4.0*x(:,2)
!       y(:,3) = 4.0*x(:,1) + x(:,2) + x(:,3) 
!       y(:,4) = 1.0    ! Necessary to make model suitable suitable
                             ! for the incompressible case.
!


      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call asithf (y, x, strl(iel:inum,:), mien(iblk)%p, fres,  &
                     shglf(lcsyst,:,1:nshl,:), &
                     shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
!
 
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)

        if (ngaussf .ne. ngauss) then
        call getstrl (y, x,      mien(iblk)%p,   &
                     strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:), &
                     shp(lcsyst,1:nshl,:))
        endif

      enddo
!
!
! must fix for abc and dynamic model
!      if(iabc==1)   !are there any axisym bc's
!     &      call rotabc(res, iBC,  'in ')
!
      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
!proc-masters = proc-masters + proc-slaves
! 
! account for periodicity in filtered variables
!
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:) ! masters = masters + slaves
        endif
      enddo
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)  !slaves get copy of the complete master
        endif
      enddo

      if (nfath .eq. 0) then
         if(numpe>1)   call commu (fres, ilwork, 24, 'out')
! again we will need a  call
!      if(iabc==1)   !are there any axisym bc's
!     &   rotabc(y, iBC,  'out')
      endif

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)  !"solve"
      enddo
!     fres(:,24) = fres(:,24) * fres(:,23)
!
!.....at this point fres is really all of our filtered quantities
!     at the nodes
!

      strnrm = sqrt(  &
        two * (fres(:,10)**2 + fres(:,11)**2 + fres(:,12)**2) &
        + four * ( fres(:,13)**2 + fres(:,14)**2 + fres(:,15)**2 ) )

      fwr = fwr1 * fres(:,22) * strnrm



!... There is a difference in xmij(:,1:3) between the compressible 
!... and incompressible forms of the model.


      if(matflg(1,1).eq.0) then ! compressible

      xmij(:,1) = -fwr &
                   * pt33 * (two*fres(:,10) - fres(:,11) - fres(:,12)) &
                   + pt33 * (two*fres(:,16) - fres(:,17) - fres(:,18))
      xmij(:,2) = -fwr &
                   * pt33 * (two*fres(:,11) - fres(:,10) - fres(:,12)) &
                   + pt33 * (two*fres(:,17) - fres(:,16) - fres(:,18))
      xmij(:,3) = -fwr &
                   * pt33 * (two*fres(:,12) - fres(:,10) - fres(:,11)) &
                   + pt33 * (two*fres(:,18) - fres(:,16) - fres(:,17))

      else

      xmij(:,1) = -fwr &
                   * fres(:,10) + fres(:,16)
      xmij(:,2) = -fwr &
                   * fres(:,11) + fres(:,17) 
      xmij(:,3) = -fwr &
                   * fres(:,12) + fres(:,18) 

      endif      

      xmij(:,4) = -fwr * fres(:,13) + fres(:,19)
      xmij(:,5) = -fwr * fres(:,14) + fres(:,20)
      xmij(:,6) = -fwr * fres(:,15) + fres(:,21)

      fres(:,22) = one / fres(:,22)

      xlij(:,1) = fres(:,4) - fres(:,1) * fres(:,1) * fres(:,22)
      xlij(:,2) = fres(:,5) - fres(:,2) * fres(:,2) * fres(:,22)
      xlij(:,3) = fres(:,6) - fres(:,3) * fres(:,3) * fres(:,22)
      xlij(:,4) = fres(:,7) - fres(:,1) * fres(:,2) * fres(:,22)
      xlij(:,5) = fres(:,8) - fres(:,1) * fres(:,3) * fres(:,22)
      xlij(:,6) = fres(:,9) - fres(:,2) * fres(:,3) * fres(:,22)

      xnum =        xlij(:,1) * xmij(:,1) + xlij(:,2) * xmij(:,2)  &
                                          + xlij(:,3) * xmij(:,3) &
           + two * (xlij(:,4) * xmij(:,4) + xlij(:,5) * xmij(:,5) &
                                          + xlij(:,6) * xmij(:,6))
      xden =        xmij(:,1) * xmij(:,1) + xmij(:,2) * xmij(:,2)  &
                                          + xmij(:,3) * xmij(:,3) &
           + two * (xmij(:,4) * xmij(:,4) + xmij(:,5) * xmij(:,5) &
                                          + xmij(:,6) * xmij(:,6))
      xden = two * xden

!  zero on processor periodic nodes so that they will not be added twice
        do j = 1,numnp
          i = iper(j)
          if (i .ne. j) then
            xnum(j) = zero
            xden(j) = zero
          endif
        enddo

      if (numpe.gt.1 .and. nsons(1).gt.1) then

         numtask = ilwork(1)
         itkbeg = 1
       
! zero the nodes that are "solved" on the other processors  
         do itask = 1, numtask

            iacc   = ilwork (itkbeg + 2)
            numseg = ilwork (itkbeg + 4)

            if (iacc .eq. 0) then
               do is = 1,numseg
                  isgbeg = ilwork (itkbeg + 3 + 2*is)
                  lenseg = ilwork (itkbeg + 4 + 2*is)
                  isgend = isgbeg + lenseg - 1
                  xnum(isgbeg:isgend) = zero
                  xden(isgbeg:isgend) = zero
               enddo
            endif
            
            itkbeg = itkbeg + 4 + 2*numseg
            
         enddo
         
      endif
!
! Description of arrays.   Each processor has an array of length equal
! to the total number of fathers times 2 xnude(nfathers,2). One to collect 
! the numerator and one to collect the denominator.  There is also an array
! of length nshg on each processor which tells the father number of each
! on processor node, ifath(nnshg).  Finally, there is an arry of length
! nfathers to tell the total (on all processors combined) number of sons
! for each father. 
!
!  Now loop over nodes and accumlate the numerator and the denominator
!  to the father nodes.  Only on processor addition at this point.
!  Note that serrogate fathers are collect some for the case where some
!  sons are on another processor
!
      ihomog=1
      if(maxval(nsons).eq.1) ihomog=0
      if(ihomog.eq.1) then
         xnude = zero
         do i = 1,nshg
            xnude(ifath(i),1) = xnude(ifath(i),1) + xnum(i)
            xnude(ifath(i),2) = xnude(ifath(i),2) + xden(i)
         enddo
      else                      ! no homogeneity assumed
         xnude(:,1)=xnum(:)
         xnude(:,2)=xden(:)
      endif

!
! Now  the true fathers and serrogates combine results and update
! each other.
!
! ONLY DO THIS IF  1) multiprocessors AND 2) homogenous directions exist 
  
      if(numpe .gt. 1 .and. ihomog.eq.1 )then
         call drvAllreduce(xnude, xnuder,2*nfath)
!
!  xnude is the sum of the sons for each father on this processor
!
!  xnuder is the sum of the sons for each father on all processor combined
!  (the same as if we had not partitioned the mesh for each processor)
!
!   For each father we have precomputed the number of sons (including
!   the sons off processor). 
!
!   Now divide by number of sons to get the average (not really necessary
!   for dynamic model since ratio will cancel nsons at each father)
!
!         xnuder(:,1) = xnuder(:,1) ! / nsons(:)
!         xnuder(:,2) = xnuder(:,2) ! / nsons(:)
!
!  the next line is c \Delta^2
!
            numNden(:,1) = whist*numNden(:,1)+wcur*xnuder(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnuder(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)
!  note that we have whist and wcur in here to allow for both time
!  averaging to be used in conjunction with spatial homogenous averaging

      else
!     
!     the next line is c \Delta^2, not nu_T but we want to save the
!     memory
!     

!$$$            write(400+myrank,555) (numNden(j*500,1),j=1,5)
!$$$            write(410+myrank,555) (numNden(j*500,2),j=1,5)
!$$$            write(500+myrank,555) (numNden(j+500,1),j=1,5)
!$$$            write(510+myrank,555) (numNden(j+500,2),j=1,5)
!$$$
!$$$            write(430+myrank,555) (xnude(j*500,1),j=1,5)
!$$$            write(440+myrank,555) (xnude(j*500,2),j=1,5)
!$$$            write(530+myrank,555) (xnude(j+500,1),j=1,5)
!$$$            write(540+myrank,555) (xnude(j+500,2),j=1,5)

            numNden(:,1) = whist*numNden(:,1)+wcur*xnude(ifath(:),1)
            numNden(:,2) = whist*numNden(:,2)+wcur*xnude(ifath(:),2)
            cdelsq(:) = numNden(:,1) / (numNden(:,2) + 1.d-09)

! 
!  to get started we hold cdelsq fixed
!
!            cdelsq(:) = 2.27e-4
!            cdelsq(:) = 0

            
!$$$            
!$$$            write(450+myrank,555) (cdelsq(j*500),j=1,5)
!$$$            write(460+myrank,555) (y(j*500,1),j=1,5)
!$$$            write(470+myrank,555) (y(j*500,2),j=1,5)
!$$$            write(480+myrank,555) (y(j*500,3),j=1,5)
!$$$            write(490+myrank,555) (strnrm(j*500),j=1,5)
!$$$
!$$$            write(550+myrank,555) (cdelsq(j+500),j=1,5)
!$$$            write(560+myrank,555) (y(j+500,1),j=1,5)
!$$$            write(570+myrank,555) (y(j+500,2),j=1,5)
!$$$            write(580+myrank,555) (y(j+500,3),j=1,5)
!$$$            write(590+myrank,555) (strnrm(j+500),j=1,5)
!$$$
!$$$            call flush(400+myrank)
!$$$            call flush(410+myrank)
!$$$            call flush(430+myrank)
!$$$            call flush(440+myrank)
!$$$            call flush(450+myrank)
!$$$            call flush(460+myrank)
!$$$            call flush(470+myrank)
!$$$            call flush(480+myrank)
!$$$            call flush(490+myrank)
!$$$            call flush(500+myrank)
!$$$            call flush(510+myrank)
!$$$            call flush(530+myrank)
!$$$            call flush(540+myrank)
!$$$            call flush(550+myrank)
!$$$            call flush(560+myrank)
!$$$            call flush(570+myrank)
!$$$            call flush(580+myrank)
!$$$            call flush(590+myrank)
      endif
 555  format(5(2x,e14.7))

! $$$$$$$$$$$$$$$$$$$$$$$$$$$
      tmp1 =  MINVAL(cdelsq)
      tmp2 =  MAXVAL(cdelsq)
      if(numpe>1) then
         call MPI_REDUCE (tmp1, tmp3, 1,MPI_DOUBLE_PRECISION, &
              MPI_MIN, master, INEWCOMM, ierr)
         call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION, &
              MPI_MAX, master, INEWCOMM, ierr)
         tmp1=tmp3
         tmp2=tmp4
      endif
      if (myrank .EQ. master) then !print CDelta^2 range
         write(34,*)lstep,tmp1,tmp2
         call flush(34)
      endif
! $$$$$$$$$$$$$$$$$$$$$$$$$$$

      if (myrank .eq. master) then
         write(*,*)'xnut=',sum(cdelsq)/nshg
         write(*,*) 'cdelsq=', cdelsq(1),cdelsq(2)        
      endif

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call scatnu (mien(iblk)%p, strl(iel:inum,:),  &
              mxmudmi(iblk)%p,cdelsq,shp(lcsyst,1:nshl,:))
      enddo
!     $$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$  tmp1 =  MINVAL(xmudmi)
!$$$  tmp2 =  MAXVAL(xmudmi)
!$$$  if(numpe>1) then
!$$$  call MPI_REDUCE (tmp1, tmp3, 1, MPI_DOUBLE_PRECISION,
!$$$  &                 MPI_MIN, master, INEWCOMM, ierr)
!$$$  call MPI_REDUCE (tmp2, tmp4, 1, MPI_DOUBLE_PRECISION,
!$$$  &                 MPI_MAX, master, INEWCOMM, ierr)
!$$$      tmp1=tmp3
!$$$  tmp2=tmp4
!$$$  endif
!$$$  if (myrank .EQ. master) then
!$$$  write(35,*) lstep,tmp1,tmp2
!$$$  call flush(35)
!$$$  endif
! $$$$$$$$$$$$$$$$$$$$$$$$$$$

!$$$c
!$$$c  if flag set, write a restart file with info (reuse xmij's memory)
!$$$c
!$$$      if(irs.eq.11) then
!$$$         lstep=999
!$$$         xmij(:,1)=xnum(:)
!$$$         xmij(:,2)=xden(:)
!$$$         xmij(:,3)=cdelsq(:)
!$$$         xmij(:,5)=xlij(:,4)    !leave M_{12} in 4 and put L_{12} here
!$$$         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
!$$$         stop
!$$$      endif
!$$$      if(lhs.eq.1) then
!$$$         lstepsafe=lstep
!$$$         lstep=999
!$$$
!$$$         xmij(:,1)=xnum(:)
!$$$         xmij(:,2)=xden(:)
!$$$         xmij(:,3)=cdelsq(:)
!$$$         xmij(:,4)=xmij(:,6)    !put M_{23} in 4 
!$$$         xmij(:,5)=xlij(:,6)    !put L_{23} here
!$$$         call restar('out ',xmij,xlij) !also dump all of L_{ij} in ac
!$$$         lstep=lstepsafe
!$$$      endif
!$$$c
!$$$c  local clipping moved to scatnu with the creation of mxmudmi pointers
!$$$c
!$$$c$$$      rmu=datmat(1,2,1)
!$$$c$$$      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
!$$$c$$$      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0
!$$$c      stop !uncomment to test dmod
!$$$c
!$$$
!$$$
!$$$c  write out the nodal values of xnut (estimate since we don't calc strain
!$$$c  there and must use the filtered strain).
!$$$c
!$$$
!$$$      if ((irs .ge. 1) .and. (mod(lstep, ntout) .eq. 0)) then
!$$$c
!$$$c  collect the average strain into xnude(2)
!$$$c
!$$$         xnude(:,2) = zero
!$$$         do i = 1,numnp
!$$$            xnude(ifath(i),2) = xnude(ifath(i),2) + strnrm(i)
!$$$         enddo
!$$$
!$$$         if(numpe .gt. 1 .and. ihomog.eq.1 )then
!$$$            call drvAllreduce(xnude(:,2), xnuder(:,2),nfath)
!$$$         else
!$$$            xnuder=xnude
!$$$         endif
!$$$c     
!$$$c          nut= cdelsq    * |S|
!$$$c 
!$$$         xnutf=xnuder(:,1)*xnuder(:,2)/nsons(:) ! off by one
!$$$c
!$$$c  collect the x and y coords into xnude
!$$$c
!$$$         xnude = zero
!$$$         do i = 1,numnp
!$$$            xnude(ifath(i),1) = xnude(ifath(i),1) + x(i,3)
!$$$            xnude(ifath(i),2) = xnude(ifath(i),2) + x(i,2)
!$$$         enddo
!$$$
!$$$         if(numpe .gt. 1 .and. nsons(1).ne.1 )then
!$$$            call drvAllreduce(xnude, xnuder,2*nfath)
!$$$            xnuder(:,1)=xnuder(:,1)/nsons(:)
!$$$            xnuder(:,2)=xnuder(:,2)/nsons(:)
!$$$         else
!$$$            xnuder=xnude
!$$$         endif
!$$$c
!$$$c  xnude is the sum of the sons for each father on this processor
!$$$c
!$$$         if((myrank.eq.master)) then
!$$$            do i=1,nfath      ! cdelsq   * |S|
!$$$               write(444,*) xnuder(i,1),xnuder(i,2),xnutf(i)
!$$$            enddo
!$$$            call flush(444)
!$$$         endif
!$$$      endif

      return
      end
