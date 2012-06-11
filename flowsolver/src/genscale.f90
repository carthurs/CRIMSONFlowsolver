      subroutine genscale(y, x, iBC)
!
!----------------------------------------------------------------------
! This subroutine calculate the y^+ and eta at inflow and internal face.
! From these generate the scaling for the inflow data.
!
! input:
!  iBC    (numnp)               : boundary condition code
!  x      (numnp,nsd)           : node coordinates
!
! output:
!  y      (numnp,ndof)          : initial values of Y variables
!
!
! Elaine Bohr december 2001
!----------------------------------------------------------------------
!
       use spebc
!       use pointer_data
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
       include "mpif.h"
       !include "auxmpi.h"
!
       dimension y(numnp,ndof),   iBC(numnp),   &
                 x(numnp,nsd), velbarR(nfint,nflow)
       dimension ifath(numnp),  velbarl(nelint,nshl,nflow), &
                 v1(nfint),       ymapped(numnp,ndof), &
      		 shapef(nshl),	shgradl(nshl,nsd), &
      		 xsi(nsd), yintl(nelint,nshl,nflow), &
      		 flucl(nelint,nshl,nflow), &
      		 ubarintl(nelint,nshl,nflow), &
      		 fluc1(npin,nflow), fluc2(npin,nflow), &
      		 ubar1(npin,nflow), ubar2(npin,nflow)
       integer   element, dir

       real*8    ymax, displTi, displTr, correction
       real*8	 freestream(nflow)
       save deltaint

       
!	return
        ymapped(:,2:4)=y(:,1:3)
	ymapped(:,1)=y(:,4)
	ymapped(:,5)=y(:,5)
	
	ubar2 = 0
	fluc2 = 0
	
	ymax = xyn(nfint)

!
! .... Localizing the solution vector on virtual plane
!

        do i = 1, nelint
          do j = 1, 3
            yintl(i,:,j+1) = y(ien2D(i,:),j)
          enddo
          yintl(i,:,1) = y(ien2D(i,:),4)
	  if(nflow.gt.4) then
            do j = 5, nflow
              yintl(i,:,j) = y(ien2D(i,:),j)
            enddo
          endif
	enddo  

!
! .... Finding averaged velocity in spanwise direction
!      for the virtual plane
!

	do i=1,nfint
	  velbarR(i,:)=0
	  do j=1,imax(i)+1
	    call shptet(ipord,xsinfin(i,j,:),shapef(:),shgradl(:,:))
	    do k=1,nshl
	      velbarR(i,:)=velbarR(i,:)  &
      		+ yintl(elcnfin(i,j),k,:)*shapef(k)
	    enddo
	  enddo
	  velbarR(i,:)=velbarR(i,:) / (imax(i)+1)
	enddo
	 
!
! .... Label the nodes that near the BL thickness
! 

       if (thetag.eq.0.0) then
         dir = 2
       else
         dir = 4
       endif

       v1(1)=10.0
       do i=2,nfint+1
          v1(i)=velbarR(i-1,dir)-0.99*vel
          if((v1(i).gt.0).and.(v1(i-1).le.0)) then
             label=i-1
             go to 200
          endif
       enddo
       label=i-1

!     
!.... Find the BL thickness by means of finding the y coord
!     

 200   continue
       dv=velbarR(label,dir)-velbarR(label-1,dir)
       dy=xyn(label)-xyn(label-1)

!     
! .... Current calculation of bl thickness at recycle plane
!

       if(istep.ne.0) then
          dlast=deltaint
          deltaint=xyn(label-1) &
      		+ dy*(0.99*vel-velbarR(label-1,dir))/dv
     
!
! .... Early transients cause jumpy delta, smooth it.
!

          deltaint=min(1.05*dlast,max(deltaint,0.95*dlast))
       else
          deltaint=xyn(label-1) &
      		+ dy*(0.99*vel-velbarR(label-1,dir))/dv
       endif

!
! .... Deltaint is now the ratio of BL thickness at the interior plane
!      to the BL thickness at the inlet plane
!

       deltaint=min(two*rbltin,max(deltaint,pt5*rbltin)) 
       rdelta = deltaint/rbltin
       
!
! .... Finding freestream solutions
!
	
	freestream = 0
	icount = 0
	do i=1, nfint
	  if (xyn(i).ge.deltaint) then
	    freestream(:) = freestream(:) + velbarR(i,:)
	    icount = icount + 1 
	  endif
	enddo
	freestream = freestream / icount
	
!
! .... Putting the freestream values into the average outside the BLT
!

	do i=1, nfint
	  if (xyn(i).ge.deltaint) then
	    velbarR(i,:) = freestream(:)
	  endif
	enddo

!
! .... Localizing the averaged velocity found above
!

	do i=1,nelint
	  do k=1,nshl
	    do j=1,nfint-1
	      if (thetag.eq.0.0) then
	        if ((x(ien2D(i,k),2).ge.xyn(j)) .and. &
      		    (x(ien2D(i,k),2).le.(xyn(j+1)+0.000001))) then
		  tmp = (x(ien2D(i,k),2) - xyn(j)) / &
      			(xyn(j+1) - xyn(j))
     		  do l=1,nflow
		     velbarl(i,k,l) =  &
      		            (velbarR(j+1,l) - velbarR(j,l)) *  &
                             tmp + velbarR(j,l)
	          enddo
		endif
	      else
	        if ((xcyl(ien2D(i,k),1).ge.xcyl(nrint(j+1),1)) .and. &
      		    (xcyl(ien2D(i,k),1).le.xcyl(nrint(j),1))) then
                  tmp = (xcyl(ien2D(i,k),1) - xcyl(nrint(j+1),1)) /  &
      	              (xcyl(nrint(j),1) - xcyl(nrint(j+1),1))
     		  do l=1,nflow
		     velbarl(i,k,l) =  &
      		            (velbarR(j,l) - velbarR(j+1,l)) *  &
                             tmp + velbarR(j+1,l)
	          enddo
		endif
     	      endif
	    enddo
	  enddo
	enddo
	
!
! --- For now only Blasius is coded ---
!
       
!
! .... Calculate fluctuations on elements of internal plane
!

!       flucl = yintl - velbarl

!
! .... Calculate mean values on elements of internal plane
!

       ubarintl = velbarl 

!
! .... Calculating the coordinates of the point from where the
!      solution will be projected to the inlet plane
!


	do i=1,npin
	
!
! .... Cartesian coodinate system
!

	  if (thetag.eq.0.0) then
	    xts1 = x(nen1(i),1) + plandist
	    if (xynin(i)*rdelta.gt.ymax) then
	      xts2 = ymax
	    else  
	      xts2 = xynin(i)*rdelta
	    endif
	    xts3 = x(nen1(i),3)

!
! .... Cylindrical coordinate system
!

	  else
	    if (xynin(i).le.0.00001) then
	      xts1 = (radcyl-xynin(i)*rdelta*sang-tolerence) &
       		*cos(xcyl(nen1(i),2))
	      xts2 = (radcyl-xynin(i)*rdelta*sang-tolerence) &
      		*sin(xcyl(nen1(i),2))
	      xts3 = (aR-(radcyl-xynin(i)*rdelta*sang-tolerence) &
      	       * (xnrml*cos(xcyl(nen1(i),2)) &
      	       +  ynrml*sin(xcyl(nen1(i),2))))/znrml
            elseif (xynin(i)*rdelta.gt.ymax) then
	      xts1 = (radcyl-ymax*sang) &
       		*cos(xcyl(nen1(i),2))
	      xts2 = (radcyl-ymax*sang) &
      		*sin(xcyl(nen1(i),2))
	      xts3 = (aR-(radcyl-ymax*sang) &
      	       * (xnrml*cos(xcyl(nen1(i),2)) &
      	       +  ynrml*sin(xcyl(nen1(i),2))))/znrml
	    else
	      xts1 = (radcyl-xynin(i)*rdelta*sang) &
       		*cos(xcyl(nen1(i),2))
	      xts2 = (radcyl-xynin(i)*rdelta*sang) &
      		*sin(xcyl(nen1(i),2))
	      xts3 = (aR-(radcyl-xynin(i)*rdelta*sang) &
      	       * (xnrml*cos(xcyl(nen1(i),2)) &
      	       +  ynrml*sin(xcyl(nen1(i),2))))/znrml
            endif
	  endif
	  
!
! .... Searching for the appropriate element
!

     	  call elem_search(xintl, xts1, xts2, xts3, &
      			   xsi(:), element, 2)
      	  call shptet(ipord,xsi(:),shapef(:),shgradl(:,:))	  

!
! .... Calculating the average velocity and fluctuations
!      for the inlet plane
!

	  do k=1,nshl
	    fluc2(i,:)= 0 !fluc2(i,:) + flucl(element,k,:)*shapef(k)
	    ubar2(i,:)=ubar2(i,:) + ubarintl(element,k,:)*shapef(k)
	  enddo
	enddo  
	        

!$$$c
!$$$c keep freestream values set through averages
!$$$c
!$$$         ubaro=0
!$$$	 tbaro=0
!$$$         icount=0
!$$$         do i=1,nfin
!$$$            if(yin(i).ge.rbltin) then
!$$$               nzl=nsons(i)  !Elaine
!$$$               nzb=ienson1(i,1)
!$$$               nze=nzb+nzl-1
!$$$	       tbaro=tbaro+2.0*ubar2(i,5)+sum(fluc2(nzb:nze,5))
!$$$               ubaro=ubaro               +sum(fluc2(nzb:nze,2))
!$$$               icount=icount+nzl
!$$$            endif
!$$$         enddo
!$$$         
!$$$c     alternative to myway
!$$$c     
!$$$         ubaro=ubaro/icount
!$$$	 tmeaninflow=0.0625097048890964
!$$$	 fact= tmeaninflow/(tbaro/icount)
!$$$	 if (fact.ge. 0.9999999 .and. fact.le.1.0000001) fact = 1.0
!$$$         
!$$$         do i=1,nfin
!$$$            if(yin(i).ge.rbltin) then
!$$$               ubar2(i,2)=1.0-ubaro
!$$$            endif
!$$$         enddo
         fact = 1.0 
	 rvscal = 1.0
	 
!
! .... Putting the freestream value outside the BLT into ubar2
!

	do i = 1, npin
	  if (xynin(i).ge.rbltin) then
	    ubar2(i,:) = freestream(:)
!	    ubar2(i,dir) = 1.0
	    fluc2(i,:) = 0
	  endif
	enddo

!$$$c
!$$$c .... For the cylindrical case the freestream velocity needs
!$$$c      to be corrected for the blockage
!$$$c
!$$$
!$$$	if (thetag.ne.0.0) then
!$$$	  displTi = 0.0
!$$$	  displTr = 0.0
!$$$	  do i = 2, nfint
!$$$
!$$$c
!$$$c .... Displacement thickness for inlet plane
!$$$c
!$$$
!$$$	    displTi = displTi + (1 - y(nrint(i),3))
!$$$     &              * (xyn(i) - xyn(i-1)) * (radcyl - xyn(i))
!$$$
!$$$c
!$$$c .... Displacement thickness for recycle plane
!$$$c
!$$$
!$$$     	    displTr = displTr + (1 - velbarR(i,4))
!$$$     &              * (xyn(i) - xyn(i-1)) * (radcyl - xyn(i))
!$$$	  enddo
!$$$c	  displTi = radcyl - sqrt(radcyl*radcyl - displTi)
!$$$c	  displTr = radcyl - sqrt(radcyl*radcyl - displTr)
!$$$	  correction = (radcyl*radcyl - displTr)
!$$$     &               / (radcyl*radcyl - displTi)
!$$$        else
	  correction = 1.0
!$$$	endif

!     
! .... Scaled plane extraction boundary condition
!     

         ymapped(nen1(1:npin),1)= correction * (ubar2(:,1)+fluc2(:,1))
         ymapped(nen1(1:npin),2)= correction * &
              (ubar2(:,2)+fluc2(:,2)) !myway *factu
         ymapped(nen1(1:npin),3)= correction * &
              (ubar2(:,3)+fluc2(:,3))*rvscal
         ymapped(nen1(1:npin),4)= correction * (ubar2(:,4)+fluc2(:,4))
      ymapped(nen1(1:npin),5)= correction * fact*(ubar2(:,5)+fluc2(:,5))



!     
! .... Ready to put the solution on the inflow plane 
!

      if(intpres.eq.1) then     !interpolating pressure at inflow
         where (btest(iBC,11))
            y(:,1) = ymapped(:,2)
            y(:,2) = ymapped(:,3)
            y(:,3) = ymapped(:,4)
	    y(:,4) = ymapped(:,1)
	    y(:,5) = ymapped(:,5)
         endwhere
      else                      ! not interpolating pressure at inflow
         where (btest(iBC,11))
            y(:,1) = ymapped(:,2)
            y(:,2) = ymapped(:,3)
            y(:,3) = ymapped(:,4)
	    y(:,5) = ymapped(:,5)
         endwhere
      endif

!     
!     debugging variables
!     
!      if(iter.eq.nitr) then
!         write(555,556)lstep+1,deltaint,label,nfint
!         write(554,556)lstep+1,yplusi(2),ypluso(2),factu,factt,gamt
!         call flush(554)
!         call flush(555)
!      endif
! 556  format(i6,5(2x,e14.7))
!
!.... return
!

      return
      end
