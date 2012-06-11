      subroutine renum_cyl(x)
!
!----------------------------------------------------------------------
! This subroutine finds all nodes that are on the inlet plane and all
! nodes that are on the recycle plane; it also blocks elements that 
! are on recycle plane; all nodes are also stored with cylindrical 
! coordinates; find all father nodes for recycle plane (i.e. for
! theta = -theta given)
!
! input:
!  x      (numnp,nsd)           : node coordinates
!
! output:
!  xcyl   (numnp,nsd)           : node cylindrical coordinates 
!  ien2D  (npro, nshl)		: connectivity array for recycle plane
!				  assuming tethraheadral elements, i.e.
!				  triangular elements on face
!
!
!  Elaine Bohr
!  June 2002
!----------------------------------------------------------------------
!
       use spebc
!       use pointer_data
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
       include "mpif.h"
       !include "auxmpi.h"
!
        dimension x(numnp,nsd), nrin(numnp), nula(numnp), &
      		  erreur(nshl), xtmp(nsd)

        integer  temp, etmp, s


	real*8,  allocatable :: xrtmp(:)
	        

!	thetag = thetag/180.0*pi


!
! .... changing to cylindrical coordinate system for nodal point
!       
	
	xcyl(:,1) = sqrt(x(:,1)*x(:,1) + x(:,2)*x(:,2))
	
	j = 0
        do i=1,numnp
	  if ((x(i,1).eq.0).and.(x(i,2).eq.0)) then
	    j = j+1
	    nula(j) = i
	  else
	    xcyl(i,2) = atan2(x(i,2),x(i,1)) 
	  endif
	enddo
        xcyl(:,3) = x(:,3)
	

!
! .... finding the minimum and maximum angles
!

	thmin = xcyl(nen1(1),2)
        thmax = xcyl(nen1(1),2)
        do i=2,npin
	  if((x(nen1(i),1).eq.0).and.(x(nen1(i),2).eq.0)) then
	    goto 10
	  else
            if(thmin.gt.xcyl(nen1(i),2)) thmin = xcyl(nen1(i),2)
            if(thmax.lt.xcyl(nen1(i),2)) thmax = xcyl(nen1(i),2)
	  endif
 10	  continue
        enddo
	
	do i=1,j
	  xcyl(nula(i),2) = thmin
	enddo
	

!
! .... Finding nodes from inlet plane with same theta given
!
       
       j = 0
       do i=1,npin
         error1 = abs(xcyl(nen1(i),2) - thmin)
	 if(error1.le.0.00001) then
	   j = j + 1
	   nrin(j)=nen1(i)
	 endif
       enddo
       nfin = j
       allocate(xrtmp(nfin))
       do i=1,nfin
	 xrtmp(i) = xcyl(nrin(i),1)
       enddo	
       
!
! .... Ordering nrin by decreasing radius
!
	j = 0
	allocate (nrint(nfin))
	
	radinl = xcyl(nrin(1),1)
	do i=1,nfin
	  if (radinl.lt.xcyl(nrin(i),1)) radinl = xcyl(nrin(i),1)
	  rmaxtemp=xrtmp(1)
	  itmp = 1
	  do k=2,nfin
	    if (rmaxtemp.le.xrtmp(k)) then
	      rmaxtemp = xrtmp(k)
	      itmp = k
	    endif
	  enddo
	  if (radcyl.ge.xrtmp(itmp)) then
	    j = j + 1
	    nrint(j)=nrin(itmp)
	  endif
	  xrtmp(itmp) = -1
	enddo
	nfint = j
	
	deallocate(xrtmp)

!
! .... off wall coordinate for the inlet plane
!
	do i=1,npin
	  xynin(i) = (radinl - xcyl(nen1(i),1))/sang
     	enddo

!
! .... off wall coordinate for the virtual points on recycle plane
!
	
	do i=1,nfint
	  xyn(i) = (radcyl - xcyl(nrint(i),1))/sang
     	enddo

!
! .... Finding corresponding elements and local coordinates
!      on recycle plane for every arc with radius from nrint
!
       s = (thmax-thmin)*radcyl/ds
       allocate (xsinfin(nfint,s+1,nsd))
       allocate (elcnfin(nfint,s+1))
       allocate (imax(nfint))

       do jj = 1, nfint
            
          xts1 = x(nrint(jj),1) !*cos(xcyl(nrint(jj),2)+tolerence)
          xts2 = x(nrint(jj),2) !*sin(xcyl(nrint(jj),2)+tolerence)
          xts3 = x(nrint(jj),3) + plandist
	  call elem_search(xintl, xts1, xts2, xts3, &
      		           xtmp(:), etmp, 1)
     	  xsinfin(jj,1,:) = xtmp(:)
	  elcnfin(jj,1) = etmp
	  imax(jj) = (thmax-thmin)*xcyl(nrint(jj),1)/ds
     	  do i=1,imax(jj)
	    if ( xcyl(nrint(jj),1) .eq. radcyl) then
	      xts1 = (xcyl(nrint(jj),1)-tolerence)  &
      		*cos(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	      xts2 = (xcyl(nrint(jj),1)-tolerence)  &
      		*sin(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	    else
!	      reel=i*ds/xcyl(nrint(jj),1)
	      xts1 = xcyl(nrint(jj),1) &
      		*cos(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	      xts2 = xcyl(nrint(jj),1) &
       		*sin(1.0*i/imax(jj)*(thmax-thmin)+thmin)
	    endif
	    xts3 = x(nrint(jj),3) + plandist
	    call elem_search(xintl, xts1, xts2, xts3, &
      		             xtmp(:), etmp, 1)
     	    xsinfin(jj,1+i,:) = xtmp(:)
	    elcnfin(jj,1+i) = etmp
          enddo
	enddo 

        return
        end


      subroutine renum_cart(x)
!
!----------------------------------------------------------------------
! This subroutine finds all nodes that are on the inlet plane and all
! nodes that are on the recycle plane; it also blocks elements that 
! are on recycle plane; all nodes are also stored with cylindrical 
! coordinates; find all father nodes for recycle plane (i.e. for
! theta = -theta given)
!
! input:
!  x      (numnp,nsd)           : node coordinates
!
! output:
!  xcyl   (numnp,nsd)           : node cylindrical coordinates 
!  ien2D  (npro, nshl)		: connectivity array for recycle plane
!				  assuming tethraheadral elements, i.e.
!				  triangular elements on face
!
!
!  Elaine Bohr
!  July 2002
!----------------------------------------------------------------------
!
       use spebc
!       use pointer_data
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
       include "mpif.h"
       !include "auxmpi.h"
!
        dimension x(numnp,nsd), nyin(numnp), &
      		  erreur(nshl), xtmp(nsd), yrtmp(numnp)

        integer  temp, etmp, s

!
! .... Finding nodes from inlet plane with minimal z value
!
       j = 0
       zmin = x(nen1(1),3)
       zmax = x(nen1(1),3)
       do i=2,npin
	 if(x(nen1(i),3).lt.zmin) zmin = x(nen1(i),3)
	 if(x(nen1(i),3).gt.zmax) zmax = x(nen1(i),3)
       enddo

       do i=1,npin
         if (x(nen1(i),3).eq.zmin) then
	   j = j + 1
	   nyin(j) = nen1(i)
	 endif  
       enddo
       nfin = j
       do i=1,nfin
	 yrtmp(i) = x(nyin(i),2)
       enddo	
       
!
! .... Ordering nyin by increasing y
!
	j = 0
	allocate (nrint(nfin))
	
	do i=1,nfin
	  rmintemp=yrtmp(1)
	  itmp = 1
	  do k=2,nfin
	    if (rmintemp.ge.yrtmp(k)) then
	      rmintemp = yrtmp(k)
	      itmp = k
	    endif
	  enddo
	  j = j + 1
	  nrint(j)=nyin(itmp)
	  yrtmp(itmp) = 10000000
	enddo
	nfint = j

!
! .... y coordinate for the inlet plane
!
	do i=1,npin
	  xynin(i) = x(nen1(i),2)
     	enddo

!
! .... y coordinate for the virtual points on recycle plane
!
	
	do i=1,nfint
	  xyn(i) = x(nrint(i),2)
     	enddo

!
! .... Finding corresponding elements and local coordinates
!      on recycle plane for every arc with radius from nrint
!
       s = (zmax - zmin) / ds
       allocate (xsinfin(nfint,s+1,nsd))
       allocate (elcnfin(nfint,s+1))
       allocate (imax(nfint))

       do jj = 1, nfint
            
          xts1 = x(nrint(jj),1) + plandist
          xts2 = x(nrint(jj),2) 
          xts3 = x(nrint(jj),3) 
	  call elem_search(xintl, xts1, xts2, xts3, &
      		           xtmp(:), etmp, 1)
     	  xsinfin(jj,1,:) = xtmp(:)
	  elcnfin(jj,1) = etmp
	  imax(jj) = s
     	  do i=1,imax(jj)
	    xts1 = x(nrint(jj),1) + plandist
	    xts2 = x(nrint(jj),2) 
	    xts3 = x(nrint(jj),3) + i*ds
	    call elem_search(xintl, xts1, xts2, xts3, &
      		             xtmp(:), etmp, 1)
     	    xsinfin(jj,1+i,:) = xtmp(:)
	    elcnfin(jj,1+i) = etmp
          enddo
	enddo 

        return
        end
