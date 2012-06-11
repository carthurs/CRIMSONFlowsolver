      subroutine eqn_plane(x, iBC)
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
       use pointer_data
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
       include "mpif.h"
       !include "auxmpi.h"
!
        dimension x(numnp,nsd),       nrin(numnp), &
        	  xA(nsd),  xB(nsd),  xC(nsd), xD(nsd), fourth(nsd), &
      		  erreur(nshl), iBC(numnp), xtmp(nsd)

        integer  temp, tempb, etmp

        integer, allocatable :: ien(:,:)
        integer, allocatable :: ienb(:,:)
		
       
!
! .... find all nodes on inlet plane
!
        j = 0
        do i=1,numnp
	  if(btest(iBC(i),11)) then
	    j = j + 1
	    nen1(j) = i
	  endif
	enddo
	npin = j
	
!
! .... find one element and its vertices on inlet plane
!
        do iblk=1,nelblb
          iel=lcblkb(1,iblk)
	  nenl=lcblkb(5,iblk)
          nenlb=lcblkb(6,iblk)
	  nshl=lcblkb(9,iblk)
          nshlb=lcblkb(10,iblk)
	  npro=lcblkb(1,iblk+1)-iel
	 
	
	  allocate (ienb(npro,nshl))
	
	  ienb(:,:)=mienb(iblk)%p(:,:)
	  do i=1,npro
	    inum=i
	    tempb=0	   
	    do k=1,nshlb
	      do j=1,npin
	        if (ienb(i,k).eq.nen1(j)) then
	          tempb=tempb+1
	        endif
	      enddo   
	    enddo
	    if (tempb.eq.3) then
	      xA(:) = x(ienb(inum,1),:)
	      xB(:) = x(ienb(inum,2),:)
	      xC(:) = x(ienb(inum,3),:)
	      xD(:) = x(ienb(inum,4),:)
	      deallocate (ienb)
	      goto 92
	    endif
	  enddo
	  deallocate (ienb) 
        enddo
 92     continue
       
!
! .... find normal to inlet plane
!
	xnrml = (xB(2) - xA(2)) * (xC(3) - xA(3)) &
              - (xB(3) - xA(3)) * (xC(2) - xA(2))
     
     	ynrml = (xB(3) - xA(3)) * (xC(1) - xA(1)) &
              - (xB(1) - xA(1)) * (xC(3) - xA(3))
     
     	znrml = (xB(1) - xA(1)) * (xC(2) - xA(2)) &
              - (xB(2) - xA(2)) * (xC(1) - xA(1))
     
     	tmp = xnrml*xnrml + ynrml*ynrml + znrml*znrml
	tmp = sqrt(tmp)
	
	xnrml = xnrml/tmp
     	ynrml = ynrml/tmp
     	znrml = znrml/tmp
	
	fourth(:) = xD(:) - xA(:)
	scal = xnrml*fourth(1) + ynrml*fourth(2) + znrml*fourth(3)
	
	if (scal .lt. 0.0) then
	  xnrml = -1.0*xnrml
	  ynrml = -1.0*ynrml
	  znrml = -1.0*znrml
	endif

!	
! .... find the equation of internal plane
!      equation of plane given by:
!      (xn)x + (yn)y + (zn)z = a
!       
        aI = xnrml*xA(1)+ynrml*xA(2)+znrml*xA(3)
!	aR = xnrml*xA(1)+ynrml*xA(2)+znrml*(xA(3)+rcydist) 
	aR = aI + plandist
	if (thetag.eq.0.0) then 
	  sang = 1.0
	else
	  angle = atan2(znrml,sqrt(xnrml*xnrml+ynrml*ynrml))
	  sang = sin(angle)
	endif

!
! .... blocking elements cutting the recycle plane
!
	itmp = lcblk(1,nelblk+1)
	allocate (ien2D(itmp,4))
       	nelint=0
       	do iblk=1,nelblk
          iel=lcblk(1,iblk)
	  nenl=lcblk(5,iblk)
	  nshl=lcblk(10,iblk)
	  npro=lcblk(1,iblk+1)-iel
	 
	  allocate (ien(npro,nshl))

	  ien(:,:)=mien(iblk)%p(:,:)

	  do i=1,npro
	    inum=iel+i-1

	    temp=0	   
	    do k=1,nshl
	      erreur(k) = aR - xnrml*x(ien(i,k),1)  &
      		- ynrml*x(ien(i,k),2) - znrml*x(ien(i,k),3) 
     	    enddo
	 
	    do j=1,nshl
	      do l=j,nshl
	        if (erreur(j)*erreur(l).le. 0.0) then
 	          nelint=nelint+1
		  ien2D(nelint,:) = ien(i,:)
		  goto 95
	        endif
	      enddo    	     
	    enddo
 95	    continue
	  enddo
	  deallocate (ien)

        enddo


! 
! .... For each node of inlet plane find the corresponding element
!      and local coordinates on recycle plane
!
	allocate (xintl(nelint,nshl,nsd))
	do i = 1, nelint
              xintl(i,:,1) = x(ien2D(i,:),1)
	      xintl(i,:,2) = x(ien2D(i,:),2)
	      xintl(i,:,3) = x(ien2D(i,:),3)
        enddo
	
	if (thetag .eq. 0.0) then
	  call renum_cart(x)
	else
	  call renum_cyl(x)
	endif

        return
        end
