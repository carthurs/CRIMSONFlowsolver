	subroutine bc3global (globMas, iBC)  
!
!----------------------------------------------------------------------
!
! This routine satisfies the BC of LHS mass matrix for a single 
! element.
!
! input:
!  iBC   (nshg) 	: boundary condition code
!  BC    (nshg,11)     : Dirichlet BC constraint parameters
!  ien   (npro,nshl)	: ien array for this element
!  EGmass(npro,nedof,nedof) : element consistent mass matrix before BC
!
! output:
!  EGmass(npro,nedof,nedof): LHS mass matrix after BC is satisfied
!
!
! Farzin Shakib, Winter 1987.
! Zdenek Johan,  Spring 1990. (Modified for general divariant gas)
!----------------------------------------------------------------------
!
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
	dimension iBC(nshg), &
                  globMas(4*nshg,4*nshg)    


	do in=1,nshg
	   i0 = (in-1)*4
!
!.... pressure
!
	   if ( btest(iBC(in),2) ) then
	      globMas(i0+1,:) = zero
	      globMas(:,i0+1) = zero
	      globMas(i0+1,i0+1) = one
	   endif
!       
!....   velocities
!       
!       
!....   x1-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 1 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	   endif
!       
!....   x2-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 2 ) then
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one

	   endif
!       
!....   x1-velocity and x2-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 3 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	   endif
!       
!....   x3-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 4 ) then
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
!       
!....   x1-velocity and x3-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 5 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
!       
!....   x2-velocity and x3-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 6 ) then
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
!       
!....   x1-velocity, x2-velocity, and x3-velocity
!       
	   if ( ibits(iBC(in),3,3) .eq. 7 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one
	   endif
	enddo
	

!       
!....   return
!       
	return
	end
