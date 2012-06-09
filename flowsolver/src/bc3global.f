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
	subroutine bc3global (globMas, iBC)  
c
c----------------------------------------------------------------------
c
c This routine satisfies the BC of LHS mass matrix for a single 
c element.
c
c input:
c  iBC   (nshg) 	: boundary condition code
c  BC    (nshg,11)     : Dirichlet BC constraint parameters
c  ien   (npro,nshl)	: ien array for this element
c  EGmass(npro,nedof,nedof) : element consistent mass matrix before BC
c
c output:
c  EGmass(npro,nedof,nedof): LHS mass matrix after BC is satisfied
c
c----------------------------------------------------------------------
c
        include "common.h"
c
	dimension iBC(nshg),
     &            globMas(4*nshg,4*nshg)    


	do in=1,nshg
	   i0 = (in-1)*4
c
c.... pressure
c
	   if ( btest(iBC(in),2) ) then
	      globMas(i0+1,:) = zero
	      globMas(:,i0+1) = zero
	      globMas(i0+1,i0+1) = one
	   endif
c       
c....   velocities
c       
c       
c....   x1-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 1 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	   endif
c       
c....   x2-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 2 ) then
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one

	   endif
c       
c....   x1-velocity and x2-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 3 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	   endif
c       
c....   x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 4 ) then
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
c       
c....   x1-velocity and x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 5 ) then
	      globMas(i0+2,:) = zero
	      globMas(:,i0+2) = zero
	      globMas(i0+2,i0+2) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
c       
c....   x2-velocity and x3-velocity
c       
	   if ( ibits(iBC(in),3,3) .eq. 6 ) then
	      globMas(i0+3,:) = zero
	      globMas(:,i0+3) = zero
	      globMas(i0+3,i0+3) = one
	      globMas(i0+4,:) = zero
	      globMas(:,i0+4) = zero
	      globMas(i0+4,i0+4) = one

	   endif
c       
c....   x1-velocity, x2-velocity, and x3-velocity
c       
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
	

c       
c....   return
c       
	return
	end
