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
        subroutine gensvb (ientmp, iBCBtmp, BCBtmp, mattmp,
     &                     ienb,   iBCB,    BCB,    materb)
c
c----------------------------------------------------------------------
c
c  This routine saves the boundary element block.
c
c input:
c  ientmp (npro,nshl)           : boundary nodal connectivity
c  iBCtmp (npro,ndiBCB)         : boundary condition codes
c  BCBtmp (npro,nshlb,ndBCB)    : boundary condition values
c  mattmp (npro)                : material type flag
c
c output:
c  ienb   (npro,nshl)           : boundary nodal connectivity
c  iBCB   (npro,ndiBCB)         : boundary condition codes
c  BCB    (npro,nshlb,ndBCB)    : boundary condition values
c  materb (npro)                : material type flag
c
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension   ientmp(npro,nshl),
     &              iBCBtmp(npro,ndiBCB),    BCBtmp(npro,ndBCB)

        dimension   mattmp(npro),           ienb(npro,nshl),
     &              iBCB(npro,ndiBCB),      BCB(npro,nshlb,ndBCB),
     &              materb(npro)
c
c.... generate the boundary element mapping
c
        do i = 1, nshl
          ienb(:,i) = ientmp(:,i)
        enddo
c
c.... save the boundary element data
c
        iBCB   = iBCBtmp
        do i = 1, nenbl ! This is NOT NSHLB as we are just copying the
                        ! piecewise constant data given by NSpre and
                        ! higher order coefficients must be zero
           do j = 1, ndBCB
              BCB(:,i,j)   = BCBtmp(:,j)
           end do
        end do
        do i = nenbl+1, nshlb
           do j = 1, ndBCB
              BCB(:,i,j)   = zero
           end do
        end do

        materb = mattmp
c
c.... return
c
        return
        end
