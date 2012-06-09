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
      subroutine rotabc (global, iBC, code)
c---------------------------------------------------------------------
c 
c This subroutine is responsible for rotating 
c the residual and solution vectors for axisymmetric BC's.
c
c input:   
c     global(nshg,n): global vector to be rotated.
c     code:            = 'in' for rotating with the residual
c                      = 'out' for rotating the solution 
c
c  note that the cos and sin of the rotation angles are preprocessed and
c  stored in acs(1 and 2) respectively.
c
c---------------------------------------------------------------------
c
      use specialBC  ! gives us acs, contains (:,1)=cos(theta) (:,2)=sin(theta)
      include "common.h"
 
      dimension global(nshg,2),             iBC(nshg),
     &          tmp(nshg)
 
      character*3 code

      if (code .eq. 'in ') then
         where( btest(iBC,10))
            tmp         =  global(:,1)*acs(:,1) - global(:,2)*acs(:,2)
            global(:,2) =  global(:,1)*acs(:,2) + global(:,2)*acs(:,1)
            global(:,1) = tmp
         endwhere
      else  if (code .eq. 'out') then
         where( btest(iBC,10))
            tmp         =  global(:,1)*acs(:,1) + global(:,2)*acs(:,2)
            global(:,2) = -global(:,1)*acs(:,2) + global(:,2)*acs(:,1)
            global(:,1) = tmp
         endwhere
      else 
         call error ('rotabc  ','code    ',0)
      endif

      return
      end
