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


c$$$      subroutine get_a_not_hex(xc,anot)
c$$$
c$$$      include "common.h"
c$$$
c$$$      dimension xc(npro,nenl,nsd), anot(npro,nenl,nsd)
c$$$
c$$$
c$$$      do i = 1, nsd
c$$$
c$$$         anot(:,1,i) = pt125*(xc(:,1,i)+xc(:,2,i)+xc(:,3,i)+xc(:,4,i)
c$$$     &        +xc(:,5,i)+xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$         
c$$$         anot(:,2,i) = pt125*(-xc(:,1,i)+xc(:,2,i)+xc(:,3,i)-xc(:,4,i)
c$$$     &        -xc(:,5,i)+xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$         anot(:,3,i) = pt125*(-xc(:,1,i)-xc(:,2,i)+xc(:,3,i)+xc(:,4,i)
c$$$     &        -xc(:,5,i)-xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$         
c$$$         anot(:,4,i) = pt125*(-xc(:,1,i)-xc(:,2,i)-xc(:,3,i)-xc(:,4,i)
c$$$     &        +xc(:,5,i)+xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$
c$$$         anot(:,5,i) = pt125*(xc(:,1,i)-xc(:,2,i)+xc(:,3,i)-xc(:,4,i)
c$$$     &        +xc(:,5,i)-xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$         anot(:,6,i) = pt125*(xc(:,1,i)+xc(:,2,i)-xc(:,3,i)-xc(:,4,i)
c$$$     &        -xc(:,5,i)-xc(:,6,i)+xc(:,7,i)+xc(:,8,i))
c$$$
c$$$         anot(:,7,i) = pt125*(xc(:,1,i)-xc(:,2,i)-xc(:,3,i)+xc(:,4,i)
c$$$     &        -xc(:,5,i)+xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$         anot(:,8,i) = pt125*(-xc(:,1,i)+xc(:,2,i)-xc(:,3,i)+xc(:,4,i)
c$$$     &        +xc(:,5,i)-xc(:,6,i)+xc(:,7,i)-xc(:,8,i))
c$$$
c$$$      enddo
c$$$
c$$$      return
c$$$      end


      subroutine get_coeff_tet(xc,anot)

      use spebc
      include "common.h"
      

      dimension xc(nelint,nenl,nsd), anot(nelint,nenl,nsd)


      do i = 1, nsd

         anot(:,1,i) = xc(:,4,i)
         anot(:,2,i) = xc(:,1,i)-xc(:,4,i)
         anot(:,3,i) = xc(:,2,i)-xc(:,4,i)
         anot(:,4,i) = xc(:,3,i)-xc(:,4,i)
         
      enddo
      
      return
      end
      
