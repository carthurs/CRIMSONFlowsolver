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
      function cname (i)

      logical beg
      CHARACTER*5 cname,cc

      ic0 = ICHAR("0")
      cc = " "
      ii = i

      i0 = mod (ii,10)
      ii = (ii - i0) / 10
      i1 = mod (ii,10)
      ii = (ii - i1) / 10
      i2 = mod (ii,10)
      ii = (ii - i2) / 10
      i3 = mod (ii,10)

      beg = .false.

      IF (i3 .ne. 0) then
        beg = .true.
        cc  = CHAR(ic0 + i3)
      ENDIF
      IF (i2 .ne. 0 .or. beg) then
        beg = .true.
        cc = TRIM(cc)//CHAR(ic0 + i2)
      ENDIF
      IF (i1 .ne. 0 .or. beg) then
        beg = .true.
        cc = TRIM(cc)//CHAR(ic0 + i1)
      ENDIF

      cc = TRIM(cc)//CHAR(ic0 + i0)
      cname = "." // cc

      return
      end


      function cname2 (i)

      logical      beg
      character*10 cname2,cc
      integer      il(0:8)

      ic0 = ICHAR("0")
      cc = " "
      ii = i

      il(0) = mod(ii,10)
      do k = 1,8
        ii = (ii - il(k-1)) / 10
        il(k) = mod (ii,10)
      enddo

      beg = .false.

      do k = 8,1,-1
        if (il(k) .ne. 0 .or. beg) then
          beg = .true.
          cc  = TRIM(cc) // CHAR(ic0 + il(k))
        endif
      enddo

      cc = TRIM(cc)//CHAR(ic0 + il(0))
      cname2 = "." // cc

      return
      end
