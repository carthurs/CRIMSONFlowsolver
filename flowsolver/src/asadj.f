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
        subroutine Asadj (   row_fill_list,
     &                    iens,        adjcnt   )
c
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        integer row_fill_list(nshg,6*nnz),
     &          ien(npro,nshl),
     &          adjcnt(nshg), ndlist(nshl)

	integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
	ien=abs(iens)

        do i=1,npro
           do j=1,nshl
              ndlist(j)=ien(i,j)
           enddo
           do j=1,nshl
              jnd=ndlist(j)
              jlngth=adjcnt(jnd) ! current length of j's list
              do k=1,nshl 
                 knd=ndlist(k)
                 ibroke=zero
                 do l= 1,jlngth
                    if(row_fill_list(jnd,l).eq. knd) then
                       ibroke=1
                       exit
                    endif
                 enddo
                 
c
c  to get here k was not in  j's list so add it
c
                 if(ibroke.eq.0) then
                    jlngth=jlngth+1 ! lenthen list
                    if(jlngth.gt.6*nnz) then
                       write(*,*) 'increase overflow factor in genadj'
                       stop
                    endif
                    row_fill_list(jnd,jlngth)=knd ! add unique entry to list
                 endif
              enddo ! finished checking all the k's for this j
              adjcnt(jnd)=jlngth  ! update the counter
           enddo                  ! done with j's
        enddo                   ! done with elements in this block
c
c
c.... end
c
        return
        end
