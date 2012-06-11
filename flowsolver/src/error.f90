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
        subroutine error (routin, variab, num)
c
c----------------------------------------------------------------------
c
c This utility routine prints out the error and stops the program.
c
c input:
c  routin       : name of the routine where the error occurred
c  variab       : an 8-character error message
c  num          : any integer number associated with the error
c
c----------------------------------------------------------------------
c
        include "common.h"
        include "mpif.h"
c
        character*8 routin, variab
c
        data ierchk /0/
c
c.... check for redundant error
c
        if (ierchk .eq. 1) stop
        ierchk = 1
c
c.... open file
c
        open (unit=ierror, file=ferror, status='unknown')
c
c.... print the error
c
        write (*,1000) title, routin, variab, num
        if (num .ne. 0) write (ierror,1000) title, routin, variab, num
        if (num .eq. 0) write (ierror,1000) title, routin, variab
c
c.... halt the process
c
        close (ierror)


        WRITE(6,'(A,G14.6)') 'Life: ',death - birth
        if (numpe > 1) then
           call MPI_ABORT(MPI_COMM_WORLD)
        endif
        
 
1000    format(' ',a80,//,
     &         ' ****** Error occurred in routine <',a8,'>',/,
     &          '  Error code :',a8,:,' : ',i8,//)
        end
