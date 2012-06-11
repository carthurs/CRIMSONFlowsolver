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
        subroutine restar (code,  q, ac)
c
c----------------------------------------------------------------------
c  This routine is the restart option.
c
c input:
c  code                 : restart option on the primitive variables
c                          eq. 'in  ', read from [restar.inp]
c                          eq. 'out ', write to  [restar.out]
c
c input or output:
c  q     (nshg,ndof)   : the variables to be read/written
c
c----------------------------------------------------------------------
c
        use readarrays          ! used to access qold, acold
        include "common.h"
        include "mpif.h"
c
        character*4 code
        character*8 mach2
        character*20 fname1,  fmt1
        character*5  cname

c
        dimension q(nshg,ndof),ac(nshg,ndof)
c arrays in the following 1 line are now dimensioned in readnblk
c        real*8    qold(nshg,ndof),acold(nshg,ndof)
c
c.... -------------------------->  'in  '  <---------------------------
c
        if (code .eq. 'in  ') then
c
c incompressible orders velocity, pressure, temperature unlike compressible
c which is what we have our files set up for
c
          q(:,1:3)=qold(:,2:4)
          q(:,4)=qold(:,1)
          if(ndof.gt.4)  q(:,5:ndof)=qold(:,5:ndof)

          ac(:,1:3)=acold(:,2:4)
          ac(:,4)=acold(:,1)
          if(ndof.gt.4)  ac(:,5:ndof)=acold(:,5:ndof)

          deallocate(qold)
          deallocate(acold)
          return
        endif
c
c.... -------------------------->  'out '  <---------------------------
c
        if (code .eq. 'out ') then

c$$$        ttim(75) = ttim(75) - tmr()
           allocate( qold(nshg,ndof) )
           allocate( acold(nshg,ndof) )
           acold=0
c           itmp = 1
c           if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
c           write (fmt1,"('(''restart.'',i',i1,',1x)')") itmp
c           write (fname1,fmt1) lstep
c
c           fname1 = trim(fname1) // cname(myrank+1)
           
c           open (unit=irstou, file=fname1, status='unknown',
c     &                                    form='unformatted', err=996)

c           write (irstou) machin, nshg, lstep
c incompressible orders velocity, pressure, temperature unlike compressible
c which is what we have our files set up for
c
           qold(:,2:4) = q(:,1:3)
           qold(:,1)   = q(:,4)
           if(ndof.gt.4) qold(:,5:ndof)   = q(:,5:ndof)
c 
           acold(:,2:4) = ac(:,1:3)
           acold(:,1)   = ac(:,4)
           if(ndof.gt.4) acold(:,5:ndof)   = ac(:,5:ndof)
c
           iqoldsiz=nshg*ndof
           call write_restart(myrank, lstep, nshg, ndof, 
     &          qold, acold)
c           write (irstou) qold
c           write (irstou) acold  ! note that this is dYdt (reordered)
c           close (irstou)

c$$$      open(unit=79,file="fort.79")
c$$$      do i=1,nshg
c$$$         write(79,78)(qold(i,j),j=1,5)
c$$$      enddo
c$$$      do i=1,nshg
c$$$         write(79,78)(acold(i,j),j=1,5)
c$$$      enddo
c$$$ 78   format(5(2x,e12.5))
c$$$      close(79)

           if (myrank.eq.master) then 
              open(unit=72,file='numstart.dat',status='old')
              write(72,*) lstep
              close(72)

           endif
           deallocate(qold)
           deallocate(acold)
c$$$          ttim(75) = ttim(75) + tmr()
           return
        endif
c
c.... ---------------------->  Error Handling  <-----------------------
c
c.... Error handling
c
        call error ('restar  ',code//'    ',0)
c
c.... file error handling
c
995     call error ('restar  ','opening ', irstin)
996     call error ('restar  ','opening ', irstou)
c
c.... end
c
        end
