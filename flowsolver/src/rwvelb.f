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
        subroutine rwvelb (code,  q, ifail)
c
        use quadfilt
        include "common.h"
        include "mpif.h"
        include "auxmpi.h"
c
        character*4 code
        character*5  cname
        character*8  mach2
        character*20 fname1,  fmt1
        character*20 fname2
        character*60 syscmd
c
        dimension q(nfath,nflow)
       logical exlog
c
c.... -------------------------->  'in  '  <---------------------------
c


       if (code .eq. 'in  ') then

          numNden=zero           ! in case the read fails 

          ifail=1
c
c.... open file
c
          fname1='bar.latest'
          fname1 = trim(fname1) // cname(myrank+1)
          inquire(file=fname1,exist=exlog)
          
          if(exlog)  then
          else
             open(unit=72,file='numstart.dat',status='old')
             read(72,*) irstart
             close(72)
             itmp = 1
             if (irstart .gt. 0) itmp = int(log10(float(irstart)))+1
             write (fmt1,"('(''bar.'',i',i1,',1x)')") itmp
             write (fname1,fmt1) irstart
             fname1 = trim(fname1) // cname(myrank+1)
          
             inquire(file=fname1,exist=exlog)
          endif

          write (*,*) 'Reading bar field file : ', fname1

          if(exlog) then        ! velb exists; open and use it
             open (unit=irstin, file=fname1, status='old',
     &            form='unformatted', err=877)
             
             read (irstin) mach2, nshg2, lstep2
             if((itwmod.gt.0) .or. (irscale.ge.0)) then             
                read (irstin,err=877,end=877) q
                write(*,*) "velb found and read properly"
                ifail=0         ! i.e. I didn't fail
             endif
             if((nsonmax.eq.1) .and. (iLES.gt.0)) then
                read (irstin,err=888,end=888) numNDen
                write(*,*) "numDen found and read properly",myrank+1
             endif
 888         continue
 877         continue    
             close (irstin)
          endif                 ! if bar field file exists

          return
       endif
c
c.... -------------------------->  'out '  <---------------------------
c
      if (code .eq. 'out ') then

         itmp = 1
         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
         write (fmt1,"('(''bar.'',i',i1,',1x)')") itmp
         write (fname1,fmt1) lstep
          fname1 = trim(fname1) // cname(myrank+1)
c     
         open (unit=irstou, file=fname1, status='unknown',
     &        form='unformatted', err=996)
              
         write (irstou) machin, nshg, lstep
         if((itwmod.gt.0) .or. (irscale.ge.0)) then             
            write (irstou) q
         endif
         if((nsonmax.eq.1) .and. (iLES.gt.0)) then
            write (irstou) numNden
         endif
         close (irstou)

         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
c
c update links of "latest"
c
         fname2='bar.latest'
         fname2 = trim(fname2) // cname(myrank+1)
c         syscmd = 'ln -sf '//trim(fname1)// ' ' //fname2
c         write(*,*) syscmd
c         call system(syscmd)

         return
      endif
c
c.... ---------------------->  Error Handling  <-----------------------
c
c.... Error handling
c
        call error ('velb  ',code//'    ',0)
c
c.... file error handling
c
995     call error ('velb  ','opening ', irstin)
996     call error ('velb  ','opening ', irstou)
c
c.... end
c
        end
