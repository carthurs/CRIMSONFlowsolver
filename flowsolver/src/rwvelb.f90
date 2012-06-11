        subroutine rwvelb (code,  q, ifail)
!
        use quadfilt
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
        !include "auxmpi.h"
!
        character*4 code
        character*5  cname
        character*8  mach2
        character*20 fname1,  fmt1
        character*20 fname2
        character*60 syscmd
!
        dimension q(nfath,nflow)
       logical exlog
!
!.... -------------------------->  'in  '  <---------------------------
!


       if (code .eq. 'in  ') then

          numNden=zero           ! in case the read fails 

          ifail=1
!
!.... open file
!
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
             open (unit=irstin, file=fname1, status='old', &
                  form='unformatted', err=877)
             
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
!
!.... -------------------------->  'out '  <---------------------------
!
      if (code .eq. 'out ') then

         itmp = 1
         if (lstep .gt. 0) itmp = int(log10(float(lstep)))+1
         write (fmt1,"('(''bar.'',i',i1,',1x)')") itmp
         write (fname1,fmt1) lstep
          fname1 = trim(fname1) // cname(myrank+1)
!     
         open (unit=irstou, file=fname1, status='unknown', &
              form='unformatted', err=996)
              
         write (irstou) machin, nshg, lstep
         if((itwmod.gt.0) .or. (irscale.ge.0)) then             
            write (irstou) q
         endif
         if((nsonmax.eq.1) .and. (iLES.gt.0)) then
            write (irstou) numNden
         endif
         close (irstou)

         call MPI_BARRIER(INEWCOMM,ierr)
!
! update links of "latest"
!
         fname2='bar.latest'
         fname2 = trim(fname2) // cname(myrank+1)
!         syscmd = 'ln -sf '//trim(fname1)// ' ' //fname2
!         write(*,*) syscmd
!         call system(syscmd)

         return
      endif
!
!.... ---------------------->  Error Handling  <-----------------------
!
!.... Error handling
!
        call error ('velb  ',code//'    ',0)
!
!.... file error handling
!
995     call error ('velb  ','opening ', irstin)
996     call error ('velb  ','opening ', irstou)
!
!.... end
!
        end
