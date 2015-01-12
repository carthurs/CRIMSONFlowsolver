      subroutine forces (y ,  ilwork)
!
!----------------------------------------------------------------------
!
! This subroutine calculates and updates the aerodynamic forces
!
!----------------------------------------------------------------------
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"
!     
      integer ilwork(nlwork)

      real*8  y(nshg,ndof)
      real*8  Forin(5), Forout(5),spmasss
      real*8  ftots(3,0:MAXSURF),ftot(3)
!
!  START OF DISTURBANCE BLOCK
!
      dimension diste(1), distesum(1)
      ! \todo-binary I removed this, together with all other lines in the code for dke or dkesum commented with the random string sdiunbgjtnz as it didnt seem to be doing anything, and was causing memory errors with strict ifort checks
!         if(iter.eq.nitr) then  ! we have completed the last iteration
!           if (numpe > 1) then
!           call MPI_REDUCE (dke,dkesum,1,MPI_DOUBLE_PRECISION, &
!                            MPI_SUM, master, INEWCOMM,ierr)
!           endif
!           if (numpe.eq.1) &
!             write (76,1000) lstep+1, dke,dke

!           if ((myrank .eq. master).and.(numpe > 1)) then 
!              write (76,1000) lstep+1, dkesum,dkesum
! !
!              call flush(76)
! !     
!           endif
!           dke=0.0               ! we must zero it back out for next step

!        endif

!
! END OF DISTURBANCE BLOCK
!

!
!.... -------------------->  Aerodynamic Forces  <----------------------
!
!.... output the forces and the heat flux
!
      if (numpe > 1) then
         Forin = (/ Force(1), Force(2), Force(3), HFlux,  &
                    entrop /)
         call MPI_ALLREDUCE (Forin(1), Forout(1),5, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
         Force = Forout(1:3)
         HFlux = Forout(4)
         entrop= Forout(5)
      endif

      if (numpe > 1) then
         call MPI_ALLREDUCE (flxID(2,isrfIM), spmasss,1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
         call MPI_ALLREDUCE (flxID(1,isrfIM), Atots,1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
         call MPI_ALLREDUCE (flxID(3,:), Ftots(1,:),MAXSURF+1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
         call MPI_ALLREDUCE (flxID(4,:), Ftots(2,:),MAXSURF+1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
         call MPI_ALLREDUCE (flxID(5,:), Ftots(3,:),MAXSURF+1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
      else
         Ftots=flxID(3:5,:)
         Atots=flxID(1,isrfIM)
         spmasss=flxID(2,isrfIM)
      endif
      call sforce()

      ftot(1)=sum(Ftots(1,1:MAXSURF))
      ftot(2)=sum(Ftots(2,1:MAXSURF))
      ftot(3)=sum(Ftots(3,1:MAXSURF))
   
      if((spmasss.ne.0) .and. (matflg(5,1).ne.0))then ! we are doing
                                                      ! force adjustment 
         tmp=Dtgl*0.025         ! .025= 1/2/tmp1 when tmp1=20
         
         select case ( matflg(5,1) )
         case ( 1 )             ! standard linear body force
            fwall=Force(1)
            vlngth=xlngth
         case ( 2 )       
         case ( 3 )     
            fwall=Force(3)
            vlngth=zlngth
         end select

         write(222,*) spmasss
         fnew=(vel-spmasss/(ro*Atots))*tmp +  &
              fwall/(vlngth*Atots)
         if(myrank.eq.0) then
            write(880,*) datmat(1,5,1),fnew
            call flush(880)
         endif
         datmat(1,5,1)=fnew     !* (one - tmp2) + datmat(1,5,1) * tmp2
      endif

      if (myrank .eq. master) then
         write (iforce,1000) lstep, (Force(i), i=1,nsd),  &
              HFlux,spmasss,(ftot(i),i=1,nsd)
         call flush(iforce)
      endif
!
      if (pzero .eq. 1) then
!

!  also adjust the pressure here so that the mean stays at zero (prevent drift)
!
!  This is only valid if there are NO Dirichlet Boundary conditions on p!!!
!
!  (method below counts partition boundary nodes twice. When we get ntopsh
!  into the method it will be easy to fix this.  Currently it would require
!  some work and it is probably not worth it 
!
!  switched to avg of the on-processor averages to get around hierarchic
!  modes problem
!
!              pave=sum(y(1:numnp,1)) 
!              call MPI_ALLREDUCE (pave, pavet, 1,
!     &             MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
!     pavet=pavet/nshgt
         pave=sum(y(1:numnp,4)) 
         if (numpe .gt. 1) then
            xnpts=numnp
            call MPI_ALLREDUCE (pave, pavet, 1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
            call MPI_ALLREDUCE (xnpts, xnptst, 1, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
            pavet=pavet/xnptst 
         else
            pavet = pave/numnp
         endif

         y(1:numnp,4)=y(1:numnp,4)-pavet
         pave=sum(y(1:numnp,4))
!         if(myrank.eq.0) then
!            write(90+myrank,*) pavet,pave
!            call flush(900)
!         endif
!
      endif

      return
!
 1000 format(1p,i6,8e13.5)
!
      end

      subroutine sforce()

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"

      real*8  SFlux(5),SFluxg(5)
      character*20 fname1
      character*5  cname

      integer icalled 
      data icalled /0/
      save icalled

      if(icalled.eq.0) then
        icalled = 1
        do isrf = 0,MAXSURF
          if ( nsrflist(isrf).ne.0 ) then
            iunit=60+isrf
            fname1 = 'forces_s'//trim(cname(isrf))// cname(myrank+1) 
            open(unit=iunit,file=trim(fname1),status='unknown')
          endif
        enddo
      endif

      do isrf = 0,MAXSURF
        if ( nsrflist(isrf).ne.0 ) then
          SFlux = (/ flxID(1,isrf),  &                                    ! Area
                      flxID(2,isrf), &                                   ! Mass Flux
                      flxID(3,isrf), &                                   ! Force 1
                      flxID(4,isrf), &                                   ! Force 2
                      flxID(5,isrf)/)                                   ! Force 3
          if ( numpe > 1 ) then
            call MPI_ALLREDUCE (SFlux, SFluxg,5, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
            SFlux = SFluxg
          endif
          iunit=60+isrf
          write(iunit,"(i7,1p5e14.5)")lstep,SFlux
          call flush(iunit)
        endif
      enddo

      return
      end
