        subroutine genBC1 (BCtmp,  iBC,  BC)
!
!----------------------------------------------------------------------
!
! This subroutine adjusts the boundary conditions to accommodate for 
! the velocity constraints in the non-axes directions. It copies the
! reduced constraint parameters in BC.
!
! input:
!  BCtmp (nshg,6+5*I3nsd) : input BC parameters (density, temperature,
!                             pressure, (nsd-1)(nsd+1) velocity params,
!                             upto 4 scalar params)
!  iBC   (nshg)           : boundary condition code
!
! output:
!  BC    (nshg,ndofBC)    : the constraint eq's parameters
!
!
! Farzin Shakib, Winter 1986.
! Zdenek Johan,  Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        
!
        dimension BCtmp(nshg,ndof+7),    iBC(nshg), &
                  BC(nshg,ndofBC),tmpbc(4)
!
        dimension tmp(nshg)
!
!.... scalars
!
        do isclr=1,nsclr
           where (btest(iBC,5+isclr)) BC(:,6+isclr) = BCtmp(:,12+isclr)
        enddo
!
!.... set up the thermodynamic properties
!
        where (btest(iBC,0)) BC(:,1) = BCtmp(:,1) ! density
        where (btest(iBC,1)) BC(:,2) = BCtmp(:,2) ! temperature
        where (btest(iBC,2)) BC(:,1) = BCtmp(:,3) ! pressure
!
!.... if the velocity in the x1-direction is specified
!
        where (ibits(iBC,3,3) .eq. 1)
          tmp     = BCtmp(:,4)**2 + BCtmp(:,5)**2 + BCtmp(:,6)**2
          BC(:,3) = tmp * BCtmp(:,7) / BCtmp(:,4)
          BC(:,4) =       BCtmp(:,5) / BCtmp(:,4)
          BC(:,5) =       BCtmp(:,6) / BCtmp(:,4)
        endwhere
!
!.... if the velocity in the x2-direction is specified
!
        where (ibits(iBC,3,3) .eq. 2)
          tmp     = BCtmp(:,4)**2 + BCtmp(:,5)**2 + BCtmp(:,6)**2
          BC(:,3) = tmp * BCtmp(:,7) / BCtmp(:,5)
          BC(:,4) =       BCtmp(:,4) / BCtmp(:,5)
          BC(:,5) =       BCtmp(:,6) / BCtmp(:,5)
        endwhere
!
!.... if the two velocities are specified (x1 & x2-direction)
!
!
!  Protect against user flipping the order of x1 and x2 in 
!  the vector 1 and vector 2.  Without this it will blow up.
!
        do i=1,nshg
          if(ibits(iBC(i),3,3) .eq. 3 .and.  &
             (BCtmp(i,4).eq.0 .or. BCtmp(i,9).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,4:7)
              BCtmp(i,4:7)=BCtmp(i,8:11)
              BCtmp(i,8:11)=tmpbc(1:4)
          endif
        enddo
        where (ibits(iBC,3,3) .eq. 3)
          tmp         = sqrt (BCtmp(:, 4)**2 + BCtmp(:, 5)**2 &
                                             + BCtmp(:, 6)**2)
          BCtmp(:, 4) = BCtmp(:, 4) / tmp
          BCtmp(:, 5) = BCtmp(:, 5) / tmp
          BCtmp(:, 6) = BCtmp(:, 6) / tmp
          BCtmp(:, 7) = BCtmp(:, 7) * tmp
!
          tmp         = sqrt (BCtmp(:, 8)**2 + BCtmp(:, 9)**2 &
                                             + BCtmp(:,10)**2)
          BCtmp(:, 8) = BCtmp(:, 8) / tmp
          BCtmp(:, 9) = BCtmp(:, 9) / tmp
          BCtmp(:,10) = BCtmp(:,10) / tmp
          BCtmp(:,11) = BCtmp(:,11) * tmp
!
          BCtmp(:, 4) = BCtmp(:, 9) * BCtmp(:, 4) &
                      - BCtmp(:, 5) * BCtmp(:, 8)
          BCtmp(:, 6) = BCtmp(:, 9) * BCtmp(:, 6) &
                      - BCtmp(:, 5) * BCtmp(:,10)
          BCtmp(:, 7) = BCtmp(:, 9) * BCtmp(:, 7) &
                      - BCtmp(:, 5) * BCtmp(:,11)
          BC(:,3)     = BCtmp(:, 7) / BCtmp(:, 4)
          BC(:,4)     = BCtmp(:, 6) / BCtmp(:, 4)
!
          BCtmp(:, 9) = BCtmp(:, 4) * BCtmp(:, 9) 
          BCtmp(:,10) = BCtmp(:, 4) * BCtmp(:,10) &
                      - BCtmp(:, 8) * BCtmp(:, 6)
          BCtmp(:,11) = BCtmp(:, 4) * BCtmp(:,11) &
                      - BCtmp(:, 8) * BCtmp(:, 7)
          BC(:,5)     = BCtmp(:,11) / BCtmp(:, 9)
          BC(:,6)     = BCtmp(:,10) / BCtmp(:, 9)
        endwhere
!
!.... if the velocity in the x3-direction is specified
!
        if (nsd .eq. 3) then
        where (ibits(iBC,3,3) .eq. 4)
          tmp     = BCtmp(:,4)**2 + BCtmp(:,5)**2 + BCtmp(:,6)**2
          BC(:,3) = tmp * BCtmp(:,7) / BCtmp(:,6)
          BC(:,4) =       BCtmp(:,4) / BCtmp(:,6)
          BC(:,5) =       BCtmp(:,5) / BCtmp(:,6)
        endwhere
        endif
!
!.... if two velocities are specified (x1 & x3-direction)
!
        if (nsd .eq. 3) then
!
!  Protect against user flipping the order of x1 and x3 in 
!  the vector 1 and vector 2.  Without this it will blow up.
!
        do i=1,nshg
          if(ibits(iBC(i),3,3) .eq. 5 .and. &
             (BCtmp(i,4).eq.0 .or. BCtmp(i,10).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,4:7)
              BCtmp(i,4:7)=BCtmp(i,8:11)
              BCtmp(i,8:11)=tmpbc(1:4)
           endif
        enddo
        where (ibits(iBC,3,3) .eq. 5)
          tmp         = sqrt (BCtmp(:, 4)**2 + BCtmp(:, 5)**2 &
                                             + BCtmp(:, 6)**2)
          BCtmp(:, 4) = BCtmp(:, 4) / tmp
          BCtmp(:, 5) = BCtmp(:, 5) / tmp
          BCtmp(:, 6) = BCtmp(:, 6) / tmp
          BCtmp(:, 7) = BCtmp(:, 7) * tmp
!
          tmp         = sqrt (BCtmp(:, 8)**2 + BCtmp(:, 9)**2 &
                                             + BCtmp(:,10)**2)
          BCtmp(:, 8) = BCtmp(:, 8) / tmp
          BCtmp(:, 9) = BCtmp(:, 9) / tmp
          BCtmp(:,10) = BCtmp(:,10) / tmp
          BCtmp(:,11) = BCtmp(:,11) * tmp
!
          BCtmp(:, 4) = BCtmp(:,10) * BCtmp(:, 4) &
                      - BCtmp(:, 6) * BCtmp(:, 8)
          BCtmp(:, 5) = BCtmp(:,10) * BCtmp(:, 5) &
                      - BCtmp(:, 6) * BCtmp(:, 9)
          BCtmp(:, 7) = BCtmp(:,10) * BCtmp(:, 7) &
                      - BCtmp(:, 6) * BCtmp(:,11)
          BC(:,3)     = BCtmp(:, 7) / BCtmp(:, 4)
          BC(:,4)     = BCtmp(:, 5) / BCtmp(:, 4)
!
          BCtmp(:, 9) = BCtmp(:, 4) * BCtmp(:, 9) &
                      - BCtmp(:, 8) * BCtmp(:, 5)
          BCtmp(:,10) = BCtmp(:, 4) * BCtmp(:,10)
          BCtmp(:,11) = BCtmp(:, 4) * BCtmp(:,11) &
                      - BCtmp(:, 8) * BCtmp(:, 7)
          BC(:,5)     = BCtmp(:,11) / BCtmp(:,10)
          BC(:,6)     = BCtmp(:, 9) / BCtmp(:,10)
        endwhere
        endif
!
!.... if two velocities are specified (x2 & x3-direction)
!
        if (nsd .eq. 3) then
!
!  Protect against user flipping the order of x2 and x3 in 
!  the vector 1 and vector 2.  Without this it will blow up.
!
        do i=1,nshg
          if(ibits(iBC(i),3,3) .eq. 6 .and. ( &
             BCtmp(i,5).eq.0 .or. BCtmp(i,10).eq.0)) then !flip them
              tmpbc(1:4)=BCtmp(i,4:7)
              BCtmp(i,4:7)=BCtmp(i,8:11)
              BCtmp(i,8:11)=tmpbc(1:4)
           endif
        enddo
        where (ibits(iBC,3,3) .eq. 6)
          tmp         = sqrt (BCtmp(:, 4)**2 + BCtmp(:, 5)**2 &
                                             + BCtmp(:, 6)**2)
          BCtmp(:, 4) = BCtmp(:, 4) / tmp
          BCtmp(:, 5) = BCtmp(:, 5) / tmp
          BCtmp(:, 6) = BCtmp(:, 6) / tmp
          BCtmp(:, 7) = BCtmp(:, 7) * tmp
!
          tmp         = sqrt (BCtmp(:, 8)**2 + BCtmp(:, 9)**2 &
                                             + BCtmp(:,10)**2)
          BCtmp(:, 8) = BCtmp(:, 8) / tmp
          BCtmp(:, 9) = BCtmp(:, 9) / tmp
          BCtmp(:,10) = BCtmp(:,10) / tmp
          BCtmp(:,11) = BCtmp(:,11) * tmp
!
          BCtmp(:, 4) = BCtmp(:,10) * BCtmp(:, 4) &
                      - BCtmp(:, 6) * BCtmp(:, 8)
          BCtmp(:, 5) = BCtmp(:,10) * BCtmp(:, 5) &
                      - BCtmp(:, 6) * BCtmp(:, 9)
          BCtmp(:, 7) = BCtmp(:,10) * BCtmp(:, 7) &
                      - BCtmp(:, 6) * BCtmp(:,11)
          BC(:,3)     = BCtmp(:, 7) / BCtmp(:, 5)
          BC(:,4)     = BCtmp(:, 4) / BCtmp(:, 5)
!
          BCtmp(:, 8) = BCtmp(:, 5) * BCtmp(:, 8) &
                      - BCtmp(:, 9) * BCtmp(:, 4) 
          BCtmp(:,10) = BCtmp(:, 5) * BCtmp(:,10)
          BCtmp(:,11) = BCtmp(:, 5) * BCtmp(:,11) &
                      - BCtmp(:, 9) * BCtmp(:, 7)
          BC(:,5)     = BCtmp(:,11) / BCtmp(:,10)
          BC(:,6)     = BCtmp(:, 8) / BCtmp(:,10)
        endwhere
        endif
!
!.... if all velocities are specified
!
        if (nsd .eq. 3) then
        where (ibits(iBC,3,3) .eq. 7)
          BC(:,3) = BCtmp(:,7) * BCtmp(:,4)
          BC(:,4) = BCtmp(:,7) * BCtmp(:,5)
          BC(:,5) = BCtmp(:,7) * BCtmp(:,6)
        endwhere
        endif
!
!.... end
!
        return
        end
