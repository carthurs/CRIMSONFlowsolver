      subroutine asbwmod (y,      ac,      x,      BC,   iBC, &
                          iper,   ilwork,  ifath,  velbar)
!
!----------------------------------------------------------------------
!
! This routine assembles traction BCs for a modeled wall
!
!----------------------------------------------------------------------
!
      use pointer_data
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension y(nshg,ndof),         x(numnp, nsd), &
                BC(nshg,ndofBC),      iBC(nshg), &
                iper(nshg),           ilwork(nlwork), &
                ifath(numnp),         velbar(nfath,nflow), &
                ac(nshg,ndof)
!
!.... compute and assemble the residuals corresponding to the 
!     boundary integral
!
              call settauw (y,              x, &
                   BC, &
                   ifath,                   velbar)
!
!.... enforce the new BC for SA variable
!
           isclr = 1
           if (iRANS.eq.-1) then ! S-A RANS
              call itrBCSclr (y, ac, iBC,  BC, iper, ilwork)
           endif
!
           return
           end
