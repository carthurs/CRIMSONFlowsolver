!-----------------------------------------------------------------------
!
!  Natural pressure boundary condition can be calculated with p, the pressure,
!  related (in some prescribed manner) to Q, the flow rate, through the same 
!  boundary.  To do this efficiently requires us to precompute the integral 
!  of N_A over the boundary for each node A and store it in a vector of length
!  nshg (a bit wasteful since only the nodes on the boundary will be non zero
!  in this vector but it is probably slower to index it than to multiply and 
!  add the extra zeros....check later).
!
!-----------------------------------------------------------------------
      module pvsQbi

      real*8, allocatable ::  NABI(:,:)
      real*8, allocatable ::  NASC(:)
      real*8, allocatable ::  PNABI(:,:)
      real*8, allocatable ::  NANBIJ(:,:,:)
      integer, allocatable :: ndsurf(:)
      
      end module
      
      subroutine DpvsQbi
      use pvsQbi
      if(allocated(NABI)) deallocate(NABI)
      if(allocated(NASC)) deallocate(NASC)
      if(allocated(PNABI)) deallocate(PNABI)
      if(allocated(NANBIJ)) deallocate(NANBIJ)
      deallocate(ndsurf)
      return
      end

!-----------------------------------------------------------------------
!
!     Initialize:
!
!-----------------------------------------------------------------------
      subroutine initNABI( x, shpb )
      
      use     pointer_data
      use     pvsQbi
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8   x(numnp,nsd)
!
! use is like
! 
!      NABI=pvsQbi -> NABI
!
        dimension   shpb(MAXTOP,maxsh,MAXQPT)
        real*8, allocatable :: tmpshpb(:,:)
        allocate ( NABI(nshg,3) ) 
        allocate ( NASC(nshg)   )
        allocate ( ndsurf(nshg) ) 

!
!....  calculate NABI
!
      NABI=zero
      NASC=zero
      ndsurf=0
        if (Lagrange .gt. 0) then
           allocate ( PNABI(nshg,3)  )
           allocate ( NANBIJ(nshg,3,3)  )
           PNABI = zero
           NANBIJ = zero
        endif
!
!.... -------------------->   boundary elements   <--------------------
!
!.... set up parameters
!
!        intrul = intg   (2,itseq)
!        intind = intptb (intrul)
!
!.... loop over the boundary elements
!
        do iblk = 1, nelblb
!
!.... set up the parameters
!
          iel    = lcblkb(1,iblk)
          lelCat = lcblkb(2,iblk)
          lcsyst = lcblkb(3,iblk)
          iorder = lcblkb(4,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          mattyp = lcblkb(7,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel 


          if(lcsyst.eq.3) lcsyst=nenbl
!
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
          
!
!.... compute and assemble the residuals corresponding to the 
!     boundary integral
!
          allocate (tmpshpb(nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)

          call AsBNABI (                       x, &
                       tmpshpb, &
                       mienb(iblk)%p, &
                       miBCB(iblk)%p)

          call AsBNASC(                       x, &
                       tmpshpb, &
                       mienb(iblk)%p, &
                       miBCB(iblk)%p)
          
          if (Lagrange .gt. 0) then 
             call AsBPNABI(             x, &
                       tmpshpb, &
                       mienb(iblk)%p, &
                       miBCB(iblk)%p) 
          endif    

          deallocate (tmpshpb)

      enddo 

!
!     note that NABI has NOT been communicated.  It
!     is the on processor contribution to this vector.  It will used to
!     build the on processor contribution to res that will then be made
!     complete via a call to commu.  Similarly the LHS usage will create
!     the on-processor contribution to the lhsK. Same for NASC
!
      return
      end


