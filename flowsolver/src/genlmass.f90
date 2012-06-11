        subroutine genlmass (x, shp,shgl)
!
        use pointer_data
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
!
        real*8 x(numnp,nsd)
!
        real*8 shp(MAXTOP,maxsh,MAXQPT),    &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT) 
!
        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        
!
! gmass came in via pointer_data and will 
! be available wherever it is included  (allocate it now).
!

        allocate (gmass(nshg))  
        gmass=zero
!
!.... loop over the element-blocks
!
        do iblk = 1, nelblk
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
!
!
!.... compute and assemble the residual and tangent matrix
!
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))
          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
         

          call AsImass (x,       tmpshp,      &
                        tmpshgl, mien(iblk)%p, &
                        gmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl )
!
!.... end of interior element loop
!
       enddo

      return
      end

        subroutine AsImass (x,      shp, &
                           shgl,    ien,      &
                           gmass)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the mass corresponding to the
! each node.
!
! Ken Jansen, Winter 2000.  (Fortran 90)
!----------------------------------------------------------------------
!
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        real*8 x(numnp,nsd),               &
               shp(nshl,maxsh),       shgl(nsd,nshl,ngauss), &
               gmass(nshg)

        integer ien(npro,nshl)

!
        real*8    xl(npro,nenl,nsd),    WdetJ(npro),  &
                  sgn(npro,nshl),       shapeVar(npro,nshl),           &
                  locmass(npro,nshl),   shg(npro,nshl,nsd), &
                  fmstot(npro),         temp(npro), &
                  dxidx(npro,nsd,nsd),  shdrv(npro,nsd,nshl)

        integer aa
!        
!
!
!.... gather the variables
!
!
!.... get the matrix of mode signs for the hierarchic basis functions. 
!
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localx(x,      xl,     ien,    nsd,    'gather  ')
!
!.... zero the matrices if they are being recalculated
!

        locmass=zero
        fmstot=zero

        do intp = 1, ngauss

           if (Qwt(lcsyst,intp) .eq. zero) cycle ! precaution
!
!.... get the hierarchic shape functions at this int point
!
           call getshp(shp,          shgl,      sgn,  &
                       shapeVar,        shdrv)

!
!.... --------------------->  Element Metrics  <-----------------------
!
           call e3metric( xl,         shdrv,       dxidx,   &
                          shg,        WdetJ)

!
!  get this quad points contribution to the integral of the square of  the 
!  shape function
!
           do aa = 1,nshl
              locmass(:,aa)= locmass(:,aa)  &
                   + shapeVar(:,aa)*shapeVar(:,aa)*WdetJ
           enddo
!
! also accumulate this quad points contribution to the integral of the element
! volume (integral Na^2 d Omega)
! 
           fmstot= fmstot + WdetJ ! intregral  d Omega
!
!.... end of integration loop
!
        enddo
!
!.... lumped mass if needed   Note that the locmass factors accumulated
!     over integration points and weighted with WdetJ already.
!

!.... scale the LHS matrix contribution with special lumping weighting
!
!  The first term we collect is the trace of integral Na^2 d Omega
!
        temp = zero
        do aa = 1, nshl
           temp = temp + locmass(:,aa) !reusing temp to save memory
        enddo

!
! scale the diagonal so that the trace will still yield Omega^e (the volume
! of the element)
!
        do aa = 1, nshl
           locmass(:,aa) = locmass(:,aa) * fmstot / temp
        enddo
!
!.... assemble the residual
!
        call local (gmass,    locmass,     ien,    1,  'scatter ')

!
!.... end
!
        return
        end

      subroutine lmassadd ( ac,       res, &
                            rowp,     colm,     &
                            lhsK,     gmass)
!     
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      real*8 ac(nshg,ndof), res(nshg,4), tmp,tmp1
      real*8 lhsK(9,nnz_tot), gmass(nshg), rho(nshg)
      integer rowp(nnz*nshg),  colm(nshg+1)
      integer	n,	k
!
      integer sparseloc
!
!
      rho=datmat(1,1,1)  ! needs to be generalized for VOF or level set
      tmp1=flmpl*almi
      if((flmpl.ne.0).and.(lhs.eq.1)) then
!
!.... Add lmass to diag of lhsK
!
         do n = 1, nshg
	    k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n ) &
             + colm(n)-1
            tmp=gmass(n)*tmp1*rho(n)
	    lhsK(1,k) = lhsK(1,k) + tmp
	    lhsK(5,k) = lhsK(5,k) + tmp
	    lhsK(9,k) = lhsK(9,k) + tmp
         enddo
      endif

      tmp1=flmpr

      if(flmpr.ne.0) then
         rho=rho*gmass*tmp1  ! reuse rho
         res(:,1)=res(:,1)-ac(:,1)*rho(:)
         res(:,2)=res(:,2)-ac(:,2)*rho(:)
         res(:,3)=res(:,3)-ac(:,3)*rho(:)
      endif
     
      return
      end

      subroutine lmassaddSclr ( ac,       res, &
                                rowp,     colm,     &
                                lhsS,     gmass)
!     
       use phcommonvars
       IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      real*8 ac(nshg),       res(nshg), tmp, tmp1
      real*8 lhsS(nnz_tot), gmass(nshg), rho(nshg)
      integer rowp(nnz*nshg),  colm(nshg+1)
      integer	n,	k
!
      integer sparseloc
!
!
      rho=datmat(1,1,1)  ! needs to be generalized for VOF or level set
      tmp1=flmpl*almi
      if((flmpl.ne.0).and.(lhs.eq.1)) then
!
!.... Add lmass to diag of lhsK
!
         do n = 1, nshg
	    k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n ) &
             + colm(n)-1
            tmp=gmass(n)*tmp1*rho(n)
	    lhsS(k) = lhsS(k) + tmp
         enddo
      endif

      tmp1=flmpr
      if(flmpr.ne.0) then
         rho=rho*gmass*tmp1  ! reuse rho
         res(:)=res(:)-ac(:)*rho(:)
      endif
     
      return
      end
