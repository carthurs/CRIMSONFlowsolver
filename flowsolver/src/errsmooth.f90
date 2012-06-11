      subroutine errsmooth(rerr,   x,     iper,   ilwork,  &
                           shp,    shgl,  iBC)
!
        use pointer_data
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
!
        dimension shp(MAXTOP,maxsh,MAXQPT),   &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                  shpb(MAXTOP,maxsh,MAXQPT), &
                  shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!
        dimension rerrsm(nshg, 10), rerr(nshg,10), rmass(nshg)
!
        dimension ilwork(nlwork), iBC(nshg), iper(nshg)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

!
! loop over element blocks for the global reconstruction
! of the smoothed error and lumped mass matrix, rmass
!
        rerrsm = zero
        rmass = zero
        
        do iblk = 1, nelblk
!
!.... set up the parameters
!
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk)   ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel
          ngauss = nint(lcsyst)
!
!.... compute and assemble diffusive flux vector residual, qres,
!     and lumped mass matrix, rmass

          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call smooth (rerr,                x,                        &
                     tmpshp,               &
                     tmpshgl, &
                     mien(iblk)%p, &
                     rerrsm,                    &
                     rmass)

          deallocate ( tmpshp )
          deallocate ( tmpshgl ) 
       enddo
!
       if (numpe > 1) then
          call commu (rerrsm , ilwork,  10   , 'in ')
          call commu (rmass  , ilwork,  1    , 'in ')
       endif       
!
!.... take care of periodic boundary conditions
!
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            rerrsm(i,:) = rerrsm(i,:) + rerrsm(j,:)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(j) = rmass(i)
            rerrsm(j,:) = rerrsm(i,:)
          endif
        enddo
!
!.... invert the diagonal mass matrix and find q
!
        rmass = one/rmass
       
       do i=1, 10
          rerrsm(:,i) = rmass*rerrsm(:,i)
       enddo
       if(numpe > 1) then
          call commu (rerrsm, ilwork, 10, 'out')    
       endif
!
!      copy the smoothed error overwriting the original error.
!

       rerr = rerrsm 

       return
       end

        subroutine smooth (rerr,       x,       shp, &
                           shgl,       ien,           &
                           rerrsm,     rmass    )
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
! interior elements for the global reconstruction of the diffusive
! flux vector.
!
! input:
!     y     (nshg,ndof)        : Y variables
!     x     (numnp,nsd)         : nodal coordinates
!     shp   (nshape,ngauss)     : element shape-functions
!     shgl  (nsd,nshape,ngauss) : element local shape-function gradients
!     ien   (npro)              : nodal connectivity array
!
! output:
!     qres  (nshg,nflow-1,nsd)  : residual vector for diffusive flux
!     rmass  (nshg)            : lumped mass matrix
!
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension rerr(nshg,10),               x(numnp,nsd),      &
                  shp(nshl,maxsh),   &
                  shgl(nsd,nshl,maxsh), &
                  ien(npro,nshl), &
                  rerrsm(nshg,10),    rmass(nshg)
!
!.... element level declarations
!
        dimension rerrl(npro,nshl,10),        xl(npro,nenl,nsd),          &
                  rerrsml(npro,nshl,10),       rmassl(npro,nshl)
!
        dimension sgn(npro,nshl),          shapeVar(npro,nshl), &
                  shdrv(npro,nsd,nshl),    WdetJ(npro), &
                  dxidx(npro,nsd,nsd),     shg(npro,nshl,nsd)
!
        dimension error(npro,10)
!
!.... create the matrix of mode signs for the hierarchic basis 
!     functions. 
!
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
!
!.... gather the variables
!

        call local(rerr,   rerrl,  ien,    10,   'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
!
!.... get the element residuals 
!
        rerrsml     = zero
        rmassl      = zero

!
!.... loop through the integration points
!
        
                
        do intp = 1, ngauss
        if (Qwt(lcsyst,intp) .eq. zero) cycle          ! precaution
!
!.... create a matrix of shape functions (and derivatives) for each
!     element at this quadrature point. These arrays will contain 
!     the correct signs for the hierarchic basis
!
        call getshp(shp,          shgl,      sgn,  &
                    shapeVar,     shdrv)
!
        call e3metric( xl,         shdrv,        dxidx,   &
                       shg,        WdetJ)
        error=zero
        do n = 1, nshl
           do i=1,10
              error(:,i)=error(:,i) + shapeVar(:,n) * rerrl(:,n,i)
           enddo
        enddo
        do i=1,nshl
           do j=1,10
              rerrsml(:,i,j)  = rerrsml(:,i,j)   &
                              + shapeVar(:,i)*WdetJ*error(:,j)
           enddo

           rmassl(:,i) = rmassl(:,i) + shapeVar(:,i)*WdetJ
        enddo
 
!.... end of the loop over integration points
!
      enddo
!
!.... assemble the diffusive flux residual 
!
        call local (rerrsm,   rerrsml,  ien,  10,'scatter ')
        call local (rmass,   rmassl,  ien,  1,  'scatter ')
!

      return
      end
