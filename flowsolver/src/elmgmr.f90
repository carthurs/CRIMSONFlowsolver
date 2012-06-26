        subroutine ElmGMR (u,         y,         ac,         &
                           x,         xdist,     xdnv,      &
                           shp,       shgl,      iBC, &
                           BC,        shpb,      shglb, &
                           res,       iper,      ilwork, &
                           rowp,      colm,      lhsK,       &
                           lhsP,      rerr)
!
!----------------------------------------------------------------------
!
! This routine computes the LHS mass matrix, the RHS residual 
! vector, and the preconditioning matrix, for use with the GMRES
! solver.
!
! Zdenek Johan, Winter 1991.      (Fortran 90)
! Chris Whiting, Winter 1998.     (Matrix EBE-GMRES)
! Alberto Figueroa, Winter 2004.  CMM-FSI
! Irene Vignon, Spring 2004.
!----------------------------------------------------------------------
!
        use pvsQbi  ! brings in NABI
        use stats   !  
        use pointer_data  ! brings in the pointers for the blocked arrays
        use local_mass
        use LagrangeMultipliers 
        use deformableWall
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndof),         ac(nshg,ndof), &
                  u(nshg,nsd), &
                  x(numnp,nsd),                &
                  xdist(numnp), &
                  xdnv(numnp,nsd), &
                  iBC(nshg),            &
                  BC(nshg,ndofBC),   &
                  res(nshg,nflow), &
                  iper(nshg)
!
        dimension shp(MAXTOP,maxsh,MAXQPT),   &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                  shpb(MAXTOP,maxsh,MAXQPT), &
                  shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!
        dimension qres(nshg,idflx),     rmass(nshg)
!
        dimension ilwork(nlwork)

        integer rowp(nshg*nnz),         colm(nshg+1)

	  real*8	lhsK(9,nnz_tot),	lhsP(4,nnz_tot)

        real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC

        real*8  rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

        real*8 spmasstot(20),  ebres(nshg)
!
!.... set up the timer
!

!AD        call timer ('Elm_Form')
!
!.... -------------------->   diffusive flux   <--------------------
!
!.... set up parameters
!
        ires   = 1

        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
!
! loop over element blocks for the global reconstruction
! of the diffusive flux vector, q, and lumped mass matrix, rmass
!
           qres = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
!     
!.... compute and assemble diffusive flux vector residual, qres,
!     and lumped mass matrix, rmass

              call AsIq (y,                x,                        &
                         shp(lcsyst,1:nshl,:),  &
                         shgl(lcsyst,:,1:nshl,:), &
                         mien(iblk)%p,     mxmudmi(iblk)%p,   &
                         qres,             rmass )
           enddo
       
!
!.... form the diffusive flux approximation
!
           call qpbc( rmass, qres, iBC, iper, ilwork )       
!
        endif 
!
!.... -------------------->   interior elements   <--------------------
!
        res    = zero
        if (stsResFlg .ne. 1) then
           flxID = zero
        endif

        if (lhs .eq. 1) then
           lhsp   = zero
           lhsk   = zero
        endif
!
!.... loop over the element-blocks
!
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iel    = lcblk(1,iblk)
          lelCat = lcblk(2,iblk)
          lcsyst = lcblk(3,iblk)
          iorder = lcblk(4,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          mattyp = lcblk(7,iblk)
          ndofl  = lcblk(8,iblk)
          nsymdl = lcblk(9,iblk)
          npro   = lcblk(1,iblk+1) - iel 
          inum   = iel + npro - 1
          ngauss = nint(lcsyst)
!
!.... allocate the element matrices
!
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
!
!..... to calculate inner product for Lagrange Multipliers
!
          if(Lagrange.gt.zero) then
             allocate(loclhsLag(npro,9,nshl,nshl,3))
          endif 
!
!.... compute and assemble the residual and tangent matrix
!
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac, &
                       x,                   mxmudmi(iblk)%p,       &
                       tmpshp,  &
                       tmpshgl, &
                       mien(iblk)%p, &
                       res, &
                       qres,                xKebe, &
                       xGoC,                rerr)
!
!.... satisfy the BC's on the implicit LHS
!     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1)  &
               call bc3lhs (iBC, BC,mien(iblk)%p, xKebe)  
             call fillsparseI (mien(iblk)%p,  &
                       xKebe,            lhsK, &
                       xGoC,             lhsP, &
                       rowp,                      colm)
          endif

          deallocate ( xKebe )
          deallocate ( xGoC  )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
!
!..... to calculate inner product for Lagrange Multipliers
!
       if(Lagrange.gt.zero) then
          deallocate(loclhsLag)
       endif 
!
!.... end of interior element loop
!
       enddo
!$$$       if(ibksiz.eq.20 .and. iwrote.ne.789) then
!$$$          do i=1,nshg
!$$$             write(789,*) 'eqn block ',i 
!$$$             do j=colm(i),colm(i+1)-1
!$$$                write(789,*) 'var block',rowp(j)
!$$$
!$$$                do ii=1,3
!$$$                   write(789,111) (lhsK((ii-1)*3+jj,j),jj=1,3)
!$$$                enddo
!$$$             enddo
!$$$          enddo
!$$$          close(789)
!$$$          iwrote=789
!$$$       endif
!$$$ 111   format(3(e14.7,2x))
!$$$c
!.... add in lumped mass contributions if needed
!
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassadd(ac,res,rowp,colm,lhsK,gmass)
       endif

       have_local_mass = 1
!
!.... time average statistics
!       
       if ( stsResFlg .eq. 1 ) then

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'in ')
          endif
          do j = 1,nshg
             if (btest(iBC(j),10)) then
                i = iper(j)
                stsVec(i,:) = stsVec(i,:) + stsVec(j,:)
             endif
          enddo
!     
          do i = 1,nshg
             stsVec(i,:) = stsVec(iper(i),:)
          enddo

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'out')
          endif
          return
          
       endif
!
!.... zero lhsLagL before adding contributions from the boundary elements
!
       if(Lagrange.gt.zero) then
          lhsLagL = zero
       endif 

!
!.... -------------------->   boundary elements   <--------------------
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
!.... allocate the element matrices
!
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
!
!..... to calculate inner product for Lagrange Multipliers
!
          if(Lagrange.gt.zero) then
             allocate(loclhsLag(npro,9,nshlb,nshlb,3))
          endif 
!
!.... compute and assemble the residuals corresponding to the 
!     boundary integral
!
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (u,                       y, &
                       ac,                       &
                       x, &
                       xdist, &
                       xdnv, &
                       tmpshpb, &
                       tmpshglb, &
                       mienb(iblk)%p,           mmatb(iblk)%p, &
                       miBCB(iblk)%p,           mBCB(iblk)%p, &
                       res,                     xKebe, &
                       mSWB(iblk)%p,            mTWB(iblk)%p, &
                       mEWB(iblk)%p, &
                       mPS_global(iblk)%p, &
                       mKwall_xKebe(iblk)%p)

!
!.... satisfy (again, for the vessel wall contributions) the BC's on the implicit LHS
!
!.... first, we need to make xGoC zero, since it doesn't have contributions from the 
!.... vessel wall elements

          xGoC = zero

          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1) &
                call bc3lhs (iBC, BC,mienb(iblk)%p, xKebe)
             call fillsparseI (mienb(iblk)%p, &
                       xKebe,           lhsK, &
                       xGoC,            lhsP, &
                       rowp,            colm)
          endif
!
!       
          deallocate ( xKebe )
          deallocate ( xGoC )
          deallocate (tmpshpb)
          deallocate (tmpshglb)
          if(Lagrange.gt.zero) then
             deallocate(loclhsLag)
          endif
!
!.... end of boundary element loop
!
       enddo
!
       if(Lagrange.gt.zero) then
          LagSwitch = 0 
          call CalcNANBLagrange(colm, rowp, y(:,1:3))
       endif
!       
       if(ipvsq.ge.1) then
!
!....  pressure vs. resistance boundary condition sets pressure at
!      outflow to linearly increase as flow through that face increases
!      (routine is at bottom of this file)
!
          call ElmpvsQ (res,y,-1.0d0)     
       endif
           
!
! before the commu we need to rotate the residual vector for axisymmetric
! boundary conditions (so that off processor periodicity is a dof add instead
! of a dof combination).  Take care of all nodes now so periodicity, like
! commu is a simple dof add.
!
       if(iabc==1) &             !are there any axisym bc's
             call rotabc(res, iBC,  'in ')
!
!
!.... -------------------->   communications <-------------------------
!

       if (numpe > 1) then
          call commu (res  , ilwork, nflow  , 'in ')
       endif

!
!.... ---------------------->   post processing  <----------------------
!
!.... satisfy the BCs on the residual
!
      call bc3Res (iBC,  BC,  res,  iper, ilwork)
!
!.... return
!
!      call timer ('Back    ')
      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!********************************************************************
!--------------------------------------------------------------------

      subroutine ElmGMRSclr (y,         ac,        x,      &
                             shp,       shgl,      iBC, &
                             BC,        shpb,      shglb, &
                             res,       iper,      ilwork, &
                             rowp,      colm,      lhsS    )
!
!----------------------------------------------------------------------
!
! This routine computes the LHS mass matrix, the RHS residual 
! vector, and the preconditioning matrix, for use with the GMRES
! solver.
!
!----------------------------------------------------------------------
!
        use pointer_data
        use local_mass
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        include "mpif.h"
!
        dimension y(nshg,ndof),         ac(nshg,ndof), &
                  x(numnp,nsd),         iBC(nshg),            &
                  BC(nshg,ndofBC),      res(nshg), &
                  iper(nshg)
!
        dimension shp(MAXTOP,maxsh,MAXQPT),   &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                  shpb(MAXTOP,maxsh,MAXQPT), &
                  shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!
        dimension qres(nshg,nsd),     rmass(nshg)
!
        integer ilwork(nlwork), rowp(nshg*nnz),   colm(nshg+1)

	real*8	lhsS(nnz_tot)

        real*8, allocatable, dimension(:,:,:) :: xSebe
!
!.... set up the timer
!

!AD        call timer ('Elm_Form')
!
!.... -------------------->   diffusive flux   <--------------------
!
        ires   = 1

        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
!
! loop over element blocks for the global reconstruction
! of the diffusive flux vector, q, and lumped mass matrix, rmass
!
           qres = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lcsyst = lcblk(3,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              
              ngauss = nint(lcsyst)
!     
!.... compute and assemble diffusive flux vector residual, qres,
!     and lumped mass matrix, rmass

              call AsIqSclr (y,                   x,                     &
                             shp(lcsyst,1:nshl,:),  &
                             shgl(lcsyst,:,1:nshl,:), &
                             mien(iblk)%p,     qres,                    &
                             rmass )
       
           enddo
       
!
!.... form the diffusive flux approximation
!
           call qpbcSclr ( rmass, qres, iBC, iper, ilwork )       
!
        endif 
!
!.... -------------------->   interior elements   <--------------------
!
        res    = zero
        spmass = zero

        if (lhs .eq. 1) then
           lhsS   = zero
        endif

        if ((impl(1)/10) .eq. 0) then   ! no flow solve so flxID was not zeroed
           flxID = zero
        endif
!
!.... loop over the element-blocks
!
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iel    = lcblk(1,iblk)
          lcsyst = lcblk(3,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          ndofl  = lcblk(8,iblk)
          npro   = lcblk(1,iblk+1) - iel

          ngauss = nint(lcsyst)
!
!.... allocate the element matrices
!
          allocate ( xSebe(npro,nshl,nshl) )
!
!.... compute and assemble the residual and tangent matrix
!
          call AsIGMRSclr(y,                   ac, &
                       x, &
                       shp(lcsyst,1:nshl,:),  &
                       shgl(lcsyst,:,1:nshl,:), &
                       mien(iblk)%p,        res, &
                       qres,                xSebe, mxmudmi(iblk)%p )
!
!.... satisfy the BC's on the implicit LHS
!     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             call fillsparseSclr (mien(iblk)%p,  &
                       xSebe,             lhsS, &
                       rowp,              colm)
          endif

          deallocate ( xSebe )
!
!.... end of interior element loop
!
       enddo

!
!.... add in lumped mass contributions if needed
!
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassaddSclr(ac(:,isclr), res,rowp,colm,lhsS,gmass)
       endif

       have_local_mass = 1
!
!
!  call DtN routine which updates the flux to be consistent with the
!  current solution values.  We will put the result in the last slot of
!  BC (we added a space in input.f).  That way we can localize this
!  value to the boundary elements.  This is important to keep from calling
!  the DtN evaluator more than once per node (it can be very expensive).
!
         if(idtn.eq.1)  call DtN(iBC,BC,y)
!
!.... -------------------->   boundary elements   <--------------------
!
!
!.... loop over the boundary elements
!
        do iblk = 1, nelblb
!
!.... set up the parameters
!
          iel    = lcblkb(1,iblk)
          lcsyst = lcblkb(3,iblk)
          nenl   = lcblkb(5,iblk)  ! no. of vertices per element
          nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
          nshl   = lcblkb(9,iblk)
          nshlb  = lcblkb(10,iblk)
          ndofl  = lcblkb(8,iblk)
          npro   = lcblkb(1,iblk+1) - iel

          if(lcsyst.eq.3) lcsyst=nenbl
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
!
! localize the dtn boundary condition
!

          if(idtn.eq.1)   call dtnl(   iBC, BC, mienb(iblk)%p, &
                    miBCB(iblk)%p,  mBCB(iblk)%p)

!
!.... compute and assemble the residuals corresponding to the 
!     boundary integral
!
          call AsBSclr (y,                       x, &
                        shpb(lcsyst,1:nshl,:), &
                        shglb(lcsyst,:,1:nshl,:), &
                        mienb(iblk)%p,           mmatb(iblk)%p, &
                        miBCB(iblk)%p,           mBCB(iblk)%p, &
                        res)
!
!.... end of boundary element loop
!
        enddo
!
!
!.... -------------------->   communications <-------------------------
!

      if (numpe > 1) then
        call commu (res  , ilwork, 1  , 'in ')
      endif

!
!.... ---------------------->   post processing  <----------------------
!
!.... satisfy the BCs on the residual
!
      call bc3ResSclr (iBC,  res,  iper, ilwork)
!
!.... return
!
!AD      call timer ('Back    ')
      return
      end
!
!
!....routine to compute and return the flow rates for coupled surfaces of a given type
!        
      subroutine GetFlowQ (qsurf, y, srfIdList, numSrfs)
        
      use pvsQbi  ! brings in NABI
!
       use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF), qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)
!
! note we only need the first three entries (u) from y

      qsurfProc=zero
!      
      do i = 1, nshg
         if(numSrfs .gt. zero) then
          do k = 1, numSrfs
            irankCoupled = 0
            if (srfIdList(k) .eq. ndsurf(i)) then
               irankCoupled = k
               do j = 1, 3              
                  qsurfProc(irankCoupled) = qsurfProc(irankCoupled) &
                                  + NABI(i,j)*y(i,j)
               enddo
            endif      
          enddo       
         endif
      enddo
!      
!     at this point, each qsurf has its "nodes" contributions to Q
!     accumulated into qsurf. Note, because NABI is on processor this
!     will NOT be Q for the surface yet
!
!.... reduce integrated Q for each surface, push on qsurf
!
      npars=MAXSURF+1
      call MPI_ALLREDUCE (qsurfProc, qsurf(:), npars, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr) 
!
!.... return
!
      return
      end 
!
!.... routine to compute and return the flow rates multiplied by a profile function
!.... for constrained surfaces
!        
      subroutine GetProfileFlowQ (qsurf, y, srfIdList, numSrfs)
        
      use pvsQbi  ! brings in PNABI, ndsurf
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"

      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF), qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)
      integer i, j, k
!
!.... clear the vectors 
!
      qsurfProc = zero
      do i = 1,nshg      
         if(numSrfs .gt. zero) then
            do k = 1, numSrfs
               irankCoupled = 0
               if (srfIdList(k) .eq. ndsurf(i)) then
                  irankCoupled=k
                  do j = 1, 3              
                     qsurfProc(irankCoupled) = qsurfProc(irankCoupled) &
                        +PNABI(i,j)*y(i,j)
                  enddo
               endif      
            enddo       
         endif      
      enddo
!      
!     at this point, each qsurf has its "nodes" contributions to Q
!     accumulated into qsurf. Note, because PNABI is on processor this
!     will NOT be Q for the surface yet
!
!.... reduce integrated Q for each surface, push on qsurf
!
      npars=MAXSURF+1
      call MPI_ALLREDUCE (qsurfProc, qsurf, npars, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
!
!.... return
!
      return
      end     
!
!....routine for computing inner products of velocity components for constrained surfaces.
!    Inner product is computed by calling GetInnerProduct after this routine.
!
	subroutine CalcNANBLagrange(col, row, y)
!
!.... Data declaration
!
      use LagrangeMultipliers
      use pvsQbi
 
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      integer	col(nshg+1),	    row(nnz_tot)
      real*8    y(nshg,3)
!
      real*8	tmp1,	tmp2,	tmp3
      integer	p,  i,	j,	k,  n,  m
!
!.... clear the vector
!
	if (LagSwitch .gt. 0) then 
         NANBLagrange(4:6,:,:) = zero
	else
	   NANBLagrange = zero
	endif
!
!....calculate NANBLagrange
!	
	do i = 1, nshg
	   do n=1, 3
   	      do p = 1, numLagrangeSrfs
	         tmp1 = 0
	         tmp2 = 0
	         tmp3 = 0
                 if (nsrflistLagrange(p).eq.ndsurf(i)) then 
	            do k = col(i), col(i+1)-1
		          j = row(k)
!
		          tmp1 = tmp1 &
      		         +lhsLagL(1,k,n)*y(j,1) &
      		         +lhsLagL(4,k,n)*y(j,2) &
      		         +lhsLagL(7,k,n)*y(j,3)
		          tmp2 = tmp2 &
      		         +lhsLagL(2,k,n)*y(j,1) &
      		         +lhsLagL(5,k,n)*y(j,2) &
      		         +lhsLagL(8,k,n)*y(j,3)
		          tmp3 = tmp3 &
      		         +lhsLagL(3,k,n)*y(j,1) &
      		         +lhsLagL(6,k,n)*y(j,2) &
      		         +lhsLagL(9,k,n)*y(j,3)
!
	            enddo
	            if (LagSwitch .gt. 0) then 
	               m = n+3
  	               NANBLagrange(m,i,1)=NANBLagrange(m,i,1)+tmp1
	               NANBLagrange(m,i,2)=NANBLagrange(m,i,2)+tmp2
	               NANBLagrange(m,i,3)=NANBLagrange(m,i,3)+tmp3
	            else
	               m=n
	               NANBLagrange(m,i,1)=NANBLagrange(m,i,1)+tmp1
	               NANBLagrange(m,i,2)=NANBLagrange(m,i,2)+tmp2
	               NANBLagrange(m,i,3)=NANBLagrange(m,i,3)+tmp3
	            endif
                 endif
 	      enddo
 	   enddo
 	enddo
!
      return 
      end	      
                        
!
!....routine to compute inner products for constrained surfaces.
!    CalcNANBLagrange should be called first
!        
      subroutine GetInnerProduct (qsurf, y, srfIdList, numSrfs)
!        
      use LagrangeMultipliers ! brings in NANBLagrange
      use pvsQbi  ! brings in ndsurf
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF,3), qsurfProc(0:MAXSURF,3)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)
      integer i, j, k, n
!
!.... clear the vector 
!
      qsurfProc = zero
      if(numSrfs.gt.zero) then
         do i = 1, nshg
            do n=1, 3
               do k = 1, numSrfs      
                  if (srfIdList(k) .eq. ndsurf(i)) then
                     do j=1, 3
                        qsurfProc(k,n)=qsurfProc(k,n) &
                           +NANBLagrange(n,i,j)*y(i,j)
                     enddo
                  endif
               enddo      
            enddo       
         enddo      
      endif
!      
!     at this point, each qsurf has its "nodes" contributions to Q
!     accumulated into qsurf. Note, because NABI is on processor this
!     will NOT be Q for the surface yet
!
!.... reduce integrated Q for each surface, push on qsurf
!
      do n=1, 3
         npars=MAXSURF+1
         call MPI_ALLREDUCE (qsurfProc(:,n), qsurf(:,n), npars, &
              MPI_DOUBLE_PRECISION,MPI_SUM, INEWCOMM,ierr)
      enddo  
!
!.... return
!
      return
      end      
!
!... routine to multiply 1/mu * L transpose matrix for Lagrange Multipliers
!
      subroutine LagMultiplyMatrixTranspose(srfIDList, numSrfs)     

      use pvsQbi  ! brings in NABI
      use LagrangeMultipliers !brings in the current part of coef for Lagrange Multipliers
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      real*8  DiagonalDelta,          DiagonalDeltaSurf
      integer  srfIDList(0:MAXSURF),  numSrfs 
     
      DiagonalDelta = -two*alfi*gami*Delt(1)
      LagAPproduct = zero

      do i=1, nshg
         do k = 1, numSrfs
            DiagonalDeltaSurf = zero
            DiagonalDeltaSurf = DiagonalDelta*LagMeanFlow(k)
            if (srfIDList(k).eq.ndsurf(i)) then 
               LagAPproduct(i,1:3)=LagAPproduct(i,1:3)+DiagonalDeltaSurf &
                  *((NANBLagrange(1,i,1:3)-PQLagrange(k,1)*NABI(i,1:3) &
                  -QLagrange(k,1)*PNABI(i,1:3)/LagProfileArea(k) &
                  +QLagrange(k,1)*NABI(i,1:3)*ProfileDelta(k)) &
                  *AddLag(k,1)+NANBLagrange(2,i,1:3)*AddLag(k,2) &
                  +NANBLagrange(3,i,1:3)*AddLag(k,3) )
            endif
         enddo
      enddo  
     
      return
      end

!
!... routine to multiply L matrix for Lagrange Multipliers
!
      subroutine LagMultiplyMatrix (Dy, CaseNumber, srfIDList, numSrfs) 

      use pvsQbi  ! brings in NABI
      use LagrangeMultipliers !brings in the current part of coef for Lagrange Multipliers
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      real*8  Dy(nshg,3),    DiagonalDeltaSurf
      real*8  DiagonalDelta, ProcAddLag(0:MAXSURF,3)
      integer CaseNumber, srfIDList(0:MAXSURF), numSrfs 

      DiagonalDelta = -two*alfi*gami*Delt(1)
      ProcAddLag = zero
      if (CaseNumber .eq. zero) then
         call GetFlowQ(ProcAddLag(:,1), Dy(:,1:3), srfIDList, numSrfs)  
         QLagrange(1:numSrfs,2)=ProcAddLag(1:numSrfs,1)
         ProcAddLag = zero
         call GetProfileFlowQ(ProcAddLag(:,1), Dy(:,1:3), srfIDList, &
            numSrfs)    
         PQLagrange(1:numSrfs,2)=ProcAddLag(1:numSrfs,1) &
            /LagProfileArea(1:numSrfs) 
         ProcAddLag = zero
         call GetInnerProduct(ProcAddLag, Dy(:,1:3), srfIDList, numSrfs)
         IPLagrange(1:numSrfs,4:6)=ProcAddLag(1:numSrfs,1:3)  
      endif
            
      do k = 1, numSrfs
         DiagonalDeltaSurf = zero
         DiagonalDeltaSurf = DiagonalDelta * LagMeanFlow(k)
         AddLag(k,1)=DiagonalDeltaSurf* &
            (IPLagrange(k,4)-PQLagrange(k,1)*QLagrange(k,2) &
            -QLagrange(k,1)*PQLagrange(k,2) &
            +QLagrange(k,1)*QLagrange(k,2)*ProfileDelta(k))
         AddLag(k,2)=DiagonalDeltaSurf*IPLagrange(k,5)
         AddLag(k,3)=DiagonalDeltaSurf*IPLagrange(k,6)
      enddo  
  
      return
      end

        
!
!... routine to couple pressure with flow rate for each coupled surface
!
      subroutine ElmpvsQ (res,y,sign)     

      use pvsQbi  ! brings in NABI
      use convolImpFlow !brings in the current part of convol coef for imp BC
      use convolRCRFlow !brings in the current part of convol coef for RCR BC
      use convolCORFlow !brings in the current park of convol coef for Cor BC
      use incpBC        !brings in the current part of coef for INCP BC
      use LagrangeMultipliers !brings in the current part of coef for Lagrange Multipliers
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
!
      real*8  res(nshg,ndof), y(nshg,3)
      real*8  p(0:MAXSURF),   q(0:MAXSURF,3)
      integer irankCoupled, i, j, k
!
!... get p for the resistance BC
!           
      if(numResistSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistResist,numResistSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        p(:)=sign*p(:)*ValueListResist(:) ! p=QR  now we have the true pressure on each
                                        ! outflow surface.  Note sign is -1
                                        ! for RHS, +1 for LHS
!
!....  multiply it by integral NA n_i
!     
       do i = 1,nshg
          do k = 1,numResistSrfs
              irankCoupled = 0
              if (nsrflistResist(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)
              endif
          enddo   
       enddo
       
      endif !end of coupling for Resistance BC     
!
!... get p for the impedance BC
!     
      if(numImpSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistImp,numImpSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numImpSrfs
            if(sign.lt.zero) then ! RHS so -1
               p(j)= sign*(poldImp(j) + p(j)*ImpConvCoef(ntimeptpT+2,j))  !pressure p=pold+ Qbeta
            elseif(sign.gt.zero) then ! LHS so sign is positive
               p(j)= sign*p(j)*ImpConvCoef(ntimeptpT+2,j)
            endif
        enddo             
!
!....  multiply it by integral NA n_i
!     
       do i = 1,nshg
          do k = 1,numImpSrfs
              irankCoupled = 0
              if (nsrflistImp(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)
              endif
          enddo   
       enddo
       
      endif !end of coupling for Impedance BC
!
!... get p for the RCR BC
!     
      if(numRCRSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistRCR,numRCRSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numRCRSrfs
            if(sign.lt.zero) then ! RHS so -1
                p(j)= sign*(poldRCR(j) + p(j)*RCRConvCoef(lstep+2,j)) !pressure p=pold+ Qbeta
                p(j)= p(j) - HopRCR(j) ! H operator contribution 
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*RCRConvCoef(lstep+2,j)
            endif
        enddo
             
!
!....  multiply it by integral NA n_i
!     
       do i = 1,nshg
          do k = 1,numRCRSrfs
              irankCoupled = 0
              if (nsrflistRCR(k).eq.ndsurf(i)) then 
                  irankCoupled=k
                  res(i,1:3)=res(i,1:3) + p(irankCoupled)*NABI(i,1:3)
              endif
          enddo   
       enddo       
      endif !end of coupling for RCR BC
!
!... get p for the Coronary BC
!    
      if(numCORSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistCOR,numCORSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numCORSrfs
            if(sign.lt.zero) then ! RHS so -1
                p(j)= sign*(poldCOR(j) +  &
                   p(j)*CORConvCoef(lstep+2,j)) !pressure p=pold+ Qbeta
                                !check lstep - need it to be integer and value n not n+1
                p(j)= p(j) +sign* HopCOR(j) ! H operator contribution 
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*CORConvCoef(lstep+2,j)
            endif
        enddo
!
!....  multiply it by integral NA n_i
!     
       do i = 1,nshg
          do k = 1,numCORSrfs
              irankCoupled = 0
              if (nsrflistCOR(k).eq.ndsurf(i)) then 
                 irankCoupled=k
                 res(i,1:3) = res(i,1:3) + p(irankCoupled)*NABI(i,1:3)
              endif
          enddo   
       enddo

      endif !end of coupling for Coronary BC
!
!.... get p for the coupled inflow BC
!  
      if (numINCPSrfs .gt. zero) then
         call GetFlowQ(p, y, nsrflistINCP, numINCPSrfs)
         do j=1, numINCPSrfs
            if (sign .lt. zero) then  !RHS so -1
               p(j)=sign*(INCPCoef(1,j)*p(j)+INCPCoef(2,j))
            elseif (sign .gt. zero) then
               p(j)=sign*p(j)*INCPCoef(1,j)
            endif
         enddo
!
!.... multiply p by integral NA*n_i
!
         do i=1, nshg
            do k=1, numINCPSrfs
               irankCoupled = 0
               if (nsrflistINCP(k) .ne. inactive(k)) then
                  if (nsrflistINCP(k) .eq. ndsurf(i)) then
                     irankCoupled = k
                     res(i,1:3)=res(i,1:3)+p(irankCoupled)*NABI(i,1:3)
                  endif
               endif
            enddo
         enddo
      endif  !end of coupling for INCP BC 
!
!... get p for the Lagrange multipliers
!           
      if(numLagrangeSrfs .gt. zero) then
         if(sign .lt. zero) then ! RHS so -1
            p = zero
            call GetFlowQ(p, y, nsrflistLagrange, &
               numLagrangeSrfs)  
            QLagrange(1:numLagrangeSrfs,1)=p(1:numLagrangeSrfs)
            p = zero
            call GetProfileFlowQ(p, y, nsrflistLagrange, &
               numLagrangeSrfs)    !flow rate multiplied by a profile function 
            PQLagrange(1:numLagrangeSrfs,1)=p(1:numLagrangeSrfs) &
               /LagProfileArea(1:numLagrangeSrfs) 
            q = zero
            call GetInnerProduct(q, y, nsrflistLagrange, &
                  numLagrangeSrfs)
            IPLagrange(1:numLagrangeSrfs,1:3)=q(1:numLagrangeSrfs,1:3)
            do k = 1,numLagrangeSrfs
               Penalty(k,1)=  &
                 abs( IPLagrange(k,1)-two*QLagrange(k,1)*PQLagrange(k,1) &
                  +QLagrange(k,1)**2*ProfileDelta(k) )*LagMeanFlow(k)
               Penalty(k,2)= abs(IPLagrange(k,2))*LagMeanFlow(k)
               Penalty(k,3)= abs(IPLagrange(k,3))*LagMeanFlow(k)
               resL(k,1)=two*ScaleFactor(k,1)*Lagalpha(k,1)-Penalty(k,1)
               resL(k,2)=two*ScaleFactor(k,2)*Lagalpha(k,2)-Penalty(k,2)
               resL(k,3)=two*ScaleFactor(k,3)*Lagalpha(k,3)-Penalty(k,3)
               if (lstep .eq. 0) then
                  do i=1, 3
                     LagErrorHist(1,(k-1)*3+i)=Penalty(k,i)
                  enddo
               endif
               if (numINCPSrfs .gt. zero) then
                  if (nsrflistLagrange(k).eq.inactive(k)) then !Lagrange surface ID should start from INCP IDs
                     resL(k,1:3) = zero
                  else
                     do i = 1,nshg
                        if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k) &
                              *(-Lagalpha(k,1)+PenaltyCoeff(k,1) &
                              *Penalty(k,1))*(NANBLagrange(1,i,1:3)- &
                              PQLagrange(k,1)*NABI(i,1:3)-QLagrange(k,1) &
                              *PNABI(i,1:3)/LagProfileArea(k)+ &
                              QLagrange(k,1)*NABI(i,1:3) &
                              *ProfileDelta(k))+sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,2)+PenaltyCoeff(k,2) &
                              *Penalty(k,2))*NANBLagrange(2,i,1:3) &
                              +sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,3)+PenaltyCoeff(k,3) &
                              *Penalty(k,3))*NANBLagrange(3,i,1:3)   
                        endif
                     enddo
                  endif
               else    
                  do i = 1,nshg
                     if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k) &
                              *(-Lagalpha(k,1)+PenaltyCoeff(k,1) &
                              *Penalty(k,1))*(NANBLagrange(1,i,1:3)- &
                              PQLagrange(k,1)*NABI(i,1:3)-QLagrange(k,1) &
                              *PNABI(i,1:3)/LagProfileArea(k)+ &
                              QLagrange(k,1)*NABI(i,1:3) &
                              *ProfileDelta(k))+sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,2)+PenaltyCoeff(k,2) &
                              *Penalty(k,2))*NANBLagrange(2,i,1:3) &
                              +sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,3)+PenaltyCoeff(k,3) &
                              *Penalty(k,3))*NANBLagrange(3,i,1:3)   
                     endif
                  enddo
               endif
            enddo
            AddLag(:,1:3) = resL(:,1:3)
            call LagMultiplyMatrixTranspose(nsrflistLagrange, &
               numLagrangeSrfs)
            res(:,1:3) = res(:,1:3) + LagAPproduct(:,1:3) &
               /ScaleFactor(1,1)/alfi/gami/two
         elseif(sign .gt. zero) then ! LHS 
            p = zero
            call GetFlowQ(p, y, nsrflistLagrange, &
               numLagrangeSrfs)  
            QLagrange(1:numLagrangeSrfs,2)=p(1:numLagrangeSrfs)
            p = zero
            call GetProfileFlowQ(p, y, nsrflistLagrange, &
               numLagrangeSrfs)    !flow rate multiplied by a profile function 
            PQLagrange(1:numLagrangeSrfs,2)=p(1:numLagrangeSrfs) &
               /LagProfileArea(1:numLagrangeSrfs) 
            q = zero
            call GetInnerProduct(q, y, nsrflistLagrange, &
               numLagrangeSrfs)
            IPLagrange(1:numLagrangeSrfs,4:6)=q(1:numLagrangeSrfs,1:3)  
            do k = 1, numLagrangeSrfs
               if (numINCPSrfs .gt. zero) then
                  if (nsrflistLagrange(k).eq.inactive(k)) then !order of INCP should be the same with the order of Lag
                  else
                     do i = 1,nshg
                        if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,1)+PenaltyCoeff(k,1) &
                              *Penalty(k,1))*(NANBLagrange(4,i,1:3)- &
                              NABI(i,1:3)*PQLagrange(k,2)-QLagrange(k,2) &
                              *PNABI(i,1:3)/LagProfileArea(k)+ &
                              QLagrange(k,2)*NABI(i,1:3)* &
                              ProfileDelta(k))   &
                              +sign*LagMeanFlow(k)*LagMeanFlow(k) &
                              *PenaltyCoeff(k,1)*(NANBLagrange(1,i,1:3)- &
                              NABI(i,1:3)*PQLagrange(k,1)-QLagrange(k,1) &
                              *PNABI(i,1:3)/LagProfileArea(k)+ &
                              QLagrange(k,1)*NABI(i,1:3)* &
                              ProfileDelta(k))*(IPLagrange(k,4) &
                              -PQLagrange(k,1)*QLagrange(k,2) &
                              -QLagrange(k,1)*PQLagrange(k,2) &
                              +QLagrange(k,1)*QLagrange(k,2) &
                              *ProfileDelta(k))+sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,2)+PenaltyCoeff(k,2) &
                              *Penalty(k,2))*NANBLagrange(5,i,1:3) &
                              +sign*LagMeanFlow(k)*(-Lagalpha(k,3)+ &
                              PenaltyCoeff(k,3)*Penalty(k,3)) &
                              *NANBLagrange(6,i,1:3)+sign*LagMeanFlow(k) &
                              *LagMeanFlow(k)*PenaltyCoeff(k,2)* &
                              NANBLagrange(2,i,1:3)*IPLagrange(k,5)+sign &
                              *LagMeanFlow(k)**2*PenaltyCoeff(k,3) &
                              *NANBLagrange(3,i,1:3)*IPLagrange(k,6)  
                        endif    
                     enddo
                  endif
               else      
                  do i = 1,nshg
                     if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,1)+PenaltyCoeff(k,1) &
                              *Penalty(k,1))*(NANBLagrange(4,i,1:3)- &
                              NABI(i,1:3)*PQLagrange(k,2)-QLagrange(k,2) &
                              *PNABI(i,1:3)/LagProfileArea(k)+ &
                              QLagrange(k,2)*NABI(i,1:3)* &
                              ProfileDelta(k))   &
                              +sign*LagMeanFlow(k)*LagMeanFlow(k) &
                              *PenaltyCoeff(k,1)*(NANBLagrange(1,i,1:3)- &
                              NABI(i,1:3)*PQLagrange(k,1)-QLagrange(k,1) &
                              *PNABI(i,1:3)/LagProfileArea(k)+ &
                              QLagrange(k,1)*NABI(i,1:3)* &
                              ProfileDelta(k))*(IPLagrange(k,4) &
                              -PQLagrange(k,1)*QLagrange(k,2) &
                              -QLagrange(k,1)*PQLagrange(k,2) &
                              +QLagrange(k,1)*QLagrange(k,2) &
                              *ProfileDelta(k))+sign*LagMeanFlow(k)* &
                              (-Lagalpha(k,2)+PenaltyCoeff(k,2) &
                              *Penalty(k,2))*NANBLagrange(5,i,1:3) &
                              +sign*LagMeanFlow(k)*(-Lagalpha(k,3)+ &
                              PenaltyCoeff(k,3)*Penalty(k,3)) &
                              *NANBLagrange(6,i,1:3)+sign*LagMeanFlow(k) &
                              *LagMeanFlow(k)*PenaltyCoeff(k,2)* &
                              NANBLagrange(2,i,1:3)*IPLagrange(k,5)+sign &
                              *LagMeanFlow(k)**2*PenaltyCoeff(k,3) &
                              *NANBLagrange(3,i,1:3)*IPLagrange(k,6)  
                     endif    
                  enddo
               endif
            enddo

            call LagMultiplyMatrix(y, 1, nsrflistLagrange, &
               numLagrangeSrfs)  
            do k = 1, numLagrangeSrfs
               if (numINCPSrfs .gt. zero) then
                  if (nsrflistLagrange(k).eq.inactive(k)) then !order of INCP should be the same with the order of Lag
                     AddLag(k,:) = zero
                  endif
               endif
            enddo  
            call LagMultiplyMatrixTranspose(nsrflistLagrange, &
               numLagrangeSrfs)  
            res(:,1:3) = res(:,1:3) - LagAPproduct(:,1:3)         &
               /ScaleFactor(1,1)/alfi/gami/two               
         endif
      endif !end of coupling for Lagrange multipliers
         
      return
      end