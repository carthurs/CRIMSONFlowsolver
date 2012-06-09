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
        subroutine ElmGMR (u,         y,         ac,        
     &                     x,         xdist,     xdnv,     
     &                     shp,       shgl,      iBC,
     &                     BC,        shpb,      shglb,
     &                     res,       iper,      ilwork,
     &                     rowp,      colm,      lhsK,      
     &                     lhsP,      rerr)
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c----------------------------------------------------------------------
c
        use pvsQbi  ! brings in NABI
        use stats   !  
        use pointer_data  ! brings in the pointers for the blocked arrays
        use local_mass
        use LagrangeMultipliers 
        use deformableWall
c
        include "common.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            u(nshg,nsd),
     &            x(numnp,nsd),               
     &            xdist(numnp),
     &            xdnv(numnp,nsd),
     &            iBC(nshg),           
     &            BC(nshg,ndofBC),  
     &            res(nshg,nflow),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,idflx),     rmass(nshg)
c
        dimension ilwork(nlwork)

        integer rowp(nshg*nnz),         colm(nshg+1)

	  real*8	lhsK(9,nnz_tot),	lhsP(4,nnz_tot)

        real*8, allocatable, dimension(:,:,:,:) :: xKebe, xGoC

        real*8  rerr(nshg,10)

        real*8, allocatable :: tmpshp(:,:), tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

        real*8 spmasstot(20),  ebres(nshg)
c
c.... set up the timer
c

CAD        call timer ('Elm_Form')
c
c.... -------------------->   diffusive flux   <--------------------
c
c.... set up parameters
c
        ires   = 1

        if (idiff==1 .or. idiff==3 .or. isurf==1) then ! global reconstruction
                                                       ! of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
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
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIq (y,                x,                       
     &                   shp(lcsyst,1:nshl,:), 
     &                   shgl(lcsyst,:,1:nshl,:),
     &                   mien(iblk)%p,     mxmudmi(iblk)%p,  
     &                   qres,             rmass )
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbc( rmass, qres, iBC, iper, ilwork )       
c
        endif 
c
c.... -------------------->   interior elements   <--------------------
c
        res    = zero
        if (stsResFlg .ne. 1) then
           flxID = zero
        endif

        if (lhs .eq. 1) then
           lhsp   = zero
           lhsk   = zero
        endif
c
c.... loop over the element-blocks
c
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
c
c.... allocate the element matrices
c
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
c
c..... to calculate inner product for Lagrange Multipliers
c
          if(Lagrange.gt.zero) then
             allocate(loclhsLag(npro,9,nshl,nshl,3))
          endif 
c
c.... compute and assemble the residual and tangent matrix
c
          allocate (tmpshp(nshl,MAXQPT))
          allocate (tmpshgl(nsd,nshl,MAXQPT))

          tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
          tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)

          call AsIGMR (y,                   ac,
     &                 x,                   mxmudmi(iblk)%p,      
     &                 tmpshp, 
     &                 tmpshgl,
     &                 mien(iblk)%p,
     &                 res,
     &                 qres,                xKebe,
     &                 xGoC,                rerr)
c
c.... satisfy the BC's on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1) 
     &         call bc3lhs (iBC, BC,mien(iblk)%p, xKebe)  
             call fillsparseI (mien(iblk)%p, 
     &                 xKebe,            lhsK,
     &                 xGoC,             lhsP,
     &                 rowp,                      colm)
          endif

          deallocate ( xKebe )
          deallocate ( xGoC  )
          deallocate ( tmpshp )
          deallocate ( tmpshgl )
c
c..... to calculate inner product for Lagrange Multipliers
c
       if(Lagrange.gt.zero) then
          deallocate(loclhsLag)
       endif 
c
c.... end of interior element loop
c
       enddo
c$$$       if(ibksiz.eq.20 .and. iwrote.ne.789) then
c$$$          do i=1,nshg
c$$$             write(789,*) 'eqn block ',i 
c$$$             do j=colm(i),colm(i+1)-1
c$$$                write(789,*) 'var block',rowp(j)
c$$$
c$$$                do ii=1,3
c$$$                   write(789,111) (lhsK((ii-1)*3+jj,j),jj=1,3)
c$$$                enddo
c$$$             enddo
c$$$          enddo
c$$$          close(789)
c$$$          iwrote=789
c$$$       endif
c$$$ 111   format(3(e14.7,2x))
c$$$c
c.... add in lumped mass contributions if needed
c
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassadd(ac,res,rowp,colm,lhsK,gmass)
       endif

       have_local_mass = 1
c
c.... time average statistics
c       
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
c     
          do i = 1,nshg
             stsVec(i,:) = stsVec(iper(i),:)
          enddo

          if (numpe > 1) then
             call commu (stsVec, ilwork, nResDims  , 'out')
          endif
          return
          
       endif
c
c.... zero lhsLagL before adding contributions from the boundary elements
c
       if(Lagrange.gt.zero) then
          lhsLagL = zero
       endif 

c
c.... -------------------->   boundary elements   <--------------------
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
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
c
          if(lcsyst.eq.3 .or. lcsyst.eq.4) then
             ngaussb = nintb(lcsyst)
          else
             ngaussb = nintb(lcsyst)
          endif
c
c.... allocate the element matrices
c
          allocate ( xKebe(npro,9,nshl,nshl) )
          allocate ( xGoC (npro,4,nshl,nshl) )
c
c..... to calculate inner product for Lagrange Multipliers
c
          if(Lagrange.gt.zero) then
             allocate(loclhsLag(npro,9,nshlb,nshlb,3))
          endif 
c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          allocate (tmpshpb(nshl,MAXQPT))
          allocate (tmpshglb(nsd,nshl,MAXQPT))
          
          tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
          tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

          call AsBMFG (u,                       y,
     &                 ac,                      
     &                 x,
     &                 xdist,
     &                 xdnv,
     &                 tmpshpb,
     &                 tmpshglb,
     &                 mienb(iblk)%p,           mmatb(iblk)%p,
     &                 miBCB(iblk)%p,           mBCB(iblk)%p,
     &                 res,                     xKebe,
     &                 mSWB(iblk)%p,            mTWB(iblk)%p,
     &                 mEWB(iblk)%p,
     &                 mPS_global(iblk)%p,
     &                 mKwall_xKebe(iblk)%p)

c
c.... satisfy (again, for the vessel wall contributions) the BC's on the implicit LHS
c
c.... first, we need to make xGoC zero, since it doesn't have contributions from the 
c.... vessel wall elements

          xGoC = zero

          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             if(ipord.eq.1)
     &          call bc3lhs (iBC, BC,mienb(iblk)%p, xKebe)
             call fillsparseI (mienb(iblk)%p,
     &                 xKebe,           lhsK,
     &                 xGoC,            lhsP,
     &                 rowp,            colm)
          endif
c
c       
          deallocate ( xKebe )
          deallocate ( xGoC )
          deallocate (tmpshpb)
          deallocate (tmpshglb)
          if(Lagrange.gt.zero) then
             deallocate(loclhsLag)
          endif
c
c.... end of boundary element loop
c
       enddo
c
       if(Lagrange.gt.zero) then
          LagSwitch = 0 
          call CalcNANBLagrange(colm, rowp, y(:,1:3))
       endif
c       
       if(ipvsq.ge.1) then
c
c....  pressure vs. resistance boundary condition sets pressure at
c      outflow to linearly increase as flow through that face increases
c      (routine is at bottom of this file)
c
          call ElmpvsQ (res,y,-1.0d0)     
       endif
           
c
c before the commu we need to rotate the residual vector for axisymmetric
c boundary conditions (so that off processor periodicity is a dof add instead
c of a dof combination).  Take care of all nodes now so periodicity, like
c commu is a simple dof add.
c
       if(iabc==1)              !are there any axisym bc's
     &       call rotabc(res, iBC,  'in ')
c
c
c.... -------------------->   communications <-------------------------
c

       if (numpe > 1) then
          call commu (res  , ilwork, nflow  , 'in ')
       endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3Res (iBC,  BC,  res,  iper, ilwork)
c
c.... return
c
c      call timer ('Back    ')
      return
      end


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!********************************************************************
!--------------------------------------------------------------------

      subroutine ElmGMRSclr (y,         ac,        x,     
     &                       shp,       shgl,      iBC,
     &                       BC,        shpb,      shglb,
     &                       res,       iper,      ilwork,
     &                       rowp,      colm,      lhsS    )
c
c----------------------------------------------------------------------
c
c This routine computes the LHS mass matrix, the RHS residual 
c vector, and the preconditioning matrix, for use with the GMRES
c solver.
c
c----------------------------------------------------------------------
c
        use pointer_data
        use local_mass
c
        include "common.h"
        include "mpif.h"
c
        dimension y(nshg,ndof),         ac(nshg,ndof),
     &            x(numnp,nsd),         iBC(nshg),           
     &            BC(nshg,ndofBC),      res(nshg),
     &            iper(nshg)
c
        dimension shp(MAXTOP,maxsh,MAXQPT),  
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT), 
     &            shpb(MAXTOP,maxsh,MAXQPT),
     &            shglb(MAXTOP,nsd,maxsh,MAXQPT) 
c
        dimension qres(nshg,nsd),     rmass(nshg)
c
        integer ilwork(nlwork), rowp(nshg*nnz),   colm(nshg+1)

	real*8	lhsS(nnz_tot)

        real*8, allocatable, dimension(:,:,:) :: xSebe
c
c.... set up the timer
c

CAD        call timer ('Elm_Form')
c
c.... -------------------->   diffusive flux   <--------------------
c
        ires   = 1

        if (idiff==1 .or. idiff==3) then ! global reconstruction of qdiff
c
c loop over element blocks for the global reconstruction
c of the diffusive flux vector, q, and lumped mass matrix, rmass
c
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
c     
c.... compute and assemble diffusive flux vector residual, qres,
c     and lumped mass matrix, rmass

              call AsIqSclr (y,                   x,                    
     &                       shp(lcsyst,1:nshl,:), 
     &                       shgl(lcsyst,:,1:nshl,:),
     &                       mien(iblk)%p,     qres,                   
     &                       rmass )
       
           enddo
       
c
c.... form the diffusive flux approximation
c
           call qpbcSclr ( rmass, qres, iBC, iper, ilwork )       
c
        endif 
c
c.... -------------------->   interior elements   <--------------------
c
        res    = zero
        spmass = zero

        if (lhs .eq. 1) then
           lhsS   = zero
        endif

        if ((impl(1)/10) .eq. 0) then   ! no flow solve so flxID was not zeroed
           flxID = zero
        endif
c
c.... loop over the element-blocks
c
        do iblk = 1, nelblk
          iblock = iblk         ! used in local mass inverse (p>2)
          iel    = lcblk(1,iblk)
          lcsyst = lcblk(3,iblk)
          nenl   = lcblk(5,iblk) ! no. of vertices per element
          nshl   = lcblk(10,iblk)
          ndofl  = lcblk(8,iblk)
          npro   = lcblk(1,iblk+1) - iel

          ngauss = nint(lcsyst)
c
c.... allocate the element matrices
c
          allocate ( xSebe(npro,nshl,nshl) )
c
c.... compute and assemble the residual and tangent matrix
c
          call AsIGMRSclr(y,                   ac,
     &                 x,
     &                 shp(lcsyst,1:nshl,:), 
     &                 shgl(lcsyst,:,1:nshl,:),
     &                 mien(iblk)%p,        res,
     &                 qres,                xSebe, mxmudmi(iblk)%p )
c
c.... satisfy the BC's on the implicit LHS
c     
          if (impl(1) .ne. 9 .and. lhs .eq. 1) then
             call fillsparseSclr (mien(iblk)%p, 
     &                 xSebe,             lhsS,
     &                 rowp,              colm)
          endif

          deallocate ( xSebe )
c
c.... end of interior element loop
c
       enddo

c
c.... add in lumped mass contributions if needed
c
       if((flmpr.ne.0).or.(flmpl.ne.0)) then
          call lmassaddSclr(ac(:,isclr), res,rowp,colm,lhsS,gmass)
       endif

       have_local_mass = 1
c
c
c  call DtN routine which updates the flux to be consistent with the
c  current solution values.  We will put the result in the last slot of
c  BC (we added a space in input.f).  That way we can localize this
c  value to the boundary elements.  This is important to keep from calling
c  the DtN evaluator more than once per node (it can be very expensive).
c
         if(idtn.eq.1)  call DtN(iBC,BC,y)
c
c.... -------------------->   boundary elements   <--------------------
c
c
c.... loop over the boundary elements
c
        do iblk = 1, nelblb
c
c.... set up the parameters
c
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
c
c localize the dtn boundary condition
c

          if(idtn.eq.1)   call dtnl(   iBC, BC, mienb(iblk)%p,
     &              miBCB(iblk)%p,  mBCB(iblk)%p)

c
c.... compute and assemble the residuals corresponding to the 
c     boundary integral
c
          call AsBSclr (y,                       x,
     &                  shpb(lcsyst,1:nshl,:),
     &                  shglb(lcsyst,:,1:nshl,:),
     &                  mienb(iblk)%p,           mmatb(iblk)%p,
     &                  miBCB(iblk)%p,           mBCB(iblk)%p,
     &                  res)
c
c.... end of boundary element loop
c
        enddo
c
c
c.... -------------------->   communications <-------------------------
c

      if (numpe > 1) then
        call commu (res  , ilwork, 1  , 'in ')
      endif

c
c.... ---------------------->   post processing  <----------------------
c
c.... satisfy the BCs on the residual
c
      call bc3ResSclr (iBC,  res,  iper, ilwork)
c
c.... return
c
CAD      call timer ('Back    ')
      return
      end
c
c
c....routine to compute and return the flow rates for coupled surfaces of a given type
c        
      subroutine GetFlowQ (qsurf, y, srfIdList, numSrfs)
        
      use pvsQbi  ! brings in NABI
c
      include "common.h"
      include "mpif.h"
c
      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF), qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)
c
c note we only need the first three entries (u) from y

      qsurfProc=zero
c      
      do i = 1, nshg
         if(numSrfs .gt. zero) then
          do k = 1, numSrfs
            irankCoupled = 0
            if (srfIdList(k) .eq. ndsurf(i)) then
               irankCoupled = k
               do j = 1, 3              
                  qsurfProc(irankCoupled) = qsurfProc(irankCoupled)
     &                            + NABI(i,j)*y(i,j)
               enddo
            endif      
          enddo       
         endif
      enddo
c      
c     at this point, each qsurf has its "nodes" contributions to Q
c     accumulated into qsurf. Note, because NABI is on processor this
c     will NOT be Q for the surface yet
c
c.... reduce integrated Q for each surface, push on qsurf
c
      npars=MAXSURF+1
      call MPI_ALLREDUCE (qsurfProc, qsurf(:), npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr) 
c
c.... return
c
      return
      end 
c
c.... routine to compute and return the flow rates multiplied by a profile function
c.... for constrained surfaces
c        
      subroutine GetProfileFlowQ (qsurf, y, srfIdList, numSrfs)
        
      use pvsQbi  ! brings in PNABI, ndsurf
c
      include "common.h"
      include "mpif.h"

      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF), qsurfProc(0:MAXSURF)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)
      integer i, j, k
c
c.... clear the vectors 
c
      qsurfProc = zero
      do i = 1,nshg      
         if(numSrfs .gt. zero) then
            do k = 1, numSrfs
               irankCoupled = 0
               if (srfIdList(k) .eq. ndsurf(i)) then
                  irankCoupled=k
                  do j = 1, 3              
                     qsurfProc(irankCoupled) = qsurfProc(irankCoupled)
     &                  +PNABI(i,j)*y(i,j)
                  enddo
               endif      
            enddo       
         endif      
      enddo
c      
c     at this point, each qsurf has its "nodes" contributions to Q
c     accumulated into qsurf. Note, because PNABI is on processor this
c     will NOT be Q for the surface yet
c
c.... reduce integrated Q for each surface, push on qsurf
c
      npars=MAXSURF+1
      call MPI_ALLREDUCE (qsurfProc, qsurf, npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
c
c.... return
c
      return
      end     
c
c....routine for computing inner products of velocity components for constrained surfaces.
c    Inner product is computed by calling GetInnerProduct after this routine.
c
	subroutine CalcNANBLagrange(col, row, y)
c
c.... Data declaration
c
      use LagrangeMultipliers
      use pvsQbi
 
      include "common.h"
c
	integer	col(nshg+1),	    row(nnz_tot)
      real*8    y(nshg,3)
c
	real*8	tmp1,	tmp2,	tmp3
	integer	p,  i,	j,	k,  n,  m
c
c.... clear the vector
c
	if (LagSwitch .gt. 0) then 
         NANBLagrange(4:6,:,:) = zero
	else
	   NANBLagrange = zero
	endif
c
c....calculate NANBLagrange
c	
	do i = 1, nshg
	   do n=1, 3
   	      do p = 1, numLagrangeSrfs
	         tmp1 = 0
	         tmp2 = 0
	         tmp3 = 0
                 if (nsrflistLagrange(p).eq.ndsurf(i)) then 
	            do k = col(i), col(i+1)-1
		          j = row(k)
c
		          tmp1 = tmp1
     1		         +lhsLagL(1,k,n)*y(j,1)
     2		         +lhsLagL(4,k,n)*y(j,2)
     3		         +lhsLagL(7,k,n)*y(j,3)
		          tmp2 = tmp2
     1		         +lhsLagL(2,k,n)*y(j,1)
     2		         +lhsLagL(5,k,n)*y(j,2)
     3		         +lhsLagL(8,k,n)*y(j,3)
		          tmp3 = tmp3
     1		         +lhsLagL(3,k,n)*y(j,1)
     2		         +lhsLagL(6,k,n)*y(j,2)
     3		         +lhsLagL(9,k,n)*y(j,3)
c
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
c
      return 
      end	      
                        
c
c....routine to compute inner products for constrained surfaces.
c    CalcNANBLagrange should be called first
c        
      subroutine GetInnerProduct (qsurf, y, srfIdList, numSrfs)
c        
      use LagrangeMultipliers ! brings in NANBLagrange
      use pvsQbi  ! brings in ndsurf
c
      include "common.h"
      include "mpif.h"
c
      real*8  y(nshg,3)
      real*8  qsurf(0:MAXSURF,3), qsurfProc(0:MAXSURF,3)
      integer numSrfs, irankCoupled, srfIdList(0:MAXSURF)
      integer i, j, k, n
c
c.... clear the vector 
c
      qsurfProc = zero
      if(numSrfs.gt.zero) then
         do i = 1, nshg
            do n=1, 3
               do k = 1, numSrfs      
                  if (srfIdList(k) .eq. ndsurf(i)) then
                     do j=1, 3
                        qsurfProc(k,n)=qsurfProc(k,n)
     &                     +NANBLagrange(n,i,j)*y(i,j)
                     enddo
                  endif
               enddo      
            enddo       
         enddo      
      endif
c      
c     at this point, each qsurf has its "nodes" contributions to Q
c     accumulated into qsurf. Note, because NABI is on processor this
c     will NOT be Q for the surface yet
c
c.... reduce integrated Q for each surface, push on qsurf
c
      do n=1, 3
         npars=MAXSURF+1
         call MPI_ALLREDUCE (qsurfProc(:,n), qsurf(:,n), npars,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
      enddo  
c
c.... return
c
      return
      end      
c
c... routine to multiply 1/mu * L transpose matrix for Lagrange Multipliers
c
      subroutine LagMultiplyMatrixTranspose(srfIDList, numSrfs)     

      use pvsQbi  ! brings in NABI
      use LagrangeMultipliers !brings in the current part of coef for Lagrange Multipliers
c
      include "common.h"
      include "mpif.h"
c
      real*8  DiagonalDelta,          DiagonalDeltaSurf
      integer  srfIDList(0:MAXSURF),  numSrfs 
     
      DiagonalDelta = -two*alfi*gami*Delt(1)
      LagAPproduct = zero

      do i=1, nshg
         do k = 1, numSrfs
            DiagonalDeltaSurf = zero
            DiagonalDeltaSurf = DiagonalDelta*LagMeanFlow(k)
            if (srfIDList(k).eq.ndsurf(i)) then 
               LagAPproduct(i,1:3)=LagAPproduct(i,1:3)+DiagonalDeltaSurf
     &            *((NANBLagrange(1,i,1:3)-PQLagrange(k,1)*NABI(i,1:3)
     &            -QLagrange(k,1)*PNABI(i,1:3)/LagProfileArea(k)
     &            +QLagrange(k,1)*NABI(i,1:3)*ProfileDelta(k))
     &            *AddLag(k,1)+NANBLagrange(2,i,1:3)*AddLag(k,2)
     &            +NANBLagrange(3,i,1:3)*AddLag(k,3) )
            endif
         enddo
      enddo  
     
      return
      end

c
c... routine to multiply L matrix for Lagrange Multipliers
c
      subroutine LagMultiplyMatrix (Dy, CaseNumber, srfIDList, numSrfs) 

      use pvsQbi  ! brings in NABI
      use LagrangeMultipliers !brings in the current part of coef for Lagrange Multipliers
c
      include "common.h"
      include "mpif.h"
c
      real*8  Dy(nshg,3),    DiagonalDeltaSurf
      real*8  DiagonalDelta, ProcAddLag(0:MAXSURF,3)
      integer CaseNumber, srfIDList(0:MAXSURF), numSrfs 

      DiagonalDelta = -two*alfi*gami*Delt(1)
      ProcAddLag = zero
      if (CaseNumber .eq. zero) then
         call GetFlowQ(ProcAddLag(:,1), Dy(:,1:3), srfIDList, numSrfs)  
         QLagrange(1:numSrfs,2)=ProcAddLag(1:numSrfs,1)
         ProcAddLag = zero
         call GetProfileFlowQ(ProcAddLag(:,1), Dy(:,1:3), srfIDList,
     &      numSrfs)    
         PQLagrange(1:numSrfs,2)=ProcAddLag(1:numSrfs,1)
     &      /LagProfileArea(1:numSrfs) 
         ProcAddLag = zero
         call GetInnerProduct(ProcAddLag, Dy(:,1:3), srfIDList, numSrfs)
         IPLagrange(1:numSrfs,4:6)=ProcAddLag(1:numSrfs,1:3)  
      endif
            
      do k = 1, numSrfs
         DiagonalDeltaSurf = zero
         DiagonalDeltaSurf = DiagonalDelta * LagMeanFlow(k)
         AddLag(k,1)=DiagonalDeltaSurf*
     &      (IPLagrange(k,4)-PQLagrange(k,1)*QLagrange(k,2)
     &      -QLagrange(k,1)*PQLagrange(k,2)
     &      +QLagrange(k,1)*QLagrange(k,2)*ProfileDelta(k))
         AddLag(k,2)=DiagonalDeltaSurf*IPLagrange(k,5)
         AddLag(k,3)=DiagonalDeltaSurf*IPLagrange(k,6)
      enddo  
  
      return
      end

        
c
c... routine to couple pressure with flow rate for each coupled surface
c
      subroutine ElmpvsQ (res,y,sign)     

      use pvsQbi  ! brings in NABI
      use convolImpFlow !brings in the current part of convol coef for imp BC
      use convolRCRFlow !brings in the current part of convol coef for RCR BC
      use convolCORFlow !brings in the current park of convol coef for Cor BC
      use incpBC        !brings in the current part of coef for INCP BC
      use LagrangeMultipliers !brings in the current part of coef for Lagrange Multipliers
c
      include "common.h"
      include "mpif.h"
c
      real*8  res(nshg,ndof), y(nshg,3)
      real*8  p(0:MAXSURF),   q(0:MAXSURF,3)
      integer irankCoupled, i, j, k
c
c... get p for the resistance BC
c           
      if(numResistSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistResist,numResistSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        p(:)=sign*p(:)*ValueListResist(:) ! p=QR  now we have the true pressure on each
                                        ! outflow surface.  Note sign is -1
                                        ! for RHS, +1 for LHS
c
c....  multiply it by integral NA n_i
c     
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
c
c... get p for the impedance BC
c     
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
c
c....  multiply it by integral NA n_i
c     
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
c
c... get p for the RCR BC
c     
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
             
c
c....  multiply it by integral NA n_i
c     
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
c
c... get p for the Coronary BC
c    
      if(numCORSrfs.gt.zero) then
        call GetFlowQ(p,y,nsrflistCOR,numCORSrfs)  !Q pushed into p but at this point 
                          ! p is just the full Q for each surface
        do j = 1,numCORSrfs
            if(sign.lt.zero) then ! RHS so -1
                p(j)= sign*(poldCOR(j) + 
     &             p(j)*CORConvCoef(lstep+2,j)) !pressure p=pold+ Qbeta
                                !check lstep - need it to be integer and value n not n+1
                p(j)= p(j) +sign* HopCOR(j) ! H operator contribution 
            elseif(sign.gt.zero) then ! LHS so sign is positive
                p(j)= sign*p(j)*CORConvCoef(lstep+2,j)
            endif
        enddo
c
c....  multiply it by integral NA n_i
c     
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
c
c.... get p for the coupled inflow BC
c  
      if (numINCPSrfs .gt. zero) then
         call GetFlowQ(p, y, nsrflistINCP, numINCPSrfs)
         do j=1, numINCPSrfs
            if (sign .lt. zero) then  !RHS so -1
               p(j)=sign*(INCPCoef(1,j)*p(j)+INCPCoef(2,j))
            elseif (sign .gt. zero) then
               p(j)=sign*p(j)*INCPCoef(1,j)
            endif
         enddo
c
c.... multiply p by integral NA*n_i
c
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
c
c... get p for the Lagrange multipliers
c           
      if(numLagrangeSrfs .gt. zero) then
         if(sign .lt. zero) then ! RHS so -1
            p = zero
            call GetFlowQ(p, y, nsrflistLagrange,
     &         numLagrangeSrfs)  
            QLagrange(1:numLagrangeSrfs,1)=p(1:numLagrangeSrfs)
            p = zero
            call GetProfileFlowQ(p, y, nsrflistLagrange,
     &         numLagrangeSrfs)    !flow rate multiplied by a profile function 
            PQLagrange(1:numLagrangeSrfs,1)=p(1:numLagrangeSrfs)
     &         /LagProfileArea(1:numLagrangeSrfs) 
            q = zero
            call GetInnerProduct(q, y, nsrflistLagrange,
     &            numLagrangeSrfs)
            IPLagrange(1:numLagrangeSrfs,1:3)=q(1:numLagrangeSrfs,1:3)
            do k = 1,numLagrangeSrfs
               Penalty(k,1)= 
     &           abs( IPLagrange(k,1)-two*QLagrange(k,1)*PQLagrange(k,1)
     &            +QLagrange(k,1)**2*ProfileDelta(k) )*LagMeanFlow(k)
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
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k)
     &                        *(-Lagalpha(k,1)+PenaltyCoeff(k,1)
     &                        *Penalty(k,1))*(NANBLagrange(1,i,1:3)-
     &                        PQLagrange(k,1)*NABI(i,1:3)-QLagrange(k,1)
     &                        *PNABI(i,1:3)/LagProfileArea(k)+
     &                        QLagrange(k,1)*NABI(i,1:3)
     &                        *ProfileDelta(k))+sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,2)+PenaltyCoeff(k,2)
     &                        *Penalty(k,2))*NANBLagrange(2,i,1:3)
     &                        +sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,3)+PenaltyCoeff(k,3)
     &                        *Penalty(k,3))*NANBLagrange(3,i,1:3)   
                        endif
                     enddo
                  endif
               else    
                  do i = 1,nshg
                     if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k)
     &                        *(-Lagalpha(k,1)+PenaltyCoeff(k,1)
     &                        *Penalty(k,1))*(NANBLagrange(1,i,1:3)-
     &                        PQLagrange(k,1)*NABI(i,1:3)-QLagrange(k,1)
     &                        *PNABI(i,1:3)/LagProfileArea(k)+
     &                        QLagrange(k,1)*NABI(i,1:3)
     &                        *ProfileDelta(k))+sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,2)+PenaltyCoeff(k,2)
     &                        *Penalty(k,2))*NANBLagrange(2,i,1:3)
     &                        +sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,3)+PenaltyCoeff(k,3)
     &                        *Penalty(k,3))*NANBLagrange(3,i,1:3)   
                     endif
                  enddo
               endif
            enddo
            AddLag(:,1:3) = resL(:,1:3)
            call LagMultiplyMatrixTranspose(nsrflistLagrange,
     &         numLagrangeSrfs)
            res(:,1:3) = res(:,1:3) + LagAPproduct(:,1:3)
     &         /ScaleFactor(1,1)/alfi/gami/two
         elseif(sign .gt. zero) then ! LHS 
            p = zero
            call GetFlowQ(p, y, nsrflistLagrange,
     &         numLagrangeSrfs)  
            QLagrange(1:numLagrangeSrfs,2)=p(1:numLagrangeSrfs)
            p = zero
            call GetProfileFlowQ(p, y, nsrflistLagrange,
     &         numLagrangeSrfs)    !flow rate multiplied by a profile function 
            PQLagrange(1:numLagrangeSrfs,2)=p(1:numLagrangeSrfs)
     &         /LagProfileArea(1:numLagrangeSrfs) 
            q = zero
            call GetInnerProduct(q, y, nsrflistLagrange,
     &         numLagrangeSrfs)
            IPLagrange(1:numLagrangeSrfs,4:6)=q(1:numLagrangeSrfs,1:3)  
            do k = 1, numLagrangeSrfs
               if (numINCPSrfs .gt. zero) then
                  if (nsrflistLagrange(k).eq.inactive(k)) then !order of INCP should be the same with the order of Lag
                  else
                     do i = 1,nshg
                        if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,1)+PenaltyCoeff(k,1)
     &                        *Penalty(k,1))*(NANBLagrange(4,i,1:3)-
     &                        NABI(i,1:3)*PQLagrange(k,2)-QLagrange(k,2)
     &                        *PNABI(i,1:3)/LagProfileArea(k)+
     &                        QLagrange(k,2)*NABI(i,1:3)*
     &                        ProfileDelta(k))  
     &                        +sign*LagMeanFlow(k)*LagMeanFlow(k)
     &                        *PenaltyCoeff(k,1)*(NANBLagrange(1,i,1:3)-
     &                        NABI(i,1:3)*PQLagrange(k,1)-QLagrange(k,1)
     &                        *PNABI(i,1:3)/LagProfileArea(k)+
     &                        QLagrange(k,1)*NABI(i,1:3)*
     &                        ProfileDelta(k))*(IPLagrange(k,4)
     &                        -PQLagrange(k,1)*QLagrange(k,2)
     &                        -QLagrange(k,1)*PQLagrange(k,2)
     &                        +QLagrange(k,1)*QLagrange(k,2)
     &                        *ProfileDelta(k))+sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,2)+PenaltyCoeff(k,2)
     &                        *Penalty(k,2))*NANBLagrange(5,i,1:3)
     &                        +sign*LagMeanFlow(k)*(-Lagalpha(k,3)+
     &                        PenaltyCoeff(k,3)*Penalty(k,3))
     &                        *NANBLagrange(6,i,1:3)+sign*LagMeanFlow(k)
     &                        *LagMeanFlow(k)*PenaltyCoeff(k,2)*
     &                        NANBLagrange(2,i,1:3)*IPLagrange(k,5)+sign
     &                        *LagMeanFlow(k)**2*PenaltyCoeff(k,3)
     &                        *NANBLagrange(3,i,1:3)*IPLagrange(k,6)  
                        endif    
                     enddo
                  endif
               else      
                  do i = 1,nshg
                     if (nsrflistLagrange(k).eq.ndsurf(i)) then 
                          res(i,1:3)=res(i,1:3)+sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,1)+PenaltyCoeff(k,1)
     &                        *Penalty(k,1))*(NANBLagrange(4,i,1:3)-
     &                        NABI(i,1:3)*PQLagrange(k,2)-QLagrange(k,2)
     &                        *PNABI(i,1:3)/LagProfileArea(k)+
     &                        QLagrange(k,2)*NABI(i,1:3)*
     &                        ProfileDelta(k))  
     &                        +sign*LagMeanFlow(k)*LagMeanFlow(k)
     &                        *PenaltyCoeff(k,1)*(NANBLagrange(1,i,1:3)-
     &                        NABI(i,1:3)*PQLagrange(k,1)-QLagrange(k,1)
     &                        *PNABI(i,1:3)/LagProfileArea(k)+
     &                        QLagrange(k,1)*NABI(i,1:3)*
     &                        ProfileDelta(k))*(IPLagrange(k,4)
     &                        -PQLagrange(k,1)*QLagrange(k,2)
     &                        -QLagrange(k,1)*PQLagrange(k,2)
     &                        +QLagrange(k,1)*QLagrange(k,2)
     &                        *ProfileDelta(k))+sign*LagMeanFlow(k)*
     &                        (-Lagalpha(k,2)+PenaltyCoeff(k,2)
     &                        *Penalty(k,2))*NANBLagrange(5,i,1:3)
     &                        +sign*LagMeanFlow(k)*(-Lagalpha(k,3)+
     &                        PenaltyCoeff(k,3)*Penalty(k,3))
     &                        *NANBLagrange(6,i,1:3)+sign*LagMeanFlow(k)
     &                        *LagMeanFlow(k)*PenaltyCoeff(k,2)*
     &                        NANBLagrange(2,i,1:3)*IPLagrange(k,5)+sign
     &                        *LagMeanFlow(k)**2*PenaltyCoeff(k,3)
     &                        *NANBLagrange(3,i,1:3)*IPLagrange(k,6)  
                     endif    
                  enddo
               endif
            enddo

            call LagMultiplyMatrix(y, 1, nsrflistLagrange,
     &         numLagrangeSrfs)  
            do k = 1, numLagrangeSrfs
               if (numINCPSrfs .gt. zero) then
                  if (nsrflistLagrange(k).eq.inactive(k)) then !order of INCP should be the same with the order of Lag
                     AddLag(k,:) = zero
                  endif
               endif
            enddo  
            call LagMultiplyMatrixTranspose(nsrflistLagrange,
     &         numLagrangeSrfs)  
            res(:,1:3) = res(:,1:3) - LagAPproduct(:,1:3)        
     &         /ScaleFactor(1,1)/alfi/gami/two               
         endif
      endif !end of coupling for Lagrange multipliers
         
      return
      end
