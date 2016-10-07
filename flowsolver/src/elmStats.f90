      subroutine elmStatsLhs( x,  iBC,   iper,  ilwork )
!-----------------------------------------------------------------------
!     compute the necessary terms for the statistics projection
!     matrices.
!-----------------------------------------------------------------------
      use     stats
      use     pointer_data
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      real*8  x(numnp,3)
      integer iBC(nshg), iper(nshg), ilwork(nlwork)
      
      real*8, allocatable :: xl(:,:,:)
      real*8, allocatable :: lStsVec(:,:,:)

!
!.... loop over element blocks
!
      stsVec = zero
      
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         lcsyst = lcblk(3,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         ndofl  = lcblk(8,iblk)
         npro   = lcblk(1,iblk+1) - iel 

         allocate ( xl(npro,nenl,3)             )
         allocate ( lStsVec(npro,nshl,nResDims) )
!
!.... localize needed data
!
         call localx ( x,    xl,  mien(iblk)%p, nsd,   'gather  ' )
!
!.... form the Lhs
!
         call e3StsLhs( xl, lStsVec )
!
!.... assemble
!
         call local (stsVec, lStsVec, mien(iblk)%p, &
                     nResDims, 'scatter ' ) 

         deallocate ( xl       )
         deallocate ( lStsVec  )
!
!.... end loop over element blocks
!
      enddo

      if (numpe > 1) then
        call commu (stsVec, ilwork, nResDims  , 'in ')
      endif
!
!.... local periodic boundary conditions (no communications)
!
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
      end
      
!-----------------------------------------------------------------------
!  Assemble the residual for the statistics
!-----------------------------------------------------------------------
      subroutine elmStatsRes( u,        y,           ac,     &
                              x,        xdist,       xdnv, &
                              shp,      shgl,  &
                              shpb,     shglb,       iBC,     BC,  &
                              iper,     ilwork,      rowp,    colm, &
                              lhsK,     lhsP, &
                              uMesh   ) !ALE variables added MAF 06/10/2016)
      
      use stats
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      
      real*8  y(nshg,ndof),             ac(nshg,ndof),  &
              u(nshg,nsd), &
              x(numnp,nsd),            &
              xdist(nshg), &
              xdnv(nshg,nsd), &
              shp(MAXTOP,maxsh,MAXQPT), shgl(MAXTOP,nsd,maxsh,MAXQPT), &
              shpb(MAXTOP,maxsh,MAXQPT), &
              shglb(MAXTOP,nsd,maxsh,MAXQPT), &
              BC(nshg,ndofBC),          lhsK(9,nnz_tot), &
              lhsP(4,nnz_tot),          res(nshg,ndof)

      integer iBC(nshg),                iper(nshg), &
              ilwork(nlwork),           rowp(nshg,nnz), &
              colm(nshg+1)

      real*8 uMesh(nshg,3)           

      lhs    = 0
      stsVec = zero
      
      stsResFlg = 1
      ierrcalctmp=ierrcalc ! save the current value of ierrcalc
      ierrcalc=0           ! set it to zero so that we don't calculate error
                           ! here (which would overflow memory around rjunk)
      call ElmGMR (u,         y,         ac,          &
                   x,         xdist,      xdnv, &
                   shp,       shgl,       iBC,        &
                   BC,        shpb,       shglb, &
                   res,       iper,       ilwork,    &
                   rowp,      colm,       lhsK,       &
                   lhsP,      rjunk, &
                   uMesh   ) !ALE variables added MAF 06/10/2016)
      stsResFlg = 0
      ierrcalc=ierrcalctmp  ! reset it back to original value
      return 
      end

