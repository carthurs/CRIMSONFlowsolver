      subroutine genadj (colm, rowp, icnt )
!     
      use pointer_data
!     
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      integer rowp(nshg*nnz),         colm(nshg+1)
      integer icnt
      
      integer adjcnt(nshg),    row_fill_list(nshg,6*nnz), mloc(1)
!                                          change ^ if overflow
!                                   also change overflow check in asadj TWICE
      integer tmprdim(1)
      real*8, allocatable, dimension(:) :: tmpr

      adjcnt=0
      
      do iblk = 1, nelblk
!     
!.... set up the parameters
!     
         iel    = lcblk(1,iblk)
         lelCat = lcblk(2,iblk)
         lcsyst = lcblk(3,iblk)
         iorder = lcblk(4,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         npro   = lcblk(1,iblk+1) - iel 
         
!     
!.... compute sparse matrix data structures
!     
         call Asadj (row_fill_list, &
                     mien(iblk)%p,  adjcnt )
         
      enddo
      
      call sumgatInt ( adjcnt, nshg, nnonzero)
      if ( myrank .eq. master) then
         write (*,*) 'Number of global nonzeros ',nnonzero
      endif

!     
!     build the colm array
!     
      colm(1)=1
      do i=1,nshg
         colm(i+1)=colm(i)+adjcnt(i)
      enddo
!     
!     sort the rowp into increasing order
!     
      ibig=10*nshg
      icnt=0
      tmprdim=maxval(adjcnt)
      allocate (tmpr(tmprdim(1)))
      do i=1,nshg
         ncol=adjcnt(i)
         tmpr(1:ncol)=row_fill_list(i,1:ncol)
         do j=1,ncol
            icnt=icnt+1
            imin=minval(tmpr(1:ncol))
            mloc=minloc(tmpr(1:ncol))
            rowp(icnt)=imin
            tmpr(mloc(1))=ibig
         enddo
      enddo

      maxfill=tmprdim(1)
      write(*,*) 'maxfill=',maxfill
      nnza=icnt/nshg +1
      if(icnt.gt.nnz*nshg) then
         write(*,*) 'increase nnz in genmat to',nnza
         stop
      else
         write(*,*) 'nnz ok  nnz=',nnz,' actually needed',nnza   
         write(*,*) myrank,' is my rank and my nnz_tot is: ',icnt   
      endif
      return
      end









