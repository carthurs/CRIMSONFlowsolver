!------------------------------------------------------------------------
!
! This file contains functions for dealing with higher order shape
! functions at the element level.
!
! Christian Whiting, Winter 1999
!------------------------------------------------------------------------

      subroutine getsgn(ien, sgn)
!------------------------------------------------------------------------
!     returns the matrix of mode signs used for negating higher order
!     basis functions. Connectivity array is assumed to have negative
!     signs on all modes to be negated.
!------------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension ien(npro,nshl),  sgn(npro,nshl)
      
      do i=nenl+1,nshl
         where ( ien(:,i) < 0 )
            sgn(:,i) = -one
         elsewhere
            sgn(:,i) = one
         endwhere
      enddo
      
      return 
      end
      
      subroutine getshp(shp, shgl, sgn, shapeVar, shdrv)
!------------------------------------------------------------------------
!     returns the matrix of element shape functions with the higher
!     order modes correctly negated at the current quadrature point.
!------------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      dimension shp(nshl,ngauss),   shgl(nsd,nshl,ngauss), &
                sgn(npro,nshl),     shapeVar(npro,nshl), &
                shdrv(npro,nsd,nshl)
      
      
      do i=1,nenl
         shapeVar(:,i) = shp(i,intp)
         do j=1,3
            shdrv(:,j,i) = shgl(j,i,intp)
         enddo
      enddo
      if ( ipord > 1 ) then
         do i=nenl+1,nshl
            shapeVar(:,i) = sgn(:,i) * shp(i,intp)
            do j=1,3
               shdrv(:,j,i) = shgl(j,i,intp)*sgn(:,i) 
            enddo
         enddo
      endif
      
      return 
      end
      
      subroutine getshpb(shp, shgl, sgn, shapeVar, shdrv)
!------------------------------------------------------------------------
!     returns the matrix of element shape functions with the higher
!     order modes correctly negated at the current quadrature point.
!------------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension shp(nshl,ngaussb),  shgl(nsd,nshl,ngaussb), &
                sgn(npro,nshl),     shapeVar(npro,nshl), &
                shdrv(npro,nsd,nshl)
      
      
      do i=1,nenl
         shapeVar(:,i) = shp(i,intp)
         do j=1,3
            shdrv(:,j,i) = shgl(j,i,intp)
         enddo
      enddo
      if ( ipord > 1 ) then
         do i=nenl+1,nshl
            shapeVar(:,i) = sgn(:,i) * shp(i,intp)
            do j=1,3
               shdrv(:,j,i) = shgl(j,i,intp)*sgn(:,i) 
            enddo
         enddo
      endif
      
      return 
      end
      
      subroutine getbnodes(lnode)
!------------------------------------------------------------------------
!     compute the higher order modes that lie on the boundary of the 
!     element.
!------------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      dimension lnode(27)

!
!.... boundary triangle of tet element
!
      if (lcsyst .eq. 1) then
         do n = 1, nenbl
            lnode(n) = n
         enddo
         if( ipord>1 ) then
            nem = ipord-1
            do n=1,3*nem
               lnode(nenbl+n) = nenbl+1+n
            enddo
         endif
         if( ipord>2 ) then
            nfm = (ipord-2)*(ipord-1)/2
            do n=1,nfm
               lnode(3+3*nem+n) = 4+6*nem+n
            enddo
         endif
!
!.....boundary quadrilateral for a hex element
!
      else if(lcsyst .eq. 2) then
         do n = 1, nenbl
            lnode(n) = n
         enddo
         if( ipord > 1) then
            nem = ipord -1
            do n=1,4*nem
               lnode(nenbl+n) = 8+n
            enddo
         endif
         
         if( ipord > 3) then
            nfm = (ipord-2)*(ipord-3)/2
            do n=1,nfm
               lnode(4+4*nem+n) = 8+12*nem+n
            enddo
         endif         
!
!
!.... This renumbers the boundary nodes for a wedge element when the
!     boundary element is a quad face.  From lnode = [1 2 3 4]
!                                        to   lnode = [1 4 5 2]
!
      else if(lcsyst .eq. 3) then
         do n = 1, nenbl
            lnode(n) = n
         enddo
!
!     Need to implement for cubic, this valid only for ipord=2
! 
         if( ipord>1 ) then
            nem = ipord-1
            do n=1,3*nem
               lnode(nenbl+n) = 6+n
            enddo
         endif
!
!     Boundary quad of wedge element
!
      else if(lcsyst .eq. 4) then
         lnode(1) = 1
         lnode(2) = 4
         lnode(3) = 5
         lnode(4) = 2
!$$$c     
!$$$c     Need to implement for cubic, this valid only for ipord=2
!$$$c     
!$$$         if( ipord > 1) then
!$$$            lnode(5) = 9
!$$$            lnode(6) = 15
!$$$            lnode(7) = 12
!$$$            lnode(8) = 13
!$$$               nem = ipord -1
!$$$               do n=1,4*nem
!$$$                  lnode(nenbl+n) = 6+n
!$$$               enddo
!$$$         endif         
!     
!     Boundary quad of pyramid element
!
      else if(lcsyst .eq. 5) then
         lnode(1) = 1
         lnode(2) = 2
         lnode(3) = 3
         lnode(4) = 4
!$$$  c
!$$$  c     Need to implement for cubic, this valid only for ipord=2
!$$$  c          
!$$$            if( ipord > 1) then
!$$$               lnode(5) = 9
!$$$               lnode(6) = 15
!$$$               lnode(7) = 12
!$$$               lnode(8) = 13
!$$$               nem = ipord -1
!$$$               do n=1,4*nem
!$$$                  lnode(nenbl+n) = 6+n
!$$$               enddo
!$$$            endif         
!     
!     Boundary triangle of pyramid element
!
      else if(lcsyst .eq. 6) then
         lnode(1) = 1
         lnode(2) = 5
         lnode(3) = 2
!$$$c
!$$$c     Need to implement for cubic, this valid only for ipord=2
!$$$c          
!$$$            if( ipord > 1) then
!$$$               lnode(5) = 9
!$$$               lnode(6) = 15
!$$$               lnode(7) = 12
!$$$               lnode(8) = 13
!$$$               nem = ipord -1
!$$$               do n=1,4*nem
!$$$                  lnode(nenbl+n) = 6+n
!$$$               enddo
!$$$            endif         
!     
!.... other element types need to be implemented
!
      else
         write (*,*) 'Boundary element not implemented for lcyst=' &
              ,lcsyst
         stop
      endif
      
      return 
      end
      
!-----------------------------------------------------------------------
!
!  Evaluate coefficient vector at its interpolation points
!
!-----------------------------------------------------------------------
      subroutine evalAtInterp( ycoeff,  yvals,  x,   nvars, npts )

      use     pointer_data
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer nvars, npts, nHits(nshg)
      
      real*8  ycoeff(nshg,ndof),   yvals(nshg,nvars), &
              shp(nshl,npts),      shgl(nsd,nshl,npts), &
              intpnt(3,npts),      x(numnp,nsd)
      
      real*8, allocatable :: ycl(:,:,:)
      real*8, allocatable :: xl(:,:,:)
      real*8, allocatable :: yvl(:,:,:)
      real*8, allocatable :: sgn(:,:)

      yvals = zero
!
!.... generate the shape functions at the interpolation points
!
      call getIntPnts(intpnt,npts)
      do i=1,npts
         call shpTet(ipord,intpnt(:,i),shp(:,i),shgl(:,:,i))
      enddo
!
!.... loop over element blocks
!
      nHits = 0
      do iblk = 1, nelblk
         iel    = lcblk(1,iblk)
         lcsyst = lcblk(3,iblk)
         nenl   = lcblk(5,iblk) ! no. of vertices per element
         nshl   = lcblk(10,iblk)
         ndofl  = lcblk(8,iblk)
         npro   = lcblk(1,iblk+1) - iel 

         allocate ( ycl(npro,nshl,ndof ) )
         allocate ( yvl(npro,nshl,nvars) )
         allocate ( xl(npro,nenl,nsd   ) )
         allocate ( sgn(npro,nshl)       )
         
         call getsgn(mien(iblk)%p,sgn)
         
         call localy( ycoeff, ycl, mien(iblk)%p, ndof,  'gather  ')
         call localx( x,      xl,  mien(iblk)%p, nsd,   'gather  ')

         call eval  ( xl,       ycl,      yvl,       &
                      shp,      shgl,     sgn,       &
                      nvars,    npts    )

!
!.... average coefficients since stresses may be discontinuous
!         
         call localSum( yvals,    yvl,    mien(iblk)%p,   &
                        nHits,    nVars)  
         
         
         deallocate ( ycl )
         deallocate ( yvl )
         deallocate ( sgn )
         deallocate ( xl  )
!
      enddo

!
!.... average the global values
!
      do i = 1, nshg
         do j = 1, nvars
            yvals(i,j) = yvals(i,j)/nHits(i) !(real(nHits(i),8))
         enddo
      enddo
      
      return
      end

!-----------------------------------------------------------------------
!
!  evaluate in element coordinate system
!
!-----------------------------------------------------------------------
      subroutine eval( xl,      ycl,     yvl,      &
                       shp,     shgl,    sgn, &
                       nvars,   npts ) 
      
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer nvars
!
      real*8  ycl(npro,nshl,ndof),   yvl(npro,nshl,nvars), &
              sgn(npro,nshl),        shapeVar(npro,nshl), &
              shdrv(npro,nsd,nshl),  shp(nshl,npts), &
              shgl(nsd,nshl,npts),   xl(npro,nenl,nsd), &
              shg(npro,nshl,nsd),    gradV(npro,nsd,nsd), &
              dxidx(npro,nsd,nsd),   tmp(npro), wtmp
      
      yvl = zero
!
!.... loop over interpolation points
!
      do intp = 1, npts
         call getshp(shp,          shgl,      sgn,  &
                     shapeVar,     shdrv)
      
!
!.... pressure and velocity
!
         do i = 1, nshl
            do j = 1, 4
               yvl(:,intp,j) = yvl(:,intp,j) + shapeVar(:,i) * ycl(:,i,j)
            enddo
         enddo
!
!.... viscous stress
!
         call e3metric( xl,         shdrv,      dxidx,   &
                        shg,        tmp)

         gradV = zero
         do n = 1, nshl
            gradV(:,1,1) = gradV(:,1,1) + shg(:,n,1) * ycl(:,n,2)
            gradV(:,2,1) = gradV(:,2,1) + shg(:,n,1) * ycl(:,n,3)
            gradV(:,3,1) = gradV(:,3,1) + shg(:,n,1) * ycl(:,n,4)
!     
            gradV(:,1,2) = gradV(:,1,2) + shg(:,n,2) * ycl(:,n,2)
            gradV(:,2,2) = gradV(:,2,2) + shg(:,n,2) * ycl(:,n,3)
            gradV(:,3,2) = gradV(:,3,2) + shg(:,n,2) * ycl(:,n,4)
!     
            gradV(:,1,3) = gradV(:,1,3) + shg(:,n,3) * ycl(:,n,2)
            gradV(:,2,3) = gradV(:,2,3) + shg(:,n,3) * ycl(:,n,3)
            gradV(:,3,3) = gradV(:,3,3) + shg(:,n,3) * ycl(:,n,4)
         enddo

         rmu = datmat(1,2,1)
            
         yvl(:,intp,6 ) = two * rmu * gradV(:,1,1)
         yvl(:,intp,7 ) = two * rmu * gradV(:,2,2)
         yvl(:,intp,8 ) = two * rmu * gradV(:,3,3)

         yvl(:,intp,9 ) = rmu * ( gradV(:,1,2) + gradV(:,2,1) )
         yvl(:,intp,10) = rmu * ( gradV(:,1,3) + gradV(:,3,1) )
         yvl(:,intp,11) = rmu * ( gradV(:,2,3) + gradV(:,3,2) )

!
!.... loop over interpolation points
!         
      enddo
      
      return
      end

         

