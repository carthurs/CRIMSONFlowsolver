      module lhsGkeep

      real*8, allocatable :: lhsG(:)

      end module

!----------------------------------------------------------------------------

      subroutine DlhsGkeep

      use lhsGkeep

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      if (allocated(lhsG)) deallocate ( lhsG )

      return
      end
!----------------------------------------------------------------------------

      subroutine keeplhsG

      use lhsGkeep

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      if (.not. allocated(lhsG)) allocate ( lhsG(nnz*nshg) )

      return

      end
      subroutine cmass (shp, shgl, xl, em)
!       
!----------------------------------------------------------------------
!     
!     This subroutine computes the consistent mass matrices
!     
!     Ken Jansen, Spring 2000
!----------------------------------------------------------------------
!     
!
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      integer ne, na, nb, nodlcla, nodlclb, iel
      dimension  & 
           shp(nshl,MAXQPT),   shgl(nsd,nshl,MAXQPT), & 
           em(npro,nshl,nshl), & 
           xl(npro,nenl,nsd)
!
      dimension shapeVar(npro,nshl),   shdrv(npro,nsd,nshl), & 
                sgn(npro,nshl),     dxidx(npro,nsd,nsd),   & 
                shg(npro,nshl,nsd),  & 
                WdetJ(npro)
!
      em = zero
!     
!.... loop through the integration points
!     
      do intp = 1, ngauss      ! (these are in common.h)
!
!.... get the hierarchic shape functions at this int point
!
         call getshp(shp,         shgl,         sgn,  & 
                     shapeVar,       shdrv,        intp)
!     
!.... calculate the determinant of the jacobian and weight it
!     
         call e3metric( xl, shdrv,dxidx,shg,WdetJ)
!     
         do iel = 1, npro
            do  na  = 1, nshl
               do  nb = 1, nshl
                  shp2 = shapeVar(iel,na) * shapeVar(iel,nb)
                  em(iel,na,nb) = em(iel,na,nb) + shp2*WdetJ(iel)
               enddo
            enddo
         enddo
      enddo
!     
!.... return
!     
      return
      end

      subroutine cmassl (shp, shgl, shpf, shglf, xl, em, Qwtf)
!       
!----------------------------------------------------------------------
!     
!     This subroutine computes the consistent mass matrices
!     
!     Ken Jansen, Spring 2000
!----------------------------------------------------------------------
!     
!
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      integer ne, na, nb, nodlcla, nodlclb, iel
      dimension  & 
           shp(nshl,MAXQPT),   shgl(nsd,nshl,MAXQPT), & 
           shpf(nshl,MAXQPT),  shglf(nsd,nshl,MAXQPT), & 
           em(npro,nshl,nshl), eml(npro,nshl), & 
           xl(npro,nenl,nsd)
!
      dimension shapeVar(npro,nshl),   shdrv(npro,nsd,nshl), & 
                sgn(npro,nshl),     dxidx(npro,nsd,nsd),   & 
                shg(npro,nshl,nsd), Qwtf(ngaussf), & 
                WdetJ(npro)
!
      em = zero
      eml= zero

      if (ifproj.eq.1)then
         nods = nshl
      else
         nods = nenl
      endif

!----------------> Get the lumped mass matrix <-----------------------

!     
!.... loop through the integration points
!     
      do intp = 1, ngaussf      ! (these are in common.h)
!
!.... get the hierarchic shape functions at this int point
!
         call getshp(shpf,         shglf,         sgn,  & 
                     shapeVar,       shdrv,        intp)
!     
!.... calculate the determinant of the jacobian and weight it
!     
         call e3metricf( xl, shdrv,dxidx,shg,WdetJ,Qwtf)
!     
         do i=1,nods !nenl !nshl
            eml(:,i) = eml(:,i) + shapeVar(:,i)*WdetJ(:)
         enddo         

      enddo ! End loop over quad points.

          
!--------------> Get the consistent mass matrix <------------------------


      shapeVar = zero
      shdrv = zero
      dxidx = zero
      WdetJ = zero
      shg   = zero

!     
!.... loop through the integration points
!     
      do intp = 1, ngauss       ! (these are in common.h)

!.... get the hierarchic shape functions at this int point
!
!
!.... for the mass matrix to be consistent shp and shgl must be
!.... evaluated with at least higher quadrature than one-pt. quad. 

         call getshp(shp,         shgl,         sgn,   & 
                     shapeVar,       shdrv,        intp)

!     
!.... calculate the determinant of the jacobian and weight it
!     
         call e3metric( xl, shdrv,dxidx,shg,WdetJ)
!     

         do iel = 1, npro
            do  na  = 1, nods !nenl !nshl
               do  nb = 1, nods !nenl !nshl
                  shp2 = shapeVar(iel,na) * shapeVar(iel,nb)
                  em(iel,na,nb) = em(iel,na,nb) + shp2*WdetJ(iel)
               enddo
            enddo
         enddo

      enddo    ! End loop over quadrature points



!----------> Obtain a mixed (lumped/consistent) mass matrix <------------

!... Different combinations of the lump and mass matrices yield
!... filters of varying widths. In the limiting case were
!... the entire matrix is lumped, we obtain the same filter as 
!... in getdmc.f. Note that in these equivalent ways of
!... filtering one-point quadrature is used for shpf and shgl..


      em = (one-flump)*em

      do iel = 1, npro
         do  na  = 1, nods !nenl !nshl
            em(iel,na,na) = em(iel,na,na)+flump*eml(iel,na)
         enddo
      enddo

!     
!.... return
!     
      return
      end



      subroutine cmasstl (shp, shgl, shpf, shglf, xl, em, Qwtf)
!       
!----------------------------------------------------------------------
!     
!     This subroutine computes the consistent mass matrices
!     
!     Ken Jansen, Spring 2000
!----------------------------------------------------------------------
!     
!
  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      integer ne, na, nb, nodlcla, nodlclb, iel
      dimension  & 
           shp(nshl,MAXQPT),   shgl(nsd,nshl,MAXQPT), & 
           shpf(nshl,MAXQPT),  shglf(nsd,nshl,MAXQPT), & 
           em(npro,nshl,nshl), eml(npro,nshl), & 
           xl(npro,nenl,nsd)
!
      dimension shapeVar(npro,nshl),   shdrv(npro,nsd,nshl), & 
                sgn(npro,nshl),     dxidx(npro,nsd,nsd),   & 
                shg(npro,nshl,nsd), Qwtf(ngaussf), & 
                WdetJ(npro)
!
      em = zero
      eml= zero

!----------------> Get the lumped mass matrix <-----------------------

!     
!.... loop through the integration points
!     
      do intp = 1, ngaussf      ! (these are in common.h)
!
!.... get the hierarchic shape functions at this int point
!
         call getshp(shpf,         shglf,         sgn,  & 
                     shapeVar,       shdrv,        intp)
!     
!.... calculate the determinant of the jacobian and weight it
!     
         call e3metricf( xl, shdrv,dxidx,shg,WdetJ,Qwtf)
!     
         do i=1,nshl
            eml(:,i) = eml(:,i) + shapeVar(:,i)*WdetJ(:)
         enddo         

      enddo ! End loop over quad points.

          
!--------------> Get the consistent mass matrix <------------------------


      shapeVar= zero
      shdrv = zero
      dxidx = zero
      WdetJ = zero
      shg   = zero

!     
!.... loop through the integration points
!     
      do intp = 1, ngauss       ! (these are in common.h)

!.... get the hierarchic shape functions at this int point
!
!
!.... for the mass matrix to be consistent shp and shgl must be
!.... evaluated with at least higher quadrature than one-pt. quad. 

         call getshp(shp,         shgl,         sgn,   & 
                     shapeVar,       shdrv,        intp)

!     
!.... calculate the determinant of the jacobian and weight it
!     
         call e3metric( xl, shdrv,dxidx,shg,WdetJ)
!     

         do iel = 1, npro
            do  na  = 1, nshl
               do  nb = 1, nshl
                  shp2 = shapeVar(iel,na) * shapeVar(iel,nb)
                  em(iel,na,nb) = em(iel,na,nb) + shp2*WdetJ(iel)
               enddo
            enddo
         enddo

      enddo    ! End loop over quadrature points



!----------> Obtain a mixed (lumped/consistent) mass matrix <------------

!... Different combinations of the lump and mass matrices yield
!... filters of varying widths. In the limiting case were
!... the entire matrix is lumped, we obtain the same filter as 
!... in getdmc.f. Note that in these equivalent ways of
!... filtering one-point quadrature is used for shpf and shgl..


      do iel = 1, npro
         do  na  = 1, nshl
            em(iel,na,na) = eml(iel,na)
         enddo
      enddo

!     
!.... return
!     
      return
      end

!-----------------------------------------------------------------------
!
!  compute the metrics of the mapping from global to local 
!  coordinates and the jacobian of the mapping (weighted by 
!  the quadrature weight
!
!-----------------------------------------------------------------------
      subroutine e3metricf(  xl,      shgl,     dxidx, & 
                            shg,     WdetJ, Qwtf)

  use phcommonvars  
  IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8     xl(npro,nenl,nsd),    shgl(npro,nsd,nshl), & 
                 dxidx(npro,nsd,nsd),  shg(npro,nshl,nsd),  & 
                 WdetJ(npro),          Qwtf(ngaussf)

      real*8     dxdxi(npro,nsd,nsd),  tmp(npro)

!
!.... compute the deformation gradient
!
      dxdxi = zero
!
       do n = 1, nenl
          dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(:,1,n)
          dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(:,2,n)
          dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(:,3,n)
          dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(:,1,n)
          dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(:,2,n)
          dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(:,3,n)
          dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(:,1,n)
          dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(:,2,n)
          dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(:,3,n)
       enddo
!
!.... compute the inverse of deformation gradient
!
       dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)  & 
                      - dxdxi(:,3,2) * dxdxi(:,2,3)
       dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)  & 
                      - dxdxi(:,1,2) * dxdxi(:,3,3)
       dxidx(:,1,3) =  dxdxi(:,1,2) * dxdxi(:,2,3)  & 
                      - dxdxi(:,1,3) * dxdxi(:,2,2)
       tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)  & 
                             + dxidx(:,1,2) * dxdxi(:,2,1)   & 
                             + dxidx(:,1,3) * dxdxi(:,3,1) )
       dxidx(:,1,1) = dxidx(:,1,1) * tmp
       dxidx(:,1,2) = dxidx(:,1,2) * tmp
       dxidx(:,1,3) = dxidx(:,1,3) * tmp
       dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)  & 
                      - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
       dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)  & 
                      - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
       dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)  & 
                      - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
       dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)  & 
                      - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
       dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)  & 
                      - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
       dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)  & 
                      - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
!
!       WdetJ = Qwt(lcsyst,intp) / tmp

       WdetJ = Qwtf(intp) / tmp
!
!.... compute the global gradient of shape-functions
!
       do n = 1, nshl
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) +  & 
                       shgl(:,2,n) * dxidx(:,2,1) + & 
                       shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) +  & 
                       shgl(:,2,n) * dxidx(:,2,2) + & 
                       shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) +  & 
                       shgl(:,2,n) * dxidx(:,2,3) + & 
                       shgl(:,3,n) * dxidx(:,3,3) 
       enddo

       return
       end



