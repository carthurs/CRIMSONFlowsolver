        subroutine AsIq (y,       x,       shp, &
                         shgl,    ien,     xmudmi, &
                         qres,    rmass    )
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
! interior elements for the global reconstruction of the diffusive
! flux vector.
!
! input:
!     y     (numnp,ndof)        : Y variables
!     x     (numnp,nsd)         : nodal coordinates
!     shp   (nen,nintg)         : element shape-functions
!     shgl  (nsd,nen,nintg)     : element local shape-function gradients
!     ien   (npro)              : nodal connectivity array
!
! output:
!     qres  (numnp,nsd,nsd)  : residual vector for diffusive flux
!     rmass  (numnp)            : lumped mass matrix
!
!----------------------------------------------------------------------
!
        use turbsa      ! access to d2wall
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndof),               x(numnp,nsd), &           
                  shp(nshl,ngauss),         shgl(nsd,nshl,ngauss), &
                  ien(npro,nshl),      dwl(npro,nenl), &
                  qres(nshg,idflx),    rmass(nshg)
!
        dimension yl(npro,nshl,ndof),          xl(npro,nenl,nsd), &
                  ql(npro,nshl,idflx),  rmassl(npro,nshl), &
                  xmudmi(npro,ngauss)
!
        dimension sgn(npro,nshl)
!
!.... create the matrix of mode signs for the hierarchic basis 
!     functions. 
!
        do i=1,nshl
           where ( ien(:,i) < 0 )
              sgn(:,i) = -one
           elsewhere
              sgn(:,i) = one
           endwhere
        enddo

!
!.... gather the variables
!

        call localy(y,      yl,     ien,    ndof,   'gather  ')
        call localx (x,      xl,     ien,    nsd,    'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
!
!.... get the element residuals 
!
        ql     = zero
        rmassl = zero

        call e3q  (yl,         dwl,      shp,      shgl, &
                   xl,         ql,       rmassl, &
                   xmudmi,     sgn  )

!
!.... assemble the diffusive flux residual 
!
        call local (qres,   ql,  ien,  idflx,  'scatter ')
        call local (rmass,  rmassl, ien,  1,          'scatter ')
!
!.... end
!
        return
        end


!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
! interior elements for the global reconstruction of the diffusive
! flux vector.
!
!----------------------------------------------------------------------
        subroutine AsIqSclr (y,       x,       shp, &
                             shgl,    ien,     qres, &   
                             rmass    )
!
        use turbsa      ! access to d2wall
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndof),             x(numnp,nsd),            &
                  shp(nshl,ngauss),         shgl(nsd,nshl,ngauss),   &
                  ien(npro,nshl),      dwl(npro,nenl),               &
                  qres(nshg,nsd),           rmass(nshg)
!
        dimension yl(npro,nshl,ndof),       xl(npro,nenl,nsd),       &
                  ql(npro,nshl,nsd),        rmassl(npro,nshl)
!
        dimension sgn(npro,nshl)

        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
!
!.... gather the variables
!
        call localy(y,      yl,     ien,    ndof,   'gather  ')
        call localx (x,      xl,     ien,    nsd,    'gather  ')
        if (iRANS .eq. -2) then ! kay-epsilon
           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
        endif
!
!.... get the element residuals 
!
        ql     = zero
        rmassl = zero

        call e3qSclr  (yl,      dwl,    shp,    shgl,    &
                       xl,      ql,     rmassl,          &
                       sgn             )

!
!.... assemble the temperature diffusive flux residual 
!
        call local (qres,   ql,  ien,  nsd,  'scatter ')
        call local (rmass,  rmassl, ien,  1, 'scatter ')
!
!.... end
!
        return
        end

