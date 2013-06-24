      subroutine AsBMFG (u,          y,          ac,          &
                         x,          xdist,      xdnv,        &
                         shpb,       shglb,                   &
                         ienb,       materb,     iBCB,        BCB, &
                         res,        xKebe,      SWB)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  boundary elements.
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
! Alberto Figueroa, Winter 2004.  CMM-FSI
! Irene Vignon, Spring 2004.
!----------------------------------------------------------------------
!
        use turbSA                ! access to d2wall
        use phcommonvars  
        use deformableWall
        use LagrangeMultipliers 
!
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndofl), &        
                  ac(nshg,ndofl),          u(nshg,nsd), &
                  shpb(nshl,ngaussb), &
                  x(numnp,nsd), &
                  xdist(nshg), &
                  xdnv(nshg,nsd), &
                  shglb(nsd,nshl,ngaussb), &
                  ienb(npro,nshl),         materb(npro), &
                  iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB), &
                  res(nshg,nflow),         dwl(npro,nenl), &
                  SWB(npro,nProps)
!
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd), &
                  rl(npro,nshl,nflow),     sgn(npro,nshl), &
                  ul(npro,nshl,nsd),       acl(npro,nshl,ndofl), &
                  xdistl(npro,nshl),       xdnvl(npro,nshl,nsd)
!
        dimension xKebe(npro,9,nshl,nshl)
     
        integer iblk
     
!
!.... get the matrix of mode signs for the hierarchic basis functions
!
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
!
!.... gather the variables
!
        call localy(y,      yl,     ienb,   ndofl,  'gather  ')
        call localy(ac,     acl,    ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
        call localx(u,      ul,     ienb,   nsd,    'gather  ')
        

        call localx(xdist,     xdistl,    ienb,   1,      'gather  ')
        call localx(xdnv,      xdnvl,     ienb,   nsd,    'gather  ')
        
        if (iRANS.eq.-2) then
           call local(d2wall, dwl, ienb, 1, 'gather  ')
        endif

!
!.... zero the matrices if they are being recalculated
!
        if (lhs .eq. 1)  then
           xKebe = zero
        endif   
        if(Lagrange.gt.zero) then
           loclhsLag = zero
        endif 
!
!.... get the boundary element residuals
!
        rl  = zero
        
        filrhsl = zero
!
!.... 3D
!
        call e3b  (ul,      yl,      acl, &
                   iBCB,    BCB,          &
                   shpb,    shglb,        &
                   xlb,     xdistl,  xdnvl, &     
                   rl,      sgn,     dwl,     xKebe, &
                   SWB)
!
!.... assemble the residual and the modified residual
!
        call local (res,    rl,     ienb,   nflow,  'scatter ')

!     
!.... end
!
        return
      end


!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  boundary elements.
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
! Alberto Figueroa, Winter 2004.  CMM-FSI
! Irene Vignon, Spring 2004.
!----------------------------------------------------------------------
!
      subroutine AsBSclr (y,       x,       shpb,    shglb, &
                         ienb,    materb,  iBCB,    BCB, &
                         res)
        use turbSA ! access to d2wall
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndofl),           x(numnp,nsd), &
                  shpb(nshl,*), &
                  shglb(nsd,nshl,*), &
                  ienb(npro,nshl),         materb(npro), &
                  iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB), &
                  res(nshg)         
!
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd), &
                  rl(npro,nshl),     sgn(npro,nshl)
        real*8 dwl(npro,nshl)
!
!.... get the matrix of mode signs for the hierarchic basis functions
!
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
!
!.... gather the variables
!
        call localy(y,      yl,     ienb,   ndofl,  'gather  ')
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
        if(iRANS.eq.-2) then
           call local(d2wall, dwl, ienb, 1, 'gather  ')
        endif
!
!.... get the boundary element residuals
!
        rl  = zero

        call e3bSclr  (yl,      iBCB,    BCB,     shpb,    shglb, &
                       xlb,     rl,      sgn,     dwl)
!
!.... assemble the residual and the modified residual
!
        call local (res,    rl,     ienb,   1,  'scatter ')
!     
!.... end
!
        return
        end


