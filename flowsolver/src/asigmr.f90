        subroutine AsIGMR (y,       ac,      x,       xmudmi, &
                           shp,     shgl,    ien,             &
                           res,     qres,                     &
                           xKebe,   xGoC,    rerr, CFLworst)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  interior elements.
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
      use stats
      !use rlssave  ! Use the resolved Leonard stresses at the nodes.
      use timedata    ! time series
      !use turbsa                ! access to d2wall
      use LagrangeMultipliers 


      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndofl),              ac(nshg,ndofl), &
                  x(numnp,nsd),                               &
                  shp(nshl,ngauss),            shgl(nsd,nshl,ngauss), &
                  ien(npro,nshl), &
                  res(nshg,nflow), &
                  qres(nshg,idflx)

!
        dimension yl(npro,nshl,ndofl),         acl(npro,nshl,ndofl), &
                  xl(npro,nenl,nsd),           dwl(npro,nenl),       &
                  rl(npro,nshl,nflow), &
                  ql(npro,nshl,idflx)
!        
        dimension xKebe(npro,9,nshl,nshl), &
                  xGoC(npro,4,nshl,nshl)
!
        dimension rlsl(npro,nshl,6) 

!
        real*8    lStsVec(npro,nshl,nResDims)
        
        dimension xmudmi(npro,ngauss)
        dimension sgn(npro,nshl)
        dimension CFLworst(npro)
!
        real*8 rerrl(npro,nshl,6), rerr(nshg,10)
!
!.... gather the variables
!
!
!.... get the matrix of mode signs for the hierarchic basis functions. 
!
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        call local (qres,   ql,     ien,    idflx,  'gather  ')
!        if (iRANS .eq. -2) then ! kay-epsilon
!           call localx (d2wall,   dwl,     ien,    1,     'gather  ')
!        endif
 
        !if( (iLES.gt.10).and.(iLES.lt.20)) then  ! bardina
        !   call local (rls, rlsl,     ien,       6, 'gather  ')
        !else
           rlsl = zero
        !endif

!
!.... zero the matrices if they are being recalculated
!
        if (lhs .eq. 1)  then
           xKebe = zero
           xGoC  = zero
        endif   
        if(Lagrange.gt.zero) then
           loclhsLag = zero
        endif 
!
!.... get the element residuals, LHS matrix, and preconditioner
!
        rl     = zero

        if(ierrcalc.eq.1) rerrl = zero

        call e3  (yl,      acl,     dwl,     shp, &
                  shgl,    xl,      rl,           &
                  ql,      xKebe,   xGoC,    xmudmi, & 
                  sgn,     rerrl,  rlsl,     CFLworst)
!
!.... assemble the statistics residual
!
        if ( stsResFlg .eq. 1 ) then
           call e3StsRes ( xl, rl, lStsVec )
           call local( stsVec, lStsVec, ien, nResDims, 'scatter ')
        else
!
!.... assemble the residual
!
           call local (res,    rl,     ien,    nflow,  'scatter ')
           
           if ( ierrcalc .eq. 1 ) then
              call local (rerr, rerrl,  ien, 6, 'scatter ')
           endif
        endif
!
!.... end
!
        if (exts) then
           if ((iter.eq.1).and.(mod(lstep,freq).eq.0)) then
              call timeseries(yl,xl,ien,sgn)
           endif
        endif
        
        return
        end



!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!-----------------------------------------------------------------------
!=======================================================================


        subroutine AsIGMRSclr(y,       ac,      x,    &   
                           shp,     shgl,    ien,     &
                           res,     qres,    xSebe, xmudmi )
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  interior elements.
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
      !use     turbSA
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndofl),              ac(nshg,ndofl), &
                  x(numnp,nsd),                               &
                  shp(nshl,ngauss),            shgl(nsd,nshl,ngauss), &
                  ien(npro,nshl),                             &
                  res(nshg),                  qres(nshg,nsd)  

!
        real*8    yl(npro,nshl,ndofl),        acl(npro,nshl,ndofl), &
                  xl(npro,nenl,nsd),                                &
                  rl(npro,nshl),              ql(npro,nshl,nsd),    &
                  dwl(npro,nenl)            
!        
        real*8    xSebe(npro,nshl,nshl),      xmudmi(npro,ngauss) 
!
!.... gather the variables
!
        real*8 sgn(npro,nshl)
!
!.... get the matrix of mode signs for the hierarchic basis functions. 
!
        if (ipord .gt. 1) then
           call getsgn(ien,sgn)
        endif
        
        call localy(y,      yl,     ien,    ndofl,  'gather  ')
        call localy(ac,    acl,     ien,    ndofl,  'gather  ')
        call localx(x,      xl,     ien,    nsd,    'gather  ')
        if(iRANS.lt. 0) &
        call localx(d2wall, dwl,    ien,    1,      'gather  ')
        call local (qres,   ql,     ien,    nsd,    'gather  ')
!
!.... zero the matrices if they are being recalculated
!
        if (lhs .eq. 1)  then
           xSebe = zero
        endif   
!
!.... get the element residuals, LHS matrix, and preconditioner
!
      rl = zero
      call e3Sclr  (yl,      acl,     shp, &
                    shgl,    xl,      dwl, &
                    rl,      ql,      xSebe,   &
                    sgn, xmudmi)
!
!.... assemble the residual
!
        call local (res,    rl,     ien,    1,  'scatter ')
!
!.... end
!
        return
        end
