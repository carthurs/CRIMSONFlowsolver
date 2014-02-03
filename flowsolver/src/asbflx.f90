      subroutine AsBFlx (u,             y,           ac, &   
                         x,             xdist,       xdnv, &
                         shpb,          shglb, &
                         ienb,          iBCB, & 
                         BCB,           invflx,      flxres, &
                         flxLHS,        flxnrm,      xKebe, &
                         SWB)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  boundary elements.
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
        !use turbSA ! access to d2wall
        use LagrangeMultipliers 
!
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndofl), &      
                  x(numnp,nsd), &
                  xdist(nshg), &
                  xdnv(nshg,nsd), &
                  ac(nshg,ndofl),          u(nshg,nsd), &
                  shpb(nshl,ngaussb), &
                  shglb(nsd,nshl,ngaussb), &
                  ienb(npro,nshl), &
                  iBCB(npro,ndiBCB),       BCB(npro,nshlb,ndBCB), &
                  invflx(nshg),            flxres(nshg,nflow), &
                  flxLHS(nshg,1),          flxnrm(nshg,nsd), &
                  SWB(npro,nProps)
!
        dimension yl(npro,nshl,ndofl),     xlb(npro,nenl,nsd), &
                  rl(npro,nshl,nflow),     sgn(npro,nshl), &
                  flhsl(npro,nshl,1),      fnrml(npro,nshl,nsd), &
                  lnflx(npro),             lnode(27), &
                  ul(npro,nshl,nsd),       acl(npro,nshl,ndofl), &
                  xdistl(npro,nshl),       xdnvl(npro,nshl,nsd)
     
        real*8 dwl(npro,nshl)
        
        dimension xKebe(npro,9,nshl,nshl) 

!
!.... compute the nodes which lie on the boundary (hierarchic)
!
        call getbnodes(lnode)
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
        
!        if(iRANS.eq.-2) then
!           call local(d2wall, dwl, ienb, 1, 'gather  ')
!        endif

        rl    = zero
        flhsl = zero
        fnrml = zero
!
        ires = 2
!
!..... to calculate inner product for Lagrange Multipliers
!
        if(Lagrange.gt.zero) then
           allocate(loclhsLag(npro,9,nshlb,nshlb,3))
           loclhsLag = zero
        endif           
!        
        call e3b  (ul,      yl,      acl, &     
                   iBCB,    BCB, &
                   shpb,    shglb, &
                   xlb,     xdistl,  xdnvl, &
                   rl,      sgn,     dwl,     xKebe, &
                   SWB)
     
        ires = 1
!
!.... assemble the residuals
!
        call local (flxres, rl,     ienb,   nflow,  'scatter ')
!
!.... compute the LHS for the flux computation (should only be done
!     once)
!
        call f3lhs (shpb,       shglb,      xlb, &
                    flhsl,      fnrml,      sgn )

!     
!.... reset the non-contributing element values
!
        lnflx = 0
        do n = 1, nshlb
          lnflx = lnflx + min(1, invflx(ienb(:,lnode(n))))
        enddo
!
        do n = 1, nshl
          where (lnflx .ne. nshlb)   flhsl(:,n,1) = zero
          do i = 1, nsd
            where (lnflx .ne. nshlb) fnrml(:,n,i) = zero
          enddo
        enddo
!
!.... assemble the boundary LHS and normal
!
        call local (flxLHS, flhsl,  ienb,   1,      'scatter ')
        call local (flxnrm, fnrml,  ienb,   nsd,    'scatter ')
!
        if(Lagrange.gt.zero) then
           deallocate(loclhsLag)
        endif
!     
!.... end
!
        return
        end
