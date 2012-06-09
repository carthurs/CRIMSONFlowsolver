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
        subroutine e3b (ul,      yl,      acl,     
     &                  iBCB,    BCB,     
     &                  shpb,    shglb,
     &                  xlb,     xdistl,  xdnvl,     
     &                  rl,      sgn,     dwl,     xKebe,
     &                  SWB,     TWB,     EWB,
     &                  PS_global,        
     &                  Kwall_xKebe)
c
c----------------------------------------------------------------------
c
c   This routine calculates the 3D RHS residual of the fluid boundary 
c   elements.
c
c input:
c  yl     (npro,nshl,ndof)      : Y variables
c  iBCB   (npro,ndiBCB)         : boundary condition code (iBCB(:,1) is
c      a bit tested boundary integral flag i.e.
c                  if set to value of BCB      if set to floating value
c      iBCB(:,1) : convective flux * 1            0  (ditto to all below)
c                  pressure   flux * 2
c                  viscous    flux * 4
c                  heat       flux * 8
c                  turbulence wall * 16
c                  scalarI   flux  * 16*2^I 
c                  (where I is the scalar number)
c
c      iBCB(:,2) is the srfID given by the user in MGI that we will
c                collect integrated fluxes for.
c
c  BCB    (npro,nshlb,ndBCB)    : Boundary Condition values
c                                  BCB (1) : mass flux 
c                                  BCB (2) : pressure 
c                                  BCB (3) : viscous flux in x1-direc.
c                                  BCB (4) : viscous flux in x2-direc.
c                                  BCB (5) : viscous flux in x3-direc.
c                                  BCB (6) : heat flux
c  shpb   (nen,ngaussb)           : boundary element shape-functions
c  shglb  (nsd,nen,ngaussb)       : boundary element grad-shape-functions
c  xlb    (npro,nenl,nsd)       : nodal coordinates at current step
c
c output:
c  rl     (npro,nshl,nflow)      : element residual
c
c Note: Always the first side of the element is on the boundary.  
c       However, note that for higher-order elements the nodes on 
c       the boundary side are not the first nshlb nodes, see the 
c       array mnodeb.
c
c----------------------------------------------------------------------
c
        use incpBC
        use LagrangeMultipliers 
c
        include "common.h"
c
        dimension yl(npro,nshl,ndof),          iBCB(npro,ndiBCB),
     &            BCB(npro,nshlb,ndBCB),       shpb(nshl,ngaussb),
     &            shglb(nsd,nshl,ngaussb),     xlb(npro,nenl,nsd),
     &            ul(npro,nshl,nsd),           acl(npro,nshl,ndof),
     &            rl(npro,nshl,nflow),         
     &            xdistl(npro,nshl),           xdnvl(npro,nshl,nsd),
     &            SWB(npro,nProps),            TWB(npro,2),
     &            EWB(npro,1)
c
        dimension g1yi(npro,ndof),             g2yi(npro,ndof),
     &            g3yi(npro,ndof),             WdetJb(npro),
     &            bnorm(npro,nsd)
c
        dimension u1(npro),                    u2(npro),
     &            u3(npro),                    rho(npro),
     &            unm(npro),                   pres(npro),
     &            vdot(npro,nsd),              rlKwall(npro,nshlb,nsd),
     &            usup1(npro,nsd),             usup2(npro,nsd), 
     &            velsup(npro,nsd)
c
        dimension rou(npro),                   rmu(npro)
c
        dimension tau1n(npro),
     &            tau2n(npro),                 tau3n(npro)
c
        dimension lnode(27),                   sgn(npro,nshl),
     &            shape(npro,nshl),            shdrv(npro,nsd,nshl),
     &            rNa(npro,4)

        real*8    xmudmi(npro,ngauss),         dwl(npro,nshl)
c
      	dimension xKebe(npro,9,nshl,nshl), 
     &            rKwall_glob(npro,9,nshl,nshl)
                 
        real*8    Kwall_global(npro,9,9),  
     &            PS_global(npro,9),        
     &            Kwall_xKebe(npro,9,nshl,nshl), 
     &            rsumall(npro,9,nshl,nshl)
     
      	integer   intp, num
      	integer   i,j,k,l,m,n

c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
                
c
c.... loop through the integration points
c
        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif
        
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shpb,        shglb,        sgn, 
     &              shape,       shdrv)
c
c     NOTE I DID NOT PASS THE lnode down.  It is not needed
c     since the shape functions are zero on the boundary
c
c     Note that xmudmi is not calculated at these quadrature 
c     points so you give it a zero.  This has implications.
c     the traction calculated by this approach will include 
c     molecular stresses ONLY.  This is why we will use the 
c     consistent flux method to obtain the forces when doing
c     effective viscosity wall modeling.  When doing slip velocity
c     this is not a problem since the traction is given from the
c     log law relation (not the viscosity).
c
        xmudmi=zero
c
c.... get necessary fluid properties (including eddy viscosity)
c
        call getdiff(dwl, yl,     shape,     xmudmi, xlb,   rmu, rho)
c
c.... calculate the integraton variables
c
        call e3bvar (yl,              acl,             ul,              
     &               shape,
     &               shdrv,           
     &               xlb,             xdistl,          xdnvl,
     &               lnode,           SWB,             TWB,
     &               EWB,
     &               WdetJb,
     &               bnorm,           pres,            
     &               u1,              u2,              u3,
     &               rmu,             unm,
     &               tau1n,           tau2n,           tau3n,
     &               vdot,            
     &               usup1,           usup2,            
     &               velsup,
     &               rlKwall,         
     &               xKebe,           !rKwall_glob,
     &               Kwall_global,    PS_global)
     

c
c....  a hack to only add the deformable wall contribution only once
c....  instead of every iteration of the integration loop
c
        if (ideformwall.eq.1 .and. intp.eq.ngaussb)  then 
c
c....  Form the residual contribution
c 
           rlKwall = zero
           
           do i = 1, nshlb
              do j = 1, nsd
                 do k = 1, nshlb
                    do m = 1, nsd
                       rlKwall(:,i,j) = rlKwall(:,i,j)
     &                 +Kwall_xKebe(:,3*(j-1)+m,i,k)*ul(:,k,m)
                    end do
                 end do
              end do
           enddo
c
c....  Add the PreStress Contribution to the residual
c           
           do i = 1, nshlb
              do j = 1, nsd
                 rlKwall(:,i,j) = rlKwall(:,i,j) 
     &              +PS_global(:,3*(i-1)+j)
              end do
           end do

        end if
        
c        
c.... -----------------> boundary conditions <-------------------
c
        do iel = 1, npro
c
c  if we have a nonzero value then
c  calculate the fluxes through this surface 
c
           iface = abs(iBCB(iel,2))
           if (iface .ne. 0 .and. ires.ne.2) then
              flxID(1,iface) =  flxID(1,iface) + WdetJb(iel)! measure area too
              flxID(2,iface) =  flxID(2,iface) - WdetJb(iel) * unm(iel)
              flxID(3,iface) = flxID(3,iface)
     &                   - ( tau1n(iel) - bnorm(iel,1)*pres(iel))
     &                   * WdetJb(iel) 
              flxID(4,iface) = flxID(4,iface)
     &                   - ( tau2n(iel) - bnorm(iel,2)*pres(iel))
     &                   * WdetJb(iel) 
              flxID(5,iface) = flxID(5,iface)
     &                   - ( tau3n(iel) - bnorm(iel,3)*pres(iel))
     &                   * WdetJb(iel) 

           endif
           
c
c.... mass flux
c

           if (btest(iBCB(iel,1),0)) then
              unm(iel)  = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 unm(iel) = unm(iel) 
     &                    + shape(iel,nodlcl) * BCB(iel,n,1)
              enddo
           endif
c
c.... pressure
c
           if (incp .gt. zero) then
              do k=1, numINCPSrfs
                 if(iBCB(iel,2).eq.inactive(k)) then
c .... do nothing in iBCB
                 else        
                    if (btest(iBCB(iel,1),1)) then
                       pres(iel) = zero
                       do n = 1, nshlb
                          nodlcl = lnode(n)
                          pres(iel) = pres(iel) 
     &                       + shape(iel,nodlcl) * BCB(iel,n,2)
                       enddo
                    endif
                 endif
              enddo
        
           else        
              if (btest(iBCB(iel,1),1)) then
                 pres(iel) = zero
                 do n = 1, nshlb
                    nodlcl = lnode(n)
                    pres(iel) = pres(iel) 
     &                        + shape(iel,nodlcl) * BCB(iel,n,2)
                 enddo
              endif
           endif
c
c.... viscous flux
c        
           if (btest(iBCB(iel,1),2)) then
              tau1n(iel) = zero
              tau2n(iel) = zero
              tau3n(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 tau1n(iel) = tau1n(iel) 
     &                      + shape(iel,nodlcl)*BCB(iel,n,3)
                 tau2n(iel) = tau2n(iel) 
     &                      + shape(iel,nodlcl)*BCB(iel,n,4)
                 tau3n(iel) = tau3n(iel) 
     &                      + shape(iel,nodlcl)*BCB(iel,n,5)
              enddo
           endif
c
c.... turbulence wall (as a way of checking for deformable wall stiffness)
c
           if (btest(iBCB(iel,1),4)) then
              
c              rlKwall(iel,:,:) = rlKwall(iel,:,:) / ngaussb ! divide by number of gauss points 
              pres(iel) = zero                              ! to avoid the gauss point loop
              tau1n(iel) = zero                             ! and make the traction contribution
              tau2n(iel) = zero                             ! zero
              tau3n(iel) = zero                              
           else
              
              if (intp.eq.ngaussb) then
                 rlKwall(iel,:,:) = zero                    ! this is not a deformable element
              end if

              vdot(iel,:) = zero                             
              usup1(iel,:) = zero
              usup2(iel,:) = zero
              velsup(iel,:) = zero
              
              ! we do this because the mass matrix contribution is already in xKebe via e3bvar
              xKebe(iel,:,:,:) = zero  

           endif
c
c..... to calculate inner products for Lagrange Multipliers
c
           if(Lagrange.gt.zero) then
              do k=1, numLagrangeSrfs
                 if (iBCB(iel,2).eq.nsrflistLagrange(k)) then
                    do n = 1, nshlb
                       nodlcl = lnode(n)
                       do m=1, nsd
                          do i=1, nshlb
                             nodlcl2 = lnode(i)
                             do l=1, nsd
                                do j=1, 3
                                   num=(m-1)*nsd+l
                                   loclhsLag(iel,num,n,i,j)=
     &                                loclhsLag(iel,num,n,i,j)+
     &                                shape(iel,nodlcl)*WdetJb(iel)*
     &                                shape(iel,nodlcl2)*
     &                                LagInplaneVectors(m,j,k)*
     &                                LagInplaneVectors(l,j,k)
                                enddo 
                             enddo
                          enddo
                       enddo
                    enddo
                 endif
              enddo
           endif  ! end of if(Lagrange.gt.zero)
c
        enddo                                               ! end of bc loop
c
c$$$c.... if we are computing the bdry for the consistent
c$$$c     boundary forces, we must not include the surface elements
c$$$c     in the computataion (could be done MUCH more efficiently!)--->
                                                                  !this
                                                                  !comment should read as for the consistent flux calculation rather than boundary forces
c$$$c
        if (ires .eq. 2) then
           do iel = 1, npro 
              if (nsrflist(iBCB(iel,2)) .ne. 0) then
                 unm(iel) = zero
                 tau1n(iel) = zero
                 tau2n(iel) = zero
                 tau3n(iel) = zero
c                 pres(iel) = zero
c
c whatever is zeroed here will beome part of the post-processed surface
c                 "traction force"
c
c uncomment the next lines to get all of the t vector coming from
c                 Alberto's wall motion model.
                 vdot(iel,:)=zero
                 usup1(iel,:)=zero
                 usup2(iel,:)=zero
                 velsup(iel,:)=zero
                 rlKwall(iel,:,:)=zero

c
c uncomment the next 8 lines to get only the tangential part
c

c                  vn=dot_product(vdot(iel,:),bnorm(iel,:))
c                  vdot(iel,:)=vn*bnorm(iel,:)
c                  walln1=dot_product(rlkwall(iel,1,:),bnorm(iel,:))
c                  walln2=dot_product(rlkwall(iel,2,:),bnorm(iel,:))
c                  walln3=dot_product(rlkwall(iel,3,:),bnorm(iel,:))
c                  rlkwall(iel,1,:)=walln1*bnorm(iel,:)
c                  rlkwall(iel,2,:)=walln2*bnorm(iel,:)
c                  rlkwall(iel,3,:)=walln3*bnorm(iel,:)
              endif
           enddo
        endif
c
c.... assemble the contributions
c
        rNa(:,1) = -tau1n + bnorm(:,1)*pres 
        rNa(:,2) = -tau2n + bnorm(:,2)*pres 
        rNa(:,3) = -tau3n + bnorm(:,3)*pres 
        rNa(:,4) = unm
        
        if (ideformwall.eq.1) then
        
           rNa(:,1) = rNa(:,1) + vdot(:,1)
           rNa(:,2) = rNa(:,2) + vdot(:,2)
           rNa(:,3) = rNa(:,3) + vdot(:,3)

           if (iwalldamp.eq.1) then  
              rNa(:,1) = rNa(:,1) + velsup(:,1)
              rNa(:,2) = rNa(:,2) + velsup(:,2)
              rNa(:,3) = rNa(:,3) + velsup(:,3)
           endif
        
           if (iwallsupp.eq.1) then  
              rNa(:,1) = rNa(:,1) + usup1(:,1)
              rNa(:,2) = rNa(:,2) + usup1(:,2)
              rNa(:,3) = rNa(:,3) + usup1(:,3)
           endif
        
           if (imeasdist.eq.1) then
              rNa(:,1) = rNa(:,1) + usup2(:,1)
              rNa(:,2) = rNa(:,2) + usup2(:,2)
              rNa(:,3) = rNa(:,3) + usup2(:,3)
           endif
           
        endif
        
c
        if(iconvflow.eq.1) then     ! conservative form was integrated
                                    ! by parts and has a convective 
                                    ! boundary integral
c
c.... assemble the contributions
c
           rou=rho*unm
           rNa(:,1) = rNa(:,1) + rou*u1 
           rNa(:,2) = rNa(:,2) + rou*u2
           rNa(:,3) = rNa(:,3) + rou*u3
        endif
        
        rNa(:,1) = rNa(:,1)*WdetJb
        rNa(:,2) = rNa(:,2)*WdetJb
        rNa(:,3) = rNa(:,3)*WdetJb
        rNa(:,4) = rNa(:,4)*WdetJb
        
c
c.... ------------------------->  Residual  <--------------------------
c
c.... add the flux to the residual
c
        do n = 1, nshlb
           nodlcl = lnode(n)

           rl(:,nodlcl,1) = rl(:,nodlcl,1) - shape(:,nodlcl) * rNa(:,1)
           rl(:,nodlcl,2) = rl(:,nodlcl,2) - shape(:,nodlcl) * rNa(:,2)
           rl(:,nodlcl,3) = rl(:,nodlcl,3) - shape(:,nodlcl) * rNa(:,3)
           rl(:,nodlcl,4) = rl(:,nodlcl,4) - shape(:,nodlcl) * rNa(:,4)

        enddo
        
        if(ideformwall.eq.1 .and. intp.eq.ngaussb) then
           rl(:,1,1) = rl(:,1,1) - rlKwall(:,1,1)
           rl(:,1,2) = rl(:,1,2) - rlKwall(:,1,2)
           rl(:,1,3) = rl(:,1,3) - rlKwall(:,1,3)
           
           rl(:,2,1) = rl(:,2,1) - rlKwall(:,2,1)
           rl(:,2,2) = rl(:,2,2) - rlKwall(:,2,2)
           rl(:,2,3) = rl(:,2,3) - rlKwall(:,2,3)
        
           rl(:,3,1) = rl(:,3,1) - rlKwall(:,3,1)
           rl(:,3,2) = rl(:,3,2) - rlKwall(:,3,2)
           rl(:,3,3) = rl(:,3,3) - rlKwall(:,3,3)
        endif 
c
c.... -------------------->  Aerodynamic Forces  <---------------------
c
        if (((ires .ne. 2) .and. (iter .eq. nitr))
     &                     .and. (abs(itwmod).eq.1)) then
c
c.... compute the forces on the body
c
           where (nsrflist(iBCB(:,2)).eq.1)
              tau1n = ( tau1n - bnorm(:,1)*pres)  * WdetJb
              tau2n = ( tau2n - bnorm(:,2)*pres)  * WdetJb
              tau3n = ( tau3n - bnorm(:,3)*pres)  * WdetJb
           elsewhere
              tau1n = zero
              tau2n = zero
              tau3n = zero
           endwhere
c
c Note that the sign has changed from the compressible code to 
c make it consistent with the way the bflux calculates the forces
c Note also that Hflux has moved to e3btemp
c
          Force(1) = Force(1) - sum(tau1n)
          Force(2) = Force(2) - sum(tau2n)
          Force(3) = Force(3) - sum(tau3n)

c
        endif
c
c.... end of integration loop
c
        enddo
        
        if(ideformwall.eq.1) then
c     
c.... -----> Wall Stiffness and Mass matrices for implicit LHS  <-----------
c     
c.... Now we simply have to add the stiffness contribution in rKwall_glob to 
c.... the mass contribution already contained in xKebe

           xKebe = xKebe*iwallmassfactor+
     &             Kwall_xKebe*(iwallstiffactor*
     &                          betai*Delt(itseq)*Delt(itseq)*alfi)

         endif
c$$$        ttim(40) = ttim(40) + tmr()
c
c.... return
c
        return
        end


c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c*********************************************************************
c*********************************************************************


        subroutine e3bSclr (yl,      iBCB,    BCB,     shpb,    shglb,
     &                      xlb,     rl,      sgn,     dwl)
        include "common.h"
c
        dimension yl(npro,nshl,ndof),          iBCB(npro,ndiBCB),
     &            BCB(npro,nshlb,ndBCB),       shpb(nshl,*),
     &            shglb(nsd,nshl,*),           
     &            xlb(npro,nenl,nsd),          
     &            rl(npro,nshl)
c
        real*8    WdetJb(npro),                bnorm(npro,nsd)
c
        dimension lnode(27),                   sgn(npro,nshl),
     &            shape(npro,nshl),            shdrv(npro,nsd,nshl),
     &            rNa(npro),                   flux(npro)
        real*8    dwl(npro,nshl)

c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... loop through the integration points
c
        if(lcsyst.eq.3.or.lcsyst.eq.4) then
           ngaussb = nintb(lcsyst)
        else
           ngaussb = nintb(lcsyst)
        endif
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
        call getshp(shpb,        shglb,        sgn, 
     &              shape,       shdrv)
c
c.... calculate the integraton variables
c
        call e3bvarSclr (yl,          shdrv,   xlb,
     &                   shape,       WdetJb,  bnorm,
     &                   flux,        dwl )
c        
c.... -----------------> boundary conditions <-------------------
c

c
c.... heat or scalar  flux
c     
        if(isclr.eq.0) then 
           iwalljump=0
        else
           iwalljump=1  !turb wall between heat and scalar flux..jump over
        endif
        ib=4+isclr+iwalljump
        ibb=6+isclr
        do iel=1, npro
c
c  if we have a nonzero value then
c  calculate the fluxes through this surface 
c
           if (iBCB(iel,2) .ne. 0 .and. ires.ne.2) then
              iface = abs(iBCB(iel,2))
              flxID(ibb,iface) =  flxID(ibb,iface) 
     &                          - WdetJb(iel) * flux(iel)
           endif

           if (btest(iBCB(iel,1),ib-1)) then
              flux(iel) = zero
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 flux(iel) = flux(iel) 
     &                     + shape(iel,nodlcl) * BCB(iel,n,ibb)
              enddo           
           endif
        enddo
c
c.... assemble the contributions
c
        rNa(:) = -WdetJb * flux
c
c.... ------------------------->  Residual  <--------------------------
c
c.... add the flux to the residual
c
        do n = 1, nshlb
           nodlcl = lnode(n)
 
           rl(:,nodlcl) = rl(:,nodlcl) - shape(:,nodlcl) * rNa(:)
        enddo
c
c.... -------------------->  Aerodynamic Forces  <---------------------
c
        if ((ires .ne. 2) .and. (iter .eq. nitr)) then
c
c.... compute the forces on the body
c
           if(isclr.eq.0)   HFlux    = sum(flux)
c
        endif
c
c.... end of integration loop
c
        enddo

c
c.... return
c
        return
        end
  
