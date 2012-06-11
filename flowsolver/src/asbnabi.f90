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
      subroutine AsBNABI ( x,       shpb,
     &                   ienb,  iBCB)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c----------------------------------------------------------------------
c
      use pvsQbi
        include "common.h"
c
        dimension xlb(npro,nenl,nsd),    bnorm(npro,nsd),
     &            rl(npro,nshl,nsd),     WdetJb(npro)

        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl),
     &            shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

c
        dimension dxdxib(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd)

c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rl  = zero
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... loop through the integration points
c
        ngaussb = nintb(lcsyst)

        
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, 
     &              shpfun,       shdrv)

c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
           if(lcsyst.ne.6) then
              ipt3=3
           else
              ipt3=5
           endif
           v1 = xlb(:,2,:) - xlb(:,1,:)
           v2 = xlb(:,ipt3,:) - xlb(:,1,:)

c
c.....The following are done in order to correct temp1..3  
c     based on the results from compressible code.  This is done only 
c     for wedges, depending on the boundary face.(tri or quad)  
c        
           if (lcsyst .eq. 1) then
              temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
              temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
              temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
           else 
              temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
              temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
              temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
           endif
c     
           temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
           bnorm(:,1) = temp1 * temp
           bnorm(:,2) = temp2 * temp
           bnorm(:,3) = temp3 * temp
c     
      
           if (lcsyst .eq. 3) then
              WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
           elseif (lcsyst .eq. 4) then
              WdetJb     = Qwtb(lcsyst,intp) / temp
           else
              WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
           endif


           if (numDirCalcSrfs.gt.zero) then
              do iel=1,npro
                 do i=1, numDirCalcSrfs
                    if (btest(iBCB(iel,1),1) .or.
     &                 iBCB(iel,2).eq.nsrflistDirCalc(i)) then
                       bnorm(iel,1:3)=bnorm(iel,1:3)*WdetJb(iel)
                    else
                        bnorm(iel,:) = zero  ! we want zeros where we are not integrating
                    endif
                 enddo
              enddo
           else
              do iel=1,npro
                 if (btest(iBCB(iel,1),1)) then 
                    bnorm(iel,1:3)=bnorm(iel,1:3)*WdetJb(iel)
                 else
                    bnorm(iel,:) = zero  ! we want zeros where we are not integrating
                 endif
              enddo
           endif

c
c  Now lets calculate Integral N_(a:e)^i n_i d Gamma
c
c
           do n = 1, nshlb
              nodlcl = lnode(n)
         rl(:,nodlcl,1) = rl(:,nodlcl,1) + shpfun(:,nodlcl) * bnorm(:,1)
         rl(:,nodlcl,2) = rl(:,nodlcl,2) + shpfun(:,nodlcl) * bnorm(:,2)
         rl(:,nodlcl,3) = rl(:,nodlcl,3) + shpfun(:,nodlcl) * bnorm(:,3) 
         
           enddo

        enddo  ! quadrature point loop
c
c.... assemble the NABI vector
c
        call local (NABI,    rl,     ienb,   3,  'scatter ')
c
c     push the surf number which we have associated with boundary
C     elements up to the global level in the array ndsurf
c
        do iel=1,npro
           if (iBCB(iel,2) .ne. 0) then
              iface = iBCB(iel,2)
              do k=1,nshlb
                 if (ndsurf(ienb(iel,k)).ne.1) then
                    ndsurf(ienb(iel,k))=iface   
                 endif
              enddo
           endif
           
           if (numVisFluxSrfs .gt. zero) then
              do k=1, numVisFluxSrfs
                 if (iBCB(iel,2) .eq. nsrflistVisFlux(k)) then
                    iBCB(iel,1) = 6
                 endif
              enddo
           endif
        enddo
c     
c.... end
c
        return
        end


      subroutine AsBNASC ( x,       shpb,
     &                   ienb,  iBCB)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles the data corresponding to the
c  boundary elements.
c
c Irene Vignon - copied from AsBNABI, Fall 2005.  (Fortran 90)
c----------------------------------------------------------------------
c
      use pvsQbi
        include "common.h"
c
        dimension xlb(npro,nenl,nsd),
     &            rl(npro,nshl),     WdetJb(npro)

        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl),
     &            shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

c
        dimension dxdxib(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd)

c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c.... get the boundary element residuals
c
        rl  = zero
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... loop through the integration points
c
        ngaussb = nintb(lcsyst)

        
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, 
     &              shpfun,       shdrv)

c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
           if(lcsyst.ne.6) then
              ipt3=3
           else
              ipt3=5
           endif
           v1 = xlb(:,2,:) - xlb(:,1,:)
           v2 = xlb(:,ipt3,:) - xlb(:,1,:)

c
c.....The following are done in order to correct temp1..3  
c     based on the results from compressible code.  This is done only 
c     for wedges, depending on the boundary face.(tri or quad)  
c        
           if (lcsyst .eq. 1) then
              temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
              temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
              temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
           else 
              temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
              temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
              temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
           endif
c     
           temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
       
           if (lcsyst .eq. 3) then
              WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
           elseif (lcsyst .eq. 4) then
              WdetJb     = Qwtb(lcsyst,intp) / temp
           else
              WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
           endif


c......here I only want the d Gamma, not the n_i
           if (numDirCalcSrfs.gt.zero) then
              do iel=1,npro
                 do i=1, numDirCalcSrfs
                    if (btest(iBCB(iel,1),1) .or.
     &                 iBCB(iel,2).eq.nsrflistDirCalc(i)) then
                    else
                       WdetJb(iel) = zero  ! we want zeros where we are not integrating
                    endif
                 enddo
              enddo
           else
              do iel=1,npro
                 if (btest(iBCB(iel,1),1)) then 
                 else
                    WdetJb(iel) = zero  ! we want zeros where we are not integrating
                 endif
              enddo
           endif

c
c  Now lets calculate Integral N_(a:e) d Gamma
c
c
           do n = 1, nshlb
              nodlcl = lnode(n)
            rl(:,nodlcl) = rl(:,nodlcl) + shpfun(:,nodlcl) * WdetJb(:)
         
           enddo

        enddo  ! quadrature point loop
c
c.... assemble the NASC vector
c
        call local (NASC,    rl,     ienb,   1,  'scatter ')
c
c     push the surf number which we have associated with boundary
C     elements up to the global level in the array ndsurf
c
        do iel=1,npro
           if (iBCB(iel,2) .ne. 0) then
              iface = iBCB(iel,2)
              do k=1,nshlb
                 if (ndsurf(ienb(iel,k)).ne.1) then
                    ndsurf(ienb(iel,k))=iface   
                 endif
              enddo
           endif
        enddo
c     
c.... end
c
        return
        end

      subroutine AsBPNABI ( x,       shpb,
     &                   ienb,  iBCB)
c
c----------------------------------------------------------------------
c
c This routine computes and assembles data required for an Augmented
c Lagrangian Method. 
c
c----------------------------------------------------------------------
c
        use pvsQbi
        use LagrangeMultipliers ! brings in face radius and center 
        include "common.h"
c
        dimension xlb(npro,nenl,nsd),    bnorm(npro,nsd),
     &            rl(npro,nshl,nsd),     WdetJb(npro),
     &            rl2(npro,nshl,nsd,3)

        dimension x(numnp,nsd),
     &            shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb),
     &            ienb(npro,nshl),
     &            iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl),
     &            shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

c
        dimension dxdxib(npro,nsd,nsd),      temp(npro),
     &            temp1(npro),               temp2(npro),
     &            temp3(npro),
     &            v1(npro,nsd),              v2(npro,nsd)
c
        real*8    intpxlb(npro,3,nsd),    intpdistance(npro,3),
     &            intpprofile(npro,3)
        real*8    tmpLagInplaneVectors(3,3,0:MAXSURF)
        real*8    tmpLagProfileArea(0:MAXSURF)
        real*8    tmpProfileDelta(0:MAXSURF)
        real*8    Inplane1, Inplane2, Inplane3, InplaneNorm
        integer   count
c
c.... get the matrix of mode signs for the hierarchic basis functions
c
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
c
c.... gather the variables
c
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
c
c....   calculate quadrature points
c        
        intpxlb = zero
        intpdistance = zero
        intpprofile = zero
        if (lcsyst.ne.6) then 
           do intp=1, 3 ! use 3 quadrature points
              do n=1, 3  
                 intpxlb(:,intp,:) = intpxlb(:,intp,:)
     &              +xlb(:,n,:)*Qptb(1,n,intp)
              enddo
           enddo
        else
           do intp=1, 3 ! use 3 quadrature points
              do n=1, 2  
                 intpxlb(:,intp,:) = intpxlb(:,intp,:)
     &              +xlb(:,n,:)*Qptb(1,n,intp)
              enddo
              intpxlb(:,intp,:) = intpxlb(:,intp,:)
     &           +xlb(:,5,:)*Qptb(1,3,intp)              
           enddo
        endif
c
c....   calculate profile functions at quadrature points
c        
        do k=1, numLagrangeSrfs
           do intp = 1, 3
              do iel=1, npro
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    intpdistance(iel,intp)=sqrt(
     &                (intpxlb(iel,intp,1)-LagCenter(1,k))**2+
     &                (intpxlb(iel,intp,2)-LagCenter(2,k))**2+
     &                (intpxlb(iel,intp,3)-LagCenter(3,k))**2)
     &                /LagRadius(k)
                    intpprofile(iel,intp)=
     &                (ProfileOrder(k)+2)/ProfileOrder(k)*
     &                (1-intpdistance(iel,intp)**ProfileOrder(k))
                 endif
              enddo
           enddo
         enddo  

c
c.... get the boundary element residuals
c
        rl  = zero
        rl2 = zero
        tmpLagProfileArea = zero
        tmpProfileDelta = zero     
        tmpLagInplaneVectors = zero
c
c.... compute the nodes which lie on the boundary (hierarchic)
c
        call getbnodes(lnode)
c
c.... loop through the integration points
c
        ngaussb = nintb(lcsyst)
c        
        do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, 
     &              shpfun,       shdrv)

c
c.... compute the normal to the boundary. This is achieved by taking
c     the cross product of two vectors in the plane of the 2-d 
c     boundary face.
c
           if(lcsyst.ne.6) then
              ipt3=3
           else
              ipt3=5
           endif
           v1 = xlb(:,2,:) - xlb(:,1,:)
           v2 = xlb(:,ipt3,:) - xlb(:,1,:)
c
c.....The following are done in order to correct temp1..3  
c     based on the results from compressible code.  This is done only 
c     for wedges, depending on the boundary face.(tri or quad)  
c        
           if (lcsyst .eq. 1) then
              temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
              temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
              temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
           else 
              temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
              temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
              temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
           endif
c     
           temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
           bnorm(:,1) = temp1 * temp
           bnorm(:,2) = temp2 * temp
           bnorm(:,3) = temp3 * temp
c 
           if (lcsyst .eq. 3) then
              WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
           elseif (lcsyst .eq. 4) then
              WdetJb     = Qwtb(lcsyst,intp) / temp
           else
              WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
           endif
c
           do iel=1,npro
              count = 0
              do k=1, numLagrangeSrfs
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    iface = iBCB(iel,2)
                    do kk=1,nshlb
                       if (ndsurf(ienb(iel,kk)).ne.1) then
                          ndsurf(ienb(iel,kk))=iface   
                       endif
                    enddo
                    count = count+1
                 endif
              enddo
              if (count .eq. 0) then
                 bnorm(iel,:) = zero  ! we want zeros where we are not integrating
                 WdetJb(iel) = zero  ! we want zeros where we are not integrating
              endif              
           enddo
c
c   Calculate two orthonormal in-plane vectors
c   |bnorm(iel,1)  bnorm(iel,2)  bnorm(iel,3) |
c   |v1(iel,1)     v1(iel,2)     v1(iel,3)    | 
c   x1 component: -v1(iel,2)*bnorm(iel,3)+v1(iel,3)*bnorm(iel,2)
c   x2 component: -v1(iel,3)*bnorm(iel,1)+v1(iel,1)*bnorm(iel,3)
c   x3 component: -v1(iel,1)*bnorm(iel,2)+v1(iel,2)*bnorm(iel,1)
c
           do k=1, numLagrangeSrfs
              do iel=1,npro
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    tmpLagInplaneVectors(1:3,1,k)=bnorm(iel,1:3)
                    tmpLagInplaneVectors(1:3,2,k)=v1(iel,1:3)
     &                 /sqrt(v1(iel,1)**2+v1(iel,2)**2+v1(iel,3)**2)
                    Inplane1=-v1(iel,2)*bnorm(iel,3)
     &                 +v1(iel,3)*bnorm(iel,2)
                    Inplane2=-v1(iel,3)*bnorm(iel,1)
     &                 +v1(iel,1)*bnorm(iel,3)
                    Inplane3=-v1(iel,1)*bnorm(iel,2)
     &                 +v1(iel,2)*bnorm(iel,1)
                    InplaneNorm=one
     &                 /sqrt(Inplane1**2+Inplane2**2+Inplane3**2)
                    tmpLagInplaneVectors(1,3,k)=Inplane1*InplaneNorm
                    tmpLagInplaneVectors(2,3,k)=Inplane2*InplaneNorm
                    tmpLagInplaneVectors(3,3,k)=Inplane3*InplaneNorm
                    exit
                 endif
              enddo
           enddo              
c
c  Now lets calculate Integral N_(a:e)^i n_i ProfileFunction  d Gamma
c
c
           do n = 1, nshlb
              nodlcl = lnode(n)
              rl(:,nodlcl,1) = rl(:,nodlcl,1) + shpfun(:,nodlcl) 
     &           * bnorm(:,1)*intpprofile(:,intp)*WdetJb(:)
              rl(:,nodlcl,2) = rl(:,nodlcl,2) + shpfun(:,nodlcl) 
     &           * bnorm(:,2)*intpprofile(:,intp)*WdetJb(:)
              rl(:,nodlcl,3) = rl(:,nodlcl,3) + shpfun(:,nodlcl) 
     &           * bnorm(:,3)*intpprofile(:,intp)*WdetJb(:)
           enddo
c
c  Now lets calculate Integral N_(a:e)^i n_i N_(b:e)^i n_i d Gamma
c
c
           do k=1, numLagrangeSrfs
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 do m=1, nsd                 
                    rl2(:,nodlcl,m,1)=rl2(:,nodlcl,m,1)+
     &                 shpfun(:,nodlcl)*shpfun(:,nodlcl)*WdetJb(:)
     &                 *tmpLagInplaneVectors(m,1,k)
     &                 *tmpLagInplaneVectors(m,1,k)
                    rl2(:,nodlcl,m,2)=rl2(:,nodlcl,m,2)+
     &                 shpfun(:,nodlcl)*shpfun(:,nodlcl)*WdetJb(:)
     &                 *tmpLagInplaneVectors(m,2,k)
     &                 *tmpLagInplaneVectors(m,2,k)
                    rl2(:,nodlcl,m,3)=rl2(:,nodlcl,m,3)+
     &                 shpfun(:,nodlcl)*shpfun(:,nodlcl)*WdetJb(:)
     &                 *tmpLagInplaneVectors(m,3,k)
     &                 *tmpLagInplaneVectors(m,3,k)
                 enddo
              enddo
           enddo
           
           do k=1, numLagrangeSrfs
              do iel=1,npro
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    tmpLagProfileArea(k)=tmpLagProfileArea(k)+
     &                 intpprofile(iel,intp)*WdetJb(iel)
                    tmpProfileDelta(k)=tmpProfileDelta(k)+
     &                 intpprofile(iel,intp)**2*WdetJb(iel)
                 endif
              enddo
           enddo
        enddo  ! quadrature point loop
c
c.... assemble the PNABI vector
c
        call local (PNABI,    rl,     ienb,   3,  'scatter ')
c
c.... assemble the NANBIJ vector
c      
        do i=1, 3
           call local (NANBIJ(:,:,i),rl2(:,:,:,i),ienb,3,'scatter ')
        enddo
        
        do k=1, numLagrangeSrfs
           LagProfileArea(k)=LagProfileArea(k)+tmpLagProfileArea(k)
           ProfileDelta(k)=ProfileDelta(k)+tmpProfileDelta(k)
           InplaneNorm=sqrt(LagInplaneVectors(1,1,k)**2+
     &        LagInplaneVectors(2,1,k)**2+LagInplaneVectors(3,1,k)**2)
           if (InplaneNorm .eq. zero) then
              LagInplaneVectors(:,:,k)=tmpLagInplaneVectors(:,:,k)
           endif
        enddo   
c                 
        return
        end
