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
      subroutine f3lhs (shpb,   shglb,  xlb,    flhsl,
     &                  fnrml,  sgn )
c
c----------------------------------------------------------------------
c
c  This subroutine computes the element LHS matrix and the normal 
c to the boundary for computation of (output) boundary fluxes. 
c 
c input:
c  shpb   (nen,nintg)           : boundary element shape-functions
c  shglb  (nsd,nen,nintg)       : boundary element grad-shape-functions
c  wghtb  (nintg)               : boundary element weight
c  xlb    (npro,nenl,nsd)       : nodal coordinates
c  sgn    (npro,nshl)           : mode signs for hierarchic basis
c
c output:
c  flhsl  (npro,nenl,1)         : element lumped lhs on flux boundary
c  fnrml  (npro,nenl,nsd)       : RHS of LS projection of normal to 
c                                  flux boundary
c
c
c Note: Special lumping technique is used to compute the LHS. 
c       See T.J.R. Hughes, "The Finite Element Method: Linear 
c       Static and Dynamic Finite Element Analysis", page 445.  
c
c Note: Least-squares projection is used to compute the normal to
c       the boundary at the nodes.  This routine provides the element
c       contribution to the RHS of the projection linear system.
c
c----------------------------------------------------------------------
c
      include "common.h"
c
      dimension shpb(nshl,ngaussb),        shglb(nsd,nshl,ngaussb),
     &          xlb(npro,nenl,nsd),
     &          flhsl(npro,nshl,1),        fnrml(npro,nshl,nsd)
c
      dimension WdetJb(npro),
     &          bnorm(npro,nsd),           fmstot(npro),
     &          temp(npro),                temp1(npro),
     &          temp2(npro),               temp3(npro)

      dimension sgn(npro,nshl),            shape(npro,nshl),
     &          shdrv(npro,nsd,nshl),
     &          v1(npro,nsd),              v2(npro,nsd)
c
c.... integrate the lumped LHS matrix and normal
c
      fmstot = zero
c
c
c.... compute the normal to the boundary 
c

      v1 = xlb(:,2,:) - xlb(:,1,:)
      v2 = xlb(:,3,:) - xlb(:,1,:)
      
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
      
      do intp = 1, ngaussb
c
c.... get the hierarchic shape functions at this int point
c
         call getshp(shpb,        shglb,        sgn, 
     &               shape,       shdrv)
c
         WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
c
c.... compute the lumped LHS and normal
c          
         do n = 1, nenl ! when changed to nshl solution degraded ipord 10
            flhsl(:,n,1) = flhsl(:,n,1) + WdetJb * shape(:,n)

c for curved geometries the below construct for the normals has to be used
            fnrml(:,n,1) = fnrml(:,n,1) + WdetJb * bnorm(:,1)
     &                                           * shape(:,n)
            fnrml(:,n,2) = fnrml(:,n,2) + WdetJb * bnorm(:,2)
     &                                           * shape(:,n)
            fnrml(:,n,3) = fnrml(:,n,3) + WdetJb * bnorm(:,3)
     &                                           * shape(:,n)
          enddo
c
c  To best represent this case it should be assigned to the vertex 
c  modes and higher entities should get zero as is done below
c
          fmstot = fmstot + WdetJb
c
        enddo
        
c$$$        do i=1,nenl
c$$$           fnrml(:,i,:)=bnorm(:,:)
c$$$        enddo
        if(ipord.gt.1)  fnrml(:,nenl:nshl,:)=zero
c
c.... scale the LHS matrix contribution
c
        temp = zero
        do n = 1, nshl
           temp = temp + flhsl(:,n,1)
        enddo
c
        do n = 1, nshl
           flhsl(:,n,1) = flhsl(:,n,1) * fmstot / temp
        enddo
c
c.... return
c
        return
        end
