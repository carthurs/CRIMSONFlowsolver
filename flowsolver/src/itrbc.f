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
        subroutine itrBC (y,ac, iBC, BC, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This program satisfies the boundary conditions on the Y-variables.
c
c input:
c  y      (nshg,nflow)   : y variables 
c  iBC    (nshg)        : Boundary Condition Code
c  BC     (nshg,ndofBC) : boundary condition constraint parameters
c
c output:
c  y      (nshg,nflow)   : Adjusted V value(s) corresponding to a 
c                           constraint d.o.f.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,nflow),             iBC(nshg),
     &            ac(nshg,nflow),            BC(nshg,ndofBC)

        dimension ilwork(nlwork),           iper(nshg)
c
c  limiting...ugly but sometimes the only way
c
        do i=1,nflow
           if(ylimit(1,i).gt.0) 
     &          y(:,i)=min(ylimit(3,i),max(ylimit(2,i),y(:,i))) 
        enddo
c
c.... ------------------------->  Velocity  <--------------------------
c.... 3D
c
c.... x1-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 1)
            y(:,1) =  BC(:,3)  - BC(:,4) * y(:,2)
     &                         - BC(:,5) * y(:,3)
          endwhere
c
c.... x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 2)
            y(:,2) = BC(:,3)  - BC(:,4) * y(:,1)
     &                        - BC(:,5) * y(:,3)
          endwhere
c
c.... x1-velocity and x2-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 3)
            y(:,1) =  BC(:,3)  - BC(:,4) * y(:,3)
            y(:,2) =  BC(:,5)  - BC(:,6) * y(:,3)
          endwhere
c
c.... x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 4)
            y(:,3) = BC(:,3) - BC(:,4) * y(:,1)
     &                       - BC(:,5) * y(:,2)
          endwhere
c
c.... x1-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 5)
            y(:,1) = BC(:,3) - BC(:,4) * y(:,2)
            y(:,3) = BC(:,5) - BC(:,6) * y(:,2)
          endwhere
c
c.... x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 6)
            y(:,2) = BC(:,3)  - BC(:,4) * y(:,1)
            y(:,3) = BC(:,5)  - BC(:,6) * y(:,1)
          endwhere
c
c.... x1-velocity, x2-velocity and x3-velocity, 3D
c
          where (ibits(iBC,3,3) .eq. 7)
            y(:,1) =  BC(:,3)
            y(:,2) =  BC(:,4)
            y(:,3) =  BC(:,5) 
          endwhere
c
c       endif
c
c.... end of velocity
c
c.... ------------------------->  Pressure  <--------------------------
c
        if (any(btest(iBC,2))) then
c
c.... pressure
c
          where (btest(iBC,2))
            y(:,4) = BC(:,1)  ! pressure here
          endwhere
c
        endif
c
c.... local periodic (and axisymmetric) boundary conditions (no communications)
c 
	do i = 1,nflow
           y(:,i) = y(iper(:),i)
           ac(:,i) = ac(iper(:),i)
	enddo
c
c.... communications
c 
        if (numpe > 1) then
           call commu (y, ilwork, nflow, 'out')
           call commu (ac, ilwork, nflow, 'out')
        endif
c
c       slave has masters value, for abc we need to rotate it
c
        if(iabc==1) then        !are there any axisym bc's
           call rotabc(y, iBC, 'out')
           call rotabc(ac, iBC, 'out')
        endif
     
c
c.... return
c
        return
        end


        subroutine itrBCSclr (y, ac, iBC, BC, iper, ilwork)
c
c----------------------------------------------------------------------
c
c This routine satisfies the boundary conditions on the isclr
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension y(nshg,ndof),             iBC(nshg),
     &            ac(nshg,ndof),            BC(nshg,ndofBC)

        dimension ilwork(nlwork),            iper(nshg)
        dimension T(nshg)

        if(isclr.eq.0) then ! this is temperature
           ib=1
           ibb=2
           id=5
        else
           ib=5+isclr
           ibb=ib+1
           id=ib
        endif
c
c  limiting...ugly but sometimes the only way
c
           if(ylimit(1,id).gt.0) 
     &          y(:,id)=min(ylimit(3,id),max(ylimit(2,id),y(:,id))) 
c
c
c.... ------------------------>  Scalar  <------------------------
c
c
        where (btest(iBC,ib))
          y(:,id) =  BC(:,ibb)
        endwhere
c
c.... local periodic (and axisymmetric) boundary conditions (no communications)
c 
	do i = 1,nshg
          y(i,id) = y(iper(i),id)
          ac(i,id) = ac(iper(i),id)
	enddo
c
c.... communications
c 
        if (numpe > 1) then
           T=y(:,id)
           call commu (T, ilwork, 1, 'out')
           y(:,id)=T
           T=ac(:,id)
           call commu (T, ilwork, 1, 'out')
           ac(:,id)=T
        endif
     
        return
        end

