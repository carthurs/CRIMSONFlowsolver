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
        subroutine genshp (shp,    shgl, nshp, nblk)  
c
c----------------------------------------------------------------------
c
c This subroutine generates shape functions for triangular,
c quadrilateral, tetrahedron, wedge and brick elements and pyramids.
c
c----------------------------------------------------------------------
c
        include "common.h"
c
        dimension shp(MAXTOP,maxsh,MAXQPT), 
     &            shgl(MAXTOP,nsd,maxsh,MAXQPT)
c
c.... loop through element blocks
c
        maxnint=1
          do iblk = 1, nblk
c
c.... get coord. system and element type 
c
            lcsyst = lcblk(3,iblk)
            nshl   = lcblk(10,iblk)
c
c.... generate the shape-functions in local coordinates
c
            select case ( lcsyst )
            case ( 1 )          ! tets
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)  
               call shpTet(ipord,Qpt(1,1:3,i),shp(1,:,i),shgl(1,:,:,i))
            enddo
            shgl(1,:,1:nshl,1:nint(lcsyst)) = 
     &      shgl(1,:,1:nshl,1:nint(lcsyst))/two
c     
            case ( 2 )          ! hexes
c     
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shphex  (ipord, Qpt(2,1:3,i),shp(2,:,i),
     &                       shgl(2,:,:,i))
            enddo
c
            case ( 3 )          ! wedges
c
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shp6w  (ipord,Qpt(3,1:3,i),shp(3,:,i),
     &                       shgl(3,:,:,i))
            enddo

         case ( 5)              ! pyramids
            
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shppyr (ipord,Qpt(5,1:3,i),shp(5,:,i),shgl(5,:,:,i))
               
            enddo
c
c.... nonexistent element
c
            case default
c
            call error ('genshp  ', 'elem Cat', lelCat)
c
            end select
c
c.... end of generation
c
          enddo
c
c.... return
c
        return
        end
