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
        subroutine genblk (IBKSZ)
c
c----------------------------------------------------------------------
c
c  This routine reads the interior elements and generates the
c  appropriate blocks.
c
c----------------------------------------------------------------------
c
        use pointer_data
c
        include "common.h"
c
        integer, allocatable :: ientp(:,:)
        integer mater(ibksz)
        integer intfromfile(50) ! integers read from headers
        character*255 fname1
c
        iel=1
        itpblk=nelblk
        nelblk=0
        mattyp = 0
        ndofl = ndof
        nsymdl = nsymdf
        do iblk = 1, itpblk
c
c           read(igeom) neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst
           iseven=7
c           call creadlist(igeom,iseven,
c     &          neltp,nenl,ipordl,nshl, ijunk, ijunk, lcsyst)
           iseven=7
           fname1='connectivity interior?'
           call readheader(igeom,fname1,intfromfile,iseven,
     &                     "integer", iotype)
           neltp  =intfromfile(1)
           nenl   =intfromfile(2)
           ipordl =intfromfile(3)
           nshl   =intfromfile(4)
           ijunk  =intfromfile(5)
           ijunk  =intfromfile(6)
           lcsyst =intfromfile(7)
           allocate (ientp(neltp,nshl))
c           read(igeom) ientp
           iientpsiz=neltp*nshl
           call readdatablock(igeom,fname1,ientp,iientpsiz,
     &                     "integer", iotype)

           do n=1,neltp,ibksz 
              nelblk=nelblk+1
              npro= min(IBKSZ, neltp - n + 1)
c
              lcblk(1,nelblk)  = iel
c              lcblk(2,nelblk)  = iopen ! available for later use
              lcblk(3,nelblk)  = lcsyst
              lcblk(4,nelblk)  = ipordl
              lcblk(5,nelblk)  = nenl
              lcblk(6,nelblk)  = nfacel
              lcblk(7,nelblk)  = mattyp
              lcblk(8,nelblk)  = ndofl
              lcblk(9,nelblk)  = nsymdl 
              lcblk(10,nelblk) = nshl ! # of shape functions per elt
c
c.... allocate memory for stack arrays
c
              allocate (mmat(nelblk)%p(npro))
c
              allocate (mien(nelblk)%p(npro,nshl))
              allocate (mxmudmi(nelblk)%p(npro,maxsh))
c
c.... save the element block
c
              n1=n
              n2=n+npro-1
              mater=1   ! all one material for now
              call gensav (ientp(n1:n2,1:nshl),
     &                     mater,           mien(nelblk)%p,
     &                     mmat(nelblk)%p)
              iel=iel+npro
c
           enddo
           deallocate(ientp)
        enddo
        lcblk(1,nelblk+1) = iel
c
c.... return
c
CAD        call timer ('Back    ')
c
        return
c
1000    format(a80,//,
     &  ' N o d a l   C o n n e c t i v i t y',//,
     &  '   Elem  ',/,
     &  '  Number  ',7x,27('Node',i2,:,2x))
1100    format(2x,i5,6x,27i8)
        end
