!
!  Copyright (c) 2000-2007, Stanford University, 
!     Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
!     Charles A. Taylor (see SimVascular Acknowledgements file 
!     for additional contributors to the source code).
!
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without 
!  modification, are permitted provided that the following conditions 
!  are met:
!
!  Redistributions of source code must retain the above copyright notice,
!  this list of conditions and the following disclaimer. 
!  Redistributions in binary form must reproduce the above copyright 
!  notice, this list of conditions and the following disclaimer in the 
!  documentation and/or other materials provided with the distribution. 
!  Neither the name of the Stanford University or Rensselaer Polytechnic
!  Institute nor the names of its contributors may be used to endorse or
!  promote products derived from this software without specific prior 
!  written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
!  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
!  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
!  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
!  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
!  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
!  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
!  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
!  DAMAGE.
!
!
        subroutine hessian ( y,         x,      &
                             shp,       shgl,      iBC, &
                             shpb,      shglb,     iper,       &
                             ilwork,    uhess,     gradu  )
        use pointer_data  ! brings in the pointers for the blocked arrays

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndof),       &
                  x(numnp,nsd),         iBC(nshg),            &
                  iper(nshg)
!
        dimension shp(MAXTOP,maxsh,MAXQPT),   &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                  shpb(MAXTOP,maxsh,MAXQPT), &
                  shglb(MAXTOP,nsd,maxsh,MAXQPT) 
!
        dimension gradu(nshg,9),     rmass(nshg), &
                  uhess(nshg,27)
!
        dimension ilwork(nlwork)

!
           gradu = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
!     
              call velocity_gradient ( y,                 &
                                       x,                        &
                                       shp(lcsyst,1:nshl,:),  &
                                       shgl(lcsyst,:,1:nshl,:), &
                                       mien(iblk)%p,      &
                                       gradu,  &
                                       rmass )

           end do

!
           call reconstruct( rmass, gradu, iBC, iper, ilwork, 9 )       

           uhess = zero
           rmass = zero
        
           do iblk = 1, nelblk
              iel    = lcblk(1,iblk)
              lelCat = lcblk(2,iblk)
              lcsyst = lcblk(3,iblk)
              iorder = lcblk(4,iblk)
              nenl   = lcblk(5,iblk) ! no. of vertices per element
              nshl   = lcblk(10,iblk)
              mattyp = lcblk(7,iblk)
              ndofl  = lcblk(8,iblk)
              nsymdl = lcblk(9,iblk)
              npro   = lcblk(1,iblk+1) - iel 
              ngauss = nint(lcsyst)
!     

              call velocity_hessian (  gradu,                 &
                                       x,                        &
                                       shp(lcsyst,1:nshl,:),  &
                                       shgl(lcsyst,:,1:nshl,:), &
                                       mien(iblk)%p,      &
                                       uhess, &
      				       rmass  )    
           end do
       

           call reconstruct( rmass, uhess, iBC, iper, ilwork, 27 )       
!
      return
      end

!-----------------------------------------------------------------------------

        subroutine velocity_gradient ( y,       x,       shp,     shgl,  &
                                       ien,     gradu,   rmass    )

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension y(nshg,ndof),               x(numnp,nsd),             &
                  shp(nshl,ngauss),           shgl(nsd,nshl,ngauss), &
                  ien(npro,nshl),             gradu(nshg,9),  &
                  shdrv(npro,nsd,nshl),       shapeVar( npro, nshl ),       &
                  gradul(npro,9) ,            rmass( nshg ) 
!
        dimension yl(npro,nshl,ndof),          xl(npro,nenl,nsd), &
                  ql(npro,nshl,9),             dxidx(npro,nsd,nsd), &
                  WdetJ(npro),		       rmassl(npro,nshl)
!
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

        call localy (y,    yl,     ien,    ndof,   'gather  ')
        call localx (x,    xl,     ien,    nsd,    'gather  ')
!
!.... get the element residuals 
!
        ql     = zero
        rmassl = zero

        do intp = 1, ngauss

            if ( Qwt( lcsyst, intp ) .eq. zero ) cycle

            gradul = zero
            call getshp( shp, shgl, sgn, shapeVar, shdrv )
            call local_gradient( yl(:,:,2:4), 3,  shdrv, xl,  &
                                 gradul , dxidx, WdetJ )

!.... assemble contribution of gradu to each element node
!     
            do i=1,nshl
                do j = 1, 9
                    ql(:,i,j) = ql(:,i,j)+shapeVar(:,i)*WdetJ*gradul(:,j)
                end do
                
                rmassl(:,i) = rmassl(:,i) + shapeVar(:,i)*WdetJ

             end do

        end do
!
!
        call local (gradu,  ql,     ien,  9,  'scatter ')
        call local (rmass,  rmassl, ien,  1,  'scatter ')
!
!.... end
!
        return
        end


!-----------------------------------------------------------------------------

        subroutine velocity_hessian ( gradu,   x,     shp,   shgl,  &
                                      ien,     uhess, rmass  )

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension gradu(nshg,9),              x(numnp,nsd),             &
                  shp(nshl,ngauss),           shgl(nsd,nshl,ngauss), &
                  ien(npro,nshl),             uhess(nshg,27),  &
                  shdrv(npro,nsd,nshl),       shapeVar( npro, nshl ), &
                  uhessl(npro,27),            rmass( nshg ) 
!
        dimension gradul(npro,nshl,9),          xl(npro,nenl,nsd),          &
                  ql(npro,nshl,27),             dxidx(npro,nsd,nsd),     &
                  WdetJ(npro),                  rmassl(npro, nshl)
!
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

        call local  (gradu,  gradul, ien,    9 ,   'gather  ')
        call localx (x,      xl,     ien,    nsd,  'gather  ')
!
!.... get the element residuals 
!
        ql     = zero
	rmassl = zero

        do intp = 1, ngauss

            if ( Qwt( lcsyst, intp ) .eq. zero ) cycle

            uhessl = zero
            call getshp( shp, shgl, sgn, shapeVar, shdrv )
            call local_gradient( gradul, 9,  shdrv, xl,  &
                                 uhessl , dxidx, WdetJ )

!.... assemble contribution of gradu .,
!     
            do i=1,nshl
                do j = 1,27 
                    ql(:,i,j)=ql(:,i,j)+shapeVar(:,i)*WdetJ*uhessl(:,j )
                end do

                rmassl(:,i) = rmassl(:,i) + shapeVar(:,i)*WdetJ
             end do

        end do
!
!
        call local (uhess,  ql,     ien,  27,     'scatter ')
        call local (rmass,  rmassl, ien,   1,     'scatter ')
!
!.... end
!
        return
        end


!--------------------------------------------------------------------
      subroutine reconstruct( rmass, qres, iBC, iper, ilwork, vsize )

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      integer vsize
      dimension rmass(nshg), qres( nshg, vsize), &
                iBC(nshg), iper(nshg)
!
!
!.... compute qi for node A, i.e., qres <-- qres/rmass
!
       if (numpe > 1) then
          call commu ( qres  , ilwork,  vsize  , 'in ')
          call commu ( rmass , ilwork,  1  , 'in ')
       endif
!
!  take care of periodic boundary conditions
!
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            qres(i,:) = qres(i,:) + qres(j,:)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(j) = rmass(i)
            qres(j,:) = qres(i,:)
          endif
        enddo
!
!.... invert the diagonal mass matrix and find q
!
        rmass = one/rmass
       
       do i=1,vsize 
          qres(:,i) = rmass*qres(:,i)
       enddo

       if(numpe > 1) then
          call commu (qres, ilwork, vsize, 'out')    
       endif

!.... return
!    
        return
        end

!-------------------------------------------------------------------------

        subroutine local_gradient ( vector,   vsize, shgl,   xl,  &
                                    gradient, dxidx,   WdetJ )
!
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
!  passed arrays

        integer vsize
!
        dimension vector(npro,nshl,vsize),  &
                  shgl(npro,nsd,nshl),        xl(npro,nenl,nsd), &
                  gradient(npro,vsize*3),     shg(npro,nshl,nsd),  &
                  dxidx(npro,nsd,nsd),        WdetJ(npro)
!
!  local arrays
!
        dimension tmp(npro),           dxdxi(npro,nsd,nsd)

!
!.... compute the deformation gradient
!
        dxdxi = zero
!
          do n = 1, nenl
            dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(:,1,n)
            dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(:,2,n)
            dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(:,3,n)
            dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(:,1,n)
            dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(:,2,n)
            dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(:,3,n)
            dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(:,1,n)
            dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(:,2,n)
            dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(:,3,n)
          enddo
!
!.... compute the inverse of deformation gradient
!
        dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)  &
                       - dxdxi(:,3,2) * dxdxi(:,2,3)
        dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)  &
                       - dxdxi(:,1,2) * dxdxi(:,3,3)
        dxidx(:,1,3) =   dxdxi(:,1,2) * dxdxi(:,2,3)  &
                       - dxdxi(:,1,3) * dxdxi(:,2,2)
        tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)  &
                             + dxidx(:,1,2) * dxdxi(:,2,1)   &
                             + dxidx(:,1,3) * dxdxi(:,3,1) )
        dxidx(:,1,1) = dxidx(:,1,1) * tmp
        dxidx(:,1,2) = dxidx(:,1,2) * tmp
        dxidx(:,1,3) = dxidx(:,1,3) * tmp
        dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)  &
                      - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
        dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)  &
                      - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
        dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)  &
                      - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
        dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)  &
                      - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
        dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)  &
                      - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
        dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)  &
                      - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
!
        WdetJ = Qwt(lcsyst,intp)/ tmp

!
!.... --------------------->  Global Gradients  <-----------------------
!
        gradient = zero
!
!
        do n = 1, nshl
!
!.... compute the global gradient of shape-function
!
!            ! N_{a,x_i}= N_{a,xi_i} xi_{i,x_j}
!
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) +  &
                       shgl(:,2,n) * dxidx(:,2,1) + &
                       shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) +  &
                       shgl(:,2,n) * dxidx(:,2,2) + &
                       shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) +  &
                       shgl(:,2,n) * dxidx(:,2,3) + &
                       shgl(:,3,n) * dxidx(:,3,3) 
!
!
!  Y_{,x_i}=SUM_{a=1}^nenl (N_{a,x_i}(int) Ya)
!
          do i = 1, 3
            do j = 1, vsize
               k = (i-1)*vsize+j
               gradient(:,k) = gradient(:,k) + shg(:,n,i)*vector(:,n,j)
            end do
          end do

       end do

!
!.... return
!
       return
       end
