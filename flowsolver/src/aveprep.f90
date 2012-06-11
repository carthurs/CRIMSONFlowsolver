      module aveall

      real*8, allocatable :: ylist(:)
      integer, allocatable :: ifathe(:,:)

      end module
      
!   destructor
      subroutine Daveall
      use aveall
      deallocate(ylist)
      deallocate(ifathe)
      return
      end

!------------------------------------------------------------------------------

      subroutine setave

      use aveall

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      idim = numel*intg(1,1)
      allocate ( ylist(idim) )
      allocate ( ifathe(numel,maxsh) )

!.... return

      return
      end

!------------------------------------------------------------------------------

      subroutine aveprep(shp,x)

      use pointer_data
      use aveall

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      dimension shp(MAXTOP,maxsh,MAXQPT)
      dimension x(numnp,3)
      !integer   nlist

      nlist  = 0 
      ylist  = zero
      lfathe = 0
      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)

        call getylist( ylist, ifathe(iel:inum,:), shp(lcsyst,1:nshl,:), &
             x, mien(iblk)%p)
        

      enddo

      return
      end















      subroutine getnu (ien, strl, xmudmi, cdelsq, lfathe)

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension  ien(npro,nshl),       strl(npro,ngauss), &
                 xmudmi(npro,ngauss),   cdelsq(nlist), &
                 lfathe(npro,ngauss)         

      do intp = 1, ngauss
      do iel  = 1, npro

         xmudmi(iel,intp) = cdelsq( lfathe(iel,intp) )* &
              strl(iel,intp)

      enddo
      enddo

!
!  local clipping
!
      rmu=datmat(1,2,1)
!c      xmudmi=min(xmudmi,1000.0*rmu) !don't let it get larger than 1000 mu
!c      xmudmi=max(xmudmi, -rmu) ! don't let (xmudmi + mu) < 0

      do int = 1, ngauss
      xmudmi(:,int) = &
           min(xmudmi(:,int),1000.0*rmu) !don't let it get larger than 1000 mu
      xmudmi(:,int) = &
           max(xmudmi(:,int), -rmu)    ! don't let (xmudmi + mu) < 0
      enddo

      return
      end
      subroutine getxnut (fres, xden, xnum, ien, shp)

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension ien(npro,nshl),          fres(nshg,22), &
                shp(nshl,maxsh),         fresli(npro,22), &
                fresl(npro,nshl,22),     strnrmi(npro), &
                xnutl(npro),             fwr(npro), &
                xmij(npro,6),            xlij(npro,6), &
                xdenli(npro),            xnumli(npro), &
                xdenl(npro),             xnuml(npro), &
                xden(1),                 xnum(1)         
     

      call local (fres, fresl, ien, 22, 'gather  ')

      xnuml  = zero
      xdenl  = zero

      do intp = 1, ngauss

         fresli = zero
         do i = 1, nenl
            do j = 1, 22
               fresli(:,j) = fresli(:,j) + shp(i,intp)*fresl(:,i,j)
            enddo
         enddo

         strnrmi(:) =  sqrt( &
         two * (fresli(:,10)**2 + fresli(:,11)**2 + fresli(:,12)**2) &
         + four * ( fresli(:,13)**2 + fresli(:,14)**2 + &
         fresli(:,15)**2 ) )

         fwr = fwr1 * fresli(:,22) * strnrmi(:)

         xmij(:,1) = -fwr  &
              * fresli(:,10) + fresli(:,16) 

         xmij(:,2) = -fwr  &
              * fresli(:,11) + fresli(:,17)

         xmij(:,3) = -fwr  &
              * fresli(:,12) + fresli(:,18)

         xmij(:,4) = -fwr * fresli(:,13) + fresli(:,19)
         xmij(:,5) = -fwr * fresli(:,14) + fresli(:,20)
         xmij(:,6) = -fwr * fresli(:,15) + fresli(:,21)
         
         fresli(:,22) = one/fresli(:,22) 

         xlij(:,1) = fresli(:,4) - fresli(:,1) * &
              fresli(:,1) * fresli(:,22)
         xlij(:,2) = fresli(:,5) - fresli(:,2) * &
              fresli(:,2) * fresli(:,22)
         xlij(:,3) = fresli(:,6) - fresli(:,3) * &
              fresli(:,3) * fresli(:,22)
         xlij(:,4) = fresli(:,7) - fresli(:,1) * &
              fresli(:,2) * fresli(:,22)
         xlij(:,5) = fresli(:,8) - fresli(:,1) * &
              fresli(:,3) * fresli(:,22)
         xlij(:,6) = fresli(:,9) - fresli(:,2) * &
              fresli(:,3) * fresli(:,22)

         xnumli(:) =    xlij(:,1) * xmij(:,1) &
              + xlij(:,2) * xmij(:,2) &
              + xlij(:,3) * xmij(:,3) &
              + two * (xlij(:,4) * xmij(:,4) &
              + xlij(:,5) * xmij(:,5) &
              + xlij(:,6) * xmij(:,6))
         xdenli(:) =        xmij(:,1) * xmij(:,1) &
              + xmij(:,2) * xmij(:,2) &
              + xmij(:,3) * xmij(:,3) &
              + two * (xmij(:,4) * xmij(:,4) &
              + xmij(:,5) * xmij(:,5) &
              + xmij(:,6) * xmij(:,6))
         
!         xdenli = xdenli

         xdenl(:) = xdenl(:) + xdenli(:)
         xnuml(:) = xnuml(:) + xnumli(:)

      enddo

      do nel = 1, npro
         xden(1) = xden(1) + xdenl(nel)
         xnum(1) = xnum(1) + xnuml(nel)
      enddo 

      return

      end









      subroutine getylist (ylist, lfathe, shp, x, ien)

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension lfathe(npro,maxsh),    shp(nshl,maxsh),  &
                x(numnp,3),            ien(npro,nshl), &
                xl(npro,nenl,nsd),     ylist(idim)
                
      
      call localx (x,      xl,     ien,    3,  'gather  ')      
                

      do nel  = 1, npro
      do intp = 1, ngauss

!... Compute the y-coordinate of the current quad. pt. and call it yint.

         yint = zero
         xint = zero
         zint = zero
         do n = 1, nenl
            yint = yint + xl(nel,n,2)*shp(n,intp)
            xint = xint + xl(nel,n,1)*shp(n,intp)
            zint = zint + xl(nel,n,3)*shp(n,intp)
         enddo

!... Check ylist to see if yint is already included in ylist.

         if(nlist .eq. 0)then ! In this case yint is definitely not in ylist.
            nlist            = nlist + 1
            ylist(nlist)     = yint
            lfathe(nel,intp) = nlist

         else           

            do ilist = 1, nlist
               imatch = ilist
               if( abs(yint-ylist(ilist)) .lt. 0.00001 ) then
                  lfathe(nel,intp) = imatch
                  goto 5
               endif
            enddo

            nlist            = nlist + 1
            ylist(nlist)     = yint
            lfathe(nel,intp) = nlist
                  
         endif

 5       yint = zero

!         if( nlist .eq. 5)then
!            write(*,*)'ylist yint=',ylist(4),ylist(5)
!            stop
!         endif

      enddo ! End loop over quad. pts.
      enddo ! End loop over elements in current block.
         
      return
      end
      subroutine projdmc (y,      shgl,      shp, &
                         iper,   ilwork,       x)

      use pointer_data

      use aveall   ! This module helps us average cdelsq computed at quad pts

      use quadfilt   ! This module gives us shglf(maxtp,nsd,maxsh,ngaussf),
!                    shpf(maxtp,maxsh,ngaussf), and Qwtf(maxtp,ngaussf). 
!                    Shpf and shglf are the shape funciotns and their 
!                    gradient evaluated using the quadrature rule desired 
!                    for computing the dmod. Qwt contains the weights of the 
!                    quad. points.  

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"

!
      dimension fres(nshg,24),         fwr(nshg), &
                strnrm(nshg),          cdelsq(nshg), &
                strl(numel,ngauss), &
                y(nshg,5),            yold(nshg,5), &
                iper(nshg), &
                ilwork(nlwork), &
                x(numnp,3), &
                shgl(MAXTOP,nsd,maxsh,MAXQPT), shp(MAXTOP,maxsh,MAXQPT), &
                xden(idim),                    xnum(idim), &
                xdent(idim),                   xnumt(idim)
     
      fres = zero
      yold(:,1)=y(:,4)
      yold(:,2:4)=y(:,1:3)
!
!  hack in an interesting velocity field (uncomment to test dmod)
!
!      yold(:,5) = 1.0  
!      yold(:,2) = 2.0*x(:,1) - 3*x(:,2) 
!      yold(:,3) = 3.0*x(:,1) + 4.0*x(:,2)
!      yold(:,4) = 4.0*x(:,1) + x(:,2) + x(:,3)
!      yold(:,1) = Rgas * yold(:,5) ! Necessary to make model suitable
!                               suitable for the


      intrul=intg(1,itseq)
      intind=intpt(intrul)

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1

        ngauss  = nint(lcsyst)
        ngaussf = nintf(lcsyst)
        
        call asithf (yold, x, strl(iel:inum,:), mien(iblk)%p, fres, &
                     shglf(lcsyst,:,1:nshl,:), &
                     shpf(lcsyst,1:nshl,:),Qwtf(lcsyst,1:ngaussf))

      enddo
!
 
      if(numpe>1) call commu (fres, ilwork, 24, 'in ')
! 
! account for periodicity in filtered variables
!
      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(i,:) = fres(i,:) + fres(j,:)
        endif
      enddo

      do j = 1,nshg
        i = iper(j)
        if (i .ne. j) then
           fres(j,:) = fres(i,:)
        endif
      enddo

      if(numpe>1)then
         call commu (fres, ilwork, 24, 'out')
      endif

      fres(:,23) = one / fres(:,23)
      do j = 1,22
        fres(:,j) = fres(:,j) * fres(:,23)
      enddo

      xden   = zero
      xdent  = zero
      xnum   = zero
      xnumt  = zero

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        inum  = iel + npro - 1
        ngauss = nint(lcsyst)
        call smagcoeff (fres(:,1:22), xden, xnum, mien(iblk)%p, &
                        shp(lcsyst,1:nshl,:), ifathe(iel:inum,:))

      enddo      	

      if(numpe.gt.1)then
         call drvAllreduce( xden, xdent, idim )
         call drvAllreduce( xnum, xnumt, idim )
      else
         xdent  = xden
         xnumt  = xnum
      endif

      do i = 1, nlist
         cdelsq(i) = pt5*xnumt(i)/xdent(i)
      enddo

!...  Uncomment for averaging over all directions

!      sumc1 = 0.0
!      sumc2 = 0.0
!      do i = 1, nlist
!         sumc1 = sumc1 + xnumt(i)
!         sumc2 = sumc2 + xdent(i)
!      enddo
!      xfact = pt5*sumc1/sumc2
!      cdelsq(:) = xfact
!      if(myrank .eq. master)then
!         write(22,*)cdelsq(1)
!      endif

!...  End of averaging over all directions.

      do iblk = 1,nelblk
        lcsyst = lcblk(3,iblk) 
        iel  = lcblk(1,iblk) !Element number where this block begins
        npro = lcblk(1,iblk+1) - iel
        lelCat = lcblk(2,iblk)
        nenl = lcblk(5,iblk)
        nshl = lcblk(10,iblk)
        inum  = iel + npro - 1
        
        ngauss = nint(lcsyst)
        ngaussf = nintf(lcsyst)
 
        if (ngauss .ne. ngaussf) then
        call getstrl (yold, x,      mien(iblk)%p, &
                     strl(iel:inum,:), shgl(lcsyst,:,1:nshl,:), &
                     shp(lcsyst,1:nshl,:))
        endif

      enddo

      do iblk = 1,nelblk
         lcsyst = lcblk(3,iblk)
         iel  = lcblk(1,iblk)
         npro = lcblk(1,iblk+1) - iel
         lelCat = lcblk(2,iblk)
         inum  = iel + npro - 1
         
         ngauss = nint(lcsyst)

         call getnu (mien(iblk)%p, strl(iel:inum,:), &
              mxmudmi(iblk)%p,cdelsq,ifathe(iel:inum,:))
      enddo


      return
      end



      subroutine smagcoeff (fres, xden, xnum, ien, shp, lfathe)

      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

      dimension ien(npro,nshl),          fres(nshg,22), &
                shp(nshl,maxsh),         fresli(npro,22), &
                fresl(npro,nshl,22),     strnrmi(npro), &
                xnutl(npro),             fwr(npro), &
                xmij(npro,6),            xlij(npro,6), &
                xdenli(npro,ngauss),      xnumli(npro,ngauss), &
                xden(idim),              xnum(idim), &
                lfathe(npro,maxsh)

      call local (fres, fresl, ien, 22, 'gather  ')

      do intp = 1, ngauss

         fresli = zero
         do i = 1, nenl
            do j = 1, 22
               fresli(:,j) = fresli(:,j) + shp(i,intp)*fresl(:,i,j)
            enddo
         enddo

         strnrmi(:) =  sqrt( &
         two * (fresli(:,10)**2 + fresli(:,11)**2 + fresli(:,12)**2) &
         + four * ( fresli(:,13)**2 + fresli(:,14)**2 + &
         fresli(:,15)**2 ) )

         fwr = fwr1 * fresli(:,22) * strnrmi(:)

         xmij(:,1) = -fwr &
              * fresli(:,10) + fresli(:,16) 

         xmij(:,2) = -fwr &
              * fresli(:,11) + fresli(:,17)

         xmij(:,3) = -fwr &
              * fresli(:,12) + fresli(:,18)

         xmij(:,4) = -fwr * fresli(:,13) + fresli(:,19)
         xmij(:,5) = -fwr * fresli(:,14) + fresli(:,20)
         xmij(:,6) = -fwr * fresli(:,15) + fresli(:,21)
         
         fresli(:,22) = one/fresli(:,22) 

         xlij(:,1) = fresli(:,4) - fresli(:,1) * &
              fresli(:,1) * fresli(:,22)
         xlij(:,2) = fresli(:,5) - fresli(:,2) * &
              fresli(:,2) * fresli(:,22)
         xlij(:,3) = fresli(:,6) - fresli(:,3) * &
              fresli(:,3) * fresli(:,22)
         xlij(:,4) = fresli(:,7) - fresli(:,1) * &
              fresli(:,2) * fresli(:,22)
         xlij(:,5) = fresli(:,8) - fresli(:,1) * &
              fresli(:,3) * fresli(:,22)
         xlij(:,6) = fresli(:,9) - fresli(:,2) * &
              fresli(:,3) * fresli(:,22)

         xnumli(:,intp) =    xlij(:,1) * xmij(:,1) &
              + xlij(:,2) * xmij(:,2) &
              + xlij(:,3) * xmij(:,3) &
              + two * (xlij(:,4) * xmij(:,4) &
              + xlij(:,5) * xmij(:,5) &
              + xlij(:,6) * xmij(:,6))
         xdenli(:,intp) =        xmij(:,1) * xmij(:,1) &
              + xmij(:,2) * xmij(:,2) &
              + xmij(:,3) * xmij(:,3) &
              + two * (xmij(:,4) * xmij(:,4) &
              + xmij(:,5) * xmij(:,5) &
              + xmij(:,6) * xmij(:,6))
         
      enddo ! End loop over quad pts

!... For debugging

!      do nel = 1, npro
!      do intp = 1, ngauss

!         if ( mod(lfathe(nel,intp),2) .eq. 0 )then
!            xnumli(nel,intp) = 2.0
!            xdenli(nel,intp) = 3.0
!         else
!            xnumli(nel,intp) = 5.0
!            xdenli(nel,intp) = 2.5
!         endif

!      enddo
!      enddo

!...  End of debugging stuff.


      do intp = 1, ngauss
      do nel  = 1, npro

         xden( lfathe(nel,intp) ) = xden( lfathe(nel,intp) ) + &
              xdenli(nel,intp)

         xnum( lfathe(nel,intp) ) = xnum( lfathe(nel,intp) ) + &
              xnumli(nel,intp)

      enddo 
      enddo

      return

      end









