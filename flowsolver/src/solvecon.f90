      subroutine solvecon(y,       x,      iBC,  BC,  &
                           iper,    ilwork, shp,  shgl)

!---------------------------------------------------------------------
! This subroutine is to calculate the constraint for the redistancing 
! scalar of the level set method. This is to prevent interface from 
! moving by applying the condition that the volume must stay constant 
! in each element when the redistance step is applied.
!--------------------------------------------------------------------
!
!
      use pointer_data
!     
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
      !include "auxmpi.h"      
!      
      dimension y(nshg,ndof),                     &
                x(numnp,nsd),            iBC(nshg), &
                BC(nshg,ndofBC),         ilwork(nlwork), &
                iper(nshg)
!
!     
      dimension shp(MAXTOP,maxsh,MAXQPT),   &
                 shgl(MAXTOP,nsd,maxsh,MAXQPT),  &
                 v_lambda(nshg),   hprime(nshg), &
                 v_lambda1(nshg), v_lambda2(nshg), &
                 rmass(nshg) 
!
        real*8, allocatable :: tmpshp(:,:),  tmpshgl(:,:,:)
        real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

!
! ... intialize
!       
      rmass  = zero
      v_lambda = zero
      v_lambda1 = zero
      v_lambda2 = zero
      hprime = zero
!
! ... loop over element blocks
!
      do iblk = 1, nelblk
!
!.... set up the parameters
!
         nenl   = lcblk(5,iblk) ! no. of vertices per element
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
!.... compute and assemble the constarint factor, and mass matrix
!     

         allocate (tmpshp(nshl,MAXQPT))
         allocate (tmpshgl(nsd,nshl,MAXQPT))

         tmpshp(1:nshl,:) = shp(lcsyst,1:nshl,:)
         tmpshgl(:,1:nshl,:) = shgl(lcsyst,:,1:nshl,:)
!   
         call volcon (y,          x,             tmpshp,               &
                      tmpshgl,    mien(iblk)%p,  rmass,      &
                      v_lambda1,  hprime,        v_lambda2)

         deallocate ( tmpshp )
         deallocate ( tmpshgl ) 
      enddo

!
! ... multiple processor communication
!
      if (numpe > 1) then
         call commu (v_lambda1  , ilwork, 1  , 'in ')
         call commu (v_lambda2  , ilwork, 1  , 'in ')
         call commu (hprime     , ilwork, 1  , 'in ')
         call commu (rmass      , ilwork, 1  , 'in ')
      endif
!
!.... take care of periodic boundary conditions
!
      do j= 1,nshg
         if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            v_lambda1(i) =  v_lambda1(i) +  v_lambda1(j)
            v_lambda2(i) =  v_lambda2(i) +  v_lambda2(j)
            hprime(i)   =  hprime(i) + hprime(j)
         endif
      enddo
!
      do j= 1,nshg
         if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(j) = rmass(i)
            v_lambda1(j) =  v_lambda1(i)
            v_lambda2(j) =  v_lambda2(i)
            hprime(j) = hprime(i)
         endif
      enddo
!
! ... calculation of constraint factor
!
      rmass  = one/rmass
      v_lambda1 = v_lambda1*rmass ! numerator of lambda
      v_lambda2 = v_lambda2*rmass ! denominator of lambda
      v_lambda  = v_lambda1/(v_lambda2+epsM**2)
      hprime    = hprime*rmass
      v_lambda  = v_lambda*hprime  
!
! ... commu out for the multiple processor
!
      if(numpe > 1) then
         call commu (v_lambda, ilwork, 1, 'out')    
      endif
!          
! ... the following commented lines are for the different way of getting 
!     the denominator of constraint (lambda) calculation 
!$$$                hprime=zero
!$$$                do kk=1, nshg
!$$$                   if (abs (y(kk,6)) .le. epsilon_ls) then
!$$$                      hprime(kk) = (0.5/epsilon_ls) * (1 
!$$$     &                   + cos(pi*y(kk,6)/epsilon_ls))
!$$$                   endif
!$$$                enddo
!$$$                y(:,7) = y(:,7)+v_lambda*hprime/dtgl
!
! ... the vlome constraint applied on the second scalar
!
      y(:,7) = y(:,7)+v_lambda/dtgl   
!  
      return
      end
!
!
!

      subroutine volcon (y,         x,      shp,       &
                         shgl,      ien,    rmass,  &
                         v_lambda1, hprime, v_lambda2)

!---------------------------------------------------------------------
!
! This subroutine is to calculate the element contribution to the 
! constraint factor and mass matrix.
!
!---------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!     
      dimension y(nshg,ndof),               x(numnp,nsd),               &
                  shp(nshl,maxsh),   &
                  shgl(nsd,nshl,maxsh), &
                  ien(npro,nshl), &
                  qres(nshg,idflx),         rmass(nshg)
!
!.... element level declarations
!
      dimension ycl(npro,nshl,ndof),      xl(npro,nenl,nsd),          &
                rmassl(npro,nshl)     
      dimension sgn(npro,nshape),         v_lambdal1(npro,nshl), &
                v_lambda1(nshg),          hprimel(npro,nshl), &
                hprime(nshg),             v_lambdal2(npro,nshl), &
                v_lambda2(nshg)
!
! local arrays
!
      dimension shg(npro,nshl,nsd), &
                dxidx(npro,nsd,nsd),      WdetJ(npro)
!
      dimension shapeVar(npro,nshl), &
                shdrv(npro,nsd,nshl)
!
!
!.... for volume constraint calculation of redistancing step
!
      dimension Sclr(npro),              Sclrtmp(npro), &
                h_prime(npro),           tmp1(npro),  &
                tmp2(npro), &
                v_lambdatmp(npro),       v_lambdal(npro,nshl)
!$$$     &          ,hprimel(npro,nshl),      v_lambdal1(npro,nshl),
!$$$     &          v_lambdal2(npro,nshl)
! above arrays must be uncommented for alternate method included below (commented)
      real epsilon_tmp

!
!.... create the matrix of mode signs for the hierarchic basis 
!     functions. 
!
      if (ipord .gt. 1) then
         call getsgn(ien,sgn)
      endif
!
!.... gather the variables
!

      call localy(y,      ycl,     ien,    ndof,   'gather  ')
      call localx(x,      xl,      ien,    nsd,    'gather  ')
!
!.... get the element contributions of the numerator and denominator 
!     of lambda 
!
      rmassl = zero
      v_lambdal1= zero
      v_lambdal2= zero
      hprimel = zero
!
!.... loop through the integration points
!
      do intp = 1, ngauss
         if (Qwt(lcsyst,intp) .eq. zero) cycle ! precaution
!
!.... create a matrix of shape functions (and derivatives) for each
!     element at this quadrature point. These arrays will contain 
!     the correct signs for the hierarchic basis
!
         call getshp(shp,          shgl,      sgn,  &
                     shapeVar,        shdrv)
!
!.... initialize
!     
         h_prime   = zero        
	 sclr      = zero
	 sclrtmp   = zero
         v_lambdatmp=zero
         tmp1       =zero
         tmp2       =zero

!
!.... --------------------->  Element Metrics  <-----------------------
!
         call e3metric( xl,         shdrv,        dxidx,   &
                        shg,        WdetJ)

!
         do i = 1, nshl 
!
!  y(intp)=SUM_{a=1}^nshl (N_a(intp) Ya)
!     
            Sclr    = Sclr    + shapeVar(:,i) * ycl(:,i,7) !d^kbar
            sclrtmp = sclrtmp + shapeVar(:,i) * ycl(:,i,6) !d^0
         enddo

         if (isclr .eq. 2) then
            epsilon_tmp = epsilon_lsd
         else
            epsilon_tmp = epsilon_ls
         endif

         do i=1,npro
            if (abs (Sclrtmp(i)) .le. epsilon_tmp) then
               h_prime(i) = (0.5/epsilon_tmp) * (1  &
                          + cos(pi*Sclrtmp(i)/epsilon_tmp))
               tmp1(i)=-h_prime(i)*(sclr(i)-sclrtmp(i))*dtgl
               tmp2(i)=h_prime(i)**2
!              v_lambdatmp(i)=tmp1(i)/(tmp2(i)+ epsM)
            endif
         enddo
!
         do i=1,nshl
            v_lambdal1(:,i)=v_lambdal1(:,i)+ shapeVar(:,i)*WdetJ*tmp1
            v_lambdal2(:,i)=v_lambdal2(:,i)+ shapeVar(:,i)*WdetJ*tmp2
!           v_lambdal(:,i)=v_lambdal(:,i)+ shapeVar(:,i)*WdetJ*v_lambdatmp
            hprimel(:,i) = hprimel(:,i) + shapeVar(:,i)*WdetJ*h_prime
            rmassl(:,i)  =rmassl(:,i) + shapeVar(:,i)*WdetJ
         enddo
!
!.... end of the loop over integration points
!
      enddo
!
      call local (v_lambda1,  v_lambdal1, ien,  1,  'scatter ')
      call local (v_lambda2,  v_lambdal2, ien,  1,  'scatter ')
      call local (hprime,     hprimel,    ien,  1,  'scatter ')
      call local (rmass,      rmassl,     ien,  1,  'scatter ') 
!
      return
      end

