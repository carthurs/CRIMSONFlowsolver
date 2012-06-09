      module deformableWall
     
      parameter ( MAXBLK2 = 4000 ) ! Note compiler was complaining 
c                                    because MAXBLK in common.h be careful
c                                    to change both places


      type r1d
         real*8, pointer :: p(:)
      end type

      type r2d
         real*8, pointer :: p(:,:)
      end type
         
      type r4d
         real*8, pointer :: p(:,:,:,:)
      end type

      type i1d
         integer, pointer :: p(:)
      end type
      
      type i2d
         integer, pointer :: p(:,:)
      end type
      
      

      ! currently dynamically allocated in genbkb
      type (r2d), dimension(MAXBLK2) ::  mSWB, mTWB, mEWB         
      
      ! currently allocated in gendat
      type (i1d) :: mWNodes, mWNodes_gtlmap 
      
      ! these are dynamically allocated in a local routine
      type (r2d), dimension(MAXBLK2) :: mPS_global
      type (r4d), dimension(MAXBLK2) :: mKwall_xKebe
      
      ! stored rhs for the DDF state filter
      type (r1d), dimension(MAXBLK2) :: mPS_sfRHS
      
      ! row numbers and column pointers
      ! allocated in itrdrv
      type (i2d) :: rowp_wall
      type (i1d) :: colm_wall
      
      ! nonzero count for the deformable wall
      integer nnz_tot_wall
      
      
      contains
      
      
      
      subroutine vlmwStTri(x,     iBC,      BC)
      
      use pointer_data
      
      include "common.h"
      
      real*8 x(numnp,nsd)
      
      real*8  BC(nshg,ndofBC)
      integer iBC(nshg)
            
c
c.... loop over the boundary elements
c
      do iblk = 1, nelblb
c
c.... set up the parameters
c
         iel    = lcblkb(1,iblk)
         lelCat = lcblkb(2,iblk)
         lcsyst = lcblkb(3,iblk)
         iorder = lcblkb(4,iblk)
         nenl   = lcblkb(5,iblk)  ! no. of vertices per element
         nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
         nshl   = lcblkb(9,iblk)
         nshlb  = lcblkb(10,iblk)
         mattyp = lcblkb(7,iblk)
         ndofl  = lcblkb(8,iblk)
         npro   = lcblkb(1,iblk+1) - iel
         
c
c.... allocate space for the arrays
c     
    
         if (.not.associated(mPS_global(iblk)%p)) then
            allocate(mPS_global(iblk)%p(npro,9))
         end if
         
         if (.not.associated(mKwall_xKebe(iblk)%p)) then
            allocate(mKwall_xKebe(iblk)%p(npro,9,nshl,nshl))
         end if

c
c.... check if SWB, TWB, or EWB files are used
c.... if no, then use the entered uniform values
c              
         
         if (iUseSWB.eq.0) then
            mSWB(iblk)%p(:,1) = thicknessvw
         end if
         
         if (iUseTWB.eq.0) then
            mTWB(iblk)%p(:,1) = tissSuppStiffCoeff
            mTWB(iblk)%p(:,2) = tissSuppDampCoeff
         end if
         
         if (iUseEWB.eq.0) then
            mEWB(iblk)%p(:,1) = stateFilterCoeff
         end if
         
c
c.... pre-assemble the element matrices for the deformable wall
c     
             
         call asbvlmwsttri(x,   
     &                     mien(iblk)%p,         iBC,           BC,
     &                     mienb(iblk)%p,        miBCB(iblk)%p,
     &                     mSWB(iblk)%p,
     &                     mPS_global(iblk)%p,
     &                     mKwall_xKebe(iblk)%p)
         
      end do
            
      end subroutine
      
      
      
      subroutine asbvlmwsttri(x,           
     &                        ien,         iBC,   BC,
     &                        ienb,        iBCB,  
     &                        SWB,
     &                        PS_global,
     &                        Kwall_xKebe)
      
      include "common.h"
      
      real*8    x(numnp,nsd)
      
      integer   iBCB(npro,ndiBCB),         ienb(npro,nshl)
      integer   ien(npro,nshape)
      
      real*8    SWB(npro,nProps)          
      
      real*8    xlb(npro,nenl,nsd)
      
      real*8    rotel(npro,9,9),
     &          temp(npro),
     &          temp1(npro),               temp2(npro),
     &          temp3(npro),
     &          v1(npro,nsd),              v2(npro,nsd),
     &          v3(npro,nsd),              
     &          x1rot(npro,nsd),           
     &          x2rot(npro,nsd),           x3rot(npro,nsd),           
     &          detJacrot(npro),
     &          B(npro,5,9),               Dmatrix(npro,5,5),
     &          DtimesB(npro,5,9),         
     &          Kwall_local(npro,9,9),     Kwall_global(npro,9,9),
     &          Kwall_temp(npro,9,9),
     &          PS_local(npro,9),          PS_global(npro,9),
     &          Kwall_xKebe(npro,9,nshl,nshl)
     
      real*8    BC(nshg,ndofBC)
      integer   iBC(nshg)
     
      integer   i,j,k,l,m,n
              
c
c.... gather the variables
c

      call localx(x,xlb,ienb,nsd,'gather  ')

c
c.... check for element coordinate system
c      
      
      if (lcsyst.eq.1) then     
         ipt2=2
         ipt3=3
      else
         write(*,*) "triangular elements required"
         stop
      end if      
      
c     
c.... Assemble the stiffness LHS matrix and residual contribution matrix
c
                  
c     
c.... B^t * D * B formulation for plane stress enhanced membrane
c

c
c.... rotation matrix
c     
      v1 = xlb(:,ipt2,:) - xlb(:,1,:)
      temp = one / sqrt ( v1(:,1)**2 + v1(:,2)**2 + v1(:,3)**2 )
      v1(:,1) = v1(:,1) * temp
      v1(:,2) = v1(:,2) * temp
      v1(:,3) = v1(:,3) * temp
         
      v2 = xlb(:,ipt3,:) - xlb(:,1,:)
         
c     compute cross product
      temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
      temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
      temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
         
      temp     = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      v3(:,1) = temp1 * temp
      v3(:,2) = temp2 * temp
      v3(:,3) = temp3 * temp
         
c     cross product again for v2
      temp1 = v3(:,2) * v1(:,3) - v1(:,2) * v3(:,3)
      temp2 = v1(:,1) * v3(:,3) - v3(:,1) * v1(:,3)
      temp3 = v3(:,1) * v1(:,2) - v1(:,1) * v3(:,2)
         
      temp     = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
      v2(:,1) = temp1 * temp
      v2(:,2) = temp2 * temp
      v2(:,3) = temp3 * temp
      
      rotel = zero   
         
      rotel(:,1,1:3) = v1(:,:)
      rotel(:,2,1:3) = v2(:,:)
      rotel(:,3,1:3) = v3(:,:)
      
      rotel(:,4,4:6) = v1(:,:)
      rotel(:,5,4:6) = v2(:,:)
      rotel(:,6,4:6) = v3(:,:)
      
      rotel(:,7,7:9) = v1(:,:)
      rotel(:,8,7:9) = v2(:,:)
      rotel(:,9,7:9) = v3(:,:)
         
c     
c.... rotated coordinates
c     
      x1rot = zero
      x2rot = zero
      x3rot = zero
         
      do i = 1, 3
         do j = 1, 3 
            x1rot(:,i) = x1rot(:,i)+rotel(:,i,j)*xlb(:,1,j)
            x2rot(:,i) = x2rot(:,i)+rotel(:,i,j)*xlb(:,ipt2,j)
            x3rot(:,i) = x3rot(:,i)+rotel(:,i,j)*xlb(:,ipt3,j)
         end do
      end do
      
c     
c.... B matrices
c     
      B = zero
      
      detJacrot = (x2rot(:,1)-x1rot(:,1)) * (x3rot(:,2)-x1rot(:,2)) - 
     &  (x3rot(:,1)-x1rot(:,1)) * (x2rot(:,2)-x1rot(:,2))
      
      B(:,1,1) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
      B(:,2,2) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
      B(:,3,1) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
      B(:,3,2) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
      B(:,4,3) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
      B(:,5,3) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
         
      B(:,1,4) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
      B(:,2,5) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
      B(:,3,4) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
      B(:,3,5) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
      B(:,4,6) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
      B(:,5,6) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
         
      B(:,1,7) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
      B(:,2,8) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
      B(:,3,7) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
      B(:,3,8) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
      B(:,4,9) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
      B(:,5,9) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
      
c     
c.... D matrix
c     
      Dmatrix = zero
      
      if (iUseSWB.gt.0) then
         if (nProps.eq.10) then 
c        This is an Isotropic Material         

            Dmatrix(:,1,1) = SWB(:,7)
            Dmatrix(:,2,2) = SWB(:,7)
            Dmatrix(:,1,2) = SWB(:,8)
            Dmatrix(:,2,1) = SWB(:,8)
            Dmatrix(:,3,3) = SWB(:,9)
            Dmatrix(:,4,4) = SWB(:,10)
            Dmatrix(:,5,5) = SWB(:,10)
         
         elseif(nProps.eq.21) then

c        This is an Orthotropic Material
            Dmatrix(:,1,1) = SWB(:,7)
            
            Dmatrix(:,2,1) = SWB(:,8)
            Dmatrix(:,1,2) = SWB(:,8)
            Dmatrix(:,2,2) = SWB(:,9)
            
            Dmatrix(:,3,1) = SWB(:,10)
            Dmatrix(:,1,3) = SWB(:,10)
            Dmatrix(:,3,2) = SWB(:,11)
            Dmatrix(:,2,3) = SWB(:,11)
            Dmatrix(:,3,3) = SWB(:,12)
            
            Dmatrix(:,4,1) = SWB(:,13)
            Dmatrix(:,1,4) = SWB(:,13)
            Dmatrix(:,4,2) = SWB(:,14)
            Dmatrix(:,2,4) = SWB(:,14)
            Dmatrix(:,4,3) = SWB(:,15)
            Dmatrix(:,3,4) = SWB(:,15)
            Dmatrix(:,4,4) = SWB(:,16)
            
            Dmatrix(:,5,1) = SWB(:,17)
            Dmatrix(:,1,5) = SWB(:,17)
            Dmatrix(:,5,2) = SWB(:,18)
            Dmatrix(:,2,5) = SWB(:,18)
            Dmatrix(:,5,3) = SWB(:,19)
            Dmatrix(:,3,5) = SWB(:,19)
            Dmatrix(:,5,4) = SWB(:,20)
            Dmatrix(:,4,5) = SWB(:,20)
            Dmatrix(:,5,5) = SWB(:,21)
         
         else
         
            write(*,*) 'Number of wall properties not set correctly!'
            stop
         
         end if
      else
      
         Dmatrix(:,1,1) = evw/(one-rnuvw**2)
         Dmatrix(:,2,2) = Dmatrix(:,1,1)
         Dmatrix(:,1,2) = evw*rnuvw/(one-rnuvw**2)
         Dmatrix(:,2,1) = Dmatrix(:,1,2)
         Dmatrix(:,3,3) = evw*pt5*(1-rnuvw)/(one-rnuvw**2)
         Dmatrix(:,4,4) = Dmatrix(:,3,3)*rshearconstantvw
         Dmatrix(:,5,5) = Dmatrix(:,4,4)
      
      end if
      
c     
c.... D * [B1|B2|B3]
c     
      DtimesB = zero
      
      do i = 1, 5
         do j = 1, 9
            do k = 1, 5
               DtimesB(:,i,j)=DtimesB(:,i,j)+Dmatrix(:,i,k)*B(:,k,j)
            end do
         end do
      end do
         
c     
c.... [B1|B2|B3]^T * D * [B1|B2|B3]
c     
      Kwall_local = zero
         
      do i = 1, 9            
         do j = 1, 9         
            do k = 1, 5
               Kwall_local(:,i,j)=
     &            Kwall_local(:,i,j)+B(:,k,i)*DtimesB(:,k,j)
            end do
         end do
      end do
      
      
c     
c.... [B1|B2|B3]^T * Sigma  Compute prestress contribution
c            

      if (iUseSWB.gt.0) then
         do i = 1, 3            ! BEWARE index 3 and 4
            PS_local(:,i) = B(:,1,i) * SWB(:,2) +
     &                      B(:,2,i) * SWB(:,4) +
     &                      B(:,3,i) * SWB(:,3) +
     &                      B(:,4,i) * SWB(:,5) +
     &                      B(:,5,i) * SWB(:,6)

            PS_local(:,i+3) = B(:,1,i+3) * SWB(:,2) +
     &                        B(:,2,i+3) * SWB(:,4) +
     &                        B(:,3,i+3) * SWB(:,3) +
     &                        B(:,4,i+3) * SWB(:,5) +
     &                        B(:,5,i+3) * SWB(:,6)

            PS_local(:,i+6) = B(:,1,i+6) * SWB(:,2) +
     &                        B(:,2,i+6) * SWB(:,4) +
     &                        B(:,3,i+6) * SWB(:,3) +
     &                        B(:,4,i+6) * SWB(:,5) +
     &                        B(:,5,i+6) * SWB(:,6)
         end do
      end if
      
c
c.... rotate element stiffness to global coordinates
c     

      Kwall_global = zero
      Kwall_temp = zero
      
      do i = 1, 9
         do j = 1, 9
            do k= 1, 9
      
               Kwall_temp(:,i,j) = Kwall_temp(:,i,j) + 
     &            Kwall_local(:,i,k) * rotel(:,k,j)     
               
            end do
         end do
      end do

      do i = 1, 9
         do j = 1, 9
            do k= 1, 9
      
               Kwall_global(:,i,j) = Kwall_global(:,i,j) + 
     &            rotel(:,k,i) * Kwall_temp(:,k,j)     
               
            end do
         end do
      end do
      
c
c.... rotate local prestress contribution to global coordinates
c      
      PS_global = zero

      if (iUseSWB.gt.0) then
         do i =1, 9
            do j = 1, 9
               PS_global(:,i) = rotel(:,j,i) * PS_local(:,j) +
     &                          PS_global(:,i)
            
            end do
         end do
      end if
      
c      
c.... integrate over element - multiply the nodal matrices by the area and the thickness      
c

      do i =1, 9
         do j = 1, 9
         
            Kwall_global(:,i,j) = Kwall_global(:,i,j)*detJacrot(:)
     &                            *pt5*SWB(:,1)
     
         enddo
         
         PS_global(:,i) = PS_global(:,i)*detJacrot(:)*pt5*SWB(:,1)
         
      enddo
      
c
c.... collect the entries of the stiffness matrix into the xKebe format
c      

      Kwall_xKebe = zero
      
      do i=1,3
         do j=1,3
            do k=1,3
               do m=1,3
                  Kwall_xKebe(:,3*(k-1)+m,i,j) = 
     &               Kwall_global(:,3*(i-1)+k,3*(j-1)+m)
               end do
            end do
         end do
      end do   
      
c
c.... turbulence wall (as a way of checking for deformable wall stiffness)
c
  
      do iel = 1, npro
         if (btest(iBCB(iel,1),4)) then      
            
         else
            
            Kwall_xKebe(iel,:,:,:) = zero
         end if
      end do
      
c      
c.... assemble lhsK with the global wall equation numbers
c.... this step is for the state filter
c      
      
c      call fillsparseW(ienb, Kwall_xKebe, Kwall_lhsK,
c     &                 row,   col)

c      
c.... apply dirichlet conditions for components
c.... this step is for the state filter
c      

      call bc3lhs(iBC, BC, ienb, Kwall_xKebe)
      
      end subroutine
      
      
c      
c.... helper routine to set up and solve
c.... a linear system corresponding to the deformable wall
c.... for the state filter
c     
      
      subroutine solveWallProb(rowp,colm,ilwork,iBC,BC,iper)
      
      use pointer_data
      
      include "common.h"

      integer rowp(nnz_tot),            colm(nshg+1),
     &        iBC(nshg),                ilwork(nlwork),
     &        iper(nshg) 
     
	real*8  lhsK(9,nnz_tot),          BC(nshg,ndofBC)
	real*8  flowDiag(nshg,3)
	
	real*8  tempK(nshg*3,nshg*3),     tempRHS(nshg,3)
	
	real*8  tempRHSPre(nshg,3),       trx(nshg,3)

      integer i,j,k
      
      integer sparseloc
      
      real*8  rand
      
      lhsK = zero;
           
c
c.... loop over the boundary elements
c
      do iblk = 1, nelblb
c
c.... set up the parameters
c
         iel    = lcblkb(1,iblk)
         nshl   = lcblkb(9,iblk)
         nshlb  = lcblkb(10,iblk)
         npro   = lcblkb(1,iblk+1) - iel
         
c
c.... assemble into lhs form
c         
         call fillsparseK(mienb(iblk)%p,mKwall_xKebe(iblk)%p,lhsK,
     &                    rowp,colm)
     
c     
c.... assemble the residuals     
c         
      end do
 
c     
c.... construct the diagonal preconditioner
c 
      
      call drvlesPrepDiagK ( flowDiag, ilwork,
     &                       iBC,      BC,      iper,
     &                       rowp,     colm,    
     &                       lhsK )
     
     
     
      tempK = zero
      do i = 1, nshg
        do k = colm(i), colm(i+1)-1
        
          j = rowp(k)
          
          tempK(3*(i-1)+1,3*(j-1)+1) = lhsK(1,k)
          tempK(3*(i-1)+1,3*(j-1)+2) = lhsK(2,k)
          tempK(3*(i-1)+1,3*(j-1)+3) = lhsK(3,k)
          
          tempK(3*(i-1)+2,3*(j-1)+1) = lhsK(4,k)
          tempK(3*(i-1)+2,3*(j-1)+2) = lhsK(5,k)
          tempK(3*(i-1)+2,3*(j-1)+3) = lhsK(6,k)
          
          tempK(3*(i-1)+3,3*(j-1)+1) = lhsK(7,k)
          tempK(3*(i-1)+3,3*(j-1)+2) = lhsK(8,k)
          tempK(3*(i-1)+3,3*(j-1)+3) = lhsK(9,k)
            
        end do
      end do
      
!      open (unit = 2, file = "lhsK.dat")
!      
!      do i=1,nshg*3
!        do j=1,nshg*3 
!          write(2,*) tempK(i,j)
!        enddo
!      enddo
!
!      close (2)
     
c
c.... apply the diagonal preconditioner to the lhs matrix
c     
      
      do i = 1, nshg
        do k = colm(i), colm(i+1)-1
          
          j = rowp(k)
          
          ! right multiply
          
          lhsK(1,k) = lhsK(1,k) * flowDiag(j,1)
          lhsK(4,k) = lhsK(4,k) * flowDiag(j,1)
          lhsK(7,k) = lhsK(7,k) * flowDiag(j,1)
          
          lhsK(2,k) = lhsK(2,k) * flowDiag(j,2)
          lhsK(5,k) = lhsK(5,k) * flowDiag(j,2)
          lhsK(8,k) = lhsK(8,k) * flowDiag(j,2)
          
          lhsK(3,k) = lhsK(3,k) * flowDiag(j,3)
          lhsK(6,k) = lhsK(6,k) * flowDiag(j,3)
          lhsK(9,k) = lhsK(9,k) * flowDiag(j,3)
          
          ! left multiply
          
          lhsK(1,k) = lhsK(1,k) * flowDiag(i,1)
          lhsK(2,k) = lhsK(2,k) * flowDiag(i,1)
          lhsK(3,k) = lhsK(3,k) * flowDiag(i,1)
          
          lhsK(4,k) = lhsK(4,k) * flowDiag(i,2)
          lhsK(5,k) = lhsK(5,k) * flowDiag(i,2)
          lhsK(6,k) = lhsK(6,k) * flowDiag(i,2)
          
          lhsK(7,k) = lhsK(7,k) * flowDiag(i,3)
          lhsK(8,k) = lhsK(8,k) * flowDiag(i,3)
          lhsK(9,k) = lhsK(9,k) * flowDiag(i,3)
          
        end do 
        
      end do
      
c
c.... ones in the diagonal positions
c      
      
      do n = 1, nshg
	  k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )
     &    + colm(n)-1
     
        lhsK(1,k) = one
        lhsK(5,k) = one
        lhsK(9,k) = one
        
      end do
            
      tempK = zero
      do i = 1, nshg
        do k = colm(i), colm(i+1)-1
        
          j = rowp(k)
          
          tempK(3*(i-1)+1,3*(j-1)+1) = lhsK(1,k)
          tempK(3*(i-1)+1,3*(j-1)+2) = lhsK(2,k)
          tempK(3*(i-1)+1,3*(j-1)+3) = lhsK(3,k)
          
          tempK(3*(i-1)+2,3*(j-1)+1) = lhsK(4,k)
          tempK(3*(i-1)+2,3*(j-1)+2) = lhsK(5,k)
          tempK(3*(i-1)+2,3*(j-1)+3) = lhsK(6,k)
          
          tempK(3*(i-1)+3,3*(j-1)+1) = lhsK(7,k)
          tempK(3*(i-1)+3,3*(j-1)+2) = lhsK(8,k)
          tempK(3*(i-1)+3,3*(j-1)+3) = lhsK(9,k)
            
        end do
      end do
      
!      open (unit = 2, file = "lhsK_pre.dat")
!      
!      do i=1,nshg*3
!        do j=1,nshg*3 
!          write(2,*) tempK(i,j)
!        enddo
!      enddo
!
!      close (2)
!      
!      open (unit = 2, file = "flowDiag.dat")
!      
!      do i=1,nshg
!        do j=1,3
!          write(2,*) flowDiag(i,j)
!        enddo
!      enddo
!
!      close (2)
      
c      
c.... create temporary right-hand side for testing
c     

      do i=1,nshg
        do j=1,3
          tempRHS(i,j) = 5.6812
        enddo
      enddo
      
      where (ibits(iBC,3,3) .eq. 7)  ! bits of iBC= xy111zab 
        tempRHS(:,1) = zero
        tempRHS(:,2) = zero
        tempRHS(:,3) = zero
      endwhere
      
      tempRHSPre = flowDiag * tempRHS

!      open (unit = 2, file = "RHS.dat")
!      
!      do i=1,nshg
!        do j=1,3
!          write(2,*) tempRHS(i,j)
!        enddo
!      enddo
!      
!      close(2)
      
      
!      open (unit = 2, file = "RHSPre.dat")
!      
!      do i=1,nshg
!        do j=1,3
!          write(2,*) tempRHSPre(i,j)
!        enddo
!      enddo
!      
!      close(2)
      
      call sparseCGK(tempRHSPre,trx,lhsK,rowp,colm,iper,ilwork,iBC,BC)
      
c
c.... recover solution by back-preconditioning
c

      where (flowDiag .ne. zero)  ! bits of iBC != xy111zab 
        trx = trx / flowDiag
      endwhere
           
      end subroutine
      
      
ccc==================================================================================
ccc==================================================================================

c INCOMPLETE ROUTINES

ccc==================================================================================
ccc================================================================================== 

      subroutine sparseCGK (rhsorig, trx, lhsM, row, col, iper,
     &                      ilwork,  iBC, BC)
c
c------------------------------------------------------------------------
c
c This subroutine uses Conjugate Gradient,
c to solve the system of equations.
c
c Farzin Shakib,  Summer 1987.
c
c Variant for 3by3 block matrix
c------------------------------------------------------------------------
c
      include "common.h"
c
      dimension rhsorig(nshg,3), trx(nshg,3)
c
      dimension d(nshg,3),     p(nshg,3),
     &          q(nshg,3),     ilwork(nlwork),
     &		    rhs(nshg,3),
     &          pp(nshg,3),    
     &          iBC(nshg),
     &          BC(nshg,ndofbc)


      integer   row(nnz_tot), col(nshg+1)
      integer   iBCdumb(nshg), iper(nshg)

      real*8    BCdumb(nshg,ndofBC)

      real*8	lhsM(9,nnz_tot)
c
      data nitercf / 100 /
c
c
      BCdumb  = one
      iBCdumb = 1
c
      rhs(:,:)=rhsorig(:,:)

c
c.... initialize
c      
      
      rr = 0
      
      do k = 1, 3
        do n = 1, nshg
          rr  = rr + rhs(n,k)*rhs(n,k)
        enddo
      enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
      dotTot=zero
      call drvAllreduce(rr,dotTot,1)
      rr=dotTot
      rr0 = rr
c
      trx = zero ! x_{0}=0
c                     ! r_{0}=b
c
c                       ! beta_{1}=0
      p = rhs ! p_{1}=r_{0}=b
c
c.... Iterate
c        
      do iter = 1, nitercf      ! 15  ! nitercf
c
c.... calculate alpha
c
         pp=p   ! make a vector that we can copy masters onto slaves
                ! and thereby make ready for an accurate Ap product

         call commOut(pp, ilwork, 3,iper,iBCdumb,BCdumb)  !slaves=master

         call fLesSparseApK( col, row, lhsM,
     &                       pp,  q,   nshg,
     &                       nnz_tot )

         call commIn(q, ilwork, 3, iper,iBC,BC) ! masters=masters+slaves
					                            ! slaves zeroed

         pap = 0
         
         do k = 1, 3
           do n = 1, nshg
             pap = pap + p(n,k)*q(n,k)
           enddo
         enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
         dotTot=zero
         call drvAllreduce(pap,dotTot,1)
         pap=dotTot
         alpha = rr / pap 
c
c.... calculate the next guess
c
         trx(:,:) = trx(:,:) + alpha * p(:,:)
c
c.... calculate the new residual
c
c
         rhs(:,:) = rhs(:,:) - alpha * q(:,:)
         tmp = rr
         rr = 0
         
         do k = 1, 3
           do n = 1, nshg
             rr = rr + rhs(n,k)*rhs(n,k)
           enddo
         enddo
c
c  if parallel the above is only this processors part of the dot product.
c  get the total result
c
         dotTot=zero
         call drvAllreduce(rr,dotTot,1)
         rr=dotTot
c
c.... check for convergence
c
         if(rr.lt.100.*epsM**2) goto 6000
c
c.... calculate a new search direction
c
         beta = rr / tmp
         p(:,:) = rhs(:,:) + beta * p(:,:)
c 
c.... end of iteration
c
      enddo
c
c.... if converged
c
6000  continue

c need a commu(out) on solution (TRX) to get slaves the correct solution AND
c on processor slaves = on processor masters

      if(numpe>1) call commu (trx, ilwork, 3, 'out ')
	
	trx(:,1) = trx(iper(:),1)
	trx(:,2) = trx(iper(:),2)
	trx(:,3) = trx(iper(:),3)

c
      write(*,9000) iter, rr / rr0
c	write(16,9000) iter, rr / rr0
c
c.... return
c
      return
9000  format(20x,'  number of iterations:', i10,/,
     &       20x,'    residual reduction:', 2x,e10.2)
     
      end subroutine
     
      
C============================================================================
C
C "fLesSparseApK":
C
C============================================================================

      subroutine fLesSparseApK( col, row, kLhs,
     1                          p,   q,   nNodes,
     2                          nnz_tot ) 
c
c.... Data declaration
c
      
      integer nNodes, nnz_tot
      integer col(nNodes+1), row(nnz_tot)
      real*8  kLhs(9,nnz_tot)
      real*8  p(nNodes,3), q(nNodes,3)
c
      real*8  tmp1, tmp2, tmp3, pisave
      integer i, j, k
c
c.... clear the vector
c
	do i = 1, nNodes
        q(i,1) = 0
        q(i,2) = 0
        q(i,3) = 0
	enddo
c
c.... Do an AP product
c
	do i = 1, nNodes
c
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        
c        pisave   = p(i,4)
cdir$ ivdep
        do k = col(i), col(i+1)-1
        
          j = row(k)
          
          tmp1 = tmp1
     1           + kLhs(1,k) * p(j,1)
     2           + kLhs(4,k) * p(j,2)
     3           + kLhs(7,k) * p(j,3)
     
          tmp2 = tmp2
     1           + kLhs(2,k) * p(j,1)
     2           + kLhs(5,k) * p(j,2)
     3           + kLhs(8,k) * p(j,3)
     
          tmp3 = tmp3
     1           + kLhs(3,k) * p(j,1)
     2           + kLhs(6,k) * p(j,2)
     3           + kLhs(9,k) * p(j,3)
c
        enddo
        
        q(i,1) = q(i,1) + tmp1
        q(i,2) = q(i,2) + tmp2
        q(i,3) = q(i,3) + tmp3
        
      enddo

c.... end
c
	return
	end subroutine
      
      
c      
c.... routine to do diagonal preconditioning
c.... specifically for the wall equations
c          
      subroutine drvlesPrepDiagK ( flowDiag, ilwork,
     &                             iBC,      BC,      iper,
     &                             rowp,     colm,    
     &                             lhsK)
c     
      use pointer_data

      include "common.h"
      include "mpif.h"
c     
      dimension flowDiag(nshg,3), ilwork(nlwork)
      dimension iBC(nshg), iper(nshg), BC(nshg,ndofBC)
      real*8    lhsK(9,nnz_tot)
      integer   rowp(nnz_tot), colm(nshg+1)
      integer   n, k
      integer   sparseloc
c
c.... Clear the flowdiag
c
      flowDiag = zero
	do n = 1, nshg
	  k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )
     &    + colm(n)-1
c     
	  flowdiag(n,1) = lhsK(1,k)
	  flowdiag(n,2) = lhsK(5,k)
	  flowdiag(n,3) = lhsK(9,k)

      enddo

c     
c.... communicate : add the slaves part to the master's part of flowDiag
c     
	if (numpe > 1) then 
        call commu (flowDiag, ilwork, 3, 'in ') 
      endif
c
c.... satisfy the boundary conditions on the diagonal
c.... assuming all three components are zeroed at the dirichlet nodes
c
      where (ibits(iBC,3,3) .eq. 7)  ! bits of iBC= xy111zab 
        flowDiag(:,1) = zero
        flowDiag(:,2) = zero
        flowDiag(:,3) = zero
      endwhere
c
c     
c.... on processor periodicity was not taken care of in the setting of the 
c     boundary conditions on the matrix.  Take care of it now.
c
      call bc3per(iBC,  flowDiag, iper, ilwork, 3)
c
c... slaves and masters have the correct values
c
c     
c.... Calculate square root
c     
      do i = 1, nshg
        do j = 1, 3
          if (flowDiag(i,j).ne.0) 
     &      flowDiag(i,j) = 1. / sqrt(abs(flowDiag(i,j)))
        enddo
      enddo
c     
      return
      end subroutine
      
      
      
      
ccc==================================================================================
ccc==================================================================================
ccc==================================================================================

c
c fillsparse routine designed for the deformable wall      
c      
      subroutine fillsparseK(iens, xKebe, lhsK,
     &                       row,  col)
c
c
c
	include "common.h"
	
	real*8  xKebe(npro,9,nshl,nshl)
	integer ien(npro,nshl), col(nshg+1), row(nshg*nnz)
	
	integer iBCB(npro,ndiBCB)
	
	real*8  lhsK(9,nnz_tot)
c
      integer aa, b, c, e, i, k, n, j, l
c
      integer sparseloc

      integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
      ien=abs(iens)
c       
c.... Accumulate the lhs
c
do e = 1, npro ! loop over the elements
      do e = 1, npro ! loop over the elements
        do aa = 1, nshl ! loop over the local equation numbers
          i = ien(e,aa) ! finds the global equation number or
                          ! block-row of our matrix
          c = col(i)    ! starting point to look for the matching column
          n = col(i+1) - c  !length of the list of entries in rowp
          
          do b = 1, nshl ! local variable number tangent respect
                           ! to
c function that searches row until it finds the match that gives the
c		   global equation number

            k = sparseloc(row(c),n,ien(e,b))+c-1

            lhsK(1,k) = lhsK(1,k) + xKebe(e,1,aa,b)
            lhsK(2,k) = lhsK(2,k) + xKebe(e,2,aa,b)
            lhsK(3,k) = lhsK(3,k) + xKebe(e,3,aa,b)
            lhsK(4,k) = lhsK(4,k) + xKebe(e,4,aa,b)
            lhsK(5,k) = lhsK(5,k) + xKebe(e,5,aa,b)
            lhsK(6,k) = lhsK(6,k) + xKebe(e,6,aa,b)
            lhsK(7,k) = lhsK(7,k) + xKebe(e,7,aa,b)
            lhsK(8,k) = lhsK(8,k) + xKebe(e,8,aa,b)
            lhsK(9,k) = lhsK(9,k) + xKebe(e,9,aa,b)

          enddo
        enddo
      enddo

c
c.... end
c
	return
	end subroutine
	
	
cc ============================================================================
cc UNUSED for now

c
c.... just like asadj, except for the ''deformable wall'
c.... part of the boundary
c      

      subroutine AsadjWall ( row_fill_list,
     &                       iens, adjcnt,
     &                       iBCB )
     
      use pointer_data

      include "common.h"
 
c
      integer row_fill_list(nwnp,6*nnz),
     &        ien(npro,nshlb),
     &        adjcnt(nwnp), ndlist(nshlb)
     
      integer iBCB(npro,ndiBCB)

      integer iens(npro,nshl)
c
c prefer to show explicit absolute value needed for cubic modes and
c higher rather than inline abs on pointer as in past versions
c iens is the signed ien array ien is unsigned
c
      ien=abs(iens)
	  
      do i=1,npro
c
c.... turbulence wall (as a way of checking for deformable wall stiffness)
c.... this is necessary due to the compression mapping
c      
         if (btest(iBCB(i,1),4)) then               
            do j=1,nshlb
               ndlist(j)=mWNodes_gtlmap%p(ien(i,j))
            enddo
            do j=1,nshlb
               jnd=ndlist(j)
               jlngth=adjcnt(jnd) ! current length of j's list
               do k=1,nshlb
                  knd=ndlist(k)
                  ibroke=zero
                  do l= 1,jlngth
                     if(row_fill_list(jnd,l).eq. knd) then
                        ibroke=1
                        exit
                     endif
                  enddo
                 
c
c.... to get here k was not in  j's list so add it
c
                  if(ibroke.eq.0) then
                     jlngth=jlngth+1 ! lenthen list
                     if(jlngth.gt.6*nnz) then
                        write(*,*) 'increase overflow factor in genadj'
                        stop
                     endif
                     row_fill_list(jnd,jlngth)=knd ! add unique entry to list
                  endif
               enddo ! finished checking all the k's for this j
               adjcnt(jnd)=jlngth  ! update the counter
            enddo                  ! done with j's
         end if
      enddo                   ! done with elements in this block
c
c
c.... end
c
      return
      
      end subroutine
      

            
c     
c.... just like genadj, instead we generate
c.... the sparse matrix pattern for the 'deformable wall' part
c.... of the boundary
c  	  
	subroutine genadjWall ( colm, rowp, icnt )

      use pointer_data
c     
      include "common.h"
c     
c     we use nwnp in place of nshg
      integer rowp(nwnp*nnz), colm(nwnp+1)
      integer adjcnt(nwnp),   row_fill_list(nwnp,6*nnz), mloc(1)
c                                           change^ if overflow

      integer tmprdim(1)
      real*8, allocatable, dimension(:) :: tmpr

      adjcnt=0
      
      do iblk = 1, nelblb
c     
c.... set up the parameters
c     
       
         iel    = lcblkb(1,iblk)
         lelCat = lcblkb(2,iblk)
         lcsyst = lcblkb(3,iblk)
         iorder = lcblkb(4,iblk)
         nenl   = lcblkb(5,iblk)  ! no. of vertices per element
         nenbl  = lcblkb(6,iblk)  ! no. of vertices per bdry. face
         nshl   = lcblkb(9,iblk)
         nshlb  = lcblkb(10,iblk)
         mattyp = lcblkb(7,iblk)
         ndofl  = lcblkb(8,iblk)
         npro   = lcblkb(1,iblk+1) - iel
         
c     
c.... compute sparse matrix data structures
c     
         
         call AsadjWall ( row_fill_list,                       
     &                    mienb(iblk)%p,  adjcnt,
     &                    miBCB(iblk)%p )
     
      enddo
     
      call sumgatInt ( adjcnt, nwnp, nnonzero )
         
      if ( myrank .eq. master) then
         write (*,*) 'Number of global wall nonzeros ',nnonzero
      endif
                  
c     
c     build the colm array
c     
      colm(1)=1
      do i=1,nwnp
         colm(i+1)=colm(i)+adjcnt(i)
      enddo
c     
c     sort the rowp into increasing order
c     
      ibig=10*nwnp
      icnt=0
      tmprdim=maxval(adjcnt)
      allocate (tmpr(tmprdim(1)))
      do i=1,nwnp
         ncol=adjcnt(i)
         tmpr(1:ncol)=row_fill_list(i,1:ncol)
         do j=1,ncol
            icnt=icnt+1
            imin=minval(tmpr(1:ncol))
            mloc=minloc(tmpr(1:ncol))
            rowp(icnt)=imin
            tmpr(mloc(1))=ibig
         enddo
      enddo

      maxfill=tmprdim(1)
      write(*,*) 'maxfill=',maxfill
      nnza=icnt/nwnp +1
      if(icnt.gt.nnz*nwnp) then
         write(*,*) 'increase nnz in genmat to',nnza
         stop
      else
         write(*,*) 'nnz ok  nnz=',nnz,' actually needed',nnza   
         write(*,*) myrank,' is my rank and my nnz_tot_wall is: ',icnt
      endif
      return
      
      end subroutine
      

      
      end module