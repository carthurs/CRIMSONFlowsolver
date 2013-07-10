module deformableWall
     
    parameter ( MAXBLK2 = 500000 )

    type r1d
        real*8, pointer :: p(:)
    end type

    type r2d
        real*8, pointer :: p(:,:)
    end type

    type r3d
        real*8, pointer :: p(:,:,:)
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
    type (r2d), dimension(MAXBLK2) ::  mSWB
      
    ! currently allocated in gendat, but should be moved here
    type (i1d) :: mWNodes, mWNodes_gtlmap

    ! these are dynamically allocated in a local routine
    ! type (r2d), dimension(MAXBLK2) :: mPS_global
    ! type (r4d), dimension(MAXBLK2) :: mKwall_xKebe


    ! the reference element local displacements
    type (r3d), dimension(MAXBLK2) :: mDisp_refl



    ! element local node tags
    type (i2d), dimension(MAXBLK2) :: mNodeTagl
      
contains
      
      
      
    subroutine vlmwStTri(u,nodetagfield,x)
      
        use pointer_data
        use measureWallDistance

        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

        real*8 u(nshg,nsd),x(numnp,nsd)
        integer nodetagfield(numnp)



        !
        !.... initialize distance evaluation (maybe move this to itrdrv)
        !
        if(ideformwall.eq.1 .and. imeasdist.eq.1) then
            ! read "observed mesh data"
            write(*,*) "reading observed mesh at rank ",myrank
            call dm_dataread("observed.dat")

            ! build octrees for observed meshes
            write(*,*) "building octrees at rank ",myrank
            call dm_initialize()
        end if

!        do i = 1, numnp
!            if (nodetagfield(i).gt.0) then
!                write(*,*) nodetagfield(i)
!            end if
!        end do

            
        !
        !.... loop over the boundary elements
        !
        do iblk = 1, nelblb
            !
            !.... set up the parameters
            !
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
         
            !
            !.... allocate space for the arrays
            !
    
!            if (.not.associated(mPS_global(iblk)%p)) then
!                allocate(mPS_global(iblk)%p(npro,9))
!            end if

!            if (.not.associated(mKwall_xKebe(iblk)%p)) then
!                allocate(mKwall_xKebe(iblk)%p(npro,9,nshl,nshl))
!            end if

            if (.not.associated(mDisp_refl(iblk)%p)) then
                allocate(mDisp_refl(iblk)%p(npro,nshl,nsd))
            end if

            !
            !.... check if SWB files are used
            !.... if no, then use the entered uniform values
            !

            !if (iUseSWB.eq.0) then
            !    mSWB(iblk)%p = zero
            !    mSWB(iblk)%p(:,1) = thicknessvw
            !end if

            ! save the initial element local displacements
            mDisp_refl(iblk)%p = zero

            ! add a flag here to say if we are using initial prestress
            if (iinitialprestress.gt.0) then
                call localx(u, mDisp_refl(iblk)%p, mienb(iblk)%p, nsd, 'gather  ')
            end if

            ! check if element nodes contain more than one ring node

            if (iwallsupp.gt.0 .or. iwalldamp.gt.0) then

                if (.not.associated(mNodeTagl(iblk)%p)) then
                    allocate(mNodeTagl(iblk)%p(npro,nenl))
                end if

                call locali(nodetagfield, mNodeTagl(iblk)%p, mienb(iblk)%p, 1, 'gather  ')

            end if

            !
            !.... pre-assemble the element matrices for the deformable wall
            !
             
!            call asbvlmwsttri(x,    &
!            mien(iblk)%p,         iBC,           BC, &
!            mienb(iblk)%p,        miBCB(iblk)%p, &
!            mSWB(iblk)%p, &
!            mPS_global(iblk)%p, &
!            mKwall_xKebe(iblk)%p)

        end do
            
    end subroutine
      
      
      
!    subroutine asbvlmwsttri(x,            &
!    ien,         iBC,   BC, &
!    ienb,        iBCB,   &
!    SWB, &
!    PS_global, &
!    Kwall_xKebe)
!
!        use phcommonvars
!        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
!        real*8    x(numnp,nsd)
!
!        integer   iBCB(npro,ndiBCB),         ienb(npro,nshl)
!        integer   ien(npro,nshape)
!
!        real*8    SWB(npro,nProps)
!
!        real*8    xlb(npro,nenl,nsd)
!
!        real*8    rotel(npro,9,9), &
!                  temp(npro), &
!                  temp1(npro),               temp2(npro), &
!                  temp3(npro), &
!                  v1(npro,nsd),              v2(npro,nsd), &
!                  v3(npro,nsd),               &
!                  x1rot(npro,nsd),            &
!                  x2rot(npro,nsd),           x3rot(npro,nsd),            &
!                  detJacrot(npro), &
!                  B(npro,5,9),               Dmatrix(npro,5,5), &
!                  DtimesB(npro,5,9),          &
!                  Kwall_local(npro,9,9),     Kwall_global(npro,9,9), &
!                  Kwall_temp(npro,9,9), &
!                  PS_local(npro,9),          PS_global(npro,9), &
!                  Kwall_xKebe(npro,9,nshl,nshl)
!
!        real*8    BC(nshg,ndofBC)
!        integer   iBC(nshg)
!
!        integer   i,j,k,l,m,n
!
!        !
!        !.... gather the variables
!        !
!
!        call localx(x,xlb,ienb,nsd,'gather  ')
!
!        !
!        !.... check for element coordinate system
!        !
!
!        if (lcsyst.eq.1) then
!            ipt2=2
!            ipt3=3
!        else
!            write(*,*) "triangular elements required"
!            stop
!        end if
!
!        !
!        !.... Assemble the stiffness LHS matrix and residual contribution matrix
!        !
!
!        !
!        !.... B^t * D * B formulation for plane stress enhanced membrane
!        !
!
!        !
!        !.... rotation matrix
!        !
!        v1 = xlb(:,ipt2,:) - xlb(:,1,:)
!        temp = one / sqrt ( v1(:,1)**2 + v1(:,2)**2 + v1(:,3)**2 )
!        v1(:,1) = v1(:,1) * temp
!        v1(:,2) = v1(:,2) * temp
!        v1(:,3) = v1(:,3) * temp
!
!        v2 = xlb(:,ipt3,:) - xlb(:,1,:)
!
!        !     compute cross product
!        temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
!        temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
!        temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
!
!        temp     = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
!        v3(:,1) = temp1 * temp
!        v3(:,2) = temp2 * temp
!        v3(:,3) = temp3 * temp
!
!        !     cross product again for v2
!        temp1 = v3(:,2) * v1(:,3) - v1(:,2) * v3(:,3)
!        temp2 = v1(:,1) * v3(:,3) - v3(:,1) * v1(:,3)
!        temp3 = v3(:,1) * v1(:,2) - v1(:,1) * v3(:,2)
!
!        temp     = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
!        v2(:,1) = temp1 * temp
!        v2(:,2) = temp2 * temp
!        v2(:,3) = temp3 * temp
!
!        rotel = zero
!
!        rotel(:,1,1:3) = v1(:,:)
!        rotel(:,2,1:3) = v2(:,:)
!        rotel(:,3,1:3) = v3(:,:)
!
!        rotel(:,4,4:6) = v1(:,:)
!        rotel(:,5,4:6) = v2(:,:)
!        rotel(:,6,4:6) = v3(:,:)
!
!        rotel(:,7,7:9) = v1(:,:)
!        rotel(:,8,7:9) = v2(:,:)
!        rotel(:,9,7:9) = v3(:,:)
!
!        !
!        !.... rotated coordinates
!        !
!        x1rot = zero
!        x2rot = zero
!        x3rot = zero
!
!        do i = 1, 3
!            do j = 1, 3
!                x1rot(:,i) = x1rot(:,i)+rotel(:,i,j)*xlb(:,1,j)
!                x2rot(:,i) = x2rot(:,i)+rotel(:,i,j)*xlb(:,ipt2,j)
!                x3rot(:,i) = x3rot(:,i)+rotel(:,i,j)*xlb(:,ipt3,j)
!            end do
!        end do
!
!        !
!        !.... B matrices
!        !
!        B = zero
!
!        detJacrot = (x2rot(:,1)-x1rot(:,1)) * (x3rot(:,2)-x1rot(:,2)) -  &
!        (x3rot(:,1)-x1rot(:,1)) * (x2rot(:,2)-x1rot(:,2))
!
!        B(:,1,1) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
!        B(:,2,2) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
!        B(:,3,1) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
!        B(:,3,2) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
!        B(:,4,3) = (x2rot(:,2)-x3rot(:,2))/detJacrot(:)
!        B(:,5,3) = (x3rot(:,1)-x2rot(:,1))/detJacrot(:)
!
!        B(:,1,4) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
!        B(:,2,5) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
!        B(:,3,4) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
!        B(:,3,5) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
!        B(:,4,6) = (x3rot(:,2)-x1rot(:,2))/detJacrot(:)
!        B(:,5,6) = (x1rot(:,1)-x3rot(:,1))/detJacrot(:)
!
!        B(:,1,7) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
!        B(:,2,8) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
!        B(:,3,7) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
!        B(:,3,8) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
!        B(:,4,9) = (x1rot(:,2)-x2rot(:,2))/detJacrot(:)
!        B(:,5,9) = (x2rot(:,1)-x1rot(:,1))/detJacrot(:)
!
!        !
!        !.... D matrix
!        !
!        Dmatrix = zero
!
!        if (iUseSWB.gt.0) then
!            if (nProps.eq.10) then
!                !        This is an Isotropic Material
!
!                Dmatrix(:,1,1) = SWB(:,7)
!                Dmatrix(:,2,2) = SWB(:,7)
!                Dmatrix(:,1,2) = SWB(:,8)
!                Dmatrix(:,2,1) = SWB(:,8)
!                Dmatrix(:,3,3) = SWB(:,9)
!                Dmatrix(:,4,4) = SWB(:,10)
!                Dmatrix(:,5,5) = SWB(:,10)
!
!            elseif(nProps.eq.21) then
!
!                !        This is an Orthotropic Material
!                Dmatrix(:,1,1) = SWB(:,7)
!
!                Dmatrix(:,2,1) = SWB(:,8)
!                Dmatrix(:,1,2) = SWB(:,8)
!                Dmatrix(:,2,2) = SWB(:,9)
!
!                Dmatrix(:,3,1) = SWB(:,10)
!                Dmatrix(:,1,3) = SWB(:,10)
!                Dmatrix(:,3,2) = SWB(:,11)
!                Dmatrix(:,2,3) = SWB(:,11)
!                Dmatrix(:,3,3) = SWB(:,12)
!
!                Dmatrix(:,4,1) = SWB(:,13)
!                Dmatrix(:,1,4) = SWB(:,13)
!                Dmatrix(:,4,2) = SWB(:,14)
!                Dmatrix(:,2,4) = SWB(:,14)
!                Dmatrix(:,4,3) = SWB(:,15)
!                Dmatrix(:,3,4) = SWB(:,15)
!                Dmatrix(:,4,4) = SWB(:,16)
!
!                Dmatrix(:,5,1) = SWB(:,17)
!                Dmatrix(:,1,5) = SWB(:,17)
!                Dmatrix(:,5,2) = SWB(:,18)
!                Dmatrix(:,2,5) = SWB(:,18)
!                Dmatrix(:,5,3) = SWB(:,19)
!                Dmatrix(:,3,5) = SWB(:,19)
!                Dmatrix(:,5,4) = SWB(:,20)
!                Dmatrix(:,4,5) = SWB(:,20)
!                Dmatrix(:,5,5) = SWB(:,21)
!
!            else
!
!                write(*,*) 'Number of wall properties not set correctly!'
!                stop
!
!            end if
!        else
!
!            Dmatrix(:,1,1) = one
!            Dmatrix(:,2,2) = Dmatrix(:,1,1)
!            Dmatrix(:,1,2) = rnuvw
!            Dmatrix(:,2,1) = Dmatrix(:,1,2)
!            Dmatrix(:,3,3) = pt5*(1-rnuvw)
!            Dmatrix(:,4,4) = pt5*(1-rnuvw)*rshearconstantvw
!            Dmatrix(:,5,5) = pt5*(1-rnuvw)*rshearconstantvw
!
!            Dmatrix(:,:,:) = one/(one-rnuvw**2)*Dmatrix(:,:,:);
!
!        end if
!
!        !
!        !.... D * [B1|B2|B3]
!        !
!        DtimesB = zero
!
!        do i = 1, 5
!            do j = 1, 9
!                do k = 1, 5
!                    DtimesB(:,i,j)=DtimesB(:,i,j)+Dmatrix(:,i,k)*B(:,k,j)
!                end do
!            end do
!        end do
!
!        !
!        !.... [B1|B2|B3]^T * D * [B1|B2|B3]
!        !
!        Kwall_local = zero
!
!        do i = 1, 9
!            do j = 1, 9
!                do k = 1, 5
!                    Kwall_local(:,i,j)= &
!                    Kwall_local(:,i,j)+B(:,k,i)*DtimesB(:,k,j)
!                end do
!            end do
!        end do
!
!
!        !
!        !.... [B1|B2|B3]^T * Sigma  Compute prestress contribution
!        !
!
!        if (iUseSWB.gt.0) then
!            do i = 1, 3            ! BEWARE index 3 and 4
!                PS_local(:,i) = B(:,1,i) * SWB(:,2) + &
!                B(:,2,i) * SWB(:,4) + &
!                B(:,3,i) * SWB(:,3) + &
!                B(:,4,i) * SWB(:,5) + &
!                B(:,5,i) * SWB(:,6)
!
!                PS_local(:,i+3) = B(:,1,i+3) * SWB(:,2) + &
!                B(:,2,i+3) * SWB(:,4) + &
!                B(:,3,i+3) * SWB(:,3) + &
!                B(:,4,i+3) * SWB(:,5) + &
!                B(:,5,i+3) * SWB(:,6)
!
!                PS_local(:,i+6) = B(:,1,i+6) * SWB(:,2) + &
!                B(:,2,i+6) * SWB(:,4) + &
!                B(:,3,i+6) * SWB(:,3) + &
!                B(:,4,i+6) * SWB(:,5) + &
!                B(:,5,i+6) * SWB(:,6)
!            end do
!        end if
!
!        !
!        !.... rotate element stiffness to global coordinates
!        !
!
!        Kwall_global = zero
!        Kwall_temp = zero
!
!        do i = 1, 9
!            do j = 1, 9
!                do k= 1, 9
!
!                    Kwall_temp(:,i,j) = Kwall_temp(:,i,j) +  &
!                    Kwall_local(:,i,k) * rotel(:,k,j)
!
!                end do
!            end do
!        end do
!
!        do i = 1, 9
!            do j = 1, 9
!                do k= 1, 9
!
!                    Kwall_global(:,i,j) = Kwall_global(:,i,j) +  &
!                    rotel(:,k,i) * Kwall_temp(:,k,j)
!
!                end do
!            end do
!        end do
!
!        !
!        !.... rotate local prestress contribution to global coordinates
!        !
!        PS_global = zero
!
!        if (iUseSWB.gt.0) then
!            do i =1, 9
!                do j = 1, 9
!                    PS_global(:,i) = rotel(:,j,i) * PS_local(:,j) + &
!                    PS_global(:,i)
!
!                end do
!            end do
!        end if
!
!        !
!        !.... integrate over element - multiply the nodal matrices by the area and the thickness
!        !
!
!        do i =1, 9
!            do j = 1, 9
!
!                Kwall_global(:,i,j) = Kwall_global(:,i,j)*detJacrot(:) &
!                *pt5*SWB(:,1)
!
!                ! hack to remove the elastic modulus factor (uniform value for now) at this point
!                ! when using the SWB option
!                ! since it will be multiplied later in e3b
!                ! for now we only use SWB for the prestress (again uniform)
!
!                if (iUseSWB.gt.0) then
!                    Kwall_global(:,i,j) = Kwall_global(:,i,j)/evw
!                end if
!
!            enddo
!
!            PS_global(:,i) = PS_global(:,i)*detJacrot(:)*pt5*SWB(:,1)
!
!        enddo
!
!        !
!        !.... collect the entries of the stiffness matrix into the xKebe format
!        !
!
!        Kwall_xKebe = zero
!
!        do i=1,3
!            do j=1,3
!                do k=1,3
!                    do m=1,3
!                        Kwall_xKebe(:,3*(k-1)+m,i,j) =  &
!                        Kwall_global(:,3*(i-1)+k,3*(j-1)+m)
!                    end do
!                end do
!            end do
!        end do
!
!        !
!        !.... boundary code (as a way of checking for which boundary elements have deformable wall stiffness)
!        !
!
!        do iel = 1, npro
!            if (btest(iBCB(iel,1),4)) then
!
!            else
!
!                Kwall_xKebe(iel,:,:,:) = zero
!            end if
!        end do
!
!        !
!        !.... assemble lhsK with the global wall equation numbers
!        !.... this step is for the state filter
!        !
!
!        !      call fillsparseW(ienb, Kwall_xKebe, Kwall_lhsK,
!        !     &                 row,   col)
!
!        !
!        !.... apply dirichlet conditions for components
!        !.... this step is for the state filter
!        !
!
!        !call bc3lhs(iBC, BC, ienb, Kwall_xKebe)
!
!    end subroutine
!
end module
      
      
!      
!.... deallocated module variables
!            
subroutine DdeformableWall
      
    use deformableWall
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      
      
    if (ideformwall.ne.0) then
        do iblk = 1, nelblb
            
            if (associated(mSWB(iblk)%p)) deallocate (mSWB(iblk)%p)

            if (associated(mDisp_refl(iblk)%p)) deallocate (mDisp_refl(iblk)%p)
            if (associated(mNodeTagl(iblk)%p)) deallocate (mNodeTagl(iblk)%p)
      
        end do
    end if
      
end subroutine
