subroutine ElmDist(u, x, xdist, xdnv, df_fem) bind(C, name="elmdist")
      
    use pointer_data  ! brings in the pointers for the blocked arrays
    use measureWallDistance
    use deformableWall
    use shapeTable
         
    !     
    ! the subroutine interface for ElmDist is defined inside phcommonvars
    ! however the fortran standard doesn't permit the subroutine to know 
    ! it's own interface, see:
    !
    ! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51268
    ! https://groups.google.com/forum/#!topic/comp.lang.fortran/7zkAqyYcXNk
    ! 
    ! the fix is to reassign to rename the function, inside the function 
    ! itself, i.e. use MODULE_NAME, LOCAL_NAME => MODULE_NAME 
    ! 
    ! KDL, June 2016

    use phcommonvars, temp => elmdist

    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
         
    real*8 x(numnp,NSD),u(nshg,NSD)
    real*8 xdist(nshg)
    real*8 xdnv(numnp,NSD)
    real*8 df_fem(nshg)
    !real*8 shpb(MAXTOP,MAXSH,MAXQPT), shglb(MAXTOP,NSD,MAXSH,MAXQPT)

    real*8 tempnv1(NSD)
    real*8 tempnv2(NSD)
         
    real*8 closestPt1(NSD),closestPt2(NSD)
         
    real*8 tempPt(NSD)
    real*8 tempDistSq1,tempDistSq2,cycleTime,obsInterval,intTime
    real*8 alphaObs
         
    real*8, allocatable :: tmpshpb(:,:), tmpshglb(:,:,:)

    real*8 df_fem_sum

    integer ii,jj,nn,nPer,nObsInt,obsFr1,obsFr2
         
    !
    !.....compute the number of periods that has passed
    !     according to the cycle length parameter that was set
    !     in observed.dat
    !
    nPer = (lstep)*Delt(1)/cycleLength
         
    !
    !.... compute the physical time
    !
                  
    cycleTime = (lstep)*Delt(1)-nPer*cycleLength
         
    !
    !.... compute the interval that we are in right now
    !

    obsInterval = cycleLength / numDataFrames ! interval between data
    nObsInt = cycleTime / obsInterval
    obsFr1 = nObsInt + 1
    obsFr2 = nObsInt + 2
    intTime = cycleTime - nObsInt*obsInterval
         
    if (obsFr2 .gt. numDataFrames) then
        obsFr2 = 1
    end if
         
    alphaObs = 1-intTime/obsInterval
         
    if (myrank .eq. 0) then
        write(*,*) "dtime",obsFr1,obsFr2,(1-alphaObs)
    end if

    !
    !        loop over the deformable wall nodes
    !

    do ii = 1, nwnp

        ! TODO: assumption is that numnp = nshg
        tempPt = x(mWNodes%p(ii),:)+ &
            u(mWNodes%p(ii),:)

        call dm_cpmeshp3(obsFr1,tempPt, &
            closestpt1,tempDistSq1, &
            tempnv1)
        call dm_cpmeshp3(obsFr2,tempPt, &
            closestpt2,tempDistSq2, &
            tempnv2)

        tempDistSq1 = sign(sqrt(abs(tempDistSq1)),tempDistSq1)

        tempDistSq2 = sign(sqrt(abs(tempDistSq2)),tempDistSq2)

        xdist(mWNodes%p(ii)) =  &
            alphaObs*tempDistSq1+(1-alphaObs)*tempDistSq2

        xdnv(mWNodes%p(ii),:) =  &
            alphaObs*tempnv1+(1-alphaObs)*tempnv2

    end do

    !
    !        loop through the boundary elements
    !

    df_fem = zero

    ! testing
    !xdist = zero
    ! testing
    !xdist(mWNodes%p) = one
    !write(*,*) sum(xdist)

    do iblk = 1, nelblb

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

        !ngauss = nint(lcsyst)

        icurrentblk = iblk  ! current block

        allocate (tmpshpb(nshl,MAXQPT))
        allocate (tmpshglb(nsd,nshl,MAXQPT))

        tmpshpb(1:nshl,:) = shpb(lcsyst,1:nshl,:)
        tmpshglb(:,1:nshl,:) = shglb(lcsyst,:,1:nshl,:)

        call AsBDist(x, xdist, xdnv, tmpshpb, tmpshglb, mienb(iblk)%p, df_fem)

        if (allocated(tmpshpb)) then
          deallocate(tmpshpb)
        endif
        if (allocated(tmpshglb)) then
          deallocate(tmpshglb)
        endif

    end do

    ! testing
!    df_fem_sum = zero
!    do n = 1,nshg
!    df_fem_sum = df_fem_sum + df_fem(n)
!    enddo
!
!    write(*,*) df_fem_sum
                  
end subroutine

subroutine AsBDist(x, xdist, xdnv, shpb, shglb, ienb, df_fem)

    use phcommonvars
    use deformableWall
    use pointer_data
    use measureWallDistance

    real*8 x(numnp,NSD)
    real*8 xdist(nshg), xdnv(nshg,nsd)
    real*8 shpb(nshl,ngaussb), shglb(nsd,nshl,ngaussb)

    integer ienb(npro,nshl)

    real*8 sgn(npro,nshl)

    real*8 xlb(npro,nenl,nsd)
    real*8 xdistl(npro,nshl), xdnvl(npro,nshl,nsd)
    real*8 shapeVar(npro,nshl), shdrv(npro,nsd,nshl)

    real*8 bnorm(npro,nsd)
    real*8 temp(npro)
    real*8 temp1(npro), temp2(npro), temp3(npro)
    real*8 v1(npro,nsd), v2(npro,nsd), v3(npro,nsd)
    real*8 WdetJb(npro)

    real*8 usup2(npro)

    real*8 dfl(npro,nshl), dfl2(npro,nshl)
    real*8 df_fem(nshg)

    integer lnode(27)

    !
    !.... get the matrix of mode signs for the hierarchic basis functions
    !
    if (ipord .gt. 1) then
        call getsgn(ienb,sgn)
    endif

    !
    !.... gather the variables
    !
    call localx(x, xlb, ienb, NSD, 'gather  ')
    call local(xdist, xdistl, ienb, NSD, 'gather  ')
    call local(xdnv, xdnvl, ienb, NSD, 'gather  ')

    !
    !.... compute the nodes which lie on the boundary (hierarchic)
    !
    call getbnodes(lnode)

    !
    !.... loop through the integration points
    !
    if(lcsyst.eq.3) lcsyst=nenbl
    !
    if(lcsyst.eq.3 .or. lcsyst.eq.4) then
        ngaussb = nintb(lcsyst)
    else
        ngaussb = nintb(lcsyst)
    endif

    if(lcsyst.eq.1) then      ! set to curl into element all others out
        ipt2=2
        ipt3=3
    elseif(lcsyst.eq.2) then
        ipt2=4
        ipt3=2
    elseif(lcsyst.eq.3) then
        ipt2=3
        ipt3=2
    elseif(lcsyst.eq.4) then
        ipt2=2
        ipt3=4
    elseif(lcsyst.eq.5) then
        ipt2=4
        ipt3=2
    elseif(lcsyst.eq.6) then
        ipt2=2
        ipt3=5
    endif

    v1 = xlb(:,ipt2,:) - xlb(:,1,:)
    v2 = xlb(:,ipt3,:) - xlb(:,1,:)
    !
    ! compute cross product
    !
    temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
    temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
    temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
    !
    ! mag is area for quads, twice area for tris
    !
    temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
    bnorm(:,1) = temp1 * temp
    bnorm(:,2) = temp2 * temp
    bnorm(:,3) = temp3 * temp

    dfl  = zero

    do intp = 1, ngaussb ! loop through quadrature points

        !
        ! .... Weighted Jacobian Determinant
        !
        if (lcsyst .eq. 1) then
            WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
        elseif (lcsyst .eq. 2) then
            WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
        elseif (lcsyst .eq. 3) then
            WdetJb     = Qwtb(lcsyst,intp) / (two*temp)
        elseif (lcsyst .eq. 4) then
            WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
        elseif (lcsyst .eq. 5) then
            WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
        elseif (lcsyst .eq. 6) then
            WdetJb     = Qwtb(lcsyst,intp) / (two*temp)
        endif

        !
        !.... get the hierarchic shapeVar functions at this int point
        !
        call getshpb(shpb, shglb, sgn, shapeVar, shdrv)

        usup2 = zero

        do n = 1, nshlb

            nodlcl = lnode(n)

            usup2(:) = usup2(:) + shapeVar(:,nodlcl) * xdistl(:,nodlcl)

        enddo

        usup2(:) = usup2(:)*WdetJb   ! multiply by the weighted Jacobian

        do n = 1, nshlb
            nodlcl = lnode(n)

            dfl(:,nodlcl) = dfl(:,nodlcl) + shapeVar(:,nodlcl) * usup2(:)

        enddo

    end do

    do n = 1, nshlb
        nodlcl = lnode(n)
        dfl(:,nodlcl) = dfl(:,nodlcl) * temp(:) * two / numwallelems_global ! this works for triangles
    enddo

    do i=1,npro
        if (.not.btest(miBCB(icurrentblk)%p(i,1),4)) then ! check element deformable
            dfl(i,:) = zero
        end if
    end do


    call local (df_fem, dfl, ienb, 1, 'scatter ')

end subroutine
