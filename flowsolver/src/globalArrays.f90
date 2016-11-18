!    module globalArrays contains the arrays that formerly were passed
!    into (or initialized in) itrdrv but now must be reformed on every 
!    adaptivity step by a call
!      to process (note some will appear to be missing because they are
!      already part of other modules that will also be included in itrdrv


module globalArrays
    use, intrinsic :: iso_c_binding
    integer, target, allocatable :: nodetagfield(:)
    real (c_double), target, allocatable :: y(:,:)
    real (c_double), target, allocatable :: ac(:,:)
    real (c_double), target, allocatable :: u(:,:)
    real (c_double), target, allocatable :: yold(:,:)
    real (c_double), target, allocatable :: acold(:,:)
    real (c_double), target, allocatable :: uold(:,:)
    real*8, allocatable :: uref(:,:)
    real*8, allocatable :: ubar(:,:)

    ! gfortran is throwing an error because xdist, xdnv, df_fem  doesn't have
    ! the TARGET attribute - KDL 14/06/2016
    real*8, target, allocatable :: xdist(:)
    real*8, target, allocatable :: xdnv(:,:)
    real*8, target, allocatable :: df_fem(:)

    ! gfortran is throwing an error because temporary_array  doesn't have the 
    ! TARGET attribute - KDL 14/06/2016
    real (c_double), target, allocatable :: temporary_array(:,:)

    integer (c_int), target, allocatable :: iBC(:)
    integer (c_int), target, allocatable :: iBC_original(:)
    real*8, allocatable :: BC(:,:)
    real*8, allocatable :: solinc(:,:)
    real*8, allocatable :: res(:,:)
    integer, allocatable :: rowp(:,:)
    integer, allocatable :: colm(:)
    real*8, allocatable :: rerr(:,:)
    real*8, allocatable :: ybar(:,:)
    real*8, allocatable :: dummyVar(:)
    real*8, allocatable :: uhess(:,:)
    real*8, allocatable :: gradu(:,:)


    ! the ones listed below were formerly in readarrays but since they are
    ! unchanged in phasta preprocessing they have been moved here to make
    ! clear that everything in readarrays can be destroyed after reading
    ! while the arrays found here live as long as the mesh is unchanged
    ! A consequence is that now readnblk must include this data module so
    ! that  these arrays can be filled as before

    real (c_double), target, allocatable :: x(:,:)
    integer (c_int), target, allocatable :: ilwork(:)
    integer (c_int), target, allocatable :: iper(:)

    integer, allocatable :: ifath(:)
    integer, allocatable :: nsons(:)

    real*8, allocatable :: velbar(:,:)

    !
    ! new arrays that are used specifically for data assimilation
    !

    integer (c_int), target, allocatable :: inodesuniq(:) ! local unique index of nodes
    integer (c_int), target, allocatable :: ilinobsfunc_sol(:,:) ! simple linear observation functions
    integer (c_int), target, allocatable :: ilinobsfunc_acc(:,:)
    integer (c_int), target, allocatable :: ilinobsfunc_disp(:,:)
    integer (c_int), target, allocatable :: obsfunc_dist(:) ! distance observation flag
 
    !linear algebra matrices

    real*8, allocatable, dimension(:,:) :: aperm,  atemp, atempS
    real*8, allocatable, dimension(:,:,:) :: apermS

    real*8, allocatable, dimension(:,:) :: lhsP, lhsK, lhsS


    !ALE global variables; MAF 06/10/2016
    real (c_double), target, allocatable :: aMesh(:,:)
    real (c_double), target, allocatable :: uMesh(:,:)
    real (c_double), target, allocatable :: dispMesh(:,:)
    real (c_double), target, allocatable :: x_iniMesh(:,:)
    real (c_double), target, allocatable :: aMeshold(:,:)
    real (c_double), target, allocatable :: uMeshold(:,:)
    real (c_double), target, allocatable :: dispMeshold(:,:)
    real (c_double), target, allocatable :: xMeshold(:,:)
    real (c_double), target, allocatable :: aMeshinc(:,:)
    ! integer, allocatable :: meshBCwallIDnodes(:)


end module


!-----------------------------------------------------------------------
!
!     Initialize:
!
!-----------------------------------------------------------------------
subroutine initGlobalArrays
    !
    use globalArrays
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    if (.not. allocated(nodetagfield)) allocate (nodetagfield(numnp))
    if (.not. allocated(y)) allocate (y(nshg,ndof))
    if (.not. allocated(ac)) allocate (ac(nshg,ndof))
    if (.not. allocated(yold)) allocate (yold(nshg,ndof))
    if (.not. allocated(acold)) allocate (acold(nshg,ndof))

    if (.not. allocated(temporary_array)) allocate(temporary_array(nshg,ndof))

    if (.not. allocated(ibc)) allocate (ibc(nshg))
    if (.not. allocated(BC)) allocate (BC(nshg,ndofBC))
    if (.not. allocated(solinc)) allocate (solinc(nshg,ndof))
    if (.not. allocated(rowp)) allocate (rowp(nshg,nnz))
    if (.not. allocated(colm)) allocate (colm(nshg+1))
    !allocate (x(nshg,nsd))
    if (.not. allocated(x)) allocate (x(numnp,nsd))
    if (.not. allocated(iper)) allocate (iper(nshg))
    if (.not. allocated(res)) allocate (res(nshg,nflow))
    !
    ! note that ilwork, ifath, nsons and velbar are all allocated in readnblk.f
    ! destructor will need to kill these.  Same is true for

    !
    ! conditional allocates
    !
    ! these are currently passed but should be moved to a
    ! module so that they can be conditionally allocated to save
    ! memory when not used
    if (.not. allocated(u)) allocate (u(nshg,nsd))
    if (.not. allocated(uold)) allocate (uold(nshg,nsd))
    if (.not. allocated(uref)) allocate (uref(nshg,nsd))
    if (.not. allocated(ubar)) allocate (ubar(nshg,nsd))
    if (.not. allocated(xdist)) allocate (xdist(nshg))
    if (.not. allocated(xdnv)) allocate (xdnv(nshg,nsd))
    if (.not. allocated(df_fem)) allocate (df_fem(nshg))
      

    if (.not. allocated(rerr)) allocate (rerr(nshg,numerr)) ! this one was pulled out because it
                                                            ! is passed down.  Soon, the error
                                                            ! fields should be put in a module and
                                                            ! "used" only where needed to save
                                                            ! this substantial memory not being
                                                            ! allocated when not needed
    if(ierrcalc.ne.0) then       ! these only live at itrdrv level and thus
                                  ! safe to conditional allocate
        if (.not. allocated(ybar)) allocate (ybar(nshg,5))
        if (.not. allocated(dummyVar)) allocate (dummyVar(nshg))
    endif
    if(ihessian.ne.0) then
        if (.not. allocated(uhess)) allocate (uhess(nshg,27))
        if (.not. allocated(gradu)) allocate (gradu(nshg,9))
    endif

    if (.not. allocated(aMesh)) allocate (aMesh(nshg,3))
    if (.not. allocated(uMesh)) allocate (uMesh(nshg,3))
    if (.not. allocated(dispMesh)) allocate (dispMesh(numnp,3))
    if (.not. allocated(x_iniMesh)) allocate (x_iniMesh(numnp,3)) ! I used numnp to be consistent with 
                                                                  ! the allocation of x, but I think that
                                                                  ! in the specific case of linear tetrahedrons
                                                                  ! numnp = nshg (to be verified) MAF 11/10/2016
    if (.not. allocated(aMeshold)) allocate (aMeshold(nshg,3)) 
    if (.not. allocated(uMeshold)) allocate (uMeshold(nshg,3))
    if (.not. allocated(dispMeshold)) allocate (dispMeshold(numnp,3))
    if (.not. allocated(xMeshold)) allocate (xMeshold(numnp,3))                                                                 
    if (.not. allocated(aMeshinc)) allocate (aMeshinc(nshg,3)) 

end subroutine initGlobalArrays

subroutine destroyGlobalArrays
    !
    use globalArrays
    use phcommonvars
    use ale
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    if (allocated(nodetagfield)) deallocate (nodetagfield)
    if (allocated(y)) deallocate (y)
    if (allocated(ac)) deallocate (ac)
    if (allocated(yold)) deallocate (yold) !
    if (allocated(acold)) deallocate (acold) !
    if (allocated(temporary_array)) deallocate(temporary_array) !
    if (allocated(ibc)) deallocate (ibc)
    if (allocated(BC)) deallocate (BC)
    if (allocated(solinc)) deallocate (solinc)
    if (allocated(rowp)) deallocate (rowp)
    if (allocated(colm)) deallocate (colm)
    if (allocated(x)) deallocate (x) !
    if (allocated(iper)) deallocate (iper)
    if (allocated(res)) deallocate (res)
    if (allocated(u)) deallocate (u)
    if (allocated(uold)) deallocate (uold)!
    if (allocated(uref)) deallocate (uref)
    if (allocated(ubar)) deallocate (ubar)
    if (allocated(xdist)) deallocate (xdist)
    if (allocated(xdnv)) deallocate (xdnv)
    if (allocated(df_fem)) deallocate (df_fem)
    if (allocated(rerr)) deallocate (rerr)
    if(ierrcalc.ne.0) then       ! these only live at itrdrv level and thus
                                  ! safe to conditional allocate
        if (allocated(ybar)) deallocate (ybar)
        if (allocated(dummyVar)) deallocate (dummyVar)
    endif
    if(ihessian.ne.0) then
        if (allocated(uhess)) deallocate (uhess)
        if (allocated(gradu)) deallocate (gradu)
    endif
    ! if(aleon.eq.1) then
    if (allocated(aMesh)) deallocate (aMesh)
    if (allocated(uMesh)) deallocate (uMesh)
    if (allocated(dispMesh)) deallocate (dispMesh)
    if (allocated(x_iniMesh)) deallocate(x_iniMesh)
    if (allocated(aMeshold)) deallocate (aMeshold)
    if (allocated(uMeshold)) deallocate (uMeshold)
    if (allocated(dispMeshold)) deallocate (dispMeshold)
    if (allocated(xMeshold)) deallocate (xMeshold)
    if (allocated(aMeshinc)) deallocate (aMeshinc)
    if (allocated(meshBCwallIDnodes)) deallocate (meshBCwallIDnodes)
    ! endif

end subroutine destroyGlobalArrays
!
    ! note that ilwork, ifath, nsons and velbar are all allocated in readnblk.f
    ! destructor will need to kill these.  Same is true for
