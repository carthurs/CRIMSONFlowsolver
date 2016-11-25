!
! ALE module
! 

module ale

implicit none

private ! everything is private by default

! public routines
public :: getmeshvelocities
public :: readUpdatedMeshVariablesFromFile
public :: addGlobalRigidVelocityToInitialSolution
public :: updateMeshVariables
public :: initialize_ALE



! public variables
real*8, public :: globalRigidVelocity(3)
! integer, public :: aleRigid = int(0)
! integer, public :: aleImposedMeshMotion = int(0)

integer, allocatable, public :: meshBCwallIDnodes(:)
integer, public :: meshBCwallnnodes
! integer, allocatable, public :: noslipnodes(:)

real*8, allocatable :: updatedMeshCoordinates(:,:)
real*8, allocatable :: updatedMeshVelocities(:,:)
real*8, allocatable :: updatedMeshAcceleration(:,:)
real*8, allocatable :: innerMeshMotionParameters(:)
character(len=50), public :: folderMeshEvolutionALE



contains

subroutine initialize_ALE(aleType)
implicit none
integer :: aleType
integer :: ierr, rerr
integer :: fnum = 145
integer :: nparam, i
character(len=100) :: filename
logical :: file_exists

! if (aleType.eq.1) then
! 	aleRigid = int(1) 
! elseif (aleType.eq.2) then
! 	aleImposedMeshMotion = int(1)
! endif

globalRigidVelocity(:) = real(0.0,8)

if (aleType.eq.1) then
    filename = 'globalRigidVelocity.dat'
    inquire(file=filename, exist=file_exists)
	if (file_exists) then
        open(fnum, file='globalRigidVelocity.dat', status='old', iostat=ierr)         
      	if (ierr .eq. int(0)) then           
            read(fnum,*,iostat=rerr) globalRigidVelocity(:)
        endif
        write(*,*) "globalRigidVelocity =",globalRigidVelocity
	else  
        write(*,*) 'Error, Rigid Body Motion. File ', filename, ' not found. Exiting.'
        stop
	endif     
    close(fnum) 

elseif ((aleType.eq.2).or.(aleType.eq.3)) then
    filename = "innerMeshMotionParameters.dat"
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
        open(fnum, file='innerMeshMotionParameters.dat', status='old', iostat=ierr)
        if (ierr .eq. int(0)) then
            read(fnum,*,iostat=rerr) nparam
            if (.not.allocated(innerMeshMotionParameters)) then
                allocate(innerMeshMotionParameters(nparam))
			endif 
	        do i=1,nparam
                read(fnum,*) innerMeshMotionParameters(i)
            enddo
            write(*,*) "inner Mesh Motion Parameters: ",innerMeshMotionParameters
        endif    
        close(fnum)

    else
        write(*,*) 'Error, Inner Mesh motion. File',filename, ' not found. Exiting.'
        stop
	endif       
endif




end subroutine initialize_ALE

    
subroutine addGlobalRigidVelocityToInitialSolution(y,nshg,ndof)
   implicit none
   integer :: nshg, ndof
   real*8, intent(inout)  :: y(nshg,ndof)

! add global mesh velocity to solution 
   y(:,1) = y(:,1) + globalRigidVelocity(1)
   y(:,2) = y(:,2) + globalRigidVelocity(2)
   y(:,3) = y(:,3) + globalRigidVelocity(3)
end subroutine addGlobalRigidVelocityToInitialSolution


subroutine getMeshVelocities(aleType,uMesh,aMesh,x_ini,nnodes,step_number,dt)


implicit none
integer :: nnodes
integer :: aleType
integer :: step_number
real*8, intent(inout)  :: uMesh(nnodes,3), aMesh(nnodes,3)
real*8  :: dt, time_current
real*8  :: x_ini(nnodes,3)
real*8  :: x_ini1(nnodes), x_ini2(nnodes), x_ini3(nnodes)
real*8  :: R_ini(nnodes), phi_ini(nnodes)
real*8  :: aR, bR, aZ, bZ, lZ, lR, t0, tini
real*8  :: r_def(nnodes), x_def1(nnodes), x_def2(nnodes)
real*8  :: x_def3(nnodes), vr_def(nnodes), v3_def(nnodes)
real*8  :: ar_def(nnodes), a3_def(nnodes)
real*8  :: sinx3(nnodes), sinr(nnodes)
real*8  :: v1_def(nnodes), v2_def(nnodes)
real*8  :: a1_def(nnodes), a2_def(nnodes)
real*8  :: pi = 3.1415926535897932384626433832795d0
real*8  :: two    = 2.0000000000000000000000000000000d0
integer :: option_analytical

option_analytical = 1


! if uniform velocity
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (aleType.eq.1) then
    uMesh(:,1) = globalRigidVelocity(1)
    uMesh(:,2) = globalRigidVelocity(2)
    uMesh(:,3) = globalRigidVelocity(3)

! if imposed inner mesh motion
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif ((aleType.eq.2).or.(aleType.eq.3)) then

    ! let's start hardwiring the mapping parameters
    aR = innerMeshMotionParameters(1)
    bR = innerMeshMotionParameters(2)
    lR = innerMeshMotionParameters(3)
   
    aZ = innerMeshMotionParameters(4)
    bZ = innerMeshMotionParameters(5)
    lZ = innerMeshMotionParameters(6)

    t0 = innerMeshMotionParameters(7)
    tini = innerMeshMotionParameters(8)
    

    !be careful because this routine changes the whole domain 
    !(not only the interior nodes). Hopefully, all the boundary nodes
    !are exactly at the analytical boundary....

    time_current = dt*step_number

    ! write(*,*) "tini = ",tini
    ! write(*,*) "tcurrent = ",time_current

    !Get initial cartesian coordinates
    !---------------------------------
    x_ini1 = x_ini(:,1)
    x_ini2 = x_ini(:,2)
    x_ini3 = x_ini(:,3)

    !Calculate initial cylindrical coordinates
    !-----------------------------------------
    R_ini = sqrt(x_ini1**2 + x_ini2**2)
    phi_ini = atan2(x_ini2,x_ini1)


    sinx3 = sin(2*pi*x_ini3*(1.0d0/lZ))
    sinr  = sin(2*pi*R_ini*(1.0d0/lR))

    !Calculate current cylindrical coordinates and velocity
    !------------------------------------------------------

    if (option_analytical.eq.1) then
        !Radial direction
        r_def = R_ini  + aR*sinx3*sinr*&
                   sin(bR*pi*(time_current-tini)*(1.0d0/t0))

        vr_def =         aR*sinx3*sinr*bR*pi*(1.0d0/t0)*&
                   cos(bR*pi*(time_current-tini)*(1.0d0/t0))

        ar_def =    (-1.0d0)*aR*sinx3*sinr*((bR*pi*(1.0d0/t0))**2)*&
                   sin(bR*pi*(time_current-tini)*(1.0d0/t0))
                 

        !Axial direction
        x_def3 = x_ini3 + aZ*sinx3*sinr*&
                   sin(bZ*pi*(time_current-tini)*(1.0d0/t0))

        v3_def =          aZ*sinx3*sinr*bZ*pi*(1.0d0/t0)*&
                   cos(bZ*pi*(time_current-tini)*(1.0d0/t0))

        a3_def =     (-1.0d0)*aZ*sinx3*sinr*((bZ*pi*(1.0d0/t0))**2)*&
                   sin(bZ*pi*(time_current-tini)*(1.0d0/t0)) 
    elseif (option_analytical.eq.2) then
        !Radial direction
        r_def = R_ini  + aR*sinx3*sinr*&
                   (cos(bR*pi*(time_current-tini)*(1.0d0/t0))-1.0d0)

        vr_def =         (-1.0d0)*aR*sinx3*sinr*bR*pi*(1.0d0/t0)*&
                   sin(bR*pi*(time_current-tini)*(1.0d0/t0))

        ar_def =    (-1.0d0)*aR*sinx3*sinr*((bR*pi*(1.0d0/t0))**2)*&
                   cos(bR*pi*(time_current-tini)*(1.0d0/t0))
                 

        !Axial direction
        x_def3 = x_ini3 + aZ*sinx3*sinr*&
                   (cos(bZ*pi*(time_current-tini)*(1.0d0/t0))-1.0d0)

        v3_def =          (-1.0d0)*aZ*sinx3*sinr*bZ*pi*(1.0d0/t0)*&
                   sin(bZ*pi*(time_current-tini)*(1.0d0/t0))

        a3_def =     (-1.0d0)*aZ*sinx3*sinr*((bZ*pi*(1.0d0/t0))**2)*&
                   cos(bZ*pi*(time_current-tini)*(1.0d0/t0)) 


    endif   



    !Calculate current cartesian coordinates and velocity
    !----------------------------------------------------					 	
    x_def1 = r_def*cos(phi_ini)
    x_def2 = r_def*sin(phi_ini) 

    v1_def = vr_def*cos(phi_ini)
    v2_def = vr_def*sin(phi_ini)

    a1_def = ar_def*cos(phi_ini)
    a2_def = ar_def*sin(phi_ini)    

    !Update uMesh and coordinates
    !-----------------------------
    if (.not.allocated(updatedMeshCoordinates)) then
        allocate(updatedMeshCoordinates(nnodes,3))
	endif

    !If t<tini, nothing
    if (time_current<tini) then
        ! write(*,*) "time_current < tini"
        uMesh(:,1) = 0.0d0
        uMesh(:,2) = 0.0d0
        uMesh(:,3) = 0.0d0

        aMesh(:,1) = 0.0d0
        aMesh(:,2) = 0.0d0
        aMesh(:,3) = 0.0d0

        updatedMeshCoordinates(:,1) = x_ini1
        updatedMeshCoordinates(:,2) = x_ini2
        updatedMeshCoordinates(:,3) = x_ini3

        ! write(*,*) "aMesh = 0"
        ! write(*,*) aMesh(1800:1810,1)
        ! write(*,*) aMesh(1800:1810,2)

    !Else, apply velocity
    else 
        ! write(*,*) "time_current not < tini"
        uMesh(:,1) = v1_def
        uMesh(:,2) = v2_def
        uMesh(:,3) = v3_def

        aMesh(:,1) = a1_def
        aMesh(:,2) = a2_def
        aMesh(:,3) = a3_def

        updatedMeshCoordinates(:,1) = x_def1
        updatedMeshCoordinates(:,2) = x_def2
        updatedMeshCoordinates(:,3) = x_def3

        ! write(*,*) "aMesh = adef"
        ! write(*,*) aMesh(1800:1810,1)
        ! write(*,*) aMesh(1800:1810,2)
    endif


! if no ALE activated
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
    uMesh(:,1) = real(0.0,8)
    uMesh(:,2) = real(0.0,8)
    uMesh(:,3) = real(0.0,8)

    aMesh(:,1) = 0.0d0
    aMesh(:,2) = 0.0d0
    aMesh(:,3) = 0.0d0
endif

end subroutine getMeshVelocities


subroutine readUpdatedMeshVariablesFromFile(nnodes,step_number)
implicit none
integer, intent(in) :: step_number
integer, intent(in) :: nnodes
integer :: i, ierr, fnum = 123
character(100) :: filename


if (.not.allocated(updatedMeshCoordinates)) then
    allocate(updatedMeshCoordinates(nnodes,3))
endif 

if (.not.allocated(updatedMeshVelocities)) then
    allocate(updatedMeshVelocities(nnodes,3))
endif

updatedMeshVelocities(:,:) = real(0.0,8)
updatedMeshCoordinates(:,:) = real(0.0,8)

! Read mesh coordinates from file

write(filename,'(a,I5.5,a)') 'coordinates_mesh',step_number,'.dat'
filename = trim(folderMeshEvolutionALE)//trim(filename)
! write(*,*) "opening ",filename
open(fnum, file=filename, status='old', iostat=ierr) 
if (ierr .eq. int(0)) then
	do i=1,nnodes
        read(fnum,*) updatedMeshCoordinates(i,:)
	enddo
else  
    write(*,*) "error reading updated mesh coordinates"
end if
close(fnum)


! Read mesh velocities from file

write(filename,'(a,I5.5,a)') 'velocities_mesh',step_number,'.dat'
filename = trim(folderMeshEvolutionALE)//trim(filename)
! write(*,*) "opening ",filename
open(fnum, file=filename, status='old', iostat=ierr) 
if (ierr .eq. int(0)) then
	do i=1,nnodes
        read(fnum,*) updatedMeshVelocities(i,:)
	enddo
else  
    write(*,*) "error reading updated mesh velocities"
end if
close(fnum)

! write(*,*) "HURRAH"

end subroutine readUpdatedMeshVariablesFromFile




    ! subroutine to update mesh variables
subroutine updateMeshVariables(x, numnp)
implicit none
integer, intent(in) :: numnp
real*8, intent(inout) :: x(numnp,3)

x(:,1:3) = updatedMeshCoordinates(:,1:3)


end subroutine updateMeshVariables



end module ale