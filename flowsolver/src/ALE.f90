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
	integer, public :: aleRigid = int(0)
	integer, public :: aleImposedMeshMotion = int(0)

	
	real*8, allocatable :: updatedMeshCoordinates(:,:)
	real*8, allocatable :: updatedMeshVelocities(:,:)



	contains

	subroutine initialize_ALE(aleType)
		implicit none
		integer :: aleType
		integer :: ierr, rerr
		integer :: fnum = 145
		character(len=50) :: filename
		logical :: file_exists


		if (aleType.eq.1) then
			aleRigid = int(1) 
		elseif (aleType.eq.2) then
			aleImposedMeshMotion = int(1)
		endif

		globalRigidVelocity(:) = real(0.0,8)

		if (aleRigid.eq.1) then
			filename = 'globalRigidVelocity.dat'
			inquire(file=filename, exist=file_exists)
			if (file_exists) then
				open(fnum, file='globalRigidVelocity.dat', status='old', iostat=ierr)         	      	
		      	if (ierr .eq. int(0)) then           	
		         	read(fnum,*,iostat=rerr) globalRigidVelocity(:)
		        endif
		        write(*,*) "globalRigidVelocity =",globalRigidVelocity
			else  
	        	write(*,*) 'Error, Rigid Body Motion and file ', filename, ' not found. Exiting.'
				stop
			endif     	
	      	close(fnum) 
	    endif




	end subroutine initialize_ALE

    ! subroutine to read global mesh velocity

	! subroutine readGlobalMeshVelocity()
	! 	implicit none
	! 	integer :: ierr, rerr      
	! 	integer :: fnum = 145, fnum2 = 146
	! 	real*8 :: testvalue

	! 	! open file and if no read error set global mesh velocity, else set zero
 !      	open(fnum, file='globalMeshVelocity.dat', status='old', iostat=ierr)         	      	
 !      	if (ierr .eq. int(0)) then           	
 !         	read(fnum,*,iostat=rerr) globalMeshVelocity(:)
 !         	uniformConstantVelocity = int(1)
	! 	else  
 !        	globalMeshVelocity(:) = real(0.0,8)
 !         	uniformConstantVelocity = int(0)      	
	! 	end if
 !      	close(fnum)     

 !  !     	open(fnum2, file='globalRigidVelocity.dat', status='old', iostat=ierr)         	      	
 !  !     	if (ierr .eq. int(0)) then           	
 !  !        	read(fnum2,*,iostat=rerr) globalRigidVelocity(:)
	! 	! else  
 !  !       	globalRigidVelocity(:) = real(0.0,8)
	! 	! end if
 !  !     	close(fnum2)  

	! end subroutine readGlobalMeshVelocity
   
    ! subroutine to add global mesh velocity to initial solution
    
	subroutine addGlobalRigidVelocityToInitialSolution(y,nshg,ndof)
		implicit none
		integer :: nshg, ndof
		real*8, intent(inout)  :: y(nshg,ndof)

		! add global mesh velocity to solution 
		y(:,1) = y(:,1) + globalRigidVelocity(1)
		y(:,2) = y(:,2) + globalRigidVelocity(2)
		y(:,3) = y(:,3) + globalRigidVelocity(3)		
	end subroutine addGlobalRigidVelocityToInitialSolution

	! subroutine to return mesh velocity 

	! subroutine getMeshVelocities(uMesh1,uMesh2,uMesh3,uMeshSize)
	! 	implicit none
	! 	integer :: uMeshSize
	! 	real*8, intent(inout)  :: uMesh1(uMeshSize), uMesh2(uMeshSize), uMesh3(uMeshSize)

 !        ! if uniform velocity
 !        if (aleRigid.eq.1) then
	!         uMesh1(:) = globalRigidVelocity(1)
	!         uMesh2(:) = globalRigidVelocity(2)
	!         uMesh3(:) = globalRigidVelocity(3)
 !    !   	endif 
 !    !     if (uniformConstantVelocity .eq. int(1)) then
 ! 			! uMesh1(:) = globalMeshVelocity(1) 
 ! 			! uMesh2(:) = globalMeshVelocity(2) 
 ! 			! uMesh3(:) = globalMeshVelocity(3) 
 ! 		else 
 !            uMesh1(:) = real(0.0,8)
 !            uMesh2(:) = real(0.0,8)
 !            uMesh3(:) = real(0.0,8)
 !        end if

	! end subroutine getMeshVelocities



	subroutine getMeshVelocities(uMesh,uMeshSize)
		implicit none
		integer :: uMeshSize
		real*8, intent(inout)  :: uMesh(uMeshSize,3)

        ! if uniform velocity
        if (aleRigid.eq.1) then
	        uMesh(:,1) = globalRigidVelocity(1)
	        uMesh(:,2) = globalRigidVelocity(2)
	        uMesh(:,3) = globalRigidVelocity(3)
    !   	endif 
    !     if (uniformConstantVelocity .eq. int(1)) then
 			! uMesh1(:) = globalMeshVelocity(1) 
 			! uMesh2(:) = globalMeshVelocity(2) 
 			! uMesh3(:) = globalMeshVelocity(3) 
 		else 
            uMesh(:,1) = real(0.0,8)
            uMesh(:,2) = real(0.0,8)
            uMesh(:,3) = real(0.0,8)
        end if

	end subroutine getMeshVelocities


	subroutine readUpdatedMeshVariablesFromFile(step_number,nnodes)
		implicit none
		integer, intent(in) :: step_number
		integer, intent(in) :: nnodes
		integer :: i, ierr, fnum = 123
		character(50) :: filename


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
		open(fnum, file=filename, status='old', iostat=ierr) 
		if (ierr .eq. int(0)) then
			do i=1,nnodes
				read(fnum,*) updatedMeshVelocities(i,:)
			enddo
		else  
        	write(*,*) "error reading updated mesh velocities"
		end if
		close(fnum)

		write(*,*) "HURRAH"

	end subroutine readUpdatedMeshVariablesFromFile




    ! subroutine to update mesh variables
	subroutine updateMeshVariables(x, y, numnp)
		implicit none
		integer, intent(in) :: numnp
		real*8, intent(inout) :: x(numnp,3), y(numnp,5)

		x(:,1:3) = updatedMeshCoordinates(:,1:3)


	end subroutine updateMeshVariables



end module ale