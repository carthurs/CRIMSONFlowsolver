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

	
	real*8, allocatable :: updatedMeshCoordinates(:,:)
	real*8, allocatable :: updatedMeshVelocities(:,:)
	character(len=50), public :: folderMeshEvolutionALE



	contains

	subroutine initialize_ALE(aleType)
		implicit none
		integer :: aleType
		integer :: ierr, rerr
		integer :: fnum = 145
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

	    elseif (aleType.eq.2) then	
	    	filename = "pathtoFolderMeshEvolution.dat"
	    	inquire(file=filename, exist=file_exists)	
	    	if (file_exists) then
	    		open(fnum, file='pathtoFolderMeshEvolution.dat', status='old', iostat=ierr)         	      	
		      	if (ierr .eq. int(0)) then           	
		         	read(fnum,*,iostat=rerr) folderMeshEvolutionALE
		        endif
		        close(fnum) 
		        write(*,*) "path to folder Mesh Evolution: ",folderMeshEvolutionALE
	    	else
	    		write(*,*) 'Error, Inner Mesh motion. File',filename, ' not found. Exiting.'
				stop
	    	endif

	    	! folderMeshEvolutionALE = "../../00.geometry/imposed_mesh_motion/"
	    	filename = "coordinates_mesh00001.dat"
	    	filename = trim(folderMeshEvolutionALE)//trim(filename)
	    	inquire(file=filename, exist=file_exists)	
	    	if (file_exists) then
				
			else  
	        	write(*,*) 'Error, Inner Mesh motion. File',filename, ' not found. Exiting.'
				stop
			endif

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


	subroutine getMeshVelocities(aleType,uMesh,uMeshSize,step_number)
		implicit none
		integer :: uMeshSize
		integer :: aleType
		integer :: step_number
		real*8, intent(inout)  :: uMesh(uMeshSize,3)

        ! if uniform velocity
	    if (aleType.eq.1) then
	        uMesh(:,1) = globalRigidVelocity(1)
	        uMesh(:,2) = globalRigidVelocity(2)
	        uMesh(:,3) = globalRigidVelocity(3)
	    elseif (aleType.eq.2) then
	    	call readUpdatedMeshVariablesFromFile(uMeshSize,step_number)
	    	uMesh = updatedMeshVelocities
 		else 
            uMesh(:,1) = real(0.0,8)
            uMesh(:,2) = real(0.0,8)
            uMesh(:,3) = real(0.0,8)
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