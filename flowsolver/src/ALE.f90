!
! ALE module
! 

module ale

	implicit none


	! prescribed velocity flags
	integer :: uniformConstantVelocity       = int(0)
	integer :: uniformTimeVaryingVelocity    = int(0)
	integer :: nonUniformConstantVelocity    = int(0)
	integer :: nonUniformTimeVaryingVelocity = int(0)

	real*8, public :: globalMeshVelocity(3) 

	contains

    ! subroutine to read global mesh velocity

	subroutine readGlobalMeshVelocity()
		implicit none
		integer :: ierr, rerr      
		integer :: fnum = 145
		real*8 :: testvalue

		! open file and if no read error set global mesh velocity, else set zero
      	open(fnum, file='globalMeshVelocity.dat', status='old', iostat=ierr)         	      	
      	if (ierr .eq. int(0)) then           	
         	read(fnum,*,iostat=rerr) globalMeshVelocity(:)
         	uniformConstantVelocity = int(1)
		else  
        	globalMeshVelocity(:) = real(0.0,8)
         	uniformConstantVelocity = int(0)      	
		end if
      	close(fnum)     
	end subroutine readGlobalMeshVelocity
   
    ! subroutine to add global mesh velocity to initial solution
    
	subroutine addGlobalMeshVelocityToSolution(y,nshg,ndof)
		implicit none
		integer :: nshg, ndof
		real*8, intent(inout)  :: y(nshg,ndof)

		! add global mesh velocity to solution 
		y(:,1) = y(:,1) + globalMeshVelocity(1)
		y(:,2) = y(:,2) + globalMeshVelocity(2)
		y(:,3) = y(:,3) + globalMeshVelocity(3)		
	end subroutine addGlobalMeshVelocityToSolution

	! subroutine to return mesh velocity 

	subroutine getMeshVelocities(uMesh1,uMesh2,uMesh3,uMeshSize)
	implicit none
		integer :: uMeshSize
		real*8, intent(inout)  :: uMesh1(uMeshSize), uMesh2(uMeshSize), uMesh3(uMeshSize)

        ! if uniform velocity 
        if (uniformConstantVelocity .eq. int(1)) then
 			uMesh1(:) = globalMeshVelocity(1) 
 			uMesh2(:) = globalMeshVelocity(2) 
 			uMesh3(:) = globalMeshVelocity(3) 
 		else 
            uMesh1(:) = real(0.0,8)
            uMesh2(:) = real(0.0,8)
            uMesh3(:) = real(0.0,8)
        end if

	end subroutine getMeshVelocities



end module ale