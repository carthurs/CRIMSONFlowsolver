!
! ALE module
! 

module ale

	implicit none


	real*8, public :: globalMeshVelocity(3) 

	contains

    ! subroutine to read global mesh velocity

	subroutine readGlobalMeshVelocity()
		implicit none
		integer :: ierr, rerr      
		integer :: fnum = 145
		real*8 :: testvalue

		! open file
      	open(fnum, file='globalMeshVelocity.dat', status='old', iostat=ierr)   
      	! if no read error set global mesh velocity, else set zero
      	if (ierr .eq. int(0)) then           	
         	read(fnum,*,iostat=rerr) globalMeshVelocity(:)
         	write(*,*) 'Global mesh velocity set to: ', globalMeshVelocity(:)
		else  
        	globalMeshVelocity(:) = real(0.0,8)
         	write(*,*) 'Global mesh velocity set to: ', globalMeshVelocity(:)        	
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



end module ale