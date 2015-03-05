module errorManagement
	implicit none

	contains

	function getProcessorRank_MPI() result(rank)

	     implicit none

	     include "mpif.h"

	     integer :: err
	     integer :: rank

	     call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)

	     if (err .ne. int(0)) then
	        write(*,*) "MPI error during call to getProcessorRank_MPI. Halting."
	        stop
	     end if

	     return

	  end function getProcessorRank_MPI

	subroutine write_to_stderr(message)
		! The Fortran standard doesn't guarantee which i/o unit is connected
		! to stderr; the intrinsic module iso_fortran_env provides named constants
		! to fix this problem. Requires a Fortran 2003 compiler.
		use, intrinsic :: iso_fortran_env, only : error_unit

		implicit none

		character(len=*), intent(in) :: message
		integer :: rank

		rank = getProcessorRank_MPI()

		if (rank .eq. 0) then
			write(error_unit,*) message
		end if

	end subroutine write_to_stderr


end module errorManagement