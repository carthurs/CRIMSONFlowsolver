module debuggingTools

	contains

	subroutine debug_write(variableName,variableValue)

		implicit none

		character(len=*), intent(in) :: variableName
		real*8, intent(in) :: variableValue
		integer :: rank

		rank = getProcessorRank_MPI()

		if (rank .eq. int(0)) then
			write(*,*) 'Debug ==> ', variableName, variableValue
		end if

	end subroutine debug_write

	subroutine debug_file(fileName,variableValue)

		implicit none

		character(len=*), intent(in) :: fileName
		real*8, intent(in) :: variableValue
		integer :: rank
		logical :: exist

		rank = getProcessorRank_MPI()

		if (rank .eq. int(0)) then
			inquire(file=fileName,exist=exist)
			if(exist) then
               open(unit=119, file=fileName, status='old',position='append',action='write')
            else
               open(119, file=fileName,status='new',action='write')
            end if
            write(119,*) variableValue
            close(119)
		end if
	end subroutine debug_file

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

end module debuggingTools