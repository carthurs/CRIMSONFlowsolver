module debuggingTools

	integer :: writeCount = 0

	contains

	subroutine debug_write(variableName,variableValue)

		implicit none

		character(len=*), intent(in) :: variableName
		real*8, intent(in) :: variableValue
		integer :: rank

		rank = getProcessorRank_MPI()

		writeCount = writeCount + 1

		if (rank .eq. int(0)) then
			write(*,*) 'Debug ==> ', variableName, variableValue, ' write-index: ', writeCount
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


	subroutine listOpenFileUnits(rangeToCheck_lowerBound,rangeToCheck_upperBound)

		implicit none

		integer, intent(in) :: rangeToCheck_lowerBound
		integer, intent(in) :: rangeToCheck_upperBound
		integer :: ii
		integer :: iostat
		logical :: opened

		do ii = rangeToCheck_lowerBound, rangeToCheck_upperBound
			inquire(unit=ii, opened=opened, iostat=iostat)
        	if ((iostat.ne.0).or.(opened)) then
        		write(*,*) 'File unit still open here: ', ii
        		write(*,*) 'Note that 0, 5 and 6 are reserved, and 100-102 are disallowed.'
        	endif
        enddo

        return

	end subroutine listOpenFileUnits

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