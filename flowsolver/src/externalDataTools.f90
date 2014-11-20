module externalDataTools

	use datatypes

	contains

	function loadFileAsLinkedList(filename,restartFlag) result(outputLinkedList)

		implicit none

		include "mpif.h"

		character(len=*), intent(in) :: filename
		integer, intent(in) :: restartFlag

		integer :: numberOfLines
		integer :: iostatus
		integer :: ii
		integer :: ierr
		integer :: rank
		logical :: file_exists
		type(linkedlist), pointer :: outputLinkedList
		type(linkedlist), pointer :: current
		integer :: linkedListNumstart
		integer :: globalNumstart
		integer :: lastGlobalNumstartPriorToLinkedlistRestartBeingWritten
		character(len=13) :: numstartWord
		character(len=100) :: fullname
		logical :: exist

		inquire(file=filename, exist=file_exists)
		if (file_exists) then
			numberOfLines = -1 ! -1 here because the code below overcounts lines by 1.
			iostatus = 1

			open(unit=74,file=filename,status='old')
			do while (iostatus .ge. int(0))
				read(74,*,IOSTAT=iostatus)
				numberOfLines = numberOfLines + 1
			end do

			call fseek(74,0,0)

			allocate(outputLinkedList)
			outputLinkedList%next => outputLinkedList
			current => outputLinkedList

			do ii=1, numberOfLines
				read(74,*) current%value
				current%lengthOfListIncludingLoopMarker = numberOfLines+1
				current%listEntryIndex = ii
				allocate(current%next)
				current%next%next => outputLinkedList
				current => current%next
			end do

			! a nonsense value to warn that the linkedlist is about to loop back to the beginning (loop marker):
			current%value = -1.0
			current%lengthOfListIncludingLoopMarker = numberOfLines + 1
			current%listEntryIndex = 7

			close(74)

			call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
			if (ierr .ne. int(0) ) then
				write(*,*) 'MPI Error in loadFileAsLinkedList. Code:', ierr
				stop
			end if
			if (rank .eq. int(0)) then
				write(*,*) 'Loaded ', numberOfLines, ' datapoints from ', filename, '.'
			end if

			numstartWord = '-numstart.dat'

			! On restarts, cycle the linkedlist to the appropriate value,
			! otherwise, write the -numstart for this file, so that if we restart
			! before a full heart period has completed, it doesnt crash:
			fullname = filename//numstartWord
			if (restartFlag .gt. int(0)) then
				inquire(file=fullname,exist=exist)
				if (exist) then
					open(unit=112,file=fullname,status='old')
					read(112,*) linkedListNumstart
					read(112,*) ! Read & ignore comment line: '# If value matches numstart.dat, must'
					read(112,*) ! Read & ignore comment line: '# subtract 1 from linkedListNumstart before use.'
					read(112,*) lastGlobalNumstartPriorToLinkedlistRestartBeingWritten
					close(112)

					open(unit=112,file='numstart.dat',status='old')
					read(112,*) globalNumstart
					close(112)
					if ((globalNumstart .eq. lastGlobalNumstartPriorToLinkedlistRestartBeingWritten) .and. &
					     (globalNumstart .gt. int(0))) then
						! In this case we must subtract 1, because at the time-step where the restart
						! is taking place, the linked list had not yet been updated to its next value.
						linkedListNumstart = linkedListNumstart - 1
					end if
				else
					write(*,*) 'Error: restart index not found: ', fullname
					stop
				end if

				do ii=1, linkedListNumstart
					outputLinkedList => outputLinkedList%next
					if (rank .eq. int(0)) then
						write(*,*) 'Cycled linked list by 1: ', filename
					end if
				end do
			else
				if (rank .eq. int(0)) then

					open(unit=112,file=fullname,status='replace')
					write(112,*) int(0)
					write(112,*) '# If value matches numstart.dat, must'
					write(112,*) '# subtract 1 from linkedListNumstart before use.'
					write(112,*) int(0)
					close(112)

				end if
			end if

		else
			write(*,*) 'Error: File ', filename, ' not found. Exiting.'
			stop
		end if

	end function loadFileAsLinkedList


	subroutine writeRestartForLinkedList(filename)
		implicit none

		include "mpif.h"

		character(len=*), intent(in) :: filename
		integer :: linkedListNumstart
		character(len=13) :: numstartWord
		character(len=100) :: fullname
		integer :: lastGlobalNumstart
		integer :: rank
		integer :: err
		integer :: exist

		call MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
		if (err .ne. int(0) ) then
			write(*,*) 'MPI Error in writeRestartForLinkedList. Code:', err
			stop
		end if

		numstartWord = '-numstart.dat'
		fullname = filename//numstartWord

		if (rank .eq. int(0)) then

			inquire(file=fullname, exist=exist)
			if (exist) then
				open(unit=112,file=fullname,status='old')
				read(112,*) linkedListNumstart
				close(112)
			else
				write(*,*) 'Error: file not found - ', fullname, 'Halting.'
				stop
			end if

			! to ensure proper alignment of the linked list data and the
			! rest of the simulation on restart, we need to make a note of the 
			! lastGlobalNumstart (from numstart.dat), so we need to know what
			! its value is - we acquire that value now.
			open(unit=113,file='numstart.dat',status='old')
			read(113,*) lastGlobalNumstart
			close(113)

			open(unit=112,file=fullname,status='replace')
			linkedListNumstart = linkedListNumstart + 1
			write(112,*) linkedListNumstart
			write(112,*) '# If value matches numstart.dat, must'
			write(112,*) '# subtract 1 from linkedListNumstart before use.'
			write(112,*) lastGlobalNumstart
			close(112)

		end if

		call MPI_BARRIER(MPI_COMM_WORLD,err) !\todo needed?

	end subroutine writeRestartForLinkedList


end module externalDataTools