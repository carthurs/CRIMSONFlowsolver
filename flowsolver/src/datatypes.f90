! 
! ***************************************************************************
! ***************************************************************************
! *** datatype module, contains data types used in the multidomain module ***
! ***************************************************************************
! ***************************************************************************
!
      module datatypes

      use pointer_data, only: r2d

      implicit none
!
! *** time data type
!
      type timedata
         private
         real*8, public, allocatable :: v(:,:)
      end type timedata   
!
! *** real*8, 1d pointers 
!
      type pntr
         private
         real*8, pointer, public :: p => null()
      end type pntr  
!
! *** surf data type
!
      type surfdata
         private
         integer, public :: num
         integer, public, allocatable :: ids(:)
      end type surfdata

      type linkedlist
         real*8 :: value
         integer :: listEntryIndex
         integer :: lengthOfListIncludingLoopMarker
         type(linkedlist), pointer :: next
      end type linkedlist

      ! 
      type modifiablestates
         real*8, pointer :: s(:,:)
         ! real*8, pointer:: s(:,:)
      end type 

      type(modifiablestates) nrcr_states
      ! real*8, allocatable :: nrcr_states


      contains





!
! *** timedata interpolation 
!
      function getvalue(curr_time,time_data)
!
      implicit none
! 
      real*8 :: getvalue
      real*8 :: curr_time
      type(timedata) :: time_data ! timedata is assumed to have size ixj
      integer :: i                ! ith index are the time points
      integer :: length           ! jth index ranges from 1 to 2   
      real*8 :: time_n            ! the 1st index is the time and the 2nd is the value
      real*8 :: time_n1
      real*8 :: d_time
      real*8 :: n_ratio
      real*8 :: n1_ratio
!
      length = size(time_data%v,1)
      time_n1 = time_data%v(length,1)
!
      if (curr_time .gt. time_n1) then
         getvalue = time_data%v(length,2)
      else
         do i = 2, length
            time_n1 = time_data%v(i,1)
            if (time_n1 .ge. curr_time) then
               time_n = time_data%v(i-1,1)
               d_time = time_n1 - time_n
               n_ratio = (time_n1 - curr_time)/d_time
               n1_ratio = (curr_time - time_n)/d_time
               getvalue = n_ratio*time_data%v(i-1,2) &
                        + n1_ratio*time_data%v(i,2)
               exit
            end if
         end do
      end if
!
      end function getvalue
!
      end module datatypes