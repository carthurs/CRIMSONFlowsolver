        function mpoint (name,  ndim1,  ndim2,  ndim3)
!
!----------------------------------------------------------------------
!
! This function dynamically allocates memory for the arrays.
!
! input:
!  name                 : name of the array
!  ndim1, ndim2, ndim3  : dimensions of the array
!
! output:
!  mpoint               : memory location
!
! Farzin Shakib, Summer 1985.
! Taken from dlearn (modified summer 1985)
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        character*8 name
!
!.... calculate the array size
!
        idim1  = ndim1
        idim2  = ndim2 * min(idim1, 1)
        idim3  = ndim3 * min(idim2, 1)
!
!.... store the array information
!
        mpoint = mbeg
!
!.... set the memory pointer 
!
        mbeg   = mpoint + max(1,idim1) * max(1,idim2) * max(1,idim3)
!
!.... if past the end of total memory allocated, set the error message
!
        if (mbeg .gt. mend) call error ('mpoint  ', name, mbeg-mend)
!
!.... return
!
        return
        end
