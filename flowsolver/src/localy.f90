        subroutine localy (global, rlocal, ientmp, n, code)
!
!----------------------------------------------------------------------
!
! This subroutine performs a vector gather/scatter operation.
!
! input:
!  global (nshg,n)              : global array
!  rlocal (npro,nshl,n)         : local array
!  ien    (npro,nshl)           : nodal connectivity
!  n                            : number of d.o.f.'s to be copied
!  code                         : the transfer code
!                                  .eq. 'gather  ', from global to local
!                                  .eq. 'scatter ', add  local to global 
!                                  .eq. 'globaliz', from local to global
!
!
! Zdenek Johan, Winter 1992.
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

        dimension global(nshg,n),           rlocal(npro,nshl,n), &
                  ien(npro,nshl),           ientmp(npro,nshl)
!
        character*8 code
        
!
!.... cubic basis has negatives in ien
!
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
!
!.... ------------------------>  'localization  '  <--------------------
!
        if (code .eq. 'gather  ') then
!
!.... set timer
!
!          call timer ('Gather  ')
!
!.... gather the data to the current block 
!

!AD      rlocal = yl={P, u, v, w, T, scalar1, ...}
!AD	 global = y = {u, v, w, P, T, scalar1, ...}

!AD      Put u,v,w in the slots 2,3,4 of yl 

          do j = 1, 3
            do i = 1, nshl
              rlocal(:,i,j+1) = global(ien(:,i),j)
            enddo
          enddo

!AD      Put Pressure in the first slot of yl

          do i = 1, nshl
             rlocal(:,i,1) = global(ien(:,i),4)
          enddo

!AD      Fill in the remaining slots with T, and additional scalars
          
          if(n.gt.4) then
             do j = 5, n
                do i = 1, nshl
                   rlocal(:,i,j) = global(ien(:,i),j)
                enddo
             enddo
          endif
!
!.... transfer count
!
          gbytes = gbytes + n*nshl*npro
!
!.... return
!
!          call timer ('Back    ')
          return
        endif
!
!.... ------------------------->  'assembling '  <----------------------
!
        if (code .eq. 'scatter ') then
           write(*,*) 'do not use localy here'
        endif
!
!.... ------------------------->  'globalizing '  <----------------------
!
        if (code .eq. 'globaliz') then
           write(*,*) 'do not use localy here'
        endif
!
!.... --------------------------->  error  <---------------------------
!
        call error ('local   ', code, 0)
!
!.... end
!
        end
