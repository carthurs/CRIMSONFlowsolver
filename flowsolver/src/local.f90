        subroutine local (global, rlocal, ientmp, n, code)
!
!----------------------------------------------------------------------
!
! This subroutine performs a vector gather/scatter operation.
!
! input:
!  global (nshg,n)             : global array
!  rlocal (npro,nshl,n)         : local array
!  ien    (npro,nshl)      : nodal connectivity
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
!.... gather the data
!

          do j = 1, n
            do i = 1, nshl
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


!
!.... transfer count
!
          gbytes = gbytes + n*nshl*npro
!
!.... return
!
          return
        endif
!
!.... ------------------------->  'assembling '  <----------------------
!
        if (code .eq. 'scatter ') then
!
!.... scatter the data (possible collisions)
!
          do j = 1, n
            do i = 1, nshl
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j)  &
                                     + rlocal(nel,i,j)
              enddo
            enddo
          enddo

!
!.... transfer and flop counts
!
          sbytes = sbytes + n*nshl*npro
          flops  = flops  + n*nshl*npro
!
!.... return
!
          return
        endif
!
!.... ------------------------->  'globalizing '  <----------------------
!
        if (code .eq. 'globaliz') then
!
!.... scatter the data (possible collisions)
!
          do j = 1, n
            do i = 1, nshl
              do nel = 1,npro
                global(ien(nel,i),j) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
!
!.... return
!
          return
        endif
!
!.... --------------------------->  error  <---------------------------
!
        call error ('local   ', code, 0)
!
!.... end
!
        end
!
        subroutine localx (global, rlocal, ien, n, code)
!
!----------------------------------------------------------------------
!
! This subroutine performs a vector gather/scatter operation for the
! nodal coordinates array.
!
! input:
!  global (numnp,n)             : global array
!  rlocal (npro,nenl,n)         : local array
!  ien    (npro,nshl)      : nodal connectivity
!  n                            : number of d.o.f.'s to be copied
!  code                         : the transfer code
!                                  .eq. 'gather  ', from global to local
!                                  .eq. 'scatter ', add  local to global 
!
!
! Zdenek Johan, Winter 1992.
!----------------------------------------------------------------------
!
         use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

        dimension global(numnp,n),           rlocal(npro,nenl,n), &
                  ien(npro,nshl)
!
        character*8 code
!
!.... ------------------------>  'localization  '  <--------------------
!
        if (code .eq. 'gather  ') then
!
!.... gather the data
!
          do j = 1, n
            do i = 1, nenl
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


!
!.... transfer count
!
          gbytes = gbytes + n*nenl*npro
!
!.... return
!
          return
        endif
!
!.... ------------------------->  'assembling '  <----------------------
!
        if (code .eq. 'scatter ') then
!
!.... scatter the data (possible collisions)
!

          do j = 1, n
            do i = 1, nenl
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j)  &
                                     + rlocal(nel,i,j)
              enddo
            enddo
          enddo


!
!.... transfer and flop counts
!
          sbytes = sbytes + n*nenl*npro
          flops  = flops  + n*nenl*npro
!
!.... return
!
          return
        endif
!
!.... --------------------------->  error  <---------------------------
!
        call error ('local   ', code, 0)
!
!.... end
!
        end
!



        subroutine locali (global, rlocal, ien, n, code)
!
!----------------------------------------------------------------------
!
! This subroutine performs a vector gather/scatter operation for an
! integer array.
!
! input:
!  global (numnp,n)             : global array
!  rlocal (npro,nenl,n)         : local array
!  ien    (npro,nshl)      : nodal connectivity
!  n                            : number of d.o.f.'s to be copied
!  code                         : the transfer code
!                                  .eq. 'gather  ', from global to local
!                                  .eq. 'scatter ', add  local to global
!
!----------------------------------------------------------------------
!
         use phcommonvars
 IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

        integer global(numnp,n),           rlocal(npro,nenl,n), &
                  ien(npro,nshl)
!
        character*8 code
!
!.... ------------------------>  'localization  '  <--------------------
!
        if (code .eq. 'gather  ') then
!
!.... gather the data
!
          do j = 1, n
            do i = 1, nenl
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


!
!.... transfer count
!
          gbytes = gbytes + n*nenl*npro
!
!.... return
!
          return
        endif
!
!.... ------------------------->  'assembling '  <----------------------
!
        if (code .eq. 'scatter ') then
!
!.... scatter the data (possible collisions)
!

          do j = 1, n
            do i = 1, nenl
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j)  &
                                     + rlocal(nel,i,j)
              enddo
            enddo
          enddo


!
!.... transfer and flop counts
!
          sbytes = sbytes + n*nenl*npro
          flops  = flops  + n*nenl*npro
!
!.... return
!
          return
        endif
!
!.... --------------------------->  error  <---------------------------
!
        call error ('local   ', code, 0)
!
!.... end
!
        end



!        subroutine localM (global, xKebe, xGoC, ien)
!c
!c----------------------------------------------------------------------
!c This routine assembles a global tangent matrix from the element
!c matrices.
!c
!c
!c 
!c
!c
!c                         |  C      G^T |
!c           globalK   =   |             |
!c                         |  G      K   |   
!c
!c
!c
!c
!c----------------------------------------------------------------------
!c
!         use phcommonvars
!
!        dimension global(nshg*4,nshg*4),xKebe(npro,3*nshl,3*nshl),
!     &            xGoC(npro,4*nshl,nshl),
!     &            ien(npro,nshape)
!c
!        character*8 code
!c
!c.... ------------------------->  'assembling '  <----------------------
!c
!
!c     
!c.... scatter the data (possible collisions)
!c
!
!c
!c.... k
!c          
!          do iel = 1, numel
!
!             do i = 1, nshl
!                i0 = (i-1)*3
!c                
!                do j = 1, nshl
!                   j0 = (j-1)*3 
!c
!                   ia = (ien(iel,i)-1)*4 + 1
!                   ib = (ien(iel,j)-1)*4 + 1 
!c                      
!                   global(ia+1,ib+1) = global(ia+1,ib+1) 
!     &                                       + xKebe(iel,i0+1,j0+1)
!                   global(ia+1,ib+2) = global(ia+1,ib+2) 
!     &                                       + xKebe(iel,i0+1,j0+2)
!                   global(ia+1,ib+3) = global(ia+1,ib+3) 
!     &                                       + xKebe(iel,i0+1,j0+3)
!                   global(ia+2,ib+1) = global(ia+2,ib+1) 
!     &                                       + xKebe(iel,i0+2,j0+1)
!                   global(ia+2,ib+2) = global(ia+2,ib+2) 
!     &                                       + xKebe(iel,i0+2,j0+2)
!                   global(ia+2,ib+3) = global(ia+2,ib+3) 
!     &                                       + xKebe(iel,i0+2,j0+3)
!                   global(ia+3,ib+1) = global(ia+3,ib+1) 
!     &                                       + xKebe(iel,i0+3,j0+1)
!                   global(ia+3,ib+2) = global(ia+3,ib+2) 
!     &                                       + xKebe(iel,i0+3,j0+2)
!                   global(ia+3,ib+3) = global(ia+3,ib+3) 
!     &                                       + xKebe(iel,i0+3,j0+3)
!c                   
!                enddo
!c
!             enddo
!c
!          enddo
!
!c
!c.... G and G^T
!c          
!          do iel = 1, numel
!
!             do i = 1, nshl
!                i0 = (i-1)*3
!                do j = 1, nshl 
!                
!                   ia = (ien(iel,i)-1)*4 + 1
!                   ib = (ien(iel,j)-1)*4 + 1 
!c                      
!                global(ia+1,ib  ) = global(ia+1,ib  )+ xGoC(iel,i0+1,j)
!                global(ia+2,ib  ) = global(ia+2,ib  )+ xGoC(iel,i0+2,j)
!                global(ia+3,ib  ) = global(ia+3,ib  )+ xGoC(iel,i0+3,j)
!                global(ia  ,ib+1) = global(ia  ,ib+1)+ xGoC(iel,i0+1,j)
!                global(ia  ,ib+2) = global(ia  ,ib+2)+ xGoC(iel,i0+2,j)
!                global(ia  ,ib+3) = global(ia  ,ib+3)+ xGoC(iel,i0+3,j)
!
!c
!             enddo
!c
!          enddo
!       enddo
!       
!c
!c.... C
!c
!          do iel = 1, numel
!             do i = 1, nshl
!                i0 = 3*nshl + i
!                do j = 1, nshl
!                   ia = (ien(iel,i)-1)*4 + 1
!                   ib = (ien(iel,j)-1)*4 + 1
!c                      
!                   global(ia,ib) = global(ia,ib) + xGoC(iel,i0,j)
!c
!                enddo
!             enddo
!             
!c
!          enddo
!       
!          
!       
!cad	  ttim(4) = ttim(4) + secs(0.0)
!
!c
!c.... transfer and flop counts
!c
!          sbytes = sbytes + nshl*nenl*npro
!          flops  = flops  + nshl*nenl*npro
!c
!c.... return
!c
!cad          call timer ('Back    ')
!          return
!c
!c.... --------------------------->  error  <---------------------------
!c
!        call error ('local   ', code, 0)
!c
!c.... end
!c
!        end
!c
!


        subroutine localSum (global, rlocal, ientmp, nHits, n)
!
!----------------------------------------------------------------------
!
!  sum the data from the local array to the global degrees of
!  freedom and keep track of the number of locals contributing
!  to each global dof. This may be used to find the average.
!
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

        dimension global(nshg,n),           rlocal(npro,nshl,n), &
                  ien(npro,nshl),           ientmp(npro,nshl), &
                  nHits(nshg)
!
!.... cubic basis has negatives in ien
!
        if (ipord > 2) then
           ien = abs(ientmp)
        else
           ien = ientmp
        endif
!
!.... ------------------------->  'assembling '  <----------------------
!
        do j = 1, n
           do i = 1, nshl
              do nel = 1,npro
                 idg = ien(nel,i)
                 global(idg,j) = global(idg,j) + rlocal(nel,i,j)
              enddo
           enddo
        enddo
        do i = 1, nshl
           do nel = 1,npro
              idg = ien(nel,i)
              nHits(idg) = nHits(idg) + 1
           enddo
        enddo
!
!.... end
!
        end
 
      subroutine localb (global, rlocal, ientmp, n, code)
!
!----------------------------------------------------------------------
!
! This subroutine performs a vector gather/scatter operation on boundary only.
!
! input:
!  global (nshg,n)             : global array
!  rlocal (npro,nshl,n)         : local array
!  ien    (npro,nshl)      : nodal connectivity
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

        dimension global(nshg,n),           rlocal(npro,nshlb,n), &
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
!ad          call timer ('Gather  ')
!
!.... gather the data
!

          do j = 1, n
            do i = 1, nshlb
              rlocal(:,i,j) = global(ien(:,i),j)
            enddo
          enddo


!
!.... transfer count
!
          gbytes = gbytes + n*nshl*npro
!
!.... return
!
!ad          call timer ('Back    ')
          return
        endif
!
!.... ------------------------->  'assembling '  <----------------------
!
        if (code .eq. 'scatter ') then
!
!.... set timer
!
!ad          call timer ('Scatter ')
!
!.... scatter the data (possible collisions)
!
          do j = 1, n
            do i = 1, nshlb
              do nel = 1,npro
                global(ien(nel,i),j) = global(ien(nel,i),j)  &
                                     + rlocal(nel,i,j)
              enddo
            enddo
          enddo

!
!.... transfer and flop counts
!
          sbytes = sbytes + n*nshlb*npro
          flops  = flops  + n*nshlb*npro
!
!.... return
!
!AD          call timer ('Back    ')
          return
        endif
!
!.... ------------------------->  'globalizing '  <----------------------
!
        if (code .eq. 'globaliz') then
!
!.... scatter the data (possible collisions)
!
          do j = 1, n
            do i = 1, nshlb
              do nel = 1,npro
                global(ien(nel,i),j) = rlocal(nel,i,j)
              enddo
            enddo
          enddo
!
!.... return
!
!ad          call timer ('Back    ')
          return
        endif
!
!.... --------------------------->  error  <---------------------------
!
        call error ('local   ', code, 0)
!
!.... end
!
        end
!




