        subroutine genshp (shp,    shgl, nshp, nblk)  
!
!----------------------------------------------------------------------
!
! This subroutine generates shape functions for triangular,
! quadrilateral, tetrahedron, wedge and brick elements and pyramids.
!
! Christian Whiting, Winter 1999
!----------------------------------------------------------------------
!
        use phcommonvars
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension shp(MAXTOP,maxsh,MAXQPT),  &
                  shgl(MAXTOP,nsd,maxsh,MAXQPT)
!
!.... loop through element blocks
!
        maxnint=1
          do iblk = 1, nblk
!
!.... get coord. system and element type 
!
            lcsyst = lcblk(3,iblk)
            nshl   = lcblk(10,iblk)
!
!.... generate the shape-functions in local coordinates
!
            select case ( lcsyst )
            case ( 1 )          ! tets
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)  
               call shpTet(ipord,Qpt(1,1:3,i),shp(1,:,i),shgl(1,:,:,i))
            enddo
            shgl(1,:,1:nshl,1:nint(lcsyst)) =  &
            shgl(1,:,1:nshl,1:nint(lcsyst))/two
!     
            case ( 2 )          ! hexes
!     
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shphex  (ipord, Qpt(2,1:3,i),shp(2,:,i), &
                             shgl(2,:,:,i))
            enddo
!
            case ( 3 )          ! wedges
!
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shp6w  (ipord,Qpt(3,1:3,i),shp(3,:,i), &
                             shgl(3,:,:,i))
            enddo

         case ( 5)              ! pyramids
            
               maxnint=max(maxnint,nint(lcsyst))
            do i=1,nint(lcsyst)
               call shppyr (ipord,Qpt(5,1:3,i),shp(5,:,i),shgl(5,:,:,i))
               
            enddo
!
!.... nonexistent element
!
            case default
!
            call error ('genshp  ', 'elem Cat', lelCat)
!
            end select
!
!.... end of generation
!
          enddo
!
!.... return
!
        return
        end
