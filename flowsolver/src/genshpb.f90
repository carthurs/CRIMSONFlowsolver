      subroutine genshpb (shpb,    shglb, nshpb, nblk)  
!
!----------------------------------------------------------------------
!
! This subroutine generates shape functions for triangular,
! quadrilateral, tetrahedron, wedge and brick elements.
!
! Christian Whiting, Winter 1999
!----------------------------------------------------------------------
!
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension shpb(MAXTOP,maxsh,MAXQPT),  &
                shglb(MAXTOP,nsd,maxsh,MAXQPT)
!
!.... loop through element blocks
!
      do iblk = 1, nblk
!
!.... get coord. system and element type 
!

         lcsyst = lcblkb(3,iblk)

!.... generate the shape-functions in local coordinates
!
         select case ( lcsyst )
         case ( 1 )             ! tets
            nshl=lcblkb(9,iblk)
            do i=1,nintb(lcsyst)  
               call shpTet(ipord,Qptb(1,1:3,i),shpb(1,:,i), &
                    shglb(1,:,:,i))
            enddo
            shglb(1,:,1:nshl,1:nintb(lcsyst)) =  &
                 shglb(1,:,1:nshl,1:nintb(lcsyst))/two
!     
         case ( 2 )             ! hexes
!     
            do i=1,nintb(lcsyst)
               call shpHex  (ipord, Qptb(2,1:3,i),shpb(2,:,i), &
                    shglb(2,:,:,i))
            enddo
!     
         case ( 3 )             ! wedges with tri bd face
            
            do i=1,nintb(lcsyst)
               call shp6w (ipord,Qptb(3,1:3,i), &
                    shpb(3,:,i),shglb(3,:,:,i))
            enddo
!     
         case ( 4 )             ! wedges with quad bd face
!     
            do i=1,nintb(lcsyst)
               call shp6w (ipord,Qptb(4,1:3,i), &
                    shpb(4,:,i),shglb(4,:,:,i))
            enddo
         case ( 5 )             ! pyramids with quad bd face
!     
            do i=1,nintb(lcsyst)
               call shppyr (ipord,Qptb(5,1:3,i), &
                    shpb(5,:,i),shglb(5,:,:,i))
            enddo
!
         case ( 6 )             ! pyramids with quad bd face
!     
            do i=1,nintb(lcsyst)
               call shppyr (ipord,Qptb(6,1:3,i), &
                    shpb(6,:,i),shglb(6,:,:,i))
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
