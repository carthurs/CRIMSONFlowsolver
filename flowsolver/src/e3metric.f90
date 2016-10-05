!-----------------------------------------------------------------------
!
!  compute the metrics of the mapping from global to local 
!  coordinates and the jacobian of the mapping (weighted by 
!  the quadrature weight
!
!-----------------------------------------------------------------------
      subroutine e3metric(  xl,      shgl,     dxidx, &
                            shg,     WdetJ)

      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      real*8     xl(npro,nenl,nsd),    shgl(npro,nsd,nshl), &
                 dxidx(npro,nsd,nsd),  shg(npro,nshl,nsd),  &
                 WdetJ(npro)

      real*8     dxdxi(npro,nsd,nsd),  tmp(npro)

      integer   iprint
      logical :: exist

!
!.... compute the deformation gradient
!
      dxdxi = zero
!
       do n = 1, nenl
          dxdxi(:,1,1) = dxdxi(:,1,1) + xl(:,n,1) * shgl(:,1,n)
          dxdxi(:,1,2) = dxdxi(:,1,2) + xl(:,n,1) * shgl(:,2,n)
          dxdxi(:,1,3) = dxdxi(:,1,3) + xl(:,n,1) * shgl(:,3,n)
          dxdxi(:,2,1) = dxdxi(:,2,1) + xl(:,n,2) * shgl(:,1,n)
          dxdxi(:,2,2) = dxdxi(:,2,2) + xl(:,n,2) * shgl(:,2,n)
          dxdxi(:,2,3) = dxdxi(:,2,3) + xl(:,n,2) * shgl(:,3,n)
          dxdxi(:,3,1) = dxdxi(:,3,1) + xl(:,n,3) * shgl(:,1,n)
          dxdxi(:,3,2) = dxdxi(:,3,2) + xl(:,n,3) * shgl(:,2,n)
          dxdxi(:,3,3) = dxdxi(:,3,3) + xl(:,n,3) * shgl(:,3,n)
       enddo
!
!.... compute the inverse of deformation gradient
!
       dxidx(:,1,1) =   dxdxi(:,2,2) * dxdxi(:,3,3)  &
                      - dxdxi(:,3,2) * dxdxi(:,2,3)
       dxidx(:,1,2) =   dxdxi(:,3,2) * dxdxi(:,1,3)  &
                      - dxdxi(:,1,2) * dxdxi(:,3,3)
       dxidx(:,1,3) =  dxdxi(:,1,2) * dxdxi(:,2,3)  &
                      - dxdxi(:,1,3) * dxdxi(:,2,2)
       tmp          = one / ( dxidx(:,1,1) * dxdxi(:,1,1)  &
                             + dxidx(:,1,2) * dxdxi(:,2,1)   &
                             + dxidx(:,1,3) * dxdxi(:,3,1) )
       dxidx(:,1,1) = dxidx(:,1,1) * tmp
       dxidx(:,1,2) = dxidx(:,1,2) * tmp
       dxidx(:,1,3) = dxidx(:,1,3) * tmp
       dxidx(:,2,1) = (dxdxi(:,2,3) * dxdxi(:,3,1)  &
                      - dxdxi(:,2,1) * dxdxi(:,3,3)) * tmp
       dxidx(:,2,2) = (dxdxi(:,1,1) * dxdxi(:,3,3)  &
                      - dxdxi(:,3,1) * dxdxi(:,1,3)) * tmp
       dxidx(:,2,3) = (dxdxi(:,2,1) * dxdxi(:,1,3)  &
                      - dxdxi(:,1,1) * dxdxi(:,2,3)) * tmp
       dxidx(:,3,1) = (dxdxi(:,2,1) * dxdxi(:,3,2)  &
                      - dxdxi(:,2,2) * dxdxi(:,3,1)) * tmp
       dxidx(:,3,2) = (dxdxi(:,3,1) * dxdxi(:,1,2)  &
                      - dxdxi(:,1,1) * dxdxi(:,3,2)) * tmp
       dxidx(:,3,3) = (dxdxi(:,1,1) * dxdxi(:,2,2)  &
                      - dxdxi(:,1,2) * dxdxi(:,2,1)) * tmp
!
       WdetJ = Qwt(lcsyst,intp) / tmp
!
!.... compute the global gradient of shape-functions
!
       do n = 1, nshl
          shg(:,n,1) = shgl(:,1,n) * dxidx(:,1,1) +  &
                       shgl(:,2,n) * dxidx(:,2,1) + &
                       shgl(:,3,n) * dxidx(:,3,1)
          shg(:,n,2) = shgl(:,1,n) * dxidx(:,1,2) +  &
                       shgl(:,2,n) * dxidx(:,2,2) + &
                       shgl(:,3,n) * dxidx(:,3,2) 
          shg(:,n,3) = shgl(:,1,n) * dxidx(:,1,3) +  &
                       shgl(:,2,n) * dxidx(:,2,3) + &
                       shgl(:,3,n) * dxidx(:,3,3) 
       enddo

#if DEBUG_ALE == 1

      ! print shgl
      !-----------   
      write(*,*) 'printing shgl inside e3metric'    
      inquire(file="shgle3metric1.dat", exist=exist)
      if (exist) then
        open(794, file="shgle3metric1.dat", status="old", position="append", action="write")
      else
        open(794, file="shgle3metric1.dat", status="new", action="write")
      end if
      do iprint = 1, npro
             write(794,'(4(e40.20))')  shgl(iprint,1,1),&
                                        shgl(iprint,1,2),&
                                        shgl(iprint,1,3),&
                                        shgl(iprint,1,4)
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      enddo 
      close(794)

      inquire(file="shgle3metric2.dat", exist=exist)
      if (exist) then
        open(794, file="shgle3metric2.dat", status="old", position="append", action="write")
      else
        open(794, file="shgle3metric2.dat", status="new", action="write")
      end if
      do iprint = 1, npro
             write(794,'(4(e40.20))')  shgl(iprint,2,1),&
                                        shgl(iprint,2,2),&
                                        shgl(iprint,2,3),&
                                        shgl(iprint,2,4)
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      enddo 
      close(794)

      inquire(file="shgle3metric3.dat", exist=exist)
      if (exist) then
        open(794, file="shgle3metric3.dat", status="old", position="append", action="write")
      else
        open(794, file="shgle3metric3.dat", status="new", action="write")
      end if
      do iprint = 1, npro
             write(794,'(4(e40.20))')  shgl(iprint,3,1),&
                                        shgl(iprint,3,2),&
                                        shgl(iprint,3,3),&
                                        shgl(iprint,3,4)
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      enddo 
      close(794)




      ! print shg
      !-----------   
      write(*,*) 'printing shg inside e3metric'    
      inquire(file="shge3metric1.dat", exist=exist)
      if (exist) then
        open(794, file="shge3metric1.dat", status="old", position="append", action="write")
      else
        open(794, file="shge3metric1.dat", status="new", action="write")
      end if
      do iprint = 1, npro
             write(794,'(4(e40.20))')  shg(iprint,1,1),&
                                        shg(iprint,2,1),&
                                        shg(iprint,3,1),&
                                        shg(iprint,4,1)
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      enddo 
      close(794)

      inquire(file="shge3metric2.dat", exist=exist)
      if (exist) then
        open(794, file="shge3metric2.dat", status="old", position="append", action="write")
      else
        open(794, file="shge3metric2.dat", status="new", action="write")
      end if
      do iprint = 1, npro
             write(794,'(4(e40.20))')  shg(iprint,1,2),&
                                        shg(iprint,2,2),&
                                        shg(iprint,3,2),&
                                        shg(iprint,4,2)
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      enddo 
      close(794)

      inquire(file="shge3metric3.dat", exist=exist)
      if (exist) then
        open(794, file="shge3metric3.dat", status="old", position="append", action="write")
      else
        open(794, file="shge3metric3.dat", status="new", action="write")
      end if
      do iprint = 1, npro
             write(794,'(4(e40.20))')  shg(iprint,1,3),&
                                        shg(iprint,2,3),&
                                        shg(iprint,3,3),&
                                        shg(iprint,4,3)
               ! write(794,'(1(e40.20))') g1yi(i,1)                                   
      enddo 
      close(794)


      
#endif




       return
       end



