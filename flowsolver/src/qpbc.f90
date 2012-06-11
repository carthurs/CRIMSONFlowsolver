      subroutine qpbc( rmass, qres, iBC, iper, ilwork )
!---------------------------------------------------------------------
!
! This routine satisfies the periodic boundary conditions
! on the diffusive flux residual and mass matrix
!
! input:
!     rmass   (nshg)              : mass matrix
!     qres    (nshg,(nflow-1)*nsd) : diffusive flux vector
! 
! output: modified qres and rmass 
!---------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      dimension rmass(nshg), qres(nshg,idflx), &
                iBC(nshg), iper(nshg),uv(nshg,2), &
                tmpvec(nshg,4), tmp(nshg)  
!
      if(iabc==1) then   !are there any axisym bc's
      do i=1,idflx/nsd
         do j=1,2
            istrt=j+(i-1)*(nflow-1)
            uv(:,j)=qres(:,istrt)
         enddo
         call rotabc(uv, iBC, 'in ')
         do j=1,2
            istrt=j+(i-1)*(nflow-1)
            qres(:,istrt)=uv(:,j)
         enddo
      enddo
      endif
!
!
!.... compute qi for node A, i.e., qres <-- qres/rmass
!
       if (numpe > 1) then
          call commu (qres  , ilwork, idflx  , 'in ')
          call commu (rmass , ilwork,  1            , 'in ')
       endif
!
!  take care of periodic boundary conditions
!  but not on surface tension terms in qres(:,10-12)
!  that are used to compute normal vector
!
        idflow = (nflow-1)*nsd
        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
!            qres(i,:) = qres(i,:) + qres(j,:)
            qres(i,1:idflow) = qres(i,1:idflow) + qres(j,1:idflow)
          endif
        enddo

        do j= 1,nshg
          if ((btest(iBC(j),10))) then
            i = iper(j)
            rmass(j) = rmass(i)
!            qres(j,:) = qres(i,:)
            qres(j,1:idflow) = qres(i,1:idflow)
          endif
        enddo
!
!.... invert the diagonal mass matrix and find q
!
        rmass = one/rmass
       
       do i=1, idflx
          qres(:,i) = rmass*qres(:,i)
       enddo
       if (isurf .eq. 1) then
          idflow=(nflow-1)*nsd
!
!.... calculation of the unit normal vector
!
           tmp =  sqrt(qres(:,idflow+1)**2 &
                     + qres(:,idflow+2)**2 &
                     + qres(:,idflow+3)**2)
           do i = 1, nshg
             if (tmp(i) .lt. 0.0001) tmp(i) = 0.0001
           end do
           tmp = one/tmp

          do i=1, nsd
             qres(:,idflow+i) = tmp*qres(:,idflow+i)
          enddo 
       endif

       if(numpe > 1) then
          call commu (qres, ilwork, idflx, 'out')    
       endif

       if(iabc==1) then         !are there any axisym bc's
!
!       slave has masters value, for abc we need to rotate it
!
          do i=1,idflx/nsd
             do j=1,2
                istrt=j+(i-1)*(nflow-1)
                uv(:,j)=qres(:,istrt)
             enddo
             call rotabc(uv, iBC, 'out')
             do j=1,2
                istrt=j+(i-1)*(nflow-1)
                qres(:,istrt)=uv(:,j)
             enddo
          enddo
       endif
       
!
!.... return
!    
        return
        end


      subroutine qpbcSclr( rmass, qres, iBC, iper, ilwork )
!---------------------------------------------------------------------
!
! This routine satisfies the periodic boundary conditions
! on the diffusive flux residual and mass matrix
!
! input:
!     rmass   (nshg)              : mass matrix
!     qres    (nshg, nsd)         : diffusive flux vector
! 
! output: modified qres and rmass 
!---------------------------------------------------------------------
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      
      dimension rmass(nshg), qres(nshg,nsd), &
                iBC(nshg), iper(nshg)

      if(iabc==1) & !are there any axisym bc's
               call rotabc(qres, iBC,  'in ')

!
!.... compute qi for node A, i.e., qres <-- qres/rmass
!
       if (numpe > 1) then
          call commu (qres  , ilwork, nsd  , 'in ')
          call commu (rmass , ilwork,  1   , 'in ')
       endif

!
!  take care of periodic boundary conditions
!
        do j= 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(i) = rmass(i) + rmass(j)
            qres(i,:) = qres(i,:) + qres(j,:)
          endif
        enddo

        do j= 1,nshg
          if (btest(iBC(j),10)) then
            i = iper(j)
            rmass(j) = rmass(i)
            qres(j,:) = qres(i,:)
          endif
        enddo
!
!.... invert the diagonal mass matrix and find q
!
        rmass = one/rmass
       
       do i=1, nsd
          qres(:,i) = rmass*qres(:,i)
       enddo

       if(numpe > 1) then
          call commu (qres, ilwork, nsd, 'out')    
       endif

      if(iabc==1) & !are there any axisym bc's
               call rotabc(qres, iBC, 'out')

!
!.... return
!    
        return
        end

