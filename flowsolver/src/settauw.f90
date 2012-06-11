        subroutine settauw (y,  x,  BC,  &
                            ifath,   velbar)
!
!----------------------------------------------------------------------
!
! This routine computes the time varying viscous flux for turbulence wall
!  boundary elements.
!
! Zdenek Johan, Winter 1991.  (Fortran 90)
!----------------------------------------------------------------------
!
      use pointer_data
      use turbSA
      use phcommonvars
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
      dimension y(nshg,ndofl),            x(numnp,nsd),  &
                BC(nshg,ndofBC), &
                ifath(numnp),             velbar(nfath,nflow), &
                ull(nsd),                 trx(numnp,nsd), &
                ullb(nsd),                dull(nsd), &
                evisc(numnp)
      real*8 outvec(nshg,nsd+nsd+nsd)
!
! calculate the traction magnitude for all nodes on the wall
!
      trx=zero
      evisc=zero
      rm=datmat(1,2,1)
      do nodw = 1, numnp ! loop over nodes
      if ( otwn(nodw).ne.nodw ) then ! wall node check
         if (iLES.gt.0) then
            ull(:)=velbar(ifath(otwn(nodw)),2:4) ! recall old #ing
         endif
         if (iRANS.lt.0) then
            ull(:)=y(otwn(nodw),1:3)
         endif
         ub=ull(1)*wnrm(nodw,1) & ! 
           +ull(2)*wnrm(nodw,2) & ! store u.n here for now
           +ull(3)*wnrm(nodw,3)   !
!     u_parallel_to_boundary=u-(u.n)n  :
         ull(:)=ull(:)-ub*wnrm(nodw,:) ! ull in flow
!     u_b is || u_ll ||
         ub=sqrt(ull(1)**2+ull(2)**2+ull(3)**2)
!     perp d2wall is (x2-x1).n   :
         dw=abs( &
               (x(otwn(nodw),1)-x(nodw,1))*wnrm(nodw,1) &
              +(x(otwn(nodw),2)-x(nodw,2))*wnrm(nodw,2) &
              +(x(otwn(nodw),3)-x(nodw,3))*wnrm(nodw,3) &
               )
         if(ub.eq.0) then
            ut=zero
            twoub=zero
         else
            ut=utau(ub,dw,rm,nodw,x(nodw,:))   ! find utau
            twoub=-ut*ut/ub     ! find -tau_w/ub
         endif 
!
! for LES ull has the mean-parallel velocity vector.  We want the 
! instantaneous-parallel velocity vector
!
         if (iLES.gt.0) then
            ullb=ull
            ull(:)=y(otwn(nodw),1:3)
            ubn=ull(1)*wnrm(nodw,1) & ! 
              +ull(2)*wnrm(nodw,2) &  ! store u.n here for now
              +ull(3)*wnrm(nodw,3)    !
!
!     u_parallel_to_boundary=u-(u.n)n  :
!
            ull(:)=ull(:)-ubn*wnrm(nodw,:) ! ull in flow
!
! hack a limiter into this fluctuating vector. Early transients have
! huge differences from mean values
!
            dull=ull-ullb ! the current vector difference
            dullm=sqrt(dull(1)*dull(1)+dull(2)*dull(2)+dull(3)*dull(3)) &
                 + 1.0e-9
            ullbm=ub ! mag of ullb still there
!
! limit the magnitude of the difference to a 40% change from the mean.
! if less than that already we will take the whole difference, otherwise
! only take a 40% change.
!

            dullmod=min(one,0.4*ullbm/dullm)
            ull=ullb+dullmod*dull
         endif

         trx(nodw,:)=twoub*ull(:)
         if(itwmod.eq.-2) then ! effective-viscosity
            tauw=ut*ut
            BC(nodw,7)=tauw*dw/ub-rm
         endif
         if(itwmod.eq.2) then ! effective-viscosity
!
!  mag of u instantaneous
!
            ullm=sqrt(ull(1)*ull(1)+ull(2)*ull(2)+ull(3)*ull(3))
            tauw=ut*ut*ullm/ub
            evisc(nodw)=tauw*dw/ub-rm
         endif
         if((itwmod.eq.-1)) then ! slip-velocity RANS
            up=sqrt( &
                 y(nodw,1)**2  & ! 
                 +y(nodw,2)**2 & ! flow is auto-|| at boundary
                 +y(nodw,3)**2 & !
                 )/ut
            BC(nodw,7)=savarw(up,rm,saCv1,nodw,x(nodw,:))
         endif
      endif ! wallnode check
      enddo ! loop over nodes
!
! Write the traction vectors to a file restart.4077.n
!
!$$$      ilstep=4077
!$$$      outvec(:,1:3)=trx(:,1:3)
!$$$      outvec(:,4:6)=0
!$$$      do i=1,numnp
!$$$         if(otwn(i).ne. i ) outvec(i,4:6)=y(otwn(i),1:3)
!$$$      enddo
!$$$      outvec(:,7:9)=wnrm(:,1:3)
!$$$      call write_restart(myrank,ilstep,numnp,nsd*3,outvec,outvec)
!$$$      write(*,*) 'Traction dumped to restart.4077.*'
!
! Put traction calculations into BCB
!
      do iblk = 1, nelblb
         iel    = lcblkb(1,iblk)
         nenl   = lcblkb(5,iblk) ! no. of vertices per element
         nenbl  = lcblkb(6,iblk) ! no. of vertices per bdry. face
         ndofl  = lcblkb(8,iblk)
         npro   = lcblkb(1,iblk+1) - iel 
! For all elements in this block that lie on a wall, assign the traction
         do i=1,npro
            if(btest(miBCB(iblk)%p(i,1),4)) then ! wall elt
               do j = 1, nenbl
                  if(itwmod.eq.-2) then ! effective-viscosity
                     mBCB(iblk)%p(i,j,3:5)=0.0
                  endif
                  if((itwmod.eq.-1).or.(itwmod.eq.1)) then ! slip-velocity
                     mBCB(iblk)%p(i,j,3:5)=trx(mienb(iblk)%p(i,j),:)
                  endif
               enddo
            endif
         enddo
      enddo


      if(itwmod.eq.2) then   !effective viscosity
!
! For the elements which touch a modeled wall,
! modify the eddy viscosity at all quadrature points to be the average
! of the "optimal" nodal values.  This is an element-wise constant.
!
         do iblk = 1,nelblk
            nenl   = lcblk(5,iblk) ! no. of vertices per element
            iel    = lcblk(1,iblk)
            npro   = lcblk(1,iblk+1) - iel
            lcsyst = lcblk(3,iblk)
            ngauss = nint(lcsyst)
            do i=1,npro
               xnave=zero
               avevisc=zero
               do j=1,nenl
                  ev = evisc(mien(iblk)%p(i,j))
                  avevisc=avevisc + ev
                  if(ev.ne.zero) xnave=xnave+1
               enddo
               if(xnave.ne.0) mxmudmi(iblk)%p(i,1:ngauss)=avevisc/xnave
            enddo
         enddo
      endif
!
!.... end
!
        return
        end

        function utaul(u,y,rm)
        real*8 utaul
        utaul=.2
        return
        end

        function utau(u,y,rm,nodw,x)
        implicit none
        real*8 u,err,utau,yrmi,y,rm,yp,up,kup,f,dfds,rat,efac
        real*8 kup2,kup3,pt5,sxth,kappa
        integer iter
        integer nodw
        real*8 x(3)
        pt5=0.5
        sxth=0.1666666666666667
        err=1.0d-6
        utau=0.04
        yrmi=y/rm
        kappa=0.4
!$$$        B=5.5
        efac=0.1108 ! exp(-kappa*B)
        do iter=1,500
           yp=yrmi*utau
           up=u/utau
           kup=kappa*up
           kup2=kup*kup
           kup3=kup*kup2
           f=   up-yp+efac*(exp(kup)-1.0-kup - kup2*pt5 - kup3*sxth)
           dfds=up+yp+efac*(exp(kup)*kup-kup - kup2 - kup3*pt5) !/-utau in rat
           rat=f*utau/dfds ! this is -f/dfds but we did not do 1/(-utau) yet.
           utau=utau+rat
           if(abs(rat).le.err) goto 20
        enddo
        write(*,*)'utau failed to converge at ',nodw,x(1),x(2),x(3)
        write(*,*) 'u,          dwallperp,      mu'
        write(*,*) u, y, rm
        write(*,*) 'dfds,         rat,         utau'
        write(*,*) dfds,rat,utau
!        stop
!        utau=0.0
!
!  if the above fails then try a simple no-slip traction
!
        utau=sqrt(u/y*rm)

 20     continue
        return
        end  

!$$$        function utau_log_layer_only(u,y,rm)
!$$$        real*8 u,err,utau,yrmi,y,rm,lnyp,f,dfds,rat
!$$$        err=1.0d-6
!$$$        utau=0.04
!$$$        yrmi=y/rm
!$$$        do iter=1,50
!$$$           lnyp=log(yrmi*utau)
!$$$           f=u-utau*(2.5*lnyp+5.2)
!$$$           dfds=-2.5*(lnyp+6.2)
!$$$           rat=-f/dfds
!$$$           utau=utau+rat
!$$$           if(utau.gt.0.5) then !flow is laminar for now give the parabolic
!$$$              utau=sqrt(2.0*rm*u)/y
!$$$              goto 20
!$$$           endif
!$$$              
!$$$           if(abs(rat).le.err) goto 20
!$$$        enddo
!$$$        write(*,*)'utau failed to converge',r,dfds,rat,utau,y,u,rm
!$$$        stop
!$$$ 20     continue
!$$$        return
!$$$        end  

        function savarw(up,rm,cv1,nodw,x)
        implicit none
        real*8 err,savarw,rm,up,cv1,f,dfds,rat,efac
        real*8 pt5,kappa,B,xmut,chi3,denom,cv1_3
        integer iter
        integer nodw
        real*8 x(3)
        pt5=0.5
        err=1.0d-6
        savarw=rm*cv1*1.2599 ! inflection point chi=cv1*cuberoot(2)
        kappa=0.4
!$$$        B=5.5
        efac=0.1108 ! exp(-kappa*B)
        xmut=rm*kappa*efac*(exp(kappa*up)-1.0-kappa*up &
             -pt5*(kappa**2)*(up**2))
        do iter=1,50
           chi3=savarw/rm
           chi3=chi3*chi3*chi3
           cv1_3=cv1**3
           denom=chi3+cv1_3

           f=savarw*chi3/denom - xmut
           dfds=chi3*(chi3+4.0*cv1_3)/(denom**2)
           rat=-f/dfds
           savarw=savarw+rat
           if(abs(rat).le.err) goto 20
        enddo
        write(*,*)'savarw failed to converge at ',nodw,x(1),x(2),x(3)
        write(*,*) 'uplus,           mu'
        write(*,*) up,rm
        write(*,*) 'dfds,        rat,        savarw'
        write(*,*) dfds,rat,savarw
        savarw=10000.0*rm
 20     continue
        return
        end  
