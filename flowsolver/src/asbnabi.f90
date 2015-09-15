      subroutine AsBNABI ( x,       shpb, &
                          ienb,  iBCB)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  boundary elements.
!
!----------------------------------------------------------------------
!
      use pvsQbi
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        real*8, allocatable :: locationsInbnormToZeroOut(:)
!
        dimension xlb(npro,nenl,nsd),    bnorm(npro,nsd), &
                  rl(npro,nshl,nsd),     WdetJb(npro)

        dimension x(numnp,nsd), &
                  shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb), &
                  ienb(npro,nshl), &
                  iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl), &
                  shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

!
        dimension dxdxib(npro,nsd,nsd),      temp(npro), &
                  temp1(npro),               temp2(npro), &
                  temp3(npro), &
                  v1(npro,nsd),              v2(npro,nsd)

!
!.... get the matrix of mode signs for the hierarchic basis functions
!
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
!
!.... gather the variables
!
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
!
!.... get the boundary element residuals
!
        rl  = zero
!
!.... compute the nodes which lie on the boundary (hierarchic)
!
        call getbnodes(lnode)
!
!.... loop through the integration points
!
        ngaussb = nintb(lcsyst)

        
        do intp = 1, ngaussb
!
!.... get the hierarchic shape functions at this int point
!
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, &
                    shpfun,       shdrv)

!
!.... compute the normal to the boundary. This is achieved by taking
!     the cross product of two vectors in the plane of the 2-d 
!     boundary face.
!
           if(lcsyst.ne.6) then
              ipt3=3
           else
              ipt3=5
           endif
           v1 = xlb(:,2,:) - xlb(:,1,:)
           v2 = xlb(:,ipt3,:) - xlb(:,1,:)

!
!.....The following are done in order to correct temp1..3  
!     based on the results from compressible code.  This is done only 
!     for wedges, depending on the boundary face.(tri or quad)  
!        
           if (lcsyst .eq. 1) then
              temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
              temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
              temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
           else 
              temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
              temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
              temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
           endif
!     
           temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
           bnorm(:,1) = temp1 * temp
           bnorm(:,2) = temp2 * temp
           bnorm(:,3) = temp3 * temp
!     
      
           if (lcsyst .eq. 3) then
              WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
           elseif (lcsyst .eq. 4) then
              WdetJb     = Qwtb(lcsyst,intp) / temp
           else
              WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
           endif

           if (numDirCalcSrfs.gt.zero) then
              ! An array to use to remember the locations that we dont want to zero out in bnorm as we find them
              allocate(locationsInbnormToZeroOut(npro))
              locationsInbnormToZeroOut = 1

              do iel=1,npro
                 do i=1, numDirCalcSrfs
                    if (btest(iBCB(iel,1),1) .or. iBCB(iel,2).eq.nsrflistDirCalc(i)) then
                       bnorm(iel,1:3)=bnorm(iel,1:3)*WdetJb(iel)
                       locationsInbnormToZeroOut(iel) = 0 ! Mark that we want to keep (i.e. dont want to zero out) bnorm for this element.
                       exit ! we've done what needs to be done on this iteration
                    ! else
                    !     bnorm(iel,:) = zero  ! we want zeros where we are not integrating
                    endif
                 enddo
              enddo

              do iel=1,npro
                if (locationsInbnormToZeroOut(iel) .eq. 1) then
                  bnorm(iel,:) = 0
                endif
              enddo

              deallocate(locationsInbnormToZeroOut)
           else
              do iel=1,npro
                 if (btest(iBCB(iel,1),1)) then 
                    bnorm(iel,1:3)=bnorm(iel,1:3)*WdetJb(iel)
                 else
                    bnorm(iel,:) = zero  ! we want zeros where we are not integrating
                 endif
              enddo
           endif

!
!  Now lets calculate Integral N_(a:e)^i n_i d Gamma
!
!
           do n = 1, nshlb
              nodlcl = lnode(n)
         rl(:,nodlcl,1) = rl(:,nodlcl,1) + shpfun(:,nodlcl) * bnorm(:,1)
         rl(:,nodlcl,2) = rl(:,nodlcl,2) + shpfun(:,nodlcl) * bnorm(:,2)
         rl(:,nodlcl,3) = rl(:,nodlcl,3) + shpfun(:,nodlcl) * bnorm(:,3) 
         
           enddo

        enddo  ! quadrature point loop
!
!.... assemble the NABI vector
!
        call local (NABI,    rl,     ienb,   3,  'scatter ')
!
!     push the surf number which we have associated with boundary
!     elements up to the global level in the array ndsurf
!
        do iel=1,npro

          if (indsurf) then        
          else
            if (iBCB(iel,2) .ne. 0) then
              iface = iBCB(iel,2)
              do k=1,nshlb
                 if (ndsurf(ienb(iel,k)).ne.1) then
                    ndsurf(ienb(iel,k))=iface   
                 endif
              enddo
            endif
          end if 
           
           if (numVisFluxSrfs .gt. zero) then
              do k=1, numVisFluxSrfs
                 if (iBCB(iel,2) .eq. nsrflistVisFlux(k)) then
                    iBCB(iel,1) = 6
                 endif
              enddo
           endif
        enddo
!     
!.... end
!
        return
        end


      subroutine AsBNASC ( x,       shpb, &
                         ienb,  iBCB)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles the data corresponding to the
!  boundary elements.
!
! Irene Vignon - copied from AsBNABI, Fall 2005.  (Fortran 90)
!----------------------------------------------------------------------
!
      use pvsQbi
      use phcommonvars  
      IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
        real*8, allocatable :: locationsInWdetJbToZeroOut(:)
!
        dimension xlb(npro,nenl,nsd), &
                  rl(npro,nshl),     WdetJb(npro)

        dimension x(numnp,nsd), &
                  shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb), &
                  ienb(npro,nshl), &
                  iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl), &
                  shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

!
        dimension dxdxib(npro,nsd,nsd),      temp(npro), &
                  temp1(npro),               temp2(npro), &
                  temp3(npro), &
                  v1(npro,nsd),              v2(npro,nsd)

!
!.... get the matrix of mode signs for the hierarchic basis functions
!
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
!
!.... gather the variables
!
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
!
!.... get the boundary element residuals
!
        rl  = zero
!
!.... compute the nodes which lie on the boundary (hierarchic)
!
        call getbnodes(lnode)
!
!.... loop through the integration points
!
        ngaussb = nintb(lcsyst)

        
        do intp = 1, ngaussb
!
!.... get the hierarchic shape functions at this int point
!
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, &
                    shpfun,       shdrv)

!
!.... compute the normal to the boundary. This is achieved by taking
!     the cross product of two vectors in the plane of the 2-d 
!     boundary face.
!
           if(lcsyst.ne.6) then
              ipt3=3
           else
              ipt3=5
           endif
           v1 = xlb(:,2,:) - xlb(:,1,:)
           v2 = xlb(:,ipt3,:) - xlb(:,1,:)

!
!.....The following are done in order to correct temp1..3  
!     based on the results from compressible code.  This is done only 
!     for wedges, depending on the boundary face.(tri or quad)  
!        
           if (lcsyst .eq. 1) then
              temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
              temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
              temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
           else 
              temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
              temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
              temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
           endif
!     
           temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
       
           if (lcsyst .eq. 3) then
              WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
           elseif (lcsyst .eq. 4) then
              WdetJb     = Qwtb(lcsyst,intp) / temp
           else
              WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
           endif


!......here I only want the d Gamma, not the n_i
           if (numDirCalcSrfs.gt.zero) then
              
              ! An array to use to remember the locations that we dont want to zero out in WdetJb as we find them
              allocate(locationsInWdetJbToZeroOut(npro))
              locationsInWdetJbToZeroOut = 1

              do iel=1,npro
                 do i=1, numDirCalcSrfs
                    if (btest(iBCB(iel,1),1) .or. iBCB(iel,2).eq.nsrflistDirCalc(i)) then
                       locationsInWdetJbToZeroOut(iel) = 0 ! Mark that we want to keep (i.e. dont want to zero out) WdetJb for this element.
                    endif
                 enddo
              enddo
              
              where(locationsInWdetJbToZeroOut .eq. 1) WdetJb = 0 ! we want zeros where we are not integrating

              deallocate(locationsInWdetJbToZeroOut)
           else
              do iel=1,npro
                 if (btest(iBCB(iel,1),1)) then 
                 else
                    WdetJb(iel) = zero  ! we want zeros where we are not integrating
                 endif
              enddo
           endif

!
!  Now lets calculate Integral N_(a:e) d Gamma
!
!
           do n = 1, nshlb
              nodlcl = lnode(n)
            rl(:,nodlcl) = rl(:,nodlcl) + shpfun(:,nodlcl) * WdetJb(:)
         
           enddo

        enddo  ! quadrature point loop
!
!.... assemble the NASC vector
!
        call local (NASC,    rl,     ienb,   1,  'scatter ')
!
!     push the surf number which we have associated with boundary
!     elements up to the global level in the array ndsurf
!
        if (indsurf) then        
        else
          do iel=1,npro
            if (iBCB(iel,2) .ne. 0) then
              iface = iBCB(iel,2)
              do k=1,nshlb
                if (ndsurf(ienb(iel,k)).ne.1) then
                  ndsurf(ienb(iel,k))=iface   
                endif
              enddo
            endif
          enddo
        end if 
!     
!.... end
!
        return
        end

      subroutine AsBPNABI ( x,       shpb, &
                         ienb,  iBCB)
!
!----------------------------------------------------------------------
!
! This routine computes and assembles data required for an Augmented
! Lagrangian Method. 
!
!----------------------------------------------------------------------
!
        use pvsQbi
        use LagrangeMultipliers ! brings in face radius and center 
        use phcommonvars  
        IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
!
        dimension xlb(npro,nenl,nsd),    bnorm(npro,nsd), &
                  rl(npro,nshl,nsd),     WdetJb(npro), &
                  rl2(npro,nshl,nsd,3)

        dimension x(numnp,nsd), &
                  shpb(nshl,ngaussb),      shglb(nsd,nshl,ngaussb), &
                  ienb(npro,nshl), &
                  iBCB(npro,ndiBCB)

        dimension lnode(27),               sgn(npro,nshl), &
                  shpfun(npro,nshl),        shdrv(npro,nsd,nshl)

!
        dimension dxdxib(npro,nsd,nsd),      temp(npro), &
                  temp1(npro),               temp2(npro), &
                  temp3(npro), &
                  v1(npro,nsd),              v2(npro,nsd)
!
        real*8    intpxlb(npro,3,nsd),    intpdistance(npro,3), &
                  intpprofile(npro,3)
        real*8    tmpLagInplaneVectors(3,3,0:MAXSURF)
        real*8    tmpLagProfileArea(0:MAXSURF)
        real*8    tmpProfileDelta(0:MAXSURF)
        real*8    Inplane1, Inplane2, Inplane3, InplaneNorm
        integer   count
!
!.... get the matrix of mode signs for the hierarchic basis functions
!
        if (ipord .gt. 1) then
           call getsgn(ienb,sgn)
        endif
!
!.... gather the variables
!
        call localx(x,      xlb,    ienb,   nsd,    'gather  ')
!
!....   calculate quadrature points
!        
        intpxlb = zero
        intpdistance = zero
        intpprofile = zero
        if (lcsyst.ne.6) then 
           do intp=1, 3 ! use 3 quadrature points
              do n=1, 3  
                 intpxlb(:,intp,:) = intpxlb(:,intp,:) &
                    +xlb(:,n,:)*Qptb(1,n,intp)
              enddo
           enddo
        else
           do intp=1, 3 ! use 3 quadrature points
              do n=1, 2  
                 intpxlb(:,intp,:) = intpxlb(:,intp,:) &
                    +xlb(:,n,:)*Qptb(1,n,intp)
              enddo
              intpxlb(:,intp,:) = intpxlb(:,intp,:) &
                 +xlb(:,5,:)*Qptb(1,3,intp)              
           enddo
        endif
!
!....   calculate profile functions at quadrature points
!        
        do k=1, numLagrangeSrfs
           do intp = 1, 3
              do iel=1, npro
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    intpdistance(iel,intp)=sqrt( &
                      (intpxlb(iel,intp,1)-LagCenter(1,k))**2+ &
                      (intpxlb(iel,intp,2)-LagCenter(2,k))**2+ &
                      (intpxlb(iel,intp,3)-LagCenter(3,k))**2) &
                      /LagRadius(k)
                    intpprofile(iel,intp)= &
                      (ProfileOrder(k)+2)/ProfileOrder(k)* &
                      (1-intpdistance(iel,intp)**ProfileOrder(k))
                 endif
              enddo
           enddo
         enddo  

!
!.... get the boundary element residuals
!
        rl  = zero
        rl2 = zero
        tmpLagProfileArea = zero
        tmpProfileDelta = zero     
        tmpLagInplaneVectors = zero
!
!.... compute the nodes which lie on the boundary (hierarchic)
!
        call getbnodes(lnode)
!
!.... loop through the integration points
!
        ngaussb = nintb(lcsyst)
!        
        do intp = 1, ngaussb
!
!.... get the hierarchic shape functions at this int point
!
           shglb=zero  ! protect debugger 
           call getshpb(shpb,        shglb,        sgn, &
                    shpfun,       shdrv)

!
!.... compute the normal to the boundary. This is achieved by taking
!     the cross product of two vectors in the plane of the 2-d 
!     boundary face.
!
           if(lcsyst.ne.6) then
              ipt3=3
           else
              ipt3=5
           endif
           v1 = xlb(:,2,:) - xlb(:,1,:)
           v2 = xlb(:,ipt3,:) - xlb(:,1,:)
!
!.....The following are done in order to correct temp1..3  
!     based on the results from compressible code.  This is done only 
!     for wedges, depending on the boundary face.(tri or quad)  
!        
           if (lcsyst .eq. 1) then
              temp1 = v1(:,2) * v2(:,3) - v2(:,2) * v1(:,3)
              temp2 = v2(:,1) * v1(:,3) - v1(:,1) * v2(:,3)
              temp3 = v1(:,1) * v2(:,2) - v2(:,1) * v1(:,2)
           else 
              temp1 = - v1(:,2) * v2(:,3) + v2(:,2) * v1(:,3)
              temp2 = - v2(:,1) * v1(:,3) + v1(:,1) * v2(:,3)
              temp3 = - v1(:,1) * v2(:,2) + v2(:,1) * v1(:,2)
           endif
!     
           temp       = one / sqrt ( temp1**2 + temp2**2 + temp3**2 )
           bnorm(:,1) = temp1 * temp
           bnorm(:,2) = temp2 * temp
           bnorm(:,3) = temp3 * temp
! 
           if (lcsyst .eq. 3) then
              WdetJb     = (1 - Qwtb(lcsyst,intp)) / (four*temp)
           elseif (lcsyst .eq. 4) then
              WdetJb     = Qwtb(lcsyst,intp) / temp
           else
              WdetJb     = Qwtb(lcsyst,intp) / (four*temp)
           endif
!
           do iel=1,npro
              count = 0
              do k=1, numLagrangeSrfs
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    iface = iBCB(iel,2)
                    do kk=1,nshlb
                      if (indsurf) then        
                      else
                        if (ndsurf(ienb(iel,kk)).ne.1) then
                          ndsurf(ienb(iel,kk))=iface   
                        endif
                      end if 
                    enddo
                    count = count+1
                 endif
              enddo
              if (count .eq. 0) then
                 bnorm(iel,:) = zero  ! we want zeros where we are not integrating
                 WdetJb(iel) = zero  ! we want zeros where we are not integrating
              endif              
           enddo
!
!   Calculate two orthonormal in-plane vectors
!   |bnorm(iel,1)  bnorm(iel,2)  bnorm(iel,3) |
!   |v1(iel,1)     v1(iel,2)     v1(iel,3)    | 
!   x1 component: -v1(iel,2)*bnorm(iel,3)+v1(iel,3)*bnorm(iel,2)
!   x2 component: -v1(iel,3)*bnorm(iel,1)+v1(iel,1)*bnorm(iel,3)
!   x3 component: -v1(iel,1)*bnorm(iel,2)+v1(iel,2)*bnorm(iel,1)
!
           do k=1, numLagrangeSrfs
              do iel=1,npro
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    tmpLagInplaneVectors(1:3,1,k)=bnorm(iel,1:3)
                    tmpLagInplaneVectors(1:3,2,k)=v1(iel,1:3) &
                       /sqrt(v1(iel,1)**2+v1(iel,2)**2+v1(iel,3)**2)
                    Inplane1=-v1(iel,2)*bnorm(iel,3) &
                       +v1(iel,3)*bnorm(iel,2)
                    Inplane2=-v1(iel,3)*bnorm(iel,1) &
                       +v1(iel,1)*bnorm(iel,3)
                    Inplane3=-v1(iel,1)*bnorm(iel,2) &
                       +v1(iel,2)*bnorm(iel,1)
                    InplaneNorm=one &
                       /sqrt(Inplane1**2+Inplane2**2+Inplane3**2)
                    tmpLagInplaneVectors(1,3,k)=Inplane1*InplaneNorm
                    tmpLagInplaneVectors(2,3,k)=Inplane2*InplaneNorm
                    tmpLagInplaneVectors(3,3,k)=Inplane3*InplaneNorm
                    exit
                 endif
              enddo
           enddo              
!
!  Now lets calculate Integral N_(a:e)^i n_i ProfileFunction  d Gamma
!
!
           do n = 1, nshlb
              nodlcl = lnode(n)
              rl(:,nodlcl,1) = rl(:,nodlcl,1) + shpfun(:,nodlcl) &
                 * bnorm(:,1)*intpprofile(:,intp)*WdetJb(:)
              rl(:,nodlcl,2) = rl(:,nodlcl,2) + shpfun(:,nodlcl) &
                 * bnorm(:,2)*intpprofile(:,intp)*WdetJb(:)
              rl(:,nodlcl,3) = rl(:,nodlcl,3) + shpfun(:,nodlcl) &
                 * bnorm(:,3)*intpprofile(:,intp)*WdetJb(:)
           enddo
!
!  Now lets calculate Integral N_(a:e)^i n_i N_(b:e)^i n_i d Gamma
!
!
           do k=1, numLagrangeSrfs
              do n = 1, nshlb
                 nodlcl = lnode(n)
                 do m=1, nsd                 
                    rl2(:,nodlcl,m,1)=rl2(:,nodlcl,m,1)+ &
                       shpfun(:,nodlcl)*shpfun(:,nodlcl)*WdetJb(:) &
                       *tmpLagInplaneVectors(m,1,k) &
                       *tmpLagInplaneVectors(m,1,k)
                    rl2(:,nodlcl,m,2)=rl2(:,nodlcl,m,2)+ &
                       shpfun(:,nodlcl)*shpfun(:,nodlcl)*WdetJb(:) &
                       *tmpLagInplaneVectors(m,2,k) &
                       *tmpLagInplaneVectors(m,2,k)
                    rl2(:,nodlcl,m,3)=rl2(:,nodlcl,m,3)+ &
                       shpfun(:,nodlcl)*shpfun(:,nodlcl)*WdetJb(:) &
                       *tmpLagInplaneVectors(m,3,k) &
                       *tmpLagInplaneVectors(m,3,k)
                 enddo
              enddo
           enddo
           
           do k=1, numLagrangeSrfs
              do iel=1,npro
                 if (iBCB(iel,2) .eq. nsrflistLagrange(k)) then
                    tmpLagProfileArea(k)=tmpLagProfileArea(k)+ &
                       intpprofile(iel,intp)*WdetJb(iel)
                    tmpProfileDelta(k)=tmpProfileDelta(k)+ &
                       intpprofile(iel,intp)**2*WdetJb(iel)
                 endif
              enddo
           enddo
        enddo  ! quadrature point loop
!
!.... assemble the PNABI vector
!
        call local (PNABI,    rl,     ienb,   3,  'scatter ')
!
!.... assemble the NANBIJ vector
!      
        do i=1, 3
           call local (NANBIJ(:,:,i),rl2(:,:,:,i),ienb,3,'scatter ')
        enddo
        
        do k=1, numLagrangeSrfs
           LagProfileArea(k)=LagProfileArea(k)+tmpLagProfileArea(k)
           ProfileDelta(k)=ProfileDelta(k)+tmpProfileDelta(k)
           InplaneNorm=sqrt(LagInplaneVectors(1,1,k)**2+ &
              LagInplaneVectors(2,1,k)**2+LagInplaneVectors(3,1,k)**2)
           if (InplaneNorm .eq. zero) then
              LagInplaneVectors(:,:,k)=tmpLagInplaneVectors(:,:,k)
           endif
        enddo   
!                 
        return
        end
