!---------------------------------------------------------------------
!     
!     drvftools.f : Bundle of Fortran driver routines for ftools.f
!     
!     Each routine is to be called by les**.c
!     
!---------------------------------------------------------------------
!     
!----------------
!     drvLesPrepDiag
!----------------
!     
subroutine drvlesPrepDiag ( flowDiag, ilwork, &
                            iBC,      BC,      iper, &
                            rowp,     colm,     &
                            lhsK,     lhsP)
    !
    use pointer_data
    use pvsQbi
    use convolImpFlow !brings in the current part of convol coef for imp BC
    use convolRCRFlow !brings in the current part of convol coef for RCR BC
    use convolTRCRFlow
    use convolCORFlow !brings in the current part of convol coef for Coronary BC
    use incpBC
    use LagrangeMultipliers

    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
    !
    dimension flowDiag(nshg,4), ilwork(nlwork)
    dimension iBC(nshg), iper(nshg), BC(nshg,ndofBC)
    real*8 lhsK(9,nnz_tot), lhsP(4,nnz_tot)
    integer rowp(nnz_tot),  colm(nshg+1)
    integer	n,	k
    !
    integer sparseloc
    !
    !
    !.... Clear the flowdiag
    !
    if((flmpl.eq.1).or.(ipord.gt.1)) then
        do n = 1, nshg
            k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n ) &
            + colm(n)-1
            !
            flowdiag(n,1) = lhsK(1,k)
            flowdiag(n,2) = lhsK(5,k)
            flowdiag(n,3) = lhsK(9,k)
            !
            flowdiag(n,4) = lhsP(4,k)
        enddo
    else
        flowDiag = zero
        do n = 1, nshg  ! rowsum put on the diagonal instead of diag entry
            do k=colm(n),colm(n+1)-1

                !
                flowdiag(n,1) = flowdiag(n,1) + abs(lhsK(1,k))
                !     &                          + lhsK(2,k) + lhsK(3,k)
                flowdiag(n,2) = flowdiag(n,2) + abs(lhsK(5,k))
                !     &                          + lhsK(4,k) + lhsK(6,k)
                flowdiag(n,3) = flowdiag(n,3) + abs(lhsK(9,k))
                !     &                          + lhsK(7,k) + lhsK(8,k)
                !
                flowdiag(n,4) = flowdiag(n,4) + abs(lhsP(4,k))
            enddo
            flowdiag(n,:)=flowdiag(n,:)*pt33
        enddo
    endif
    if(ipvsq.ge.3) then ! for first cut only do diagonal extraction
        ! this is not yet correct for multi procs I suspect if partition
        ! boundary cuts a p=QR face
        tfact=alfi * gami * Delt(1)
        do n=1,nshg
            if(numResistSrfs.gt.zero) then
                do k = 1,numResistSrfs
                    if (nsrflistResist(k).eq.ndsurf(n)) then
                        irankCoupled=k
                        flowdiag(n,1:3) = flowdiag(n,1:3) &
                        + tfact*ValueListResist(irankCoupled)* &
                        NABI(n,:)*NABI(n,:)
                    endif
                enddo
            elseif(numImpSrfs.gt.zero) then
                do k = 1,numImpSrfs
                    if (nsrflistImp(k).eq.ndsurf(n)) then
                        irankCoupled=k
                        flowdiag(n,1:3) = flowdiag(n,1:3) &
                        + tfact*ImpConvCoef(ntimeptpT+2,irankCoupled)* &
                        NABI(n,:)*NABI(n,:)
                    endif
                enddo
            elseif(numRCRSrfs.gt.zero) then
                do k = 1,numRCRSrfs
                    if (nsrflistRCR(k).eq.ndsurf(n)) then
                        irankCoupled=k
                        flowdiag(n,1:3) = flowdiag(n,1:3) &
                        + tfact*RCRConvCoef(currentTimestepIndex+2,irankCoupled)* & !check currentTimestepIndex+2 if restart from t.ne.0 &
                        NABI(n,:)*NABI(n,:)
                    endif
                enddo
            elseif(numTRCRSrfs.gt.zero) then
                do k = 1,numTRCRSrfs
                    if (nsrflistTRCR(k).eq.ndsurf(n)) then
                        irankCoupled=k
                        flowdiag(n,1:3)=flowdiag(n,1:3) &
                        + tfact*TRCRConvCoef(currentTimestepIndex+2,irankCoupled)* &
                        NABI(n,:)*NABI(n,:)
                  endif
               enddo
            elseif(numCORSrfs.gt.zero) then
                do k = 1,numCORSrfs
                    if (nsrflistCOR(k).eq.ndsurf(n)) then
                        irankCoupled=k
                        flowdiag(n,1:3) = flowdiag(n,1:3) &
                        + tfact*CORConvCoef(currentTimestepIndex+2,irankCoupled)* &
                        NABI(n,:)*NABI(n,:)
                    endif
                enddo
            elseif(numINCPSrfs.gt.zero) then
                do k = 1,numINCPSrfs
                    if (nsrflistINCP(k) .ne. inactive(k)) then
                        if (nsrflistINCP(k).eq.ndsurf(n)) then
                            irankCoupled=k
                            flowdiag(n,1:3) = flowdiag(n,1:3) &
                            + tfact*INCPCoef(1,irankCoupled)*  &
                            NABI(n,:)*NABI(n,:)
                        endif
                    endif
                enddo
            endif
        enddo
    endif
    !
    !.... Now I am adding contributions from the Lagrange multipliers to preconditioning
    !
    !
    if(Lagrange .gt. zero) then
        tfact=alfi * gami * Delt(1)
        call LagAddDiag (flowDiag, tfact)
    endif

    !
    if(iabc==1) &   !are there any axisym bc's &
    call rotabc(flowdiag, iBC, 'in ')
    !

    !
    !.... communicate : add the slaves part to the master's part of flowDiag
    !
    if (numpe > 1) then
        call commu (flowDiag, ilwork, nflow, 'in ')
    endif
    !
    !.... satisfy the boundary conditions on the diagonal
    !
    call bc3diag(iBC, BC,  flowDiag)
    !
    !
    !.... on processor periodicity was not taken care of in the setting of the
    !     boundary conditions on the matrix.  Take care of it now.
    !
    call bc3per(iBC,  flowDiag, iper, ilwork, 4)
    !
    !... slaves and masters have the correct values
    !
    !
    !.... Calculate square root
    !
    do i = 1, nshg
        do j = 1, nflow
            if (flowDiag(i,j).ne.0)  &
            flowDiag(i,j) = 1. / sqrt(abs(flowDiag(i,j)))
        enddo
    enddo
    !
    return
end

!     
!-------------
!     drvsclrDiag
!-------------
!     
subroutine drvsclrDiag ( sclrDiag, ilwork, iBC, BC, iper,  &
                         rowp,     colm,   lhsS )
    !
    use pointer_data
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
      include "mpif.h"
    !
    integer  ilwork(nlwork),    iBC(nshg),     iper(nshg), &
    rowp(nnz_tot),    colm(nshg+1)

    real*8   sclrDiag(nshg),    lhsS(nnz_tot), BC(nshg,ndofBC)
    integer sparseloc

    sclrDiag = zero
    do n = 1, nshg
        k = sparseloc( rowp(colm(n)), colm(n+1)-colm(n), n )  &
        + colm(n)-1
        !
        sclrDiag(n) = lhsS(k)
    enddo
    !
    !.... communicate : add the slaves part to the master's part of sclrDiag
    !
    if (numpe > 1) then
        call commu (sclrDiag, ilwork, 1, 'in ')
    endif
    !
    !.... satisfy the boundary conditions on the diagonal
    !
    call bc3SclrDiag(iBC,  sclrDiag)
    !
    !
    !.... on processor periodicity was not taken care of in the setting of the
    !     boundary conditions on the matrix.  Take care of it now.
    !
    call bc3per(iBC,  sclrDiag, iper, ilwork, 1)
    !
    !... slaves and masters have the correct values
    !
    !
    !.... Calculate square root
    !
    do i = 1, nshg
        if (sclrDiag(i).ne.0) then
            sclrDiag(i) = 1. / sqrt(abs(sclrDiag(i)))
        endif
    enddo
    !
    return
end

!============================================================================
!
! "fLesSparseApG":
!
!============================================================================
subroutine fLesSparseApG(	col,	row,	pLhs,	 &
                            p,      q,      nNodes, &
                            nnz_tot )
    !
    !.... Data declaration
    !
    implicit none
    integer	nNodes, nnz_tot
    integer	col(nNodes+1),	row(nnz_tot)
    real*8	pLhs(4,nnz_tot),	p(nNodes),	q(nNodes,3)
    !
    real*8	pisave
    integer	i,	j,	k
    !
    !.... clear the vector
    !
    do i = 1, nNodes
        q(i,1) = 0
        q(i,2) = 0
        q(i,3) = 0
    enddo

    !
    !.... Do an AP product
    !
    do i = 1, nNodes

        pisave = p(i)
        !dir$ ivdep
        do k = col(i), col(i+1)-1
            j = row(k)
            !
            q(j,1) = q(j,1) - pLhs(1,k) * pisave
            q(j,2) = q(j,2) - pLhs(2,k) * pisave
            q(j,3) = q(j,3) - pLhs(3,k) * pisave
        enddo
    enddo

    !
    !.... end
    !
    return
end

!============================================================================
!
! "fLesSparseApKG":
!
!============================================================================

subroutine fLesSparseApKG(	col,    row,    kLhs,   pLhs, &
                            p,      q,      nNodes, &
                            nnz_tot_hide )
    !
    !.... Data declaration
    !
    !	implicit none
    use pvsQbi
    use LagrangeMultipliers

    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    integer	nNodes !, nnz_tot
    integer	col(nNodes+1),	row(nnz_tot)
    real*8	kLhs(9,nnz_tot),	pLhs(4,nnz_tot)
    real*8 	p(nNodes,4),	q(nNodes,3)
    real*8	tmp1,	tmp2,	tmp3,	pisave
    integer	i,	j,	k
    !
    !.... clear the vector
    !
    do i = 1, nNodes
        q(i,1) = 0
        q(i,2) = 0
        q(i,3) = 0
    enddo

    !
    !.... Do an AP product
    !
    do i = 1, nNodes
        !
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        pisave   = p(i,4)
        !dir$ ivdep
        do k = col(i), col(i+1)-1
            j = row(k)
            tmp1 = tmp1 &
            + kLhs(1,k) * p(j,1) &
            + kLhs(4,k) * p(j,2) &
            + kLhs(7,k) * p(j,3)
            tmp2 = tmp2 &
            + kLhs(2,k) * p(j,1) &
            + kLhs(5,k) * p(j,2) &
            + kLhs(8,k) * p(j,3)
            tmp3 = tmp3 &
            + kLhs(3,k) * p(j,1) &
            + kLhs(6,k) * p(j,2) &
            + kLhs(9,k) * p(j,3)
            !
            q(j,1) = q(j,1) - pLhs(1,k) * pisave
            q(j,2) = q(j,2) - pLhs(2,k) * pisave
            q(j,3) = q(j,3) - pLhs(3,k) * pisave
        enddo
        q(i,1) = q(i,1) + tmp1
        q(i,2) = q(i,2) + tmp2
        q(i,3) = q(i,3) + tmp3
    enddo

    if (Lagrange .gt. 0) then
        LagSwitch = 1
        call CalcNANBLagrange(col, row, p(:,1:3))
    endif

    if(ipvsq.ge.2) then
        tfact=alfi * gami * Delt(1)
        call ElmpvsQ(q,p,tfact)
    endif
    !
    !.... end
    !
    return
end


!============================================================================
!
! "fLesSparseApNGt":
!
!============================================================================

subroutine fLesSparseApNGt(	col,	row,	pLhs,	&
                            p,	q,	nNodes, &
                            nnz_tot   )
    !
    !.... Data declaration
    !
    implicit none
    integer	nNodes, nnz_tot
    integer	col(nNodes+1),	row(nnz_tot)
    real*8	pLhs(4,nnz_tot),	p(nNodes,3),	q(nNodes)
    !
    real*8	tmp
    integer	i,	j,	k
    !
    !.... Do an AP product
    !
    do i = nNodes, 1, -1
        !
        tmp = 0
        do k = col(i), col(i+1)-1
            j = row(k)
            !
            tmp = tmp &
            + pLhs(1,k) * p(j,1) &
            + pLhs(2,k) * p(j,2) &
            + pLhs(3,k) * p(j,3)
        enddo
        q(i) = tmp
    enddo
    !
    !.... end
    !
    return
end

!============================================================================
!
! "fLesSparseApNGtC":
!
!============================================================================

subroutine fLesSparseApNGtC(	col,	row,	pLhs,	&
                                p,	q,	nNodes, &
                                nnz_tot )
    !
    !.... Data declaration
    !
    implicit none
    integer	nNodes, nnz_tot
    integer	col(nNodes+1),	row(20*nNodes)
    real*8	pLhs(4,nnz_tot),	p(nNodes,4),	q(nNodes)
    !
    real*8	tmp
    integer	i,	j,	k
    !
    !.... Do an AP product
    !
    do i = nNodes, 1, -1
        !
        tmp = 0
        do k = col(i), col(i+1)-1
            j = row(k)
            !
            tmp = tmp &
            + pLhs(1,k) * p(j,1) &
            + pLhs(2,k) * p(j,2) &
            + pLhs(3,k) * p(j,3) &
            + pLhs(4,k) * p(j,4)
        enddo
        q(i) = tmp
    enddo
    !
    !.... end
    !
    return
end

!============================================================================
!
! "fLesSparseApFull":
!
!============================================================================

subroutine fLesSparseApFull(	col,	row,	kLhs,	pLhs, &
                                p,	q,	nNodes, &
                                nnz_tot_hide )
    !
    !.... Data declaration
    !
    !	implicit none
    use pvsQbi
    use LagrangeMultipliers
        
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision

    integer	nNodes !, nnz_tot
    integer	col(nNodes+1),	row(nnz_tot)
    real*8	kLhs(9,nnz_tot),	pLhs(4,nnz_tot)
    real*8    p(nNodes,4),	q(nNodes,4)
    real*8	tmp1,	tmp2,	tmp3,	tmp4,	pisave
    integer	i,	j,	k
    !
    !.... clear the vector
    !
    do i = 1, nNodes
        q(i,1) = 0
        q(i,2) = 0
        q(i,3) = 0
    enddo
    !
    !.... Do an AP product
    !
    do i = 1, nNodes
        !
        tmp1 = 0
        tmp2 = 0
        tmp3 = 0
        tmp4 = 0
        pisave   = p(i,4)
        !dir$ ivdep
        do k = col(i), col(i+1)-1
            j = row(k)
            !
            tmp1 = tmp1 &
            + kLhs(1,k) * p(j,1) &
            + kLhs(4,k) * p(j,2) &
            + kLhs(7,k) * p(j,3)
            tmp2 = tmp2 &
            + kLhs(2,k) * p(j,1) &
            + kLhs(5,k) * p(j,2) &
            + kLhs(8,k) * p(j,3)
            tmp3 = tmp3 &
            + kLhs(3,k) * p(j,1) &
            + kLhs(6,k) * p(j,2) &
            + kLhs(9,k) * p(j,3)
            !
            tmp4 = tmp4 &
            + pLhs(1,k) * p(j,1) &
            + pLhs(2,k) * p(j,2) &
            + pLhs(3,k) * p(j,3) &
            + pLhs(4,k) * p(j,4)
            !
            q(j,1) = q(j,1) - pLhs(1,k) * pisave
            q(j,2) = q(j,2) - pLhs(2,k) * pisave
            q(j,3) = q(j,3) - pLhs(3,k) * pisave
        enddo
        q(i,1) = q(i,1) + tmp1
        q(i,2) = q(i,2) + tmp2
        q(i,3) = q(i,3) + tmp3
        q(i,4) = tmp4
    enddo
	
    if (Lagrange .gt. 0) then
        LagSwitch = 1
        call CalcNANBLagrange(col, row, p(:,1:3))
    endif

    if(ipvsq.ge.2) then
        tfact=alfi * gami * Delt(1)
        call ElmpvsQ(q,p,tfact)
    endif
    !
    !.... end
    !
    return
end

!============================================================================
!
! "fLesSparseApSclr":
!
!============================================================================

subroutine fLesSparseApSclr(	col,	row,	lhs,	&
                                p,	q,	nNodes, &
                                nnz_tot)
    !
    !.... Data declaration
    !
    implicit none
    integer	nNodes, nnz_tot
    integer	col(nNodes+1),	row(nnz_tot)
    real*8	lhs(nnz_tot),	p(nNodes),	q(nNodes)
    !
    real*8	tmp
    integer	i,	j,	k
    !
    !.... Do an AP product
    !
    do i = nNodes, 1, -1
        !
        tmp = 0
        do k = col(i), col(i+1)-1
            tmp = tmp + lhs(k) * p(row(k))
        enddo
        q(i) = tmp
    enddo
    !
    !.... end
    !
    return
end

!============================================================================
subroutine commOut(  global,  ilwork,  n,  &
                     iper,    iBC, BC  )
	
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
	
    real*8  global(nshg,n), BC(nshg,ndofBC)
    integer ilwork(nlwork), iper(nshg), iBC(nshg)
    !
    if ( numpe .gt. 1) then
        call commu ( global, ilwork, n, 'out')
    endif
    !
    !     before doing AP product P must be made periodic
    !     on processor slaves did not get updated with the
    !     commu (out) so do it here
    !
    do i=1,n
        global(:,i) = global(iper(:),i)  ! iper(i)=i if non-slave so no danger
    enddo
    !
    !       slave has masters value, for abc we need to rotate it
    !        (if this is a vector only no SCALARS)
    if((iabc==1) .and. (n.gt.1)) & !are there any axisym bc's &
    call rotabc(global, iBC,  'out')


    !$$$        do j = 1,nshg
    !$$$           if (btest(iBC(j),10)) then
    !$$$              i = iper(j)
    !$$$              res(j,:) = res(i,:)
    !$$$           endif
    !$$$        enddo
	
    return
end

!============================================================================
subroutine commIn(  global,  ilwork,  n,  &
                    iper,    iBC, BC )
	
    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
	
    real*8  global(nshg,n), BC(nshg,ndofBC)
    integer ilwork(nlwork), iper(nshg), iBC(nshg)
    !
    if((iabc==1) .and. (n.gt.1)) & !are there any axisym bc's &
    call rotabc(global, iBC, 'in ')
    !

    if ( numpe .gt. 1 ) then
        call commu ( global, ilwork, n, 'in ')
    endif
		
    call bc3per ( iBC, global, iper, ilwork, n)
	
    return
end


!============================================================================
subroutine LagAddDiag (flowdiag, tfact)

    use incpBC
    use LagrangeMultipliers
    use pvsQbi
    use multidomain, only: hrt

    use phcommonvars
    IMPLICIT REAL*8 (a-h,o-z)  ! change default real type to be double precision
    real*8    flowdiag(nshg,4), tfact, tfactSurf

    !
    !.... Now I am adding contributions from the Lagrange multipliers to preconditioning
    !
    !
    do n=1,nshg
        do k = 1,numLagrangeSrfs
            tfactSurf = zero
            tfactSurf = tfact * LagMeanFlow(k)
!!            if (numINCPSrfs .gt. zero) then
!!                if (nsrflistLagrange(k).eq.inactive(k)) then ! nsrflistLagrange should begin with nsrflistINCP
            if (iheart .gt. int(0)) then
                if (hrt%hassurfid(nsrflistLagrange(k)) .and. hrt%isavopen() .eq. int(0)) then    
                    Lag(k,:) = zero
                    Lagalpha(k,:) = zero
                else
                    if (nsrflistLagrange(k).eq.ndsurf(n)) then
                        do i=1, 3
                            flowdiag(n,i)=flowdiag(n,i)+ abs(  &
                            tfactSurf*(-Lagalpha(k,1) &
                            +PenaltyCoeff(k,1)*Penalty(k,1)) &
                            *(NANBIJ(n,i,1)-two*NABI(n,i)*PNABI(n,i) &
                            /LagProfileArea(k)+NABI(n,i)*NABI(n,i)* &
                            ProfileDelta(k))+tfactSurf &
                            *LagMeanFlow(k)*PenaltyCoeff(k,1) &
                            *(NANBLagrange(1,n,i)-PQLagrange(k,1)* &
                            NABI(n,i)-QLagrange(k,1)*PNABI(n,i) &
                            /LagProfileArea(k)+QLagrange(k,1)* &
                            NABI(n,i)*ProfileDelta(k))**2 &
                            +tfactSurf*(-Lagalpha(k,2) &
                            +PenaltyCoeff(k,2)*Penalty(k,2)) &
                            *NANBIJ(n,i,2)+tfactSurf*(-Lagalpha(k,3) &
                            +PenaltyCoeff(k,3)*Penalty(k,3)) &
                            *NANBIJ(n,i,3)+tfactSurf &
                            *LagMeanFlow(k)*(PenaltyCoeff(k,2)* &
                            NANBLagrange(2,n,i)**2+PenaltyCoeff(k,3) &
                            *NANBLagrange(3,n,i)**2) )
                            flowdiag(n,i)=flowdiag(n,i)+abs(tfactSurf**2 &
                            *( (NANBLagrange(1,n,i)- &
                            PQLagrange(k,1)*NABI(n,i)-QLagrange(k,1) &
                            *PNABI(n,i)/LagProfileArea(k)+ &
                            QLagrange(k,1)*NABI(n,i)*ProfileDelta(k))**2+ &
                            NANBLagrange(2,n,i)**2 &
                            +NANBLagrange(3,n,i)**2 ) &
                            /ScaleFactor(k,1)/two/gami/alfi)
                        enddo
                    endif
                endif
            else
                if (nsrflistLagrange(k).eq.ndsurf(n)) then
                    do i=1, 3
                        flowdiag(n,i)=flowdiag(n,i)+ abs(  &
                        tfactSurf*(-Lagalpha(k,1) &
                        +PenaltyCoeff(k,1)*Penalty(k,1)) &
                        *(NANBIJ(n,i,1)-two*NABI(n,i)*PNABI(n,i) &
                        /LagProfileArea(k)+NABI(n,i)*NABI(n,i)* &
                        ProfileDelta(k))+tfactSurf &
                        *LagMeanFlow(k)*PenaltyCoeff(k,1) &
                        *(NANBLagrange(1,n,i)-PQLagrange(k,1)* &
                        NABI(n,i)-QLagrange(k,1)*PNABI(n,i) &
                        /LagProfileArea(k)+QLagrange(k,1)* &
                        NABI(n,i)*ProfileDelta(k))**2 &
                        +tfactSurf*(-Lagalpha(k,2) &
                        +PenaltyCoeff(k,2)*Penalty(k,2)) &
                        *NANBIJ(n,i,2)+tfactSurf*(-Lagalpha(k,3) &
                        +PenaltyCoeff(k,3)*Penalty(k,3)) &
                        *NANBIJ(n,i,3)+tfactSurf &
                        *LagMeanFlow(k)*(PenaltyCoeff(k,2)* &
                        NANBLagrange(2,n,i)**2+PenaltyCoeff(k,3) &
                        *NANBLagrange(3,n,i)**2) )
                        flowdiag(n,i)=flowdiag(n,i)+abs(tfactSurf**2 &
                        *( (NANBLagrange(1,n,i)- &
                        PQLagrange(k,1)*NABI(n,i)-QLagrange(k,1) &
                        *PNABI(n,i)/LagProfileArea(k)+ &
                        QLagrange(k,1)*NABI(n,i)*ProfileDelta(k))**2+ &
                        NANBLagrange(2,n,i)**2 &
                        +NANBLagrange(3,n,i)**2 ) &
                        /ScaleFactor(k,1)/two/gami/alfi)
                    enddo
                endif
            endif
        enddo
    enddo

    return
end
