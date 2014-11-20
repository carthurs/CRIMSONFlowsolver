!--------------------------------------------------------------------
!     Created by Mahdi Esmaily Moghadam
!     contact memt63@gmail.com for reporting the bugs.
!--------------------------------------------------------------------
!
!     UC Copyright Notice
!     This software is Copyright ©2012 The Regents of the University of
!     California. All Rights Reserved.
!
!     Permission to copy and modify this software and its documentation
!     for educational, research and non-profit purposes, without fee,
!     and without a written agreement is hereby granted, provided that
!     the above copyright notice, this paragraph and the following three
!     paragraphs appear in all copies.
!
!     Permission to make commercial use of this software may be obtained
!     by contacting:
!     Technology Transfer Office
!     9500 Gilman Drive, Mail Code 0910
!     University of California
!     La Jolla, CA 92093-0910
!     (858) 534-5815
!     invent@ucsd.edu
!
!     This software program and documentation are copyrighted by The
!     Regents of the University of California. The software program and
!     documentation are supplied "as is", without any accompanying
!     services from The Regents. The Regents does not warrant that the
!     operation of the program will be uninterrupted or error-free. The
!     end-user understands that the program was developed for research
!     purposes and is advised not to rely exclusively on the program for
!     any reason.
!
!     IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
!     PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
!     DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS
!     SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
!     CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY
!     WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
!     OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
!     SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE
!     UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE
!     MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

module memLS

    USE ISO_FORTRAN_ENV
    IMPLICIT NONE
    INCLUDE "mpif.h"

    !     Communication parameters
    INTEGER, PARAMETER :: stdout = OUTPUT_UNIT
    INTEGER, PARAMETER :: mplog  = MPI_LOGICAL
    INTEGER, PARAMETER :: mpint  = MPI_INTEGER
    INTEGER, PARAMETER :: mpreal = MPI_DOUBLE_PRECISION
    INTEGER, PARAMETER :: mpchar = MPI_CHARACTER

    !     Some defenitions
    INTEGER, PARAMETER :: LS_TYPE_CG=1, LS_TYPE_GMRES=2, LS_TYPE_NS=3
    INTEGER, PARAMETER :: BC_TYPE_Dir = 0, BC_TYPE_Neu = 1
    INTEGER, PARAMETER :: BCOP_TYPE_ADD = 0, BCOP_TYPE_PRE = 1

    !     Communication structure
    TYPE memLS_commuType
        !        Free of created          (USE)
        LOGICAL :: foC = .FALSE.
        !        If this the master       (USE)
        LOGICAL masF
        !        Master ID                (USE)
        INTEGER master
        !        ID of this proc.         (USE)
        INTEGER task
        !        Task in FORTRAN indexing (USE)
        INTEGER tF
        !        Total number of tasks    (USE)
        INTEGER nTasks
        !        MPI communicator         (IN)
        INTEGER comm
    END TYPE memLS_commuType

    !     LHS matrix related data
    TYPE memLS_faceType
        !        ||Sai||**2D0                   (USE)
        REAL*8 nS
        !        Neu: P = res*Q                 (IN)
        REAL*8 :: res = 0D0
        !        nodal Sai for Neu              (IN)
        REAL*8, ALLOCATABLE :: val(:,:)
        !        Neu W*Sai                      (TMP)
        REAL*8, ALLOCATABLE :: valM(:,:)
        !        Free or created                (USE)
        LOGICAL :: foC = .FALSE.
        !        Neu: P/Q coupling              (USE)
        LOGICAL coupledFlag
        !        Neu: shared between proces     (USE)
        LOGICAL :: sharedFlag = .FALSE.
        !        Included in the computations   (IN)
        LOGICAL incFlag
        !        Number of nodes                (IN)
        INTEGER :: nNo = 0
        !        Degrees of freedom for val     (IN)
        INTEGER dof
        !        Dir/Neu                        (IN)
        INTEGER :: bGrp = BC_TYPE_Dir
        !        Only for data alignment
        INTEGER reserved
        !        Global node number             (IN)
        INTEGER, ALLOCATABLE :: glob(:)
    END TYPE memLS_faceType

    !     All following are in (USE)
    TYPE memLS_cSType
        !        Pointer to start of data for commu (only 2 proc shared points)
        INTEGER ptr
        !        Number of data to be commu  (only 2 proc shared points)
        INTEGER n
        !        Commu tag
        INTEGER tag
        !        Commu req
        INTEGER req
        !        Number of blocks for commu  (for 3 < proc shared points)
        INTEGER nBl
        !        Only for data alignment
        INTEGER reserved
        !        Pointer to beggining of each block (for 3 < proc shared points)
        INTEGER, ALLOCATABLE :: blPtr(:)
        !        Length of each block (for 3 < proc shared points)
        INTEGER, ALLOCATABLE :: blN(:)
    END TYPE memLS_cSType

    TYPE memLS_lhsType
        TYPE(memLS_commuType) commu
        TYPE(memLS_cSType), ALLOCATABLE :: cS(:)
        TYPE(memLS_faceType), ALLOCATABLE :: face(:)
        !        Free of created                     (USE)
        LOGICAL :: foC = .FALSE.
        !        Global number of nodes              (IN)
        INTEGER :: gnNo = 0
        !        Number of nodes                     (IN)
        INTEGER :: nNo = 0
        !        Number of non-zero in lhs           (IN)
        INTEGER :: nnz = 0
        !        Number of faces                     (IN)
        INTEGER :: nFaces = 0
        !        nNo of this proc                    (USE)
        INTEGER mynNo
        !        Column pointer                      (USE)
        INTEGER, ALLOCATABLE :: colPtr(:)
        !        Row pointer                         (USE)
        INTEGER, ALLOCATABLE :: rowPtr(:,:)
        !        Diagonal pointer                    (USE)
        INTEGER, ALLOCATABLE :: diagPtr(:)
        !        Mapping of nodes                    (USE)
        INTEGER, ALLOCATABLE :: map(:)
    END TYPE memLS_lhsType

    !     LS related structures
    TYPE memLS_subLsType
        !        Absolute tolerance            (IN)
        REAL*8 absTol
        !        Relative tolerance            (IN)
        REAL*8 relTol
        !        Initial norm of residual      (OUT)
        REAL*8 iNorm
        !        Final norm of residual        (OUT)
        REAL*8 fNorm
        !        Res. rduction in last itr.    (OUT)
        REAL*8 dB
        !        Calling duration              (OUT)
        REAL*8 callD
        !        Successful solving            (OUT)
        LOGICAL suc
        !        Maximum iteration             (IN)
        INTEGER mItr
        !        Space dimension               (IN)
        INTEGER sD
        !        Number of iteration           (OUT)
        INTEGER itr
        !        Number of Ax multiply         (OUT)
        INTEGER cM
        !        Number of |x| norms           (OUT)
        INTEGER cN
        !        Number of <x.y> dot products  (OUT)
        INTEGER cD
        !        Only for data alignment       (-)
        INTEGER reserve
    END TYPE memLS_subLsType

    TYPE memLS_lsType
        TYPE(memLS_subLsType) GM
        TYPE(memLS_subLsType) CG
        TYPE(memLS_subLsType) RI
        !        Free of created             (USE)
        LOGICAL :: foC = .FALSE.
        !        Which one of LS             (IN)
        INTEGER LS_type
        !        Contribution of mom. res.   (OUT)
        INTEGER Resm
        !        Contribution of cont. res.  (OUT)
        INTEGER Resc
    END TYPE memLS_lsType

contains

    SUBROUTINE ADDBCMUL(op_Type, nFaces, dof, nNo, mynNo, commu, face, X, Y)

        INTEGER, INTENT(IN) :: op_type, nFaces, dof, nNo, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_faceType), INTENT(IN) :: face(nFaces)
        REAL*8, INTENT(IN) :: X(dof, nNo)
        REAL*8, INTENT(INOUT) :: Y(dof, nNo)

        INTEGER faIn, i, a, Ac, nsd
        REAL*8 S!, DOTV
        REAL*8, ALLOCATABLE :: v(:,:), coef(:)

        ALLOCATE(coef(nFaces))

        IF (op_Type .EQ. BCOP_TYPE_ADD) THEN
            coef = face%res
        ELSE IF(op_Type .EQ. BCOP_TYPE_PRE) THEN
            coef = -face%res/(1D0 + face%res*face%nS)
        ELSE
            PRINT *, "op_Type is not defined"
            STOP
        END IF

        DO faIn=1, nFaces
            nsd = MIN(face(faIn)%dof,dof)
            IF (face(faIn)%coupledFlag) THEN
                IF (face(faIn)%sharedFlag) THEN
                    IF (.NOT.ALLOCATED(v)) ALLOCATE(v(dof,nNo))
                    v = 0D0
                    DO a=1, face(faIn)%nNo
                        Ac = face(faIn)%glob(a)
                        DO i=1, nsd
                            v(i,Ac) = face(faIn)%valM(i,a)
                        END DO
                    END DO
                    S = coef(faIn)*DOTV(dof, mynNo, commu, v, X)
                    DO a=1, face(faIn)%nNo
                        Ac = face(faIn)%glob(a)
                        DO i=1, nsd
                            Y(i,Ac) = Y(i,Ac) + v(i,Ac)*S
                        END DO
                    END DO
                ELSE
                    S = 0D0
                    DO a=1, face(faIn)%nNo
                        Ac = face(faIn)%glob(a)
                        DO i=1, nsd
                            S = S + face(faIn)%valM(i,a)*X(i,Ac)
                        END DO
                    END DO
                    S = coef(faIn)*S
                    DO a=1, face(faIn)%nNo
                        Ac = face(faIn)%glob(a)
                        DO i=1, nsd
                            Y(i,Ac) = Y(i,Ac) + face(faIn)%valM(i,a)*S
                        END DO
                    END DO
                END IF
            END IF
        END DO

        RETURN
    END SUBROUTINE ADDBCMUL

    !====================================================================

    SUBROUTINE memLS_BC_CREATE (lhs, faIn, nNo, dof, BC_type, gNodes, Val)

        TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
        INTEGER, INTENT(IN) :: faIn, nNo, dof
        INTEGER, INTENT(IN) :: BC_type
        INTEGER, INTENT(IN) :: gNodes(nNo)
        REAL*8, INTENT(IN), OPTIONAL :: Val(dof,nNo)

        INTEGER a, Ac, i
        REAL*8, ALLOCATABLE :: v(:,:)

        IF (faIn .GT. lhs%nFaces) THEN
            PRINT *, "faIn is exceeding lhs structure maximum number of", "face:", lhs%nFaces, ">", faIn
            STOP
        END IF
        IF (faIn .LT. 0) THEN
            PRINT *, "faIn is should be greater than zero"
            STOP
        END IF

        IF (lhs%face(faIn)%foC) THEN
            PRINT *, "BC(", faIn,") is not free"
            PRINT *, "You may use memLS_BC_FREE to free this structure"
        END IF

        lhs%face(faIn)%nNo  = nNo
        lhs%face(faIn)%dof  = dof
        lhs%face(faIn)%bGrp = BC_type

        if ((.not.allocated(lhs%face(faIn)%glob)).and.(.not.allocated(lhs%face(faIn)%val)).and.(.not.allocated(lhs%face(faIn)%valM))) then
          ALLOCATE(lhs%face(faIn)%glob(nNo), lhs%face(faIn)%val(dof,nNo), lhs%face(faIn)%valM(dof,nNo))
        endif

        DO a=1, nNo
            Ac = lhs%map(gNodes(a))
            lhs%face(faIn)%glob(a) = Ac
        END DO

        IF (PRESENT(Val)) THEN
            DO a=1, nNo
                lhs%face(faIn)%val(:,a) = Val(:,a)
            END DO
        ELSE
            lhs%face(faIn)%val = 0D0
        END IF

        IF (lhs%commu%nTasks .GT. 1) THEN
            a = 0
            IF (lhs%face(faIn)%nNo .NE. 0) a = 1
            CALL MPI_ALLREDUCE(a, Ac, 1, mpint, MPI_SUM, lhs%commu%comm, i)
            IF (Ac .GT. 1) THEN
                lhs%face(faIn)%sharedFlag = .TRUE.
                IF (.NOT.ALLOCATED(v)) ALLOCATE(v(dof,lhs%nNo))
                v = 0D0
                DO a=1, nNo
                    Ac = lhs%face(faIn)%glob(a)
                    v(:,Ac) = lhs%face(faIn)%val(:,a)
                END DO
                CALL COMMUV(dof, lhs%nNo, lhs%commu, lhs%cS, v)

                DO a=1, nNo
                    Ac = lhs%face(faIn)%glob(a)
                    lhs%face(faIn)%val(:,a) = v(:,Ac)
                END DO
            END IF
        END IF

        RETURN
    END SUBROUTINE memLS_BC_CREATE

    !====================================================================

    SUBROUTINE memLS_BC_FREE (lhs, faIn)

        TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
        INTEGER, INTENT(IN) :: faIn

        IF (.NOT.lhs%face(faIn)%foC) THEN
            PRINT *, 'Cannot free face:', faIn
            PRINT *, 'It is not created yet'
            STOP
        END IF
        lhs%face(faIn)%foC        = .FALSE.
        lhs%face(faIn)%nNo        = 0
        lhs%face(faIn)%bGrp       = BC_TYPE_Dir
        lhs%face(faIn)%res        = 0D0
        lhs%face(faIn)%sharedFlag = .FALSE.

        DEALLOCATE(lhs%face(faIn)%glob, lhs%face(faIn)%val,lhs%face(faIn)%valM)

        RETURN
    END SUBROUTINE memLS_BC_FREE

    !====================================================================

    SUBROUTINE CGRAD(nFaces, dof, nNo, nnz, mynNo, commu, cS, face, ls, rowPtr, colPtr, D, G, L, R)

        INTEGER, INTENT(IN) :: nFaces, dof, nNo, nnz, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_faceType), INTENT(IN) :: face(nFaces)
        TYPE(memLS_subLsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: D(dof,nnz), G(dof,nnz), L(nnz)
        REAL*8, INTENT(INOUT) :: R(nNo)

        INTEGER i
        REAL*8 errO, err, alpha, eps, time
        !REAL*8 CPUT, NORMS, DOTS
        REAL*8, ALLOCATABLE :: X(:), P(:), SP(:), DGP(:), GP(:,:), unCondU(:,:)

        ALLOCATE(X(nNo), P(nNo), SP(nNo), DGP(nNo), GP(dof,nNo), unCondU(dof,nNo))

        time     = CPUT()
        ls%suc   = .FALSE.
        errO     = NORMS(mynNo, commu, R)
        ls%iNorm = errO
        eps      = MAX(ls%absTol,ls%relTol*errO)
        eps      = eps*eps
        X        = 0D0

        err      = errO
        err      = err*err
        P        = R
        DO i=1, ls%mItr
            IF (err .LT. eps) THEN
                ls%suc = .TRUE.
                EXIT
            END IF
            errO = err
            CALL SPARMULSV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, G, P, GP)

            IF (ANY(face%coupledFlag)) THEN
                unCondU = GP
                CALL ADDBCMUL(BCOP_TYPE_PRE, nFaces, dof, nNo, mynNo, commu, face, unCondU, GP)
            END IF

            CALL SPARMULVS(dof, nNo, nnz, commu, cS, rowPtr, colPtr, D, GP, DGP)

            CALL SPARMULSS(     nNo, nnz, commu, cS, rowPtr, colPtr, L, P, SP)

            SP    = SP - DGP
            alpha = errO/DOTS(mynNo, commu, P, SP)
            X     = X + alpha*P
            R     = R - alpha*SP
            err   = NORMS(mynNo, commu, R)
            err   = err*err
            P     = R + err/errO*P
        END DO
        R        = X
        ls%fNorm = SQRT(err)
        ls%callD = CPUT() - time + ls%callD
        ls%dB    = 5D0*LOG(err/errO)
        ls%itr   = ls%itr + i - 1

        DEALLOCATE(X, P, SP, DGP, GP)

        RETURN
    END SUBROUTINE CGRAD

    !====================================================================

    SUBROUTINE CGRADV(dof, nNo, nnz, mynNo, commu, cS, ls, rowPtr, colPtr, K, R)

        INTEGER, INTENT(IN) :: dof, nNo, nnz, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_subLsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: K(dof*dof,nnz)
        REAL*8, INTENT(INOUT) :: R(dof,nNo)

        INTEGER i
        REAL*8 errO, err, alpha, eps
        !REAL*8 CPUT, NORMV, DOTV
        REAL*8, ALLOCATABLE :: P(:,:), KP(:,:), X(:,:)

        ALLOCATE(P(dof,nNo), KP(dof,nNo), X(dof,nNo))

        ls%callD = CPUT()
        ls%suc   = .FALSE.
        err      = NORMV(dof, mynNo, commu, R)
        ls%iNorm = err
        eps      = MAX(ls%absTol,ls%relTol*err)
        eps      = eps*eps
        err      = err*err
        X        = 0D0
        P        = R

        DO i=1, ls%mItr
            IF (err .LT. eps) THEN
                ls%suc = .TRUE.
                EXIT
            END IF
            errO = err
            CALL SPARMULVV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, K, P, KP)

            alpha = errO/DOTV(dof, mynNo, commu, P, KP)
            X     = X + alpha*P
            R     = R - alpha*KP
            err   = NORMV(dof, mynNo, commu, R)
            err   = err*err
            P = R + err/errO*P
        END DO

        R        = X
        ls%itr   = i - 1
        ls%fNorm = SQRT(err)
        ls%callD = CPUT() - ls%callD
        ls%dB    = 5D0*LOG(err/errO)
        IF (i .GT. ls%mItr) ls%itr = ls%mItr

        RETURN
    END SUBROUTINE CGRADV

    !====================================================================

    SUBROUTINE CGRADS(nNo, nnz, mynNo, commu, cS, ls, rowPtr, colPtr, K, R)

        INTEGER, INTENT(IN) :: nNo, nnz, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_subLsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: K(nnz)
        REAL*8, INTENT(INOUT) :: R(nNo)

        INTEGER i
        REAL*8 errO, err, alpha, eps
        !REAL*8 CPUT, NORMS, DOTS
        REAL*8, ALLOCATABLE :: P(:), KP(:), X(:)

        ALLOCATE(P(nNo), KP(nNo), X(nNo))

        ls%callD = CPUT()
        ls%suc   = .FALSE.
        err      = NORMS(mynNo, commu, R)
        ls%iNorm = err
        eps      = MAX(ls%absTol,ls%relTol*err)
        eps      = eps*eps
        err      = err*err
        X        = 0D0
        P        = R

        DO i=1, ls%mItr
            IF (err .LT. eps) THEN
                ls%suc = .TRUE.
                EXIT
            END IF
            errO = err
            CALL SPARMULSS(nNo, nnz, commu, cS, rowPtr, colPtr, K, P, KP)
            alpha = errO/DOTS(mynNo, commu, P, KP)
            X     = X + alpha*P
            R     = R - alpha*KP
            err   = NORMS(mynNo, commu, R)
            err   = err*err
            P = R + err/errO*P
        END DO

        R        = X
        ls%itr   = i - 1
        ls%fNorm = SQRT(err)
        ls%callD = CPUT() - ls%callD
        ls%dB    = 5D0*LOG(err/errO)
        IF (i .GT. ls%mItr) ls%itr = ls%mItr

        RETURN
    END SUBROUTINE CGRADS

    !====================================================================

    ! memLS_COMMU_CREATE
    !
    ! creation subroutine for memLS_commuType,
    ! the datatype for an object containing MPI communicator info
    !
    ! commu : Handle to the communicator structure
    ! commi : Handle to the MPI communicator

    SUBROUTINE memLS_COMMU_CREATE(commu, commi)

        TYPE(memLS_commuType), INTENT(INOUT) :: commu
        INTEGER, INTENT(IN) :: commi

        INTEGER ierr
        INTEGER comm

        IF (commu%foC) THEN
            PRINT *, "COMMU is not free"
            PRINT *, "You may use memLS_COMMU_FREE to free this structure"
        END IF

        commu%foC  = .TRUE.
        commu%comm = commi

        comm       = commi

        CALL MPI_COMM_RANK(comm, commu%task, ierr)
        CALL MPI_COMM_SIZE(comm, commu%nTasks, ierr)

        CALL MPI_ALLREDUCE(commu%task, commu%master, 1, mpint, MPI_MIN, comm, ierr)

        IF (commu%master .NE. 0) THEN
            PRINT *, "master is not zero"
            CALL MPI_FINALIZE(comm, ierr)
            STOP
        END IF

        commu%masF  = .FALSE.
        commu%tF    = commu%task + 1
        IF (commu%task .EQ. commu%master) THEN
            commu%masF = .TRUE.
        END IF

        RETURN
    END SUBROUTINE memLS_COMMU_CREATE

    !====================================================================

    SUBROUTINE memLS_COMMU_FREE(commu)

        TYPE(memLS_commuType), INTENT(INOUT) :: commu

        IF (.NOT.commu%foC) THEN
            PRINT *, 'Cannot free commu'
            PRINT *, 'It is not created yet'
            STOP
        END IF
        commu%foC  = .FALSE.

        RETURN
    END SUBROUTINE memLS_COMMU_FREE

    !====================================================================

    SUBROUTINE GE (N, A, B)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N
        REAL*8, INTENT(IN) :: A(N,N)
        REAL*8, INTENT(INOUT) :: B(N)

        INTEGER m, ipv, i, j
        REAL*8 pivot, saveEl
        REAL*8, ALLOCATABLE :: C(:,:)

        IF (N .EQ. 2) THEN
            pivot  = A(1,1)*A(2,2) - A(2,1)*A(1,2)
            saveEl = (B(1)*A(2,2) - B(2)*A(1,2))/pivot
            B(2)   = (B(2)*A(1,1) - B(1)*A(2,1))/pivot
            B(1)   = saveEl
            RETURN
        END IF
        ALLOCATE(C(N,N+1))

        C(:N,:N) = A
        C(:,N+1) = B

        DO m=1,N-1
            ipv = m
            pivot = ABS(C(m,m))
            DO i=m+1,N
                IF (ABS(C(i,m)) .GT. pivot) THEN
                    ipv = i
                    pivot = ABS(C(i,m))
                END IF
            END DO
            IF (pivot .LT. 2D0*EPSILON(pivot)) THEN
                PRINT *, 'Singular matrix'
                STOP
            END IF
            IF (ipv .NE. m) THEN
                DO j=m, N+1
                    saveEl = C(m,j)
                    C(m,j) = C(ipv,j)
                    C(ipv,j) = saveEl
                END DO
                DO j=1, m-1
                END DO
            END IF

            DO i=m+1,N
                saveEl = C(i,m)/C(m,m)
                C(i,m) = 0D0
                DO j=m+1,N+1
                    C(i,j) = C(i,j) - saveEl*C(m,j)
                END DO
            END DO
        END DO

        DO j=N,1,-1
            DO i=j+1,N
                C(j,N+1) = C(j,N+1) - C(j,i)*C(i,N+1)
            END DO
            C(j,N+1) = C(j,N+1)/C(j,j)
        END DO

        B = C(:,N+1)

        if (allocated(C)) then
          deallocate(C)
        endif

        RETURN
    END SUBROUTINE GE

    !====================================================================

    SUBROUTINE GMRES(nFaces, dof, nNo, nnz, mynNo, commu, cS, face, ls, rowPtr, colPtr, Val, R, X)

        INTEGER, INTENT(IN) :: nFaces, dof, nNo, nnz, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_faceType), INTENT(IN) :: face(nFaces)
        TYPE(memLS_subLsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: Val(dof*dof,nnz), R(dof,nNo)
        REAL*8, INTENT(OUT) :: X(dof,nNo)

        INTEGER i, j, k, l
        !REAL*8 CPUT, NORMV, DOTV
        REAL*8 eps, tmp, time, y(ls%sD), c(ls%sD), s(ls%sD), err(ls%sD+1)
        REAL*8, ALLOCATABLE :: u(:,:,:), h(:,:), unCondU(:,:)

        ALLOCATE(h(ls%sD+1,ls%sD), u(dof,nNo,ls%sD+1), unCondU(dof,nNo))

        time   = CPUT()
        ls%suc = .FALSE.

        X = 0D0
        DO l=1, ls%mItr
            IF (l .EQ. 1) THEN
                u(:,:,1) = R
            ELSE
                ls%itr = ls%itr + 1
                CALL SPARMULVV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, Val, X, u(:,:,1))

                CALL ADDBCMUL(BCOP_TYPE_ADD, nFaces, dof, nNo, mynNo, commu, face, X, u(:,:,1))

                u(:,:,1) = R - u(:,:,1)
            END IF
            IF (ANY(face%coupledFlag)) THEN
                unCondU = u(:,:,1)
                CALL ADDBCMUL(BCOP_TYPE_PRE, nFaces, dof, nNo, mynNo, commu, face, unCondU, u(:,:,1))
            END IF

            err(1)   = NORMV(dof, mynNo, commu, u(:,:,1))
            IF (l .EQ. 1) THEN
                eps       = err(1)
                IF (eps .LE. ls%absTol) THEN
                    ls%callD = 0D0
                    ls%dB    = 0D0
                    RETURN
                END IF
                ls%iNorm  = eps
                ls%fNorm  = eps
                eps       = MAX(ls%absTol,ls%relTol*eps)
            END IF
            ls%dB = ls%fNorm
            u(:,:,1) = u(:,:,1)/err(1)
            DO i=1, ls%sD
                ls%itr = ls%itr + 1
                CALL SPARMULVV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, Val, u(:,:,i), u(:,:,i+1))

                CALL ADDBCMUL(BCOP_TYPE_ADD, nFaces, dof, nNo, mynNo, commu, face, u(:,:,i), u(:,:,i+1))

                IF (ANY(face%coupledFlag)) THEN
                    unCondU = u(:,:,i+1)
                    CALL ADDBCMUL(BCOP_TYPE_PRE, nFaces, dof, nNo, mynNo, commu, face, unCondU, u(:,:,i+1))
                END IF
                DO j=1, i
                    h(j,i) = DOTV(dof, mynNo, commu, u(:,:,i+1), u(:,:,j))
                    u(:,:,i+1) = u(:,:,i+1) - h(j,i)*u(:,:,j)
                END DO
                h(i+1,i)   = NORMV(dof, mynNo, commu, u(:,:,i+1))

                u(:,:,i+1) = u(:,:,i+1)/h(i+1,i)
                DO j=1, i-1
                    tmp      =  c(j)*h(j,i) + s(j)*h(j+1,i)
                    h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i)
                    h(j,i)   =  tmp
                END DO
                tmp      = SQRT(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
                c(i)     = h(i,i)/tmp
                s(i)     = h(i+1,i)/tmp
                h(i,i)   = tmp
                h(i+1,i) = 0D0
                err(i+1) = -s(i)*err(i)
                err(i)   =  c(i)*err(i)
                IF (ABS(err(i+1)) .LT. eps) THEN
                    ls%suc = .TRUE.
                    EXIT
                END IF
            END DO
            IF (i .GT. ls%sD) i = ls%sD

            y = err(1:i)
            DO j=i, 1, -1
                DO k=j+1, i
                    y(j) = y(j) - h(j,k)*y(k)
                END DO
                y(j) = y(j)/h(j,j)
            END DO

            DO j=1, i
                X = X + u(:,:,j)*y(j)
            END DO
            ls%fNorm = ABS(err(i+1))
            IF (ls%suc) EXIT
        END DO

        ls%callD = CPUT() - time + ls%callD
        ls%dB    = 1D1*LOG(ls%fNorm/ls%dB)

        RETURN
    END SUBROUTINE GMRES

    !====================================================================

    SUBROUTINE GMRESV(nFaces, dof, nNo, nnz, mynNo, commu, cS, face, ls, rowPtr, colPtr, Val, R)

        INTEGER, INTENT(IN) :: nFaces, dof, nNo, nnz, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_faceType), INTENT(IN) :: face(nFaces)
        TYPE(memLS_subLsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: Val(dof*dof,nnz)
        REAL*8, INTENT(INOUT) :: R(dof,nNo)

        INTEGER i, j, k, l
        !REAL*8 CPUT, NORMV, DOTV
        REAL*8 eps, tmp, y(ls%sD), c(ls%sD), s(ls%sD), err(ls%sD+1)
        REAL*8, ALLOCATABLE :: u(:,:,:), h(:,:), X(:,:)


        ALLOCATE(h(ls%sD+1,ls%sD), u(dof,nNo,ls%sD+1), X(dof,nNo))

        ls%callD  = CPUT()
        ls%suc    = .FALSE.
        eps       = NORMV(dof, mynNo, commu, R)
        ls%iNorm  = eps
        ls%fNorm  = eps
        eps       = MAX(ls%absTol,ls%relTol*eps)
        ls%itr    = 0
        X         = 0D0

        IF (ls%iNorm .LE. ls%absTol) THEN
            ls%callD = 0D0
            ls%dB    = 0D0
            RETURN
        END IF
        DO l=1, ls%mItr
            ls%dB = ls%fNorm
            ls%itr = ls%itr + 1
            CALL SPARMULVV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, Val, X, u(:,:,1))
            CALL ADDBCMUL(BCOP_TYPE_ADD, nFaces, dof, nNo, mynNo, commu, face, X, u(:,:,1))

            u(:,:,1) = R - u(:,:,1)
            err(1)   = NORMV(dof, mynNo, commu, u(:,:,1))
            u(:,:,1) = u(:,:,1)/err(1)
            DO i=1, ls%sD
                ls%itr = ls%itr + 1
                CALL SPARMULVV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, Val, u(:,:,i), u(:,:,i+1))
                CALL ADDBCMUL(BCOP_TYPE_ADD, nFaces, dof, nNo, mynNo, commu, face, u(:,:,i), u(:,:,i+1))

                DO j=1, i
                    h(j,i) = DOTV(dof, mynNo, commu, u(:,:,i+1), u(:,:,j))
                    u(:,:,i+1) = u(:,:,i+1) - h(j,i)*u(:,:,j)
                END DO
                h(i+1,i)   = NORMV(dof, mynNo, commu, u(:,:,i+1))
                u(:,:,i+1) = u(:,:,i+1)/h(i+1,i)
                DO j=1, i-1
                    tmp      =  c(j)*h(j,i) + s(j)*h(j+1,i)
                    h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i)
                    h(j,i)   =  tmp
                END DO
                tmp      = SQRT(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
                c(i)     = h(i,i)/tmp
                s(i)     = h(i+1,i)/tmp
                h(i,i)   = tmp
                h(i+1,i) = 0D0
                err(i+1) = -s(i)*err(i)
                err(i)   =  c(i)*err(i)
                IF (ABS(err(i+1)) .LT. eps) THEN
                    ls%suc = .TRUE.
                    EXIT
                END IF
            END DO
            IF (i .GT. ls%sD) i = ls%sD

            y = err(1:i)
            DO j=i, 1, -1
                DO k=j+1, i
                    y(j) = y(j) - h(j,k)*y(k)
                END DO
                y(j) = y(j)/h(j,j)
            END DO

            DO j=1, i
                X = X + u(:,:,j)*y(j)
            END DO
            ls%fNorm = ABS(err(i+1))
            IF (ls%suc) EXIT
        END DO
        R = X
        ls%callD = CPUT() - ls%callD
        ls%dB    = 1D1*LOG(ls%fNorm/ls%dB)

        RETURN
    END SUBROUTINE GMRESV

    !====================================================================

    SUBROUTINE COMMUV(dof, nNo, commu, cS, R)

        INTEGER, INTENT(IN) :: dof, nNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        REAL*8, INTENT(INOUT) :: R(dof,nNo)

        INTEGER i, j, k, s, e, ierr, nTasks, tF, stat(MPI_STATUS_SIZE)
        INTEGER comm
        REAL*8, ALLOCATABLE :: rTmp(:,:)

        IF (commu%nTasks .EQ. 1) RETURN

        nTasks = commu%nTasks
        tF     = commu%tF
        comm   = commu%comm

        IF (tF .NE. 1) THEN
            i = tF - 1
            i = cS(i)%ptr + cS(i)%n - 1
            ALLOCATE(rTmp(dof,i))
        END IF

        DO i=1, nTasks
            IF (cS(i)%tag .NE. 0) THEN
                s = cS(i)%ptr
                e = s + cS(i)%n - 1
                IF (i .LT. tF) THEN
                    CALL MPI_IRECV(rTmp(:,s:e), cS(i)%n*dof, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                ELSE
                    CALL MPI_ISEND(R(:,s:e), cS(i)%n*dof, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                END IF
            END IF
        END DO

        k = 1
        DO i=1, tF - 1
            IF (cS(i)%tag .NE. 0) THEN
                CALL MPI_WAIT(cS(i)%req, stat, ierr)
                DO j=1, cS(i)%nBl
                    s = cS(i)%blPtr(j)
                    e = s + cS(i)%blN(j) - 1
                    R(:,s:e) = R(:,s:e) + rTmp(:,k:k+e-s)
                    k = k + cS(i)%blN(j)
                END DO
            END IF
        END DO

        k = 1
        DO i=1, tF - 1
            DO j=1, cS(i)%nBl
                s = cS(i)%blPtr(j)
                e = s + cS(i)%blN(j) - 1
                rTmp(:,k:k+e-s) = R(:,s:e)
                k = k + cS(i)%blN(j)
            END DO
        END DO

        DO i=1, nTasks
            IF (cS(i)%tag .NE. 0) THEN
                s = cS(i)%ptr
                e = s + cS(i)%n - 1
                IF (i .GT. tF) THEN
                    CALL MPI_WAIT(cS(i)%req, stat, ierr)
                    CALL MPI_IRECV(R(:,s:e), cS(i)%n*dof, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                ELSE
                    CALL MPI_ISEND(rTmp(:,s:e), cS(i)%n*dof, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                END IF
            END IF
        END DO

        DO i=1, nTasks
            IF (cS(i)%tag .NE. 0) THEN
                CALL MPI_WAIT(cS(i)%req, stat, ierr)
            END IF
        END DO

        RETURN
    END SUBROUTINE COMMUV

    !====================================================================

    SUBROUTINE COMMUS(nNo, commu, cS, R)

        INTEGER, INTENT(IN) :: nNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        REAL*8, INTENT(INOUT) :: R(nNo)

        INTEGER i, j, k, s, e, ierr, nTasks, tF, stat(MPI_STATUS_SIZE)
        INTEGER comm
        REAL*8, ALLOCATABLE :: rTmp(:)

        IF (commu%nTasks .EQ. 1) RETURN

        nTasks = commu%nTasks
        tF     = commu%tF
        comm   = commu%comm

        IF (tF .NE. 1) THEN
            i = tF - 1
            i = cS(i)%ptr + cS(i)%n - 1
            ALLOCATE(rTmp(i))
        END IF

        DO i=1, nTasks
            IF (cS(i)%tag .NE. 0) THEN
                s = cS(i)%ptr
                e = s + cS(i)%n - 1
                IF (i .LT. tF) THEN
                    CALL MPI_IRECV(rTmp(s:e), cS(i)%n, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                ELSE
                    CALL MPI_ISEND(R(s:e), cS(i)%n, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                END IF
            END IF
        END DO

        k = 1
        DO i=1, tF - 1
            IF (cS(i)%tag .NE. 0) THEN
                CALL MPI_WAIT(cS(i)%req, stat, ierr)
                DO j=1, cS(i)%nBl
                    s = cS(i)%blPtr(j)
                    e = s + cS(i)%blN(j) - 1
                    R(s:e) = R(s:e) + rTmp(k:k+e-s)
                    k = k + cS(i)%blN(j)
                END DO
            END IF
        END DO

        k = 1
        DO i=1, tF - 1
            DO j=1, cS(i)%nBl
                s = cS(i)%blPtr(j)
                e = s + cS(i)%blN(j) - 1
                rTmp(k:k+e-s) = R(s:e)
                k = k + cS(i)%blN(j)
            END DO
        END DO

        DO i=1, nTasks
            IF (cS(i)%tag .NE. 0) THEN
                s = cS(i)%ptr
                e = s + cS(i)%n - 1
                IF (i .GT. tF) THEN
                    CALL MPI_WAIT(cS(i)%req, stat, ierr)
                    CALL MPI_IRECV(R(s:e), cS(i)%n, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                ELSE
                    CALL MPI_ISEND(rTmp(s:e), cS(i)%n, mpreal, i-1, cS(i)%tag, comm, cS(i)%req, ierr)
                END IF
            END IF
        END DO

        DO i=1, nTasks
            IF (cS(i)%tag .NE. 0) THEN
                CALL MPI_WAIT(cS(i)%req, stat, ierr)
            END IF
        END DO

        RETURN
    END SUBROUTINE COMMUS

    !====================================================================

    ! --------------------------------------------------
    ! memLS_COMMU_CREATE
    !
    ! creation subroutine for memLS_lhsType, dataype for
    ! an object that containing information about
    ! the Left Hand Side (LHS) matrix
    ! (sparse, compressed column storage).
    !
    ! lhs    : Handle to the LHS structure
    ! commu  : Handle to the communicator structure
    ! gnNo   : Total number of nodes (global node No.)
    ! nNo    : Number of nodes on this processor
    ! nnz    : Number of nonzero entries in the LHS
    !          (of calling processor)
    ! ltg    : Local to global index
    ! rowPtr : Pointer to the beginning of each row
    ! colPtr : Column number of LHS sparse array
    ! nFaces : Number of boundary conditions supposed to be added
    !          to this structure
    !
    ! The local to global pointer vector, ltg(a),
    ! should return the node number of the node a in
    ! the entire mesh.
    ! The row pointer, rowPtr(i),
    ! should return the LHS entries position that
    ! corresponds to the first element of row i.
    ! The column pointer, colPtr(j), should return the
    ! column number of jth entire of LHS sparse vector.

    SUBROUTINE memLS_LHS_CREATE(lhs, commu, gnNo, nNo, nnz, gNodes, rowPtr, colPtr, nFaces)

        TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
        TYPE(memLS_commuType), INTENT(IN) :: commu
        INTEGER, INTENT(IN) :: gnNo, nNo, nnz
        INTEGER, INTENT(IN) :: gNodes(nNo), rowPtr(nNo+1), colPtr(nnz)
        INTEGER, INTENT(IN) :: nFaces

        INTEGER i, j, k, a, Ac, ai, s, e, nTasks, tF, maxnNo, ierr, stat(MPI_STATUS_SIZE)

        INTEGER comm
        INTEGER, ALLOCATABLE :: aNodes(:,:), gtlPtr(:), ltg(:), part(:), sCount(:), disp(:)

        IF (lhs%foC) THEN
            PRINT *, "LHS is not free"
            PRINT *, "You may use memLS_LHS_FREE to free this structure"
        END IF

        lhs%foC    = .TRUE.
        lhs%gnNo   = gnNo
        lhs%nNo    = nNo
        lhs%nnz    = nnz
        lhs%commu  = commu
        lhs%nFaces = nFaces

        nTasks = commu%nTasks
        comm   = commu%comm
        tF     = commu%tF

        if ((.not.allocated(lhs%colPtr)).and.(.not.allocated(lhs%rowPtr)).and.(.not.allocated(lhs%diagPtr)).and.(.not.allocated(lhs%map)).and.(.not.allocated(lhs%cS)).and.(.not.allocated(lhs%face))) then
          ALLOCATE (lhs%colPtr(nnz), lhs%rowPtr(2,nNo), lhs%diagPtr(nNo), lhs%map(nNo), lhs%cS(nTasks), lhs%face(nFaces))
        endif

        IF (nTasks .EQ. 1) THEN
            DO i=1, nnz
                lhs%colPtr(i) = colPtr(i)
            END DO
            DO Ac=1, nNo
                s = rowPtr(Ac)
                e = rowPtr(Ac+1) - 1
                DO i=s, e
                    a = colPtr(i)
                    IF (Ac .EQ. a) THEN
                        lhs%diagPtr(Ac) = i
                        EXIT
                    END IF
                END DO

                lhs%rowPtr(1,Ac) = s
                lhs%rowPtr(2,Ac) = e

                lhs%map(Ac) = Ac
            END DO

            lhs%mynNo = nNo
            RETURN
        END IF

        CALL MPI_ALLREDUCE (nNo, maxnNo, 1, mpint, MPI_MAX, comm, ierr)

        if ((.not.allocated(aNodes)).and.(.not.allocated(part)).and.(.not.allocated(sCount)).and.(.not.allocated(disp)).and.(.not.allocated(gtlPtr)).and.(.not.allocated(ltg))) then
          ALLOCATE(aNodes(maxnNo,nTasks), part(maxnNo), sCount(nTasks), disp(nTasks), gtlPtr(gnNo), ltg(nNo))
        endif

        part = 0
        part(1:nNo) = gNodes

        DO i=1, nTasks
            disp(i)   = (i-1)*maxnNo
            sCount(i) = maxnNo
        END DO
        CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp, mpint, comm, ierr)

        gtlPtr = 0
        DO a=1, nNo
            Ac = gNodes(a)
            gtlPtr(Ac) = a
        END DO

        DO i=nTasks, 1, -1
            IF (i .EQ. tF) CYCLE

            DO a=1, maxnNo
                Ac = aNodes(a,i)
                IF (Ac .EQ. 0) EXIT
                ai = gtlPtr(Ac)
                IF (ai .NE. 0) THEN
                    IF (aNodes(ai,tF) .NE. 0) THEN
                        IF (i .GT. tF) aNodes(ai,tF) = 0
                    ELSE
                        aNodes(a,i) = 0
                    END IF
                ELSE
                    aNodes(a,i) = 0
                END IF
            END DO
        END DO

        j = 1
        lhs%cS(1)%ptr = 1
        DO i=1, nTasks
            lhs%cS(i)%n = 0
            lhs%cS(i)%ptr = j
            IF (i.NE.tF .AND. i.NE.1)  THEN
                lhs%cS(i)%ptr = lhs%cS(i-1)%ptr + lhs%cS(i-1)%n
            END IF

            DO a=1, maxnNo
                Ac = aNodes(a,i)
                IF (Ac .NE. 0) THEN
                    lhs%cS(i)%n = lhs%cS(i)%n + 1
                    ai = gtlPtr(Ac)
                    IF (i.GT.tF .OR. aNodes(ai,tF).NE.0) THEN
                        ltg(j) = Ac
                        j = j + 1
                        aNodes(ai,tF) = 0
                    END IF
                END IF
            END DO

            IF (i .LT. tF) THEN
                lhs%cS(i)%tag = nTasks*i  + tF
            ELSE
                lhs%cS(i)%tag = nTasks*tF + i
            END IF
            IF (lhs%cS(i)%n .EQ. 0) lhs%cS(i)%tag = 0
        END DO

        lhs%cS(tF)%tag = 0
        lhs%mynNo = lhs%cS(tF)%ptr + lhs%cS(tF)%n - 1

        gtlPtr = 0
        DO a=1, nNo
            Ac = ltg(a)
            gtlPtr(Ac) = a
        END DO
        DO a=1, nNo
            Ac = gNodes(a)
            lhs%map(a) = gtlPtr(Ac)
        END DO

        DEALLOCATE(aNodes, part, gtlPtr, sCount, disp)

        DO a=1, nNo
            Ac = lhs%map(a)
            lhs%rowPtr(1,Ac) = rowPtr(a)
            lhs%rowPtr(2,Ac) = rowPtr(a+1) - 1
        END DO

        DO i=1, nnz
            lhs%colPtr(i) = lhs%map(colPtr(i))
        END DO

        DO Ac=1, nNo
            DO i=lhs%rowPtr(1,Ac), lhs%rowPtr(2,Ac)
                a = lhs%colPtr(i)
                IF (Ac .EQ. a) THEN
                    lhs%diagPtr(Ac) = i
                    EXIT
                END IF
            END DO
        END DO

        IF (tF .NE. 1) THEN
            i = tF - 1
            i = lhs%cS(i)%ptr + lhs%cS(i)%n - 1
            if (.not.allocated(part)) then
              ALLOCATE(part(i))
            endif
        END IF

        DO i=1, nTasks
            lhs%cS(i)%nBl = 0
            IF (lhs%cS(i)%tag .NE. 0) THEN
                s = lhs%cS(i)%ptr
                e = s + lhs%cS(i)%n - 1
                IF (i .LT. tF) THEN
                    CALL MPI_RECV(part(s:e), lhs%cS(i)%n, mpint, i-1, lhs%cS(i)%tag, comm, stat, ierr)

                    k = 0
                    DO j=s, e
                        k = k + 1
                        IF (part(j).NE.ltg(k) .OR. j.EQ.s) THEN
                            lhs%cS(i)%nBl = lhs%cS(i)%nBl + 1
                            DO k=1, lhs%cS(tF)%ptr
                                IF (part(j) .EQ. ltg(k)) EXIT
                            END DO
                        END IF
                    END DO
                    a = lhs%cS(i)%nBl
                    if ((.not.allocated(lhs%cS(i)%blPtr)).and.(.not.allocated(lhs%cS(i)%blN))) then
                      ALLOCATE(lhs%cS(i)%blPtr(a), lhs%cS(i)%blN(a))
                    endif

                    k = 0
                    a = 0
                    DO j=s, e
                        k = k + 1
                        IF (part(j).NE.ltg(k) .OR. j.EQ.s) THEN
                            a = a + 1
                            lhs%cS(i)%blN(a) = 1
                            DO k=1, lhs%cS(tF)%ptr
                                IF (part(j) .EQ. ltg(k)) THEN
                                    lhs%cS(i)%blPtr(a) = k
                                    EXIT
                                END IF
                            END DO
                        ELSE
                            lhs%cS(i)%blN(a) = lhs%cS(i)%blN(a) + 1
                        END IF
                    END DO
                ELSE
                    CALL MPI_SEND(ltg(s:e), lhs%cS(i)%n, mpint, i-1, lhs%cS(i)%tag, comm, stat, ierr)
                END IF
            END IF
        END DO

        IF (ALLOCATED(part)) DEALLOCATE(part)
        if (allocated(ltg)) then
          deallocate(ltg)
        endif

        RETURN
    END SUBROUTINE memLS_LHS_CREATE

    !====================================================================

    SUBROUTINE memLS_LHS_FREE(lhs)

        TYPE(memLS_lhsType), INTENT(INOUT) :: lhs

        INTEGER faIn, i

        IF (.NOT.lhs%foC) THEN
            PRINT *, 'Cannot free LHS'
            PRINT *, 'It is not created yet'
            STOP
        END IF

        DO faIn = 1, lhs%nFaces
            IF (lhs%face(faIn)%foC) CALL memLS_BC_FREE(lhs, faIn)
        END DO

        DO i=1, lhs%commu%nTasks
            IF (ALLOCATED(lhs%cS(i)%blPtr)) THEN
                DEALLOCATE(lhs%cS(i)%blPtr, lhs%cS(i)%blN)
            END IF
        END DO

        lhs%foC    = .FALSE.
        lhs%gnNo   = 0
        lhs%nNo    = 0
        lhs%nnz    = 0
        lhs%nFaces = 0

        DEALLOCATE (lhs%colPtr, lhs%rowPtr, lhs%diagPtr, lhs%map, lhs%cS, lhs%face)

        RETURN
    END SUBROUTINE memLS_LHS_FREE

    !====================================================================

    SUBROUTINE memLS_LHS_CREATE_C(pLHS, commu, gnNo, nNo, nnz, gNodes, rowPtr, colPtr, nFaces)

        TYPE(memLS_lhsType), POINTER, INTENT(OUT) :: pLHS
        TYPE(memLS_commuType), INTENT(IN) :: commu
        INTEGER, INTENT(IN) :: gnNo, nNo, nnz
        INTEGER, INTENT(IN) :: gNodes(nNo), rowPtr(nNo+1), colPtr(nnz)
        INTEGER, INTENT(IN) :: nFaces

        TYPE(memLS_lhsType), TARGET, SAVE :: lhs

        CALL memLS_LHS_CREATE(lhs, commu, gnNo, nNo, nnz, gNodes, rowPtr, colPtr, nFaces)

        pLHS => lhs

        RETURN
    END SUBROUTINE memLS_LHS_CREATE_C

    !====================================================================

    ! --------------------------------------------------
    ! memLS_LS_CREATE
    !
    ! creation subroutine for memLS_lsType, the
    ! dataype for an object containing information about the type of solver.
    !
    ! (also contains tolerance parameters)
    !
    ! ls       : Handle to LS object
    ! LS_type  : LS algorithm LS TYPE GMRES to be used in memLS
    ! relTol   : Relative tolerance (ratio between final & initial residuals)
    ! absTole  : Absolute tolerance (final residual)
    ! maxItr   : Maximum number of iterations
    ! dimKry   : Krylov space dimension (for GMRES and NS solvers)
    ! relTolIn : Inner loop relative tolerances (for NS solver)
    ! absTolIn : Inner loop absolute tolerances (for NS solver)
    ! maxItrIn : Inner loop maximum number of iterations (for NS solver)

    SUBROUTINE memLS_LS_CREATE(ls, LS_type, relTol, absTol, maxItr, dimKry, relTolIn, absTolIn, maxItrIn)

        TYPE(memLS_lsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: LS_type
        REAL*8, INTENT(IN), OPTIONAL :: relTol, absTol, relTolIn(2), absTolIn(2)
        INTEGER, INTENT(IN), OPTIONAL :: maxItr, dimKry, maxItrIn(2)

        IF (ls%foC) THEN
            PRINT *, "LS is not free"
            PRINT *, "You may use memLS_LS_FREE to free this structure"
        END IF

        ls%foC     = .TRUE.
        ls%LS_type = LS_type

        SELECT CASE (LS_type)
            CASE (LS_TYPE_NS)
                ls%RI%relTol = 4D-1
                ls%GM%relTol = 1D-2
                ls%CG%relTol = 1D-1
                ls%RI%mItr = 10
                ls%GM%mItr = 3
                ls%CG%mItr = 500
                ls%GM%sD   = 50
            CASE (LS_TYPE_GMRES)
                ls%RI%relTol = 1D-2
                ls%RI%mItr   = 2
                ls%RI%sD     = 150
            CASE (LS_TYPE_CG)
                ls%RI%reltol = 1D-4
                ls%RI%mItr   = 1000
            CASE DEFAULT
                PRINT *, 'Solver type LS_TYPE is not defined'
                STOP
        END SELECT
        ls%RI%absTol = 1D-10
        ls%GM%absTol = 1D-10
        ls%CG%absTol = 1D-10

        IF (PRESENT(relTol)) ls%RI%relTol = relTol
        IF (PRESENT(absTol)) ls%RI%absTol = absTol
        IF (PRESENT(maxItr)) ls%RI%mItr   = maxItr

        IF (PRESENT(dimKry)) THEN
            ls%RI%sD = dimKry
            ls%GM%sD = dimKry
        END IF
        IF (PRESENT(relTolIn)) THEN
            ls%GM%relTol = relTolIn(1)
            ls%CG%relTol = relTolIn(2)
        END IF
        IF (PRESENT(absTolIn)) THEN
            ls%GM%absTol = absTolIn(1)
            ls%CG%absTol = absTolIn(2)
        END IF
        IF (PRESENT(maxItrIn)) THEN
            ls%GM%mItr = maxItrIn(1)
            ls%CG%mItr = maxItrIn(2)
        END IF

        RETURN
    END SUBROUTINE memLS_LS_CREATE

    !====================================================================

    SUBROUTINE memLS_LS_FREE (ls)

        TYPE(memLS_lsType), INTENT(INOUT) :: ls

        IF (.NOT.ls%foC) THEN
            PRINT *, 'Cannot free LS'
            PRINT *, 'It is not created yet'
            STOP
        END IF
        ls%foC  = .FALSE.

        RETURN
    END SUBROUTINE memLS_LS_FREE

    !====================================================================

    SUBROUTINE NSSOLVER(nFaces, gnNo, dof, nNo, nnz, mynNo, commu, cS, face, ls, rowPtr, colPtr, Val, Ri)

        INTEGER, INTENT(IN) :: nFaces, gnNo, dof, nNo, nnz, mynNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_faceType), INTENT(INOUT) :: face(nFaces)
        TYPE(memLS_lsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: Val(dof*dof,nnz)
        REAL*8, INTENT(INOUT) :: Ri(dof,nNo)

        INTEGER i, j, k, iB, iBB, nB, nsd
        REAL*8 eps!REAL*8 CPUT, NORMS, NORMV, DOTS, DOTV, eps
        REAL*8, ALLOCATABLE :: U(:,:,:), P(:,:), &
            MU(:,:,:), MP(:,:), A(:,:), B(:), xB(:), mK(:,:), mG(:,:), &
            mD(:,:), mL(:), Gt(:,:), Rm(:,:), Rc(:), Rmi(:,:), Rci(:)

        nsd = dof - 1
        iB = ls%RI%mItr
        nB = 2*iB
        ALLOCATE(Rm(nsd,nNo), Rc(nNo), Rmi(nsd,nNo), Rci(nNo), &
            U(nsd,nNo,iB), P(nNo,iB), MU(nsd,nNo,nB), MP(nNo,nB), &
            A(nB,nB), B(nB), xB(nB))

        Rmi = Ri(1:nsd,:)
        Rci = Ri(dof,:)

        xB          = 0D0
        B           = 0D0
        Rm          = Rmi
        Rc          = Rci
        eps         = SQRT(NORMV(nsd, mynNo, commu, Rm)**2D0 &
            +      NORMS(     mynNo, commu, Rc)**2D0)
        ls%RI%iNorm = eps
        ls%RI%fNorm = eps
        ls%CG%callD = 0D0
        ls%GM%callD = 0D0
        ls%CG%itr   = 0
        ls%GM%itr   = 0
        ls%RI%callD = CPUT()
        ls%RI%suc   = .FALSE.
        eps         = MAX(ls%RI%absTol,ls%RI%relTol*eps)

        CALL DEPART
        CALL BCPRE

        DO i=1, ls%RI%mItr
            iB  = 2*i - 1
            iBB = 2*i
            ls%RI%dB = ls%RI%fNorm
            CALL GMRES(nFaces, nsd, nNo, nnz, mynNo, commu, cS, face, ls%GM, rowPtr, colPtr, mK, Rm, U(:,:,i))

            CALL SPARMULVS(nsd, nNo, nnz, commu, cS, rowPtr, colPtr, mD, U(:,:,i), P(:,i))

            P(:,i) = Rc - P(:,i)
            CALL CGRAD(nFaces, nsd, nNo, nnz, mynNo, commu, cS, face, ls%CG, rowPtr, colPtr, Gt, mG, mL, P(:,i))

            CALL SPARMULSV(nsd, nNo, nnz, commu, cS, rowPtr, colPtr, mG, P(:,i), MU(:,:,iB))

            MU(:,:,iBB) = Rm - MU(:,:,iB)
            CALL GMRES(nFaces, nsd, nNo, nnz, mynNo, commu, cS, face, ls%GM, rowPtr, colPtr, mK, MU(:,:,iBB), U(:,:,i))

            CALL SPARMULVV(nsd, nNo, nnz, commu, cS, rowPtr, colPtr, mK, U(:,:,i), MU(:,:,iBB))

            CALL ADDBCMUL(BCOP_TYPE_ADD, nFaces, nsd, nNo, mynNo, commu, face, U(:,:,i), MU(:,:,iBB))

            CALL SPARMULSS(nNo, nnz, commu, cS, rowPtr, colPtr, mL, P(:,i), MP(:,iB))

            CALL SPARMULVS(nsd, nNo, nnz, commu, cS, rowPtr, colPtr, mD, U(:,:,i), MP(:,iBB))

            DO k=iB, iBB
                DO j=1, k - 1
                    A(j,k) = DOTV(nsd, mynNo, commu, MU(:,:,j), MU(:,:,k)) &
                        + DOTS(     mynNo, commu, MP(:,j),   MP(:,k))
                    A(k,j) = A(j,k)
                END DO
                A(k,k) = NORMV(nsd, mynNo, commu, MU(:,:,k))**2D0 &
                    + NORMS(     mynNo, commu, MP(:,k))**2D0
                B(k)   = DOTV (nsd, mynNo, commu, MU(:,:,k), Rmi) &
                    + DOTS (     mynNo, commu, MP(:,k), Rci)
            END DO

            xB = B
            CALL GE(iBB, A(1:iBB,1:iBB), xB(1:iBB))

            ls%RI%fNorm = SQRT(ls%RI%iNorm**2D0 - SUM(xB(1:iBB)*B(1:iBB)))
            IF(ls%RI%fNorm .LT. eps) THEN
                ls%RI%suc = .TRUE.
                EXIT
            END IF

            Rm = Rmi - xB(1)*MU(:,:,1)
            Rc = Rci - xB(1)*MP(:,1)
            DO j=2, iBB
                Rm = Rm - xB(j)*MU(:,:,j)
                Rc = Rc - xB(j)*MP(:,j)
            END DO
        END DO
        IF (i .GT. ls%RI%mItr) THEN
            ls%RI%itr = ls%RI%mItr
        ELSE
            ls%RI%itr = i

            Rc = Rci - xB(1)*MP(:,1)
            DO j=2, iBB
                Rc = Rc - xB(j)*MP(:,j)
            END DO
        END IF
        ls%Resc = NINT(1D2*(NORMS(mynNo, commu, Rc)/ls%RI%fNorm)**2D0)
        ls%Resm = 100 - ls%Resc

        Rmi = xB(2)*U(:,:,1)
        Rci = xB(1)*P(:,1)
        DO i=2, ls%RI%itr
            iB  = 2*i - 1
            iBB = 2*i

            Rmi = Rmi + xB(iBB)*U(:,:,i)
            Rci = Rci + xB(iB)*P(:,i)
        END DO

        ls%RI%callD = CPUT() - ls%RI%callD
        ls%RI%dB    = 1D1*LOG(ls%RI%fNorm/ls%RI%dB)

        Ri(1:nsd,:) = Rmi
        Ri(dof,:) = Rci

        if ((allocated(Rm)).and. &
            (allocated(Rc)).and. &
            (allocated(Rmi)).and. &
            (allocated(Rci)).and. &
            (allocated(U)).and. &
            (allocated(P)).and. &
            (allocated(MU)).and. &
            (allocated(MP)).and. &
            (allocated(A)).and. &
            (allocated(B)).and. &
            (allocated(mK)).and. &
            (allocated(mD)).and. &
            (allocated(mG)).and. &
            (allocated(mL)).and. &
            (allocated(Gt))) then
              DEALLOCATE (Rm, Rc, Rmi, Rci, U, P, MU, MP, A, B, mK, mD, mG, mL, Gt)
        endif

        IF (commu%masF) CALL LOGFILE

        RETURN
    CONTAINS

        !====================================================================

        SUBROUTINE DEPART

            IMPLICIT NONE

            INTEGER i, j, k, l
            REAL*8 tmp((nsd+1)*(nsd+1))

            ALLOCATE(mK(nsd*nsd,nnz), mG(nsd,nnz), mD(nsd,nnz), mL(nnz), Gt(nsd,nnz))

            IF (nsd .EQ. 2) THEN
                DO i=1, nnz
                    tmp = Val(:,i)

                    mK(1,i) = tmp(1)
                    mK(2,i) = tmp(2)
                    mK(3,i) = tmp(4)
                    mK(4,i) = tmp(5)

                    mG(1,i) = tmp(3)
                    mG(2,i) = tmp(6)

                    mD(1,i) = tmp(7)
                    mD(2,i) = tmp(8)

                    mL(i)   = tmp(9)
                END DO
            ELSE IF(nsd .EQ. 3) THEN
                DO i=1, nnz
                    tmp = Val(:,i)

                    mK(1,i) = tmp(1)
                    mK(2,i) = tmp(2)
                    mK(3,i) = tmp(3)
                    mK(4,i) = tmp(5)
                    mK(5,i) = tmp(6)
                    mK(6,i) = tmp(7)
                    mK(7,i) = tmp(9)
                    mK(8,i) = tmp(10)
                    mK(9,i) = tmp(11)

                    mG(1,i) = tmp(4)
                    mG(2,i) = tmp(8)
                    mG(3,i) = tmp(12)

                    mD(1,i) = tmp(13)
                    mD(2,i) = tmp(14)
                    mD(3,i) = tmp(15)

                    mL(i)   = tmp(16)
                END DO
            ELSE
                PRINT *, "Not defined nsd for DEPART", nsd
            END IF

            DO i=1, nNo
                Do j=rowPtr(1,i), rowPtr(2,i)
                    k = colPtr(j)
                    DO l=rowPtr(1,k), rowPtr(2,k)
                        IF (colPtr(l) .EQ. i) THEN
                            Gt(:,l) = -mG(:,j)
                            EXIT
                        END IF
                    END DO
                END DO
            END DO

            RETURN
        END SUBROUTINE DEPART

        !====================================================================

        SUBROUTINE BCPRE

            IMPLICIT NONE

            INTEGER faIn, i, a, Ac
            !REAL*8 NORMV
            REAL*8, ALLOCATABLE :: v(:,:)

            DO faIn=1, nFaces
                IF (face(faIn)%coupledFlag) THEN
                    IF (face(faIn)%sharedFlag) THEN
                        IF (.NOT.ALLOCATED(v)) ALLOCATE(v(nsd,nNo))
                        v = 0D0
                        DO a=1, face(faIn)%nNo
                            Ac = face(faIn)%glob(a)
                            DO i=1, nsd
                                v(i,Ac) = face(faIn)%valM(i,a)
                            END DO
                        END DO
                        face(faIn)%nS = NORMV(nsd, mynNo, commu, v)**2D0
                    ELSE
                        face(faIn)%nS = 0D0
                        DO a=1, face(faIn)%nNo
                            Ac = face(faIn)%glob(a)
                            DO i=1, nsd
                                face(faIn)%nS = face(faIn)%nS + face(faIn)%valM(i,a)**2D0
                            END DO
                        END DO
                    END IF
                END IF
            END DO

            RETURN
        END SUBROUTINE BCPRE

        !====================================================================

        SUBROUTINE LOGFILE

            IMPLICIT NONE

            LOGICAL flag
            INTEGER fid, i, j
            CHARACTER*16, PARAMETER :: fName = 'memLS_NS.log'

            INQUIRE(FILE=fName, EXIST=flag)

            fid = 11232
            OPEN(fid, FILE=fName, POSITION='APPEND')

            IF (.NOT.flag) THEN
                i = 0
                DO j=1, nFaces
                    IF (face(j)%coupledFlag) i = i + 1
                END DO
                WRITE(fid,*) gnNo, commu%nTasks, i
            END IF

            i = 0
            IF (ls%RI%suc) i = i + 100
            IF (ls%GM%suc) i = i + 10
            IF (ls%CG%suc) i = i + 1

            WRITE(fid,"(I4.3,I3,I4,I5,3I4,3ES9.2E2,3I4)") &
                i, ls%RI%itr, ls%GM%itr, ls%CG%itr, &
                NINT((ls%RI%CallD-ls%GM%CallD-ls%CG%CallD)/ls%RI%CallD*1D2), &
                NINT(ls%GM%callD/ls%RI%CallD*1D2), &
                NINT(ls%CG%callD/ls%RI%CallD*1D2), &
                ls%RI%iNorm, ls%RI%fNorm/ls%RI%iNorm, ls%RI%CallD, &
                ls%Resm, ls%Resc, NINT(ls%RI%dB)

            CLOSE(fid)

            RETURN
        END SUBROUTINE LOGFILE

    END SUBROUTINE NSSOLVER

    !====================================================================

    SUBROUTINE PRECOND(nFaces, dof, nNo, nnz, commu, cS, face, rowPtr, colPtr, diagPtr, Val, R, W)

        INTEGER, INTENT(IN) :: nFaces, dof, nNo, nnz
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        TYPE(memLS_faceType), INTENT(INOUT) :: face(nFaces)
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz), diagPtr(nNo)
        REAL*8, INTENT(INOUT) :: Val(dof*dof,nnz), R(dof,nNo)
        REAL*8, INTENT(OUT) :: W(dof,nNo)

        INTEGER i, j, a, b, d, Ac, faIn

        SELECT CASE (dof)
            CASE (1)
                DO Ac=1, nNo
                    W(1,Ac) = Val(1,diagPtr(Ac))
                END DO
            CASE(2)
                DO Ac=1, nNo
                    d       = diagPtr(Ac)
                    W(1,Ac) = Val(1,d)
                    W(2,Ac) = Val(4,d)
                END DO
            CASE(3)
                DO Ac=1, nNo
                    d       = diagPtr(Ac)
                    W(1,Ac) = Val(1,d)
                    W(2,Ac) = Val(5,d)
                    W(3,Ac) = Val(9,d)
                END DO
            CASE(4)
                DO Ac=1, nNo
                    d       = diagPtr(Ac)
                    W(1,Ac) = Val(1,d)
                    W(2,Ac) = Val(6,d)
                    W(3,Ac) = Val(11,d)
                    W(4,Ac) = Val(16,d)
                END DO
            CASE DEFAULT
                DO Ac=1, nNo
                    d = diagPtr(Ac)
                    DO i=1, dof
                        W(i,Ac) = Val(i*dof-dof+i,d)
                    END DO
                END DO
        END SELECT

        CALL COMMUV(dof, nNo, commu, cS, W)

        DO Ac=1, nNo
            d = diagPtr(Ac)
            DO i=1, dof
                IF (W(i,Ac) .EQ. 0D0) THEN
                    W(i,Ac) = 1D0
                    Val(i*dof-dof+i,d) = 1D0
                END IF
            END DO
        END DO

        W = 1D0/SQRT(ABS(W))
        DO faIn=1, nFaces
            IF (.NOT.face(faIn)%incFlag) CYCLE
            i = MIN(face(faIn)%dof,dof)
            IF (face(faIn)%bGrp .EQ. BC_TYPE_Dir) THEN
                DO a=1, face(faIn)%nNo
                    Ac = face(faIn)%glob(a)
                    W(1:i,Ac) = W(1:i,Ac)*face(faIn)%val(1:i,a)
                END DO
            END IF
        END DO

        SELECT CASE (dof)
            CASE (1)
                DO Ac=1, nNo
                    a          = rowPtr(1,Ac)
                    b          = rowPtr(2,Ac)
                    Val(1,a:b) = Val(1,a:b)*W(1,Ac)
                END DO
            CASE(2)
                DO Ac=1, nNo
                    a            = rowPtr(1,Ac)
                    b            = rowPtr(2,Ac)
                    Val(1:2,a:b) = Val(1:2,a:b)*W(1,Ac)
                    Val(3:4,a:b) = Val(3:4,a:b)*W(2,Ac)
                END DO
            CASE(3)
                DO Ac=1, nNo
                    a            = rowPtr(1,Ac)
                    b            = rowPtr(2,Ac)
                    Val(1:3,a:b) = Val(1:3,a:b)*W(1,Ac)
                    Val(4:6,a:b) = Val(4:6,a:b)*W(2,Ac)
                    Val(7:9,a:b) = Val(7:9,a:b)*W(3,Ac)
                END DO
            CASE(4)
                DO Ac=1, nNo
                    a              = rowPtr(1,Ac)
                    b              = rowPtr(2,Ac)
                    Val(1:4,a:b)   = Val(1:4,a:b)*W(1,Ac)
                    Val(5:8,a:b)   = Val(5:8,a:b)*W(2,Ac)
                    Val(9:12,a:b)  = Val(9:12,a:b)*W(3,Ac)
                    Val(13:16,a:b) = Val(13:16,a:b)*W(4,Ac)
                END DO
            CASE DEFAULT
                DO Ac=1, nNo
                    a = rowPtr(1,Ac)
                    b = rowPtr(2,Ac)
                    DO i=1, dof
                        j = i*dof - dof + 1
                        Val(j:i*dof,a:b) = Val(j:i*dof,a:b)*W(i,Ac)
                    END DO
                END DO
        END SELECT
        R = W*R
        SELECT CASE (dof)
            CASE (1)
                DO Ac=1, nNo
                    DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                        a = colPtr(i)
                        Val(1,i) = Val(1,i)*W(1,a)
                    END DO
                END DO
            CASE (2)
                DO Ac=1, nNo
                    DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                        a = colPtr(i)
                        Val(1:3:2,i) = Val(1:3:2,i)*W(1,a)
                        Val(2:4:2,i) = Val(2:4:2,i)*W(2,a)
                    END DO
                END DO
            CASE (3)
                DO Ac=1, nNo
                    DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                        a = colPtr(i)
                        Val(1:7:3,i) = Val(1:7:3,i)*W(1,a)
                        Val(2:8:3,i) = Val(2:8:3,i)*W(2,a)
                        Val(3:9:3,i) = Val(3:9:3,i)*W(3,a)
                    END DO
                END DO
            CASE (4)
                DO Ac=1, nNo
                    DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                        a = colPtr(i)
                        Val(1:13:4,i) = Val(1:13:4,i)*W(1,a)
                        Val(2:14:4,i) = Val(2:14:4,i)*W(2,a)
                        Val(3:15:4,i) = Val(3:15:4,i)*W(3,a)
                        Val(4:16:4,i) = Val(4:16:4,i)*W(4,a)
                    END DO
                END DO
            CASE DEFAULT
                DO Ac=1, nNo
                    DO i=rowPtr(1,Ac), rowPtr(2,Ac)
                        a = colPtr(i)
                        DO b=1, dof
                            j = dof*(dof-1) + b
                            Val(b:j:dof,i) = Val(b:j:dof,i)*W(b,a)
                        END DO
                    END DO
                END DO
        END SELECT

        DO Ac=1, nNo
            i = diagPtr(Ac)
            DO a=1, dof
                b = (a-1)*dof + a
                IF (Val(b,i) .EQ. 0D0) Val(b,i) = 1D0
            END DO
        END DO

        DO faIn=1, nFaces
            IF (face(faIn)%coupledFlag) THEN
                DO a=1, face(faIn)%nNo
                    Ac = face(faIn)%glob(a)
                    DO i=1, MIN(face(faIn)%dof,dof)
                        face(faIn)%valM(i,a) = face(faIn)%val(i,a)*W(i,Ac)
                    END DO
                END DO
            END IF
        END DO

        RETURN
    END SUBROUTINE PRECOND

    !====================================================================

    SUBROUTINE memLS_SOLVE (lhs, ls, dof, Ri, Val, incL, res)

        TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
        TYPE(memLS_lsType), INTENT(INOUT) :: ls
        INTEGER, INTENT(IN) :: dof
        REAL*8, INTENT(INOUT) :: Ri(dof,lhs%nNo)
        REAL*8, INTENT(INOUT) :: Val(dof*dof,lhs%nnz)
        INTEGER, INTENT(IN), OPTIONAL :: incL(lhs%nFaces)
        REAL*8, INTENT(IN), OPTIONAL :: res(lhs%nFaces)

        LOGICAL flag
        INTEGER faIn, a, nNo, nnz, nFaces
        REAL*8, ALLOCATABLE :: R(:,:), W(:,:)

        nNo    = lhs%nNo
        nnz    = lhs%nnz
        nFaces = lhs%nFaces

        IF (lhs%nFaces .NE. 0) THEN
            lhs%face%incFlag = .TRUE.
            IF (PRESENT(incL)) THEN
                DO faIn=1, lhs%nFaces
                    IF (incL(faIn) .EQ. 0) lhs%face(faIn)%incFlag = .FALSE.
                END DO
            END IF

            flag = ANY(lhs%face%bGrp.EQ.BC_TYPE_Neu)
            IF (.NOT.PRESENT(res) .AND. flag) THEN
                PRINT *, "res is required when there is a Neu surface"
            END IF
            DO faIn=1, lhs%nFaces
                lhs%face(faIn)%coupledFlag = .FALSE.
                IF (.NOT.lhs%face(faIn)%incFlag) CYCLE
                flag = lhs%face(faIn)%bGrp .EQ. BC_TYPE_Neu
                IF (flag .AND. res(faIn).NE.0D0) THEN
                    lhs%face(faIn)%res = res(faIn)
                    lhs%face(faIn)%coupledFlag = .TRUE.
                END IF
            END DO
        END IF

        ALLOCATE(R(dof,nNo), W(dof,nNo))
        DO a=1, nNo
            R(:,lhs%map(a)) = Ri(:,a)
        END DO

        CALL COMMUV(dof, nNo, lhs%commu, lhs%cS, R)
        CALL PRECOND(nFaces, dof, nNo, nnz, lhs%commu, lhs%cS, lhs%face, lhs%rowPtr, lhs%colPtr, lhs%diagPtr, Val, R, W)

        SELECT CASE (ls%LS_type)
            CASE (LS_TYPE_NS)
                CALL NSSOLVER(nFaces, lhs%gnNo, dof, nNo, nnz, lhs%mynNo, lhs%commu, lhs%cS, lhs%face, ls, lhs%rowPtr, lhs%colPtr, Val, R)
            CASE (LS_TYPE_GMRES)
                CALL GMRESV(nFaces, dof, nNo, nnz, lhs%mynNo, lhs%commu, lhs%cS, lhs%face, ls%RI, lhs%rowPtr, lhs%colPtr, Val, R)
            CASE (LS_TYPE_CG)
                IF (dof .EQ. 1) THEN
                    CALL CGRADS(nNo, nnz, lhs%mynNo, lhs%commu, lhs%cS, ls%RI, lhs%rowPtr, lhs%colPtr, Val, R)
                ELSE
                    CALL CGRADV(dof, nNo, nnz, lhs%mynNo, lhs%commu, lhs%cS, ls%RI, lhs%rowPtr, lhs%colPtr, Val, R)
                END IF
            CASE DEFAULT
                PRINT *, 'LS_type not defined'
                STOP
        END SELECT
        R = R*W

        DO a=1, nNo
            Ri(:,a) = R(:,lhs%map(a))
        END DO

        DEALLOCATE(R, W)

        RETURN
    END SUBROUTINE memLS_SOLVE

    !====================================================================

    SUBROUTINE SPARMULVV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, K, U, KU)

        INTEGER, INTENT(IN) :: dof, nNo, nnz
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: K(dof*dof,nnz), U(dof,nNo)
        REAL*8, INTENT(OUT) :: KU(dof,nNo)

        INTEGER i, j, l, col, s, e

        KU = 0D0
        SELECT CASE (dof)
            CASE (1)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        KU(1,i) = KU(1,i) + K(1,j)*U(1,colPtr(j))
                    END DO
                END DO
            CASE(2)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(1,i) = KU(1,i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)
                        KU(2,i) = KU(2,i) + K(3,j)*U(1,col) + K(4,j)*U(2,col)
                    END DO
                END DO
            CASE(3)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(1,i) = KU(1,i) + K(1,j)*U(1,col) + K(2,j)*U(2,col) &
                            + K(3,j)*U(3,col)
                        KU(2,i) = KU(2,i) + K(4,j)*U(1,col) + K(5,j)*U(2,col) &
                            + K(6,j)*U(3,col)
                        KU(3,i) = KU(3,i) + K(7,j)*U(1,col) + K(8,j)*U(2,col) &
                            + K(9,j)*U(3,col)
                    END DO
                END DO
            CASE(4)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(1,i) = KU(1,i) + K(1 ,j)*U(1,col) + K(2 ,j)*U(2,col) &
                            + K(3 ,j)*U(3,col) + K(4 ,j)*U(4,col)
                        KU(2,i) = KU(2,i) + K(5 ,j)*U(1,col) + K(6 ,j)*U(2,col) &
                            + K(7 ,j)*U(3,col) + K(8 ,j)*U(4,col)
                        KU(3,i) = KU(3,i) + K(9 ,j)*U(1,col) + K(10,j)*U(2,col) &
                            + K(11,j)*U(3,col) + K(12,j)*U(4,col)
                        KU(4,i) = KU(4,i) + K(13,j)*U(1,col) + K(14,j)*U(2,col) &
                            + K(15,j)*U(3,col) + K(16,j)*U(4,col)
                    END DO
                END DO
            CASE DEFAULT
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        DO l=1, dof
                            e = l*dof
                            s = e - dof + 1
                            KU(l,i) = KU(l,i) + SUM(K(s:e,j)*U(:,col))
                        END DO
                    END DO
                END DO
        END SELECT

        CALL COMMUV(dof, nNo, commu, cS, KU)

        RETURN
    END SUBROUTINE SPARMULVV

    !====================================================================

    SUBROUTINE SPARMULVS(dof, nNo, nnz, commu, cS, rowPtr, colPtr, K, U, KU)

        INTEGER, INTENT(IN) :: dof, nNo, nnz
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: K(dof,nnz), U(dof,nNo)
        REAL*8, INTENT(OUT) :: KU(nNo)

        INTEGER i, j, col

        KU = 0D0
        SELECT CASE (dof)
            CASE (1)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        KU(i) = KU(i) + K(1,j)*U(1,colPtr(j))
                    END DO
                END DO
            CASE(2)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(i) = KU(i) + K(1,j)*U(1,col) + K(2,j)*U(2,col)
                    END DO
                END DO
            CASE(3)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(i) = KU(i) + K(1,j)*U(1,col) + K(2,j)*U(2,col) &
                            + K(3,j)*U(3,col)
                    END DO
                END DO
            CASE(4)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(i) = KU(i) + K(1,j)*U(1,col) + K(2,j)*U(2,col) &
                            + K(3,j)*U(3,col) + K(4,j)*U(4,col)
                    END DO
                END DO
            CASE DEFAULT
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        KU(i) = KU(i) + SUM(K(:,j)*U(:,colPtr(j)))
                    END DO
                END DO
        END SELECT

        CALL COMMUS(nNo, commu, cS, KU)

        RETURN
    END SUBROUTINE SPARMULVS

    !====================================================================

    SUBROUTINE SPARMULSV(dof, nNo, nnz, commu, cS, rowPtr, colPtr, K, U, KU)

        INTEGER, INTENT(IN) :: dof, nNo, nnz
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: K(dof,nnz), U(nNo)
        REAL*8, INTENT(OUT) :: KU(dof,nNo)

        INTEGER i, j, col

        KU = 0D0
        SELECT CASE (dof)
            CASE (1)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        KU(1,i) = KU(1,i) + K(1,j)*U(colPtr(j))
                    END DO
                END DO
            CASE(2)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(1,i) = KU(1,i) + K(1,j)*U(col)
                        KU(2,i) = KU(2,i) + K(2,j)*U(col)
                    END DO
                END DO
            CASE(3)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(1,i) = KU(1,i) + K(1,j)*U(col)
                        KU(2,i) = KU(2,i) + K(2,j)*U(col)
                        KU(3,i) = KU(3,i) + K(3,j)*U(col)
                    END DO
                END DO
            CASE(4)
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        col = colPtr(j)
                        KU(1,i) = KU(1,i) + K(1,j)*U(col)
                        KU(2,i) = KU(2,i) + K(2,j)*U(col)
                        KU(3,i) = KU(3,i) + K(3,j)*U(col)
                        KU(4,i) = KU(4,i) + K(4,j)*U(col)
                    END DO
                END DO
            CASE DEFAULT
                DO i=1, nNo
                    DO j=rowPtr(1,i), rowPtr(2,i)
                        KU(:,i) = KU(:,i) + K(:,j)*U(colPtr(j))
                    END DO
                END DO
        END SELECT

        CALL COMMUV(dof, nNo, commu, cS, KU)

        RETURN
    END SUBROUTINE SPARMULSV

    !====================================================================

    SUBROUTINE SPARMULSS(nNo, nnz, commu, cS, rowPtr, colPtr, K, U, KU)

        INTEGER, INTENT(IN) :: nNo, nnz
        TYPE(memLS_commuType), INTENT(IN) :: commu
        TYPE(memLS_cSType), INTENT(IN) :: cS(commu%nTasks)
        INTEGER, INTENT(IN) :: rowPtr(2,nNo), colPtr(nnz)
        REAL*8, INTENT(IN) :: K(nnz), U(nNo)
        REAL*8, INTENT(OUT) :: KU(nNo)

        INTEGER i, j

        KU = 0D0
        DO i=1, nNo
            DO j=rowPtr(1,i), rowPtr(2,i)
                KU(i) = KU(i) + K(j)*U(colPtr(j))
            END DO
        END DO

        CALL COMMUS(nNo, commu, cS, KU)

        RETURN
    END SUBROUTINE SPARMULSS

    !====================================================================

    FUNCTION NORMV(dof, nNo, commu, U)

        INTEGER, INTENT(IN) :: dof, nNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        REAL*8, INTENT(IN) :: U(dof,nNo)

        INTEGER i, ierr
        REAL*8 tmp, NORMV

        NORMV = 0D0
        SELECT CASE(dof)
            CASE(1)
                DO i=1, nNo
                    NORMV = NORMV + U(1,i)*U(1,i)
                END DO
            CASE(2)
                DO i=1, nNo
                    NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i)
                END DO
            CASE(3)
                DO i=1, nNo
                    NORMV = NORMV+ U(1,i)*U(1,i) + U(2,i)*U(2,i) + U(3,i)*U(3,i)
                END DO
            CASE(4)
                DO i=1, nNo
                    NORMV = NORMV + U(1,i)*U(1,i) + U(2,i)*U(2,i) &
                        + U(3,i)*U(3,i) + U(4,i)*U(4,i)
                END DO
            CASE DEFAULT
                DO i=1, nNo
                    NORMV = NORMV + SUM(U(:,i)*U(:,i))
                END DO
        END SELECT

        IF (commu%nTasks .NE. 1) THEN
            CALL MPI_ALLREDUCE(NORMV, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)
            NORMV = tmp
        END IF

        NORMV = SQRT(NORMV)

        RETURN
    END FUNCTION NORMV

    !====================================================================

    FUNCTION NORMS(nNo, commu, U)

        INTEGER, INTENT(IN) :: nNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        REAL*8, INTENT(IN) :: U(nNo)

        INTEGER i, ierr
        REAL*8 tmp, NORMS

        NORMS = 0D0
        DO i=1, nNo
            NORMS = NORMS + U(i)*U(i)
        END DO

        IF (commu%nTasks .NE. 1) THEN
            CALL MPI_ALLREDUCE(NORMS, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)
            NORMS = tmp
        END IF

        NORMS = SQRT(NORMS)

        RETURN
    END FUNCTION NORMS

    !====================================================================

    FUNCTION DOTV(dof, nNo, commu, U, V)

        INTEGER, INTENT(IN) :: dof, nNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        REAL*8, INTENT(IN) :: V(dof,nNo), U(dof,nNo)

        INTEGER i, ierr
        REAL*8 tmp, DOTV

        DOTV = 0D0
        SELECT CASE(dof)
            CASE(1)
                DO i=1, nNo
                    DOTV = DOTV + U(1,i)*V(1,i)
                END DO
            CASE(2)
                DO i=1, nNo
                    DOTV = DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i)
                END DO
            CASE(3)
                DO i=1, nNo
                    DOTV = DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i) +U(3,i)*V(3,i)
                END DO
            CASE(4)
                DO i=1, nNo
                    DOTV = DOTV + U(1,i)*V(1,i) + U(2,i)*V(2,i) &
                        + U(3,i)*V(3,i) + U(4,i)*V(4,i)
                END DO
            CASE DEFAULT
                DO i=1, nNo
                    DOTV = DOTV + SUM(U(:,i)*V(:,i))
                END DO
        END SELECT

        IF (commu%nTasks .EQ. 1) RETURN

        CALL MPI_ALLREDUCE(DOTV, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)

        DOTV = tmp

        RETURN
    END FUNCTION DOTV

    !====================================================================

    FUNCTION DOTS(nNo, commu, U, V)

        INTEGER, INTENT(IN) :: nNo
        TYPE(memLS_commuType), INTENT(IN) :: commu
        REAL*8, INTENT(IN) :: V(nNo), U(nNo)

        INTEGER i, ierr
        REAL*8 tmp, DOTS

        DOTS = 0D0
        DO i=1, nNo
            DOTS = DOTS + U(i)*V(i)
        END DO

        IF (commu%nTasks .EQ. 1) RETURN

        CALL MPI_ALLREDUCE(DOTS, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)

        DOTS = tmp

        RETURN
    END FUNCTION DOTS

    !====================================================================

    FUNCTION CPUT()

        IMPLICIT NONE

        INTEGER timeArray(8), i
        INTEGER, PARAMETER::nD(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)

        REAL*8 CPUT

        CALL DATE_AND_TIME (VALUES=timeArray)

        timeArray(3) = timeArray(3) + (timeArray(1) - 2010)*365
        DO i=1, timeArray(2) - 1
            timeArray(3) = timeArray(3) + nD(i)
        END DO
        CPUT = timeArray(3)*8.64D4 + timeArray(5)*3.6D3 + timeArray(6)*6D1 + timeArray(7)*1D0 + timeArray(8)*1D-3

        RETURN
    END FUNCTION CPUT

end module
