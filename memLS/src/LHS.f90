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

      INCLUDE "STD.h"

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

      ALLOCATE (lhs%colPtr(nnz), lhs%rowPtr(2,nNo), lhs%diagPtr(nNo), lhs%map(nNo), lhs%cS(nTasks), lhs%face(nFaces))

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

      ALLOCATE(aNodes(maxnNo,nTasks), part(maxnNo), sCount(nTasks), disp(nTasks), gtlPtr(gnNo), ltg(nNo))

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
         ALLOCATE(part(i))
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
               ALLOCATE(lhs%cS(i)%blPtr(a), lhs%cS(i)%blN(a))

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
      DEALLOCATE(ltg)

      RETURN
      END SUBROUTINE memLS_LHS_CREATE

!====================================================================

      SUBROUTINE memLS_LHS_FREE(lhs)

      INCLUDE "STD.h"

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

      INCLUDE "STD.h"

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
