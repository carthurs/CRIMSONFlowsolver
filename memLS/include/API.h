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

      INTERFACE
         SUBROUTINE memLS_COMMU_CREATE(commu, commi)
            INCLUDE "STD.h"
            TYPE(memLS_commuType), INTENT(INOUT) :: commu
            INTEGER, INTENT(IN) :: commi
         END SUBROUTINE memLS_COMMU_CREATE
         SUBROUTINE memLS_COMMU_FREE(commu)
            INCLUDE "STD.h"
            TYPE(memLS_commuType), INTENT(INOUT) :: commu
         END SUBROUTINE memLS_COMMU_FREE

         SUBROUTINE memLS_LHS_CREATE(lhs, commu, gnNo, nNo, nnz, &
            gNodes, rowPtr, colPtr, nFaces)
            INCLUDE "STD.h"
            TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
            TYPE(memLS_commuType), INTENT(IN) :: commu
            INTEGER, INTENT(IN) :: gnNo, nNo, nnz
            INTEGER, INTENT(IN) :: gNodes(nNo), rowPtr(nNo+1), &
               colPtr(nnz)
            INTEGER, INTENT(IN) :: nFaces
         END SUBROUTINE memLS_LHS_CREATE
         SUBROUTINE memLS_LHS_FREE(lhs)
            INCLUDE "STD.h"
            TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
         END SUBROUTINE memLS_LHS_FREE

         SUBROUTINE memLS_LS_CREATE(ls, LS_type, relTol, absTol, maxItr, &
            dimKry, relTolIn, absTolIn, maxItrIn)
            INCLUDE "STD.h"
            TYPE(memLS_lsType), INTENT(INOUT) :: ls
            INTEGER, INTENT(IN) :: LS_type
            REAL*8, INTENT(IN), OPTIONAL :: relTol, absTol, relTolIn(2), &
               absTolIn(2)
            INTEGER, INTENT(IN), OPTIONAL :: maxItr, dimKry, maxItrIn(2)
         END SUBROUTINE memLS_LS_CREATE
         SUBROUTINE memLS_LS_FREE(ls)
            INCLUDE "STD.h"
            TYPE(memLS_lsType), INTENT(INOUT) :: ls
         END SUBROUTINE memLS_LS_FREE

         SUBROUTINE memLS_BC_CREATE(lhs, faIn, nNo, dof, BC_type, &
               gNodes, Val)
            INCLUDE "STD.h"
            TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: faIn, nNo, dof
            INTEGER, INTENT(IN) :: BC_type
            INTEGER, INTENT(IN) :: gNodes(nNo)
            REAL*8, INTENT(IN), OPTIONAL :: Val(dof,nNo)
         END SUBROUTINE memLS_BC_CREATE
         SUBROUTINE memLS_BC_FREE(lhs, faIn)
            INCLUDE "STD.h"
            TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
            INTEGER, INTENT(IN) :: faIn
         END SUBROUTINE memLS_BC_FREE

         SUBROUTINE memLS_SOLVE (lhs, ls, dof, Ri, Val, incL, res)
            INCLUDE "STD.h"
            TYPE(memLS_lhsType), INTENT(INOUT) :: lhs
            TYPE(memLS_lsType), INTENT(INOUT) :: ls
            INTEGER, INTENT(IN) :: dof
            REAL*8, INTENT(INOUT) :: Ri(dof,lhs%nNo)
            REAL*8, INTENT(IN) :: Val(dof*dof,lhs%nnz)
            INTEGER, INTENT(IN), OPTIONAL :: incL(lhs%nFaces)
            REAL*8, INTENT(IN), OPTIONAL :: res(lhs%nFaces)
         END SUBROUTINE memLS_SOLVE
      END INTERFACE





