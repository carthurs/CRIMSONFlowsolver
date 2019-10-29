// Copyright (C) 2011 Marc Duruflé
// Copyright (C) 2010-2011 Vivien Mallet
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_COMPUTATION_SPARSESOLVER_CXX
#define SELDON_FILE_COMPUTATION_SPARSESOLVER_CXX


#include "Ordering.cxx"
#include "SparseSolver.hxx"


#ifdef SELDON_WITH_UMFPACK
#include "camd.h"
#include "colamd.h"
#endif

namespace Seldon
{


  /*************************
   * Default Seldon solver *
   *************************/


  template<class T, class Allocator>
  SparseSeldonSolver<T, Allocator>::SparseSeldonSolver()
  {
    print_level = -1;
    symmetric_matrix = false;
    permtol = 0.1;
  }


  template<class T, class Allocator>
  void SparseSeldonSolver<T, Allocator>::Clear()
  {
    mat_sym.Clear();
    mat_unsym.Clear();
  }


  template<class T, class Allocator>
  void SparseSeldonSolver<T, Allocator>::HideMessages()
  {
    print_level = -1;
  }


  template<class T, class Allocator>
  void SparseSeldonSolver<T, Allocator>::ShowMessages()
  {
    print_level = 1;
  }


  template<class T, class Allocator>
  double SparseSeldonSolver<T, Allocator>::GetPivotThreshold() const
  {
    return permtol;
  }


  template<class T, class Allocator>
  void SparseSeldonSolver<T, Allocator>::SetPivotThreshold(const double& a)
  {
    permtol = a;
  }


  template<class T, class Allocator>
  template<class T0, class Storage0, class Allocator0>
  void SparseSeldonSolver<T, Allocator>::
  FactorizeMatrix(const IVect& perm,
		  Matrix<T0, General, Storage0, Allocator0>& mat,
		  bool keep_matrix)
  {
    IVect inv_permutation;

    // We convert matrix to unsymmetric format.
    Copy(mat, mat_unsym);

    // Old matrix is erased if needed.
    if (!keep_matrix)
      mat.Clear();

    // We keep permutation array in memory, and check it.
    int n = mat_unsym.GetM();
    if (perm.GetM() != n)
      throw WrongArgument("FactorizeMatrix(IVect&, Matrix&, bool)",
                          "Numbering array is of size "
                          + to_str(perm.GetM())
                          + " while the matrix is of size "
                          + to_str(mat.GetM()) + " x "
                          + to_str(mat.GetN()) + ".");

    permutation_row.Reallocate(n);
    permutation_col.Reallocate(n);
    inv_permutation.Reallocate(n);
    inv_permutation.Fill(-1);
    for (int i = 0; i < n; i++)
      {
        permutation_row(i) = i;
        permutation_col(i) = i;
        inv_permutation(perm(i)) = i;
      }

    for (int i = 0; i < n; i++)
      if (inv_permutation(i) == -1)
        throw WrongArgument("FactorizeMatrix(IVect&, Matrix&, bool)",
                            "The numbering array is invalid.");

    IVect iperm = inv_permutation;

    // Rows of matrix are permuted.
    ApplyInversePermutation(mat_unsym, perm, perm);

    // Temporary vector used for solving.
    xtmp.Reallocate(n);

    // Factorization is performed.
    // Columns are permuted during the factorization.
    symmetric_matrix = false;
    inv_permutation.Fill();
    GetLU(mat_unsym, permutation_col, inv_permutation, permtol, print_level);

    // Combining permutations.
    IVect itmp = permutation_col;
    for (int i = 0; i < n; i++)
      permutation_col(i) = iperm(itmp(i));

    permutation_row = perm;

  }


  template<class T, class Allocator>
  template<class T0, class Storage0, class Allocator0>
  void SparseSeldonSolver<T, Allocator>::
  FactorizeMatrix(const IVect& perm,
		  Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
		  bool keep_matrix)
  {
    IVect inv_permutation;

    // We convert matrix to symmetric format.
    Copy(mat, mat_sym);

    // Old matrix is erased if needed.
    if (!keep_matrix)
      mat.Clear();

    // We keep permutation array in memory, and check it.
    int n = mat_sym.GetM();
    if (perm.GetM() != n)
      throw WrongArgument("FactorizeMatrix(IVect&, Matrix&, bool)",
                          "Numbering array is of size "
                          + to_str(perm.GetM())
                          + " while the matrix is of size "
                          + to_str(mat.GetM()) + " x "
                          + to_str(mat.GetN()) + ".");

    permutation_row.Reallocate(n);
    inv_permutation.Reallocate(n);
    inv_permutation.Fill(-1);
    for (int i = 0; i < n; i++)
      {
        permutation_row(i) = perm(i);
        inv_permutation(perm(i)) = i;
      }

    for (int i = 0; i < n; i++)
      if (inv_permutation(i) == -1)
        throw WrongArgument("FactorizeMatrix(IVect&, Matrix&, bool)",
                            "The numbering array is invalid.");

    // Matrix is permuted.
    ApplyInversePermutation(mat_sym, perm, perm);

    // Temporary vector used for solving.
    xtmp.Reallocate(n);

    // Factorization is performed.
    symmetric_matrix = true;
    GetLU(mat_sym, print_level);
  }


  template<class T, class Allocator> template<class Vector1>
  void SparseSeldonSolver<T, Allocator>::Solve(Vector1& z)
  {
    if (symmetric_matrix)
      {
	for (int i = 0; i < z.GetM(); i++)
	  xtmp(permutation_row(i)) = z(i);

	SolveLU(mat_sym, xtmp);

	for (int i = 0; i < z.GetM(); i++)
	  z(i) = xtmp(permutation_row(i));
      }
    else
      {
	for (int i = 0; i < z.GetM(); i++)
	  xtmp(permutation_row(i)) = z(i);

	SolveLU(mat_unsym, xtmp);

	for (int i = 0; i < z.GetM(); i++)
	  z(permutation_col(i)) = xtmp(i);
      }
  }


  template<class T, class Allocator>
  template<class TransStatus, class Vector1>
  void SparseSeldonSolver<T, Allocator>
  ::Solve(const TransStatus& TransA, Vector1& z)
  {
    if (symmetric_matrix)
      Solve(z);
    else
      if (TransA.Trans())
        {
          for (int i = 0; i < z.GetM(); i++)
            xtmp(i) = z(permutation_col(i));

          SolveLU(SeldonTrans, mat_unsym, xtmp);

          for (int i = 0; i < z.GetM(); i++)
            z(i) = xtmp(permutation_row(i));
        }
      else
        Solve(z);
  }


  /************************************************
   * GetLU and SolveLU used by SeldonSparseSolver *
   ************************************************/


  //! LU factorization with pivot of unsymmetric matrix
  template<class T, class Allocator>
  void GetLU(Matrix<T, General, ArrayRowSparse, Allocator>& A,
	     IVect& iperm, IVect& rperm,
	     double permtol, int print_level)
  {
    int n = A.GetN();

    T fact, s, t;
    double tnorm, zero = 0.0;
    int length_lower, length_upper, jpos, jrow, i_row, j_col;
    int i, j, k, length, size_row, index_lu;

    T czero, cone;
    SetComplexZero(czero);
    SetComplexOne(cone);
    Vector<T, VectFull, Allocator> Row_Val(n);
    IVect Index(n), Row_Ind(n), Index_Diag(n);
    Row_Val.Fill(czero);
    Row_Ind.Fill(-1);
    Index_Diag.Fill(-1);

    Index.Fill(-1);
    // main loop
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row) / double(n-1) * 78.);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#";
                cout.flush();
              }
            old_percent = new_percent;
          }

	size_row = A.GetRowSize(i_row);
	tnorm = zero;

	// tnorm is the sum of absolute value of coefficients of row i_row.
	for (k = 0 ; k < size_row; k++)
          if (A.Value(i_row, k) != czero)
            tnorm += abs(A.Value(i_row, k));

	if (tnorm == zero)
          throw WrongArgument("GetLU(Matrix<ArrayRowSparse>&, IVect&, "
                              "IVect&, double, int)",
                              "The matrix is structurally singular. "
                              "The norm of row " + to_str(i_row)
                              + " is equal to 0.");

        // Unpack L-part and U-part of row of A.
	length_upper = 1;
	length_lower = 0;
	Row_Ind(i_row) = i_row;
	Row_Val(i_row) = czero;
	Index(i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
            k = rperm(A.Index(i_row, j));

	    t = A.Value(i_row,j);
	    if (k < i_row)
	      {
		Row_Ind(length_lower) = k;
		Row_Val(length_lower) = t;
		Index(k) = length_lower;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		Row_Val(i_row) = t;
	      }
	    else
	      {
		jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
	      }
          }

	j_col = 0;
	length = 0;

	// Eliminates previous rows.
	while (j_col < length_lower)
	  {
	    // In order to do the elimination in the correct order, we must
            // select the smallest column index.
	    jrow = Row_Ind(j_col);
	    k = j_col;

	    // Determine smallest column index.
	    for (j = j_col + 1; j < length_lower; j++)
              if (Row_Ind(j) < jrow)
                {
                  jrow = Row_Ind(j);
                  k = j;
                }

	    if (k != j_col)
	      {
		// Exchanging column coefficients.
		j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;

		Index(jrow) = j_col;
		Index(j) = k;

		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }

	    // Zero out element in row.
	    Index(jrow) = -1;

            // Gets the multiplier for row to be eliminated (jrow).
	    // first_index_upper points now on the diagonal coefficient.
	    fact = Row_Val(j_col) * A.Value(jrow, Index_Diag(jrow));

	    // Combines current row and row jrow.
	    for (k = (Index_Diag(jrow)+1); k < A.GetRowSize(jrow); k++)
	      {
		s = fact * A.Value(jrow,k);
		j = rperm(A.Index(jrow,k));

		jpos = Index(j);

		if (j >= i_row)
                  // Dealing with upper part.
                  if (jpos == -1)
                    {
                      // This is a fill-in element.
                      i = i_row + length_upper;
                      Row_Ind(i) = j;
                      Index(j) = i;
                      Row_Val(i) = -s;
                      ++length_upper;
                    }
                  else
                    // This is not a fill-in element.
                    Row_Val(jpos) -= s;
		else
                  // Dealing  with lower part.
                  if (jpos == -1)
                    {
                      // this is a fill-in element
                      Row_Ind(length_lower) = j;
                      Index(j) = length_lower;
                      Row_Val(length_lower) = -s;
                      ++length_lower;
                    }
                  else
                    // This is not a fill-in element.
                    Row_Val(jpos) -= s;
	      }

	    // Stores this pivot element from left to right -- no danger
	    // of overlap with the working elements in L (pivots).
	    Row_Val(length) = fact;
	    Row_Ind(length) = jrow;
	    ++length;
	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row + k)) = -1;

	size_row = length;
	A.ReallocateRow(i_row,size_row);

        // store L-part
	index_lu = 0;
        for (k = 0 ; k < length ; k++)
          {
            A.Value(i_row,index_lu) = Row_Val(k);
            A.Index(i_row,index_lu) = iperm(Row_Ind(k));
            ++index_lu;
          }

	// Saves pointer to beginning of row i_row of U.
	Index_Diag(i_row) = index_lu;

        // Updates. U-matrix -- first apply dropping strategy.
	length = 0;
	for (k = 1; k <= (length_upper-1); k++)
	  {
	    ++length;
	    Row_Val(i_row+length) = Row_Val(i_row+k);
	    Row_Ind(i_row+length) = Row_Ind(i_row+k);
	  }

	length++;

	// Determines next pivot.
        int imax = i_row;
        double xmax = abs(Row_Val(imax));
        double xmax0 = xmax;
        for (k = i_row + 1; k <= i_row + length - 1; k++)
          {
            tnorm = abs(Row_Val(k));
            if (tnorm > xmax && tnorm * permtol > xmax0)
              {
                imax = k;
                xmax = tnorm;
              }
          }

        // Exchanges Row_Val.
        s = Row_Val(i_row);
        Row_Val(i_row) = Row_Val(imax);
        Row_Val(imax) = s;

        // Updates iperm and reverses iperm.
        j = Row_Ind(imax);
        i = iperm(i_row);
        iperm(i_row) = iperm(j);
        iperm(j) = i;

        // Reverses iperm.
        rperm(iperm(i_row)) = i_row;
        rperm(iperm(j)) = j;

	// Copies U-part in original coordinates.
        int index_diag = index_lu;
	A.ResizeRow(i_row, size_row+length);

        for (k = i_row ; k <= i_row + length - 1; k++)
          {
            A.Index(i_row,index_lu) = iperm(Row_Ind(k));
            A.Value(i_row,index_lu) = Row_Val(k);
            ++index_lu;
          }

        // Stores inverse of diagonal element of u.
	A.Value(i_row, index_diag) = cone / Row_Val(i_row);

      } // end main loop

    if (print_level > 0)
      cout << endl;

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize() * (sizeof(T) + 4)) / (1024. * 1024.))
           << " MB" << endl;

    for (i = 0; i < n; i++ )
      for (j = 0; j < A.GetRowSize(i); j++)
        A.Index(i,j) = rperm(A.Index(i,j));
  }


  template<class cplx,
	   class Allocator, class Storage2, class Allocator2>
  void SolveLU(const Matrix<cplx, General, ArrayRowSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x)
  {
    SolveLU(SeldonNoTrans, A, x);
  }


  //! Resolution of LU y = x (x is overwritten with y)
  /*! L and U are assumed to be stored in A. The diagonal of A contains the
    inverse of the diagonal of U.
  */
  template<class cplx, class TransStatus,
	   class Allocator, class Storage2, class Allocator2>
  void SolveLU(const TransStatus& transA,
               const Matrix<cplx, General, ArrayRowSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x)
  {
    int i, k, n, k_;
    cplx inv_diag;
    n = A.GetM();

    if (transA.Trans())
      {
	// Forward solve (with U^T).
	for (i = 0 ; i < n ; i++)
	  {
	    k_ = 0; k = A.Index(i, k_);
	    while ( k < i)
	      {
		k_++;
		k = A.Index(i, k_);
	      }

	    x(i) *= A.Value(i, k_);
	    for (k = k_ + 1; k < A.GetRowSize(i); k++)
		x(A.Index(i, k)) -= A.Value(i, k) * x(i);
	  }

	// Backward solve (with L^T).
	for (i = n-1; i >= 0; i--)
	  {
	    k_ = 0; k = A.Index(i, k_);
	    while (k < i)
	      {
		x(k) -= A.Value(i, k_) * x(i);
		k_++;
		k = A.Index(i, k_);
	      }
	  }
      }
    else
      {
	// Forward solve.
	for (i = 0; i < n; i++)
	  {
	    k_ = 0;
            k = A.Index(i, k_);
	    while (k < i)
	      {
		x(i) -= A.Value(i, k_) * x(k);
		k_++;
		k = A.Index(i, k_);
	      }
	  }

	// Backward solve.
	for (i = n-1; i >= 0; i--)
	  {
	    k_ = 0;
            k = A.Index(i, k_);
	    while (k < i)
	      {
		k_++;
		k = A.Index(i, k_);
	      }

	    inv_diag = A.Value(i, k_);
	    for (k = k_ + 1; k < A.GetRowSize(i); k++)
	      x(i) -= A.Value(i, k) * x(A.Index(i, k));

	    x(i) *= inv_diag;
	  }
      }
  }


  //! LDLt factorization without pivot for symmetric matrix.
  template<class T, class Allocator>
  void GetLU(Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A,
             int print_level)
  {
    int size_row;
    int n = A.GetN();
    double zero = 0.0;

    T fact, s, t;
    double tnorm;
    int length_lower, length_upper, jpos, jrow;
    int i_row, j_col, index_lu, length;
    int i, j, k;

    Vector<T, VectFull, Allocator> Row_Val(n);
    IVect Index(n), Row_Ind(n);
    Row_Val.Zero();
    Row_Ind.Fill(-1);

    Index.Fill(-1);

    // We convert A into an unsymmetric matrix.
    Matrix<T, General, ArrayRowSparse, Allocator> B;
    Seldon::Copy(A, B);

    A.Clear();
    A.Reallocate(n, n);

    // Main loop.
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row) / double(n-1) * 78.);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#";
                cout.flush();
              }
            old_percent = new_percent;
          }

        // 1-norm of the row of initial matrix.
	size_row = B.GetRowSize(i_row);
	tnorm = zero;
	for (k = 0 ; k < size_row; k++)
          tnorm += abs(B.Value(i_row, k));

	if (tnorm == zero)
          throw WrongArgument("GetLU(Matrix<ArrayRowSymSparse>&, int)",
                              "The matrix is structurally singular. "
                              "The norm of row " + to_str(i_row)
                              + " is equal to 0.");

        // Separating lower part from upper part for this row.
	length_upper = 1;
	length_lower = 0;
	Row_Ind(i_row) = i_row;
	Row_Val(i_row) = 0.0;
	Index(i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
	    k = B.Index(i_row,j);
            t = B.Value(i_row,j);
	    if (k < i_row)
	      {
		Row_Ind(length_lower) = k;
		Row_Val(length_lower) = t;
		Index(k) = length_lower;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		Row_Val(i_row) = t;
	      }
	    else
	      {
		jpos = i_row + length_upper;
		Row_Ind(jpos) = k;
		Row_Val(jpos) = t;
		Index(k) = jpos;
		length_upper++;
	      }
          }

        // This row of B is cleared.
        B.ClearRow(i_row);

	j_col = 0;
	length = 0;

        // We eliminate previous rows.
	while (j_col <length_lower)
	  {
	    // In order to do the elimination in the correct order, we must
            // select the smallest column index.
	    jrow = Row_Ind(j_col);
	    k = j_col;

	    // We determine smallest column index.
	    for (j = (j_col+1) ; j < length_lower; j++)
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
	      }

            // If needed, we exchange positions of this element in
            // Row_Ind/Row_Val so that it appears first.
	    if (k != j_col)
	      {

		j = Row_Ind(j_col);
		Row_Ind(j_col) = Row_Ind(k);
		Row_Ind(k) = j;

		Index(jrow) = j_col;
		Index(j) = k;

		s = Row_Val(j_col);
		Row_Val(j_col) = Row_Val(k);
	        Row_Val(k) = s;
	      }

            // Zero out element in row by setting Index to -1.
	    Index(jrow) = -1;

	    // Gets the multiplier for row to be eliminated.
	    fact = Row_Val(j_col) * A.Value(jrow, 0);

	    // Combines current row and row jrow.
	    for (k = 1; k < A.GetRowSize(jrow); k++)
	      {
		s = fact * A.Value(jrow, k);
		j = A.Index(jrow, k);

		jpos = Index(j);
		if (j >= i_row)
		  {

		    // Dealing with upper part.
		    if (jpos == -1)
		      {

			// This is a fill-in element.
			i = i_row + length_upper;
			Row_Ind(i) = j;
			Index(j) = i;
			Row_Val(i) = -s;
			++length_upper;
		      }
		    else
		      {
			// This is not a fill-in element.
			Row_Val(jpos) -= s;
		      }
		  }
		else
		  {
		    // Dealing  with lower part.
		    if (jpos == -1)
		      {
			// This is a fill-in element.
			Row_Ind(length_lower) = j;
			Index(j) = length_lower;
			Row_Val(length_lower) = -s;
			++length_lower;
		      }
		    else
		      {
			// This is not a fill-in element.
			Row_Val(jpos) -= s;
		      }
		  }
	      }

	    // We store this pivot element (from left to right -- no
	    // danger of overlap with the working elements in L (pivots).
	    Row_Val(length) = fact;
	    Row_Ind(length) = jrow;
	    ++length;
	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row + k)) = -1;

	// Updating U-matrix
	length = 0;
	for (k = 1; k <= length_upper - 1; k++)
	  {
	    ++length;
	    Row_Val(i_row + length) = Row_Val(i_row + k);
	    Row_Ind(i_row + length) = Row_Ind(i_row + k);
	  }

	length++;

	// Copies U-part in matrix A.
	A.ReallocateRow(i_row, length);
	index_lu = 1;
	for (k = i_row + 1 ; k <= i_row + length - 1 ; k++)
	  {
	    A.Index(i_row, index_lu) = Row_Ind(k);
	    A.Value(i_row, index_lu) = Row_Val(k);
	    ++index_lu;
	  }

	// Stores the inverse of the diagonal element of u.
	A.Value(i_row,0) = 1.0 / Row_Val(i_row);

      } // end main loop.

    if (print_level > 0)
      cout << endl;

    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i, j) *= A.Value(i, 0);

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize() * (sizeof(T) + 4)) / (1024. * 1024.))
           << " MB" << endl;

  }


  //! Resolution of L D L^t y = x (result is overwritten in x)
  /*! The factor L^t is assumed to be stored in matrix A. The diagonal of A is
    equal to the inverse of diagonal D.
  */
  template<class real, class cplx, class Allocator,
           class Storage2, class Allocator2>
  void SolveLU(const Matrix<real, Symmetric, ArrayRowSymSparse, Allocator>& A,
               Vector<cplx, Storage2, Allocator2>& x)
  {
    int n = A.GetM(), j_row;
    cplx tmp;

    // We solve first L y = b.
    for (int i_col = 0; i_col < n ; i_col++)
      for (int k = 1; k < A.GetRowSize(i_col) ; k++)
        {
          j_row = A.Index(i_col, k);
          x(j_row) -= A.Value(i_col, k)*x(i_col);
        }

    // Inverting by diagonal D.
    for (int i_col = 0; i_col < n ; i_col++)
      x(i_col) *= A.Value(i_col, 0);

    // Then we solve L^t x = y.
    for (int i_col = n-1; i_col >=0; i_col--)
      {
	tmp = x(i_col);
	for (int k = 1; k < A.GetRowSize(i_col); k++)
	  {
	    j_row = A.Index(i_col,k);
	    tmp -= A.Value(i_col,k) * x(j_row);
	  }
	x(i_col) = tmp;
      }
  }


  /********************************************
   * GetLU and SolveLU for SeldonSparseSolver *
   ********************************************/


  template<class T, class Storage, class Allocator, class Alloc2>
  void GetLU(Matrix<T, Symmetric, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     bool keep_matrix = false)
  {
    // identity ordering
    IVect perm(A.GetM());
    perm.Fill();

    // factorization
    mat_lu.FactorizeMatrix(perm, A, keep_matrix);
  }


  template<class T, class Storage, class Allocator, class Alloc2>
  void GetLU(Matrix<T, General, Storage, Allocator>& A,
	     SparseSeldonSolver<T, Alloc2>& mat_lu,
	     bool keep_matrix = false)
  {
    // identity ordering
    IVect perm(A.GetM());
    perm.Fill();

    // factorization
    mat_lu.FactorizeMatrix(perm, A, keep_matrix);
  }


  template<class T, class Alloc2, class Allocator>
  void SolveLU(SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  template<class T, class Alloc2, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       SparseSeldonSolver<T, Alloc2>& mat_lu,
	       Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  /*************************
   * Choice of best solver *
   *************************/


  //! Default constructor
  template<class T>
  SparseDirectSolver<T>::SparseDirectSolver()
  {
    n = 0;
    type_ordering = SparseMatrixOrdering::AUTO;

    type_solver = SELDON_SOLVER;

    // We try to use an other available solver.
    // The order of preference is Pastix, Mumps, UmfPack and SuperLU.
#ifdef SELDON_WITH_SUPERLU
    type_solver = SUPERLU;
#endif
#ifdef SELDON_WITH_UMFPACK
    type_solver = UMFPACK;
#endif
#ifdef SELDON_WITH_MUMPS
    type_solver = MUMPS;
#endif
#ifdef SELDON_WITH_PASTIX
    type_solver = PASTIX;
#endif

    number_threads_per_node = 1;
    threshold_matrix = 0;
    enforce_unsym_ilut = false;
  }


  //! Hiding all messages.
  template<class T>
  void SparseDirectSolver<T>::HideMessages()
  {
    mat_seldon.HideMessages();

#ifdef SELDON_WITH_MUMPS
    mat_mumps.HideMessages();
#endif

#ifdef SELDON_WITH_SUPERLU
    mat_superlu.HideMessages();
#endif

#ifdef SELDON_WITH_UMFPACK
    mat_umf.HideMessages();
#endif

#ifdef SELDON_WITH_PASTIX
    mat_pastix.HideMessages();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
    mat_ilut.SetPrintLevel(0);
#endif
  }


  //! Displaying basic messages.
  template<class T>
  void SparseDirectSolver<T>::ShowMessages()
  {
    mat_seldon.ShowMessages();

#ifdef SELDON_WITH_MUMPS
    mat_mumps.ShowMessages();
#endif

#ifdef SELDON_WITH_SUPERLU
    mat_superlu.ShowMessages();
#endif

#ifdef SELDON_WITH_UMFPACK
    mat_umf.ShowMessages();
#endif

#ifdef SELDON_WITH_PASTIX
    mat_pastix.ShowMessages();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
    mat_ilut.SetPrintLevel(1);
#endif
  }


  //! Displaying all the messages.
  template<class T>
  void SparseDirectSolver<T>::ShowFullHistory()
  {
    mat_seldon.ShowMessages();

#ifdef SELDON_WITH_MUMPS
    mat_mumps.ShowMessages();
#endif

#ifdef SELDON_WITH_SUPERLU
    mat_superlu.ShowMessages();
#endif

#ifdef SELDON_WITH_UMFPACK
    mat_umf.ShowMessages();
#endif

#ifdef SELDON_WITH_PASTIX
    mat_pastix.ShowFullHistory();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
    mat_ilut.SetPrintLevel(1);
#endif
  }


  //! Clearing factorization.
  template<class T>
  void SparseDirectSolver<T>::Clear()
  {
    if (n > 0)
      {
	n = 0;
	mat_seldon.Clear();

#ifdef SELDON_WITH_MUMPS
	mat_mumps.Clear();
#endif

#ifdef SELDON_WITH_SUPERLU
	mat_superlu.Clear();
#endif

#ifdef SELDON_WITH_UMFPACK
	mat_umf.Clear();
#endif

#ifdef SELDON_WITH_PASTIX
	mat_pastix.Clear();
#endif

#ifdef SELDON_WITH_PRECONDITIONING
        mat_ilut.Clear();
#endif
      }
  }


  //! Returns the number of rows of the factorized matrix.
  template<class T>
  int SparseDirectSolver<T>::GetM() const
  {
    return n;
  }


  //! Returns the number of rows of the factorized matrix.
  template<class T>
  int SparseDirectSolver<T>::GetN() const
  {
    return n;
  }


  //! Returns the ordering algorithm to use.
  template<class T>
  int SparseDirectSolver<T>::GetTypeOrdering() const
  {
    return type_ordering;
  }


  //! Sets directly the new ordering (by giving a permutation vector).
  template<class T>
  void SparseDirectSolver<T>::SetPermutation(const IVect& num)
  {
    type_ordering = SparseMatrixOrdering::USER;
    permut = num;
  }


  //! Modifies the ordering algorithm to use.
  template<class T>
  void SparseDirectSolver<T>::SelectOrdering(int type)
  {
    type_ordering = type;
  }


  //! Modifies the number of threads per node (for Pastix only).
  template<class T>
  void SparseDirectSolver<T>::SetNumberThreadPerNode(int p)
  {
    number_threads_per_node = p;
  }


  //! Modifies the direct solver to use.
  template<class T>
  void SparseDirectSolver<T>::SelectDirectSolver(int type)
  {
    type_solver = type;
  }


  //! Enforces the use of unsymmetric algorithm for ilut solver.
  template<class T>
  void SparseDirectSolver<T>::SetNonSymmetricIlut()
  {
    enforce_unsym_ilut = true;
  }


  //! Returns the direct solver to use.
  template<class T>
  int SparseDirectSolver<T>::GetDirectSolver()
  {
    return type_solver;
  }


  //! Returns threshold used for ilut (if this solver is selected).
  template<class T>
  double SparseDirectSolver<T>::GetThresholdMatrix() const
  {
    return threshold_matrix;
  }


  //! Computation of the permutation vector in order to reduce fill-in.
  template<class T> template<class MatrixSparse>
  void SparseDirectSolver<T>::ComputeOrdering(MatrixSparse& A)
  {
    bool user_ordering = false;
    // We set the ordering for each direct solver interfaced.
    switch (type_ordering)
      {
      case SparseMatrixOrdering::AUTO :
	{
	  // We choose the default strategy proposed by the direct solver that
	  // will be called.
          if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      mat_mumps.SelectOrdering(7);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_PASTIX
	      mat_pastix.SelectOrdering(API_ORDER_SCOTCH);
#endif
	    }
	  else if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      mat_umf.SelectOrdering(UMFPACK_ORDERING_AMD);
#endif
	    }
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      mat_superlu.SelectOrdering(COLAMD);
#endif
	    }
	  else
	    {

              type_ordering = SparseMatrixOrdering::IDENTITY;

#ifdef SELDON_WITH_UMFPACK
              type_ordering = SparseMatrixOrdering::AMD;
#endif

#ifdef SELDON_WITH_MUMPS
              type_ordering = SparseMatrixOrdering::METIS;
#endif

#ifdef SELDON_WITH_PASTIX
              type_ordering = SparseMatrixOrdering::SCOTCH;
#endif

	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::IDENTITY :
      case SparseMatrixOrdering::REVERSE_CUTHILL_MCKEE :
      case SparseMatrixOrdering::USER :
	{
	  user_ordering = true;
	}
	break;
      case SparseMatrixOrdering::PORD :
      case SparseMatrixOrdering::QAMD :
	{
	  // Mumps ordering.
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      if (type_ordering == SparseMatrixOrdering::PORD)
		mat_mumps.SelectOrdering(4);
	      else
		mat_mumps.SelectOrdering(6);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::SCOTCH :
	{
	  // Available for Mumps and Pastix.
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      mat_mumps.SelectOrdering(3);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_SCOTCH
	      mat_pastix.SelectOrdering(API_ORDER_SCOTCH);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::METIS :
	{
	  // Available for Mumps and Pastix.
	  if (type_solver == MUMPS)
	    {
#ifdef SELDON_WITH_MUMPS
	      mat_mumps.SelectOrdering(5);
#endif
	    }
	  else if (type_solver == PASTIX)
	    {
#ifdef SELDON_WITH_SCOTCH
	      mat_pastix.SelectOrdering(API_ORDER_METIS);
#endif
	    }
	  else if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_SCOTCH
	      mat_umf.SelectOrdering(UMFPACK_ORDERING_METIS);
#endif
	    }

	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      case SparseMatrixOrdering::AMD :
      case SparseMatrixOrdering::COLAMD :
	{
	  // Default ordering for UmfPack.
	  if (type_solver == UMFPACK)
	    {
#ifdef SELDON_WITH_UMFPACK
	      mat_umf.SelectOrdering(UMFPACK_ORDERING_AMD);
#endif
	    }
	  else if (type_solver == SUPERLU)
	    {
#ifdef SELDON_WITH_SUPERLU
	      mat_superlu.SelectOrdering(COLAMD);
#endif
	    }
	  else
	    {
	      user_ordering = true;
	    }
	}
	break;
      }

    if (user_ordering)
      {
	// Case where the ordering is not natively available in the direct
	// solver computing separetely the ordering.
	FindSparseOrdering(A, permut, type_ordering);

        if (type_solver == MUMPS)
	  {
#ifdef SELDON_WITH_MUMPS
	    mat_mumps.SetPermutation(permut);
#endif
	  }
	else if (type_solver == PASTIX)
	  {
#ifdef SELDON_WITH_PASTIX
	    mat_pastix.SetPermutation(permut);
#endif
	  }
	else if (type_solver == UMFPACK)
	  {
#ifdef SELDON_WITH_UMFPACK
	    mat_umf.SetPermutation(permut);
#endif
	  }
	else if (type_solver == SUPERLU)
	  {
#ifdef SELDON_WITH_SUPERLU
	    mat_superlu.SetPermutation(permut);
#endif
	  }
	else
	  {
	  }
      }

  }


  //! Factorization of matrix A.
  /*!
    LU factorization is stored in the current object.
    You can ask to clear the matrix given on input (to spare memory).
  */
  template<class T> template<class MatrixSparse>
  void SparseDirectSolver<T>::Factorize(MatrixSparse& A, bool keep_matrix)
  {
    ComputeOrdering(A);

    n = A.GetM();
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	mat_umf.FactorizeMatrix(A, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with UmfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	mat_superlu.FactorizeMatrix(A, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	mat_mumps.FactorizeMatrix(A, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        mat_pastix.SetNumberThreadPerNode(number_threads_per_node);
	mat_pastix.FactorizeMatrix(A, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING
        // Setting some parameters.
        if (enforce_unsym_ilut || !IsSymmetricMatrix(A))
          mat_ilut.SetUnsymmetricAlgorithm();
        else
          mat_ilut.SetSymmetricAlgorithm();

        // Then performing factorization.
        mat_ilut.FactorizeMatrix(permut, A, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::Factorize(MatrixSparse&, bool)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	mat_seldon.FactorizeMatrix(permut, A);
      }

  }


  //! Returns error code of the direct solver (for Mumps only).
  template <class T>
  int SparseDirectSolver<T>::GetInfoFactorization(int& ierr) const
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
        ierr = mat_umf.GetInfoFactorization();
        switch (ierr)
          {
          case UMFPACK_OK :
            return FACTO_OK;
          case UMFPACK_WARNING_singular_matrix :
            return NUMERICALLY_SINGULAR_MATRIX;
          case UMFPACK_ERROR_out_of_memory :
            return OUT_OF_MEMORY;
          case UMFPACK_ERROR_invalid_Numeric_object :
          case UMFPACK_ERROR_invalid_Symbolic_object :
          case UMFPACK_ERROR_argument_missing :
          case UMFPACK_ERROR_different_pattern :
          case UMFPACK_ERROR_invalid_system :
            return INVALID_ARGUMENT;
          case UMFPACK_ERROR_n_nonpositive :
            return INCORRECT_NUMBER_OF_ROWS;
          case UMFPACK_ERROR_invalid_matrix :
            return MATRIX_INDICES_INCORRECT;
          case UMFPACK_ERROR_invalid_permutation :
            return INVALID_PERMUTATION;
          case UMFPACK_ERROR_ordering_failed :
            return ORDERING_FAILED;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
        ierr = mat_superlu.GetInfoFactorization();
        if (ierr > 0)
          return INTERNAL_ERROR;
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
        ierr = mat_mumps.GetInfoFactorization();
        switch (ierr)
          {
          case -2 :
            // nz out of range
            return MATRIX_INDICES_INCORRECT;
          case -3 :
            // invalid job number
            return INVALID_ARGUMENT;
          case -4 :
            // invalid permutation
            return INVALID_PERMUTATION;
          case -5 :
            // problem of real workspace allocation
            return OUT_OF_MEMORY;
          case -6 :
            // structurally singular matrix
            return STRUCTURALLY_SINGULAR_MATRIX;
          case -7 :
            // problem of integer workspace allocation
            return OUT_OF_MEMORY;
          case -10 :
            // numerically singular matrix
            return NUMERICALLY_SINGULAR_MATRIX;
          case -13 :
            // allocate failed
            return OUT_OF_MEMORY;
          case -16 :
            // N out of range
            return INCORRECT_NUMBER_OF_ROWS;
          case -22 :
            // invalid pointer
            return INVALID_ARGUMENT;
          case 1 :
            // index out of range
            return MATRIX_INDICES_INCORRECT;
	  case 0 :
	    return FACTO_OK;
          default :
            return INTERNAL_ERROR;
          }
#endif
      }

    return FACTO_OK;
  }


  //! x_solution is overwritten by solution of A x = b.
  /*!
    We assume that Factorize has been called previously.
  */
  template<class T> template<class Vector1>
  void SparseDirectSolver<T>::Solve(Vector1& x_solution)
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	Seldon::SolveLU(mat_umf, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with UmfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	Seldon::SolveLU(mat_superlu, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	Seldon::SolveLU(mat_mumps, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	Seldon::SolveLU(mat_pastix, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING
	mat_ilut.Solve(x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(Vector&)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	Seldon::SolveLU(mat_seldon, x_solution);
      }
  }


  //! x_solution is overwritten with solution of A x = b or A^T x = b.
  template<class T> template<class TransStatus, class Vector1>
  void SparseDirectSolver<T>
  ::Solve(const TransStatus& TransA, Vector1& x_solution)
  {
    if (type_solver == UMFPACK)
      {
#ifdef SELDON_WITH_UMFPACK
	Seldon::SolveLU(TransA, mat_umf, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with UmpfPack support.");
#endif
      }
    else if (type_solver == SUPERLU)
      {
#ifdef SELDON_WITH_SUPERLU
	Seldon::SolveLU(TransA, mat_superlu, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with SuperLU support.");
#endif
      }
    else if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	Seldon::SolveLU(TransA, mat_mumps, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	Seldon::SolveLU(TransA, mat_pastix, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else if (type_solver == ILUT)
      {
#ifdef SELDON_WITH_PRECONDITIONING
	mat_ilut.Solve(TransA, x_solution);
#else
        throw Undefined("SparseDirectSolver::Solve(TransStatus, Vector&)",
                        "Seldon was not compiled with the preconditioners.");
#endif
      }
    else
      {
	Seldon::SolveLU(TransA, mat_seldon, x_solution);
      }
  }


#ifdef SELDON_WITH_MPI
  //! Factorization of a matrix.
  /*! The matrix is given on each processor of the communicator in CSC
    format. If the matrix is assumed to be symmetric, you provide only the
    lower part of the matrix.
    \param[in] comm_facto communicator grouping processors involved in the
    factorization.
    \param[in,out] Ptr start indices
    \param[in,out] Row column indices
    \param[in,out] Val values of non-zero entries
    \param[in] glob_num global column numbers, each column of the global
    matrix is associated with one processor and only one.
    \param[in] sym if true, the matrix is assumed to be symmetric
    \param[in] keep_matrix if false the input matrix is erased
  */
  template<class T> template<class Tint>
  void SparseDirectSolver<T>::
  FactorizeDistributed(MPI::Comm& comm_facto,
                       Vector<Tint>& Ptr, Vector<Tint>& Row, Vector<T>& Val,
                       const IVect& glob_num, bool sym, bool keep_matrix)
  {
    n = Ptr.GetM()-1;
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	mat_mumps.FactorizeDistributedMatrix(comm_facto, Ptr, Row, Val,
                                             glob_num, sym, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        mat_pastix.FactorizeDistributedMatrix(comm_facto, Ptr, Row, Val,
                                              glob_num, sym, keep_matrix);
#else
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else
      {
        throw Undefined("SparseDirectSolver::FactorizeDistributed(MPI::Comm&,"
                        " IVect&, IVect&, Vector<T>&, IVect&, bool, bool)",
                        "The method is defined for Mumps and Pastix only.");
      }
  }


  //! Solution of linear system Ax = b by using LU factorization.
  /*!
    \param[in,out] x_solution on input right hand side, on output solution
  */
  template<class T> template<class Vector1>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI::Comm& comm_facto, Vector1& x_solution,
                   const IVect& glob_number)
  {
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	mat_mumps.SolveDistributed(comm_facto, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
        mat_pastix.SolveDistributed(comm_facto, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else
      {
        throw Undefined("SparseDirectSolver::SolveDistributed(MPI::Comm&,"
                        " Vector&, IVect&)",
                        "The method is defined for Mumps and Pastix only.");
      }

  }


  /*! \brief Solution of linear system A^T x = b by using LU factorization
    (without scaling). */
  /*!
    \param[in,out] x_solution on input right hand side, on output solution.
  */
  template<class T> template<class TransStatus, class Vector1>
  void SparseDirectSolver<T>::
  SolveDistributed(MPI::Comm& comm_facto, const TransStatus& TransA,
                   Vector1& x_solution, const IVect& glob_number)
  {
    if (type_solver == MUMPS)
      {
#ifdef SELDON_WITH_MUMPS
	mat_mumps.SolveDistributed(comm_facto, TransA, x_solution,
                                   glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "Seldon was not compiled with Mumps support.");
#endif
      }
    else if (type_solver == PASTIX)
      {
#ifdef SELDON_WITH_PASTIX
	mat_pastix.SolveDistributed(comm_facto, TransA, x_solution, glob_number);
#else
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "Seldon was not compiled with Pastix support.");
#endif
      }
    else
      {
        throw Undefined("SparseDirectSolver::SolveDistributed(TransStatus, "
                        "MPI::Comm&, Vector&, IVect&)",
                        "The method is defined for Mumps and Pastix only.");
      }
  }
#endif


  /*************************
   * Solve and SparseSolve *
   *************************/


  //! Solves a sparse linear system using LU factorization.
  /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
    \f$ X \f$ and \f$ Y \f$ are vectors.
    \param[in] M the sparse matrix of the linear system, to be factorized in
    LU form by UMFPACK, SuperLU or Mumps. On exit, \a M is cleared.
    \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
    solution \f$ X \f$ of the system.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                   Vector<T1, Storage1, Allocator1>& Y)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(M, Y, "SparseSolve(M, Y)");
#endif

    SparseDirectSolver<T0> matrix_lu;

    matrix_lu.Factorize(M);
    matrix_lu.Solve(Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M, Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  void Solve(Matrix<T, Prop0, ColSparse, Allocator0>& M,
             Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


  //! \copydoc SparseSolve(Matrix<T0, Prop0, Storage0, Allocator0>& M, Vector<T1, Storage1, Allocator1>& Y)
  template <class T, class Prop0, class Allocator0, class Allocator1>
  void Solve(Matrix<T, Prop0, RowSparse, Allocator0>& M,
             Vector<T, VectFull, Allocator1>& Y)
  {
    SparseSolve(M, Y);
  }


}  // namespace Seldon.


#endif
