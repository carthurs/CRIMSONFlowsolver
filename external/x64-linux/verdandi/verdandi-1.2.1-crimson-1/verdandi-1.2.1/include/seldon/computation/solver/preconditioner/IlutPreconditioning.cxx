// Copyright (C) 2010 Marc Duruflé
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


#ifndef SELDON_FILE_ILUT_PRECONDITIONING_CXX

#include "SymmetricIlutPreconditioning.cxx"

namespace Seldon
{

  template<class real, class cplx, class Allocator>
  IlutPreconditioning<real, cplx, Allocator>::IlutPreconditioning()
  {
    print_level = 0;
    symmetric_algorithm = false;
    type_ilu = ILUT;
    fill_level = 1000000;
    additional_fill = 1000000;
    mbloc = 1000000;
    alpha = 1.0;
    droptol = 0.01;
    permtol = 0.1;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::Clear()
  {
    permutation_row.Clear();
    mat_sym.Clear();
    mat_unsym.Clear();
    xtmp.Clear();
  }


  template<class real, class cplx, class Allocator>
  int IlutPreconditioning<real, cplx, Allocator>::GetFactorisationType() const
  {
    return type_ilu;
  }


  template<class real, class cplx, class Allocator>
  int IlutPreconditioning<real, cplx, Allocator>::GetFillLevel() const
  {
    return fill_level;
  }


  template<class real, class cplx, class Allocator>
  int IlutPreconditioning<real, cplx, Allocator>::GetAdditionalFillNumber()
    const
  {
    return additional_fill;
  }


  template<class real, class cplx, class Allocator>
  int IlutPreconditioning<real, cplx, Allocator>::GetPrintLevel() const
  {
    return print_level;
  }


  template<class real, class cplx, class Allocator>
  int IlutPreconditioning<real, cplx, Allocator>::GetPivotBlockInteger() const
  {
    return mbloc;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>
  ::SetFactorisationType(int type)
  {
    type_ilu = type;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::SetFillLevel(int level)
  {
    fill_level = level;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>
  ::SetAdditionalFillNumber(int level)
  {
    additional_fill = level;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::SetPrintLevel(int level)
  {
    print_level = level;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::SetPivotBlockInteger(int i)
  {
    mbloc = i;
  }


  template<class real, class cplx, class Allocator>
  real IlutPreconditioning<real, cplx, Allocator>
  ::GetDroppingThreshold() const
  {
    return droptol;
  }


  template<class real, class cplx, class Allocator>
  real IlutPreconditioning<real, cplx, Allocator>
  ::GetDiagonalCoefficient() const
  {
    return alpha;
  }


  template<class real, class cplx, class Allocator>
  real IlutPreconditioning<real, cplx, Allocator>::GetPivotThreshold() const
  {
    return permtol;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>
  ::SetDroppingThreshold(real tol)
  {
    droptol = tol;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>
  ::SetDiagonalCoefficient(real beta)
  {
    alpha = beta;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::SetPivotThreshold(real tol)
  {
    permtol = tol;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::SetSymmetricAlgorithm()
  {
    symmetric_algorithm = true;
  }


  template<class real, class cplx, class Allocator>
  void IlutPreconditioning<real, cplx, Allocator>::SetUnsymmetricAlgorithm()
  {
    symmetric_algorithm = false;
  }


  template<class real, class cplx, class Allocator>
  template<class T0, class Storage0, class Allocator0>
  void IlutPreconditioning<real, cplx, Allocator>::
  FactorizeMatrix(const IVect& perm,
                  Matrix<T0, General, Storage0, Allocator0>& mat,
                  bool keep_matrix)
  {
    if (symmetric_algorithm)
      {
        cout << "Conversion to symmetric matrices not implemented." << endl;
        abort();
      }
    else
      FactorizeUnsymMatrix(perm, mat, keep_matrix);
  }


  template<class real, class cplx, class Allocator>
  template<class T0, class Storage0, class Allocator0>
  void IlutPreconditioning<real, cplx, Allocator>::
  FactorizeMatrix(const IVect& perm,
                  Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
                  bool keep_matrix)
  {
    if (symmetric_algorithm)
      FactorizeSymMatrix(perm, mat, keep_matrix);
    else
      FactorizeUnsymMatrix(perm, mat, keep_matrix);
  }


  template<class real, class cplx, class Allocator>
  template<class MatrixSparse>
  void IlutPreconditioning<real, cplx, Allocator>::
  FactorizeSymMatrix(const IVect& perm, MatrixSparse& mat, bool keep_matrix)
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
      {
        cout << "Numbering array should have the same size as matrix.";
        cout << endl;
        abort();
      }

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
        {
          cout << "Error in the numbering array." << endl;
          abort();
        }

    // Matrix is permuted.
    ApplyInversePermutation(mat_sym, perm, perm);

    // Temporary vector used for solving.
    xtmp.Reallocate(n);

    // Factorization is performed.
    GetIlut(*this, mat_sym);
  }


  template<class real, class cplx, class Allocator>
  template<class MatrixSparse>
  void IlutPreconditioning<real, cplx, Allocator>::
  FactorizeUnsymMatrix(const IVect& perm, MatrixSparse& mat, bool keep_matrix)
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
      {
        cout << "Numbering array should have the same size as matrix.";
        cout << endl;
        abort();
      }

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
        {
          cout << "Error in the numbering array." << endl;
          abort();
        }

    IVect iperm = inv_permutation;

    // Rows of matrix are permuted.
    ApplyInversePermutation(mat_unsym, perm, perm);

    // Temporary vector used for solving.
    xtmp.Reallocate(n);

    // Factorization is performed.
    // Columns are permuted during the factorization.
    inv_permutation.Fill();
    GetIlut(*this, mat_unsym, permutation_col, inv_permutation);

    // Combining permutations.
    IVect itmp = permutation_col;
    for (int i = 0; i < n; i++)
      permutation_col(i) = iperm(itmp(i));

    permutation_row = perm;
  }


  template<class real, class cplx, class Allocator>
  template<class Matrix1, class Vector1>
  void IlutPreconditioning<real, cplx, Allocator>::
  Solve(const Matrix1& A, const Vector1& r, Vector1& z)
  {
    if (symmetric_algorithm)
      {
        for (int i = 0; i < r.GetM(); i++)
          xtmp(permutation_row(i)) = r(i);

        SolveLU(mat_sym, xtmp);

        for (int i = 0; i < r.GetM(); i++)
          z(i) = xtmp(permutation_row(i));
      }
    else
      {
        for (int i = 0; i < r.GetM(); i++)
          xtmp(permutation_row(i)) = r(i);

        SolveLU(mat_unsym, xtmp);

        for (int i = 0; i < r.GetM(); i++)
          z(permutation_col(i)) = xtmp(i);
      }
  }


  template<class real, class cplx, class Allocator>
  template<class Matrix1, class Vector1>
  void IlutPreconditioning<real, cplx, Allocator>::
  TransSolve(const Matrix1& A, const Vector1& r, Vector1& z)
  {
    if (symmetric_algorithm)
      Solve(A, r, z);
    else
      {
        for (int i = 0; i < r.GetM(); i++)
          xtmp(i) = r(permutation_col(i));

        SolveLU(SeldonTrans, mat_unsym, xtmp);

        for (int i = 0; i < r.GetM(); i++)
          z(i) = xtmp(permutation_row(i));
      }
  }


  template<class real, class cplx, class Allocator>
  template<class Vector1>
  void IlutPreconditioning<real, cplx, Allocator>::Solve(Vector1& z)
  {
    if (symmetric_algorithm)
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


  template<class real, class cplx, class Allocator>
  template<class Vector1>
  void IlutPreconditioning<real, cplx, Allocator>::TransSolve(Vector1& z)
  {
    if (symmetric_algorithm)
      Solve(z);
    else
      {
        for (int i = 0; i < z.GetM(); i++)
          xtmp(i) = z(permutation_col(i));

        SolveLU(SeldonTrans, mat_unsym, xtmp);

        for (int i = 0; i < z.GetM(); i++)
          z(i) = xtmp(permutation_row(i));
      }
  }


  template<class real, class cplx, class Allocator>
  template<class TransStatus, class Vector1>
  void IlutPreconditioning<real, cplx, Allocator>
  ::Solve(const TransStatus& transA, Vector1& z)
  {
    if (transA.Trans())
      TransSolve(z);
    else
      Solve(z);
  }


  template<class real, class cplx, class Storage, class Allocator>
  void qsplit_ilut(Vector<cplx, Storage, Allocator>& a, IVect& ind, int first,
                   int n, int ncut, const real& abs_ncut)
  {
    //-----------------------------------------------------------------------
    //     does a quick-sort split of a real array.
    //     on input a(1:n). is a real array
    //     on output a(1:n) is permuted such that its elements satisfy:
    //
    //     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
    //     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
    //
    //     ind is an integer array which permuted in the same way as a
    //-----------------------------------------------------------------------
    int last = n-1;
    int ncut_ = ncut-1;
    int first_ = first;

    if ((ncut_ < first_) || (ncut_ > last))
      return;

    cplx tmp; real abskey;
    int mid, itmp;
    bool test_loop = true;

    //     outer loop -- while mid .ne. ncut do
    while (test_loop)
      {
	mid = first_;
        abskey = abs(a(mid));
        for (int j = (first_+1); j <= last; j++)
	  {
	    if (abs(a(j)) > abskey)
	      {
		mid++;
		// Interchange.
		tmp = a(mid);
		itmp = ind(mid);
		a(mid) = a(j);
		ind(mid) = ind(j);
		a(j)  = tmp;
		ind(j) = itmp;
              }
	  }

	// Interchange.
        tmp = a(mid);
        a(mid) = a(first_);
        a(first_)  = tmp;

        itmp = ind(mid);
        ind(mid) = ind(first_);
        ind(first_) = itmp;

	// Test for while loop.
        if (mid == ncut_)
	  return;

        if (mid > ncut_)
	  last = mid-1;
        else
	  first_ = mid+1;
      }
  }


  //! Incomplete factorization with pivot for unsymmetric matrix.
  template<class real, class cplx, class Allocator1, class Allocator2>
  void GetIlut(const IlutPreconditioning<real, cplx, Allocator1>& param,
               Matrix<cplx, General, ArrayRowSparse, Allocator2>& A,
               IVect& iperm, IVect& rperm)
  {
    int size_row;
    int n = A.GetN();
    int type_factorization = param.GetFactorizationType();
    int lfil = param.GetFillLevel();
    real zero(0);
    real droptol = param.GetDroppingThreshold();
    real alpha = param.GetDiagonalCoefficient();
    bool variable_fill = false;
    bool standard_dropping = true;
    int additional_fill = param.GetAdditionalFillNumber();
    int mbloc = param.GetPivotBlockInteger();
    real permtol = param.GetPivotThreshold();
    int print_level = param.GetPrintLevel();

    if (type_factorization == param.ILUT)
      standard_dropping = false;
    else if (type_factorization == param.ILU_D)
      standard_dropping = true;
    else if (type_factorization == param.ILUT_K)
      {
	variable_fill = true;   // We use a variable lfil
	standard_dropping = false;
      }
    else if (type_factorization == param.ILU_0)
      {
	GetIlu0(A);
        return;
      }
    else if (type_factorization == param.MILU_0)
      {
	GetMilu0(A);
        return;
      }
    else if (type_factorization == param.ILU_K)
      {
	GetIluk(lfil, A);
        return;
      }

    cplx fact, s, t;
    real tnorm;
    int length_lower, length_upper, jpos, jrow, i_row, j_col;
    int i, j, k, index_lu, length;


    if (lfil < 0)
      {
        cout << "Incorrect fill level." << endl;
        abort();
      }

    cplx czero, cone;
    SetComplexZero(czero);
    SetComplexOne(cone);
    typedef Vector<cplx, VectFull, Allocator2> VectCplx;
    VectCplx Row_Val(n);
    IVect Index(n), Row_Ind(n), Index_Diag(n);
    Row_Val.Fill(czero);
    Row_Ind.Fill(-1);
    Index_Diag.Fill(-1);

    Index.Fill(-1);

    bool element_dropped; cplx dropsum;

    // main loop
    int new_percent = 0, old_percent = 0;
    for (i_row = 0; i_row < n; i_row++)
      {
        // Progress bar if print level is high enough.
        if (print_level > 0)
          {
            new_percent = int(double(i_row)/(n-1)*80);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#"; cout.flush();
              }

            old_percent = new_percent;
          }

	size_row = A.GetRowSize(i_row);
	tnorm = zero;

	dropsum = czero;
	for (k = 0 ; k < size_row; k++)
          if (A.Value(i_row, k) != czero)
            tnorm += abs(A.Value(i_row, k));

	if (tnorm == zero)
	  {
            cout << "Structurally singular matrix." << endl;
            cout << "Norm of row " << i_row << " is equal to 0." << endl;
            abort();
          }

        // tnorm is the sum of absolute value of coefficients of row i_row.
	tnorm /= real(size_row);
	if (variable_fill)
	  lfil = size_row + additional_fill;

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
	      {
		if (Row_Ind(j) < jrow)
		  {
		    jrow = Row_Ind(j);
		    k = j;
		  }
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

	    element_dropped = false;
	    if (standard_dropping)
	      if (abs(Row_Val(j_col)) <= droptol*tnorm)
		{
		  dropsum += Row_Val(j_col);
		  element_dropped = true;
		}

            // Gets the multiplier for row to be eliminated (jrow).
	    if (!element_dropped)
	      {
		// first_index_upper points now on the diagonal coefficient.
		fact = Row_Val(j_col) * A.Value(jrow, Index_Diag(jrow));

		if (!standard_dropping)
		  {
		    if (abs(fact) <= droptol)
		      element_dropped = true;
		  }
	      }

	    if (!element_dropped)
	      {
		// Combines current row and row jrow.
		for (k = (Index_Diag(jrow)+1); k < A.GetRowSize(jrow); k++)
		  {
		    s = fact * A.Value(jrow,k);
                    j = rperm(A.Index(jrow,k));

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
                            // this is a fill-in element
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

		// Stores this pivot element from left to right -- no danger
		// of overlap with the working elements in L (pivots).
                Row_Val(length) = fact;
		Row_Ind(length) = jrow;
		++length;
	      }

	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row+k )) = -1;

	// Updates L-matrix.
	if (!standard_dropping)
	  {
	    length_lower = length;
	    length = min(length_lower,lfil);

	    // Sorts by quick-split.
            qsplit_ilut(Row_Val, Row_Ind, 0, length_lower, length,tnorm);
          }

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
	    if (abs(Row_Val(i_row+k)) > droptol * tnorm)
	      {
		++length;
		Row_Val(i_row+length) = Row_Val(i_row+k);
		Row_Ind(i_row+length) = Row_Ind(i_row+k);
	      }
	    else
	      dropsum += Row_Val(i_row+k);
	  }

        if (!standard_dropping)
	  {
	    length_upper = length + 1;
	    length = min(length_upper,lfil);

            qsplit_ilut(Row_Val, Row_Ind, i_row+1,
                        i_row+length_upper, i_row+length+1, tnorm);
	  }
	else
	  length++;

	// Determines next pivot.
        int imax = i_row;
        real xmax = abs(Row_Val(imax));
        real xmax0 = xmax;
        int icut = i_row + mbloc - 1 - i_row%mbloc;
        for ( k = i_row + 1; k <= i_row + length - 1; k++)
          {
            tnorm = abs(Row_Val(k));
            if ((tnorm > xmax) && (tnorm*permtol > xmax0)
                && (Row_Ind(k)<= icut))
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
	if (standard_dropping)
	  Row_Val(i_row) += alpha*dropsum;

	if (Row_Val(i_row) == czero)
          Row_Val(i_row) = (droptol + 1e-4) * tnorm;

	A.Value(i_row, index_diag) = cone / Row_Val(i_row);

      } // end main loop

    if (print_level > 0)
      cout << endl;

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize()*(sizeof(cplx)+4))/(1024*1024)) << " MB" << endl;

    for (i = 0; i < n; i++ )
      for (j = 0; j < A.GetRowSize(i); j++)
        A.Index(i,j) = rperm(A.Index(i,j));
  }


  template<class cplx, class Allocator>
  void GetIluk(int lfil, Matrix<cplx, General, ArrayRowSparse, Allocator>& A)
  {
    int n = A.GetM();
    Vector<cplx, VectFull, CallocAlloc<cplx> > w;
    w.Reallocate(n+1);
    IVect jw(3*n), Index_Diag(n);
    Vector<IVect, VectFull, NewAlloc<IVect> > levs(n);

    cplx czero, cone;
    SetComplexZero(czero);
    SetComplexOne(cone);

    // Local variables
    cplx fact, s, t;
    int length_lower,length_upper, jpos, jrow, i_row, j_col;
    int i, j, k, index_lu, length;
    bool element_dropped;

    int n2 = 2*n, jlev, k_, size_upper;
    jw.Fill(-1);

    // Main loop.
    for (i_row = 0; i_row < n; i_row++)
      {
	int size_row = A.GetRowSize(i_row);

	// Unpacks L-part and U-part of row of A in arrays w, jw.
	length_upper = 1;
	length_lower = 0;
	jw(i_row) = i_row;
	w(i_row) = 0.0;
	jw(n + i_row) = i_row;

	for (j = 0; j < size_row; j++)
	  {
	    k = A.Index(i_row,j);
	    t = A.Value(i_row,j);
	    if (k < i_row)
	      {
		jw(length_lower) = k;
		w(length_lower) = t;
		jw(n + k) = length_lower;
		jw(n2+length_lower) = -1;
		++length_lower;
	      }
	    else if (k == i_row)
	      {
		w(i_row) = t;
		jw(n2+length_lower) = -1;
	      }
	    else
	      {
		jpos = i_row + length_upper;
		jw(jpos) = k;
		w(jpos) = t;
		jw(n + k) = jpos;
		length_upper++;
	      }
	  }

	j_col = 0;
        length = 0;


        // Eliminates previous rows.
	while (j_col <length_lower)
	  {
            // In order to do the elimination in the correct order, we must
            // select the smallest column index among jw(k); k = j_col + 1,
            // ..., length_lower.
	    jrow = jw(j_col);
	    k = j_col;

	    // Determines smallest column index.
	    for (j = (j_col+1) ; j < length_lower; j++)
	      {
		if (jw(j) < jrow)
		  {
		    jrow = jw(j);
		    k = j;
		  }
              }

	    if (k != j_col)
	      {
		// Exchanges in jw.
		j = jw(j_col);
		jw(j_col) = jw(k);
		jw(k) = j;

                // Exchanges in jw(n+  (pointers/ nonzero indicator).
		jw(n+jrow) = j_col;
		jw(n+j) = k;

		// Exchanges in jw(n2+  (levels).
		j = jw(n2+j_col);
		jw(n2+j_col)  = jw(n2+k);
		jw(n2+k) = j;

                // Exchanges in w.
		s = w(j_col);
		w(j_col) = w(k);
		w(k) = s;
	      }

	    // Zero out element in row by setting jw(n+jrow) to zero.
	    jw(n + jrow) = -1;

	    element_dropped = false;

	    // Gets the multiplier for row to be eliminated (jrow).
	    fact = w(j_col) * A.Value(jrow,Index_Diag(jrow));

	    jlev = jw(n2+j_col) + 1;
	    if (jlev > lfil)
	      element_dropped = true;

	    if (!element_dropped)
	      {
		// Combines current row and row jrow.
		k_ = 0;
		for (k = (Index_Diag(jrow)+1); k < A.GetRowSize(jrow) ; k++)
		  {
		    s = fact * A.Value(jrow,k);
		    j = A.Index(jrow,k);

		    jpos = jw(n + j);
		    if (j >= i_row)
		      {
                        // Dealing with upper part.
                        if (jpos == -1)
			  {
                            // This is a fill-in element.
			    i = i_row + length_upper;
			    jw(i) = j;
			    jw(n + j) = i;
			    w(i) = -s;

			    jw(n2+i) = jlev + levs(jrow)(k_) + 1;
			    ++length_upper;
			  }
			else
			  {
                            // This is not a fill-in element.
                            w(jpos) -= s;
			    jw(n2+jpos) = min(jw(n2+jpos),
                                              jlev + levs(jrow)(k_)+1);
			  }
		      }
		    else
		      {
			// Dealing  with lower part.
			if (jpos == -1)
			  {
                            // This is a fill-in element.
			    jw(length_lower) = j;
			    jw(n + j) = length_lower;
			    w(length_lower) = -s;
			    jw(n2+length_lower) = jlev + levs(jrow)(k_) + 1;
			    ++length_lower;
			  }
			else
			  {
                            // This is not a fill-in element.
			    w(jpos) -= s;
			    jw(n2+jpos) = min(jw(n2 + jpos),
                                              jlev + levs(jrow)(k_) + 1);
			  }
		      }

                    k_++;
		  }

	      }

            // Stores this pivot element from left to right -- no danger of
            // overlap with the working elements in L (pivots).
	    w(j_col) = fact;
	    jw(j_col) = jrow;

	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          jw(n + jw(i_row + k )) = -1;

	// Updates L-matrix.
	size_row = 1; // we have the diagonal value.
	// Size of L-matrix.
	for (k = 0; k < length_lower; k++)
	  if (jw(n2+k) < lfil)
	    size_row++;

	// Size of U-matrix.
	size_upper = 0;
	for (k = (i_row+1) ; k <= (i_row+length_upper-1) ; k++)
	  if (jw(n2+k) < lfil)
	    size_upper++;

	size_row += size_upper;
	A.ReallocateRow(i_row,size_row);
	levs(i_row).Reallocate(size_upper);

	index_lu = 0;
	for (k = 0; k < length_lower; k++)
	  {
	    if (jw(n2+k) < lfil)
	      {
		A.Value(i_row,index_lu) = w(k);
		A.Index(i_row,index_lu) = jw(k);
                ++index_lu;
	      }
	  }

	// Saves pointer to beginning of row i_row of U.
	Index_Diag(i_row) = index_lu;
	A.Value(i_row,index_lu) = cone / w(i_row);
	A.Index(i_row,index_lu++) = i_row;

	for (k = (i_row+1) ; k <= (i_row+length_upper-1) ; k++)
	  {
	    if (jw(n2+k) < lfil)
	      {
		A.Index(i_row,index_lu) = jw(k);
		A.Value(i_row,index_lu) = w(k);
		levs(i_row)(index_lu-Index_Diag(i_row)-1) = jw(n2+k);
                ++index_lu;
	      }
	  }
      }
  }


  template<class cplx, class Allocator>
  void GetIlu0(Matrix<cplx, General, ArrayRowSparse, Allocator>& A)
  {
    int j_col, jrow, jw, n = A.GetM();
    IVect Index(n), ju(n);

    cplx czero, cone;
    SetComplexZero(czero);
    SetComplexOne(cone);

    // Initializes work vector to zero's.
    Index.Fill(-1); ju.Fill(-1);
    cplx tl;

    // Main loop.
    for (int i_row = 0 ; i_row < n ; i_row++)
      {
	// Generating row number i_row of L and U.
        for (int j = 0 ; j < A.GetRowSize(i_row) ; j++ )
	  {
	    j_col = A.Index(i_row, j);
	    if (j_col == i_row)
	      ju(i_row) = j;

	    Index(j_col) = j;
	  }

	int jm = ju(i_row)-1;

	// Exit if diagonal element is reached.
	for (int j = 0; j <= jm; j++)
	  {
	    jrow = A.Index(i_row, j);
	    tl = A.Value(i_row, j)*A.Value(jrow, ju(jrow));
	    A.Value(i_row, j) = tl;

	    // Performs linear combination.
            for ( j_col = (ju(jrow)+1); j_col < A.GetRowSize(jrow); j_col++)
	      {
		jw = Index(A.Index(jrow,j_col));
		if (jw != -1)
		  A.Value(i_row, jw) -= tl*A.Value(jrow, j_col);
	      }
	  }


        // Inverts and stores diagonal element.
        if (A.Value(i_row, ju(i_row)) == czero)
          {
            cout << "Factorization fails because we found a null coefficient"
                 << " on diagonal " << i_row << endl;
            abort();
          }

        A.Value(i_row,ju(i_row)) = cone / A.Value(i_row,ju(i_row));

        // Resets pointer Index to zero.
        Index(i_row) = -1;
        for (int i = 0; i < A.GetRowSize(i_row); i++)
          Index(A.Index(i_row, i)) = -1;
      }

  }


  template<class cplx, class Allocator>
  void GetMilu0(Matrix<cplx, General, ArrayRowSparse, Allocator>& A)
  {
    int j_col, jrow, jw, n = A.GetM();
    IVect Index(n), ju(n);

    cplx czero, cone;
    SetComplexZero(czero);
    SetComplexOne(cone);

    // Initializes work vector to zero's.
    Index.Fill(-1); ju.Fill(-1);
    cplx tl;

    // Main loop.
    for (int i_row = 0 ; i_row < n ; i_row++)
      {
        // Generating row number i_row of L and U.
        for (int j = 0; j < A.GetRowSize(i_row); j++ )
	  {
	    j_col = A.Index(i_row, j);
	    if (j_col == i_row)
	      ju(i_row) = j;

	    Index(j_col) = j;
	  }

	int jm = ju(i_row)-1;
	// Exit if diagonal element is reached.
	// s accumulates fill-in values.
	cplx s(0);
	for (int j = 0; j <= jm; j++)
	  {
	    jrow = A.Index(i_row, j);
	    tl = A.Value(i_row, j)*A.Value(jrow, ju(jrow));
	    A.Value(i_row, j) = tl;

            // Performs linear combination.
            for ( j_col = (ju(jrow)+1); j_col < A.GetRowSize(jrow); j_col++ )
	      {
		jw = Index(A.Index(jrow, j_col));
		if (jw != -1)
		  A.Value(i_row, jw) -= tl*A.Value(jrow, j_col);
		else
		  s += tl*A.Value(jrow, j_col);
	      }
	  }

	// Inverts and stores diagonal element.
	A.Value(i_row, ju(i_row)) -= s;
	if (A.Value(i_row, ju(i_row)) == czero)
          {
            cout << "Factorization fails because we found a null coefficient"
                 << " on diagonal " << i_row << endl;
            abort();
          }

	A.Value(i_row, ju(i_row)) = cone /A.Value(i_row, ju(i_row));

        // Resets pointer Index to zero.
        Index(i_row) = -1;
	for (int i = 0; i < A.GetRowSize(i_row); i++)
	  Index(A.Index(i_row, i)) = -1;
      }
  }

}

#define SELDON_FILE_ILUT_PRECONDITIONING_CXX
#endif
