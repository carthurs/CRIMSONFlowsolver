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


#ifndef SELDON_FILE_SYMMETRIC_ILUT_PRECONDITIONING_CXX

namespace Seldon
{

  //! Incomplete Factorization without pivot for symmetric matrix.
  template<class real, class cplx, class Allocator1, class Allocator2>
  void GetIlut(const IlutPreconditioning<real, cplx, Allocator1>& param,
               Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator2>& A)
  {
    int size_row;
    int n = A.GetN();
    int lfil = param.GetFillLevel();
    real zero = 0.0;
    real droptol = param.GetDroppingThreshold();
    real alpha = param.GetDiagonalCoefficient();
    bool variable_fill = false;
    bool standard_dropping = true;
    int type_factorization = param.GetFactorizationType();
    int additional_fill = param.GetAdditionalFillNumber();
    int print_level = param.GetPrintLevel();
    if (type_factorization == param.ILUT)
      standard_dropping = false;
    else if (type_factorization == param.ILU_D)
      standard_dropping = true;
    else if (type_factorization == param.ILUT_K)
      {
        variable_fill = true;
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

    cplx fact, s, t; real tnorm;
    int length_lower, length_upper, jpos, jrow;
    int i_row, j_col, index_lu, length;
    int i, j, k;

    if (lfil < 0)
      {
        cout << "Incorrect fill level." << endl;
        abort();
      }

    typedef Vector<cplx, VectFull, Allocator2> VectCplx;
    VectCplx Row_Val(n);
    IVect Index(n), Row_Ind(n);
    Row_Val.Zero(); Row_Ind.Fill(-1);

    Index.Fill(-1);

    bool element_dropped; cplx dropsum;

    // We convert A into an unsymmetric matrix.
    Matrix<cplx, General, ArrayRowSparse, Allocator2> B;
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
            new_percent = int(real(i_row)/(n-1)*80);
            for (int percent = old_percent; percent < new_percent; percent++)
              {
                cout << "#"; cout.flush();
              }

            old_percent = new_percent;
          }

        // 1-norm of the row of initial matrix.
	size_row = B.GetRowSize(i_row);
	tnorm = zero;
	dropsum = zero;
	for (k = 0 ; k < size_row; k++)
          tnorm += abs(B.Value(i_row, k));

	if (tnorm == zero)
	  {
            cout << "Structurally singular matrix." << endl;
            cout << "Norm of row " << i_row << " is equal to 0." << endl;
            abort();
          }

        // tnorm is the sum of absolute value of coefficients of row i_row
	tnorm /= real(size_row);
	if (variable_fill)
	  lfil = size_row + additional_fill;


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

	    element_dropped = false;
	    if (standard_dropping)
	      if (abs(Row_Val(j_col)) <= droptol*tnorm)
		{
		  dropsum += Row_Val(j_col);
		  element_dropped = true;
		}

            // Gets the multiplier for row to be eliminated.
	    if (!element_dropped)
	      {
                fact = Row_Val(j_col) * A.Value(jrow, 0);

		if (!standard_dropping)
		  {
		    if (abs(fact) <= droptol)
		      element_dropped = true;
		  }
	      }

	    if (!element_dropped)
	      {
		// Combines current row and row jrow.
		for (k = 1; k < A.GetRowSize(jrow); k++)
		  {
		    s = fact * A.Value(jrow,k);
		    j = A.Index(jrow,k);

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
	      }

	    j_col++;
	  }

	// Resets double-pointer to zero (U-part).
	for (k = 0; k < length_upper; k++)
          Index(Row_Ind(i_row+k )) = -1;

	// Updating U-matrix -- first apply dropping strategy.

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
	    length = min(length_upper, lfil);

            qsplit_ilut(Row_Val, Row_Ind, i_row+1,
                        i_row+length_upper, i_row+length+1, tnorm);
          }
	else
	  length++;

	// Copies U-part in matrix A.
	A.ReallocateRow(i_row, length);
	index_lu = 1;
	for (k = (i_row+1) ; k <= (i_row+length-1) ; k++)
	  {
	    A.Index(i_row,index_lu) = Row_Ind(k);
	    A.Value(i_row,index_lu) = Row_Val(k);
	    ++index_lu;
	  }

	// Stores the inverse of the diagonal element of u.
	if (standard_dropping)
	  Row_Val(i_row) += alpha*dropsum;

	if (Row_Val(i_row) == zero)
          Row_Val(i_row) = (droptol + 1e-4) * tnorm;

	A.Value(i_row,0) = 1.0 / Row_Val(i_row);

      } // end main loop.

    if (print_level > 0)
      cout<<endl;

    // for each row of A, we divide by diagonal value
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= A.Value(i,0);

    if (print_level > 0)
      cout << "The matrix takes " <<
        int((A.GetDataSize()*(sizeof(cplx)+4))/(1024*1024)) << " MB" << endl;

  }


  template<class cplx, class Allocator>
  void GetIluk(int lfil,
               Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    int n = A.GetM();
    Vector<cplx, VectFull, CallocAlloc<cplx> > w;
    w.Reallocate(n+1);
    IVect jw(3*n);
    Vector<IVect, VectFull, NewAlloc<IVect> > levs(n);

    // Local variables.
    cplx fact, s, t;
    int length_lower, length_upper, jpos, jrow, i_row, j_col;
    int i, j, k, index_lu, length;
    bool element_dropped;

    // Initializes nonzero indicator array.
    int n2 = 2*n, jlev, k_, size_upper;
    jw.Fill(-1);

    // We convert A into an unsymmetric matrix.
    Matrix<cplx, General, ArrayRowSparse, Allocator> B;
    Seldon::Copy(A,B);

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
	    k = B.Index(i_row,j);
	    t = B.Value(i_row,j);
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

        B.ClearRow(i_row);

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

	    if (k!=j_col)
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
	    fact = w(j_col) * A.Value(jrow,0);

	    jlev = jw(n2+j_col) + 1;
	    if (jlev > lfil)
	      element_dropped = true;

	    if (!element_dropped)
	      {
		// Combines current row and row jrow.
		k_ = 0;
		for (k = 1; k < A.GetRowSize(jrow) ; k++)
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
                                              jlev+levs(jrow)(k_)+1);
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
			    jw(n2+jpos) = min(jw(n2+jpos),
                                              jlev+levs(jrow)(k_)+1);
			  }
		      }
                    k_++;
		  }

	      }

            // Stores this pivot element from left to right : no danger of
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

	// Size of U-matrix.
	size_upper = 0;
	for (k = (i_row+1) ; k <= (i_row+length_upper-1) ; k++)
	  if (jw(n2+k) < lfil)
	    size_upper++;

	size_row += size_upper;
	A.ReallocateRow(i_row, size_row);
	levs(i_row).Reallocate(size_upper);

	index_lu = 0;

	A.Value(i_row,index_lu) = 1.0 / w(i_row);
	A.Index(i_row,index_lu++) = i_row;

        // Updates U-matrix -- first apply dropping strategy.
	for (k = (i_row+1) ; k <= (i_row+length_upper-1) ; k++)
	  {
	    if (jw(n2+k) < lfil)
	      {
		A.Index(i_row,index_lu) = jw(k);
		A.Value(i_row,index_lu) = w(k);
		levs(i_row)(index_lu-1) = jw(n2+k);
                ++index_lu;
	      }
	  }
      } // End main loop.

    // For each row of A, we divide by diagonal value.
    for (int i = 0; i < n; i++)
      for (int j = 1; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= A.Value(i,0);

  }


  template<class cplx, class Allocator>
  void GetIlu0(Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    cout << "Not implemented." << endl;
    abort();
  }


  template<class cplx, class Allocator>
  void GetMilu0(Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    cout << "Not implemented." << endl;
    abort();
  }


} // end namespace

#define SELDON_SYMMETRIC_ILUT_PRECONDITIONING_CXX
#endif
