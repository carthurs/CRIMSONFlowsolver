// Copyright (C) 2001-2010 INRIA, Vivien Mallet
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_CXX

#include "Functions.hxx"

#include "../computation/basic_functions/Functions_Vector.cxx"

namespace Seldon
{


  //! Extracts a row from a sparse matrix.
  /*!
    \param[in] M sparse matrix.
    \param[in] i row index.
    \param[out] X extracted row.
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetRow(const Matrix<T0, General, RowSparse, Allocator0>& M,
	      int i, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int m = M.GetM();
    if (i < 0 || i >= m)
      throw WrongIndex("GetRow()",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    double* data = M.GetData();
    int size_row = ptr[i+1] - ptr[i];

    X.Reallocate(size_row);
    int shift = ptr[i];
    for (int j = 0; j < size_row; j++)
      {
	X.Index(j) = ind[shift + j];
	X.Value(j) = data[shift + j];
      }
  }


  //! Extracts a row from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] i row index.
    \param[out] X extracted row.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetRow(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int i, Vector<T1, Storage1, Allocator1>& X)
  {
    X.Reallocate(M.GetN());
    for (int j = 0; j < M.GetN(); j++)
      X(j) = M(i, j);
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void GetCol(const Matrix<T0, General, RowSparse, Allocator0>& M,
	      int j, Vector<T1, VectSparse, Allocator1>& X)
  {
#ifdef SELDON_CHECK_BOUNDS
    int n = M.GetN();
    if (j < 0 || j >= n)
      throw WrongIndex("GetCol()",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    double* data = M.GetData();
    int m = M.GetM();
    Vector<int> index;
    Vector<double> value;

    for (int i = 0; i < m; i++)
      for (int k = ptr[i]; k < ptr[i+1]; k++)
	if (ind[k] == j)
	  {
	    index.PushBack(i);
	    value.PushBack(data[k]);
	  }

    X.SetData(value, index);
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M,
              int j, Vector<T1, PETScSeq, Allocator1>& X)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M,
              int j, Vector<T1, PETScPar, Allocator1>& X)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      X.SetBuffer(i, M(i, j));
    X.Flush();
  }


  //! Extracts a column from a matrix.
  /*!
    \param[in] M matrix.
    \param[in] j column index.
    \param[out] X extracted column.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int j, Vector<T1, Storage1, Allocator1>& X)
  {
    X.Reallocate(M.GetM());
    for (int i = 0; i < M.GetM(); i++)
      X(i) = M(i, j);
  }


  //! Extracts columns of a matrix.
  /*! Columns [\a begin, \a end[ of \a M_in are returned in \a M_out.
    \param[in] M_in input matrix.
    \param[in] begin first column of \a M_in to extract.
    \param[in] end the last column to be extracted from \a M_in is \a end - 1.
    \param[out] M_out on exit, matrix composed of the columns \a begin to
    \a end - 1 of \a M_in. \a M_out is reallocated if necessary.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M_in,
	      int begin, int end,
              Matrix<T1, Prop1, Storage1, Allocator1>& M_out)
  {
    M_out.Reallocate(M_in.GetM(), end - begin);
    for (int i = 0; i < M_in.GetM(); i++)
      for (int j = begin, k = 0; j < end; j++, k++)
        M_out(i, k) = M_in(i, j);
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, Storage0, Allocator0>& M)
  {
    for (int j = 0; j < M.GetN(); j++)
      M.Set(i, j, X(j));
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScSeq, Allocator1>& X,
              int i, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M)
  {
    for (int j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X(j));
    M.Flush();
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetRow(const Vector<T1, PETScPar, Allocator1>& X,
              int i, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    for (int j = 0; j < M.GetN(); j++)
      M.SetBuffer(i, j, X(j));
    M.Flush();
  }


  //! Sets a row of a matrix.
  /*!
    \param[in] X new row \a i of \a M.
    \param[in] i row index.
    \param[out] M matrix.
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
	      int i, Matrix<T0, General, RowSparse, Allocator0>& M)
  {
    int m = M.GetM();
    int n = M.GetN();
    int nnz = M.GetDataSize();
    int Nx = X.GetSize();

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= m)
      throw WrongIndex("SetRow(Vector, int, Matrix<RowSparse>)",
                       string("Index should be in [0, ") + to_str(m - 1)
                       + "], but is equal to " + to_str(i) + ".");
#endif

    int *ptr_vector =  M.GetPtr();
    int ptr_i0 =  ptr_vector[i], ptr_i1 = ptr_vector[i + 1];
    int row_size_difference = Nx - ptr_i1 + ptr_i0;

    if (row_size_difference == 0)
      {
	for (int k = 0; k < Nx; k++)
	  M.GetInd()[k + ptr_i0] = X.Index(k);
	for (int k = 0; k < Nx; k++)
	  M.GetData()[k + ptr_i0] = X.Value(k);
	return;
      }

    Vector<int, VectFull, CallocAlloc<int> >
      new_ind_vector(nnz + row_size_difference);
    for (int k = 0; k <  ptr_i0; k++)
      new_ind_vector(k) = M.GetInd()[k];
    for (int k = 0; k < Nx; k++)
      new_ind_vector(k + ptr_i0) = X.Index(k);
    for (int k = 0; k < nnz - ptr_i1; k++)
      new_ind_vector(k + ptr_i0 + Nx) =  M.GetInd()[k + ptr_i1];

    Vector<T1, VectFull, Allocator0 >
      new_data_vector(nnz + row_size_difference);
    for (int k = 0; k <  ptr_i0; k++)
      new_data_vector(k) = M.GetData()[k];
    for (int k = 0; k < Nx; k++)
      new_data_vector(k + ptr_i0) = X.Value(k);
    for (int k = 0; k < nnz - ptr_i1; k++)
      new_data_vector(k + ptr_i0 + Nx) =  M.GetData()[k + ptr_i1];

    Vector<int, VectFull, CallocAlloc<int> > new_ptr_vector(m + 1);
    for (int j = 0; j < i + 1; j++)
      new_ptr_vector(j) = ptr_vector[j];
    for (int j = i + 1; j < m+1; j++)
      new_ptr_vector(j) =  ptr_vector[j] + row_size_difference;

    M.SetData(m, n, new_data_vector, new_ptr_vector, new_ind_vector);
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, Storage0, Allocator0>& M)
  {
    for (int i = 0; i < M.GetM(); i++)
      M.Set(i, j, X(i));
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScSeq, Allocator1>& X,
              int j, Matrix<T0, Prop0, PETScSeqDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Prop0, class Allocator0,
            class T1, class Allocator1>
  void SetCol(const Vector<T1, PETScPar, Allocator1>& X,
              int j, Matrix<T0, Prop0, PETScMPIDense, Allocator0>& M)
  {
    int a, b;
    M.GetProcessorRowRange(a, b);
    for (int i = a; i < b; i++)
      M.SetBuffer(i, j, X(i));
    M.Flush();
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Allocator0,
	    class T1, class Allocator1>
  void SetCol(const Vector<T1, VectFull, Allocator1>& X,
	      int j, Matrix<T0, General, RowSparse, Allocator0>& M)
  {
    Vector<T1, VectSparse, Allocator1> X_sparse;
    for (int k = 0; k < X.GetLength(); k++)
      {
        T1 value = X(k);
        if (value != T1(0.))
          X_sparse.AddInteraction(k, value);
      }

    SetCol(X_sparse, j, M);
  }


  //! Sets a column of a matrix.
  /*!
    \param[in] X new column \a j of \a M.
    \param[in] j column index.
    \param[out] M matrix.
  */
  template <class T0, class Allocator0, class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
	      int j, Matrix<T0, General, RowSparse, Allocator0>& M)
  {
    int m = M.GetM();
    int n = M.GetN();
    int nnz = M.GetDataSize();
    int Nx = X.GetSize();

#ifdef SELDON_CHECK_BOUNDS
    if (j < 0 || j >= n)
      throw WrongIndex("SetCol(Vector, int, Matrix<RowSparse>)",
                       string("Index should be in [0, ") + to_str(n - 1)
                       + "], but is equal to " + to_str(j) + ".");
#endif

    // The column to be changed.
    Vector<T1, VectSparse, Allocator1> column_j;
    GetCol(M, j, column_j);
    int Ncolumn_j = column_j.GetSize();
    int column_size_difference = Nx - Ncolumn_j;

    // Built a vector indexed with the rows of column_j and X.
    Vector<int, VectSparse> column_j_mask;
    Vector<int> index_j(Ncolumn_j);
    Vector<int> value_j(Ncolumn_j);
    for (int p = 0; p < Ncolumn_j; p++)
      index_j(p) = column_j.Index(p);
    value_j.Fill(-1);
    column_j_mask.SetData(value_j, index_j);
    value_j.Nullify();
    index_j.Nullify();
    Vector<int, VectSparse> X_mask;
    Vector<int> index_x(Nx);
    Vector<int> value_x(Nx);
    for (int p = 0; p < Nx; p++)
      index_x(p) = X.Index(p);
    value_x.Fill(1);
    X_mask.SetData(value_x, index_x);
    value_x.Nullify();
    index_x.Nullify();
    X_mask.AddInteractionRow(column_j_mask.GetSize(),
			     column_j_mask.GetIndex(),
			     column_j_mask.GetData(), true);

    // Built the new pointer vector.
    Vector<int, VectFull, CallocAlloc<int> > ptr_vector;
    ptr_vector.SetData(m + 1, M.GetPtr());
    Vector<int, VectFull, CallocAlloc<int> > new_ptr_vector(m + 1);
    new_ptr_vector.Zero();
    for (int p = 0; p < X_mask.GetSize(); p++)
      new_ptr_vector(X_mask.Index(p) + 1) = X_mask.Value(p);
    for (int p = 0; p < m; p++)
      new_ptr_vector(p + 1) += new_ptr_vector(p);

    Add(1, ptr_vector, new_ptr_vector);

    // Built the new index and the new data vectors row by row.
    Vector<int, VectFull, CallocAlloc<int> >
      new_ind_vector(nnz + column_size_difference);
    Vector<T0, VectFull, Allocator0>
      new_data_vector(nnz + column_size_difference);

    Vector<T0, VectSparse, Allocator0> working_vector;
    int Nworking_vector;

    int line = 0;
    for (int interaction = 0; interaction < X_mask.GetSize(); interaction++)
      {
	int ind_x =  X_mask.Index(interaction);
	for (int k = 0; k < ptr_vector(ind_x) -  ptr_vector(line); k++)
	  new_ind_vector.GetData()[k + new_ptr_vector(line)] =
	    M.GetInd()[k + ptr_vector(line)];
	for (int k = 0; k < ptr_vector(ind_x) -  ptr_vector(line); k++)
	  new_data_vector.GetData()[k + new_ptr_vector(line)] =
	    M.GetData()[k + ptr_vector(line)];

	int ind_j;
	Nworking_vector = ptr_vector(ind_x + 1) - ptr_vector(ind_x);
	working_vector.SetData(Nworking_vector,
			       M.GetData() + ptr_vector(ind_x),
			       M.GetInd() + ptr_vector(ind_x));
	switch(X_mask.Value(interaction))
	  {
	    // Collision.
	  case 0:
	    working_vector.Get(j) = X(ind_x);
	    for (int k = 0; k < Nworking_vector; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetIndex()[k];
	    for (int k = 0; k < Nworking_vector; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetData()[k];
	    break;

	    // Suppression.
	  case -1:
	    ind_j = 0;
	    while (ind_j < Nworking_vector &&
		   working_vector.Index(ind_j) != j)
	      ind_j++;

	    for (int k = 0; k < ind_j; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetIndex()[k];
	    for (int k = 0; k < Nworking_vector - ind_j - 1; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j] =
		working_vector.GetIndex()[k + ind_j + 1];

	    for (int k = 0; k < ind_j; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetData()[k];
	    for (int k = 0; k < Nworking_vector - ind_j - 1; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j] =
		working_vector.GetData()[k + ind_j + 1];
	    break;

	    // Addition.
	  case 1:
	    ind_j = 0;
	    while (ind_j < Nworking_vector &&
		   working_vector.Index(ind_j) < j)
	      ind_j++;
	    for (int k = 0; k < ind_j; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetIndex()[k];
	    new_ind_vector.GetData()[new_ptr_vector(ind_x) + ind_j] = j;
	    for (int k = 0; k < Nworking_vector - ind_j; k++)
	      new_ind_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j + 1]
		= working_vector.GetIndex()[k + ind_j];

	    for (int k = 0; k < ind_j; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x)] =
		working_vector.GetData()[k];
	    new_data_vector.GetData()[new_ptr_vector(ind_x)  + ind_j]
	      = X(ind_x);
	    for (int k = 0; k < Nworking_vector - ind_j; k++)
	      new_data_vector.GetData()[k + new_ptr_vector(ind_x) + ind_j + 1]
		= working_vector.GetData()[k + ind_j];
	  }

	line = ind_x + 1;
	working_vector.Nullify();
      }
    for (int k = 0; k < ptr_vector(m) -  ptr_vector(line); k++)
      new_ind_vector.GetData()[k + new_ptr_vector(line)] =
	M.GetInd()[k + ptr_vector(line)];
    for (int k = 0; k < ptr_vector(m) -  ptr_vector(line); k++)
      new_data_vector.GetData()[k + new_ptr_vector(line)] =
	M.GetData()[k + ptr_vector(line)];

    M.SetData(m, n, new_data_vector, new_ptr_vector, new_ind_vector);
    ptr_vector.Nullify();
    new_data_vector.Nullify();
    new_ind_vector.Nullify();
    new_ptr_vector.Nullify();
  }


  //! Permutation of a general matrix stored by rows.
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, RowMajor, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, RowMajor, Allocator> A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(i, j) = A_copy(row_perm(i) - starting_index,
                         col_perm(j) - starting_index);
  }


  //! Permutation of a general matrix stored by columns.
  /*!
    B(i, j) = A(row_perm(i), col_perm(j)) and A = B.
  */
  template<class T, class Prop, class Allocator>
  void ApplyPermutation(Matrix<T, Prop, ColMajor, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm,
                        int starting_index)
  {
    Matrix<T, Prop, ColMajor, Allocator> A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i < A.GetM(); i++)
        A(i, j) = A_copy(row_perm(i) - starting_index,
                         col_perm(j) - starting_index);
  }


  //! Inverse permutation of a general matrix stored by rows.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, RowMajor, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, RowMajor, Allocator> A_copy = A;

    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        A(row_perm(i) - starting_index, col_perm(j) - starting_index)
          = A_copy(i, j);
  }


  //! Inverse permutation of a general matrix stored by columns.
  /*!
    B(row_perm(i), col_perm(j)) = A(i,j) and A = B.

    Equivalent Matlab operation: A(row_perm, col_perm) = A.
  */
  template<class T, class Prop, class Allocator>
  void ApplyInversePermutation(Matrix<T, Prop, ColMajor, Allocator>& A,
                               const Vector<int>& row_perm,
                               const Vector<int>& col_perm,
                               int starting_index)
  {
    Matrix<T, Prop, ColMajor, Allocator> A_copy = A;

    for (int j = 0; j < A.GetN(); j++)
      for (int i = 0; i < A.GetM(); i++)
        A(row_perm(i) - starting_index, col_perm(j) - starting_index)
          = A_copy(i, j);
  }


} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_CXX
#endif
