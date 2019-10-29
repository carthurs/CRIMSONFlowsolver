// Copyright (C) 2001-2009 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_CXX
#define SELDON_FILE_FUNCTIONS_MATRIX_CXX


#include "Functions_Matrix.hxx"


/*
  Function defined in this file:

  alpha A -> A
  Mlt(alpha, A)

  A B -> C
  Mlt(A, B, C)

  alpha A B -> C
  Mlt(alpha, A, B, C)

  alpha A B + beta C -> C
  MltAdd(alpha, A, B, beta, C)

  alpha A + B -> B
  Add(alpha, A, B)

  LU factorization of matrix A without pivoting.
  GetLU(A)

  Highest absolute value of A.
  MaxAbs(A)

  1-norm of matrix A.
  Norm1(A)

  infinity norm of matrix A.
  NormInf(A)

  Transpose(A)
*/


namespace Seldon
{


  /////////
  // MLT //


  //! Multiplies a matrix by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void Mlt(const T0 alpha, Matrix<T1, Prop1, Storage1, Allocator1>& A)
  {
    T1 alpha_ = alpha;

    typename Matrix<T1, Prop1, Storage1, Allocator1>::pointer
      data = A.GetData();

    for (int i = 0; i < A.GetDataSize(); i++)
      data[i] = alpha_ * data[i];
  }


  //! Multiplies a matrix by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1>
  void Mlt(const T0 alpha,
	   Matrix<T1, Prop1, ColMajorCollection, Allocator1>& A)
  {
    typename T1::value_type alpha_ = alpha;
    for (int i = 0; i < A.GetMmatrix(); i++)
      for (int j = 0; j < A.GetNmatrix(); j++)
	Mlt(alpha, A.GetMatrix(i, j));
  }


  //! Multiplies a matrix by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1>
  void Mlt(const T0 alpha,
	   Matrix<T1, Prop1, RowMajorCollection, Allocator1>& A)
  {
    typename T1::value_type alpha_ = alpha;
    for (int i = 0; i < A.GetMmatrix(); i++)
      for (int j = 0; j < A.GetNmatrix(); j++)
	Mlt(alpha_, A.GetMatrix(i, j));
  }


  //! Multiplies a FloatDouble collection by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T0, class Allocator>
  void Mlt(const T0 alpha,
	   Matrix<FloatDouble, General, DenseSparseCollection, Allocator>& A)
  {
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator>
      ::float_dense_m m0;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator>
      ::float_sparse_m m1;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator>
      ::double_dense_m m2;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator>
      ::double_sparse_m m3;

    for (int i = 0; i < A.GetMmatrix(); i++)
      for (int j = 0; j < A.GetNmatrix(); j++)
	{
	  switch (A.GetType(i, j))
	    {
	    case 0:
	      A.GetMatrix(i, j, m0);
	      Mlt(float(alpha), m0);
	      A.SetMatrix(i, j, m0);
	      m0.Nullify();
	      break;
	    case 1:
	      A.GetMatrix(i, j, m1);
	      Mlt(float(alpha), m1);
	      A.SetMatrix(i, j, m1);
	      m1.Nullify();
	      break;
	    case 2:
	      A.GetMatrix(i, j, m2);
	      Mlt(double(alpha), m2);
	      A.SetMatrix(i, j, m2);
	      m2.Nullify();
	      break;
	    case 3:
	      A.GetMatrix(i, j, m3);
	      Mlt(double(alpha), m3);
	      A.SetMatrix(i, j, m3);
	      m3.Nullify();
	      break;
	    default:
	      throw WrongArgument("Mlt(alpha, Matrix<FloatDouble, "
				  "DenseSparseCollection)",
				  "Underlying matrix (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}
  }


  //! Multiplies two matrices.
  /*! It performs the operation \f$ C = \alpha A B \f$ where \f$ \alpha \f$ is
    a scalar, and \f$ A \f$, \f$ B \f$ and \f$ C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[out] C matrix, result of the product of \a A with \a B, times \a
    alpha.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3, class Prop3, class Storage3, class Allocator3>
  void Mlt(const T0 alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	   const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	   Matrix<T3, Prop3, Storage3, Allocator3>& C)
  {
    C.Fill(T3(0));
    MltAdd(alpha, A, B, T3(0), C);
  }


  //! Multiplies two matrices.
  /*! It performs the operation \f$ C = A B \f$ where \f$ A \f$, \f$ B \f$ and
    \f$ C \f$ are matrices.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[out] C matrix, result of the product of \a A with \a B.
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& B,
	   Matrix<T2, Prop2, Storage2, Allocator2>& C)
  {
    C.Fill(T2(0));
    MltAdd(T0(1), A, B, T2(0), C);
  }


  //! Multiplies two row-major sparse matrices in Harwell-Boeing format.
  /*! It performs the operation \f$ C = A B \f$ where \f$ A \f$, \f$ B \f$ and
    \f$ C \f$ are row-major sparse matrices in Harwell-Boeing format.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[out] C row-major sparse matrix in Harwell-Boeing format, result of
    the product of \a A with \a B. It does not need to have the right non-zero
    entries.
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	   const Matrix<T1, Prop1, RowSparse, Allocator1>& B,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, "Mlt(const Matrix<RowSparse>& A, const "
             "Matrix<RowSparse>& B, Matrix<RowSparse>& C)");
#endif

    int h, i, k, l, col;
    int Nnonzero, Nnonzero_row, Nnonzero_row_max;
    IVect column_index;
    Vector<T2> row_value;
    T1 value;
    int m = A.GetM();

    int* c_ptr = NULL;
    int* c_ind = NULL;
    T2* c_data = NULL;
    C.Clear();

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	c_ptr = reinterpret_cast<int*>(calloc(m + 1, sizeof(int)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        c_ptr = NULL;
      }

    if (c_ptr == NULL)
      throw NoMemory("Mlt(const Matrix<RowSparse>& A, const "
                     "Matrix<RowSparse>& B, Matrix<RowSparse>& C)",
		     "Unable to allocate memory for an array of "
		     + to_str(m + 1) + " integers.");
#endif

    c_ptr[0] = 0;

    // Number of non-zero elements in C.
    Nnonzero = 0;
    for (i = 0; i < m; i++)
      {
        c_ptr[i + 1] = c_ptr[i];

        if (A.GetPtr()[i + 1] != A.GetPtr()[i])
          // There are elements in the i-th row of A, so there can be non-zero
          // entries in C as well. Checks whether any column in B has an
          // element whose row index matches a column index of a non-zero in
          // the i-th row of A.
          {
            // Maximum number of non-zero entry on the i-th row of C.
            Nnonzero_row_max = 0;
            // For every element in the i-th row.
            for (k = A.GetPtr()[i]; k < A.GetPtr()[i + 1]; k++)
              {
                col = A.GetInd()[k];
                Nnonzero_row_max += B.GetPtr()[col + 1] - B.GetPtr()[col];
              }
            // Now gets the column indexes.
            column_index.Reallocate(Nnonzero_row_max);
            row_value.Reallocate(Nnonzero_row_max);
            h = 0;
            // For every element in the i-th row.
            for (k = A.GetPtr()[i]; k < A.GetPtr()[i + 1]; k++)
              {
                // The k-th column index (among the nonzero entries) on the
                // i-th row, and the corresponding value.
                col = A.GetInd()[k];
                value = A.GetData()[k];
                // Loop on all elements in the col-th row in B. These elements
                // are multiplied with the element (i, col) of A.
                for (l = B.GetPtr()[col]; l < B.GetPtr()[col + 1]; l++)
                  {
                    column_index(h) = B.GetInd()[l];
                    row_value(h) = value * B.GetData()[l];
                    h++;
                  }
              }
            // Now gathers and sorts all elements on the i-th row of C.
            Nnonzero_row = column_index.GetLength();
            Assemble(Nnonzero_row, column_index, row_value);

#ifdef SELDON_CHECK_MEMORY
            try
              {
#endif

                // Reallocates 'c_ind' and 'c_data' in order to append the
                // elements of the i-th row of C.
                c_ind = reinterpret_cast<int*>
                  (realloc(reinterpret_cast<void*>(c_ind),
                           (Nnonzero + Nnonzero_row) * sizeof(int)));
                c_data = reinterpret_cast<T2*>
                  (C.GetAllocator().reallocate(c_data,
                                               Nnonzero + Nnonzero_row));

#ifdef SELDON_CHECK_MEMORY
              }
            catch (...)
              {
                c_ind = NULL;
                c_data = NULL;
              }

            if ((c_ind == NULL || c_data == NULL)
                && Nnonzero + Nnonzero_row != 0)
              throw NoMemory("Mlt(const Matrix<RowSparse>& A, const "
                             "Matrix<RowSparse>& B, Matrix<RowSparse>& C)",
                             "Unable to allocate memory for an array of "
                             + to_str(Nnonzero + Nnonzero_row) + " integers "
                             "and for an array of "
                             + to_str(sizeof(T2) * (Nnonzero + Nnonzero_row))
                             + " bytes.");
#endif

            c_ptr[i + 1] += Nnonzero_row;
            for (h = 0; h < Nnonzero_row; h++)
              {
                c_ind[Nnonzero + h] = column_index(h);
                c_data[Nnonzero + h] = row_value(h);
              }
            Nnonzero += Nnonzero_row;
          }
      }

    C.SetData(A.GetM(), B.GetN(), Nnonzero, c_data, c_ptr, c_ind);
  }


  //! Multiplies two row-major sparse matrices in Harwell-Boeing format.
  /*! It performs the operation \f$ C = A B \f$ where \f$ A \f$, \f$ B \f$ and
    \f$ C \f$ are row-major sparse matrices in Harwell-Boeing format.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[out] C row-major sparse matrix in Harwell-Boeing format, result of
    the product of \a A with \a B. It does not need to have the right non-zero
    entries.
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, RowMajor, Allocator0>& A,
	   const Matrix<T1, Prop1, RowSparse, Allocator1>& B,
	   Matrix<T2, Prop2, RowMajor, Allocator2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, "Mlt(const Matrix<RowMajor>& A, const "
             "Matrix<RowSparse>& B, Matrix<RowMajor>& C)");
#endif

    int m = A.GetM();
    int n = A.GetN();

    C.Reallocate(A.GetM(), B.GetN());
    C.Fill(T2(0));

    for (int i = 0; i < m; i++)
      {
        for (int k = 0; k < n; k++)
          {
            // Loop on all elements in the k-th row in B. These elements
            // are multiplied with the element (i, k) of A.
            for (int l = B.GetPtr()[k]; l < B.GetPtr()[k + 1]; l++)
              C(i, B.GetInd()[l]) += A(i, k) * B.GetData()[l];
          }
      }
  }


  //! Multiplies two row-major sparse matrices in Harwell-Boeing format.
  /*! It performs the operation \f$ C = A B \f$ where \f$ A \f$, \f$ B \f$ and
    \f$ C \f$ are row-major sparse matrices in Harwell-Boeing format.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[out] C row-major sparse matrix in Harwell-Boeing format,result of
    the product of \a A with \a B transposed. It does not need to have the
    right non-zero entries.
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2>
  void MltNoTransTrans(const Matrix<T0, Prop0, RowMajor, Allocator0>& A,
                       const Matrix<T1, Prop1, RowSparse, Allocator1>& B,
                       Matrix<T2, Prop2, RowMajor, Allocator2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, SeldonTrans, B,
             "MltNoTransTrans(const Matrix<RowMajor>& A, const "
             "Matrix<RowSparse>& B, Matrix<RowMajor>& C)");
#endif

    int m = A.GetM();
    C.Reallocate(A.GetM(), B.GetM());
    C.Fill(T2(0));
    for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < B.GetM(); j++)
          {
            // Loop on all elements in the i-th row in B. These elements
            // are multiplied with the element (i, k) of A.
            for (int l = B.GetPtr()[j]; l < B.GetPtr()[j + 1]; l++)
              C(i, j) += A(i, B.GetInd()[l]) * B.GetData()[l];
          }
      }
  }


  //! Multiplies two row-major sparse matrices in Harwell-Boeing format.
  /*! It performs the operation \f$ C = A B^T \f$ where \f$ A \f$, \f$ B \f$
    and \f$ C \f$ are row-major sparse matrices in Harwell-Boeing format.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[out] C row-major sparse matrix in Harwell-Boeing format, result of
    the product of \a A with \a B transposed. It does not need to have the
    right non-zero entries.
  */
  template <class T0, class Prop0, class Allocator0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2>
  void MltNoTransTrans(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
                       const Matrix<T1, Prop1, RowSparse, Allocator1>& B,
                       Matrix<T2, Prop2, RowSparse, Allocator2>& C)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, SeldonTrans, B,
             "MltNoTransTrans(const Matrix<RowSparse>& A, "
             "const Matrix<RowSparse>& B, Matrix<RowSparse>& C)");
#endif

    int h, i, k, col;
    int ib, kb;
    int Nnonzero_row;
    int Nnonzero;

    // 'MallocAlloc' is specified so that reallocations may be efficient.
    // There will be no need for 'Resize': 'Reallocate' will do the job.
    Vector<int, VectFull, MallocAlloc<int> > column_index;
    Vector<T2, VectFull, MallocAlloc<T2> > row_value;
    T2 value = 0;

    int m = A.GetM();
    int n = B.GetM();

    int* c_ptr = NULL;
    int* c_ind = NULL;
    T2* c_data = NULL;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	c_ptr = reinterpret_cast<int*>(calloc(m + 1, sizeof(int)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        c_ptr = NULL;
      }

    if (c_ptr == NULL)
      throw NoMemory("MltNoTransTrans(const Matrix<RowSparse>& A, "
                     "const Matrix<RowSparse>& B, Matrix<RowSparse>& C)",
		     "Unable to allocate memory for an array of "
		     + to_str(m + 1) + " integers.");
#endif

    c_ptr[0] = 0;

    // Number of non-zero elements in C.
    Nnonzero = 0;

    for (i = 0; i < m; i++)
      {
        c_ptr[i + 1] = c_ptr[i];

        if (A.GetPtr()[i + 1] != A.GetPtr()[i])
          // There are elements in the i-th row of A, so there can be non-zero
          // entries in C as well. It is checked below whether any row in B
          // has an element whose row index matches a column index of a
          // non-zero in the i-th row of A.
          {
            // For every element in the i-th row.
            for (k = A.GetPtr()[i]; k < A.GetPtr()[i + 1]; k++)
              {
                col = A.GetInd()[k];
                // For every row in B.
                for (ib = 0; ib < n; ib++)
                  {
                    for (kb = B.GetPtr()[ib]; kb < B.GetPtr()[ib + 1]; kb++)
                      if (col == B.GetInd()[kb])
                        value += A.GetData()[k] * B.GetData()[kb];
                    if (value != T2(0))
                      {
                        row_value.Append(value);
                        column_index.Append(ib);
                        value = T2(0);
                      }
                  }
              }

            Nnonzero_row = column_index.GetLength();
            Assemble(Nnonzero_row, column_index, row_value);

#ifdef SELDON_CHECK_MEMORY
            try
              {
#endif

                // Reallocates 'c_ind' and 'c_data' in order to append the
                // elements of the i-th row of C.
                c_ind = reinterpret_cast<int*>
                  (realloc(reinterpret_cast<void*>(c_ind),
                           (Nnonzero + Nnonzero_row) * sizeof(int)));
                c_data = reinterpret_cast<T2*>
                  (C.GetAllocator().reallocate(c_data,
                                               Nnonzero + Nnonzero_row));

#ifdef SELDON_CHECK_MEMORY
              }
            catch (...)
              {
                c_ind = NULL;
                c_data = NULL;
              }

            if ((c_ind == NULL || c_data == NULL)
                && Nnonzero_row != 0)
              throw NoMemory("MltNoTransTrans(const Matrix<RowSparse>& A, "
                             "const Matrix<RowSparse>& B, "
                             "Matrix<RowSparse>& C)",
                             "Unable to allocate memory for an array of "
                             + to_str(Nnonzero + Nnonzero_row) + " integers "
                             "and for an array of "
                             + to_str(sizeof(T2) * (Nnonzero + Nnonzero_row))
                             + " bytes.");
#endif

            c_ptr[i + 1] += Nnonzero_row;
            for (h = 0; h < Nnonzero_row; h++)
              {
                c_ind[Nnonzero + h] = column_index(h);
                c_data[Nnonzero + h] = row_value(h);
              }
            Nnonzero += Nnonzero_row;
          }

        column_index.Clear();
        row_value.Clear();
      }

    C.SetData(A.GetM(), B.GetM(), Nnonzero, c_data, c_ptr, c_ind);
  }


  // MLT //
  /////////


  ////////////
  // MLTADD //


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	      const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    int na = A.GetN();
    int mc = C.GetM();
    int nc = C.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    T4 temp;
    T4 alpha_(alpha);
    T4 beta_(beta);

    if (beta_ != T4(0))
      for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
	for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
	  C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j))
	    *= beta_;
    else
      for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
	for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
	  C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j)) = T4(0);

    for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
      for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
	{
	  temp = T4(0);
	  for (int k = 0; k < na; k++)
	    temp += A(Storage4::GetFirst(i, j), k)
	      * B(k, Storage4::GetSecond(i, j));
	  C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j))
	    += alpha_ * temp;
	}
  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$
    or \f$C = \alpha A B^T + \beta C \f$, where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] TransA status of A: it must be SeldonNoTrans. This argument
    is required for consistency with the interface for full matrices.
    \param[in] A matrix.
    \param[in] TransB status of B: SeldonNoTrans or SeldonTrans.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B or \a
    B^T, times \a alpha, plus \a beta times \a C.
  */
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	      const SeldonTranspose& TransB,
	      const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    if (!TransA.NoTrans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransA' must be equal to 'SeldonNoTrans'.");
    if (!TransB.NoTrans() && !TransB.Trans())
      throw WrongArgument("MltAdd(const T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<T1, Prop1, Storage1, Allocator1>& A"
                          "SeldonTranspose TransB,"
                          "const Matrix<T2, Prop2, Storage2, Allocator2>& B,"
                          "const T3 beta,"
                          "Matrix<T4, Prop4, Storage4, Allocator4>& C)",
                          "'TransB' must be equal to 'SeldonNoTrans' or "
                          "'SeldonTrans'.");

    if (!TransB.Trans())
      MltAdd(alpha, A, B, beta, C);
    else
      {

        int na = A.GetN();
        int mc = C.GetM();
        int nc = C.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
        CheckDim(SeldonNoTrans, A, SeldonTrans, B, C,
                 "MltAdd(alpha, TransA, A, TransB, B, beta, C)");
#endif

        T4 temp;
        T4 alpha_(alpha);
        T4 beta_(beta);

        if (beta_ != T4(0))
          for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
            for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
              C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j))
                *= beta_;
        else
          for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
            for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
              C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j)) = T4(0);

        for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
          for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
            {
              temp = T4(0);
              for (int k = 0; k < na; k++)
                temp += A(Storage4::GetFirst(i, j), k)
                  * B(Storage4::GetSecond(i, j), k);
              C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j))
                += alpha_ * temp;
            }
      }
  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& A,
              const Matrix<T2, General, RowMajor, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, PETScMPIDense, Allocator4>& C)
  {
    int na = A.GetN();
    int mc = C.GetM();
    int nc = C.GetN();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif
    T1 *local_a;
    MatGetArray(A.GetPetscMatrix(), &local_a);
    int nlocal_A;
    int mlocal_A;
    MatGetLocalSize(A.GetPetscMatrix(), &mlocal_A, &nlocal_A);
    Matrix<T1, Prop1, ColMajor, Allocator1> local_A;
    local_A.SetData(mlocal_A, na, local_a);

    T4 *local_c;
    MatGetArray(C.GetPetscMatrix(), &local_c);
    int nlocal_C;
    int mlocal_C;
    MatGetLocalSize(C.GetPetscMatrix(), &mlocal_C, &nlocal_C);
    Matrix<T4, Prop4, ColMajor, Allocator4> local_C;
    local_C.SetData(mlocal_C, nc, local_c);

    MltAdd(alpha, local_A, B, beta, local_C);

    local_A.Nullify();
    MatRestoreArray(A.GetPetscMatrix(), &local_a);
    A.Flush();

    local_C.Nullify();
    MatRestoreArray(C.GetPetscMatrix(), &local_c);
    C.Flush();
  }


   //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& A,
	      const Matrix<T2, Prop2, RowMajorCollection, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, RowMajorCollection, Allocator4>& C)
  {
    int na = A.GetNmatrix();
    int mc = C.GetMmatrix();
    int nc = C.GetNmatrix();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    typedef typename T4::value_type value_type;

    Mlt(value_type(beta), C);

    for (int i = 0; i < mc; i++ )
      for (int j = 0; j < nc; j++)
	for (int k = 0; k < na; k++)
	  MltAdd(value_type(alpha), A.GetMatrix(i, k), B.GetMatrix(k, j),
		 value_type(1), C.GetMatrix(i, j));
  }


  //! Multiplies two matrices, and adds the result to a third matrix.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$
    \alpha \f$ and \f$ \beta \f$ are scalars, and \f$ A \f$, \f$ B \f$ and \f$
    C \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in] B matrix.
    \param[in] beta scalar.
    \param[in,out] C matrix, result of the product of \a A with \a B, times \a
    alpha, plus \a beta times \a C.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& A,
	      const Matrix<T2, Prop2, ColMajorCollection, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, ColMajorCollection, Allocator4>& C)
  {
    int na = A.GetNmatrix();
    int mc = C.GetMmatrix();
    int nc = C.GetNmatrix();

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    typedef typename T4::value_type value_type;

    Mlt(value_type(beta), C);

    for (int i = 0; i < mc; i++ )
      for (int j = 0; j < nc; j++)
	for (int k = 0; k < na; k++)
	  MltAdd(value_type(alpha), A.GetMatrix(i, k), B.GetMatrix(k, j),
		 value_type(1), C.GetMatrix(i, j));
  }


  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, General, RowMajor, Allocator1>& A,
	      const Matrix<T2, General, RowMajor, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, General, RowSparse, Allocator4>& C)
  {
    throw Undefined("void MltAdd(const T0 alpha,"
		    "const Matrix<T1, General, RowMajor, Allocator1>& A,"
		    "const Matrix<T2, General, RowMajor, Allocator2>& B,"
		    "const T3 beta,"
		    "Matrix<T4, General, RowSparse, Allocator4>& C)");
  }


  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, General, RowMajor, Allocator1>& A,
	      const Matrix<T2, General, RowSparse, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, General, RowSparse, Allocator4>& C)
  {
    throw Undefined("void MltAdd(const T0 alpha,"
		    "const Matrix<T1, General, RowMajor, Allocator1>& A,"
		    "const Matrix<T2, General, RowSparse, Allocator2>& B,"
		    "const T3 beta,"
		    "Matrix<T4, General, RowSparse, Allocator4>& C)");
  }


  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, General, RowSparse, Allocator1>& A,
	      const Matrix<T2, General, RowMajor, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, General, RowSparse, Allocator4>& C)
  {
    throw Undefined("void MltAdd(const T0 alpha,"
		    "const Matrix<T1, General, RowSparse, Allocator1>& A,"
		    "const Matrix<T2, General, RowMajor, Allocator2>& B,"
		    "const T3 beta,"
		    "Matrix<T4, General, RowSparse, Allocator4>& C)");
  }


  template <class T0,
            class Allocator1,
            class Allocator2,
            class Allocator3,
            class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd_heterogeneous(const T0 alpha,
			    const Matrix<FloatDouble, General,
			    DenseSparseCollection, Allocator1>& A,
			    const Matrix<FloatDouble, General,
			    DenseSparseCollection, Allocator2>& B,
			    Matrix<FloatDouble, General,
			    DenseSparseCollection, Allocator3>& C,
			    Matrix<T4, Prop4, Storage4, Allocator4>& mc,
			    int i, int j)
  {
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator1>
      ::float_dense_m m0a;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator1>
      ::float_sparse_m m1a;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator1>
      ::double_dense_m m2a;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator1>
      ::double_sparse_m m3a;

    int na = A.GetNmatrix();
    for (int k = 0; k < na; k++)
      {
	switch (A.GetType(i, k))
	  {
	  case 0:
	    A.GetMatrix(i, k, m0a);
	    MltAdd_heterogeneous2(alpha, m0a, B, C, mc, j, k);
	    m0a.Nullify();
	    break;
	  case 1:
	    A.GetMatrix(i, k, m1a);
	    MltAdd_heterogeneous2(alpha, m1a, B, C, mc, j, k);
	    m1a.Nullify();
	    break;
	  case 2:
	    A.GetMatrix(i, k, m2a);
	    MltAdd_heterogeneous2(alpha, m2a, B, C, mc, j, k);
	    m2a.Nullify();
	    break;
	  case 3:
	    A.GetMatrix(i, k, m3a);
	    MltAdd_heterogeneous2(alpha, m3a, B, C, mc, j, k);
	    m3a.Nullify();
	    break;
	  default:
	    throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>::"
				"MltAdd_heterogeneous(alpha, A, B, beta, C) ",
				"Underlying matrix  A (" + to_str(i) + " ,"
				+ to_str(k) + " ) not defined.");
	  }
      }
  }


  template<class T0,
	   class T1, class Prop1, class Storage1, class Allocator1,
	   class Allocator2,
	   class Allocator3,
	   class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd_heterogeneous2(const T0 alpha,
			     const Matrix<T1, Prop1,
                             Storage1, Allocator1>& ma,
			     const Matrix<FloatDouble, General,
			     DenseSparseCollection, Allocator2>& B,
			     Matrix<FloatDouble, General,
			     DenseSparseCollection, Allocator3>& C,
			     Matrix<T4, Prop4, Storage4, Allocator4>& mc,
			     int j, int k)
  {
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::float_dense_m m0b;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::float_sparse_m m1b;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::double_dense_m m2b;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::double_sparse_m m3b;

    switch (B.GetType(k, j))
      {
      case 0:
	B.GetMatrix(k, j, m0b);
	MltAdd(alpha, ma, m0b, 1., mc);
	m0b.Nullify();
	break;
      case 1:
	B.GetMatrix(k, j, m1b);
	MltAdd(alpha, ma, m1b, 1., mc);
	m1b.Nullify();
	break;
      case 2:
	B.GetMatrix(k, j, m2b);
	MltAdd(alpha, ma, m2b, 1., mc);
	m2b.Nullify();
	break;
      case 3:
	B.GetMatrix(k, j, m3b);
	MltAdd(alpha, ma, m3b, 1., mc);
	m3b.Nullify();
	break;
      default:
	throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
			    "::MltAdd_heterogeneous2(alpha, A, B, beta, C)",
			    "Underlying matrix  B (" + to_str(k) + " ,"
			    + to_str(j) + " ) not defined.");
      }
  }



  //! Multiplies a FloatDouble collection by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T0, class Allocator1, class Allocator2, class T3,
	    class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<FloatDouble, General, DenseSparseCollection,
	      Allocator1>& A,
	      const Matrix<FloatDouble, General, DenseSparseCollection,
	      Allocator2>& B,
	      const T3 beta,
	      Matrix<FloatDouble, General, DenseSparseCollection,
	      Allocator4>& C)
  {
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator4>
      ::float_dense_m m0c;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator4>
      ::float_sparse_m m1c;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator4>
      ::double_dense_m m2c;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator4>
      ::double_sparse_m m3c;

    Mlt(beta, C);

    int mc = C.GetMmatrix();
    int nc = C.GetNmatrix();
    for (int i = 0; i < mc; i++ )
      for (int j = 0; j < nc; j++)
	{
	  switch (C.GetType(i, j))
	    {
	    case 0:
	      C.GetMatrix(i, j, m0c);
	      MltAdd_heterogeneous(float(alpha), A, B, C, m0c, i, j);
	      C.SetMatrix(i, j, m0c);
	      m0c.Nullify();
	      break;
	    case 1:
	      C.GetMatrix(i, j, m1c);
	      MltAdd_heterogeneous(float(alpha), A, B, C, m1c, i, j);
	      C.SetMatrix(i, j, m1c);
	      m1c.Nullify();
	      break;
	    case 2:
	      C.GetMatrix(i, j, m2c);
	      MltAdd_heterogeneous(double(alpha), A, B, C, m2c, i, j);
	      C.SetMatrix(i, j, m2c);
	      m2c.Nullify();
	      break;
	    case 3:
	      C.GetMatrix(i, j, m3c);
	      MltAdd_heterogeneous(double(alpha), A, B, C, m3c, i, j);
	      C.SetMatrix(i, j, m3c);
	      m3c.Nullify();
	      break;
	    default:
	      throw WrongArgument("Matrix<FloatDouble, DenseSparseCollection>"
				  "::MltAdd(alpha, A, B, beta, C) ",
				  "Underlying matrix  C (" + to_str(i) + " ,"
				  + to_str(j) + " ) not defined.");
	    }
	}
  }


  //! Multiplies two row-major sparse matrices and adds the result to a third.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ where \f$ A
    \f$, \f$ B \f$ and \f$ C \f$ are row-major sparse matrices in
    Harwell-Boeing format, and \f$ \alpha \f$ and \f$ \beta \f$ are scalars.
    \param[in] alpha scalar.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[in] beta scalar.
    \param[in,out] C row-major sparse matrix in Harwell-Boeing format. On
    exit, it is equal to \f$ \alpha A B + \beta C \f$.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
              const Matrix<T2, Prop2, RowSparse, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, RowSparse, Allocator4>& C)
  {
    if (beta == T3(0))
      {
        if (alpha == T0(0))
          Mlt(alpha, C);
        else
          {
            Mlt(A, B, C);
            if (alpha != T0(1))
              Mlt(alpha, C);
          }
      }
    else
      {
        if (alpha == T0(0))
          Mlt(beta, C);
        else
          {
            Matrix<T4, Prop4, RowSparse, Allocator4> tmp;
            Mlt(A, B, tmp);
            if (beta != T0(1))
              Mlt(beta, C);
            Add(alpha, tmp, C);
          }
      }
  }


  //! Multiplies two row-major sparse matrices and adds the result to a third.
  /*! It performs the operation \f$ C = \alpha A B^T + \beta C \f$ where \f$ A
    \f$, \f$ B \f$ and \f$ C \f$ are row-major sparse matrices in
    Harwell-Boeing format, and \f$ \alpha \f$ and \f$ \beta \f$ are scalars.
    \param[in] alpha scalar.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[in] beta scalar.
    \param[in,out] C row-major sparse matrix in Harwell-Boeing format. On
    exit, it is equal to \f$ \alpha A B^T + \beta C \f$.
    \warning It is not recommended to call that function directly: it might be
    removed in next versions of Seldon. Use MltAdd instead.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
            class T4, class Prop4, class Allocator4>
  void MltNoTransTransAdd(const T0 alpha,
                          const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
                          const Matrix<T2, Prop2, RowSparse, Allocator2>& B,
                          const T3 beta,
                          Matrix<T4, Prop4, RowSparse, Allocator4>& C)
  {
    if (beta == T3(0))
      {
        if (alpha == T0(0))
          Mlt(alpha, C);
        else
          {
            MltNoTransTrans(A, B, C);
            if (alpha != T0(1))
              Mlt(alpha, C);
          }
      }
    else
      {
        if (alpha == T0(0))
          Mlt(beta, C);
        else
          {
            Matrix<T4, Prop4, RowSparse, Allocator4> tmp;
            MltNoTransTrans(A, B, tmp);
            if (beta != T0(1))
              Mlt(beta, C);
            Add(alpha, tmp, C);
          }
      }
  }


  //! Multiplies two row-major sparse matrices and adds the result to a third.
  /*! It performs the operation \f$ C = \alpha A B + \beta C \f$ or \f$ C =
    \alpha A B^T + \beta C \f$ where \f$ A \f$, \f$ B \f$ and \f$ C \f$ are
    row-major sparse matrices in Harwell-Boeing format, and \f$ \alpha \f$ and
    \f$ \beta \f$ are scalars.
    \param[in] alpha scalar.
    \param[in] TransA status of A: it must be SeldonNoTrans. This argument
    is required for consistency with the interface for full matrices.
    \param[in] A row-major sparse matrix in Harwell-Boeing format.
    \param[in] TransB status of B: SeldonNoTrans or SeldonTrans.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[in] beta scalar.
    \param[in,out] C row-major sparse matrix in Harwell-Boeing format. On
    exit, it is equal to \f$ \alpha A B + \beta C \f$ or \f$ \alpha A B^T +
    \beta C \f$.
    \note If \a TransA is not SeldonNoTrans, an exception will be thrown.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const SeldonTranspose& TransA,
              const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
	      const SeldonTranspose& TransB,
              const Matrix<T2, Prop2, RowSparse, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, RowSparse, Allocator4>& C)
  {
    if (!TransA.NoTrans())
      throw WrongArgument("MltAdd(T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<RowSparse>& A, SeldonTranspose "
                          "TransB, const Matrix<RowSparse>& B, T3 beta, "
                          "Matrix<RowSparse>& C)",
                          "'TransA' must be equal to 'SeldonNoTrans'.");
    if (!TransB.NoTrans() && !TransB.Trans())
      throw WrongArgument("MltAdd(T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<RowSparse>& A, SeldonTranspose "
                          "TransB, const Matrix<RowSparse>& B, T3 beta, "
                          "Matrix<RowSparse>& C)",
                          "'TransB' must be equal to 'SeldonNoTrans' or "
                          "'SeldonTrans'.");

    if (TransB.Trans())
      MltNoTransTransAdd(alpha, A, B, beta, C);
    else
      MltAdd(alpha, A, B, beta, C);
  }


  /*! \brief Multiplies a row-major matrix with a row-major sparse matrix
    and adds it to a third.*/
  /*! It performs the operation \f$ C = \alpha A B^T + \beta C \f$ where \f$ A
    \f$ and \f$ C \f$ are row-major dense matrices and \f$ B \f$ a row-major
    sparse matrix in Harwell-Boeing format.
    \param[in] A row-major dense matrix.
    \param[in] B row-major sparse matrix in Harwell-Boeing format.
    \param[out] C resulting row-major dense matrix.
  */
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
              const SeldonTranspose& TransA,
              const Matrix<T1, Prop1, RowMajor, Allocator1>& A,
              const SeldonTranspose& TransB,
              const Matrix<T2, Prop2, RowSparse, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, RowMajor, Allocator4>& C)
  {
    if (!TransA.NoTrans())
      throw WrongArgument("MltAdd(T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<RowMajor>& A, SeldonTranspose "
                          "TransB, const Matrix<RowSparse>& B, T3 beta, "
                          "Matrix<RowSparse>& C)",
                          "'TransA' must be equal to 'SeldonNoTrans'.");
    if (!TransB.NoTrans() && !TransB.Trans())
      throw WrongArgument("MltAdd(T0 alpha, SeldonTranspose TransA, "
                          "const Matrix<RowMajor>& A, SeldonTranspose "
                          "TransB, const Matrix<RowSparse>& B, T3 beta, "
                          "Matrix<RowMajor>& C)",
                          "'TransB' must be equal to 'SeldonNoTrans' or "
                          "'SeldonTrans'.");

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, SeldonTrans, B, C,
             "MltAdd(T0 alpha, const Matrix<RowMajor>& A, const "
             "Matrix<RowSparse>& B, T3 beta, Matrix<RowMajor>& C)");
#endif

    if (TransB.Trans())
      {
        int m = A.GetM();
        if (beta == 0)
          C.Fill(T2(0));
        else
          Mlt(beta, C);
        if (alpha != 0)
          for (int i = 0; i < m; i++)
            for (int j = 0; j < B.GetM(); j++)
              {
                // Loop on all elements in the i-th row in B. These elements
                // are multiplied with the element (i, k) of A.
                for (int l = B.GetPtr()[j]; l < B.GetPtr()[j + 1]; l++)
                  C(i, j) += alpha * A(i, B.GetInd()[l]) * B.GetData()[l];
              }
      }
    else
      MltAdd(alpha, A, B, beta, C);
  }


  // MLTADD //
  ////////////



  /////////
  // ADD //


  //! Adds two matrices.
  /*! It performs the operation \f$ B = \alpha A + B \f$ where \f$ \alpha \f$
    is a scalar, and \f$ A \f$ and \f$ B \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in,out] B matrix, result of the addition of \a B (on entry) and \a
    A times \a alpha.
  */
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
	   class T2, class Prop2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	   Matrix<T2, Prop2, Storage2, Allocator2>& B)
  {
    int i, j;
    for (i = 0; i < A.GetM(); i++)
      for (j = 0; j < A.GetN(); j++)
	B(i, j) += alpha * A(i, j);
  }


  //! Adds two matrices.
  /*! It performs the operation \f$ B = \alpha A + B \f$ where \f$ \alpha \f$
    is a scalar, and \f$ A \f$ and \f$ B \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in,out] B matrix, result of the addition of \a B (on entry) and \a
    A times \a alpha.
  */
  template<class T0,
	   class T1, class Prop1, class Allocator1,
	   class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& A,
	   Matrix<T2, Prop2, RowMajorCollection, Allocator2>& B)
  {
    int na = A.GetNmatrix();
    int ma = A.GetMmatrix();

    typedef typename T2::value_type value_type;

    for (int i = 0; i < ma; i++ )
      for (int j = 0; j < na; j++)
	Add(value_type(alpha), A.GetMatrix(i, j), B.GetMatrix(i, j));
  }


  //! Adds two matrices.
  /*! It performs the operation \f$ B = \alpha A + B \f$ where \f$ \alpha \f$
    is a scalar, and \f$ A \f$ and \f$ B \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in,out] B matrix, result of the addition of \a B (on entry) and \a
    A times \a alpha.
  */
  template<class T0,
	   class T1, class Prop1, class Allocator1,
	   class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& A,
	   Matrix<T2, Prop2, ColMajorCollection, Allocator2>& B)
  {
    int na = A.GetNmatrix();
    int ma = A.GetMmatrix();

    typedef typename T2::value_type value_type;

    for (int i = 0; i < ma; i++ )
      for (int j = 0; j < na; j++)
	Add(value_type(alpha), A.GetMatrix(i, j), B.GetMatrix(i, j));
  }


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class Allocator2>
  void Add_heterogeneous(const T0 alpha,
			 const  Matrix<T1, Prop1, Storage1, Allocator1 >& ma,
			 Matrix<FloatDouble, General,
                         DenseSparseCollection, Allocator2>& B,
			 int i, int j)
  {
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::float_dense_m m0b;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::float_sparse_m m1b;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::double_dense_m m2b;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::double_sparse_m m3b;

    T0 alpha_ = alpha;

    switch (B.GetType(i, j))
      {
      case 0:
	B.GetMatrix(i, j, m0b);
	Add(float(alpha_), ma, m0b);
	B.SetMatrix(i, j, m0b);
	m0b.Nullify();
	break;
      case 1:
	B.GetMatrix(i, j, m1b);
	Add(float(alpha_), ma, m1b);
	B.SetMatrix(i, j, m1b);
	m1b.Nullify();
	break;
      case 2:
	B.GetMatrix(i, j, m2b);
	Add(double(alpha_), ma, m2b);
	B.SetMatrix(i, j, m2b);
	m2b.Nullify();
	break;
      case 3:
	B.GetMatrix(i, j, m3b);
	Add(double(alpha_), ma, m3b);
	B.SetMatrix(i, j, m3b);
	m3b.Nullify();
	break;
      default:
	throw WrongArgument("Add_heterogeneous(alpha, Matrix<FloatDouble, "
			    "DenseSparseCollection>, Matrix<FloatDouble,"
			    "DenseSparseCollection> )",
			    "Underlying matrix (" + to_str(i) + " ,"
			    + to_str(j) + " ) not defined.");
      }
  }


  //! Multiplies a FloatDouble collection by a scalar.
  /*!
    \param[in] alpha scalar.
    \param[in,out] M matrix to be multiplied.
  */
  template <class T0, class Allocator1, class Allocator2>
  void Add(const T0 alpha,
	   const Matrix<FloatDouble, General,
           DenseSparseCollection, Allocator1>& A,
	   Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>& B)
  {
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::float_dense_m m0a;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::float_sparse_m m1a;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::double_dense_m m2a;
    typename Matrix<FloatDouble, General, DenseSparseCollection, Allocator2>
      ::double_sparse_m m3a;

    T0 alpha_ = alpha;

    for (int i = 0; i < B.GetMmatrix(); i++)
      for (int j = 0; j < B.GetNmatrix(); j++)
	{
	  switch (B.GetType(i, j))
	    {
	    case 0:
	      A.GetMatrix(i, j, m0a);
	      Add_heterogeneous(float(alpha_), m0a, B, i, j);
	      m0a.Nullify();
	      break;
	    case 1:
	      A.GetMatrix(i, j, m1a);
	      Add_heterogeneous(float(alpha_), m1a, B, i, j);
	      m1a.Nullify();
	      break;
	    case 2:
	      A.GetMatrix(i, j, m2a);
	      Add_heterogeneous(double(alpha_), m2a, B, i, j);
	      m2a.Nullify();
	      break;
	    case 3:
	      A.GetMatrix(i, j, m3a);
	      Add_heterogeneous(double(alpha_), m3a, B, i, j);
	      m3a.Nullify();
	      break;
	    default:
	      throw
                WrongArgument("Add(alpha, Matrix<FloatDouble, "
                              "DenseSparseCollection>, Matrix<FloatDouble,"
                              "DenseSparseCollection> )",
                              "Underlying matrix (" + to_str(i) + " ,"
                              + to_str(j) + " ) not defined.");
	    }
	}
  }


  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, RowMajor, Allocator1>& A,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& B)
  {
    throw Undefined("void Add(const T0& alpha,"
		    "const Matrix<T1, Prop1, RowMajor, Allocator1>& A,"
		    "Matrix<T2, Prop2, RowSparse, Allocator2>& B)");
  }


  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ColMajor, Allocator1>& A,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& B)
  {
    throw Undefined("void Add(const T0& alpha,"
		    "const Matrix<T1, Prop1, RowMajor, Allocator1>& A,"
		    "Matrix<T2, Prop2, RowSparse, Allocator2>& B)");
  }


  template<class T0, class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, Storage1, Allocator1>& A,
	   Matrix<T2, Symmetric, Storage2, Allocator2>& B)
  {
    int i, j;
    for (i = 0; i < A.GetM(); i++)
      for (j = i; j < A.GetN(); j++)
	B(i, j) += alpha * A(i, j);
  }


  //! Adds two matrices.
  /*! It performs the operation \f$ B = \alpha A + B \f$ where \f$ \alpha \f$
    is a scalar, and \f$ A \f$ and \f$ B \f$ are matrices.
    \param[in] alpha scalar.
    \param[in] A matrix.
    \param[in,out] B matrix, result of the addition of \a B (on entry) and \a
    A times \a alpha.
  */
  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& B)
  {
#ifdef SELDON_CHECK_BOUNDS
    if (A.GetM() != B.GetM() || A.GetN() != B.GetN())
      throw WrongDim("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
                     "Unable to add a " + to_str(A.GetM()) + " x "
                     + to_str(A.GetN()) + " matrix with a "
                     + to_str(B.GetM()) + " x " + to_str(B.GetN())
                     + " matrix.");
#endif

    int i = 0;
    int j = 0;
    int k;

    if (A.GetNonZeros() == B.GetNonZeros())
      // A and B might have the same structure.
      {
        // Loop over all non-zeros. If the structures of A and B differ at any
        // time, the loop is broken and a different strategy is undertaken.
        for (i = 0; i < A.GetM(); i++)
          if (A.GetPtr()[i + 1] == B.GetPtr()[i + 1])
            {
              for (j = A.GetPtr()[i]; j < A.GetPtr()[i + 1]; j++)
                if (A.GetInd()[j] == B.GetInd()[j])
                  B.GetData()[j] += alpha * A.GetData()[j];
                else
                  break;
              if (j != A.GetPtr()[i + 1])
                break;
            }
          else
            break;
        // Success: A and B have the same structure.
        if (i == A.GetM())
          return;
      }

    // The addition is performed row by row in the following lines. Thus the
    // additions already performed in the current line, if started, should be
    // canceled.
    for (k = A.GetPtr()[i]; k < j; k++)
      if (A.GetInd()[k] == B.GetInd()[k])
        B.GetData()[k] -= alpha * A.GetData()[k];

    // Number of non zero entries currently found.
    int Nnonzero = A.GetPtr()[i];

    // A and B do not have the same structure. An intermediate matrix will be
    // needed. The first i rows have already been added. These computations
    // should be preserved.
    int* c_ptr = NULL;
    int* c_ind = NULL;
    T2* c_data = NULL;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	c_ptr = reinterpret_cast<int*>(calloc(A.GetM() + 1, sizeof(int)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        c_ptr = NULL;
      }

    if (c_ptr == NULL)
      throw NoMemory("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
		     "Unable to allocate memory for an array of "
		     + to_str(A.GetM() + 1) + " integers.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

        // Reallocates 'c_ind' and 'c_data' for the first i rows.
        c_ind = reinterpret_cast<int*>
          (realloc(reinterpret_cast<void*>(c_ind), Nnonzero * sizeof(int)));
        c_data = reinterpret_cast<T2*>
          (B.GetAllocator().reallocate(c_data, Nnonzero));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        c_ind = NULL;
        c_data = NULL;
      }

    if ((c_ind == NULL || c_data == NULL)
        && Nnonzero != 0)
      throw NoMemory("Add(alpha, const Matrix<RowSparse>& A, "
                     "Matrix<RowSparse>& B)",
                     "Unable to allocate memory for an array of "
                     + to_str(Nnonzero) + " integers and for an array of "
                     + to_str(sizeof(T2) * Nnonzero) + " bytes.");
#endif

    // The pointers of the first i rows are correct.
    for (k = 0; k < i + 1; k++)
      c_ptr[k] = B.GetPtr()[k];

    // Copies the elements from the first i rows, as they were already added.
    for (k = 0; k < Nnonzero; k++)
      {
        c_ind[k] = B.GetInd()[k];
        c_data[k] = B.GetData()[k];
      }

    int Nnonzero_row_max;
    int Nnonzero_max;
    int ja, jb(0), ka, kb;
    // Now deals with the remaining lines.
    for (; i < A.GetM(); i++)
      {
        Nnonzero_row_max = A.GetPtr()[i + 1] - A.GetPtr()[i]
          + B.GetPtr()[i + 1] - B.GetPtr()[i];
        // Maximum number of non zero entries up to row i.
        Nnonzero_max = Nnonzero + Nnonzero_row_max;

#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif

            c_ind = reinterpret_cast<int*>
              (realloc(reinterpret_cast<void*>(c_ind),
                       Nnonzero_max * sizeof(int)));
            c_data = reinterpret_cast<T2*>
              (B.GetAllocator().reallocate(c_data, Nnonzero_max));

#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            c_ind = NULL;
            c_data = NULL;
          }

        if ((c_ind == NULL || c_data == NULL)
            && Nnonzero_max != 0)
          throw NoMemory("Add(alpha, const Matrix<RowSparse>& A, "
                         "Matrix<RowSparse>& B)",
                         "Unable to allocate memory for an array of "
                         + to_str(Nnonzero_max) + " integers and for an "
                         "array of "
                         + to_str(sizeof(T2) * Nnonzero_max) + " bytes.");
#endif

        kb = B.GetPtr()[i];
        if (kb < B.GetPtr()[i + 1])
          jb = B.GetInd()[kb];
        for (ka = A.GetPtr()[i]; ka < A.GetPtr()[i + 1]; ka++)
          {
            ja = A.GetInd()[ka];
            while (kb < B.GetPtr()[i + 1] && jb < ja)
              // For all elements in B that are before the ka-th element of A.
              {
                c_ind[Nnonzero] = jb;
                c_data[Nnonzero] = B.GetData()[kb];
                kb++;
                if (kb < B.GetPtr()[i + 1])
                  jb = B.GetInd()[kb];
                Nnonzero++;
              }

            if (kb < B.GetPtr()[i + 1] && ja == jb)
              // The element in A is also in B.
              {
                c_ind[Nnonzero] = jb;
                c_data[Nnonzero] = B.GetData()[kb] + alpha * A.GetData()[ka];
                kb++;
                if (kb < B.GetPtr()[i + 1])
                  jb = B.GetInd()[kb];
              }
            else
              {
                c_ind[Nnonzero] = ja;
                c_data[Nnonzero] = alpha * A.GetData()[ka];
              }
            Nnonzero++;
          }

        // The remaining elements from B.
        while (kb < B.GetPtr()[i + 1])
          {
            c_ind[Nnonzero] = jb;
            c_data[Nnonzero] = B.GetData()[kb];
            kb++;
            if (kb < B.GetPtr()[i + 1])
              jb = B.GetInd()[kb];
            Nnonzero++;
          }

        c_ptr[i + 1] = Nnonzero;
      }

    B.SetData(B.GetM(), B.GetN(), Nnonzero, c_data, c_ptr, c_ind);
  }


  // ADD //
  /////////


  ///////////
  // GETLU //


  //! Returns the LU factorization of a matrix.
  /*! It factorizes the matrix \a A into \a L and \a U, so that \f$ A = L U
    \f$, \a L is a lower triangular matrix with ones on the diagonal, and \a U
    is an upper triangular matrix. On exit, the LU factorization is stored
    inside \a A: \a L in the lower part and \a U in the upper part. The
    diagonal elements are those of \a U. The diagonal elements of \a L are
    known to be ones.
    \param[in,out] A on entry, the matrix to be factorized; on exit, the LU
    factorization.
    \sa Seldon::SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& M, Vector<T1, Storage1, Allocator1>& Y)
  */
  template <class T0, class Prop0, class Storage0, class Allocator0>
  void GetLU(Matrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    int i, p, q, k;
    T0 temp;

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("GetLU(A)", "The matrix must be squared.");
#endif

    for (i = 0; i < ma; i++)
      {
	for (p = i; p < ma; p++)
	  {
	    temp = 0;
	    for (k = 0; k < i; k++)
	      temp += A(p, k) * A(k, i);
	    A(p, i) -= temp;
	  }
	for (q = i+1; q < ma; q++)
	  {
	    temp = 0;
	    for (k = 0; k < i; k++)
	      temp += A(i, k) * A(k, q);
	    A(i, q) = (A(i,q ) - temp) / A(i, i);
	  }
      }
  }


  // GETLU //
  ///////////


  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that A B + C -> C is possible according to the dimensions of
    the matrices A, B and C. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param A matrix.
    \param B matrix.
    \param C matrix.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function)
  {
    if (B.GetM() != A.GetN() || C.GetM() != A.GetM() || B.GetN() != C.GetN())
      throw WrongDim(function, string("Operation A B + C -> C not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that A B + C -> C or B A + C -> C is possible according to the
    dimensions of the matrices A, B and C. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param side side by which A is multiplied by B.
    \param A matrix.
    \param B matrix.
    \param C matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const SeldonSide& side,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function)
  {
    if ( SeldonSide(side).Left() &&
	 (B.GetM() != A.GetN() || C.GetM() != A.GetM()
	  || B.GetN() != C.GetN()) )
      throw WrongDim(function, string("Operation A B + C -> C not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
    else if ( SeldonSide(side).Right() &&
	      (B.GetN() != A.GetM() || C.GetM() != B.GetM()
	       || A.GetN() != C.GetN()) )
      throw WrongDim(function, string("Operation B A + C -> C not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that A B -> C is possible according to the dimensions of the
    matrices A and B. If the dimensions are incompatible, an exception is
    raised (a WrongDim object is thrown).
    \param TransA status of A, e.g. transposed.
    \param A matrix.
    \param TransB status of B, e.g. transposed.
    \param B matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const SeldonTranspose& TransA,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const SeldonTranspose& TransB,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function)
  {
    SeldonTranspose status_A(TransA);
    SeldonTranspose status_B(TransB);
    string op;
    if (status_A.Trans())
      op = string("A'");
    else if (status_A.ConjTrans())
      op = string("A*");
    else
      op = string("A");
    if (status_B.Trans())
      op += string(" B'");
    else if (status_B.ConjTrans())
      op += string(" B*");
    else
      op += string(" B");
    op = string("Operation ") + op + string(" not permitted:");
    if (B.GetM(status_B) != A.GetN(status_A))
      throw WrongDim(function, op
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that A B + C -> C is possible according to the dimensions of the
    matrices A, B and C. If the dimensions are incompatible, an exception is
    raised (a WrongDim object is thrown).
    \param TransA status of A, e.g. transposed.
    \param A matrix.
    \param TransB status of B, e.g. transposed.
    \param B matrix.
    \param C matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const SeldonTranspose& TransA,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const SeldonTranspose& TransB,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function)
  {
    string op;
    if (TransA.Trans())
      op = string("A'");
    else if (TransA.ConjTrans())
      op = string("A*");
    else
      op = string("A");
    if (TransB.Trans())
      op += string(" B' + C");
    else if (TransB.ConjTrans())
      op += string(" B* + C");
    else
      op += string(" B + C");
    op = string("Operation ") + op + string(" not permitted:");
    if (B.GetM(TransB) != A.GetN(TransA) || C.GetM() != A.GetM(TransA)
	|| B.GetN(TransB) != C.GetN())
      throw WrongDim(function, op
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that A B is possible according to the dimensions of the matrices
    A and B. If the dimensions are incompatible, an exception is raised (a
    WrongDim object is thrown).
    \param A matrix.
    \param B matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function)
  {
    CheckDim(SeldonLeft, A, B, function);
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that A B or B A is possible according to the dimensions of
    the matrices A and B. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param side side by which A is multiplied by B.
    \param A matrix.
    \param B matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const SeldonSide& side,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function)
  {
    if (side.Left() && B.GetM() != A.GetN())
      throw WrongDim(function, string("Operation A B not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix."));
    else if (side.Right() && B.GetN() != A.GetM())
      throw WrongDim(function, string("Operation B A not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix."));
  }


  // CHECKDIM //
  //////////////


  ///////////
  // NORMS //


  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in absolute value) of all elements of \a A.
  */
  template <class T, class Prop, class Storage, class Allocator>
  T MaxAbs(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    T res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
	res = max(res, abs(A(i, j)) );

    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T Norm1(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int j = 0; j < A.GetN(); j++)
      {
	sum = T(0);
	for (int i = 0; i < A.GetM(); i++)
	  sum += abs( A(i, j) );

	res = max(res, sum);
      }

    return res;
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T NormInf(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
	sum = T(0);
	for (int j = 0; j < A.GetN(); j++)
	  sum += abs( A(i, j) );

	res = max(res, sum);
      }

    return res;
  }


  //! Returns the maximum (in modulus) of a matrix.
  /*!
    \param[in] A matrix.
    \return The maximum (in modulus) of all elements of \a A.
  */
  template <class T, class Prop, class Storage, class Allocator>
  T MaxAbs(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    T res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
	{
	  res = max(res, abs(A(i, j)) );
	}

    return res;
  }


  //! Returns the 1-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T Norm1(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int j = 0; j < A.GetN(); j++)
      {
	sum = T(0);
	for (int i = 0; i < A.GetM(); i++)
	  sum += abs( A(i, j) );

	res = max(res, sum);
      }

    return res;
  }


  //! Returns the infinity-norm of a matrix.
  /*!
    \param[in] A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T NormInf(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
	sum = T(0);
	for (int j = 0; j < A.GetN(); j++)
	  sum += abs( A(i, j) );

	res = max(res, sum);
      }

    return res;
  }


  // NORMS //
  ///////////


  ///////////////
  // TRANSPOSE //


  //! Matrix transposition.
  template<class T, class Prop, class Storage, class Allocator>
  void Transpose(Matrix<T, Prop, Storage, Allocator>& A)
  {
    int m = A.GetM();
    int n = A.GetN();

    if (m == n)
      {
	T tmp;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < i; j++)
	    {
	      tmp = A(i,j);
	      A.Get(i,j) = A(j,i);
	      A.Get(j,i) = tmp;
	    }
      }
    else
      {
	Matrix<T, Prop, Storage, Allocator> B;
	B = A;
	A.Reallocate(n,m);
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < n; j++)
	    A.Get(j,i) = B(i,j);
      }
  }


  //! Matrix transposition.
  template<class T, class Allocator>
  void Transpose(Matrix<T, General, RowSparse, Allocator>& A)
  {
    int m = A.GetM();
    int n = A.GetN();
    int nnz = A.GetDataSize();
    Vector<int, VectFull, CallocAlloc<int> > ptr_T(n+1), ptr;
    Vector<int, VectFull, CallocAlloc<int> > ind_T(nnz), ind;
    Vector<T, VectFull, Allocator> data_T(nnz), data;

    ptr.SetData(m+1, A.GetPtr());
    ind.SetData(nnz,  A.GetInd());
    data.SetData(nnz, A.GetData());

    ptr_T.Zero();
    ind_T.Zero();
    data_T.Zero();

    // For each column j, computes number of its non-zeroes and stores it in
    // ptr_T[j].
    for (int i = 0; i < nnz; i++)
      ptr_T(ind(i) + 1)++;

    // Computes the required number of non-zeroes ptr_T(j).
    for (int j = 1; j < n + 1; j++)
      ptr_T(j) += ptr_T(j - 1);

    Vector<int, VectFull, CallocAlloc<int> > row_ind(n+1);
    row_ind.Fill(0);
    for (int i = 0; i < m; i++)
      for (int jp = ptr(i); jp < ptr(i+1); jp++)
    	{
	  int j = ind(jp);
	  int k = ptr_T(j) + row_ind(j);
	  ++row_ind(j);
	  data_T(k) = data(jp);
	  ind_T(k) = i;
    	}

    A.SetData(n, m, data_T, ptr_T, ind_T);

    data.Nullify();
    ptr.Nullify();
    ind.Nullify();
  }


  //! Matrix transposition and conjugation.
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(Matrix<T, Prop, Storage, Allocator>& A)
  {
    int i, j;

    int m = A.GetM();
    int n = A.GetN();

    if (m == n)
      {
	T tmp;
	for (i = 0; i < m; i++)
	  {
            // Extra-diagonal part.
            for (j = 0; j < i; j++)
              {
                tmp = A(i, j);
                A(i, j) = conj(A(j, i));
                A(j, i) = conj(tmp);
              }

            // Diagonal part.
            A(i, i) = conj(A(i, i));
          }
      }
    else
      {
	Matrix<T, Prop, Storage, Allocator> B;
	B = A;
	A.Reallocate(n, m);
	for (i = 0; i < m; i++)
	  for (j = 0; j < n; j++)
	    A(j, i) = conj(B(i, j));
      }
  }


  // TRANSPOSE //
  ///////////////


  ///////////////////////
  // ISSYMMETRICMATRIX //


  //! returns true if the matrix is symmetric
  template<class T, class Prop, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    return false;
  }


  //! returns true if the matrix is symmetric
  template<class T, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Symmetric, Storage, Allocator>& A)
  {
    return true;
  }


  // ISSYMMETRICMATRIX //
  ///////////////////////


} // namespace Seldon.


#endif
