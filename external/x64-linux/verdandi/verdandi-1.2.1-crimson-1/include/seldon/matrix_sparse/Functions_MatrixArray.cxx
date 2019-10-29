// Copyright (C) 2003-2009 Marc Duruflé
// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_CXX

/*
  Functions defined in this file:
  (storage ArrayRowSparse, ArrayColSparse, etc)

  alpha.M*X + beta.Y -> Y
  MltAdd(alpha, M, X, beta, Y)

  alpha.A + B -> B
  Add(alpha, A, B)

  alpha.M -> M
  Mlt(alpha, M)

*/

namespace Seldon
{


  ////////////
  // MltAdd //


  /*** ArrayRowSymComplexSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) += complex<T1>(0, val) * B(i);
		else
		  {
		    C(i) += complex<T1>(0, val) * B(p);
		    C(p) += complex<T1>(0, val) * B(i);
		  }
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonConjTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(),n,p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		if (p == i)
		  C(i) -= complex<T1>(0, val) * B(i);
		else
		  {
		    C(i) -= complex<T1>(0, val) * B(p);
		    C(p) -= complex<T1>(0, val) * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		if (p == i)
		  C(i) += val_cplx * B(i);
		else
		  {
		    C(i) += val_cplx * B(p);
		    C(p) += val_cplx * B(i);
		  }
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		if (p == i)
		  C(i) -= val_cplx * B(i);
		else
		  {
		    C(i) -= val_cplx * B(p);
		    C(p) -= val_cplx * B(i);
		  }
	      }
	  }
      }
  }


  /*** ArrayRowComplexSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(), n, p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1, 0))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(i) += val * B(p);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(i) += complex<T1>(0, val) * B(p);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		C(i) += val_cplx * B(p);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		C(i) += val_cplx * B(p);
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(),n,p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) += complex<T1>(0, val) * B(i);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		C(p) += val_cplx * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx = alpha * complex<T1>(0, A.ValueImag(i, k));
		C(p) += val_cplx * B(i);
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonConjTrans& Trans, const Matrix<T1, General,
	      ArrayRowComplexSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(),n,p;
    T1 val;
    T3 val_cplx;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val = A.ValueReal(i, k);
		C(p) += val * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val = A.ValueImag(i, k);
		C(p) -= complex<T1>(0, val) * B(i);
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRealRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexReal(i, k);
		val_cplx = alpha * A.ValueReal(i, k);
		C(p) += val_cplx * B(i);
	      }
	    n = A.GetImagRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.IndexImag(i, k);
		val_cplx -= alpha * complex<T1>(0, A.ValueImag(i, k));
		C(p) += val_cplx * B(i);
	      }
	  }
      }
  }


  /*** ArrayRowSymSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha, const Matrix<T1, Symmetric,
	      ArrayRowSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (B.GetM() <= 0)
      return;

    T3 zero(B(0));
    zero *= 0;

    if (beta == zero)
      C.Fill(zero);
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T3 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);

		if (p == i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
    else
      {
	// alpha different from 1
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = alpha * A.Value(i, k);

		if (p==i)
		  C(i) += val * B(i);
		else
		  {
		    C(i) += val * B(p);
		    C(p) += val * B(i);
		  }
	      }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans, const Matrix<T1, Symmetric,
	      ArrayRowSymSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  /*** ArrayRowSparse ***/


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);

    int m = A.GetM(), n, p;
    T1 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(i) += val * B(p);
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(i) += alpha * val * B(p);
              }
	  }
      }
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T4, VectFull, Allocator3>& C)
  {
    MltAdd(alpha, A, B, beta, C);
  }


  template<class T0, class T1, class T2, class T3, class T4,
	   class Allocator1, class Allocator2, class Allocator3>
  void MltAdd(const T0& alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	      const Vector<T2, VectFull, Allocator2>& B,
	      const T4& beta,
	      Vector<T3, VectFull, Allocator3>& C)
  {
    if (beta == T4(0))
      C.Fill(T3(0));
    else
      Mlt(beta, C);
    int m = A.GetM(), n, p;
    T1 val;
    if (alpha == T0(1))
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(p) += val * B(i);
	      }
	  }
      }
    else // alpha != 1.
      {
	for (int i = 0 ; i < m ; i++)
	  {
	    n = A.GetRowSize(i);
	    for (int k = 0; k < n ; k++)
	      {
		p = A.Index(i, k);
		val = A.Value(i, k);
		C(p) += alpha * val * B(i);
	      }
	  }
      }
  }


  // MltAdd //
  ////////////



  /////////
  // Add //


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = T2(A.Value(i, j));

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.Value(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetIndex(i), value.GetData());
      }
  }


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric,
	   ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T2> > value;
    IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
	value.Reallocate(n + ni);
	index.Reallocate(n + ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.ValueReal(i, j);
	    index(j) = A.IndexReal(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(j+n) = complex<T2>(0, 1) * A.ValueImag(i, j);
	    index(j+n) = A.IndexImag(i, j);
	  }

	Mlt(alpha, value);
	B.AddInteractionRow(i, n+ni, index, value);
      }
  }


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, Symmetric,
	   ArrayRowSymComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B)
  {
    int m = B.GetM(), n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.ValueReal(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
	n = A.GetImagRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = complex<T1>(0, 1) * A.ValueImag(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  template<class T0, class T1, class T2, class Allocator1, class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, Symmetric, ArrayRowSparse, Allocator2>& B)
  {
    int m = B.GetM(),n;
    Vector<T2, VectFull, Allocator2> value;
    for (int i = 0; i < m; i++)
      {
	n = A.GetRealRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = A.ValueReal(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetRealInd(i), value.GetData());
	n = A.GetImagRowSize(i);
	value.Reallocate(n);
	for (int j = 0; j < n; j++)
	  value(j) = complex<T1>(0, 1) * A.ValueImag(i, j);

	Mlt(alpha, value);
	B.AddInteractionRow(i, n, A.GetImagInd(i), value.GetData());
      }
  }


  template<class T0, class T1, class T2, class Allocator1,class Allocator2>
  void Add(const T0& alpha, const Matrix<T1, General,
	   ArrayRowComplexSparse, Allocator1>& A,
	   Matrix<T2, General, ArrayRowComplexSparse, Allocator2>& B)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T2> > value; IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRealRowSize(i);
	ni = A.GetImagRowSize(i);
	value.Reallocate(n + ni);
	index.Reallocate(n + ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.ValueReal(i, j);
	    index(j) = A.IndexReal(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = complex<T2>(0, 1) * A.ValueImag(i, j);
	    index(n+j) = A.IndexImag(i, j);
	  }

	Mlt(alpha, value);
	B.AddInteractionRow(i, n+ni, index, value);
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3, class Allocator1,
	   class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, General, ArrayRowSparse, Allocator2>& B,
	   Matrix<complex<T3>, General, ArrayRowSparse, Allocator3>& C)
  {
    int m = B.GetM(),n1,n2,size_row;;
    Vector<complex<T3>, VectFull, Allocator3> val_row;
    IVect ind_row;
    for (int i = 0 ; i < m ; i++)
      {
	n1 = A.GetRowSize(i);
	n2 = B.GetRowSize(i);
	size_row = n1 + n2;
	val_row.Reallocate(size_row);
	ind_row.Reallocate(size_row);
	for (int j = 0 ; j < n1 ; j++)
	  {
	    ind_row(j) = A.Index(i, j);
	    val_row(j) = alpha*complex<T3>(A.Value(i, j), 0);
	  }

	for (int j = 0 ; j < n2 ; j++)
	  {
	    ind_row(j+n1) = B.Index(i, j);
	    val_row(j+n1) = alpha * complex<T3>(B.Value(i, j));
	  }

	C.AddInteractionRow(i, size_row, ind_row, val_row);
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
	   const Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B,
	   Matrix<complex<T3>, Symmetric, ArrayRowSymSparse, Allocator3>& C)
  {
    int m = B.GetM(), n1, n2, size_row;
    Vector<complex<T3>, VectFull, Allocator3> val_row;
    IVect ind_row;
    for (int i = 0 ; i < m ; i++)
      {
	n1 = A.GetRowSize(i);
	n2 = B.GetRowSize(i);
	size_row = n1 + n2;
	val_row.Reallocate(size_row);
	ind_row.Reallocate(size_row);
	for (int j = 0 ; j < n1 ; j++)
	  {
	    ind_row(j) = A.Index(i, j);
	    val_row(j) = alpha * complex<T3>(A.Value(i, j), 0);
	  }

	for (int j = 0 ; j < n2 ; j++)
	  {
	    ind_row(j+n1) = B.Index(i, j);
	    val_row(j+n1) = alpha * complex<T3>(B.Value(i, j));
	  }

	C.AddInteractionRow(i, size_row, ind_row, val_row);
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, ArrayRowSymSparse, Allocator1>& A,
	   const Matrix<T2, Symmetric, ArrayRowSymSparse, Allocator2>& B,
	   Matrix<T3, Symmetric, ArrayRowSymComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T3> > value; IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	value.Reallocate(n+ni);
	index.Reallocate(n+ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.Value(i, j);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = complex<T3>(0, 1) * B.Value(i, j);
	    index(n+j) = B.Index(i, j);
	  }

	Mlt(alpha, value);
	C.AddInteractionRow(i, n+ni, index, value);
      }
  }


  // C = C + complex(A,B)
  template<class T0, class T1, class T2, class T3,
	   class Allocator1, class Allocator2, class Allocator3>
  void Add(const T0& alpha,
	   const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, General, ArrayRowSparse, Allocator2>& B,
	   Matrix<T3, General, ArrayRowComplexSparse, Allocator3>& C)
  {
    int m = B.GetM(), n, ni;
    Vector<complex<T3> > value;
    IVect index;
    for (int i = 0 ; i < m ; i++)
      {
	n = A.GetRowSize(i);
	ni = B.GetRowSize(i);
	value.Reallocate(n + ni);
	index.Reallocate(n + ni);
	for (int j = 0; j < n; j++)
	  {
	    value(j) = A.Value(i, j);
	    index(j) = A.Index(i, j);
	  }

	for (int j = 0; j < ni; j++)
	  {
	    value(n+j) = complex<T3>(0, 1) * B.Value(i, j);
	    index(n+j) = B.Index(i, j);
	  }

	Mlt(alpha, value);
	C.AddInteractionRow(i, n+ni, index, value);
      }
  }


  // Add //
  /////////



  /////////
  // Mlt //


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, General, ArrayRowSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, General, ArrayColSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Symmetric, ArrayRowSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, Symmetric, ArrayColSymSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetN(); i++)
      for (int j = 0; j < A.GetColSize(i); j++)
	A.Value(i,j) *= alpha;
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha,
	   Matrix<T, General, ArrayRowComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= alpha;
	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= alpha;
      }
  }


  template<class T0, class T, class Allocator>
  void Mlt(const T0& alpha, Matrix<T, Symmetric,
	   ArrayRowSymComplexSparse, Allocator>& A)
  {
    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetRealRowSize(i); j++)
	  A.ValueReal(i,j) *= alpha;

	for (int j = 0; j < A.GetImagRowSize(i); j++)
	  A.ValueImag(i,j) *= alpha;
      }
  }


  // Matrix-matrix product (sparse matrix against full matrix)
  template<class T1, class Allocator1, class T2, class Prop2,
	   class Storage2, class Allocator2, class T3, class Prop3,
	   class Storage3, class Allocator3>
  void Mlt(const Matrix<T1, General, ArrayRowSparse, Allocator1>& A,
	   const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	   Matrix<T3, Prop3, Storage3, Allocator3>& C)
  {
    int m = A.GetM();
    int n = B.GetN();
    C.Reallocate(m,n);
    T3 val;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
	{
	  val = T3(0);
	  for (int ind = 0; ind < A.GetRowSize(i); ind++)
	    {
	      int k = A.Index(i, ind);
	      val += A.Value(i, ind) * B(k, j);
	    }
	  C(i, j) = val;
	}
  }


  // Mlt //
  /////////


} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_MATRIX_ARRAY_CXX
#endif
