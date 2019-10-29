// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_HXX
#define SELDON_FILE_FUNCTIONS_MATRIX_HXX


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


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void Mlt(const T0 alpha, Matrix<T1, Prop1, Storage1, Allocator1>& A);

  template <class T0,
	    class T1, class Prop1, class Allocator1>
  void Mlt(const T0 alpha,
           Matrix<T1, Prop1, ColMajorCollection, Allocator1>& A);

  template <class T0,
	    class T1, class Prop1, class Allocator1>
  void Mlt(const T0 alpha,
           Matrix<T1, Prop1, RowMajorCollection, Allocator1>& A);

  template <class T0, class Allocator>
  void Mlt(const T0 alpha,
           Matrix<FloatDouble, General, DenseSparseCollection, Allocator>& A);

  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3, class Prop3, class Storage3, class Allocator3>
  void Mlt(const T0 alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	   const Matrix<T2, Prop2, Storage2, Allocator2>& B,
           Matrix<T3, Prop3, Storage3, Allocator3>& C);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& B,
           Matrix<T2, Prop2, Storage2, Allocator2>& C);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
	   const Matrix<T1, Prop1, RowSparse, Allocator1>& B,
           Matrix<T2, Prop2, RowSparse, Allocator2>& C);

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2>
  void MltNoTransTrans(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
                       const Matrix<T1, Prop1, RowSparse, Allocator1>& B,
                       Matrix<T2, Prop2, RowSparse, Allocator2>& C);

  // MLT //
  /////////


  ////////////
  // MLTADD //


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	      const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, Storage4, Allocator4>& C);

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
	      Matrix<T4, Prop4, Storage4, Allocator4>& C);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& A,
              const Matrix<T2, General, RowMajor, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, PETScMPIDense, Allocator4>& C);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& A,
	      const Matrix<T2, Prop2, RowMajorCollection, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, RowMajorCollection, Allocator4>& C);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& A,
	      const Matrix<T2, Prop2, ColMajorCollection, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, ColMajorCollection, Allocator4>& C);

  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, General, RowMajor, Allocator1>& A,
	      const Matrix<T2, General, RowMajor, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, General, RowSparse, Allocator4>& C);

  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, General, RowMajor, Allocator1>& A,
	      const Matrix<T2, General, RowSparse, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, General, RowSparse, Allocator4>& C);

  template <class T0,
	    class T1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, General, RowSparse, Allocator1>& A,
	      const Matrix<T2, General, RowMajor, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, General, RowSparse, Allocator4>& C);

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
			    int i, int j);

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
			     int j, int k);

  template <class T0, class Allocator1, class Allocator2, class T3,
	    class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<FloatDouble, General, DenseSparseCollection,
	      Allocator1>& A,
	      const Matrix<FloatDouble, General, DenseSparseCollection,
	      Allocator2>& B,
	      const T3 beta,
	      Matrix<FloatDouble, General, DenseSparseCollection,
	      Allocator4>& C);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
            class T4, class Prop4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
              const Matrix<T2, Prop2, RowSparse, Allocator2>& B,
              const T3 beta,
              Matrix<T4, Prop4, RowSparse, Allocator4>& C);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Prop2, class Allocator2,
	    class T3,
            class T4, class Prop4, class Allocator4>
  void MltNoTransTransAdd(const T0 alpha,
                          const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
                          const Matrix<T2, Prop2, RowSparse, Allocator2>& B,
                          const T3 beta,
                          Matrix<T4, Prop4, RowSparse, Allocator4>& C);

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
              Matrix<T4, Prop4, RowSparse, Allocator4>& C);


  // MLTADD //
  ////////////


  /////////
  // ADD //


  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
	   class T2, class Prop2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	   Matrix<T2, Prop2, Storage2, Allocator2>& B);

  template<class T0,
	   class T1, class Prop1, class Allocator1,
	   class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& A,
	   Matrix<T2, Prop2, RowMajorCollection, Allocator2>& B);

  template<class T0,
	   class T1, class Prop1, class Allocator1,
	   class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& A,
	   Matrix<T2, Prop2, ColMajorCollection, Allocator2>& B);

  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class Allocator2>
  void Add_heterogeneous(const T0 alpha,
			 const  Matrix<T1, Prop1, Storage1, Allocator1 >& ma,
			 Matrix<FloatDouble, General,
                         DenseSparseCollection, Allocator2>& B,
			 int i, int j);

  template <class T0, class Allocator1, class Allocator2>
  void Add(const T0 alpha,
	   const Matrix<FloatDouble, General,
           DenseSparseCollection, Allocator1>& A,
	   Matrix<FloatDouble, General,
           DenseSparseCollection, Allocator2>& B);

  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, RowMajor, Allocator1>& A,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& B);

  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, ColMajor, Allocator1>& A,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& B);

  template<class T0, class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, Storage1, Allocator1>& A,
	   Matrix<T2, Symmetric, Storage2, Allocator2>& B);

  template<class T0, class T1, class Prop1, class Allocator1,
           class T2, class Prop2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Prop1, RowSparse, Allocator1>& A,
	   Matrix<T2, Prop2, RowSparse, Allocator2>& B);


  // ADD //
  /////////


  ///////////
  // GETLU //


  template <class T0, class Prop0, class Storage0, class Allocator0>
  void GetLU(Matrix<T0, Prop0, Storage0, Allocator0>& A);


  // GETLU //
  ///////////


  //////////////
  // CHECKDIM //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const SeldonSide& side,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const SeldonTranspose& TransA,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const SeldonTranspose& TransB,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const SeldonTranspose& TransA,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const SeldonTranspose& TransB,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const SeldonSide& side,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function = "");


  // CHECKDIM //
  //////////////


  ///////////
  // NORMS //


  template <class T, class Prop, class Storage, class Allocator>
  T MaxAbs(const Matrix<T, Prop, Storage, Allocator>& A);

  template <class T, class Prop, class Storage, class Allocator>
  T Norm1(const Matrix<T, Prop, Storage, Allocator>& A);

  template <class T, class Prop, class Storage, class Allocator>
  T NormInf(const Matrix<T, Prop, Storage, Allocator>& A);

  template <class T, class Prop, class Storage, class Allocator>
  T MaxAbs(const Matrix<complex<T>, Prop, Storage, Allocator>& A);

  template <class T, class Prop, class Storage, class Allocator>
  T Norm1(const Matrix<complex<T>, Prop, Storage, Allocator>& A);

  template <class T, class Prop, class Storage, class Allocator>
  T NormInf(const Matrix<complex<T>, Prop, Storage, Allocator>& A);


  // NORMS //
  ///////////


  ///////////////
  // TRANSPOSE //


  template<class T, class Prop, class Storage, class Allocator>
  void Transpose(Matrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Allocator>
  void Transpose(Matrix<T, General, RowSparse, Allocator>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(Matrix<T, Prop, Storage, Allocator>& A);


  // TRANSPOSE //
  ///////////////


  ///////////////////////
  // ISSYMMETRICMATRIX //


  template<class T, class Prop, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Symmetric, Storage, Allocator>& A);


  // ISSYMMETRICMATRIX //
  ///////////////////////


} // namespace Seldon.


#endif
