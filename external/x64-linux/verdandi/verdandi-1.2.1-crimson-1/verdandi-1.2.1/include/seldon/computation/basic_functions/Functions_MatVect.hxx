// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_HXX
#define SELDON_FILE_FUNCTIONS_MATVECT_HXX


/*
  Functions defined in this file:

  M X -> Y
  Mlt(M, X, Y)

  alpha M X -> Y
  Mlt(alpha, M, X, Y)

  M X -> Y
  M^T X -> Y
  Mlt(Trans, M, X, Y)

  alpha M X + beta Y -> Y
  MltAdd(alpha, M, X, beta, Y)

  Gauss(M, X)

  GaussSeidel(M, X, Y, iter)

  SOR(M, X, Y, omega, iter)

  SolveLU(M, Y)

  Solve(M, Y)
*/

namespace Seldon
{


  /////////
  // MLT //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T3& alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);

  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const SeldonTranspose& Trans,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y);


  // MLT //
  /////////


  ////////////
  // MLTADD //


  /*** PETSc matrices ***/

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, PETScSeq, Allocator2>& X,
              const T3 beta, Vector<T4, PETScSeq, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, PETScPar, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScSeq, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIAIJ, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  template <class T0,
            class T1, class Prop1, class Allocator1,
            class T2, class Allocator2,
            class T3,
            class T4, class Allocator4>
  void MltAdd(const T0 alpha,
              const Matrix<T1, Prop1, PETScMPIDense, Allocator1>& M,
              const Vector<T2, VectFull, Allocator2>& X,
              const T3 beta, Vector<T4, PETScPar, Allocator4>& Y);

  /*** Sparse matrices ***/

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Collection, Allocator4>& Y);

  /*** Complex sparse matrices ***/

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  /*** Symmetric sparse matrices ***/

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Collection, Allocator2>& X,
	      const T3 beta, Vector<T4, Collection, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Collection, Allocator4>& Y);

  /*** Symmetric complex sparse matrices ***/

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  /*** Sparse matrices, *Trans ***/

  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  /*** Complex sparse matrices, *Trans ***/

  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  /*** Symmetric sparse matrices, *Trans ***/

  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  /*** Symmetric complex sparse matrices, *Trans ***/

  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);

  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y);


  // MLTADD //
  ////////////


  ////////////
  // MLTADD //


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Collection, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowMajorCollection, Allocator1>& M,
	      const Vector<T2, Collection, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Collection, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Allocator2,
	    class T3,
	    class T4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, ColMajorCollection, Allocator1>& M,
	      const Vector<T2, Collection, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Collection, Allocator4>& Y);

  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y);


  // MLTADD //
  ////////////


  ///////////
  // GAUSS //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  inline void Gauss(Matrix<T0, Prop0, Storage0, Allocator0>& M,
		    Vector<T1, Storage1, Allocator1>& X);

  // GAUSS //
  ///////////


  //////////////////
  // GAUSS-SEIDEL //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  inline void GaussSeidel(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
			  const Vector<T1, Storage1, Allocator1>& X,
			  Vector<T2, Storage2, Allocator2>& Y,
			  int iter);


  // GAUSS-SEIDEL //
  //////////////////



  ///////////////////
  // S.O.R. METHOD //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  inline void SOR(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		  const Vector<T1, Storage1, Allocator1>& X,
		  Vector<T2, Storage2, Allocator2>& Y,
		  T3 omega,
		  int iter);


  // S.O.R. METHOD //
  ///////////////////



  /////////////
  // SOLVELU //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	       Vector<T1, Storage1, Allocator1>& Y);

  // SOLVELU //
  /////////////


  ///////////
  // SOLVE //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetAndSolveLU(Matrix<T0, Prop0, Storage0, Allocator0>& M,
                     Vector<T1, Storage1, Allocator1>& Y);

  // SOLVE //
  ///////////


  //////////////
  // CHECKDIM //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "");

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, RowMajorCollection, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y,
		string function = "");

  template <class T0, class Prop0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, ColMajorCollection, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Collection, Allocator2>& Y,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Collection, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const SeldonTranspose& trans,
		const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "", string op = "M X + Y -> Y");

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		string function = "", string op = "M X");


  // CHECKDIM //
  //////////////


}  // namespace Seldon.


#endif
