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


#ifndef SELDON_FILE_BLAS_2_CXX

/*
  Functions included in this file:

  xTRMV   (Mlt)
  xTPMV   (Mlt)
  xGEMV   (MltAdd)
  xHEMV   (MltAdd)
  xHPMV   (MltAdd)
  xSYMV   (MltAdd)
  xSPMV   (MltAdd)
  xGER    (Rank1Update)
  xGERU   (Rank1Update)
  xGERC   (Rank1Update)
  xSPR    (Rank1Update)
  xHPR    (Rank1Update)
  xSPR2   (Rank2Update)
  xHPR2   (Rank2Update)
  xTRSV   (Solve)
  xTPSV   (Solve)
*/

extern "C"
{
#include "cblas.h"
}

namespace Seldon
{


  /////////
  // MLT //


  /*** ColUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_strmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_strmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtrmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctrmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztrmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_strmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtrmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctrmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztrmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_strmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtrmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctrmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztrmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_strmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_strmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtrmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctrmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowUpTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztrmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_strmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_strmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtrmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctrmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowLoTriang, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztrmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColUpTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_stpmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtpmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctpmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztpmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColUpTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_stpmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtpmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctpmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztpmv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_stpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_stpmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtpmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctpmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, ColLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztpmv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_stpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_stpmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtpmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctpmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowUpTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztpmv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_stpmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_dtpmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ctpmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const Matrix<complex<double>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(M, X)");
#endif

    cblas_ztpmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_stpmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Mlt(const SeldonTranspose& TransA,
	   const SeldonDiag& DiagA,
	   const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	   Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_dtpmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ctpmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Mlt(const SeldonTranspose& TransA,
      const SeldonDiag& DiagA,
      const Matrix<complex<double>, Prop0, RowLoTriangPacked, Allocator0>& A,
      Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Mlt(status, diag, M, X)");
#endif

    cblas_ztpmv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  // MLT //
  /////////



  ////////////
  // MLTADD //


  // General //

  /*** ColMajor and NoTrans ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_sgemv(CblasColMajor, CblasNoTrans,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_dgemv(CblasColMajor, CblasNoTrans,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_cgemv(CblasColMajor, CblasNoTrans,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_zgemv(CblasColMajor, CblasNoTrans,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  /*** ColMajor and TransA ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_sgemv(CblasColMajor, TransA,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_dgemv(CblasColMajor, TransA,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_cgemv(CblasColMajor, TransA,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_zgemv(CblasColMajor, TransA,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  /*** RowMajor and NoTrans ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_sgemv(CblasRowMajor, CblasNoTrans,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_dgemv(CblasRowMajor, CblasNoTrans,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_cgemv(CblasRowMajor, CblasNoTrans,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_zgemv(CblasRowMajor, CblasNoTrans,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  /*** RowMajor and TransA ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_sgemv(CblasRowMajor, TransA,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_dgemv(CblasRowMajor, TransA,
		A.GetM(), A.GetN(), alpha, A.GetData(), A.GetLD(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_cgemv(CblasRowMajor, TransA,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonTranspose& TransA,
	      const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(TransA, A, X, Y, "MltAdd(alpha, status, M, X, beta, Y)");
#endif

    cblas_zgemv(CblasRowMajor, TransA,
		A.GetM(), A.GetN(), reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(A.GetData()), A.GetLD(),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(&beta),
		reinterpret_cast<void*>(Y.GetData()), 1);
  }


  // Hermitian //

  /*** ColHerm and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_chemv(CblasColMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_zhemv(CblasColMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  /*** ColHerm and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_chemv(CblasColMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, ColHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_zhemv(CblasColMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  /*** RowHerm and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_chemv(CblasRowMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_zhemv(CblasRowMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  /*** RowHerm and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_chemv(CblasRowMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0, RowHerm, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_zhemv(CblasRowMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(), A.GetM(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  // HermPacked //

  /*** ColHermPacked and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_chpmv(CblasColMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_zhpmv(CblasColMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  /*** ColHermPacked and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_chpmv(CblasColMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0,
	      ColHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_zhpmv(CblasColMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  /*** RowHermPacked and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const Matrix<complex<float>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_chpmv(CblasRowMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const Matrix<complex<double>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_zhpmv(CblasRowMajor, CblasUpper,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  /*** RowHermPacked and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<float> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<float>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<float>, VectFull, Allocator1>& X,
	      const complex<float> beta,
	      Vector<complex<float>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_chpmv(CblasRowMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const complex<double> alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<complex<double>, Prop0,
	      RowHermPacked, Allocator0>& A,
	      const Vector<complex<double>, VectFull, Allocator1>& X,
	      const complex<double> beta,
	      Vector<complex<double>, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_zhpmv(CblasRowMajor, Uplo,
		A.GetM(), reinterpret_cast<const void*>(&alpha),
		A.GetDataConstVoid(),
		X.GetDataConstVoid(), 1, reinterpret_cast<const void*>(&beta),
		Y.GetDataVoid(), 1);
  }


  // Symmetric //

  /*** ColSym and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_ssymv(CblasColMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_dsymv(CblasColMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  /*** ColSym and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, ColSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_ssymv(CblasColMajor, Uplo,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, ColSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_dsymv(CblasColMajor, Uplo,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  /*** RowSym and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_ssymv(CblasRowMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_dsymv(CblasRowMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  /*** RowSym and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, RowSym, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_ssymv(CblasRowMajor, Uplo,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, RowSym, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_dsymv(CblasRowMajor, Uplo,
		A.GetM(), alpha, A.GetData(), A.GetM(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  // SymPacked //

  /*** ColSymPacked and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_sspmv(CblasColMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_dspmv(CblasColMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  /*** ColSymPacked and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_sspmv(CblasColMajor, Uplo,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, ColSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_dspmv(CblasColMajor, Uplo,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  /*** RowSymPacked and Upper ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_sspmv(CblasRowMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    cblas_dspmv(CblasRowMajor, CblasUpper,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  /*** RowSymPacked and Uplo ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const float alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<float, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<float, VectFull, Allocator1>& X,
	      const float beta,
	      Vector<float, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_sspmv(CblasRowMajor, Uplo,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void MltAdd(const double alpha,
	      const SeldonUplo& Uplo,
	      const Matrix<double, Prop0, RowSymPacked, Allocator0>& A,
	      const Vector<double, VectFull, Allocator1>& X,
	      const double beta,
	      Vector<double, VectFull, Allocator2>& Y)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, Y, "MltAdd(alpha, uplo, M, X, beta, Y)");
#endif

    cblas_dspmv(CblasRowMajor, Uplo,
		A.GetM(), alpha, A.GetData(),
		X.GetData(), 1, beta, Y.GetData(), 1);
  }


  // MLTADD //
  ////////////



  /////////////////
  // RANK1UPDATE //


  /*** ColMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator1>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop0, ColMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_sger(CblasColMajor, A.GetM(), A.GetN(), alpha, X.GetData(), 1,
	       Y.GetData(), 1, A.GetData(), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator1>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop0, ColMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_dger(CblasColMajor, A.GetM(), A.GetN(), alpha, X.GetData(), 1,
	       Y.GetData(), 1, A.GetData(), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_cgeru(CblasColMajor, A.GetM(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_zgeru(CblasColMajor, A.GetM(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  /*** ColMajor and ConjY ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, ColMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, status, Y, M)", "X.Y' + M");
#endif

    if (ConjY.Conj())
      cblas_cgerc(CblasColMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
    else
      cblas_cgeru(CblasColMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, status, Y, M)", "X.Y' + M");
#endif

    if (ConjY.Conj())
      cblas_zgerc(CblasColMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
    else
      cblas_zgeru(CblasColMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator1>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop0, RowMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_sger(CblasRowMajor, A.GetM(), A.GetN(), alpha, X.GetData(), 1,
	       Y.GetData(), 1, A.GetData(), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator1>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop0, RowMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_dger(CblasRowMajor, A.GetM(), A.GetN(), alpha, X.GetData(), 1,
	       Y.GetData(), 1, A.GetData(), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_cgeru(CblasRowMajor, A.GetM(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, Y, M)", "X.Y' + M");
#endif

    cblas_zgeru(CblasRowMajor, A.GetM(), A.GetN(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  /*** RowMajor and ConjY ***/


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<float> alpha,
		   const Vector<complex<float>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<float>, VectFull, Allocator2>& Y,
		   Matrix<complex<float>, Prop0, RowMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, status, Y, M)", "X.Y' + M");
#endif

    if (ConjY.Conj())
      cblas_cgerc(CblasRowMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
    else
      cblas_cgeru(CblasRowMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1, class Allocator2>
  void Rank1Update(const complex<double> alpha,
		   const Vector<complex<double>, VectFull, Allocator1>& X,
		   const SeldonConjugate& ConjY,
		   const Vector<complex<double>, VectFull, Allocator2>& Y,
		   Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(SeldonNoTrans, A, Y, X,
	     "Rank1Update(alpha, X, status, Y, M)", "X.Y' + M");
#endif

    if (ConjY.Conj())
      cblas_zgerc(CblasRowMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
    else
      cblas_zgeru(CblasRowMajor, A.GetM(), A.GetN(),
		  reinterpret_cast<const void*>(&alpha),
		  reinterpret_cast<const void*>(X.GetData()), 1,
		  reinterpret_cast<const void*>(Y.GetData()), 1,
		  reinterpret_cast<void*>(A.GetData()), A.GetLD());
  }


  /*** ColSymPacked and Upper ***/


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_sspr(CblasColMajor, CblasUpper, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_dspr(CblasColMajor, CblasUpper, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_chpr(CblasColMajor, CblasUpper, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_zhpr(CblasColMajor, CblasUpper, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  /*** ColSymPacked and Uplo ***/


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_sspr(CblasColMajor, Uplo, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_dspr(CblasColMajor, Uplo, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_chpr(CblasColMajor, Uplo, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_zhpr(CblasColMajor, Uplo, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  /*** RowSymPacked and Upper ***/


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_sspr(CblasRowMajor, CblasUpper, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_dspr(CblasRowMajor, CblasUpper, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_chpr(CblasRowMajor, CblasUpper, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, M)", "X.X' + M");
#endif

    cblas_zhpr(CblasRowMajor, CblasUpper, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  /*** RowSymPacked and Uplo ***/


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_sspr(CblasRowMajor, Uplo, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void Rank1Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_dspr(CblasRowMajor, Uplo, A.GetM(), alpha,
	       X.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const float alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_chpr(CblasRowMajor, Uplo, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1>
  void
  Rank1Update(const double alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(A, X, "Rank1Update(alpha, X, uplo, M)", "X.X' + M");
#endif

    cblas_zhpr(CblasRowMajor, Uplo, A.GetM(), alpha,
	       reinterpret_cast<const void*>(X.GetData()), 1,
	       reinterpret_cast<void*>(A.GetData()));
  }


  // RANK1UPDATE //
  /////////////////



  /////////////////
  // RANK2UPDATE //


  /*** ColSymPacked and Upper ***/


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_sspr2(CblasColMajor, CblasUpper, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_dspr2(CblasColMajor, CblasUpper, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_chpr2(CblasColMajor, CblasUpper, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_zhpr2(CblasColMajor, CblasUpper, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  /*** ColSymPacked and Uplo ***/


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_sspr2(CblasColMajor, Uplo, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, ColSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_dspr2(CblasColMajor, Uplo, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_chpr2(CblasColMajor, Uplo, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, ColHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_zhpr2(CblasColMajor, Uplo, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  /*** RowSymPacked and Upper ***/


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_sspr2(CblasRowMajor, CblasUpper, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_dspr2(CblasRowMajor, CblasUpper, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_chpr2(CblasRowMajor, CblasUpper, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, M)", "X.Y' + Y.X' + M");
#endif

    cblas_zhpr2(CblasRowMajor, CblasUpper, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  /*** RowSymPacked and Uplo ***/


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const float alpha,
		   const Vector<float, VectFull, Allocator0>& X,
		   const Vector<float, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<float, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_sspr2(CblasRowMajor, Uplo, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void Rank2Update(const double alpha,
		   const Vector<double, VectFull, Allocator0>& X,
		   const Vector<double, VectFull, Allocator2>& Y,
		   const SeldonUplo& Uplo,
		   Matrix<double, Prop1, RowSymPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_dspr2(CblasRowMajor, Uplo, A.GetM(), alpha,
		X.GetData(), 1, Y.GetData(), 1, A.GetData());
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<float> alpha,
	      const Vector<complex<float>, VectFull, Allocator0>& X,
	      const Vector<complex<float>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<float>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_chpr2(CblasRowMajor, Uplo, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  template <class Allocator0,
	    class Prop1, class Allocator1,
	    class Allocator2>
  void
  Rank2Update(const complex<double> alpha,
	      const Vector<complex<double>, VectFull, Allocator0>& X,
	      const Vector<complex<double>, VectFull, Allocator2>& Y,
	      const SeldonUplo& Uplo,
	      Matrix<complex<double>, Prop1, RowHermPacked, Allocator1>& A)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
    CheckDim(A, Y, "Rank2Update(alpha, X, Y, uplo, M)", "X.Y' + Y.X' + M");
#endif

    cblas_zhpr2(CblasRowMajor, Uplo, A.GetM(),
		reinterpret_cast<const void*>(&alpha),
		reinterpret_cast<const void*>(X.GetData()), 1,
		reinterpret_cast<const void*>(Y.GetData()), 1,
		reinterpret_cast<void*>(A.GetData()));
  }


  // RANK2UPDATE //
  /////////////////



  ///////////
  // SOLVE //


  /*** ColUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_strsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_strsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtrsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctrsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztrsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_strsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, uplo, M, X)");
#endif

    cblas_strsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, uplo, M, X)");
#endif

    cblas_dtrsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, uplo, M, X)");
#endif

    cblas_ctrsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, uplo, M, X)");
#endif

    cblas_ztrsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_strsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztrsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_strsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowUpTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtrsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowUpTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctrsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowUpTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztrsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriang, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_strsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriang ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_strsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowLoTriang, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtrsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), A.GetM(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowLoTriang, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctrsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowLoTriang, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztrsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		A.GetM(), reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColUpTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_stpsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtpsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctpsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztpsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColUpTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_stpsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtpsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctpsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztpsv(CblasColMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_stpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	ColLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** ColLoTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_stpsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, ColLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtpsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, ColLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctpsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	ColLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztpsv(CblasColMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_stpsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtpsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctpsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztpsv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowUpTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_stpsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowUpTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtpsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowUpTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctpsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowUpTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztpsv(CblasRowMajor, CblasUpper, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriangPacked, NoTrans and NonUnit ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_stpsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_dtpsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ctpsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const Matrix<complex<double>, Prop0,
	RowLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(M, X)");
#endif

    cblas_ztpsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  /*** RowLoTriangPacked ***/


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<float, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_stpsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void Solve(const SeldonTranspose& TransA,
	     const SeldonDiag& DiagA,
	     const Matrix<double, Prop0, RowLoTriangPacked, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_dtpsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), A.GetData(), X.GetData(), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<float>, Prop0, RowLoTriangPacked, Allocator0>& A,
	Vector<complex<float>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ctpsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  template <class Prop0, class Allocator0,
	    class Allocator1>
  void
  Solve(const SeldonTranspose& TransA,
	const SeldonDiag& DiagA,
	const Matrix<complex<double>, Prop0,
	RowLoTriangPacked, Allocator0>& A,
	Vector<complex<double>, VectFull, Allocator1>& X)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    CheckDim(A, X, "Solve(status, diag, M, X)");
#endif

    cblas_ztpsv(CblasRowMajor, CblasLower, TransA, DiagA,
		A.GetN(), reinterpret_cast<const void*>(A.GetData()),
		reinterpret_cast<void*>(X.GetData()), 1);
  }


  // SOLVE //
  ///////////


} // namespace Seldon.

#define SELDON_FILE_BLAS_2_CXX
#endif
