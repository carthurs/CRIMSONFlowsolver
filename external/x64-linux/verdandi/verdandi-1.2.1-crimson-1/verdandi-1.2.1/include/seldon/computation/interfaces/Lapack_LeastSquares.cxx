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


#ifndef SELDON_FILE_LAPACK_LEAST_SQUARES_CXX

/*
  Functions included in this file:

  xGEQRF   (GetQR, GetLQ)
  xGELQF   (GetQR, GetLQ)
  xGEQP3   (GetQR_Pivot)
  xORGQR   (GetQ_FromQR)
  xUNGQR   (GetQ_FromQR)
  xUNMQR   (MltQ_FromQR)
  xORMQR   (MltQ_FromQR)
  xORMQR + xTRSM   (SolveQR)
  ZUNMQR + ZTRSM   (SolveQR)
  xORMLQ + xTRSM   (SolveQR)
  ZUNMLQ + ZTRSM   (SolveQR)
  xTRSM + xORMLQ   (SolveLQ)
  ZTRSM + ZUNMLQ   (SolveLQ)
  xTRSM + xORMQR   (SolveLQ)
  ZTRSM + ZUNMQR   (SolveLQ)
*/

namespace Seldon
{


  ///////////
  // GETQR //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    sgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    dgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    zgeqrf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  /*** RowMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    sgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    dgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    zgelqf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  // GETQR //
  ///////////


  /////////////////
  // GETQR_PIVOT //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQR_Pivot(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   Vector<int>& ipivot, LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 4 * max(m, n);
    ipivot.Fill(0);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    dgeqp3_(&m, &n, A.GetData(), &m, ipivot.GetData(), tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  // GETQR_PIVOT //
  /////////////////


  /////////////////
  // GETQ_FROMQR //


  /*** ColMajor ***/

  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   Vector<double, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    dorgqr_(&m, &m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetQ_FromQR(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = 2 * max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    zungqr_(&m, &m, &n, A.GetDataVoid(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2, class Side, class Trans>
  void MltQ_FromQR(const Side& side, const Trans& trans,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
		   Vector<complex<double>, VectFull, Allocator1>& tau,
		   Matrix<complex<double>, Prop0, ColMajor, Allocator2>& C,
		   LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    char side_ = side.Char();
    char trans_ = trans.Char();
    int k = m;
    if (side_ == 'R')
      k = n;

    zunmqr_(&side, &trans, &m, &n, &k, A.GetDataVoid(), &m, tau.GetDataVoid(),
	    C.GetDataVoid(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  // GETQ_FROMQR //
  /////////////////


  ///////////
  // GETLQ //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<float, Prop0, ColMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    sgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<double, Prop0, ColMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    dgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    zgelqf_(&m, &n, A.GetData(), &m, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  /*** RowMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<float, Prop0, RowMajor, Allocator0>& A,
	     Vector<float, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<float, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization QR of A^t.
    sgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<double, Prop0, RowMajor, Allocator0>& A,
	     Vector<double, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    dgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  template<class Prop0, class Allocator0,
	   class Allocator1>
  void GetLQ(Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	     Vector<complex<double>, VectFull, Allocator1>& tau,
	     LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int lwork = max(m,n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    tau.Reallocate(min(m, n));
    // Factorization LQ of A^t.
    zgeqrf_(&n, &m, A.GetData(), &n, tau.GetData(),
	    work.GetData(), &lwork, &info.GetInfoRef());
  }


  // GETLQ //
  ///////////


  /////////////////
  // MLTQ_FROMQR //


  /*** ColMajor ***/


  template<class Prop0, class Allocator0,
	   class Allocator1, class Allocator2, class IsTranspose>
  void MltQ_FromQR(Matrix<double, Prop0, ColMajor, Allocator0>& A,
		   const IsTranspose& trans,
		   Vector<double, VectFull, Allocator1>& tau,
		   Vector<double, VectFull, Allocator2>& b,
		   LapackInfo& info = lapack_info)
  {
    int m = b.GetM();
    int n = 1;
    int k = tau.GetM();
    int lwork = max(m,n);
    Vector<double, VectFull, Allocator1> work(lwork);
    char side('L');
    char trans_(trans);
    dormqr_(&side, &trans_, &m, &n, &k, A.GetData(), &m, tau.GetData(),
	    b.GetData(), &m, work.GetData(), &lwork,
	    &info.GetInfoRef());
  }


  // MLTQ_FROMQR //
  /////////////////


  /////////////
  // SOLVEQR //


  /*** ColMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    sormqr_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Then solves R x = Q^t b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    dormqr_(&side, &trans, &lwork, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &lwork, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Then solves R x = Q^t b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Computes Q^t b.
    zunmqr_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Then solves R x = Q^t b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
		CblasNonUnit, b.GetM(), nrhs,
		&alpha, A.GetData(), A.GetM(), b.GetData(), b.GetM());
  }


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Computes Q b.
    sormlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Solves L^t y = b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, n, nrhs,
		alpha, A.GetData(), n, b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Computes Q b.
    dormlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Solves L^t y = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, n, nrhs,
		alpha, A.GetData(), n, b.GetData(), b.GetM());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveQR(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Computes Q b.
    zunmlq_(&side, &trans, &m, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &m, work.GetData(), &lwork, &info.GetInfoRef());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Solves L^t y = b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans,
		CblasNonUnit, n, nrhs,
		&alpha, A.GetData(), n, b.GetData(), b.GetM());
  }


  // SOLVEQR //
  /////////////


  /////////////
  // SOLVELQ //


  /*** ColMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<float, Prop0, ColMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Solves L y = b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		alpha, A.GetData(), m, b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Computes Q^t b.
    sormlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<double, Prop0, ColMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('T');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Solves L y = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		alpha, A.GetData(), m, b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Computes Q^t b.
    dormlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<double>, Prop0, ColMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Solve L y = b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans,
		CblasNonUnit, m, nrhs,
		&alpha, A.GetData(), m, b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Computes Q^t.
    zunmlq_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &m, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  /*** RowMajor ***/


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<float, Prop0, RowMajor, Allocator0>& A,
	       const Vector<float, VectFull, Allocator1>& tau,
	       Vector<float, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<float, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    float alpha(1);
    cblas_strsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit, b.GetM(), nrhs,
		alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Multiplies by Q.
    sormqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<double, Prop0, RowMajor, Allocator0>& A,
	       const Vector<double, VectFull, Allocator1>& tau,
	       Vector<double, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('N');
    int lwork = max(m, n);
    Vector<double, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    double alpha(1);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit, b.GetM(), nrhs,
		alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Multiplies by Q.
    dormqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  template <class Prop0, class Allocator0,
	    class Allocator1,class Allocator2>
  void SolveLQ(const Matrix<complex<double>, Prop0, RowMajor, Allocator0>& A,
	       const Vector<complex<double>, VectFull, Allocator1>& tau,
	       Vector<complex<double>, VectFull, Allocator2>& b,
	       LapackInfo& info = lapack_info)
  {
    int m = A.GetM();
    int n = A.GetN();
    int k = tau.GetM();
    int nrhs = 1, nb = b.GetM();
    char side('L');
    char trans('C');
    int lwork = max(m, n);
    Vector<complex<double>, VectFull, Allocator1> work(lwork);
    // Solves R^t x = b.
    complex<double> alpha(1);
    cblas_ztrsm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans,
		CblasNonUnit, b.GetM(), nrhs,
		&alpha, A.GetData(), A.GetN(), b.GetData(), b.GetM());

    b.Resize(n);
    for (int i = nb; i < n; i++)
      b(i) = 0;

    // Computes Q b.
    zunmqr_(&side, &trans, &n, &nrhs, &k, A.GetData(),
	    &n, tau.GetData(), b.GetData(),
	    &n, work.GetData(), &lwork, &info.GetInfoRef());
  }


  // SOLVELQ //
  /////////////


} // end namespace

#define SELDON_FILE_LAPACK_LEAST_SQUARES_CXX
#endif

