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


#ifndef SELDON_FILE_CHOLMOD_CXX

#include "Cholmod.hxx"

namespace Seldon
{

  MatrixCholmod::MatrixCholmod()
  {
    // Sets default parameters.
    cholmod_start(&param_chol);
    Lsparse = NULL;
    n = 0;
  }



  MatrixCholmod::~MatrixCholmod()
  {
    if (n > 0)
      {
        Clear();
        cholmod_finish(&param_chol);
      }
  }


  void MatrixCholmod::HideMessages()
  {
  }


  void MatrixCholmod::ShowMessages()
  {
  }


  void MatrixCholmod::ShowFullHistory()
  {
  }


  void MatrixCholmod::Clear()
  {
    if (n > 0)
      {
        n = 0;
        cholmod_free_factor(&L, &param_chol);
        if (Lsparse != NULL)
          cholmod_free_sparse(&Lsparse, &param_chol);

        Lsparse = NULL;
      }
  }


  template<class Prop, class Storage, class Allocator>
  void MatrixCholmod::
  FactorizeMatrix(Matrix<double, Prop, Storage, Allocator> & mat,
                  bool keep_matrix)
  {
    Clear();

    n = mat.GetM();
    Matrix<double, Symmetric, RowSymSparse, MallocAlloc<double> > Acsc;
    Copy(mat, Acsc);
    if (!keep_matrix)
      mat.Clear();

    // Initialization of sparse matrix.
    cholmod_sparse A;

    A.nrow = n;
    A.ncol = n;
    A.nzmax = Acsc.GetDataSize();
    A.nz = NULL;
    A.p = Acsc.GetPtr();
    A.i = Acsc.GetInd();
    A.x = Acsc.GetData();
    A.z = NULL;
    A.stype = -1;
    A.xtype = CHOLMOD_REAL;
    A.dtype = CHOLMOD_DOUBLE;
    A.sorted = true;
    A.packed = true;
    L = cholmod_analyze(&A, &param_chol);

    // Cholesky factorization.
    cholmod_factorize(&A, L, &param_chol);

    cholmod_change_factor(CHOLMOD_REAL, true, L->is_super,
                          true, L->is_monotonic, L, &param_chol);

    // We convert the factorization to column sparse row format.
    cholmod_factor* B = cholmod_copy_factor(L, &param_chol);
    Lsparse = cholmod_factor_to_sparse(B, &param_chol);

    cholmod_free_factor(&B, &param_chol);
  }


  //! Solves L x = b or L^T x = b.
  template<class Transpose_status, class Allocator>
  void MatrixCholmod::Solve(const Transpose_status& TransA,
                            Vector<double, VectFull, Allocator>& x)
  {
    // Dense right hand side.
    cholmod_dense b_rhs;
    b_rhs.nrow = x.GetM();
    b_rhs.ncol = 1;
    b_rhs.nzmax = b_rhs.nrow;
    b_rhs.d = b_rhs.nrow;
    b_rhs.x = x.GetData();
    b_rhs.z = NULL;
    b_rhs.xtype = CHOLMOD_REAL;
    b_rhs.dtype = CHOLMOD_DOUBLE;

    cholmod_dense* x_sol, *y;
    if (TransA.Trans())
      {
        y = cholmod_solve(CHOLMOD_Lt, L, &b_rhs, &param_chol);
        x_sol = cholmod_solve(CHOLMOD_Pt, L, y, &param_chol);
      }
    else
      {
        y = cholmod_solve(CHOLMOD_P, L, &b_rhs, &param_chol);
        x_sol = cholmod_solve(CHOLMOD_L, L, y, &param_chol);
      }

    double* data = reinterpret_cast<double*>(x_sol->x);
    for (int i = 0; i < x.GetM(); i++)
      x(i) = data[i];

    cholmod_free_dense(&x_sol, &param_chol);
    cholmod_free_dense(&y, &param_chol);
  }


  //! Performs the matrix vector product y = L X or y = L^T X.
  template<class Transpose_status, class Allocator>
  void MatrixCholmod::Mlt(const Transpose_status& TransA,
                          Vector<double, VectFull, Allocator>& X)
  {
    int trans = 1;
    if(TransA.Trans())
      trans = 0;

    Vector<double, VectFull, Allocator> Y = X;

    cholmod_dense Xchol,Ychol;

    Xchol.nrow = X.GetM();
    Xchol.ncol = 1;
    Xchol.nzmax = Xchol.nrow;
    Xchol.d = Xchol.nrow;
    Xchol.x = X.GetData();
    Xchol.z = NULL;
    Xchol.xtype = CHOLMOD_REAL;
    Xchol.dtype = CHOLMOD_DOUBLE;

    Ychol.nrow = X.GetM();
    Ychol.ncol = 1;
    Ychol.nzmax = Ychol.nrow;
    Ychol.d = Ychol.nrow;
    Ychol.x = Y.GetData();
    Ychol.z = NULL;
    Ychol.xtype = CHOLMOD_REAL;
    Ychol.dtype = CHOLMOD_DOUBLE;

    // Conversion from Lsparse to Seldon structure.
    Matrix<double, General, RowSparse, MallocAlloc<double> > Lcsr;
    Lcsr.SetData(n, n, Lsparse->nzmax, reinterpret_cast<double*>(Lsparse->x),
                 reinterpret_cast<int*>(Lsparse->p),
                 reinterpret_cast<int*>(Lsparse->i));

    if (TransA.Trans())
      {
        // Computing L^T P^T x.
        cholmod_dense* x_sol;
        x_sol = cholmod_solve(CHOLMOD_P, L, &Xchol, &param_chol);

        double* data = reinterpret_cast<double*>(x_sol->x);
        for (int i = 0; i < n; i++)
          Y(i) = data[i];

        Seldon::Mlt(Lcsr, Y, X);

        cholmod_free_dense(&x_sol, &param_chol);
      }
    else
      {
        // Computing P L x.
        Seldon::MltAdd(1.0, SeldonTrans, Lcsr, X, 0.0, Y);

        cholmod_dense* x_sol;
        x_sol = cholmod_solve(CHOLMOD_Pt, L, &Ychol, &param_chol);

        double* data = reinterpret_cast<double*>(x_sol->x);
        for (int i = 0; i < X.GetM(); i++)
          X(i) = data[i];

        cholmod_free_dense(&x_sol, &param_chol);
      }

    Lcsr.Nullify();
  }


  template<class T, class Prop, class Storage, class Allocator>
  void GetCholesky(Matrix<T, Prop, Storage, Allocator>& A,
                   MatrixCholmod& mat_chol, bool keep_matrix = false)
  {
    mat_chol.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Allocator, class Transpose_status>
  void
  SolveCholesky(const Transpose_status& TransA,
                MatrixCholmod& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Solve(TransA, x);
  }


  template<class T, class Allocator, class Transpose_status>
  void
  MltCholesky(const Transpose_status& TransA,
              MatrixCholmod& mat_chol, Vector<T, VectFull, Allocator>& x)
  {
    mat_chol.Mlt(TransA, x);
  }

}

#define SELDON_FILE_CHOLMOD_CXX
#endif
