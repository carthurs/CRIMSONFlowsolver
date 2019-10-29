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


#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS

#include "Seldon.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

typedef double Real_wp;
typedef complex<double> Complex_wp;

// We test different iterative solvers on the matrix
// A = S B S^{-1}
//
// where
//
// S = |1 \beta  0   ...  0  |
//     |0   1  \beta ...  0  |
//     |      ...            |
//     |0 ...      1   \beta |
//     |0 ...      0     1   |
//
// and
//
// B = diag(\lambda_1, ... , \lambda_n)


// Class for a new type of matrix.
template<class T>
class BlackBoxMatrix
{
public :
  int n;
  T beta;
  Seldon::Vector<T> lambda;

  BlackBoxMatrix(int n_)
  {
    beta = 0.9;
    n = n_;
    lambda.Reallocate(n);
    for (int i = 0; i < n; i++)
      lambda(i) = i+1;
  }

  BlackBoxMatrix(int n_, const T& beta_)
  {
    beta = beta_; n = n_;
    lambda.Reallocate(n);
    for (int i = 0; i < n; i++)
      lambda(i) = i+1;
  }

  int GetM() const
  {
    return n;
  }

  int GetN() const
  {
    return n;
  }

};


// We need to define MltAdd(alpha, A, X, beta, Y) : Y <- beta Y + alpha A X
template<class T, class Vector>
void MltAdd(const T& alpha, const BlackBoxMatrix<T>& A,
	    const Vector& X, const T& beta, Vector& Y)
{
  if (beta == T(0))
    Y.Fill(T(0));
  else
    Mlt(beta, Y);

  int n = A.GetM();
  Vector C(X);
  // C = S^{-1} C
  for (int i = (n-2); i >= 0; i--)
    C(i) -= A.beta*C(i+1);

  // C = B C
  for (int i=0 ; i<n ; i++)
    C(i) *= A.lambda(i);

  // C = S C
  for (int i=0 ; i<(n-1) ; i++)
    C(i) += A.beta*C(i+1);

  // Y = Y + alpha C
  Add(alpha, C, Y);
}


// and the transpose...
template<class T, class Vector>
void MltAdd(const T& alpha, const class_SeldonTrans& Trans,
	    const BlackBoxMatrix<T>& A,
	    const Vector& X, const T& beta, Vector& Y)
{
  if (beta == T(0))
    Y.Fill(T(0));
  else
    Mlt(beta, Y);

  int n = A.GetM();
  // Transpose of S B S^{-1} is S^{-t} B S^t
  Vector C(X);
  // Y = S^t Y
  for (int i = (n-1); i >= 1; i--)
    C(i) += A.beta*C(i-1);

  // Y = B Y
  for (int i = 0; i < n; i++)
    C(i) *= A.lambda(i);

  // Y = S^{-t} Y
  for (int i = 1; i < n; i++)
    C(i) -= A.beta*C(i-1);

  // Y = Y + alpha C
  Add(alpha, C, Y);
}

// Y = A*X
template<class T, class Vector>
void Mlt(const BlackBoxMatrix<T>& A, const Vector& X, Vector& Y)
{
  Y.Zero();
  MltAdd(T(1), A, X, T(0), Y);
}

// Y = transpose(A)*X
template<class T, class Vector>
void Mlt(const class_SeldonTrans& Trans, const BlackBoxMatrix<T>& A,
	 const Vector& X, Vector& Y)
{
  Y.Zero();
  MltAdd(T(1), Trans, A, X, T(0), Y);
}


int main(int argc, char **argv)
{
  cout.precision(5);
  int restart = 3;
  int n = 5;
  int nb_max_iter = 5*n;
  double stopping_criterion = 1e-6;

  srand(time(NULL));

  // Resolution of unsymmetric real system.
  cout << "Resolution of an unsymmetric real system " << endl << endl;
  {
    double beta = 0.9;
    BlackBoxMatrix<double> A(n, beta);
    DVect b_rhs(n), x_sol(n);
    x_sol.Fill();
    Mlt(A, x_sol, b_rhs);
    DISP(b_rhs);
    x_sol.Zero();
    Preconditioner_Base prec;

    Iteration<double> iter(nb_max_iter, stopping_criterion);
    iter.ShowFullHistory();
    iter.SetRestart(restart);

    // You can use iterative solvers for unsymmetric real systems not using
    // matrix-vector product with transpose matrix.
    x_sol.Zero(); cout<<"BiCgStab"<<endl;
    BiCgStab(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"BiCgStabl"<<endl;
    BiCgStabl(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Cgs"<<endl;
    Cgs(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Gcr"<<endl;
    Gcr(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Gmres"<<endl;
    Gmres(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"QCgs"<<endl;
    QCgs(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"TfQmr"<<endl;
    TfQmr(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Lsqr"<<endl;
    Lsqr(A, x_sol, b_rhs, prec, iter);

    // And iterative solvers using transpose matrix.
    x_sol.Zero(); cout<<"BiCg"<<endl;
    BiCg(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Cgne"<<endl;
    Cgne(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Qmr"<<endl;
    Qmr(A, x_sol, b_rhs, prec, iter);
  }

  cout << endl << endl
       << "Resolution of an unsymmetric complex system " << endl;
  {
    complex<double> beta(0.9, 0.3);
    BlackBoxMatrix<complex<double> > A(n, beta);
    ZVect b_rhs(n), x_sol(n);
    x_sol.Fill();
    Mlt(A, x_sol, b_rhs);
    DISP(b_rhs);
    x_sol.Zero();
    Preconditioner_Base prec;

    Iteration<double> iter(nb_max_iter, stopping_criterion);
    iter.ShowFullHistory();
    iter.SetRestart(restart);

    // You can use iterative solvers for unsymmetric complex systems not using
    // matrix-vector product with transpose matrix.
    x_sol.Zero(); cout<<"BiCgStab"<<endl;
    BiCgStab(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"BiCgStabl"<<endl;
    BiCgStabl(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Cgs"<<endl;
    Cgs(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Gcr"<<endl;
    Gcr(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Gmres"<<endl;
    Gmres(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"QCgs"<<endl;
    QCgs(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"TfQmr"<<endl;
    TfQmr(A, x_sol, b_rhs, prec, iter);

    // And iterative solvers using transpose matrix.
    x_sol.Zero(); cout<<"BiCg"<<endl;
    BiCg(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Cgne"<<endl;
    Cgne(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero(); cout<<"Qmr"<<endl;
    Qmr(A, x_sol, b_rhs, prec, iter);
  }

  // Resolution of symmetric real system.
  cout << "Resolution of a symmetric real system " << endl << endl;
  {
    Matrix<double, Symmetric, RowSymPacked> A(n,n);
    A.FillRand();
    DISP(A);

    DVect b_rhs(n), x_sol(n);
    x_sol.Fill();
    Mlt(A, x_sol, b_rhs);
    DISP(b_rhs);
    x_sol.Zero();
    Preconditioner_Base prec;

    Iteration<double> iter(nb_max_iter, stopping_criterion);
    iter.ShowFullHistory();
    iter.SetRestart(restart);

    x_sol.Zero();
    cout << "BiCgcr" << endl;
    BiCgcr(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero();
    cout << "Cg" << endl;
    Cg(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero();
    cout << "Symmlq" << endl;
    Symmlq(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero();
    cout << "MinRes" << endl;
    MinRes(A, x_sol, b_rhs, prec, iter);
  }

  // Resolution of symmetric complex system.
  cout << "Resolution of a symmetric complex system " << endl << endl;
  {
    Matrix<double, Symmetric, RowSymPacked> Ar(n,n), Ai(n,n);
    Matrix<complex<double>, Symmetric, RowSymPacked> A(n,n);
    Ar.FillRand(); Ai.FillRand();
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
	A(i,j) = complex<double>(Ar(i,j), Ai(i,j));

    DISP(A);

    ZVect b_rhs(n), x_sol(n);
    x_sol.Fill();
    Mlt(A, x_sol, b_rhs);
    DISP(b_rhs);
    x_sol.Zero();
    Preconditioner_Base prec;

    Iteration<double> iter(nb_max_iter, stopping_criterion);
    iter.ShowFullHistory();
    iter.SetRestart(restart);

    x_sol.Zero();
    cout << "BiCgcr" << endl;
    BiCgcr(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero();
    cout << "CoCg" << endl;
    CoCg(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero();
    cout << "Symmlq" << endl;
    Symmlq(A, x_sol, b_rhs, prec, iter);
    x_sol.Zero();
    cout << "MinRes" << endl;
    MinRes(A, x_sol, b_rhs, prec, iter);
  }

  return 0;
}
