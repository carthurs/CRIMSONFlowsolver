// Copyright (C) 2003-2009 Marc Duruflé
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
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
#define SELDON_WITH_BLAS

#include "Seldon.hxx"

using namespace Seldon;

// testing class Vector without Blas

template<class T>
void Fill(Vector<T>& x, Vector<T>& y, T& coef, T& coef2)
{
  for (int i = 0; i < x.GetM(); i++)
    x(i) = 1.3*i + 0.5;

  for (int i = 0; i < y.GetM(); i++)
    y(i) = -0.2*i*i + 0.6*i + 1.3;

  coef = 0.9;
  coef2 = -1.6;
}

template<class T>
void Fill(Vector<complex<T> >& x, Vector<complex<T> >& y, complex<T>& coef, complex<T>& coef2)
{
  for (int i = 0; i < x.GetM(); i++)
    x(i) = complex<T>(1.3*i + 0.5, 0.4*log(0.2+i) -0.3);

  for (int i = 0; i < y.GetM(); i++)
    y(i) = complex<T>(-0.2*i*i + 0.6*i + 1.3, 1.2*cos(0.31*i+0.22));

  coef = complex<T>(0.9, 1.1);
  coef2 = complex<T>(-1.6, 0.7);
}

template<class T>
void TestGivensRotation(T x)
{
  Vector<T> U(3), V(3);
  T a, b;
  Fill(U, V, a, b);

  // givens rotation
  T cos_, sin_;
  DISP(a); DISP(b);
  GenRot(a, b, cos_, sin_);
  DISP(cos_); DISP(sin_);

  cout << "vectors before applyrot " << endl;
  DISP(U); DISP(V);
  ApplyRot(U, V, cos_, sin_);
  cout << "vectors after applyrot " << endl;
  DISP(U); DISP(V);
}


template<class T>
void TestBlas1(T x)
{
  Vector<T> U(3), V(3), W(3);
  T coef, coef2;
  Fill(U, V, coef, coef2);
  W.Zero();

  // swaps two vectors
  cout << "vectors before swap " << endl;
  DISP(U); DISP(V);
  Swap(U, V);
  cout << "vectors after swap " << endl;
  DISP(U); DISP(V);

  // W = V
  cout << "vectors before copy " << endl;
  DISP(V); DISP(W);
  Copy(V, W);
  cout << "vectors after copy " << endl;
  DISP(V); DISP(W);

  cout << "vector before Mlt " << endl;
  DISP(W);
  Mlt(coef, W);
  cout << "vector after Mlt " << endl;
  DISP(W);

  cout << "vectors before Add " << endl;
  DISP(U); DISP(W);
  Add(coef2, U, W);
  cout << "vectors after Add " << endl;
  DISP(U); DISP(W);

  coef = DotProd(U, W);
  cout << "scalar product U.W " << coef << endl;

  coef = Norm1(W);
  cout << "1-norm of W " << coef << endl;

  coef = Norm2(W);
  cout << "2-norm of W " << coef << endl;

  int imax = GetMaxAbsIndex(W);
  cout << "Index where W is maximal " << imax << endl;
}

template<class T>
void TestBlas1(complex<T> x)
{
  Vector<complex<T> > U(3), V(3), W(3);
  complex<T> coef, coef2; T norme;
  Fill(U, V, coef, coef2);
  W.Zero();

  // swaps two vectors
  cout << "vectors before swap " << endl;
  DISP(U); DISP(V);
  Swap(U, V);
  cout << "vectors after swap " << endl;
  DISP(U); DISP(V);

  // W = V
  cout << "vectors before copy " << endl;
  DISP(V); DISP(W);
  Copy(V, W);
  DISP(W);
  cout << "vectors after copy " << endl;
  DISP(V); DISP(W);

  cout << "vector before Mlt " << endl;
  DISP(W);
  Mlt(coef, W);
  cout << "vector after Mlt " << endl;
  DISP(W);

  cout << "vectors before Add " << endl;
  DISP(U); DISP(W);
  Add(coef2, U, W);
  cout << "vectors after Add " << endl;
  DISP(U); DISP(W);

  coef = DotProd(U, W);
  cout << "scalar product U.W " << coef << endl;

  coef = DotProdConj(U, W);
  cout << "scalar product U'.W " << coef << endl;

  norme = Norm1(W);
  cout << "1-norm of W " << norme << endl;

  norme = Norm2(W);
  cout << "2-norm of W " << norme << endl;

  Conjugate(W);
  cout << "conjugate of W " << W << endl;

  int imax = GetMaxAbsIndex(W);
  cout << "Index where W is maximal " << imax << endl;
}

template<class T>
void TestBlas23(T x)
{
  int m = 4, n = 3, k = 5;
  Vector<T> U(m*n), V(n*k), X(k), Y(n);
  T coef, coef2;
  Fill(U, V, coef, coef2);

  Matrix<T> A(m, n), B(n, k), C(m, k);

  for (int i = 0; i < k; i++)
    X(i) = U(i);

  int nb = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A(i, j) = U(nb++);

  nb = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++)
      B(i, j) = V(nb++);

  C.SetIdentity();
  Y.Fill();

  cout << "matrices before MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);
  MltAdd(coef, A, B, coef2, C);
  cout << "matrices after MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);

  cout << "matrices before MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);
  Transpose(A); Transpose(B);
  MltAdd(coef, SeldonTrans, A, SeldonTrans, B, coef2, C);
  cout << "matrices after MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);

  cout << "vectors before MltAdd " << endl;
  DISP(X); DISP(Y);
  Transpose(B);
  MltAdd(coef, B, X, coef2, Y);
  cout << "vectors after MltAdd " << endl;
  DISP(X); DISP(Y);

  cout << "vectors before MltAdd " << endl;
  DISP(X); DISP(Y);
  Transpose(B);
  MltAdd(coef, SeldonTrans, B, X, coef2, Y);
  cout << "vectors after MltAdd " << endl;
  DISP(X); DISP(Y);

  cout << "matrix/vectors before Rank1Update " << endl;
  DISP(B); DISP(X); DISP(Y);
  Rank1Update(coef, X, Y, B);
  cout << "matrix/vectors after Rank1Update " << endl;
  DISP(B); DISP(X); DISP(Y);

  Matrix<T, Symmetric, RowSymPacked> As(n, n);
  nb = 0;
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      As(i, j) = V(nb++);

  cout << "matrix/vectors before Rank1Update " << endl;
  DISP(As); DISP(Y);
  Rank1Update(coef, Y, As);
  cout << "matrix/vectors after Rank1Update " << endl;
  DISP(As); DISP(Y);

  X.Resize(n);
  cout << "matrix/vectors before Rank2Update " << endl;
  DISP(As); DISP(X); DISP(Y);
  Rank2Update(coef, X, Y, As);
  cout << "matrix/vectors after Rank2Update " << endl;
  DISP(As); DISP(X); DISP(Y);

}


template<class T>
void TestBlas23(complex<T> x)
{
  int m = 4, n = 3, k = 5;
  Vector<complex<T> > U(m*n), V(n*k), X(k), Y(n);
  complex<T> coef, coef2;
  Fill(U, V, coef, coef2);

  Matrix<complex<T> > A(m, n), B(n, k), C(m, k);

  for (int i = 0; i < k; i++)
    X(i) = U(i);

  int nb = 0;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A(i, j) = U(nb++);

  nb = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++)
      B(i, j) = V(nb++);

  C.SetIdentity();
  for (int i = 0; i < n; i++)
    Y(i) = V(i);

  cout << "matrices before MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);
  MltAdd(coef, A, B, coef2, C);
  cout << "matrices after MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);

  cout << "matrices before MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);
  Transpose(A); Transpose(B);
  MltAdd(coef, SeldonTrans, A, SeldonTrans, B, coef2, C);
  cout << "matrices after MltAdd " << endl;
  DISP(A); DISP(B); DISP(C);

  cout << "vectors before MltAdd " << endl;
  DISP(X); DISP(Y);
  Transpose(B);
  MltAdd(coef, B, X, coef2, Y);
  cout << "vectors after MltAdd " << endl;
  DISP(X); DISP(Y);

  cout << "vectors before MltAdd " << endl;
  DISP(X); DISP(Y);
  Transpose(B);
  MltAdd(coef, SeldonTrans, B, X, coef2, Y);
  cout << "vectors after MltAdd " << endl;
  DISP(X); DISP(Y);

  cout << "matrix/vectors before Rank1Update " << endl;
  DISP(B); DISP(X); DISP(Y);
  Rank1Update(coef, X, Y, B);
  cout << "matrix/vectors after Rank1Update " << endl;
  DISP(B); DISP(X); DISP(Y);

  Matrix<complex<T>, General, RowHermPacked> As(n, n);
  nb = 0; T norme = 0.5;
  for (int i = 0; i < n; i++)
    {
      As(i, i) = real(V(nb++));
      for (int j = i+1; j < n; j++)
	As(i, j) = V(nb++);
    }

  cout << "matrix/vectors before Rank1Update " << endl;
  DISP(As); DISP(Y);
  Rank1Update(norme, Y, As);
  cout << "matrix/vectors after Rank1Update " << endl;
  DISP(As); DISP(Y);

  X.Resize(n);
  cout << "matrix/vectors before Rank2Update " << endl;
  DISP(As); DISP(X); DISP(Y);
  Rank2Update(norme, X, Y, As);
  cout << "matrix/vectors after Rank2Update " << endl;
  DISP(As); DISP(X); DISP(Y);

}



int main()
{
  cout << "Seldon: compilation test of class Vector with Blas" << endl;

  TestGivensRotation(float(0));

  TestGivensRotation(double(0));

  TestBlas1(float(0));

  TestBlas1(double(0));

  TestBlas1(complex<float>(0, 0));

  TestBlas1(complex<double>(0, 0));

  TestBlas23(float(0));

  TestBlas23(double(0));

  TestBlas23(complex<float>(0));

  TestBlas23(complex<double>(0));

  return 0;

}
