// Copyright (C) 2010, INRIA
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


#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <ctime>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
#include "matrix/MatrixCollection.cxx"
using namespace Seldon;


class MatrixCollectionTest: public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(MatrixCollectionTest);
  CPPUNIT_TEST(test_mlt);
  CPPUNIT_TEST(test_add);
  CPPUNIT_TEST(test_mlt_add);
  CPPUNIT_TEST(test_copy);
  CPPUNIT_TEST(test_mat_vect);
  CPPUNIT_TEST(test_write_read);
  CPPUNIT_TEST_SUITE_END();

protected:
  int Nloop_;
  int M_;
  int N_;
  int Mlocal_max_;
  int Nlocal_max_;


public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_mlt()
  {
    Nloop_ = 10;
    M_ = 2;
    N_ = 2;
    Mlocal_max_ = 3;
    Nlocal_max_ = 3;
    mlt();
  }


  void test_add()
  {
    Nloop_ = 10;
    M_ = 15;
    N_ = 35;
    Mlocal_max_ = 100;
    Nlocal_max_ = 80;
    add();
  }


  void test_mlt_add()
  {
    Nloop_ = 10;
    M_ = 15;
    N_ = 35;
    Mlocal_max_ = 100;
    Nlocal_max_ = 80;
    mlt_add();
  }


  void test_copy()
  {
    Nloop_ = 10;
    M_ = 15;
    N_ = 35;
    Mlocal_max_ = 100;
    Nlocal_max_ = 80;
    copy();
  }


  void test_mat_vect()
  {
    Nloop_ = 10;
    M_ = 15;
    N_ = 35;
    Mlocal_max_ = 100;
    Nlocal_max_ = 80;
    mat_vect();
  }


  void test_write_read()
  {
    Nloop_ = 10;
    M_ = 15;
    N_ = 35;
    Mlocal_max_ = 100;
    Nlocal_max_ = 80;
    write_read();
  }


  void mlt()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal(M_), Nlocal(N_);
    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection> mlt(M_, N_),
            A(M_, N_);
	  matrix_real U, W;
	  real alpha;
	  alpha = rand();

	  for (int k = 0; k < M_; k++)
	    Mlocal(k) = rand() % Mlocal_max_ + 1;

	  for (int k = 0; k < N_; k++)
	    Nlocal(k) = rand() % Nlocal_max_ + 1;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal(i), Nlocal(j));
		U.FillRand();
		A.SetMatrix(i, j, U);

		W.Reallocate(Mlocal(i), Nlocal(j));
		W.Copy(U);

		Mlt(alpha, W);

		mlt.SetMatrix(i, j, W);

		U.Nullify();
		W.Nullify();
	      }

	  Mlt(alpha, A);

	  for (int i = 0; i < Norm1(Mlocal); i++)
	    for (int j = 0; j < Norm1(Nlocal); j++)
	      CPPUNIT_ASSERT(A(i, j) == mlt(i, j));

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      A.GetMatrix(i, j).Clear();

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      mlt.GetMatrix(i, j).Clear();
	}

	{  // For array row sparse matrix collection.
	  typedef Matrix<double, General, RowSparse> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection> A(M_, N_);
	  matrix_real U;
	  Matrix<double, General, ArrayRowSparse> U_array;
	  int r, s;
	  real value;
	  real alpha;
	  alpha = rand();

	  for (int k = 0; k < M_; k++)
	    Mlocal(k) = rand() % Mlocal_max_ + 1;

	  for (int k = 0; k < N_; k++)
	    Nlocal(k) = rand() % Nlocal_max_ + 1;

	  int M_full = Norm1(Mlocal);
	  int N_full = Norm1(Nlocal);

	  Matrix<real> A_full(M_full, N_full);

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
              {
                U_array.Reallocate(Mlocal(i), Nlocal(j));

                for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
                  {
                    r = rand() % Mlocal(i);
                    s = rand() % Nlocal(j);
                    value = double(rand());
                    U_array.AddInteraction(r, s, value);
                  }

                Copy(U_array, U);

                A.SetMatrix(i, j, U);

                U.Nullify();
              }

	  for (int i = 0; i < M_full; i++)
	    for (int j = 0; j < N_full; j++)
	      A_full(i, j) = A(i, j);

	  Mlt(alpha, A);
	  Mlt(alpha, A_full);

	  for (int i = 0; i < Norm1(Mlocal); i++)
	    for (int j = 0; j < Norm1(Nlocal); j++)
	      CPPUNIT_ASSERT(A(i, j) == A_full(i, j));

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      A.GetMatrix(i, j).Clear();

	}
      }
  }


  void add()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal_a(M_), Nlocal_a(N_);

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection>
	    A(M_, N_), B(M_, N_);

	  matrix_real U;
	  real alpha;
	  alpha = rand();

	  for (int k = 0; k < M_; k++)
	    Mlocal_a(k) = rand() % Mlocal_max_ + 1;

	  for (int k = 0; k < N_; k++)
	    Nlocal_a(k) = rand() % Nlocal_max_ + 1;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal_a(i), Nlocal_a(j));
		U.FillRand();
		A.SetMatrix(i, j, U);
		U.Nullify();
	      }

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal_a(i), Nlocal_a(j));
		U.FillRand();
		B.SetMatrix(i, j, U);
		U.Nullify();
	      }

	  int Ma_full, Na_full;
	  Ma_full = Norm1(Mlocal_a);
	  Na_full = Norm1(Nlocal_a);

	  matrix_real A_Full(Ma_full, Na_full), B_Full(Ma_full, Na_full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      A_Full(i, j) = A(i, j);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      B_Full(i, j) = B(i, j);

	  Add(alpha, A, B);

	  Add(alpha, A_Full, B_Full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      CPPUNIT_ASSERT(B_Full(i, j) == B(i, j));

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      A.GetMatrix(i, j).Clear();

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      B.GetMatrix(i, j).Clear();
	}

	{   // For array row sparse matrix collection.
	  typedef Matrix<double, General, RowSparse> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection> A(M_, N_),
            B(M_, N_);
	  matrix_real U;
	  Matrix<double, General, ArrayRowSparse> U_array;
	  int r, s;
	  real value;
	  real alpha;
	  alpha = rand();

	  for (int k = 0; k < M_; k++)
	    Mlocal_a(k) = rand() % Mlocal_max_ + 1;

	  for (int k = 0; k < N_; k++)
	    Nlocal_a(k) = rand() % Nlocal_max_ + 1;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
              {
                U_array.Reallocate(Mlocal_a(i), Nlocal_a(j));

                for (int l = 0; l < int (Mlocal_a(i) * Nlocal_a(j) / 3); l++)
                  {
                    r = rand() % Mlocal_a(i);
                    s = rand() % Nlocal_a(j);
                    value = double(rand());
                    U_array.AddInteraction(r, s, value);
                  }

                Copy(U_array, U);
                A.SetMatrix(i, j, U);
                U.Nullify();


                U_array.Reallocate(Mlocal_a(i), Nlocal_a(j));

                for (int l = 0; l < int (Mlocal_a(i) * Nlocal_a(j) / 3); l++)
                  {
                    r = rand() % Mlocal_a(i);
                    s = rand() % Nlocal_a(j);
                    value = double(rand());
                    U_array.AddInteraction(r, s, value);
                  }

                Copy(U_array, U);
                B.SetMatrix(i, j, U);
                U.Nullify();
              }

	  int Ma_full, Na_full;
	  Ma_full = Norm1(Mlocal_a);
	  Na_full = Norm1(Nlocal_a);

	  Matrix<double> A_Full(Ma_full, Na_full), B_Full(Ma_full, Na_full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      A_Full(i, j) = A(i, j);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      B_Full(i, j) = B(i, j);

	  Add(alpha, A, B);

	  Add(alpha, A_Full, B_Full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      CPPUNIT_ASSERT(B_Full(i, j) == B(i, j));

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      A.GetMatrix(i, j).Clear();

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      B.GetMatrix(i, j).Clear();
	}
      }
  }

  void mlt_add()
  {
    srand(time(NULL));

    typedef int real;
    Vector<int> Mlocal_a(M_), Nlocal_a(N_);
    Vector<int> Mlocal_b(N_), Nlocal_b(M_);
    Vector<int> Mlocal_c(M_), Nlocal_c(M_);
    int tmp;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection>
	    A(M_, N_), B(N_, M_), C(M_, M_);

	  matrix_real U, V, W;
	  real alpha, beta;
	  alpha = rand();
	  beta = rand();

	  for (int k = 0; k < M_; k++)
	    Mlocal_a(k) = rand() % Mlocal_max_ + 1;

	  tmp = rand() % Nlocal_max_ + 1;
	  for (int k = 0; k < N_; k++)
	    Nlocal_a(k) = tmp;

	  for (int k = 0; k < N_; k++)
	    Mlocal_b(k) = tmp;

	  for (int k = 0; k < M_; k++)
	    Nlocal_b(k) = rand() % Nlocal_max_ + 1;

	  for (int k = 0; k < M_; k++)
	    Mlocal_c(k) =  Mlocal_a(k);

	  for (int k = 0; k < M_; k++)
	    Nlocal_c(k) =  Nlocal_b(k);

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal_a(i), Nlocal_a(j));
		U.FillRand();
		A.SetMatrix(i, j, U);
		U.Nullify();
	      }

	  for (int i = 0; i < N_; i++)
	    for (int j = 0; j < M_; j++)
	      {
		V.Reallocate(Mlocal_b(i), Nlocal_b(j));
		V.FillRand();
		B.SetMatrix(i, j, V);
		V.Nullify();
	      }

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < M_; j++)
	      {
		W.Reallocate(Mlocal_c(i), Nlocal_c(j));
		W.FillRand();
		C.SetMatrix(i, j, W);
		W.Nullify();
	      }

	  int Ma_full, Na_full, Mb_full, Nb_full, Mc_full, Nc_full;
	  Ma_full = Norm1(Mlocal_a);
	  Na_full = Norm1(Nlocal_a);
	  Mb_full = Norm1(Mlocal_b);
	  Nb_full = Norm1(Nlocal_b);
	  Mc_full = Norm1(Mlocal_c);
	  Nc_full = Norm1(Nlocal_c);

	  matrix_real A_Full(Ma_full, Na_full), B_Full(Mb_full, Nb_full),
	    C_Full(Mc_full, Nc_full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      A_Full(i, j) = A(i, j);

	  for (int i = 0; i < Mb_full; i++)
	    for (int j = 0; j < Nb_full; j++)
	      B_Full(i, j) = B(i, j);

	  for (int i = 0; i < Mc_full; i++)
	    for (int j = 0; j < Nc_full; j++)
	      C_Full(i, j) = C(i, j);


	  MltAdd(alpha, A, B, beta, C);

	  MltAdd(alpha, A_Full, B_Full, beta, C_Full);

	  for (int i = 0; i < Norm1(Mlocal_c); i++)
	    for (int j = 0; j < Norm1(Nlocal_c); j++)
	      CPPUNIT_ASSERT(C_Full(i, j) == C(i, j));

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      A.GetMatrix(i, j).Clear();

	  for (int i = 0; i < N_; i++)
	    for (int j = 0; j < M_; j++)
	      B.GetMatrix(i, j).Clear();

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < M_; j++)
	      C.GetMatrix(i, j).Clear();
	}
      }
  }


  void copy()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal(M_), Nlocal(N_);

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection> A(M_, N_),
            B(M_, N_);
	  matrix_real U;

	  for (int k = 0; k < M_; k++)
	    Mlocal(k) = rand() % Mlocal_max_ + 1;

	  for (int k = 0; k < N_; k++)
	    Nlocal(k) = rand() % Nlocal_max_ + 1;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal(i), Nlocal(j));
		U.FillRand();
		A.SetMatrix(i, j, U);
		U.Nullify();
	      }

	  B.Copy(A);

	  for (int i = 0; i < Norm1(Mlocal); i++)
	    for (int j = 0; j < Norm1(Nlocal); j++)
	      CPPUNIT_ASSERT(A(i, j) == B(i, j));

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      A.GetMatrix(i, j).Clear();
	}
      }
  }


  void mat_vect()
  {
    srand(time(NULL));

    typedef int real;
    Vector<int> Mlocal(M_), Nlocal(N_);

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real> matrix_real;
	  typedef Vector<real> vector_real;

	  Matrix<matrix_real, General, RowMajorCollection> A(M_, N_);

	  matrix_real U;
	  real alpha, beta;
	  alpha = rand();
	  beta = rand();

	  int tmp;

	  for (int k = 0; k < M_; k++)
	    Mlocal(k) = rand() % Mlocal_max_ + 1;

	  tmp = rand() % Nlocal_max_ + 1;
	  for (int k = 0; k < N_; k++)
	    Nlocal(k) = tmp;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal(i), Nlocal(j));
		U.FillRand();
		A.SetMatrix(i, j, U);
		U.Nullify();
	      }

	  int Ma_full, Na_full;
	  Ma_full = Norm1(Mlocal);
	  Na_full = Norm1(Nlocal);


	  matrix_real A_Full(Ma_full, Na_full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      A_Full(i, j) = A(i, j);

	  Vector<vector_real, Collection> B, X;

	  vector_real V, W, B_Full, X_Full;
	  for (int k = 0; k < N_; k++)
	    {

	      V.Reallocate(Nlocal(k));
	      V.FillRand();
	      X.AddVector(V);

	      W.Copy(V);
	      X_Full.PushBack(W);

	      V.Nullify();
	    }

	  for (int k = 0; k < M_; k++)
	    {

	      V.Reallocate(Mlocal(k));
	      V.FillRand();
	      B.AddVector(V);

	      W.Copy(V);
	      B_Full.PushBack(W);

	      V.Nullify();
	    }


	  MltAdd(alpha, A, X, beta, B);

	  MltAdd(alpha, A_Full, X_Full, beta, B_Full);

	  for (int i = 0; i < B.GetM(); i++)
	    CPPUNIT_ASSERT(B_Full(i) == B(i));

	  A.Deallocate();

	  X.Deallocate();
	  B.Deallocate();

	  X_Full.Clear();
	  B_Full.Clear();
	  A_Full.Clear();
	}
      }
  }


  void write_read()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal(M_), Nlocal(N_);

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real> matrix_real;

	  Matrix<matrix_real, General, RowMajorCollection> A(M_, N_),
            B(M_, N_);
	  matrix_real U;

	  for (int k = 0; k < M_; k++)
	    Mlocal(k) = rand() % Mlocal_max_ + 1;

	  for (int k = 0; k < N_; k++)
	    Nlocal(k) = rand() % Nlocal_max_ + 1;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		U.Reallocate(Mlocal(i), Nlocal(j));
		U.FillRand();
		A.SetMatrix(i, j, U);
		U.Nullify();
	      }

	  A.Write("test.bin", true);

	  B.Read("test.bin");

	  for (int i = 0; i < Norm1(Mlocal); i++)
	    for (int j = 0; j < Norm1(Nlocal); j++)
	      CPPUNIT_ASSERT(A(i, j) == B(i, j));
	  A.Deallocate();
	  B.Deallocate();
	}
      }
  }

};
