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
#include "matrix/HeterogeneousMatrixCollection.cxx"
using namespace Seldon;


class HeterogeneousMatrixCollectionTest: public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(HeterogeneousMatrixCollectionTest);
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
    M_ = 15;
    N_ = 35;
    Mlocal_max_ = 100;
    Nlocal_max_ = 80;
    mlt();
  }


  void test_add()
  {
    Nloop_ = 10;
    M_ = 2;
    N_ = 2;
    Mlocal_max_ = 3;
    Nlocal_max_ = 3;
    add();
  }


  void test_mlt_add()
  {
    Nloop_ = 10;
    M_ = 10;
    N_ = 10;
    Mlocal_max_ = 30;
    Nlocal_max_ = 30;
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
    bool dense;

    for (int N = 0; N < Nloop_; N++)
      {
	typedef Matrix<real, General, RowMajor, NewAlloc<real> > matrix_dense;
	typedef Matrix<real, General, RowSparse, NewAlloc<real> >
          matrix_sparse;

	Matrix<FloatDouble, General, DenseSparseCollection, NewAlloc<double> >
          A(M_, N_);
	matrix_dense U0;

	matrix_sparse U1;
	Matrix<real, General, ArrayRowSparse, NewAlloc<real> > U1_array;

	int r, s;
	real value;
	real alpha;
	alpha = 2.0;

	for (int k = 0; k < M_; k++)
	  Mlocal(k) = rand() % Mlocal_max_ + 1;
	for (int k = 0; k < N_; k++)
	  Nlocal(k) = rand() % Nlocal_max_ + 1;


	for (int i = 0; i < M_; i++)
	  for (int j = 0; j < N_; j++)
	    {
	      dense = bool(rand() % 2);
	      if (dense)
		{
		  U0.Reallocate(Mlocal(i), Nlocal(j));
		  U0.FillRand();
		  A.SetMatrix(i, j, U0);
		  U0.Nullify();
		}
	      else
		{
		  U1_array.Reallocate(Mlocal(i), Nlocal(j));
		  for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
		    {
		      r = rand() % Mlocal(i);
		      s = rand() % Nlocal(j);
		      value = double(rand());
		      U1_array.AddInteraction(r, s, value);
		    }

		  Copy(U1_array, U1);

		  A.SetMatrix(i, j, U1);
		  U1.Nullify();
		}
	    }

	int M_full = Norm1(Mlocal);
	int N_full = Norm1(Nlocal);

	matrix_dense A_full(M_full, N_full);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    A_full(i, j) = A(i, j);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    CPPUNIT_ASSERT(A(i, j) == A_full(i, j));

	Mlt(alpha, A_full);


	Mlt(alpha, A);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    CPPUNIT_ASSERT(A(i, j) == A_full(i, j));

	A.Deallocate();

      }
  }


  void add()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal(M_), Nlocal(N_);
    bool dense;

    for (int N = 0; N < Nloop_; N++)
      {
	typedef Matrix<real, General, RowMajor, MallocAlloc<real> >
          matrix_dense;
	typedef Matrix<real, General, RowSparse, MallocAlloc<real> >
          matrix_sparse;

	Matrix<FloatDouble, General, DenseSparseCollection,
          MallocAlloc<double> >
	  A(M_, N_), B(M_, N_);
	matrix_dense U0;

	matrix_sparse U1;
	Matrix<real, General, ArrayRowSparse, MallocAlloc<real> > U1_array;

	int r, s;
	real value;

	for (int k = 0; k < M_; k++)
	  Mlocal(k) = rand() % Mlocal_max_ + 1;
	for (int k = 0; k < N_; k++)
	  Nlocal(k) = rand() % Nlocal_max_ + 1;

	for (int i = 0; i < M_; i++)
	  for (int j = 0; j < N_; j++)
	    {
	      dense = bool(rand() % 2);
	      if (dense)
		{
		  U0.Reallocate(Mlocal(i), Nlocal(j));
		  U0.FillRand();
		  A.SetMatrix(i, j, U0);
		  U0.Nullify();

		  U0.Reallocate(Mlocal(i), Nlocal(j));
		  U0.FillRand();
		  B.SetMatrix(i, j, U0);
		  U0.Nullify();
		}
	      else
		{
		  U1_array.Reallocate(Mlocal(i), Nlocal(j));
		  for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
		    {
		      r = rand() % Mlocal(i);
		      s = rand() % Nlocal(j);
		      value = double(rand());
		      U1_array.AddInteraction(r, s, value);
		    }

		  Copy(U1_array, U1);

		  A.SetMatrix(i, j, U1);
		  U1.Nullify();
		  U1_array.Clear();

		  U1_array.Reallocate(Mlocal(i), Nlocal(j));
		  for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
		    {
		      r = rand() % Mlocal(i);
		      s = rand() % Nlocal(j);
		      value = double(rand());
		      U1_array.AddInteraction(r, s, value);
		    }

		  Copy(U1_array, U1);

		  B.SetMatrix(i, j, U1);
		  U1.Nullify();
		}
	    }

	int M_full = Norm1(Mlocal);
	int N_full = Norm1(Nlocal);

	matrix_dense A_full(M_full, N_full), B_full(M_full, N_full);

	real beta = 1.;

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    A_full(i, j) = A(i, j);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    B_full(i, j) = B(i, j);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    CPPUNIT_ASSERT(A(i, j) == A_full(i, j));

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    CPPUNIT_ASSERT(B(i, j) == B_full(i, j));

	Add(beta, B_full, A_full);

	Add(beta, B, A);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    CPPUNIT_ASSERT(A(i, j) == A_full(i, j));

	A.Deallocate();
	B.Deallocate();
      }
  }


  void mlt_add()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal_a(M_), Nlocal_a(N_);
    Vector<int> Mlocal_b(N_), Nlocal_b(M_);
    Vector<int> Mlocal_c(M_), Nlocal_c(M_);
    int tmp;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real, General, RowMajor, MallocAlloc<real> >
            matrix_real;

	  Matrix<FloatDouble, General, DenseSparseCollection,
            MallocAlloc<double> >
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
	      CPPUNIT_ASSERT_DOUBLES_EQUAL(C_Full(i, j), C(i, j),
					   1.e-14 * C_Full(i, j));

	  A.Deallocate();
	  B.Deallocate();
	  C.Deallocate();
	}
      }
  }


  void copy()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal(M_), Nlocal(N_);
    bool dense;

    for (int N = 0; N < Nloop_; N++)
      {
	typedef Matrix<real, General, RowMajor, NewAlloc<real> > matrix_dense;
	typedef Matrix<real, General, RowSparse, NewAlloc<real> >
          matrix_sparse;

	Matrix<FloatDouble, General, DenseSparseCollection, NewAlloc<double> >
          A(M_, N_), B;
	matrix_dense U0;

	matrix_sparse U1;
	Matrix<real, General, ArrayRowSparse, NewAlloc<real> > U1_array;

	int r, s;
	real value;

	for (int k = 0; k < M_; k++)
	  Mlocal(k) = rand() % Mlocal_max_ + 1;
	for (int k = 0; k < N_; k++)
	  Nlocal(k) = rand() % Nlocal_max_ + 1;

	for (int i = 0; i < M_; i++)
	  for (int j = 0; j < N_; j++)
	    {
	      dense = bool(rand() % 2);
	      if (dense)
		{
		  U0.Reallocate(Mlocal(i), Nlocal(j));
		  U0.FillRand();
		  A.SetMatrix(i, j, U0);
		  U0.Nullify();
		}
	      else
		{
		  U1_array.Reallocate(Mlocal(i), Nlocal(j));
		  for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
		    {
		      r = rand() % Mlocal(i);
		      s = rand() % Nlocal(j);
		      value = double(rand());
		      U1_array.AddInteraction(r, s, value);
		    }

		  Copy(U1_array, U1);

		  A.SetMatrix(i, j, U1);
		  U1.Nullify();
		}
	    }

	B.Copy(A);

	int M_full = Norm1(Mlocal);
	int N_full = Norm1(Nlocal);

	for (int i = 0; i < M_full; i++)
	  for (int j = 0; j < N_full; j++)
	    CPPUNIT_ASSERT(A(i, j) == B(i, j));

	A.Deallocate();
	B.Nullify();
      }
  }


  void mat_vect()
  {
    srand(time(NULL));

    typedef double real;
    bool dense;
    Vector<int> Mlocal(M_), Nlocal(N_);
    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense matrix collection.
	  typedef Matrix<real, General, RowMajor, NewAlloc<real> >
            matrix_dense;
	  typedef Matrix<real, General, RowSparse, NewAlloc<real> >
            matrix_sparse;

	  Matrix<FloatDouble, General, DenseSparseCollection,
            NewAlloc<double> > A(M_, N_);

	  real alpha, beta;
	  alpha = rand();
	  beta = rand();

	  matrix_dense U0;
	  matrix_sparse U1;
	  Matrix<real, General, ArrayRowSparse, NewAlloc<real> > U1_array;

	  int r, s;
	  real value;

	  for (int k = 0; k < M_; k++)
	    Mlocal(k) = rand() % Mlocal_max_ + 1;
	  for (int k = 0; k < N_; k++)
	    Nlocal(k) = rand() % Nlocal_max_ + 1;

	  for (int i = 0; i < M_; i++)
	    for (int j = 0; j < N_; j++)
	      {
		dense = bool(rand() % 2);
		if (dense)
		  {
		    U0.Reallocate(Mlocal(i), Nlocal(j));
		    U0.FillRand();
		    A.SetMatrix(i, j, U0);
		    U0.Nullify();
		  }
		else
		  {
		    U1_array.Reallocate(Mlocal(i), Nlocal(j));
		    for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
		      {
			r = rand() % Mlocal(i);
			s = rand() % Nlocal(j);
			value = double(rand());
			U1_array.AddInteraction(r, s, value);
		      }

		    Copy(U1_array, U1);

		    A.SetMatrix(i, j, U1);
		    U1.Nullify();
		  }
	      }

	  int Ma_full, Na_full;
	  Ma_full = Norm1(Mlocal);
	  Na_full = Norm1(Nlocal);

	  matrix_dense A_Full(Ma_full, Na_full);

	  for (int i = 0; i < Ma_full; i++)
	    for (int j = 0; j < Na_full; j++)
	      A_Full(i, j) = A(i, j);

	  Vector<real> X(Na_full), B1(Ma_full), B2;
	  X.FillRand();
	  B1.FillRand();
	  B2.Copy(B1);

	  MltAdd(alpha, A, X, beta, B1);

	  MltAdd(alpha, A_Full, X, beta, B2);

	  for (int i = 0; i < B1.GetM(); i++)
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(B1(i), B2(i),
					 1.e-14 *B1(i));

	  A.Deallocate();

	}
      }
  }


  void write_read()
  {
    srand(time(NULL));

    typedef double real;
    Vector<int> Mlocal(M_), Nlocal(N_);
    bool dense;

    for (int N = 0; N < Nloop_; N++)
      {
	typedef Matrix<real, General, RowMajor, NewAlloc<real> > matrix_dense;
	typedef Matrix<real, General, RowSparse, NewAlloc<real> >
          matrix_sparse;

	Matrix<FloatDouble, General, DenseSparseCollection, NewAlloc<double> >
          A(M_, N_), B;
	matrix_dense U0;

	matrix_sparse U1;
	Matrix<real, General, ArrayRowSparse, NewAlloc<real> > U1_array;

	int r, s;
	real value;

	for (int k = 0; k < M_; k++)
	  Mlocal(k) = rand() % Mlocal_max_ + 1;
	for (int k = 0; k < N_; k++)
	  Nlocal(k) = rand() % Nlocal_max_ + 1;

	for (int i = 0; i < M_; i++)
	  for (int j = 0; j < N_; j++)
	    {
	      // dense = bool(rand() % 2);
	      // The method Write is not defined for RowSparse matrix so we
	      // fill the heterogeneous container with dense matrices.
	      dense = true;
	      if (dense)
		{
		  U0.Reallocate(Mlocal(i), Nlocal(j));
		  U0.FillRand();
		  A.SetMatrix(i, j, U0);
		  U0.Nullify();
		}
	      else
		{
		  U1_array.Reallocate(Mlocal(i), Nlocal(j));
		  for (int l = 0; l < int (Mlocal(i) * Nlocal(j) / 3); l++)
		    {
		      r = rand() % Mlocal(i);
		      s = rand() % Nlocal(j);
		      value = double(rand());
		      U1_array.AddInteraction(r, s, value);
		    }

		  Copy(U1_array, U1);

		  A.SetMatrix(i, j, U1);
		  U1.Nullify();
		}
	    }

	A.Write("test.bin", true);

	B.Read("test.bin");

	for (int i = 0; i < A.GetM(); i++)
	  for (int j = 0; j < A.GetN(); j++)
	    CPPUNIT_ASSERT(A(i, j) == B(i, j));

	A.Deallocate();
	B.Deallocate();

      }
  }

};
