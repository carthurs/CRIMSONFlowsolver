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


#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <ctime>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
#include "vector/HeterogeneousCollection.cxx"
using namespace Seldon;


class HeterogeneousCollectionTest: public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(HeterogeneousCollectionTest);
  CPPUNIT_TEST(test_mlt);
  CPPUNIT_TEST(test_add);
  CPPUNIT_TEST(test_copy);
  CPPUNIT_TEST(test_dot_product);
  CPPUNIT_TEST(test_mlt_add);
  CPPUNIT_TEST(test_write_read);
  CPPUNIT_TEST(test_label);
  CPPUNIT_TEST_SUITE_END();

protected:
  int Nloop_;
  int Nvector_;
  int Nsub_vector_max_;
  int m_;

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
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    mlt();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    mlt();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    mlt();
  }


  void test_add()
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    add();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    add();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    add();
  }


  void test_copy()
  {
    Nloop_ = 1;
    Nvector_ = 10;
    Nsub_vector_max_ = 100;
    copy();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    copy();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    copy();
  }


  void test_dot_product()
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    dot_product();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    dot_product();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    dot_product();
  }


  void test_mlt_add()
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    m_ = 50;
    mlt_add();

    Nloop_ = 10;
    Nvector_ = 10;
    Nsub_vector_max_ = 10000;
    m_ = 10;
    mlt_add();

    Nloop_ = 10;
    Nvector_ = 100;
    Nsub_vector_max_ = 100;
    m_ = 100;
    mlt_add();
  }


  void test_write_read()
  {
    Nloop_ = 1;
    Nvector_ = 20;
    Nsub_vector_max_ = 10;
    m_ = 50;
    write_read();
  }


  void test_label()
  {
    Nloop_ = 10;
    Nvector_ = 1000;
    Nsub_vector_max_ = 10;
    label();
  }


  void mlt()
  {
    srand(time(NULL));

    int length;

    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > mlt;
    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A;

    for (int N = 0; N < Nloop_; N++)
      {
	double alpha;
	alpha = rand();

	typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;
	vector_float_dense X0, Y0;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;

	    X0.Reallocate(length);
	    X0.FillRand();

	    A.AddVector(X0);

	    Y0.Reallocate(length);
	    Y0.Copy(X0);

	    Mlt(alpha, Y0);

	    mlt.AddVector(Y0);

	    X0.Nullify();
	    Y0.Nullify();
	  }

	typedef Vector<float, VectSparse, NewAlloc<float> >
	  vector_float_sparse;
	vector_float_sparse X1, Y1;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;

	    X1.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X1.Index(l) = rand() % length;
		X1.Value(l) = rand();
	      }
	    X1.Assemble();

	    Y1.Reallocate(length);
	    Y1.Copy(X1);

	    A.AddVector(X1);

	    Mlt(alpha, Y1);

	    mlt.AddVector(Y1);

	    X1.Nullify();
	    Y1.Nullify();
	  }

	typedef Vector<double, VectFull, NewAlloc<double> > vector_double_dense;
	vector_double_dense X2, Y2;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;

	    X2.Reallocate(length);
	    X2.FillRand();

	    A.AddVector(X2);

	    Y2.Reallocate(length);
	    Y2.Copy(X2);

	    Mlt(alpha, Y2);

	    mlt.AddVector(Y2);

	    X2.Nullify();
	    Y2.Nullify();
	  }

	typedef Vector<double, VectSparse, NewAlloc<double> >
	  vector_double_sparse;
	vector_double_sparse X3, Y3;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;

	    X3.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X3.Index(l) = rand() % length;
		X3.Value(l) = rand();
	      }
	    X3.Assemble();

	    Y3.Reallocate(length);
	    Y3.Copy(X3);

	    A.AddVector(X3);

	    Mlt(alpha, Y3);

	    mlt.AddVector(Y3);

	    X3.Nullify();
	    Y3.Nullify();
	  }

	Mlt(alpha, A);

	for (int j = 0; j < A.GetFloatDense().GetNvector(); j++)
	  for (int l = 0; l < A.GetFloatDense().GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetFloatDense().GetVector(j)(l) ==
			   mlt.GetFloatDense().GetVector(j)(l));

	for (int j = 0; j < A.GetDoubleDense().GetNvector(); j++)
	  for (int l = 0; l < A.GetDoubleDense().GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetDoubleDense().GetVector(j)(l) ==
			   mlt.GetDoubleDense().GetVector(j)(l));

	for (int j = 0; j < A.GetFloatSparse().GetNvector(); j++)
	  for (int l = 0; l < A.GetFloatSparse().GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Index(l) ==
			       mlt.GetFloatSparse().GetVector(j).Index(l));
		CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Value(l) ==
			       mlt.GetFloatSparse().GetVector(j).Value(l));
	      }

	for (int j = 0; j < A.GetDoubleSparse().GetNvector(); j++)
	  for (int l = 0; l < A.GetDoubleSparse().GetVectorLength()(j); l++)
	    {
	      CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Index(l) ==
			     mlt.GetDoubleSparse().GetVector(j).Index(l));
	      CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Value(l) ==
			     mlt.GetDoubleSparse().GetVector(j).Value(l));
	    }

	A.Deallocate();
	mlt.Deallocate();
      }
  }


  void add()
  {

    srand(time(NULL));

    int length;

    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > sum;
    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A, B;
    Vector<float, VectFull, NewAlloc<float> > C;

    for (int N = 0; N < Nloop_; N++)
      {
	// Added two collections.
	{
	  typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;
	  vector_float_dense X0, Y0, V0, W0;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      X0.Reallocate(length);
	      X0.FillRand();
	      A.AddVector(X0);

	      Y0.Reallocate(length);
	      Y0.Copy(X0);

	      V0.Reallocate(length);
	      V0.FillRand();
	      B.AddVector(V0);

	      W0.Reallocate(length);
	      W0.Copy(V0);

	      Add(1.0, Y0, W0);

	      sum.AddVector(W0);

	      X0.Nullify();
	      V0.Nullify();
	      W0.Nullify();
	    }

	  typedef Vector<float, VectSparse, NewAlloc<float> > vector_float_sparse;
	  vector_float_sparse X1, Y1, V1, W1;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      X1.Reallocate(length);
	      for (int l = 0; l < length; l++)
	      {
		X1.Index(l) = rand() % length;
		X1.Value(l) = rand();
	      }
	      X1.Assemble();
	      A.AddVector(X1);

	      Y1.Reallocate(length);
	      Y1.Copy(X1);

	      V1.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  V1.Index(l) = rand() % length;
		  V1.Value(l) = rand();
		}
	      V1.Assemble();
	      B.AddVector(V1);

	      W1.Reallocate(length);
	      W1.Copy(V1);

	      Add(1.0, Y1, W1);

	      sum.AddVector(W1);

	      X1.Nullify();
	      V1.Nullify();
	      W1.Nullify();
	    }

	  typedef Vector<double, VectFull, NewAlloc<double> > vector_double_dense;
	  vector_double_dense X2, Y2, V2, W2;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      X2.Reallocate(length);
	      X2.FillRand();
	      A.AddVector(X2);

	      Y2.Reallocate(length);
	      Y2.Copy(X2);

	      V2.Reallocate(length);
	      V2.FillRand();
	      B.AddVector(V2);

	      W2.Reallocate(length);
	      W2.Copy(V2);

	      Add(1.0, Y2, W2);

	      sum.AddVector(W2);

	      X2.Nullify();
	      V2.Nullify();
	      W2.Nullify();
	    }

	  typedef Vector<double, VectSparse, NewAlloc<double> > vector_double_sparse;
	  vector_double_sparse X3, Y3, V3, W3;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      X3.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  X3.Index(l) = rand() % length;
		  X3.Value(l) = rand();
		}
	      X3.Assemble();
	      A.AddVector(X3);

	      Y3.Reallocate(length);
	      Y3.Copy(X3);

	      V3.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  V3.Index(l) = rand() % length;
		  V3.Value(l) = rand();
		}
	      V3.Assemble();
	      B.AddVector(V3);

	      W3.Reallocate(length);
	      W3.Copy(V3);

	      Add(1.0, Y3, W3);

	      sum.AddVector(W3);

	      X3.Nullify();
	      V3.Nullify();
	      W3.Nullify();
	    }

	  Add(1.0, B, A);

	  for (int j = 0; j < A.GetFloatDense().GetNvector(); j++)
	    for (int l = 0; l < A.GetFloatDense().GetVectorLength()(j); l++)
	      CPPUNIT_ASSERT(A.GetFloatDense().GetVector(j)(l) ==
			     sum.GetFloatDense().GetVector(j)(l));

	  for (int j = 0; j < A.GetDoubleDense().GetNvector(); j++)
	    for (int l = 0; l < A.GetDoubleDense().GetVectorLength()(j); l++)
	      CPPUNIT_ASSERT(A.GetDoubleDense().GetVector(j)(l) ==
			     sum.GetDoubleDense().GetVector(j)(l));

	  for (int j = 0; j < A.GetFloatSparse().GetNvector(); j++)
	    for (int l = 0; l < A.GetFloatSparse().GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Index(l) ==
			       sum.GetFloatSparse().GetVector(j).Index(l));
		CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Value(l) ==
			       sum.GetFloatSparse().GetVector(j).Value(l));
	      }

	  for (int j = 0; j < A.GetDoubleSparse().GetNvector(); j++)
	    for (int l = 0; l < A.GetDoubleSparse().GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Index(l) ==
			       sum.GetDoubleSparse().GetVector(j).Index(l));
		CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Value(l) ==
			       sum.GetDoubleSparse().GetVector(j).Value(l));
	      }

	  A.Deallocate();
	  B.Deallocate();
	  sum.Deallocate();
	}

	// Added a collection to a dense vector.
	{
	  typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;
	  vector_float_dense X0, Y0, V0, W0;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      X0.Reallocate(length);
	      X0.FillRand();
	      A.AddVector(X0);

	      Y0.Reallocate(length);
	      Y0.Copy(X0);

	      V0.Reallocate(length);
	      V0.FillRand();

	      W0.Reallocate(length);
	      W0.Copy(V0);

	      C.PushBack(V0);

	      Add(1.0, Y0, W0);

	      sum.AddVector(W0);

	      X0.Nullify();
	      V0.Nullify();
	      W0.Nullify();
	    }

	  Add(1.0, A, C);

	  for (int j = 0; j < C.GetM(); j++)
	    CPPUNIT_ASSERT(C(j) == sum(j));

	  A.Deallocate();
	  C.Clear();
	  sum.Deallocate();

	}

      }
  }

  void copy()
  {
    srand(time(NULL));

    int length;

    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A, B;

    for (int N = 0; N < Nloop_; N++)
      {

	typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;
	vector_float_dense X0;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X0.Reallocate(length);
	    X0.FillRand();
	    A.AddVector(X0);
	    X0.Nullify();
	  }

	typedef Vector<float, VectSparse, NewAlloc<float> >
	  vector_float_sparse;
	vector_float_sparse X1, Y1;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X1.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X1.Index(l) = rand() % length;
		X1.Value(l) = rand();
	      }
	    X1.Assemble();
	    A.AddVector(X1);
	    X1.Nullify();
	  }

	typedef Vector<double, VectFull, NewAlloc<double> > vector_double_dense;
	vector_double_dense X2;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X2.Reallocate(length);
	    X2.FillRand();
	    A.AddVector(X2);
	    X2.Nullify();
	  }

	typedef Vector<double, VectSparse, NewAlloc<double> >
	  vector_double_sparse;
	vector_double_sparse X3;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X3.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X3.Index(l) = rand() % length;
		X3.Value(l) = rand();
	      }
	    X3.Assemble();
	    A.AddVector(X3);
	    X3.Nullify();
	  }

	Copy(A, B);

	for (int j = 0; j < A.GetFloatDense().GetNvector(); j++)
	  for (int l = 0; l < A.GetFloatDense().GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetFloatDense().GetVector(j)(l) ==
			   B.GetFloatDense().GetVector(j)(l));

	for (int j = 0; j < A.GetDoubleDense().GetNvector(); j++)
	  for (int l = 0; l < A.GetDoubleDense().GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetDoubleDense().GetVector(j)(l) ==
			   B.GetDoubleDense().GetVector(j)(l));

	for (int j = 0; j < A.GetFloatSparse().GetNvector(); j++)
	  for (int l = 0; l < A.GetFloatSparse().GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Index(l) ==
			       B.GetFloatSparse().GetVector(j).Index(l));
		CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Value(l) ==
			       B.GetFloatSparse().GetVector(j).Value(l));
	      }

	for (int j = 0; j < A.GetDoubleSparse().GetNvector(); j++)
	  for (int l = 0; l < A.GetDoubleSparse().GetVectorLength()(j); l++)
	    {
	      CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Index(l) ==
			     B.GetDoubleSparse().GetVector(j).Index(l));
	      CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Value(l) ==
			     B.GetDoubleSparse().GetVector(j).Value(l));
	    }

	for (int j = 0; j < A.GetNvector(); j++)
	  CPPUNIT_ASSERT(A.GetVectorLength()(j) ==
			 A.GetVectorLength()(j));


	A.Deallocate();
      }
  }


  void dot_product()
  {
    srand(time(NULL));

    int length;

    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A, B;

    for (int N = 0; N < Nloop_; N++)
      {
	double dot1(0.), dot2;

	typedef Vector<double, VectFull, NewAlloc<double> > vector_double_dense;
	vector_double_dense X2, V2;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;

	    X2.Reallocate(length);
	    X2.FillRand();
	    A.AddVector(X2);

	    V2.Reallocate(length);
	    V2.FillRand();
	    B.AddVector(V2);

	    dot1 += DotProd(X2, V2);

	    X2.Nullify();
	    V2.Nullify();
	  }

	typedef Vector<double, VectSparse, NewAlloc<double> > vector_double_sparse;
	vector_double_sparse X3, V3;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;

	    X3.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X3.Index(l) = rand() % length;
		X3.Value(l) = rand();
	      }
	    X3.Assemble();
	    A.AddVector(X3);

	    V3.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		V3.Index(l) = rand() % length;
		V3.Value(l) = rand();
	      }
	    V3.Assemble();
	    B.AddVector(V3);

	    dot1 += DotProd(X3, V3);

	    X3.Nullify();
	    V3.Nullify();
	  }

	dot2 = DotProd(A, B);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(dot1, dot2,
				     1.e-6 * dot2);

	A.Deallocate();
	B.Deallocate();
      }
  }


  void mlt_add()
  {
    srand(time(NULL));

    int length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // Add two dense vector collections.
	  typedef Vector<double, VectFull, NewAlloc<double> > vector_real_dense;
	  typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;

	  Matrix<real> M;
	  real alpha, beta;
	  alpha = rand();
	  beta = rand();

	  Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A;
	  vector_real_dense A_dense, U, W;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);
	      A_dense.PushBack(U);

	      U.Nullify();
	    }

	  vector_float_dense u;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      u.Reallocate(length);
	      u.FillRand();
	      A.AddVector(u);

	      W.Reallocate(length);
	      for (int i = 0; i < length; i++)
		W(i) = static_cast<double>(u(i));
	      A_dense.PushBack(W);

	      u.Nullify();
	    }


	  vector_real_dense Y1(m_), Y2(m_);
	  Y1.FillRand();
	  Copy(Y1, Y2);

	  M.Reallocate(m_, A.GetM());
	  M.FillRand();

	  MltAdd(alpha, M, A, beta, Y1);

	  MltAdd(alpha, M, A_dense, beta, Y2);

	  for (int l = 0; l < m_; l++)
	    CPPUNIT_ASSERT_DOUBLES_EQUAL(Y1(l), Y2(l),
					 1.e-6 * Y2(l));

	  A.Deallocate();
	}
      }
  }


  void write_read()
  {
     srand(time(NULL));

    int length;
    typedef double real;
    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A, B;
    for (int N = 0; N < Nloop_; N++)
      {
	typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;
	vector_float_dense X0;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X0.Reallocate(length);
	    X0.FillRand();
	    A.AddVector(X0);
	    X0.Nullify();
	  }

	typedef Vector<float, VectSparse, NewAlloc<float> >
	  vector_float_sparse;
	vector_float_sparse X1, Y1;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X1.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X1.Index(l) = rand() % length;
		X1.Value(l) = rand();
	      }
	    X1.Assemble();
	    A.AddVector(X1);
	    X1.Nullify();
	  }

	typedef Vector<double, VectFull, NewAlloc<double> > vector_double_dense;
	vector_double_dense X2;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X2.Reallocate(length);
	    X2.FillRand();
	    A.AddVector(X2);
	    X2.Nullify();
	  }

	typedef Vector<double, VectSparse, NewAlloc<double> >
	  vector_double_sparse;
	vector_double_sparse X3;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    X3.Reallocate(length);
	    for (int l = 0; l < length; l++)
	      {
		X3.Index(l) = rand() % length;
		X3.Value(l) = rand();
	      }
	    X3.Assemble();
	    A.AddVector(X3);
	    X3.Nullify();
	  }

	A.Write("test.bin");

	B.Read("test.bin");

	for (int j = 0; j < A.GetFloatDense().GetNvector(); j++)
	  for (int l = 0; l < A.GetFloatDense().GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetFloatDense().GetVector(j)(l) ==
			   B.GetFloatDense().GetVector(j)(l));

	for (int j = 0; j < A.GetDoubleDense().GetNvector(); j++)
	  for (int l = 0; l < A.GetDoubleDense().GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetDoubleDense().GetVector(j)(l) ==
			   B.GetDoubleDense().GetVector(j)(l));

	for (int j = 0; j < A.GetFloatSparse().GetNvector(); j++)
	  for (int l = 0; l < A.GetFloatSparse().GetVectorLength()(j); l++)
	    {
	      CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Index(l) ==
			       B.GetFloatSparse().GetVector(j).Index(l));
	      CPPUNIT_ASSERT(A.GetFloatSparse().GetVector(j).Value(l) ==
			     B.GetFloatSparse().GetVector(j).Value(l));
	    }

	for (int j = 0; j < A.GetDoubleSparse().GetNvector(); j++)
	  for (int l = 0; l < A.GetDoubleSparse().GetVectorLength()(j); l++)
	    {
	      CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Index(l) ==
			     B.GetDoubleSparse().GetVector(j).Index(l));
	      CPPUNIT_ASSERT(A.GetDoubleSparse().GetVector(j).Value(l) ==
			     B.GetDoubleSparse().GetVector(j).Value(l));
	    }

	for (int j = 0; j < A.GetNvector(); j++)
	  CPPUNIT_ASSERT(A.GetVectorLength()(j) ==
			 A.GetVectorLength()(j));

	A.Deallocate();
	B.Deallocate();
      }
  }


  void label()
  {
    srand(time(NULL));

    int length;

    Vector<FloatDouble, DenseSparseCollection, NewAlloc<float> > A;

    typedef Vector<float, VectFull, NewAlloc<float> > vector_float_dense;
    vector_float_dense X0, Y0;
    length = rand() % Nsub_vector_max_ + 1;
    X0.Reallocate(length);
    X0.FillRand();
    A.AddVector(X0, "X0");

    typedef Vector<float, VectSparse, NewAlloc<float> >
      vector_float_sparse;
    vector_float_sparse X1, Y1;
    length = rand() % Nsub_vector_max_ + 1;
    X1.Reallocate(length);
    for (int l = 0; l < length; l++)
      {
	X1.Index(l) = rand() % length;
	X1.Value(l) = rand();
      }
    X1.Assemble();
    A.AddVector(X1, "X1");

    typedef Vector<double, VectFull, NewAlloc<double> > vector_double_dense;
    vector_double_dense X2, Y2;
    length = rand() % Nsub_vector_max_ + 1;
    X2.Reallocate(length);
    X2.FillRand();
    A.AddVector(X2, "X2");

    typedef Vector<double, VectSparse, NewAlloc<double> >
      vector_double_sparse;
    vector_double_sparse X3, Y3;
    X3.Reallocate(length);
    for (int l = 0; l < length; l++)
      {
	X3.Index(l) = rand() % length;
	X3.Value(l) = rand();
      }
    X3.Assemble();
    A.AddVector(X3, "X3");

    A.GetVector("X0", Y0);
    A.GetVector("X1", Y1);
    A.GetVector("X2", Y2);
    A.GetVector("X3", Y3);

    for (int l = 0; l < X0.GetM(); l++)
      CPPUNIT_ASSERT(X0(l) == Y0(l));

    for (int l = 0; l < X1.GetM(); l++)
      {
	CPPUNIT_ASSERT(X1.Index(l) == Y1.Index(l));
	CPPUNIT_ASSERT(X1.Value(l) == Y1.Value(l));
      }

    for (int l = 0; l < X2.GetM(); l++)
      CPPUNIT_ASSERT(X2(l) == Y2(l));


    for (int l = 0; l < X3.GetM(); l++)
      {
	CPPUNIT_ASSERT(X3.Index(l) == Y3.Index(l));
	CPPUNIT_ASSERT(X3.Value(l) == Y3.Value(l));
      }


    A.Nullify();
    Y0.Nullify();
    Y1.Nullify();
    Y2.Nullify();
    Y3.Nullify();
  }

};
