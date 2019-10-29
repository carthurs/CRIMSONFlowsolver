// Copyright (C) 2001-2009 Marc Fragu
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
#include "vector/VectorCollection.cxx"
using namespace Seldon;


class VectorCollectionTest: public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE(VectorCollectionTest);
  CPPUNIT_TEST(test_mlt);
  CPPUNIT_TEST(test_add);
  CPPUNIT_TEST(test_copy);
  CPPUNIT_TEST(test_dot_product);
  CPPUNIT_TEST(test_mlt_add);
  CPPUNIT_TEST(test_write_read);
  CPPUNIT_TEST(test_label);
  CPPUNIT_TEST(test_collection);
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
    Nvector_ = 5;
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


  void test_collection()
  {
    Nvector_ = 2;
    Nsub_vector_max_ = 4;
    collection();
  }


  void mlt()
  {
    srand(time(NULL));

    int length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense vector collection.
	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> mlt;
	  Vector<vector_real_dense, Collection> A;
	  vector_real_dense U, W;
	  real alpha;
	  alpha = rand();
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();

	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);

	      Mlt(alpha, W);

	      mlt.AddVector(W);

	      U.Nullify();
	      W.Nullify();
	    }

	  Mlt(alpha, A);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      CPPUNIT_ASSERT(A.GetVector(j)(l) == mlt.GetVector(j)(l));

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    mlt.GetVector()(k).~Vector();
	}

	{
	  // For sparse vector collection.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> mlt;
	  Vector<vector_real_sparse, Collection> A;
	  vector_real_sparse U, W;
	  real alpha;
	  alpha = rand();
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      W.Reallocate(length);
	      W.Copy(U);

	      A.AddVector(U);

	      Mlt(alpha, W);

	      mlt.AddVector(W);

	      U.Nullify();
	      W.Nullify();
	    }

	  Mlt(alpha, A);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(A.GetVector(j).Index(l) ==
			       mlt.GetVector(j).Index(l));
		CPPUNIT_ASSERT(A.GetVector(j).Value(l) ==
			       mlt.GetVector(j).Value(l));
	      }

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    mlt.GetVector(k).~Vector();
	}

      }
  }


  void add()
  {
    srand(time(NULL));

    int length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // Add two dense vector collections.
	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> sum;
	  Vector<vector_real_dense, Collection> A, B;
	  vector_real_dense U, V, W, X;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.AddVector(V);

	      X.Reallocate(length);
	      X.Copy(V);

	      Add(1.0, W, X);

	      sum.AddVector(X);

	      U.Nullify();
	      V.Nullify();
	      X.Nullify();
	    }

	  Add(1.0, A, B);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      CPPUNIT_ASSERT(B.GetVector(j)(l) == sum.GetVector(j)(l));

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    B.GetVector()(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    sum.GetVector()(k).~Vector();
	}

	{
	  // Add two sparse vector collections.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> sum;
	  Vector<vector_real_sparse, Collection> A, B;
	  vector_real_sparse U, V, W, X;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      W.Reallocate(length);
	      W.Copy(U);

	      V.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  V.Index(l) = rand() % length;
		  V.Value(l) = rand();
		}
	      V.Assemble();

	      X.Reallocate(length);
	      X.Copy(V);

	      A.AddVector(U);
	      B.AddVector(V);

	      Add(1.0, W, X);

	      sum.AddVector(X);

	      U.Nullify();
	      V.Nullify();
	      X.Nullify();
	    }

	  Add(1.0, A, B);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(B.GetVector(j).Index(l) ==
			       sum.GetVector(j).Index(l));
		CPPUNIT_ASSERT(B.GetVector(j).Value(l) ==
			       sum.GetVector(j).Value(l));
	      }

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    B.GetVector(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
            sum.GetVector(k).~Vector();
	}

	{  // Add a dense vector to a dense vector collections.
	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> sum;
	  Vector<vector_real_dense, Collection> A;
	  vector_real_dense B, U, V, W, X;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      W.Reallocate(length);
	      W.Copy(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.PushBack(V);

	      X.Reallocate(length);
	      X.Copy(V);

	      Add(1.0, W, X);

	      sum.AddVector(X);

	      U.Nullify();
	      X.Nullify();
	    }

	  Add(1.0, B, A);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      CPPUNIT_ASSERT(A.GetVector(j)(l) == sum.GetVector(j)(l));

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    sum.GetVector()(k).~Vector();
	}

      }
  }


  void copy()
  {
    srand(time(NULL));

    int length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense vector collection.

	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> A, B;
	  vector_real_dense U, V;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.AddVector(V);

	      U.Nullify();
	      V.Nullify();
	    }

	  B.Deallocate();

	  Copy(A, B);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      CPPUNIT_ASSERT(A.GetVector(j)(l) == B.GetVector(j)(l));

	  A.Deallocate();

	}

	{
	  // For sparse vector collection.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> A, B;
	  vector_real_sparse U, V;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      V.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  V.Index(l) = rand() % length;
		  V.Value(l) = rand();
		}
	      V.Assemble();

	      A.AddVector(U);
	      B.AddVector(V);

	      U.Nullify();
	      V.Nullify();
	    }

	  B.Deallocate();

	  Copy(A, B);

	  for (int j = 0; j < Nvector_; j++)
	    for (int l = 0; l < A.GetVectorLength()(j); l++)
	      {
		CPPUNIT_ASSERT(B.GetVector(j).Index(l) ==
			       A.GetVector(j).Index(l));
		CPPUNIT_ASSERT(B.GetVector(j).Value(l) ==
			       A.GetVector(j).Value(l));
	      }

	  A.Deallocate();
	}
      }
  }


  void dot_product()
  {
    srand(time(NULL));

    int length;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	{  // For dense vector collection.

	  typedef Vector<real> vector_real_dense;

	  Vector<vector_real_dense, Collection> A, B;
	  vector_real_dense U, V;
	  real dot1(0.), dot2;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      U.FillRand();
	      A.AddVector(U);

	      V.Reallocate(length);
	      V.FillRand();
	      B.AddVector(V);

	      dot1 += DotProd(U, V);

	      U.Nullify();
	      V.Nullify();
	    }

	  dot2 = DotProd(A, B);

	  CPPUNIT_ASSERT(dot1 == dot2);

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    B.GetVector()(k).~Vector();
	}

	{
	  // For sparse vector collection.

	  typedef Vector<real, VectSparse, NewAlloc<real> >
	    vector_real_sparse;

	  Vector<vector_real_sparse, Collection> A, B;
	  real dot1 = real(0.), dot2;
	  vector_real_sparse U, V;
	  for (int k = 0; k < Nvector_; k++)
	    {
	      length = rand() % Nsub_vector_max_ + 1;

	      U.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  U.Index(l) = rand() % length;
		  U.Value(l) = rand();
		}
	      U.Assemble();

	      V.Reallocate(length);
	      for (int l = 0; l < length; l++)
		{
		  V.Index(l) = rand() % length;
		  V.Value(l) = rand();
		}
	      V.Assemble();

	      A.AddVector(U);
	      B.AddVector(V);

	      dot1 += DotProd(U, V);

	      U.Nullify();
	      V.Nullify();
	    }

	  dot2 = DotProd(A, B);

	  CPPUNIT_ASSERT(dot1 == dot2);

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector(k).~Vector();
	  for (int k = 0; k < Nvector_; k++)
	    B.GetVector(k).~Vector();
	}
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
	  typedef Vector<int> vector_real_dense;

	  Matrix<real> M;
	  real alpha, beta;
	  alpha = rand();
	  beta = rand();

	  Vector<vector_real_dense, Collection> A;
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
	  vector_real_dense Y1(m_), Y2(m_);
	  Y1.FillRand();
	  Copy(Y1, Y2);

	  M.Reallocate(m_, A.GetM());
	  M.FillRand();

	  MltAdd(alpha, M, A, beta, Y1);

	  MltAdd(alpha, M, A_dense, beta, Y2);

	  for (int l = 0; l < m_; l++)
	    CPPUNIT_ASSERT(Y1(l) == Y2(l));

	  for (int k = 0; k < Nvector_; k++)
	    A.GetVector()(k).~Vector();
	}
      }
  }


  void write_read()
  {
    srand(time(NULL));

    int length;
    Vector<int, VectFull, MallocAlloc<int> > length_vector;
    typedef double real;

    for (int N = 0; N < Nloop_; N++)
      {
	typedef Vector<real> vector_real_dense;

	Vector<vector_real_dense, Collection> A, B;

	length_vector.Clear();
	vector_real_dense U;
	for (int k = 0; k < Nvector_; k++)
	  {
	    length = rand() % Nsub_vector_max_ + 1;
	    length_vector.PushBack(length);
	    U.Reallocate(length);
	    U.FillRand();
	    A.AddVector(U);
	    U.Nullify();
	  }

	A.Write("test.bin");

	B.Read("test.bin", length_vector);

	for (int j = 0; j < Nvector_; j++)
	  for (int l = 0; l < A.GetVectorLength()(j); l++)
	    CPPUNIT_ASSERT(A.GetVector(j)(l) == B.GetVector(j)(l));

	A.Deallocate();
	B.Deallocate();
      }
  }


  void label()
  {
    typedef Vector<double> vector_real_dense;
    int length;

    srand(time(NULL));

    Vector<vector_real_dense, Collection> A;
    vector_real_dense vector1, vector2;

    length = rand() % Nsub_vector_max_ + 1;
    vector1.Reallocate(length);
    vector1.FillRand();
    A.AddVector(vector1, "vector1");

    length = rand() % Nsub_vector_max_ + 1;
    vector2.Reallocate(length);
    vector2.FillRand();
    A.AddVector(vector2, "vector2");

    vector_real_dense U, V;

    U.Copy(A.GetVector("vector1"));
    V.Copy(A.GetVector("vector2"));

    for (int l = 0; l < U.GetM(); l++)
      CPPUNIT_ASSERT(vector1(l) == U(l));
    for (int l = 0; l < V.GetM(); l++)
      CPPUNIT_ASSERT(vector2(l) == V(l));

    A.Nullify();
  }


  void collection()
  {
    typedef double real;
    typedef Vector<real> vector_real_dense;
    typedef Vector<vector_real_dense, Collection>
      collection_real_dense;

    int length;
    Vector<collection_real_dense, Collection> C;
    collection_real_dense A, B, D;
    vector_real_dense U;

    for (int k = 0; k < Nvector_; k++)
      {
	length = rand() % Nsub_vector_max_ + 1;
	U.Reallocate(length);
	U.FillRand();
	A.AddVector(U);
	U.Nullify();
      }

    for (int k = 0; k < Nvector_; k++)
      {
	length = rand() % Nsub_vector_max_ + 1;
	U.Reallocate(length);
	U.FillRand();
	B.AddVector(U);
	U.Nullify();
      }

    D.Copy(A);

    C.AddVector(A);
    C.AddVector(B);
    C.AddVector(D);

    // Checking adress.
    for (int i = 0; i < A.GetNvector(); i++)
      CPPUNIT_ASSERT(&A.GetVector(i).GetData()[0] ==
		     &C.GetVector(0).GetVector(i).GetData()[0]);

    // Checking adress.
    for (int i = 0; i < B.GetNvector(); i++)
      CPPUNIT_ASSERT(&B.GetVector(i).GetData()[0] ==
		     &C.GetVector(1).GetVector(i).GetData()[0]);

    A.Deallocate();
    B.Deallocate();
  }

};
