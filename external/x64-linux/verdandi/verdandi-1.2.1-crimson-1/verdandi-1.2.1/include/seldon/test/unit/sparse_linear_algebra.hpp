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


#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <ctime>

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;


class SparseLinearAlgebraTest: public CppUnit::TestFixture
{


  CPPUNIT_TEST_SUITE(SparseLinearAlgebraTest);
  CPPUNIT_TEST(test_add);
  CPPUNIT_TEST(test_mlt);
  CPPUNIT_TEST(test_mlt_trans);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;
  int p_;
  int Nelement_;
  int Nloop_;

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_add()
  {
    Nloop_ = 1;

    m_ = 5;
    n_ = 10;
    Nelement_ = 0;
    add();

    m_ = 1;
    n_ = 1;
    Nelement_ = 1;
    add();

    Nloop_ = 5;

    m_ = 10;
    n_ = 1;
    Nelement_ = 20;
    add();

    m_ = 1;
    n_ = 10;
    Nelement_ = 20;
    add();

    Nloop_ = 100;

    m_ = 10;
    n_ = 1;
    Nelement_ = 5;
    add();

    m_ = 1;
    n_ = 10;
    Nelement_ = 5;
    add();

    Nloop_ = 200;

    m_ = 10;
    n_ = 25;
    Nelement_ = 1;
    add();

    m_ = 5;
    n_ = 5;
    Nelement_ = 2;
    add();

    m_ = 4;
    n_ = 4;
    Nelement_ = 2;
    add();

    m_ = 4;
    n_ = 4;
    Nelement_ = 10;
    add();

    m_ = 10;
    n_ = 25;
    Nelement_ = 3;
    add();

    m_ = 10;
    n_ = 25;
    Nelement_ = 20;
    add();

    m_ = 100;
    n_ = 250;
    Nelement_ = 20;
    add();

    m_ = 10;
    n_ = 10;
    Nelement_ = 200;
    add();
  }


  void add()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        Matrix<double> B_full(m_, n_);
        A_full.Zero();
        B_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        Matrix<double, General, ArrayRowSparse> B_array(m_, n_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            B_array.AddInteraction(i, j, value);
            B_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;
        Matrix<double, General, RowSparse> B;

        Copy(A_array, A);
        Copy(B_array, B);

        Add(1., A, B);
        Add(1., A_full, B_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < n_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(B_full(i, j), B(i, j),
                                         1.e-14 * B_full(i, j));
      }
  }


  void test_mlt()
  {
    Nloop_ = 1;

    m_ = 10;
    n_ = 25;
    p_ = 9;
    Nelement_ = 0;
    mlt();

    Nloop_ = 20;

    m_ = 10;
    n_ = 25;
    p_ = 9;
    Nelement_ = 20;
    mlt();

    m_ = 100;
    n_ = 250;
    p_ = 9;
    Nelement_ = 20;
    mlt();

    m_ = 10;
    n_ = 10;
    p_ = 2000;
    Nelement_ = 200;
    mlt();
  }


  void mlt()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        Matrix<double> B_full(n_, p_);
        Matrix<double> C_full(m_, p_);
        A_full.Zero();
        B_full.Zero();
        C_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        Matrix<double, General, ArrayRowSparse> B_array(n_, p_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % n_;
            j = rand() % p_;
            value = double(rand());
            B_array.AddInteraction(i, j, value);
            B_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;
        Matrix<double, General, RowSparse> B;
        Matrix<double, General, RowSparse> C;

        Copy(A_array, A);
        Copy(B_array, B);

        Mlt(A, B, C);
        Mlt(A_full, B_full, C_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < p_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(C_full(i, j), C(i, j),
                                         1.e-14 * C_full(i, j));

        MltAdd(-0.3, SeldonNoTrans, A, SeldonNoTrans, B, 0.9, C);
        MltAdd(-0.3, SeldonNoTrans, A_full,
               SeldonNoTrans, B_full, 0.9, C_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < p_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(C_full(i, j), C(i, j),
                                         1.e-14 * C_full(i, j));

        Matrix<double> C_full_2(m_, p_);

        Mlt(A_full, B, C_full_2);
        Mlt(A_full, B_full, C_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < p_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(C_full(i, j), C_full_2(i, j),
                                         1.e-14 * C_full(i, j));
        Transpose(B);
        Transpose(B_full);

        C_full.Fill(0);
        C_full_2.Fill(0);

        MltNoTransTrans(A_full, B, C_full_2);
        MltAdd(1., SeldonNoTrans, A_full, SeldonTrans, B_full, 0., C_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < p_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(C_full(i, j), C_full_2(i, j),
                                         1.e-14 * C_full(i, j));
      }
  }


  void test_mlt_trans()
  {
    Nloop_ = 1;

    m_ = 10;
    n_ = 25;
    p_ = 9;
    Nelement_ = 0;
    mlt_trans();

    Nloop_ = 20;

    m_ = 10;
    n_ = 25;
    p_ = 9;
    Nelement_ = 20;
    mlt_trans();

    m_ = 100;
    n_ = 250;
    p_ = 9;
    Nelement_ = 20;
    mlt_trans();

    m_ = 10;
    n_ = 10;
    p_ = 2000;
    Nelement_ = 200;
    mlt_trans();
  }


  void mlt_trans()
  {
    srand(time(NULL));

    int i, j;
    double value;

    for (int k = 0; k < Nloop_; k++)
      {
        Matrix<double> A_full(m_, n_);
        Matrix<double> B_full(p_, n_);
        Matrix<double> C_full(m_, p_);
        A_full.Zero();
        B_full.Zero();
        C_full.Zero();

        Matrix<double, General, ArrayRowSparse> A_array(m_, n_);
        Matrix<double, General, ArrayRowSparse> B_array(p_, n_);
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % m_;
            j = rand() % n_;
            value = double(rand());
            A_array.AddInteraction(i, j, value);
            A_full(i, j) += value;
          }
        for (int l = 0; l < Nelement_; l++)
          {
            i = rand() % p_;
            j = rand() % n_;
            value = double(rand());
            B_array.AddInteraction(i, j, value);
            B_full(i, j) += value;
          }

        Matrix<double, General, RowSparse> A;
        Matrix<double, General, RowSparse> B;
        Matrix<double, General, RowSparse> C;

        Copy(A_array, A);
        Copy(B_array, B);

        MltAdd(-0.1, SeldonNoTrans, A, SeldonTrans, B, 0., C);
        MltAdd(-0.1, SeldonNoTrans, A_full, SeldonTrans, B_full, 0., C_full);

        MltAdd(2.1, SeldonNoTrans, A, SeldonTrans, B, -1.5, C);
        MltAdd(2.1, SeldonNoTrans, A_full, SeldonTrans, B_full, -1.5, C_full);

        for (int i = 0; i < m_; i++)
          for (int j = 0; j < p_; j++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(C_full(i, j), C(i, j),
                                         1.e-14 * C_full(i, j));
      }
  }


};
