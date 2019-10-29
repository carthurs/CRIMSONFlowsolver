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

#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;

typedef complex<float> complexfloat;
typedef complex<double> complexdouble;


class SubMatrixTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(SubMatrixTest);
  CPPUNIT_TEST(test_constructor);
  CPPUNIT_TEST(test_access);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;
  int m_sub_;
  int n_sub_;
  Vector<int> row_;
  Vector<int> column_;

public:
  void setUp()
  {
    m_ = 25;
    n_ = 10;
    m_sub_ = 7;
    n_sub_ = 3;
  }


  void tearDown()
  {
  }


  void test_constructor()
  {
    row_.Reallocate(m_sub_);
    column_.Reallocate(n_sub_);

    row_.Fill();
    column_.Fill();

    Matrix<@real_complex, General, @storage_full_real_complex> M@storage_full_real_complex_@real_complex(m_, n_);
    {
      SubMatrix<Matrix<@real_complex, General, @storage_full_real_complex> > SubM@storage_full_real_complex_@real_complex(M@storage_full_real_complex_@real_complex, row_, column_);
      CPPUNIT_ASSERT(SubM@storage_full_real_complex_@real_complex.GetM() == m_sub_);
      CPPUNIT_ASSERT(SubM@storage_full_real_complex_@real_complex.GetN() == n_sub_);
    }

    Matrix<@complex, General, @storage_full_complex> M@storage_full_complex_@complex(m_, n_);
    {
      SubMatrix<Matrix<@complex, General, @storage_full_complex> > SubM@storage_full_complex_@complex(M@storage_full_complex_@complex, row_, column_);
      CPPUNIT_ASSERT(SubM@storage_full_complex_@complex.GetM() == m_sub_);
      CPPUNIT_ASSERT(SubM@storage_full_complex_@complex.GetN() == n_sub_);
    }
  }


  void test_access()
  {

    /*** Checks with the whole matrix extracted  ***/

    row_.Reallocate(m_);
    column_.Reallocate(n_);

    row_.Fill();
    column_.Fill();

    {
      Matrix<@real_complex, General, @storage_rectangular_full> M@storage_rectangular_full_@real_complex(m_, n_);
      M@storage_rectangular_full_@real_complex.Fill();
      SubMatrix<Matrix<@real_complex, General, @storage_rectangular_full> > SubM@storage_rectangular_full_@real_complex(M@storage_rectangular_full_@real_complex, row_, column_);
      for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
          {
            CPPUNIT_ASSERT(SubM@storage_rectangular_full_@real_complex(i, j) == M@storage_rectangular_full_@real_complex(i, j));
          }
    }

    {
      Matrix<@real, General, @storage_rectangular_sparse> M@storage_rectangular_sparse_@real(m_, n_);
      SubMatrix<Matrix<@real, General, @storage_rectangular_sparse> > SubM@storage_rectangular_sparse_@real(M@storage_rectangular_sparse_@real, row_, column_);
      for (int i = 0; i < m_; i++)
        for (int j = 0; j < n_; j++)
          {
            CPPUNIT_ASSERT(SubM@storage_rectangular_sparse_@real(i, j) == M@storage_rectangular_sparse_@real(i, j));
          }
    }

    /*** Checks with a submatrix  ***/

    row_.Reallocate(m_sub_);
    column_.Reallocate(n_sub_);

    for (int i = 0; i < m_sub_; i++)
      row_(i) = (i * m_) / m_sub_;
    for (int j = 0; j < n_sub_; j++)
      column_(j) = (j * n_) / n_sub_;

    {
      Matrix<@real_complex, General, @storage_rectangular_full> M@storage_rectangular_full_@real_complex(m_, n_);
      M@storage_rectangular_full_@real_complex.Fill();
      SubMatrix<Matrix<@real_complex, General, @storage_rectangular_full> > SubM@storage_rectangular_full_@real_complex(M@storage_rectangular_full_@real_complex, row_, column_);
      for (int i = 0; i < m_sub_; i++)
        for (int j = 0; j < n_sub_; j++)
          {
            CPPUNIT_ASSERT(SubM@storage_rectangular_full_@real_complex(i, j) == M@storage_rectangular_full_@real_complex(row_(i), column_(j)));
          }
    }

    {
      Matrix<@real, General, @storage_rectangular_sparse> M@storage_rectangular_sparse_@real(m_, n_);
      SubMatrix<Matrix<@real, General, @storage_rectangular_sparse> > SubM@storage_rectangular_sparse_@real(M@storage_rectangular_sparse_@real, row_, column_);
      for (int i = 0; i < m_sub_; i++)
        for (int j = 0; j < n_sub_; j++)
          {
            CPPUNIT_ASSERT(SubM@storage_rectangular_sparse_@real(i, j) == M@storage_rectangular_sparse_@real(row_(i), column_(j)));
          }
    }

  }
};
