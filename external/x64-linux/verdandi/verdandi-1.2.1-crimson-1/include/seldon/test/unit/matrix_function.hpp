// Copyright (C) 2010 INRIA
// Author(s): Vivien Mallet
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
using namespace Seldon;


class MatrixFunctionTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(MatrixFunctionTest);
  CPPUNIT_TEST(test_permute);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_permute()
  {
    m_ = 10;
    n_ = 10;
    permute();

    m_ = 10;
    n_ = 5;
    permute();
  }


  void permute()
  {
    Matrix<double, General, RowMajor> Ar(m_, n_), Ar_copy;
    Matrix<double, General, ColMajor> Ac(m_, n_);

    Ar.FillRand();
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        Ac(i, j) = Ar(i, j);
    Copy(Ar, Ar_copy);

    Vector<int> row_perm(m_), col_perm(n_);

    row_perm.Fill();
    col_perm.Fill();
    ApplyPermutation(Ar, row_perm, col_perm);
    ApplyPermutation(Ac, row_perm, col_perm);
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        {
          CPPUNIT_ASSERT(Ar(i, j) == Ac(i, j));
          CPPUNIT_ASSERT(Ar(i, j) == Ar_copy(i, j));
        }

    for (int k = 0; k < 2 * max(m_, n_); k++)
      {
        int is = rand() % m_;
        int it = rand() % m_;
        int tmp = row_perm(it);
        row_perm(it) = row_perm(is);
        row_perm(is) = tmp;

        int js = rand() % n_;
        int jt = rand() % n_;
        tmp = col_perm(jt);
        col_perm(jt) = col_perm(js);
        col_perm(js) = tmp;
      }

    ApplyPermutation(Ar, row_perm, col_perm);
    ApplyPermutation(Ac, row_perm, col_perm);
    ApplyInversePermutation(Ar, row_perm, col_perm);
    ApplyInversePermutation(Ac, row_perm, col_perm);
    for (int i = 0; i < m_; i++)
      for (int j = 0; j < n_; j++)
        {
          CPPUNIT_ASSERT(Ar(i, j) == Ac(i, j));
          CPPUNIT_ASSERT(Ar(i, j) == Ar_copy(i, j));
        }
  }
};
