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
using namespace Seldon;

typedef complex<float> complexfloat;
typedef complex<double> complexdouble;


class LapackTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(LapackTest);
  CPPUNIT_TEST(test_svd);
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_svd()
  {
    svd(1, 1, 2);
    svd(1, 5, 10);
    svd(5, 1, 10);
    svd(5, 5, 10);
    svd(5, 15, 10);
    svd(15, 5, 30);
  }


  void svd(int m, int n, int Nloop)
  {
    Matrix<@real_complex, General, @storage_blasGE> M@storage_blasGE_@real_complex(m, n), A@storage_blasGE_@real_complex, U@storage_blasGE_@real_complex, V@storage_blasGE_@real_complex, S@storage_blasGE_@real_complex(m, n), tmp@storage_blasGE_@real_complex(m, n);
    Vector<@real> sigma@storage_blasGE_@real;
    Vector<@real> sigma@storage_blasGE_complex@real;

    tmp@storage_blasGE_@real_complex.Zero();

    for (int l = 0; l < Nloop; l++)
      {
        M@storage_blasGE_@real_complex.FillRand();
        A@storage_blasGE_@real_complex = M@storage_blasGE_@real_complex;
        GetSVD(M@storage_blasGE_@real, sigma@storage_blasGE_@real, U@storage_blasGE_@real, V@storage_blasGE_@real);
        GetSVD(M@storage_blasGE_complex@real, sigma@storage_blasGE_complex@real, U@storage_blasGE_complex@real, V@storage_blasGE_complex@real);

        S@storage_blasGE_@real_complex.Zero();
        for (int i = 0; i < min(m, n); i++)
          {
            S@storage_blasGE_@real_complex(i, i) = sigma@storage_blasGE_@real_complex(i);
          }

        MltAdd(1., U@storage_blasGE_@real_complex, S@storage_blasGE_@real_complex, 0., tmp@storage_blasGE_@real_complex);
        S@storage_blasGE_@real_complex = tmp@storage_blasGE_@real_complex;
        MltAdd(1., S@storage_blasGE_@real_complex, V@storage_blasGE_@real_complex, 0., tmp@storage_blasGE_@real_complex);

        for (int i = 0; i < m; i++)
          for (int j = 0; j < n; j++)
            {
              CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp@storage_blasGE_@real(i, j), A@storage_blasGE_@real(i, j), 1.e-2 * A@storage_blasGE_@real(i, j));
              CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp@storage_blasGE_complex@real(i, j).real(), A@storage_blasGE_complex@real(i, j).real(), 1.e-2 * A@storage_blasGE_complex@real(i, j).real());
              CPPUNIT_ASSERT_DOUBLES_EQUAL(tmp@storage_blasGE_complex@real(i, j).imag(), A@storage_blasGE_complex@real(i, j).imag(), 1.e-2 * A@storage_blasGE_complex@real(i, j).imag());
            }
      }
  }
};
