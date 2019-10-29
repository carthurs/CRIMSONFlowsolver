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


class Array3DTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Array3DTest);
  CPPUNIT_TEST(test_fill);
  CPPUNIT_TEST(test_reallocate);
  CPPUNIT_TEST_SUITE_END();

protected:
  int length1_;
  int length2_;
  int length3_;

public:
  void setUp()
  {
    length1_ = 10;
    length2_ = 25;
    length3_ = 9;
  }


  void tearDown()
  {
  }


  void test_fill()
  {
    int i, j, k;
    Array3D<double> A(length1_, length2_, length3_);
    A.Fill();
    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
	for (k = 0; k < length3_; k++)
	  CPPUNIT_ASSERT(A(i, j, k) == double(i * length2_ * length3_
					      + j * length3_ + k));
  }


  void test_reallocate()
  {
    Array3D<double> A(length1_, length2_, length3_);
    A.Reallocate(length2_, length1_, length3_);
    CPPUNIT_ASSERT(A.GetLength1() == length2_);
    CPPUNIT_ASSERT(A.GetLength2() == length1_);
    CPPUNIT_ASSERT(A.GetLength3() == length3_);
    A.Reallocate(0, length2_, 0);
    CPPUNIT_ASSERT(A.GetLength1() == 0);
    CPPUNIT_ASSERT(A.GetLength2() == length2_);
    CPPUNIT_ASSERT(A.GetLength3() == 0);
    A.Reallocate(0, 0, 0);
    CPPUNIT_ASSERT(A.GetLength1() == 0);
    CPPUNIT_ASSERT(A.GetLength2() == 0);
    CPPUNIT_ASSERT(A.GetLength3() == 0);
  }
};
