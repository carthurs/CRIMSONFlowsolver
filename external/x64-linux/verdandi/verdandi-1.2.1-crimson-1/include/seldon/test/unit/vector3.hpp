// Copyright (C) 2010 Marc Fragu
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

#include "vector/Vector3.cxx"


class Vector3Test: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Vector3Test);
  CPPUNIT_TEST(test_reallocate);
  CPPUNIT_TEST(test_fill);
  CPPUNIT_TEST_SUITE_END();

protected:
  int length1_;
  int length2_;
  int length3_;

public:
  void setUp()
  {
  }


  void tearDown()
  {
  }


  void test_reallocate()
  {
    length1_ = 10;
    reallocate();
  }


  void test_fill()
  {
    length1_ = 10;
    length2_ = 25;
    length3_ = 87;
    fill();
  }


  void reallocate()
  {
    int i;
    Vector3<double> V;

    V.Reallocate(0);
    CPPUNIT_ASSERT(V.GetSize() == 0);
    V.Reallocate(length1_);
    CPPUNIT_ASSERT(V.GetSize() == length1_);

    V.Reallocate(0);
    CPPUNIT_ASSERT(V.GetLength() == 0);
    V.Reallocate(length1_);
    CPPUNIT_ASSERT(V.GetLength() == length1_);

    for (i = 0; i < length1_; i++)
      V.Reallocate(i, 0);
    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == 0);

    V.Reallocate(0);
    CPPUNIT_ASSERT(V.GetLength() == 0);
    V.Reallocate(length1_);
    CPPUNIT_ASSERT(V.GetLength() == length1_);

    for (i = 0; i < length1_; i++)
      V.Reallocate(i, i);
    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == i);

    for (i = 0; i < length1_; i++)
      V(i).Clear();
    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == 0);

    for (i = 0; i < length1_; i++)
      V(i).Reallocate(i);
    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == i);

    V.Reallocate(0);
    CPPUNIT_ASSERT(V.GetLength() == 0);
    V.Reallocate(length1_);
    CPPUNIT_ASSERT(V.GetLength() == length1_);

    for (i = 0; i < length1_; i++)
      V.Reallocate(i, length1_);

    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == length1_);

    for (i = 0; i < length1_; i++)
      V.Reallocate(i, i, i);

    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i, i) == i);

    V.Reallocate(0);
    CPPUNIT_ASSERT(V.GetLength() == 0);
  }


  void fill()
  {
    int i, j, k;
    Vector3<double> V(length1_);

    for (i = 0; i < length1_; i++)
      V.Reallocate(i, length2_);

    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
	V.Reallocate(i, j, length3_);

    V.Fill(-3.);
    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == length2_);

    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
	CPPUNIT_ASSERT(V.GetLength(i, j) == length3_);

    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
	for (k = 0; k < length3_; k++)
          {
            CPPUNIT_ASSERT(V(i, j, k) == -3.);
            CPPUNIT_ASSERT(V(i)(j)(k) == -3.);
            CPPUNIT_ASSERT(V.GetVector(i, j)(k) == -3.);
          }
  }
};
