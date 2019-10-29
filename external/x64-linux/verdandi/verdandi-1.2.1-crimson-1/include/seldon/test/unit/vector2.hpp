// Copyright (C) 2010 Vivien Mallet
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

#include "vector/Vector2.cxx"


class Vector2Test: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(Vector2Test);
  CPPUNIT_TEST(test_reallocate);
  CPPUNIT_TEST(test_fill);
  CPPUNIT_TEST_SUITE_END();

protected:
  int length1_;
  int length2_;

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
    fill();
  }


  void reallocate()
  {
    int i;
    Vector2<double> V;

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
  }


  void fill()
  {
    int i, j;
    Vector2<double> V(length1_);

    for (i = 0; i < length1_; i++)
      V.Reallocate(i, length2_);
    V.Fill(-3.);
    for (i = 0; i < length1_; i++)
      CPPUNIT_ASSERT(V.GetLength(i) == length2_);
    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
        {
          CPPUNIT_ASSERT(V(i, j) == -3.);
          CPPUNIT_ASSERT(V(i)(j) == -3.);
          CPPUNIT_ASSERT(V.GetVector(i)(j) == -3.);
        }

    for (i = 0; i < length1_; i++)
      V(i).Fill();
    for (i = 0; i < length1_; i++)
      for (j = 0; j < length2_; j++)
        {
          CPPUNIT_ASSERT(V(i, j) == double(j));
          CPPUNIT_ASSERT(V(i)(j) == double(j));
          CPPUNIT_ASSERT(V.GetVector(i)(j) == double(j));
        }
  }
};
