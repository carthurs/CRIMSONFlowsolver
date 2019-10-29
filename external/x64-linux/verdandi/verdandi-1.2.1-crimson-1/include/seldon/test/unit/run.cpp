// Copyright (C) 2001-2010 Vivien Mallet, Marc Fragu
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

#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_ABORT

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#include "array3d.hpp"
#include "matrix.hpp"
#include "matrix_function.hpp"
#include "submatrix.hpp"
#include "lapack.hpp"
#include "sparse_linear_algebra.hpp"
#include "sparse_matrix.hpp"
#include "vector2.hpp"
#include "vector3.hpp"
#include "vector_collection.hpp"
#include "heterogeneous_collection.hpp"
#include "matrix_collection.hpp"
#include "heterogeneous_matrix_collection.hpp"
#include <cppunit/TestResult.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
using namespace CppUnit;


CPPUNIT_TEST_SUITE_REGISTRATION(MatrixTest);
CPPUNIT_TEST_SUITE_REGISTRATION(MatrixFunctionTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SubMatrixTest);
CPPUNIT_TEST_SUITE_REGISTRATION(Array3DTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SparseLinearAlgebraTest);
CPPUNIT_TEST_SUITE_REGISTRATION(SparseMatrixTest);
CPPUNIT_TEST_SUITE_REGISTRATION(LapackTest);
CPPUNIT_TEST_SUITE_REGISTRATION(Vector2Test);
CPPUNIT_TEST_SUITE_REGISTRATION(Vector3Test);
CPPUNIT_TEST_SUITE_REGISTRATION(VectorCollectionTest);
CPPUNIT_TEST_SUITE_REGISTRATION(HeterogeneousCollectionTest);
CPPUNIT_TEST_SUITE_REGISTRATION(MatrixCollectionTest);
CPPUNIT_TEST_SUITE_REGISTRATION(HeterogeneousMatrixCollectionTest);

int main()
{
  TRY;

  TextUi::TestRunner runner;

  TestFactoryRegistry &registry = TestFactoryRegistry::getRegistry();

  runner.addTest(registry.makeTest());

  return runner.run("", false);

  END;
}
