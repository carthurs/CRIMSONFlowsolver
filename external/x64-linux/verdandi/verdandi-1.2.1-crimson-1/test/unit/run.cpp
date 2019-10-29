// Copyright (C) 2010 INRIA
// Author(s): Vivien Mallet
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#define VERDANDI_DEBUG_LEVEL_4
#define VERDANDI_WITH_ABORT

#include "blue.hpp"
#include "chi_2.hpp"
#include "cholesky.hpp"
#include "useful_function.hpp"

#include <cppunit/TestResult.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
using namespace CppUnit;

CPPUNIT_TEST_SUITE_REGISTRATION(BLUETest);
CPPUNIT_TEST_SUITE_REGISTRATION(Chi2Test);
CPPUNIT_TEST_SUITE_REGISTRATION(CholeskyTest);
CPPUNIT_TEST_SUITE_REGISTRATION(UsefulFunctionTest);

int main()
{
  VERDANDI_TRY;

  TextUi::TestRunner runner;

  TestFactoryRegistry &registry = TestFactoryRegistry::getRegistry();

  runner.addTest(registry.makeTest());

  return runner.run("", false);

  VERDANDI_END;
}
