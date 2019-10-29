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


// These tests modify Seldon calls. They cannot be run along with regular
// tests.

#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_ABORT

#include "blas_call.hpp"

#include <cppunit/TestResult.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/HelperMacros.h>
using namespace CppUnit;

CPPUNIT_TEST_SUITE_REGISTRATION(BlasCallTest);

int main()
{
  TextUi::TestRunner runner;

  TestFactoryRegistry &registry = TestFactoryRegistry::getRegistry();

  runner.addTest(registry.makeTest());

  return runner.run("", false);
}
