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

// If SELDON_UNIT_TEST_INITIALIZE is defined, then the vector and matrices are
// filled with zeros. Otherwise, they are left with uninitialized values,
// which is recommended as Valgrind can detect calls to generic
// functions. Several generic functions call Blas: this activates the Blas
// flag, although the wrong function is actually called.
// #define SELDON_UNIT_TEST_INITIALIZE

// In order to load the call-testing Blas functions:
#define CBLAS_H
int blas_called = 0;
#include "cblas_call.h"
#define IS_CALL_FLAG_ACTIVATED CPPUNIT_ASSERT(blas_called == 1); blas_called = 0

#define SELDON_WITH_BLAS
#include "Seldon.hxx"
#include "SeldonSolver.hxx"
using namespace Seldon;

typedef complex<float> complexfloat;
typedef complex<double> complexdouble;


class BlasCallTest: public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(BlasCallTest);
  CPPUNIT_TEST(test_blas_call_1);
  CPPUNIT_TEST(test_blas_call_2);
  CPPUNIT_TEST(test_blas_call_3);
  CPPUNIT_TEST_SUITE_END();

protected:
  int m_;
  int n_;

  @real a@real, b@real, c@real, d@real;
  @complex a@complex, b@complex, c@complex, d@complex;
  Vector<@real> V@real;
  Vector<@complex> V@complex;

  Matrix<@real_complex, General, @storage_blasGE> M@storage_blasGE_@real_complex;
  Matrix<@complex, General, @storage_blasHE> M@storage_blasHE_@complex;
  Matrix<@complex, General, @storage_blasHP> M@storage_blasHP_@complex;
  Matrix<@real_complex, General, @storage_blasSY> M@storage_blasSY_@real_complex;
  Matrix<@real_complex, General, @storage_blasSP> M@storage_blasSP_@real_complex;
  Matrix<@real_complex, General, @storage_blasTR> M@storage_blasTR_@real_complex;
  Matrix<@real_complex, General, @storage_blasTP> M@storage_blasTP_@real_complex;

public:
  void setUp()
  {
    m_ = 25;

    a@real = @real(1);
    b@real = @real(1);
    c@real = @real(0);
    d@real = @real(0);
    a@complex = @complex(1, 0);
    b@complex = @complex(1, 0);
    c@complex = @complex(0, 0);
    d@complex = @complex(0, 0);

    V@real_complex.Reallocate(m_);

    M@storage_blasGE_@real_complex.Reallocate(m_, m_);
    M@storage_blasHE_@complex.Reallocate(m_, m_);
    M@storage_blasHP_@complex.Reallocate(m_, m_);
    M@storage_blasSY_@real_complex.Reallocate(m_, m_);
    M@storage_blasSP_@real_complex.Reallocate(m_, m_);
    M@storage_blasTR_@real_complex.Reallocate(m_, m_);
    M@storage_blasTP_@real_complex.Reallocate(m_, m_);

#ifdef SELDON_UNIT_TEST_INITIALIZE
    V@real_complex.Zero();

    M@storage_blasGE_@real_complex.Zero();
    M@storage_blasHE_@complex.Zero();
    M@storage_blasHP_@complex.Zero();
    M@storage_blasSY_@real_complex.Zero();
    M@storage_blasSP_@real_complex.Zero();
    M@storage_blasTR_@real_complex.Zero();
    M@storage_blasTP_@real_complex.Zero();
#endif
  }


  void tearDown()
  {
  }


  void test_blas_call_1()
  {
    GenRot(a@real, b@real, c@real, d@real); IS_CALL_FLAG_ACTIVATED;
    ApplyRot(V@real, V@real, c@real, d@real); IS_CALL_FLAG_ACTIVATED;
    ApplyModifRot(V@real, V@real, &c@real); IS_CALL_FLAG_ACTIVATED;
    Swap(V@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    Mlt(a@real, Vcomplex@real); IS_CALL_FLAG_ACTIVATED;
    Mlt(a@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    Copy(V@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    Add(a@real_complex, V@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    DotProd(V@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    {
      float a = 0;
      Vector<float> V;
      ScaledDotProd(a, V, V); IS_CALL_FLAG_ACTIVATED;
    }
    DotProdConj(V@complex, V@complex); IS_CALL_FLAG_ACTIVATED;
    Norm1(V@real_complex); IS_CALL_FLAG_ACTIVATED;
    Norm2(V@real_complex); IS_CALL_FLAG_ACTIVATED;
    GetMaxAbsIndex(V@real_complex); IS_CALL_FLAG_ACTIVATED;
  }


  void test_blas_call_2()
  {
    Mlt(M@storage_blasTR_@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    Mlt(M@storage_blasTP_@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;

    MltAdd(a@real_complex, M@storage_blasGE_@real_complex, V@real_complex, b@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(a@real_complex, @trans, M@storage_blasGE_@real_complex, V@real_complex, b@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(a@complex, M@storage_blasHE_@complex, V@complex, b@complex, V@complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(a@complex, M@storage_blasHP_@complex, V@complex, b@complex, V@complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(a@real, M@storage_blasSY_@real, V@real, b@real, V@real); IS_CALL_FLAG_ACTIVATED;
    MltAdd(a@real, M@storage_blasSP_@real, V@real, b@real, V@real); IS_CALL_FLAG_ACTIVATED;

    Rank1Update(a@real_complex, V@real_complex, V@real_complex, M@storage_blasGE_@real_complex); IS_CALL_FLAG_ACTIVATED;
    Rank1Update(a@complex, V@complex, @conj, V@complex, M@storage_blasGE_@complex); IS_CALL_FLAG_ACTIVATED;
    Rank1Update(a@real, V@real, M@storage_blasSP_@real); IS_CALL_FLAG_ACTIVATED;
    Rank1Update(a@real, V@real, @uplo, M@storage_blasSP_@real); IS_CALL_FLAG_ACTIVATED;
    Rank1Update(a@real, Vcomplex@real, M@storage_blasHP_complex@real); IS_CALL_FLAG_ACTIVATED;
    Rank1Update(a@real, Vcomplex@real, @uplo, M@storage_blasHP_complex@real); IS_CALL_FLAG_ACTIVATED;

    Rank2Update(a@real, V@real, V@real, M@storage_blasSP_@real); IS_CALL_FLAG_ACTIVATED;
    Rank2Update(a@real, V@real, V@real, @uplo, M@storage_blasSP_@real); IS_CALL_FLAG_ACTIVATED;
    Rank2Update(a@real, Vcomplex@real, Vcomplex@real, M@storage_blasHP_complex@real); IS_CALL_FLAG_ACTIVATED;
    Rank2Update(a@real, Vcomplex@real, Vcomplex@real, @uplo, M@storage_blasHP_complex@real); IS_CALL_FLAG_ACTIVATED;

    Solve(M@storage_blasT_@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
    Solve(@trans, @diag, M@storage_blasT_@real_complex, V@real_complex); IS_CALL_FLAG_ACTIVATED;
  }


  void test_blas_call_3()
  {
    MltAdd(a@real_complex, M@storage_blasGE_@real_complex, M@storage_blasGE_@real_complex, b@real_complex, M@storage_blasGE_@real_complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(a@real_complex, @trans, M@storage_blasGE_@real_complex, @trans+, M@storage_blasGE_@real_complex, b@real_complex, M@storage_blasGE_@real_complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(@side, a@real_complex, M@colrowSym_@real_complex, M@colrowMajor_@real_complex, b@real_complex, M@colrowMajor_@real_complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(@side, a@real_complex, @uplo, M@colrowSym_@real_complex, M@colrowMajor_@real_complex, b@real_complex, M@colrowMajor_@real_complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(@side, a@complex, M@colrowHerm_@complex, M@colrowMajor_@complex, b@complex, M@colrowMajor_@complex); IS_CALL_FLAG_ACTIVATED;
    MltAdd(@side, a@complex, @uplo, M@colrowHerm_@complex, M@colrowMajor_@complex, b@complex, M@colrowMajor_@complex); IS_CALL_FLAG_ACTIVATED;

    Mlt(@side, a@real_complex, M@colrow@ulTriang_@real_complex, M@colrowMajor_@real_complex); IS_CALL_FLAG_ACTIVATED;
    Mlt(@side, a@real_complex, @trans, @diag, M@colrow@ulTriang_@real_complex, M@colrowMajor_@real_complex); IS_CALL_FLAG_ACTIVATED;

    Solve(@side, a@real_complex, M@colrow@ulTriang_@real_complex, M@colrowMajor_@real_complex); IS_CALL_FLAG_ACTIVATED;
    Solve(@side, a@real_complex, @trans, @diag, M@colrow@ulTriang_@real_complex, M@colrowMajor_@real_complex); IS_CALL_FLAG_ACTIVATED;
  }
};
