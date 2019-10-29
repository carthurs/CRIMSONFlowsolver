// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
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



#include <cppunit/extensions/HelperMacros.h>

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU
#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"

using namespace Verdandi;


class UsefulFunctionTest: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(UsefulFunctionTest);
    CPPUNIT_TEST(test_get_inverse);
    CPPUNIT_TEST_SUITE_END();

protected:
    int Nloop_;
    int N_;

public:
    void setUp()
    {
        srand(time(NULL));
    }


    void tearDown()
    {
    }


    void test_get_inverse()
    {
        Nloop_ = 10;
        N_ = 15;
        get_inverse();
    }


    void get_inverse()
    {
        typedef double real;
        Matrix<real> A_dense, working_matrice;
        Matrix<real, General, RowSparse> A;
        Matrix<real, General, ArrayRowSparse> A_array;

        for (int k = 0; k < Nloop_; k++)
        {
            A_array.Reallocate(N_, N_);
            A_dense.Reallocate(N_, N_);
            A_dense.Fill(real(0));

            int r, s;
            double value;
            for (int l = 0; l < int(N_ * N_ / 2); l++)
            {
                r = rand() % N_;
                s = rand() % N_;
                value = real(rand());
                A_array.AddInteraction(r, s, value);
            }

            Copy(A_array, A);

            Copy(A, A_dense);

            Vector<real> wr, wi;
            working_matrice.Copy(A_dense);
            GetEigenvalues(working_matrice, wr, wi);
            bool is_inversible = true;
            for (int i = 0; i < wr.GetLength(); i++)
            {
                if (wr(i) == real(0) && wi(i) == real(0))
                {
                    is_inversible = false;
                    break;
                }
            }

            if (is_inversible)
            {
                GetInverse(A_dense);
                GetInverse(A);

                for (int i = 0; i < N_; i++)
                    for (int j = 0; j < N_; j++)
                        CPPUNIT_ASSERT(A(i, j) == A_dense(i, j));
            }
        }
    }


};
