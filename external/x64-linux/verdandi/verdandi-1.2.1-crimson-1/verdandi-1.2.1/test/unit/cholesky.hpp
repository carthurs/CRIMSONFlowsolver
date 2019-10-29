// Copyright (C) 2010 INRIA
// Author(s): Anne Tilloy
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
#include "method/BLUE.cxx"
using namespace Verdandi;


class CholeskyTest: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(CholeskyTest);
    CPPUNIT_TEST(test_dense);
    CPPUNIT_TEST_SUITE_END();

protected:
    int m_;
    int n_;
    int Nelement_;
    int Nloop_;

public:
    void setUp()
    {
        srand(time(NULL));
    }


    void tearDown()
    {
    }


    void test_dense()
    {
        dense_8x8();

        m_ = n_ = 3;
        Nloop_ = 10;
        dense();

        m_ = n_ = 7;
        Nloop_ = 10;
        dense();

        m_ = n_ = 30;
        Nloop_ = 10;
        dense();
    }


    void dense_8x8()
    {
        double data[8][8] = {{97.7354, 1.62757,	-0.452557, 0.163853, -111.94,
                                79.9274, -21.7045, 10.4545},
                            {1.62757, 97.2829, 1.85054, -0.67367, 79.9273,
                            -133.644, 90.4937, -43.6304},
                            {-0.452557, 1.85054, 96.8236, 2.60474, -21.7041,
                             90.4929, -156.021, 171.425},
                            {0.163853, -0.67367, 2.60474, 97.3002, 4.72342,
                            -20.4546, 83.4412, -131.921},
                            {-111.94, 79.9273, -21.7041, 4.72342, 305.574,
                            -231.258, 146.142, -99.9832},
                            {79.9274, -133.644, 90.4929, -20.4546, -231.258,
                             451.716, -320.1, 270.35},
                            {-21.7045, 90.4937, -156.021, 83.4412, 146.142,
                            -320.1, 530.871, -451.518},
                            {10.4545, -43.6304, 171.425, -131.921, -99.9832,
                            270.35,	-451.518, 623.277}};
        Matrix<double> A(8, 8), A_copy;
        A.SetData(8, 8, &data[0][0]);
        Copy(A, A_copy);

        GetCholesky(A);

        Matrix<double, General, RowMajor> B(8, 8);
        B.Zero();
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                for (int h = 0; h < 8; h++)
                    B(i, j) += A(i, h) * A(j, h);

        A.Nullify();

        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++)
                CPPUNIT_ASSERT_DOUBLES_EQUAL(A_copy(i, j),
                                             B(i, j),
                                             1.e-10
                                             + 1.e-10 * abs(A_copy(i, j)));
    }


    void dense()
    {
        for (int k = 0; k < Nloop_; k++)
        {

            Matrix<double, Symmetric, RowSymPacked> A(m_, n_), A_copy(m_, n_);
            double total = 0.;
            for (int i = 0; i < m_; i++)
                for (int j = i; j < n_; j++)
                {
                    A(i, j) = abs(double(rand()));
                    total += A(i, j);
                }
            for (int i = 0; i < m_; i++)
                A(i, i) += total;
            A_copy = A;

            GetCholesky(A);

            Matrix<double, Symmetric, RowSymPacked> B(m_, n_);
            B.Zero();
            for (int i = 0; i < m_; i++)
                for (int j = i; j < n_; j++)
                    for (int h = 0; h < i + 1; h++)
                        B(i, j) += A(i, h) * A(j, h);

            for (int i = 0; i < m_; i++)
                for (int j = i; j < n_; j++)
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(A_copy(i, j),
                                                 B(i, j),
                                                 1.e-10
                                                 + 1.e-10 * B(i, j));
        }

        for (int k = 0; k < Nloop_; k++)
        {

            Matrix<double, General, RowMajor> A(m_, n_), A_copy(m_, n_);
            double total = 0.;
            for (int i = 0; i < m_; i++)
                for (int j = i; j < n_; j++)
                {
                    A(i, j) = abs(double(rand()));
                    total += A(i, j);
                }
            for (int i = 0; i < m_; i++)
                A(i, i) += total;
            for (int i = 0; i < m_; i++)
                for (int j = 0; j < i; j++)
                    A(i, j) = A(j, i);
            A_copy = A;

            GetCholesky(A);

            Matrix<double, General, RowMajor> B(m_, n_);
            B.Zero();
            for (int i = 0; i < m_; i++)
                for (int j = 0; j < n_; j++)
                    for (int h = 0; h < n_; h++)
                        B(i, j) += A(i, h) * A(j, h);
            for (int i = 0; i < m_; i++)
                for (int j = 0; j < n_; j++)
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(A_copy(i, j),
                                                 B(i, j),
                                                 1.e-10
                                                 + 1.e-10 * B(i, j));
        }

    }

};
