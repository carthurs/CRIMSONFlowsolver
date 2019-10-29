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



#include <cppunit/extensions/HelperMacros.h>

#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define SELDON_WITH_SUPERLU
#include "Verdandi.hxx"
#include "seldon/SeldonSolver.hxx"
#include "method/chi_2.cxx"
using namespace Verdandi;


class Chi2Test: public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE(Chi2Test);
    CPPUNIT_TEST(test_chi_2);
    CPPUNIT_TEST_SUITE_END();

protected:
    int Nx_;
    int Ny_;

public:
    void setUp()
    {
    }


    void tearDown()
    {
    }


    // Fills a matrix with ones.
    template <class MatrixType>
    void fill_matrix(int m, int n, MatrixType& M)
    {
        typedef typename MatrixType::value_type T;
        M.Reallocate(m, n);
        M.Fill(T(1));
    }


    // Fills a matrix with ones.
    template <class T, class Prop, class Allocator>
    void fill_matrix(int m, int n, Matrix<T, Prop, RowSparse, Allocator>& M)
    {
        Matrix<T, Prop, ArrayRowSparse, Allocator> Marray(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                Marray.Get(i, j) = T(1);
        Copy(Marray, M);
    }


    void test_chi_2()
    {
        int Nx[5] = {10, 10, 2, 2, 1};
        int Ny[5] = { 2, 10, 2, 1, 1};

        for (int i = 0; i < 5; i++)
        {
            Nx_ = Nx[i];
            Ny_ = Ny[i];

            compute_chi_2<Matrix<double, General, RowSparse>,
                Matrix<double, General, RowSparse>,
                Matrix<double, General, RowSparse>,
                Vector<double>,
                Matrix<double, General, RowSparse>,
                Vector<double> >();

            compute_chi_2<Matrix<double>, Matrix<double>, Matrix<double>,
                Vector<double>, Matrix<double>, Vector<double> >();

            compute_chi_2<Matrix<float>, Matrix<float>, Matrix<float>,
                Vector<float>, Matrix<float>, Vector<float> >();
        }
    }


    template <class StateErrorVariance,
              class ObservationOperator, class CrossedMatrix,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    void compute_chi_2()
    {
        {
            StateErrorVariance B(Nx_, Nx_);
            ObservationOperator H(Ny_, Nx_);
            CrossedMatrix tmp;
            ObservationVector y(Ny_);
            ObservationErrorVariance R(Ny_, Ny_);
            StateVector x(Nx_), analysis;

            B.SetIdentity();
            H.SetIdentity();
            R.SetIdentity();
            y.Fill();
            x.Fill(typename StateVector::value_type(1));

            double chi2 = 0;
            if (Ny_ == 1 || Ny_ == 2)
                chi2 = 1.;
            else
                chi2 = 1.
                    + double(((Ny_ - 2) * (Ny_ - 1) * (2 * Ny_ - 3)) / 6);
            CPPUNIT_ASSERT(chi_2(x, B, H, y, R) == chi2 / 2.);
        }
        {
            typedef typename ObservationVector::value_type T;

            StateErrorVariance B(Nx_, Nx_);
            ObservationOperator H;
            CrossedMatrix tmp;
            ObservationVector y(Ny_);
            ObservationErrorVariance R(Ny_, Ny_);
            StateVector x(Nx_), analysis;

            B.SetIdentity();
            fill_matrix(Ny_, Nx_, H);
            R.SetIdentity();
            y.Fill();
            T bg_value = 5;
            x.Fill(bg_value);

            ObservationVector innovation(Ny_), aux;
            for (int i = 0; i < Ny_; i++)
                innovation(i) = T(i) - bg_value * T(Nx_);
            aux = innovation;
            ObservationErrorVariance S;
            fill_matrix(Ny_, Ny_, S);
            for (int i = 0; i < Ny_; i++)
                for (int j = 0; j < Ny_; j++)
                    S.Val(i, j) *= double(Nx_);
            for (int i = 0; i < Ny_; i++)
                S.Val(i, i) += 1;
            GetAndSolveLU(S, aux);

            CPPUNIT_ASSERT(chi_2(x, B, H, y, R) == DotProd(innovation, aux));
        }
    }
};
