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


#ifndef VERDANDI_FILE_METHOD_SIGMA_POINT_CXX


#include "SigmaPoint.hxx"

#include "seldon/vector/VectorCollection.cxx"

namespace Verdandi
{


    //! Computes 'canonical' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'canonical' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
      \param[out] alpha_constant boolean to indicate if the coefficients alpha
      are constants.
    */
    template <class SigmaPoint>
    void
    ComputeCanonicalSigmaPoint(int p,
                               Vector<SigmaPoint, Collection>& sigma_point,
                               SigmaPoint& alpha, bool& alpha_constant)
    {
        typedef typename SigmaPoint::value_type value_type;

        int r;
        SigmaPoint working_vector;

        r = 2 * p;
        value_type sqrt_p = sqrt(value_type(p));

        for (int i = 0; i < p; i++)
        {
            working_vector.Reallocate(p);
            working_vector.Fill(value_type(0));
            working_vector(i) = sqrt_p;
            sigma_point.AddVector(working_vector);
            working_vector.Nullify();
        }

        for (int i = p; i < r; i++)
        {
            working_vector.Reallocate(p);
            working_vector.Fill(value_type(0));
            working_vector(r - i - 1) = -sqrt_p;
            sigma_point.AddVector(working_vector);
            working_vector.Nullify();
        }

        alpha.Reallocate(r);
        alpha.Fill(value_type(1) / value_type(r));

        alpha_constant = true;
    }


    //! Computes 'canonical' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'canonical' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
      \param[out] alpha_constant boolean to indicate if the coefficients alpha
      are constants.
    */
    template <class T,  template <class U> class Allocator>
    void ComputeCanonicalSigmaPoint(int p,
                                    Matrix<T, General, RowMajor,
                                    Allocator<T> >& sigma_point,
                                    Vector<T, VectFull, Allocator<T> >&
                                    alpha, bool& alpha_constant)
    {
        typedef Vector<T, VectFull, Allocator<T> > SigmaPoint;
        Vector<SigmaPoint, Collection> sigma_point_collection;
        int r;

        ComputeCanonicalSigmaPoint(p, sigma_point_collection,
                                   alpha, alpha_constant);

        r = sigma_point_collection.GetNvector();

        sigma_point.Reallocate(r, p);
        for (int i = 0; i < r; i++)
            SetRow(sigma_point_collection.GetVector(i), i, sigma_point);

        sigma_point_collection.Deallocate();
    }


    //! Computes 'canonical' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'canonical' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
      \param[out] alpha_constant boolean to indicate if the coefficients alpha
      are constants.
    */
    template <class T,  template <class U> class Allocator>
    void ComputeCanonicalSigmaPoint(int p,
                                    Matrix<T, General, RowSparse,
                                    Allocator<T> >& sigma_point,
                                    Vector<T, VectFull, Allocator<T> >&
                                    alpha, bool& alpha_constant)
    {
        Matrix<T, General, RowMajor, Allocator<T> > sigma_point_dense;
        ComputeCanonicalSigmaPoint(p, sigma_point_dense, alpha,
                                   alpha_constant);
        Matrix<T, General, ArrayRowSparse, Allocator<T> > sigma_point_array;
        ConvertDenseToArrayRowSparse(sigma_point_dense, sigma_point_array);
        Copy(sigma_point_array, sigma_point);
    }


    //! Computes 'star' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'star' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
      \param[out] alpha_constant boolean to indicate if the coefficients alpha
      are constants.
    */
    template <class SigmaPoint>
    void ComputeStarSigmaPoint(int p,
                               Vector<SigmaPoint, Collection>& sigma_point,
                               SigmaPoint& alpha, bool& alpha_constant)
    {
        typedef typename SigmaPoint::value_type value_type;

        ComputeCanonicalSigmaPoint(p, sigma_point, alpha, alpha_constant);

        value_type scal = sqrt(value_type(2 * p + 1) / value_type(2 * p));
        for (int i = 0; i < sigma_point.GetNvector(); i++)
            Mlt(scal, sigma_point.GetVector(i));

        SigmaPoint working_vector;
        working_vector.Reallocate(p);
        working_vector.Fill(value_type(0));
        sigma_point.AddVector(working_vector);
        working_vector.Nullify();

        alpha.Reallocate(2 * p + 1);
        alpha.Fill(value_type(1) / value_type(2 * p + 1));

        alpha_constant = true;
    }


    //! Computes 'star' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'star' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
      \param[out] alpha_constant boolean to indicate if the coefficients alpha
      are constants.
    */
    template <class T,  template <class U> class Allocator>
    void ComputeStarSigmaPoint(int p,
                               Matrix<T, General, RowMajor, Allocator<T> >&
                               sigma_point,
                               Vector<T, VectFull, Allocator<T> >&
                               alpha, bool& alpha_constant)
    {
        typedef Vector<T, VectFull, Allocator<T> > SigmaPoint;
        Vector<SigmaPoint, Collection> sigma_point_collection;
        int r;

        ComputeStarSigmaPoint(p, sigma_point_collection,
                              alpha, alpha_constant);

        r = sigma_point_collection.GetNvector();

        sigma_point.Reallocate(r, p);
        for (int i = 0; i < r; i++)
            SetRow(sigma_point_collection.GetVector(i), i, sigma_point);

        sigma_point_collection.Deallocate();
    }


    //! Computes 'star' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'star' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
      \param[out] alpha_constant boolean to indicate if the coefficients alpha
      are constants.
    */
    template <class T,  template <class U> class Allocator>
    void ComputeStarSigmaPoint(int p,
                               Matrix<T, General, RowSparse, Allocator<T> >&
                               sigma_point,
                               Vector<T, VectFull, Allocator<T> >&
                               alpha, bool& alpha_constant)
    {
        Matrix<T, General, RowMajor, Allocator<T> > sigma_point_dense;
        ComputeStarSigmaPoint(p, sigma_point_dense, alpha,
                              alpha_constant);
        Matrix<T, General, ArrayRowSparse, Allocator<T> > sigma_point_array;
        ConvertDenseToArrayRowSparse(sigma_point_dense, sigma_point_array);
        Copy(sigma_point_array, sigma_point);
    }


    //! Computes 'simplex' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'simplex' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
    */
    template <class SigmaPoint>
    void ComputeSimplexSigmaPoint(int p,
                                  Vector<SigmaPoint, Collection>& sigma_point,
                                  SigmaPoint& alpha, bool& alpha_constant)
    {
        typedef typename SigmaPoint::value_type value_type;
        SigmaPoint working_vector;
        int r;
        r = p + 1;
        Matrix<value_type> working_matrix(p, r);
        working_matrix.Fill(value_type(0));

        value_type scal, beta, sqrt_p;
        beta = value_type(p) / (value_type(p) + 1);
        sqrt_p = sqrt(value_type(p));

        scal = value_type(1) / sqrt(value_type(2) * beta);
        working_matrix(0, 0) = -scal;
        working_matrix(0, 1) = scal;

        for (int i = 1; i < p; i++)
        {
            scal = value_type(1)
                / sqrt(beta * value_type(i + 1) * value_type(i + 2));
            for (int j = 0; j < i + 1; j++)
                working_matrix(i, j) = scal;
            working_matrix(i, i + 1) = - value_type(i + 1) * scal;
        }

        for (int i = 0; i < r; i++)
        {
            GetCol(working_matrix, i, working_vector);
            Mlt(sqrt_p, working_vector);
            sigma_point.AddVector(working_vector);
            working_vector.Nullify();
        }

        alpha.Reallocate(r);
        alpha.Fill(value_type(1) / value_type(r));

        alpha_constant = true;
    }


    //! Computes 'simplex' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'simplex' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
    */
    template <class T,  template <class U> class Allocator>
    void ComputeSimplexSigmaPoint(int p,
                                  Matrix<T, General, RowMajor, Allocator<T> >&
                                  sigma_point,
                                  Vector<T, VectFull, Allocator<T> >&
                                  alpha, bool& alpha_constant)
    {
        typedef Vector<T, VectFull, Allocator<T> > SigmaPoint;
        Vector<SigmaPoint, Collection> sigma_point_collection;
        int r;

        ComputeSimplexSigmaPoint(p, sigma_point_collection,
                                 alpha, alpha_constant);

        r = sigma_point_collection.GetNvector();

        sigma_point.Reallocate(r, p);
        for (int i = 0; i < r; i++)
            SetRow(sigma_point_collection.GetVector(i), i, sigma_point);

        sigma_point_collection.Deallocate();
    }


    //! Computes 'simplex' sigma-points.
    /*!
      \param[in] p dimension of the model state.
      \param[out] sigma_point 'simplex' sigma-points.
      \param[out] alpha coefficient vector associated with sigma-points.
    */
    template <class T,  template <class U> class Allocator>
    void ComputeSimplexSigmaPoint(int p,
                                  Matrix<T, General, RowSparse, Allocator<T> >
                                  &sigma_point,
                                  Vector<T, VectFull, Allocator<T> >&
                                  alpha, bool& alpha_constant)
    {
        Matrix<T, General, RowMajor, Allocator<T> > sigma_point_dense;
        ComputeSimplexSigmaPoint(p, sigma_point_dense, alpha,
                                 alpha_constant);
        Matrix<T, General, ArrayRowSparse, Allocator<T> > sigma_point_array;
        ConvertDenseToArrayRowSparse(sigma_point_dense, sigma_point_array);
        Copy(sigma_point_array, sigma_point);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_SIGMA_POINT_CXX
#endif
