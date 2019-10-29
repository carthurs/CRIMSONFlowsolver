// Copyright (C) 2008, INRIA
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


#ifndef VERDANDI_FILE_ERROR_ISOTROPICBALGOVINDMATRIX_CXX


#include "IsotropicBalgovindMatrix.hxx"


namespace Verdandi
{


    //////////////////
    // CONSTRUCTORS //
    //////////////////


    //! Constructor for 2D cases.
    /*! This constructors builds an isotropic Balgovind matrix for 2D regular
      grids.
      \param x_min abscissa of the center of the lower-left grid cell.
      \param delta_x step along x.
      \param Nx number of cells along x.
      \param y_min ordinate of the center of the lower-left grid cell.
      \param delta_y step along x.
      \param Ny number of cells along y.
      \param length_scale decorrelation length.
      \param variance variance.
    */
    template <class T>
    IsotropicBalgovindMatrix<T>
    ::IsotropicBalgovindMatrix(T x_min, T delta_x, int Nx,
                               T y_min, T delta_y, int Ny,
                               T length_scale, T variance):
        min_(2), step_(2), N_(2), length_scale_(length_scale),
        variance_(variance)
    {
        min_(0) = x_min;
        step_(0) = delta_x;
        N_(0) = Nx;
        min_(1) = y_min;
        step_(1) = delta_y;
        N_(1) = Ny;

        dimension_ = Nx * Ny;
    }


    //! Constructor for 3D cases.
    /*! This constructors builds an isotropic Balgovind matrix for 3D regular
      grids.
      \param x_min abscissa of the center of the lower-left grid cell.
      \param delta_x step along x.
      \param Nx number of cells along x.
      \param y_min ordinate of the center of the lower-left grid cell.
      \param delta_y step along x.
      \param Ny number of cells along y.
      \param z_min applicate of the center of the bottom grid cells.
      \param delta_z step along z.
      \param Nz number of cells along z.
      \param length_scale decorrelation length.
      \param variance variance.
    */
    template <class T>
    IsotropicBalgovindMatrix<T>
    ::IsotropicBalgovindMatrix(T x_min, T delta_x, int Nx,
                               T y_min, T delta_y, int Ny,
                               T z_min, T delta_z, int Nz,
                               T length_scale, T variance):
        min_(3), step_(3), N_(3), length_scale_(length_scale),
        variance_(variance)
    {
        min_(0) = x_min;
        step_(0) = delta_x;
        N_(0) = Nx;
        min_(1) = y_min;
        step_(1) = delta_y;
        N_(1) = Ny;
        min_(2) = z_min;
        step_(2) = delta_z;
        N_(2) = Nz;

        dimension_ = Nx * Ny * Nz;
    }


    ////////////
    // ACCESS //
    ////////////


    //! Returns the number of rows.
    /*!
      \return The number of rows.
    */
    template <class T>
    inline int IsotropicBalgovindMatrix<T>::GetM() const
    {
        return dimension_;
    }


    //! Returns the number of columns.
    /*!
      \return The number of columns.
    */
    template <class T>
    inline int IsotropicBalgovindMatrix<T>::GetN() const
    {
        return dimension_;
    }


    //! Access to one matrix entry.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The value of the entry (i, j) of the matrix.
    */
    template <class T>
    T IsotropicBalgovindMatrix<T>::operator()(int i, int j) const
    {
        Vector<int> position_i, position_j;
        get_position(i, N_, position_i);
        get_position(j, N_, position_j);
        T ratio, distance(0), tmp;
        for (int d = 0; d < N_.GetLength(); d++)
        {
            tmp = T(position_i(d) - position_j(d)) * step_(d);
            distance += tmp * tmp;
        }
        ratio = sqrt(distance) / length_scale_;
        return variance_ * (1. + ratio) * exp(-ratio);
    }


    //! Access to a row.
    /*!
      \param[in] i row index.
      \param[out] row the \a i-th row.
    */
    template <class T>
    void IsotropicBalgovindMatrix<T>::GetRow(int i, Vector<T>& row) const
    {
        row.Reallocate(dimension_);
        for (int d = 0; d < dimension_; d++)
            row(d) = (*this)(i, d);
    }


    //! Access to a column.
    /*!
      \param[in] j column index.
      \param[out] column the \a j-th row.
    */
    template <class T>
    void IsotropicBalgovindMatrix<T>::GetCol(int j, Vector<T>& column) const
    {
        GetRow(j, column);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_ERROR_ISOTROPICBALGOVINDMATRIX_CXX
#endif
