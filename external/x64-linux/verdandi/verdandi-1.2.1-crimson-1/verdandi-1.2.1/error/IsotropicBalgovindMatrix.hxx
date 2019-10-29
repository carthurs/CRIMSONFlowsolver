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


#ifndef VERDANDI_FILE_ERROR_ISOTROPICBALGOVINDMATRIX_HXX


namespace Verdandi
{


    //! This class defines a covariance matrix in isotropic Balgovind form.
    /*! In a covariance matrix in isotropic Balgovind form, the covariance
      between two points only depends on the Euclidean distances between the
      points. For example, in the 2D case, if the entry \f$(i, j)\f$ of the
      matrix is the covariance between the values at \f$(x_i, y_i)\f$ and
      \f$(x_j, y_j)\f$, its value will be:

      \f[B_{i, j} = v \left(1 + \frac{\sqrt{(x_j - x_i)^2 + (y_j -
      y_i)^2}}{L}\right) \exp\left(-\frac{\sqrt{(x_j - x_i)^2 + (y_j -
      y_i)^2}}{L}\right)\f]

      where \f$v\f$ is a variance, and \f$L\f$ is the decorrelation length.
    */
    template <class T>
    class IsotropicBalgovindMatrix
    {
    protected:
        //! Coordinates of the domain first cell center.
        Vector<T> min_;
        //! Space steps.
        Vector<T> step_;
        //! Number of points along each dimension.
        Vector<int> N_;

        //! Length scale.
        T length_scale_;

        //! Variance.
        T variance_;

        //! Matrix dimension.
        int dimension_;

    public:
        // Constructors.
        IsotropicBalgovindMatrix(T x_min, T delta_x, int Nx,
                                 T y_min, T delta_y, int Ny,
                                 T length_scale, T variance);
        IsotropicBalgovindMatrix(T x_min, T delta_x, int Nx,
                                 T y_min, T delta_y, int Ny,
                                 T z_min, T delta_z, int Nz,
                                 T length_scale, T variance);

        // Access.
        int GetM() const;
        int GetN() const;
        T operator()(int i, int j) const;
        void GetRow(int i, Vector<T>& row) const;
        void GetCol(int i, Vector<T>& column) const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_ERROR_ISOTROPICBALGOVINDMATRIX_HXX
#endif
