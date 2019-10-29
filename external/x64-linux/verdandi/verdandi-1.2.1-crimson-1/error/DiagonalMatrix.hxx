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


#ifndef VERDANDI_FILE_ERROR_DIAGONALMATRIX_HXX


namespace Verdandi
{


    //! Diagonal covariance error matrix.
    /*! \brief This class implements a covariance error matrix that is
      diagonal.
    */
    template <class T>
    class DiagonalMatrix
    {
    protected:
        //! Size of the matrix.
        int dimension_;

        //! Variance.
        T variance_;

    public:
        // Constructor.
        DiagonalMatrix(int dimension, T variance);

        // Access.
        int GetM() const;
        int GetN() const;
        T operator()(int i, int j) const;
        void GetRow(int i, Vector<T>& row) const;
        void GetCol(int i, Vector<T>& column) const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_ERROR_DIAGONALMATRIX_HXX
#endif
