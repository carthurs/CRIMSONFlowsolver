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


#ifndef VERDANDI_FILE_ERROR_DIAGONALMATRIX_CXX


#include "DiagonalMatrix.hxx"


namespace Verdandi
{


    /////////////////
    // CONSTRUCTOR //
    /////////////////


    //! Default constructor.
    /*!
      \param[in] dimension size of the matrix.
      \param[in] variance value on the diagonal.
    */
    template <class T>
    DiagonalMatrix<T>
    ::DiagonalMatrix(int dimension, T variance):
        dimension_(dimension), variance_(variance)
    {
    }


    ////////////
    // ACCESS //
    ////////////


    //! Returns the number of rows.
    /*!
      \return The number of rows.
    */
    template <class T>
    inline int DiagonalMatrix<T>::GetM() const
    {
        return dimension_;
    }


    //! Returns the number of columns.
    /*!
      \return The number of columns.
    */
    template <class T>
    inline int DiagonalMatrix<T>::GetN() const
    {
        return dimension_;
    }


    //! Access to one matrix entry.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The value of element (\a i, \a j).
    */
    template <class T>
    T DiagonalMatrix<T>::operator()(int i, int j) const
    {
        if (i == j)
            return variance_;
        else
            return T(0);
    }


    //! Access to a row.
    /*!
      \param[in] i row index.
      \param[out] row the row #\a i.
    */
    template <class T>
    void DiagonalMatrix<T>::GetRow(int i, Vector<T>& row) const
    {
        row.Reallocate(dimension_);
        row.Zero();
        row(i) = variance_;
    }


    //! Access to a column.
    /*!
      \param[in] j column index.
      \param[out] column the column #\a j.
    */
    template <class T>
    void DiagonalMatrix<T>::GetCol(int j, Vector<T>& column) const
    {
        GetRow(j, column);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_ERROR_DIAGONALMATRIX_CXX
#endif
