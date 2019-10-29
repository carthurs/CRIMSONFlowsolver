// Copyright (C) 2009 Vivien Mallet
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


#ifndef SELDON_FILE_SUBMATRIX_BASE_CXX


#include "SubMatrix_Base.hxx"


namespace Seldon
{


  ////////////////////
  // SUBMATRIX_BASE //
  ////////////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Main constructor.
  template <class T, class Prop, class M, class Allocator>
  inline SubMatrix_Base<T, Prop, M, Allocator>
  ::SubMatrix_Base(M& A, Vector<int>& row_list, Vector<int>& column_list):
    Matrix_Base<T, Allocator>(row_list.GetLength(), column_list.GetLength()),
    matrix_(&A), row_list_(row_list), column_list_(column_list)
  {
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class M, class Allocator>
  inline SubMatrix_Base<T, Prop, M, Allocator>::~SubMatrix_Base()
  {
  }


  /************************
   * ACCESS AND AFFECTION *
   ************************/


  //! Access operator.
  /*!
    Returns the value of element (\a i, \a j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the sub-matrix.
  */
  template <class T, class Prop, class M, class Allocator>
  inline typename SubMatrix_Base<T, Prop, M, Allocator>::access_type
  SubMatrix_Base<T, Prop, M, Allocator>::operator() (int i, int j)
  {
    return (*this->matrix_)(this->row_list_(i), this->column_list_(j));
  }


  //! Access operator.
  /*!
    Returns the value of element (\a i, \a j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the sub-matrix.
  */
  template <class T, class Prop, class M, class Allocator>
  inline typename SubMatrix_Base<T, Prop, M, Allocator>::const_access_type
  SubMatrix_Base<T, Prop, M, Allocator>::operator() (int i, int j) const
  {
    return (*this->matrix_)(this->row_list_(i), this->column_list_(j));
  }


  //! Access operator.
  /*!
    Returns the value of element (\a i, \a j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the sub-matrix.
  */
  template <class T, class Prop, class M, class Allocator>
  inline typename SubMatrix_Base<T, Prop, M, Allocator>::entry_type&
  SubMatrix_Base<T, Prop, M, Allocator>::Val(int i, int j)
  {
    return this->matrix_->Val(this->row_list_(i), this->column_list_(j));
  }


  //! Access operator.
  /*!
    Returns the value of element (\a i, \a j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the sub-matrix.
  */
  template <class T, class Prop, class M, class Allocator>
  inline const typename SubMatrix_Base<T, Prop, M, Allocator>::entry_type&
  SubMatrix_Base<T, Prop, M, Allocator>::Val(int i, int j) const
  {
    return this->matrix_->Val(this->row_list_(i), this->column_list_(j));
  }


  /*****************
   * BASIC METHODS *
   *****************/


  //! Returns the number of rows.
  /*!
    \return The number of rows.
  */
  template <class T, class Prop, class M, class Allocator>
  inline int SubMatrix_Base<T, Prop, M, Allocator>::GetM() const
  {
    return row_list_.GetLength();
  }


  //! Returns the number of columns.
  /*!
    \return The number of columns.
  */
  template <class T, class Prop, class M, class Allocator>
  inline int SubMatrix_Base<T, Prop, M, Allocator>::GetN() const
  {
    return column_list_.GetLength();
  }


  //! Returns the number of rows of the sub-matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the sub-matrix.
    \return The number of rows of the possibly-transposed sub-matrix.
  */
  template <class T, class Prop, class M, class Allocator>
  inline int SubMatrix_Base<T, Prop, M, Allocator>
  ::GetM(const SeldonTranspose& status) const
  {
    if (status.NoTrans())
      return row_list_.GetLength();
    else
      return column_list_.GetLength();
  }


  //! Returns the number of columns of the sub-matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the sub-matrix.
    \return The number of columns of the possibly-transposed sub-matrix.
  */
  template <class T, class Prop, class M, class Allocator>
  inline int SubMatrix_Base<T, Prop, M, Allocator>
  ::GetN(const SeldonTranspose& status) const
  {
    if (status.NoTrans())
      return column_list_.GetLength();
    else
      return row_list_.GetLength();
  }


  //! Prints a matrix on screen.
  /*! Displays all elements on the standard output, in text format. Each row
    is displayed on a single line, and the elements of a row are delimited by
    tabulations.
  */
  template <class T, class Prop, class M, class Allocator>
  void SubMatrix_Base<T, Prop, M, Allocator>::Print() const
  {
    for (int i = 0; i < this->GetM(); i++)
      {
	for (int j = 0; j < this->GetN(); j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


} // namespace Seldon.


#define SELDON_FILE_SUBMATRIX_BASE_CXX
#endif
