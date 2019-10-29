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


#ifndef SELDON_FILE_SUBMATRIX_CXX


#include "SubMatrix.hxx"


namespace Seldon
{


  ////////////////////////
  // MATRIX<SUBSTORAGE> //
  ////////////////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Main constructor.
  template <class T, class Prop, class M, class Allocator>
  inline Matrix<T, Prop, SubStorage<M>, Allocator>
  ::Matrix(M& A, Vector<int>& row_list, Vector<int>& column_list):
    SubMatrix_Base<T, Prop, M, Allocator>(A, row_list, column_list)
  {
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class M, class Allocator>
  inline Matrix<T, Prop, SubStorage<M>, Allocator>::~Matrix()
  {
  }


  ///////////////
  // SUBMATRIX //
  ///////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Main constructor.
  template <class M>
  inline SubMatrix<M>
  ::SubMatrix(M& A, Vector<int>& row_list, Vector<int>& column_list):
    Matrix<typename M::value_type, typename M::property,
    SubStorage<M>, typename M::allocator>(A, row_list, column_list)
  {
  }


} // namespace Seldon.


#define SELDON_FILE_SUBMATRIX_CXX
#endif
