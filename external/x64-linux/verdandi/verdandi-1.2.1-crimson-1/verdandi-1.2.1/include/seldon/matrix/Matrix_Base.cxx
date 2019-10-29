// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_BASE_CXX

#include "Matrix_Base.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Allocator>
  inline Matrix_Base<T, Allocator>::Matrix_Base()
  {
    m_ = 0;
    n_ = 0;
    data_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j matrix, but data array is not initialized.
    \param i number of rows.
    \param j number of columns.
    \warning the data array is not allocated.
  */
  template <class T, class Allocator>
  inline Matrix_Base<T, Allocator>::Matrix_Base(int i, int j)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (i < 0 || j < 0)
      throw WrongDim("Matrix_Base::Matrix_Base(int, int)",
                     "Unable to define a matrix with size "
                     + to_str(i) + " by " + to_str(j) + ".");
#endif

    m_ = i;
    n_ = j;
    data_ = NULL;
  }


  //! Copy constructor.
  /*!
    \param A base matrix to be copied.
    \warning Only the length is copied.
  */
  template <class T, class Allocator>
  inline Matrix_Base<T, Allocator>::
  Matrix_Base(const Matrix_Base<T, Allocator>& A)
  {
    m_ = A.GetM();
    n_ = A.GetN();
    data_ = NULL;
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  /*!
    \note Memory for the data array is not freed.
  */
  template <class T, class Allocator>
  inline Matrix_Base<T, Allocator>::~Matrix_Base()
  {

  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of rows.
  /*!
    \return The number of rows.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetM() const
  {
    return m_;
  }


  //! Returns the number of columns.
  /*!
    \return The number of columns.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetN() const
  {
    return n_;
  }


  //! Returns the number of rows of the matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the matrix.
    \return The number of rows of the possibly-transposed matrix.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetM(const SeldonTranspose& status) const
  {
    if (status.NoTrans())
      return m_;
    else
      return n_;
  }


  //! Returns the number of columns of the matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the matrix.
    \return The number of columns of the possibly-transposed matrix.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetN(const SeldonTranspose& status) const
  {
    if (status.NoTrans())
      return n_;
    else
      return m_;
  }


#ifdef SELDON_WITH_BLAS
  //! Returns the number of rows of the matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the matrix.
    \return The number of rows of the possibly-transposed matrix.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetM(const CBLAS_TRANSPOSE& status) const
  {
    if (status == CblasNoTrans)
      return m_;
    else
      return n_;
  }
#endif


#ifdef SELDON_WITH_BLAS
  //! Returns the number of columns of the matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the matrix.
    \return The number of columns of the possibly-transposed matrix.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetN(const CBLAS_TRANSPOSE& status) const
  {
    if (status == CblasNoTrans)
      return n_;
    else
      return m_;
  }
#endif


  //! Returns the number of elements in the matrix.
  /*!
    Returns the number of elements in the matrix, i.e.
    the number of rows multiplied by the number of columns.
    \return The number of elements in the matrix.
  */
  template <class T, class Allocator>
  int Matrix_Base<T, Allocator>::GetSize() const
  {
    return m_ * n_;
  }


  //! Returns a pointer to the data array.
  /*!
    Returns a pointer to data, i.e. the data array 'data_'.
    \return A pointer to the data array.
  */
  template <class T, class Allocator>
  typename Matrix_Base<T, Allocator>::pointer
  Matrix_Base<T, Allocator>::GetData() const
  {
    return data_;
  }


  //! Returns a const pointer to the data array.
  /*!
    Returns a const pointer to data, i.e. the data array 'data_'.
    \return A const pointer to the data array.
  */
  template <class T, class Allocator>
  typename Matrix_Base<T, Allocator>::const_pointer
  Matrix_Base<T, Allocator>::GetDataConst() const
  {
    return reinterpret_cast<typename Matrix_Base<T,
      Allocator>::const_pointer>(data_);
  }


  //! Returns a pointer of type "void*" to the data array.
  /*!
    Returns a pointer of type "void*" to data, i.e. the data array 'data_'.
    \return A pointer of type "void*" to the data array.
  */
  template <class T, class Allocator>
  void* Matrix_Base<T, Allocator>::GetDataVoid() const
  {
    return reinterpret_cast<void*>(data_);
  }


  //! Returns a pointer of type "const void*" to the data array.
  /*!
    Returns a pointer of type "const void*" to data, i.e.
    the data array 'data_'.
    \return A const pointer of type "void*" to the data array.
  */
  template <class T, class Allocator>
  const void* Matrix_Base<T, Allocator>::GetDataConstVoid() const
  {
    return reinterpret_cast<const void*>(data_);
  }


  //! Returns the allocator of the matrix.
  /*!
    \return The allocator.
  */
  template <class T, class Allocator>
  Allocator& Matrix_Base<T, Allocator>::GetAllocator()
  {
    return allocator_;
  }


  // operator<< overloaded for matrices.
  /*!
    \param out output stream.
    \param A matrix to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Matrix<T, Prop, Storage, Allocator>& A)
  {

    for (int i = 0; i < A.GetM(); i++)
      {
	for (int j = 0; j < A.GetN(); j++)
	  out << A(i, j) << '\t';
	out << endl;
      }

    return out;

  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_BASE_CXX
#endif
