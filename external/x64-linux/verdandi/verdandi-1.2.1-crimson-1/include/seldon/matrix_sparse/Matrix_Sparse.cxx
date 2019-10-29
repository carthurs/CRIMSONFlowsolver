// Copyright (C) 2001-2009 Vivien Mallet
// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_SPARSE_CXX

#include "Matrix_Sparse.hxx"

#include "../vector/Functions_Arrays.cxx"

#include <set>

#if defined(SELDON_WITH_MKL)
#include "mkl_service.h"
#endif

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>::Matrix_Sparse():
    Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>::Matrix_Sparse(int i,
								   int j):
    Matrix_Base<T, Allocator>()
  {
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;

    Reallocate(i, j);
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i by j , with nz non-zero elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>::
  Matrix_Sparse(int i, int j, int nz):
    Matrix_Base<T, Allocator>()
  {
    this->nz_ = 0;
    ind_ = NULL;
    ptr_ = NULL;

    Reallocate(i, j, nz);
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>::
  Matrix_Sparse(int i, int j,
		Vector<T, Storage0, Allocator0>& values,
		Vector<int, Storage1, Allocator1>& ptr,
		Vector<int, Storage2, Allocator2>& ind):
    Matrix_Base<T, Allocator>(i, j)
  {
    nz_ = values.GetLength();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (ind.GetLength() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_Sparse::Matrix_Sparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_Sparse::Matrix_Sparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column")
		       + string(" start indices (plus the number")
		       + " of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns (" + to_str(i) + " by "
		       + to_str(j) + " matrix).");
      }

    if (nz_ > 0
        && (j == 0
            || static_cast<long int>(nz_-1) / static_cast<long int>(j)
            >= static_cast<long int>(i)))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_Sparse::Matrix_Sparse(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }


  //! Copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>::
  Matrix_Sparse(const Matrix_Sparse<T, Prop, Storage, Allocator>& A):
    Matrix_Base<T, Allocator>(A)
  {
    this->m_ = 0;
    this->n_ = 0;
    this->nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>::~Matrix_Sparse()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (ptr_ != NULL)
	  {
#if defined(SELDON_WITH_MKL)
	    mkl_free(ptr_);
#else
	    free(ptr_);
#endif
	    ptr_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (ind_ != NULL)
	  {
#if defined(SELDON_WITH_MKL)
		mkl_free(ind_);
#else
	    free(ind_);
#endif
	    ind_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_, nz_);
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->nz_ = 0;
	this->data_ = NULL;
      }
#endif

    this->nz_ = 0;
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix
    is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Sparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_Sparse();
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by
    'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  Vector<T, Storage0, Allocator0>& values,
	  Vector<int, Storage1, Allocator1>& ptr,
	  Vector<int, Storage2, Allocator2>& ind)
  {
    this->Clear();
    this->m_ = i;
    this->n_ = j;
    this->nz_ = values.GetLength();

#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.

    if (ind.GetLength() != nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_Sparse::SetData(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(nz_) + " values but "
		       + to_str(ind.GetLength()) + " row or column indices.");
      }

    if (ptr.GetLength()-1 != Storage::GetFirst(i, j))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_Sparse::SetData(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices contains ")
		       + to_str(ptr.GetLength()-1) + string(" row or column")
		       + string(" start indices (plus the number")
		       + " of non-zero entries) but there are "
		       + to_str(Storage::GetFirst(i, j))
		       + " rows or columns (" + to_str(i) + " by "
		       + to_str(j) + " matrix).");
      }

    if (nz_ > 0
        && (j == 0
            || static_cast<long int>(nz_-1) / static_cast<long int>(j)
            >= static_cast<long int>(i)))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim(string("Matrix_Sparse::SetData(int, int, ")
		       + "const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(values.GetLength())
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

    this->ptr_ = ptr.GetData();
    this->ind_ = ind.GetData();
    this->data_ = values.GetData();

    ptr.Nullify();
    ind.Nullify();
    values.Nullify();
  }


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by arrays
    'values' (values), 'ptr' (pointers) and 'ind' (indices).
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero entries.
    \param values values of non-zero entries.
    \param ptr row or column start indices.
    \param ind row or column indices.
    \warning On exit, arrays 'values', 'ptr' and 'ind' are managed by the
    matrix.
    For example, it means that the destructor will released those arrays;
    therefore, the user mustn't release those arrays.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::SetData(int i, int j, int nz,
	    typename Matrix_Sparse<T, Prop, Storage, Allocator>
	    ::pointer values,
	    int* ptr, int* ind)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = j;

    this->nz_ = nz;

    this->data_ = values;
    ind_ = ind;
    ptr_ = ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Sparse<T, Prop, Storage, Allocator>::Nullify()
  {
    this->data_ = NULL;
    this->m_ = 0;
    this->n_ = 0;
    nz_ = 0;
    ptr_ = NULL;
    ind_ = NULL;
  }


  //! Initialization of an empty sparse matrix with i rows and j columns
  /*!
    \param i number of rows
    \param j number of columns
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Reallocate(int i, int j)
  {
    // clearing previous entries
    Clear();

    this->m_ = i;
    this->n_ = j;

    // we try to allocate ptr_
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<int*>( calloc(Storage::GetFirst(i, j)+1,
					      sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Reallocate(int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (Storage::GetFirst(i, j)+1) )
		     + " bytes to store " + to_str(Storage::GetFirst(i, j)+1)
		     + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

    // then filing ptr_ with 0
    for (int k = 0; k <= Storage::GetFirst(i, j); k++)
      ptr_[k] = 0;

  }


  //! Initialization of a sparse matrix with i rows and j columns
  /*!
    \param i number of rows
    \param j number of columns
    \param nz number of non-zero entries
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Reallocate(int i, int j, int nz)
  {
    // clearing previous entries
    Clear();

    this->nz_ = nz;
    this->m_ = i;
    this->n_ = j;

#ifdef SELDON_CHECK_DIMENSIONS
    if (nz_ < 0)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		       "Invalid number of non-zero elements: " + to_str(nz)
                       + ".");
      }
    if (nz_ > 0
        && (j == 0
            || static_cast<long int>(nz_-1) / static_cast<long int>(j)
            >= static_cast<long int>(i)))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<int*>( calloc(Storage::GetFirst(i, j)+1,
					      sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (Storage::GetFirst(i, j)+1) )
		     + " bytes to store " + to_str(Storage::GetFirst(i, j)+1)
		     + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ind_ = reinterpret_cast<int*>( calloc(nz_, sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
    free(ptr_);
#endif
	ptr_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz)
		     + " row or column indices, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate(nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
#if defined(SELDON_WITH_MKL)
	mkl_free(ind_);
#else
	free(ind_);
#endif
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
#if defined(SELDON_WITH_MKL)
	mkl_free(ind_);
#else
	free(ind_);
#endif
	ind_ = NULL;
      }
    if (this->data_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz) + " values, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif
  }


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    if (Storage::GetFirst(i, j) < Storage::GetFirst(this->m_, this->n_))
      Resize(i, j, ptr_[Storage::GetFirst(i, j)]);
    else
      Resize(i, j, nz_);
  }


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \param nz number of non-zero elements.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::Resize(int i, int j, int nz)
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (nz < 0)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_Sparse::Resize(int, int, int)",
		       "Invalid number of non-zero elements: " + to_str(nz)
                       + ".");
      }
    if (nz > 0
        && (j == 0
            || static_cast<long int>(nz_-1) / static_cast<long int>(j)
            >= static_cast<long int>(i)))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_Sparse::Resize(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

    if (nz != nz_)
      {
        // trying to resize ind_ and data_
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif

            ind_
              = reinterpret_cast<int*>(realloc(reinterpret_cast<void*>(ind_),
                                               nz*sizeof(int)));

#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
#if defined(SELDON_WITH_MKL)
	        mkl_free(ptr_);
#else
	        free(ptr_);
#endif
            ptr_ = NULL;
            ind_ = NULL;
            this->data_ = NULL;
          }
        if (ind_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
#if defined(SELDON_WITH_MKL)
	        mkl_free(ptr_);
#else
	        free(ptr_);
#endif
            ptr_ = NULL;
            this->data_ = NULL;
          }
        if (ind_ == NULL && i != 0 && j != 0)
          throw NoMemory("Matrix_Sparse::Resize(int, int, int)",
                         string("Unable to allocate ") + to_str(sizeof(int) * nz)
                         + " bytes to store " + to_str(nz)
                         + " row or column indices, for a "
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

        Vector<T, VectFull, Allocator> val;
        val.SetData(nz_, this->data_);
        val.Resize(nz);

        this->data_ = val.GetData();
        nz_ = nz;
        val.Nullify();
      }


    if (Storage::GetFirst(this->m_, this->n_) != Storage::GetFirst(i, j))
      {
#ifdef SELDON_CHECK_MEMORY
        try
          {
#endif
            // trying to resize ptr_
            ptr_ = reinterpret_cast<int*>( realloc(ptr_, (Storage::GetFirst(i, j)+1)*
                                                   sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
          }
        catch (...)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
            ptr_ = NULL;
            ind_ = NULL;
            this->data_ = NULL;
          }
        if (ptr_ == NULL)
          {
            this->m_ = 0;
            this->n_ = 0;
            nz_ = 0;
            ind_ = NULL;
            this->data_ = NULL;
          }
        if (ptr_ == NULL && i != 0 && j != 0)
          throw NoMemory("Matrix_Sparse::Resize(int, int)",
                         string("Unable to allocate ")
                         + to_str(sizeof(int) * (Storage::GetFirst(i, j)+1) )
                         + " bytes to store " + to_str(Storage::GetFirst(i, j)+1)
                         + " row or column start indices, for a "
                         + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

        // then filing last values of ptr_ with nz_
        for (int k = Storage::GetFirst(this->m_, this->n_);
             k <= Storage::GetFirst(i, j); k++)
          ptr_[k] = this->nz_;
      }

    this->m_ = i;
    this->n_ = j;
  }


  //! Copies a matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Sparse<T, Prop, Storage, Allocator>::
  Copy(const Matrix_Sparse<T, Prop, Storage, Allocator>& A)
  {
    this->Clear();
    int nz = A.nz_;
    int i = A.m_;
    int j = A.n_;
    this->nz_ = nz;
    this->m_ = i;
    this->n_ = j;
    if ((i == 0)||(j == 0))
      {
	this->m_ = 0;
	this->n_ = 0;
	this->nz_ = 0;
	return;
      }

#ifdef SELDON_CHECK_DIMENSIONS
    if (nz_ > 0
        && (j == 0
            || static_cast<long int>(nz_-1) / static_cast<long int>(j)
            >= static_cast<long int>(i)))
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
	throw WrongDim("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		       string("There are more values (") + to_str(nz)
		       + " values) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ptr_ = reinterpret_cast<int*>( calloc(Storage::GetFirst(i, j)+1,
					      sizeof(int)) );
	memcpy(this->ptr_, A.ptr_,
               (Storage::GetFirst(i, j) + 1) * sizeof(int));
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ptr_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (Storage::GetFirst(i, j)+1) )
		     + " bytes to store " + to_str(Storage::GetFirst(i, j)+1)
		     + " row or column start indices, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	ind_ = reinterpret_cast<int*>( calloc(nz_, sizeof(int)) );
	memcpy(this->ind_, A.ind_, nz_ * sizeof(int));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	nz_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
	this->data_ = NULL;
      }
    if (ind_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz)
		     + " row or column indices, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate(nz_, this);
	this->allocator_.memorycpy(this->data_, A.data_, nz_);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
#if defined(SELDON_WITH_MKL)
	mkl_free(ind_);
#else
	free(ind_);
#endif
	ind_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
#if defined(SELDON_WITH_MKL)
	mkl_free(ptr_);
#else
	free(ptr_);
#endif
	ptr_ = NULL;
#if defined(SELDON_WITH_MKL)
	mkl_free(ind_);
#else
	free(ind_);
#endif
	ind_ = NULL;
      }
    if (this->data_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Sparse::Matrix_Sparse(int, int, int)",
		     string("Unable to allocate ") + to_str(sizeof(int) * nz)
		     + " bytes to store " + to_str(nz) + " values, for a "
		     + to_str(i) + " by " + to_str(j) + " matrix.");
#endif

  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of non-zero elements.
  /*!
    \return The number of non-zero elements.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_Sparse<T, Prop, Storage, Allocator>::GetNonZeros() const
  {
    return nz_;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_Sparse<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return nz_;
  }


  //! Returns (row or column) start indices.
  /*!
    Returns the array ('ptr_') of start indices.
    \return The array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_Sparse<T, Prop, Storage, Allocator>::GetPtr() const
  {
    return ptr_;
  }


  //! Returns (row or column) indices of non-zero entries.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries. This array defines non-zero entries
    indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_Sparse<T, Prop, Storage, Allocator>::GetInd() const
  {
    return ind_;
  }


  //! Returns the length of the array of start indices.
  /*!
    \return The length of the array of start indices.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_Sparse<T, Prop, Storage, Allocator>::GetPtrSize() const
  {
    return (Storage::GetFirst(this->m_, this->n_) + 1);
  }


  //! Returns the length of the array of (column or row) indices.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries. This array defines non-zero entries indices
    if coupled with (column or row) start indices.
    \return The length of the array of (column or row) indices.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_Sparse<T, Prop, Storage, Allocator>::GetIndSize() const
  {
    return nz_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Sparse<T, Prop, Storage, Allocator>::value_type
  Matrix_Sparse<T, Prop, Storage, Allocator>::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Sparse::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Sparse::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int k,l;
    int a,b;

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      return T(0);

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (ind_[k] < l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      return T(0);
  }


  //! Access method.
  /*! Returns the value of element (\a i, \a j) if it can be returned as a
    reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Sparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_Sparse<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Sparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Sparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int k, l;
    int a, b;

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_Sparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (ind_[k] < l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      throw WrongArgument("Matrix_Sparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");
  }


  //! Access method.
  /*! Returns the value of element (\a i, \a j) if it can be returned as a
    reference.
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    \throw WrongArgument No reference can be returned because the element is a
    zero entry (not stored in the matrix).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline
  const typename Matrix_Sparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_Sparse<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Sparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Sparse::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int k, l;
    int a, b;

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a == b)
      throw WrongArgument("Matrix_Sparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");

    l = Storage::GetSecond(i, j);

    for (k = a; (k < b-1) && (ind_[k] < l); k++);

    if (ind_[k] == l)
      return this->data_[k];
    else
      throw WrongArgument("Matrix_Sparse::Val(int, int)",
                          "No reference to element (" + to_str(i) + ", "
                          + to_str(j)
                          + ") can be returned: it is a zero entry.");
  }


  //! Access method.
  /*! Returns reference to element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
    If the element does not belong to sparsity pattern of the matrix,
    the matrix is resized.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Sparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_Sparse<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Sparse::Get(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Sparse::Get(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int k, l;
    int a, b;

    a = ptr_[Storage::GetFirst(i, j)];
    b = ptr_[Storage::GetFirst(i, j) + 1];

    if (a < b)
      {
        l = Storage::GetSecond(i, j);

        for (k = a; (k < b) && (ind_[k] < l); k++);

        if ( (k < b) && (ind_[k] == l))
          return this->data_[k];
      }
    else
      k = a;

    // adding a non-zero entry
    Resize(this->m_, this->n_, nz_+1);

    for (int m = Storage::GetFirst(i, j)+1;
         m <= Storage::GetFirst(this->m_, this->n_); m++)
      ptr_[m]++;

    for (int m = nz_-1; m >= k+1; m--)
      {
        ind_[m] = ind_[m-1];
        this->data_[m] = this->data_[m-1];
      }

    ind_[k] = Storage::GetSecond(i, j);

    // value of new non-zero entry is set to 0
    SetComplexZero(this->data_[k]);

    return this->data_[k];
  }


  //! Access method.
  /*! Returns reference to element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \return Element (\a i, \a j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename Matrix_Sparse<T, Prop, Storage, Allocator>::value_type&
  Matrix_Sparse<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return Val(i, j);
  }


  //! Add a value to a non-zero entry.
  /*! This function adds \a val to the element (\a i, \a j), provided that
    this element is already a non-zero entry. Otherwise
    a non-zero entry is inserted equal to \a val.
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val value to be added to the element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    Get(i, j) += val;
  }


  //! Sets an element (i, j) to a value
  /*! This function sets \a val to the element (\a i, \a j)
    \param[in] i row index.
    \param[in] j column index.
    \param[in] val A(i, j) = val
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    Get(i, j) = val;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Sparse<T, Prop, Storage, Allocator>&
  Matrix_Sparse<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_Sparse<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Resets all non-zero entries to 0-value.
  /*! The sparsity pattern remains unchanged. */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Zero()
  {
    this->allocator_.memoryset(this->data_, char(0),
			       this->nz_ * sizeof(value_type));
  }


  //! Sets the matrix to identity.
  /*! This method fills the diagonal of the matrix with ones. It can be
    applied to non square matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::SetIdentity()
  {
    int m = this->m_;
    int n = this->n_;
    int nz = min(m, n);

    if (nz == 0)
      return;

    Clear();

    Vector<T, VectFull, Allocator> values(nz);
    Vector<int, VectFull, CallocAlloc<int> > ptr(Storage::GetFirst(m, n) + 1);
    Vector<int, VectFull, CallocAlloc<int> > ind(nz);

    T one; SetComplexOne(one);
    values.Fill(one);
    ind.Fill();
    int i;
    for (i = 0; i < nz + 1; i++)
      ptr(i) = i;
    for (i = nz + 1; i < ptr.GetLength(); i++)
      ptr(i) = nz;

    SetData(m, n, values, ptr, ind);
  }


  //! Fills the non-zero entries with 0, 1, 2, ...
  /*! On exit, the non-zero entries are 0, 1, 2, 3, ... The order of the
    numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Fill()
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = i;
  }


  //! Fills the non-zero entries with a given value.
  /*!
    \param x the value to set the non-zero entries to.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x;
  }


  //! Fills the non-zero entries randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = rand();
  }


  //! Fills the matrix with random elements.
  /*! The matrix is cleared and then filled with \a n random elements. Both
    the position of the elements and their values are randomly generated. On
    exit, the matrix may not have \a n non-zero elements: it is possible that
    the randomly-generated positions of two elements are the same.
    \param[in] Nelement the number of random elements to be inserted in the
    matrix.
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::FillRand(int Nelement)
  {
    if (this->m_ == 0 || this->n_ == 0)
      return;

    Vector<int> i(Nelement), j(Nelement);
    Vector<T> value(Nelement);

    set<pair<int, int> > skeleton;
    set<pair<int, int> >::iterator it;

    srand(time(NULL));

    // generation of triplet (i, j, value)
    while (static_cast<int>(skeleton.size()) != Nelement)
      skeleton.insert(make_pair(rand() % this->m_, rand() % this->n_));

    int l = 0;
    for (it = skeleton.begin(); it != skeleton.end(); it++)
      {
	i(l) = it->first;
	j(l) = it->second;
	value(l) = double(rand());
	l++;
      }

    // then conversion to current sparse matrix
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);

    ConvertMatrix_from_Coordinates(i, j, value, leaf_class, 0);
  }


  //! Fills the matrix with one value inserted at random positions.
  /*! The matrix is cleared and then filled with \a n random elements. Only
    the position of the elements is randomly generated. Their value will
    always be \a x. On exit, the matrix may not have \a n non-zero elements:
    it is possible that the randomly-generated positions of two elements are
    the same.
    \param[in] Nelement the number of random elements to be inserted in the
    matrix.
    \param[in] x the value to be inserted.
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::FillRand(int Nelement, const T& x)
  {
    if (this->m_ == 0 || this->n_ == 0)
      return;

    Vector<int> i(Nelement), j(Nelement);
    Vector<T> value(Nelement);
    value.Fill(x);

    srand(time(NULL));
    for (int l = 0; l < Nelement; l++)
      {
        i(l) = rand() % this->m_;
        j(l) = rand() % this->n_;
      }

    // then conversion to current sparse matrix
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);

    ConvertMatrix_from_Coordinates(i, j, value, leaf_class, 0);
  }


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Sparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Stores the matrix in an output stream in binary format.
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Sparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->nz_)),
		     sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->ptr_),
		     sizeof(int)*(Storage::GetFirst(this->m_, this->n_)+1));
    FileStream.write(reinterpret_cast<char*>(this->ind_),
		     sizeof(int)*this->nz_);
    FileStream.write(reinterpret_cast<char*>(this->data_),
		     sizeof(T)*this->nz_);
  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in ascii format.
    The entries are written in coordinate format (row column value)
    1-index convention is used
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  WriteText(string FileName) const
  {
    ofstream FileStream; FileStream.precision(14);
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Sparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Stores the matrix in a file in ascii format.
    The entries are written in coordinate format (row column value)
    1-index convention is used
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Sparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    // conversion in coordinate format (1-index convention)
    IVect IndRow, IndCol; Vector<T> Value;
    const Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this);

    ConvertMatrix_to_Coordinates(leaf_class, IndRow, IndCol,
				 Value, 1, true);

    for (int i = 0; i < IndRow.GetM(); i++)
      FileStream << IndRow(i) << " " << IndCol(i) << " " << Value(i) << '\n';

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Sparse::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    \param FileStream input stream
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Sparse::Read(istream& FileStream)",
		    "Stream is not ready.");
#endif

    int m, n, nz;
    FileStream.read(reinterpret_cast<char*>(&m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&n), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&nz), sizeof(int));

    Reallocate(m, n, nz);

    FileStream.read(reinterpret_cast<char*>(ptr_),
                    sizeof(int)*(Storage::GetFirst(m, n)+1));
    FileStream.read(reinterpret_cast<char*>(ind_), sizeof(int)*nz);
    FileStream.read(reinterpret_cast<char*>(this->data_), sizeof(T)*nz);

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_Sparse::Read(istream& FileStream)",
                    string("Input operation failed.")
		    + string(" The input file may have been removed")
		    + " or may not contain enough data.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads the matrix from a file in text format.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Sparse::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix from a stream in text format.
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Sparse<T, Prop, Storage, Allocator>::
  ReadText(istream& FileStream)
  {
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);

    T zero; int index = 1;
    ReadCoordinateMatrix(leaf_class, FileStream, zero, index);
  }


  ///////////////////////
  // MATRIX<COLSPARSE> //
  ///////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSparse, Allocator>::Matrix():
    Matrix_Sparse<T, Prop, ColSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSparse, Allocator>::Matrix(int i, int j):
    Matrix_Sparse<T, Prop, ColSparse, Allocator>(i, j, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSparse, Allocator>::Matrix(int i, int j, int nz):
    Matrix_Sparse<T, Prop, ColSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr row start indices.
    \param ind column indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, ColSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<int, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_Sparse<T, Prop, ColSparse, Allocator>(i, j, values, ptr, ind)
  {
  }



  ///////////////////////
  // MATRIX<ROWSPARSE> //
  ///////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSparse, Allocator>::Matrix():
    Matrix_Sparse<T, Prop, RowSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSparse, Allocator>::Matrix(int i, int j):
    Matrix_Sparse<T, Prop, RowSparse, Allocator>(i, j, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with nz non-zero elements.
    \param i number of rows.
    \param j number of columns.
    \param nz number of non-zero elements.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSparse, Allocator>::Matrix(int i, int j, int nz):
    Matrix_Sparse<T, Prop, RowSparse, Allocator>(i, j, nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'values' (values), 'ptr' (pointers) and 'ind' (indices).
    Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param values values of non-zero entries.
    \param ptr column start indices.
    \param ind row indices.
    \warning Input vectors 'values', 'ptr' and 'ind' are empty on exit.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, RowSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& values,
	 Vector<int, Storage1, Allocator1>& ptr,
	 Vector<int, Storage2, Allocator2>& ind):
    Matrix_Sparse<T, Prop, RowSparse, Allocator>(i, j, values, ptr, ind)
  {
  }


  //////////////////////////
  // ReadCoordinateMatrix //
  //////////////////////////


  //! Reading of matrix in coordinate format
  /*!
    \param FileStream stream where the matrix is read
    \param row_numbers row indices
    \param col_numbers column indices
    \param values
  */
  template<class Tint, class AllocInt, class T, class Allocator>
  void ReadCoordinateMatrix(istream& FileStream,
                            Vector<Tint, VectFull, AllocInt>& row_numbers,
                            Vector<Tint, VectFull, AllocInt>& col_numbers,
                            Vector<T, VectFull, Allocator>& values)
  {
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("ReadCoordinateMatrix", "Stream is not ready.");
#endif

    T entry; int row = 0, col = 0;
    int nb_elt = 0;
    SetComplexZero(entry);

    while (!FileStream.eof())
      {
	// new entry is read (1-index)
	FileStream >> row >> col >> entry;

	if (FileStream.fail())
	  break;
	else
	  {
#ifdef SELDON_CHECK_IO
	    if (row < 1)
	      throw IOError(string("ReadCoordinateMatrix"),
			    string("Error : Row number should be greater ")
			    + "than 0 but is equal to " + to_str(row));

	    if (col < 1)
	      throw IOError(string("ReadCoordinateMatrix"),
			    string("Error : Column number should be greater")
			    + " than 0 but is equal to " + to_str(col));
#endif

	    nb_elt++;

	    // inserting new element
	    if (nb_elt > values.GetM())
	      {
		values.Resize(2*nb_elt);
		row_numbers.Resize(2*nb_elt);
		col_numbers.Resize(2*nb_elt);
	      }

	    values(nb_elt-1) = entry;
	    row_numbers(nb_elt-1) = row;
	    col_numbers(nb_elt-1) = col;
	  }
      }

    if (nb_elt > 0)
      {
	row_numbers.Resize(nb_elt);
	col_numbers.Resize(nb_elt);
	values.Resize(nb_elt);
      }
  }


  //! Reading of matrix in coordinate format
  /*!
    \param A output matrix
    \param FileStream stream where the matrix is read
    \param zero type of value to read (double or complex<double>)
    \param index starting index (usually 0 or 1)
    \param nnz number of non-zero entries
    If nnz is equal to -1, we consider that the number of non-zero entries
    is unknown and is deduced from the number of lines present in the stream
  */
  template<class Matrix1, class T>
  void ReadCoordinateMatrix(Matrix1& A, istream& FileStream, T& zero,
                            int index, int nnz)
  {
    // previous elements are removed
    A.Clear();

    Vector<int> row_numbers, col_numbers;
    Vector<T> values;
    if (nnz >= 0)
      {
        values.Reallocate(nnz);
        row_numbers.Reallocate(nnz);
        col_numbers.Reallocate(nnz);
      }

    ReadCoordinateMatrix(FileStream, row_numbers, col_numbers, values);

    if (row_numbers.GetM() > 0)
      ConvertMatrix_from_Coordinates(row_numbers, col_numbers, values,
                                     A, index);

  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_SPARSE_CXX
#endif
