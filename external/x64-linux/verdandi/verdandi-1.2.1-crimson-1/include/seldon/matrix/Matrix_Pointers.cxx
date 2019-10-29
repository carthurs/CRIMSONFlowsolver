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


#ifndef SELDON_FILE_MATRIX_POINTERS_CXX

#include "Matrix_Pointers.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>::Matrix_Pointers():
    Matrix_Base<T, Allocator>()
  {
    me_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j full matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Matrix_Pointers(int i, int j): Matrix_Base<T, Allocator>(i, j)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>( calloc(Storage::GetFirst(i, j),
						 sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
      }
    if (me_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Pointers::Matrix_Pointers(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate(i * j, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(me_);
	me_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Pointers::Matrix_Pointers(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " elements).");
#endif

    pointer ptr = this->data_;
    int lgth = Storage::GetSecond(i, j);
    for (int k = 0; k < Storage::GetFirst(i, j); k++, ptr += lgth)
      me_[k] = ptr;

  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Matrix_Pointers(const Matrix_Pointers<T, Prop, Storage, Allocator>& A):
    Matrix_Base<T, Allocator>(A)
  {
    this->m_ = 0;
    this->n_ = 0;
    this->data_ = NULL;
    this->me_ = NULL;

    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>::~Matrix_Pointers()
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_, this->m_ * this->n_);
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (me_ != NULL)
	  {
	    free(me_);
	    me_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
      }
#endif

  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_Pointers();
    this->m_ = 0;
    this->n_ = 0;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of rows multiplied by the number of columns
    because the matrix is full.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_Pointers<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return this->m_ * this->n_;
  }


  //! Returns the pointer 'me_'.
  /*! Returns the pointer 'me_' that defines an array pointing to the first
    row or column elements, so that 'me_[1]' points to the first element of
    the second row or column.
    \return The pointer 'me_'.
  */
  template <class T, class Prop, class Storage, class Allocator>
  typename
  Matrix_Pointers<T, Prop, Storage, Allocator>::pointer*
  Matrix_Pointers<T, Prop, Storage, Allocator>::GetMe() const
  {
    return me_;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {

    if (i != this->m_ || j != this->n_)
      {
	this->m_ = i;
	this->n_ = j;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    me_ = reinterpret_cast<pointer*>( realloc(me_,
						      Storage::GetFirst(i, j)
						      * sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (me_ == NULL && i != 0 && j != 0)
	  throw NoMemory("Matrix_Pointers::Reallocate(int, int)",
			 string("Unable to reallocate memory for")
			 + " a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->allocator_
					.reallocate(this->data_, i * j,
						    this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    free(me_);
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (this->data_ == NULL && i != 0 && j != 0)
	  throw NoMemory("Matrix_Pointers::Reallocate(int, int)",
			 string("Unable to reallocate memory")
			 + " for a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " elements).");
#endif

	pointer ptr = this->data_;
	int lgth = Storage::GetSecond(i, j);
	for (int k = 0; k < Storage::GetFirst(i, j); k++, ptr += lgth)
	  me_[k] = ptr;
      }
  }


  //! Changes the size of the matrix and sets its data array
  //! (low level method).
  /*!
    The matrix is first cleared (memory is freed). The matrix is then resized
    to a i x j matrix, and the data array of the matrix is set to 'data'.
    'data' elements are not duplicated: the new data array of the matrix is
    the 'data' array. It is useful to create a matrix from pre-existing data.
    \param i new number of rows.
    \param j new number of columns.
    \param data new array storing elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The matrix
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::SetData(int i, int j,
	    typename Matrix_Pointers<T, Prop, Storage, Allocator>
	    ::pointer data)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = j;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>( calloc(Storage::GetFirst(i, j),
						 sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
	return;
      }
    if (me_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
	return;
      }
#endif

    this->data_ = data;

    pointer ptr = this->data_;
    int lgth = Storage::GetSecond(i, j);
    for (int k = 0; k < Storage::GetFirst(i, j); k++, ptr += lgth)
      me_[k] = ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
    \warning Memory is not released except for me_.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (me_ != NULL)
	  {
	    free(me_);
	    me_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
      }
#endif

    this->data_ = NULL;
  }


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Resize(int i, int j)
  {

    if (i == this->m_ && j == this->n_)
      return;

    // Storing the old values of the matrix.
    int iold = Storage::GetFirst(this->m_, this->n_);
    int jold = Storage::GetSecond(this->m_, this->n_);
    Vector<value_type, VectFull, Allocator> xold(this->GetDataSize());
    for (int k = 0; k < this->GetDataSize(); k++)
      xold(k) = this->data_[k];

    // Reallocation.
    int inew = Storage::GetFirst(i, j);
    int jnew = Storage::GetSecond(i, j);
    this->Reallocate(i,j);

    // Filling the matrix with its old values.
    int imin = min(iold, inew), jmin = min(jold, jnew);
    for (int k = 0; k < imin; k++)
      for (int l = 0; l < jmin; l++)
	this->data_[k*jnew+l] = xold(l+jold*k);
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
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator() (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Pointers::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Pointers::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference Matrix_Pointers<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Pointers::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Pointers::operator()",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Pointers::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Pointers::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    return Val(i, j);
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Pointers::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Pointers::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return Val(i, j);
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>::reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_Pointers::operator[] (int)",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Pointers<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_Pointers::operator[] (int) const",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& val)
  {
    this->Val(i, j) = val;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Pointers<T, Prop, Storage, Allocator>&
  Matrix_Pointers<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_Pointers<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_Pointers<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    this->allocator_.memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Returns the leading dimension.
  /*!
    \return The leading dimension.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_Pointers<T, Prop, Storage, Allocator>::GetLD() const
  {
    return Storage::GetSecond(this->m_, this->n_);
  }


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Zero()
  {
    this->allocator_.memoryset(this->data_, char(0),
			       this->GetDataSize() * sizeof(value_type));
  }


  //! Sets the current matrix to the identity.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::SetIdentity()
  {
    Fill(T(0));

    T one(1);
    for (int i = 0; i < min(this->m_, this->n_); i++)
      (*this)(i,i) = one;
  }


  //! Fills the matrix the matrix with 0, 1, 2, ...
  /*!
    On exit, the matrix is filled with 0, 1, 2, 3, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Fill()
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = i;
  }


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x;
  }


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_Pointers<T, Prop, Storage, Allocator>&
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the matrix randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = rand();
  }


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


  //! Displays a sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its upper-left corner (a, b)
    and its bottom-right corner (m, n). So, elements with indices
    in [a, m] x [b, n] are displayed on the standard output,
    in text format. Each row is displayed on a single line and
    elements of a row are delimited by tabulations.
    \param a row index of the upper-left corner.
    \param b column index of the upper-left corner.
    \param m row index of the bottom-right corner.
    \param n column index of the bottom-right corner.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Print(int a, int b,
							   int m, int n) const
  {
    for (int i = a; i < min(this->m_, a+m); i++)
      {
	for (int j = b; j < min(this->n_, b+n); j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


  //! Displays a square sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its bottom-right corner (l, l).
    So, elements with indices in [0, 0] x [l, l] are displayed
    on the standard output, in text format. Each row is displayed
    on a single line and elements of a row are delimited
    by tabulations.
    \param l dimension of the square matrix to be displayed.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Print(int l) const
  {
    Print(0, 0, l, l);
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileName output file name.
    \param with_size if set to 'false', the dimensions of the matrix are not
    saved.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output stream.
    \param with_size if set to 'false', the dimensions of the matrix are not
    saved.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Write(ostream& FileStream)",
		    "The stream is not ready.");
#endif

    if (with_size)
      {
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
			 sizeof(int));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
			 sizeof(int));
      }

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * this->n_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Write(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    int i, j;
    for (i = 0; i < this->GetM(); i++)
      {
	for (j = 0; j < this->GetN(); j++)
	  FileStream << (*this)(i, j) << '\t';
	FileStream << endl;
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName input file name.
    \param with_size if set to 'false', the dimensions of the matrix are not
    available in the file. In this case, the current dimensions (M, N) of the
    matrix are unchanged, and MxN elements are read in the file.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream, with_size);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream input stream.
    \param with_size if set to 'false', the dimensions of the matrix are not
    available in the stream. In this case, the current dimensions (M, N) of
    the matrix are unchanged, and MxN elements are read in the stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(istream& FileStream)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int new_m, new_n;
        FileStream.read(reinterpret_cast<char*>(&new_m), sizeof(int));
        FileStream.read(reinterpret_cast<char*>(&new_n), sizeof(int));
        this->Reallocate(new_m, new_n);
      }

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    this->GetM() * this->GetN() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(istream& FileStream)",
                    "Input operation failed.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix from a file in text format.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in text format from an input stream.
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream)
  {
    // Clears the matrix.
    Clear();

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::ReadText(istream& FileStream)",
                    "The stream is not ready.");
#endif

    // Reads the first line.
    string line;
    getline(FileStream, line);
    if (FileStream.fail())
      // Is the file empty?
      return;

    // Converts the first line into a vector.
    istringstream line_stream(line);
    Vector<T> first_row;
    first_row.ReadText(line_stream);

    // Now reads all other rows, and puts them in a single vector.
    Vector<T> other_row;
    other_row.ReadText(FileStream);

    // Number of rows and columns.
    int n = first_row.GetM();
    int m = 1 + other_row.GetM() / n;

#ifdef SELDON_CHECK_IO
    // Checks that enough elements were read.
    if (other_row.GetM() != (m - 1) * n)
      throw IOError("Matrix_Pointers::ReadText(istream& FileStream)",
                    "Not all rows have the same number of columns.");
#endif

    this->Reallocate(m, n);
    // Fills the matrix.
    for (int j = 0; j < n; j++)
      this->Val(0, j) = first_row(j);

    int k = 0;
    for (int i = 1; i < m; i++)
      for (int j = 0; j < n; j++)
	this->Val(i, j) = other_row(k++);
  }


  //////////////////////
  // MATRIX<COLMAJOR> //
  //////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColMajor, Allocator>::Matrix():
    Matrix_Pointers<T, Prop, ColMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColMajor, Allocator>::Matrix(int i, int j):
    Matrix_Pointers<T, Prop, ColMajor, Allocator>(i, j)
  {
  }


  //! Copy constructor.
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColMajor, Allocator>
  ::Matrix(const Matrix<T, Prop, ColMajor, Allocator>& A):
    Matrix_Pointers<T, Prop, ColMajor, Allocator>(A)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, ColMajor, Allocator>&
  Matrix<T, Prop, ColMajor, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColMajor, Allocator>&
  Matrix<T, Prop, ColMajor, Allocator>
  ::operator= (const Matrix<T, Prop, ColMajor, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Multiplies the matrix by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Prop, class Allocator> template<class T0>
  inline Matrix<T, Prop, ColMajor, Allocator>&
  Matrix<T, Prop, ColMajor, Allocator>::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->m_ * this->n_; i++)
      this->data_[i] *= alpha;

    return *this;
  }


  //////////////////////
  // MATRIX<ROWMAJOR> //
  //////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds a empty matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowMajor, Allocator>::Matrix():
    Matrix_Pointers<T, Prop, RowMajor, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowMajor, Allocator>::Matrix(int i, int j):
    Matrix_Pointers<T, Prop, RowMajor, Allocator>(i, j)
  {
  }


  //! Copy constructor.
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowMajor, Allocator>
  ::Matrix(const Matrix<T, Prop, RowMajor, Allocator>& A):
    Matrix_Pointers<T, Prop, RowMajor, Allocator>(A)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, RowMajor, Allocator>&
  Matrix<T, Prop, RowMajor, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowMajor, Allocator>&
  Matrix<T, Prop, RowMajor, Allocator>
  ::operator= (const Matrix<T, Prop, RowMajor, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Multiplies the matrix by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Prop, class Allocator> template<class T0>
  inline Matrix<T, Prop, RowMajor, Allocator>&
  Matrix<T, Prop, RowMajor, Allocator>::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->m_*this->n_; i++)
      this->data_[i] *= alpha;

    return *this;
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_POINTERS_CXX
#endif
