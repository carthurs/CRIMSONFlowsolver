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


#ifndef SELDON_FILE_MATRIX_HERMITIAN_CXX

#include "Matrix_Hermitian.hxx"

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
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>::Matrix_Hermitian():
    Matrix_Base<T, Allocator>()
  {
    me_ = NULL;
  }


  //! Main constructor.
  /*! Builds a i x j full matrix.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Matrix_Hermitian(int i, int j):
    Matrix_Base<T, Allocator>(i, i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>(calloc(i, sizeof(pointer)));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
      }
    if (me_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
      }
    if ((me_ == NULL) && (i != 0))
      throw NoMemory("Matrix_Hermitian::Matrix_Hermitian(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(i)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(i)
		     + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate(i * i, this);

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
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(me_);
	me_ = NULL;
      }
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Matrix_Hermitian::Matrix_Hermitian(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(i)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(i)
		     + " elements).");
#endif

    pointer ptr = this->data_;
    int lgth = i;
    for (int k = 0; k < i; k++, ptr += lgth)
      me_[k] = ptr;
  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Matrix_Hermitian(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A):
    Matrix_Base<T, Allocator>()
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
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>::~Matrix_Hermitian()
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
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_Hermitian();
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
  int Matrix_Hermitian<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return this->m_ * this->n_;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x i matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {

    if (i != this->m_)
      {
	this->m_ = i;
	this->n_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    me_ = reinterpret_cast<pointer*>( realloc(me_,
						      i * sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (me_ == NULL)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    this->data_ = NULL;
	  }
	if (me_ == NULL && i != 0)
	  throw NoMemory("Matrix_Hermitian::Reallocate(int, int)",
			 string("Unable to reallocate memory")
			 + " for a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(i)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(i)
			 + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->allocator_.reallocate(this->data_,
								    i * i,
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
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    free(me_);
	    me_ = NULL;
	  }
	if (this->data_ == NULL && i != 0)
	  throw NoMemory("Matrix_Hermitian::Reallocate(int, int)",
			 string("Unable to reallocate memory")
			 + " for a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(i)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(i)
			 + " elements).");
#endif

	pointer ptr = this->data_;
	int lgth = Storage::GetSecond(i, i);
	for (int k = 0; k < Storage::GetFirst(i, i); k++, ptr += lgth)
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
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::SetData(int i, int j,
	    typename Matrix_Hermitian<T, Prop, Storage, Allocator>
	    ::pointer data)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>(calloc(i, sizeof(pointer)));

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
    int lgth = i;
    for (int k = 0; k < i; k++, ptr += lgth)
      me_[k] = ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
    \warning Memory is not released except for me_.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>::Nullify()
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
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Resize(int i, int j)
  {

    // Storing the old values of the matrix.
    int iold = Storage::GetFirst(this->m_, this->n_);
    int jold = Storage::GetSecond(this->m_, this->n_);
    Vector<value_type, VectFull, Allocator> xold(this->GetDataSize());
    for (int k = 0; k < this->GetDataSize(); k++)
      xold(k) = this->data_[k];

    // Reallocation.
    int inew = Storage::GetFirst(i, j);
    int jnew = Storage::GetSecond(i, j);
    this->Reallocate(i, j);

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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::value_type
  Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::operator() const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::operator() const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    if (i > j)
      return conj(me_[Storage::GetSecond(i, j)][Storage::GetFirst(i, j)]);
    else
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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Returns the element (i, j).
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Get(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Returns the element (i, j).
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_Hermitian::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (i > j)
      throw WrongRow("Matrix_Hermitian::Val(int, int)",
		     string("Attempted to access to element (")
		     + to_str(i) + ", " + to_str(j)
		     + ") but row index should not be strictly"
		     + " greater than column index.");
#endif

    return me_[Storage::GetFirst(i, j)][Storage::GetSecond(i, j)];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>::reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_Hermitian::operator[] (int)",
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
  inline typename Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_Hermitian::operator[] (int) const",
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
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    if (i > j)
      this->Val(j, i) = conj(x);
    else
      this->Val(i, j) = x;
  }


  //! Duplicates a matrix (assignement operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_Hermitian<T, Prop, Storage, Allocator>&
  Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::operator=(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A)
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
  inline void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_Hermitian<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    this->allocator_.memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::Zero()
  {
    this->allocator_.memoryset(this->data_, char(0),
			       this->GetDataSize() * sizeof(value_type));
  }


  //! Sets the matrix to the identity.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, discard this method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::SetIdentity()
  {
    this->Fill(T(0));

    T one;
    SetComplexOne(one);
    for (int i = 0; i < min(this->m_, this->n_); i++)
      this->Val(i, i) = one;
  }


  //! Fills the matrix with 0, 1, 2, ...
  /*!
    On exit, the matrix is filled with 0, 1, 2, 3, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::Fill()
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = i;
  }


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
    \note If the imaginary part of x is non-null, the upper part will
    contain x, whereas lower part will contain conj(x).
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x;
  }


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
    \note If the imaginary part of x is non-null, the upper part will
    contain x, whereas lower part will contain conj(x).
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_Hermitian<T, Prop, Storage, Allocator>&
  Matrix_Hermitian<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the matrix randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::FillRand()
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
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::Print() const
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
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Print(int a, int b, int m, int n) const
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
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::Print(int l) const
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
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Hermitian::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Hermitian::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * this->n_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_Hermitian::Write(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + " or there is no space left on device.");
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
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Hermitian::WriteText(string FileName)",
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
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Hermitian::WriteText(ofstream& FileStream)",
                    "Stream is not ready.");
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
      throw IOError("Matrix_Hermitian::WriteText(ofstream& FileStream)",
                    string("Output operation failed.")
		    + string(" The output file may have been removed")
		    + " or there is no space left on device.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Hermitian::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Hermitian::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    int new_m, new_n;
    FileStream.read(reinterpret_cast<char*>(&new_m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&new_n), sizeof(int));
    this->Reallocate(new_m, new_n);

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    new_m * new_n * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_Hermitian::Read(ifstream& FileStream)",
                    string("Input operation failed.")
		    + string(" The input file may have been removed")
		    + " or may not contain enough data.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in text format in a file.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Hermitian<T, Prop, Storage, Allocator>::ReadText(string FileName)
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
  void Matrix_Hermitian<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream)
  {
    // clears previous matrix
    Clear();

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::ReadText(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    // we read first line
    string line;
    getline(FileStream, line);

    if (FileStream.fail())
      {
	// empty file ?
	return;
      }

    // converting first line into a vector
    istringstream line_stream(line);
    Vector<T> first_row;
    first_row.ReadText(line_stream);

    // and now the other rows
    Vector<T> other_rows;
    other_rows.ReadText(FileStream);

    // number of rows and columns
    int n = first_row.GetM();
    int m = 1 + other_rows.GetM()/n;

#ifdef SELDON_CHECK_IO
    // Checking number of elements
    if (other_rows.GetM() != (m-1)*n)
      throw IOError("Matrix_Pointers::ReadText(ifstream& FileStream)",
                    "The file should contain same number of columns.");
#endif

    this->Reallocate(m,n);
    // filling matrix
    for (int j = 0; j < n; j++)
      this->Val(0, j) = first_row(j);

    int nb = 0;
    for (int i = 1; i < m; i++)
      {
	for (int j = 0; j < i; j++)
	  nb++;

	for (int j = i; j < n; j++)
	  this->Val(i, j) = other_rows(nb++);
      }
  }



  ////////////////////
  // MATRIX<COLHERM> //
  ////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColHerm, Allocator>::Matrix():
    Matrix_Hermitian<T, Prop, ColHerm, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j full column-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColHerm, Allocator>::Matrix(int i, int j):
    Matrix_Hermitian<T, Prop, ColHerm, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator= (const Matrix<T, Prop,
                                                  ColHerm,
                                                  Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
    \warning The imaginary part of x should be null to keep an Hermitian
    matrix.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, ColHerm, Allocator>&
  Matrix<T, Prop, ColHerm, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  ////////////////////
  // MATRIX<ROWHERM> //
  ////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowHerm, Allocator>::Matrix():
    Matrix_Hermitian<T, Prop, RowHerm, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i by j full row-major matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowHerm, Allocator>::Matrix(int i, int j):
    Matrix_Hermitian<T, Prop, RowHerm, Allocator>(i, j)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Fills a matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator= (const Matrix<T, Prop,
                                                       RowHerm,
                                                       Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }


  //! Multiplies the matrix by a given value.
  /*!
    \param x multiplication coefficient
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, RowHerm, Allocator>&
  Matrix<T, Prop, RowHerm, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_HERMITIAN_CXX
#endif
