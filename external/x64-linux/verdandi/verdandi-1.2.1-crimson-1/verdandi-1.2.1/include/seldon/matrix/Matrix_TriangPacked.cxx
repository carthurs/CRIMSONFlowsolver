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


#ifndef SELDON_FILE_MATRIX_TRIANGPACKED_CXX

#include "Matrix_TriangPacked.hxx"

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
  inline Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Matrix_TriangPacked(): Matrix_Base<T, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j triangular matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Matrix_TriangPacked(int i, int j): Matrix_Base<T, Allocator>(i, i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->allocator_.allocate((i * (i + 1)) / 2, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
	return;
      }
    if (this->data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	return;
      }
#endif

  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Matrix_TriangPacked(const Matrix_TriangPacked<T, Prop,
			Storage, Allocator>& A)
    : Matrix_Base<T, Allocator>()
  {
    this->m_ = 0;
    this->n_ = 0;
    this->data_ = NULL;

    this->Copy(A);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::~Matrix_TriangPacked()
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    this->allocator_.deallocate(this->data_,
					(this->m_ * (this->m_ + 1)) / 2);
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
      }
#endif

  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_TriangPacked();
    this->m_ = 0;
    this->n_ = 0;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_TriangPacked<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return (this->m_ * (this->m_ + 1)) / 2;
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
  inline void Matrix_TriangPacked<T, Prop, Storage, Allocator>
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

	    this->data_ =
	      reinterpret_cast<pointer>(this->allocator_.reallocate(this->data_,
								    (i*(i+1))/2,
								    this) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    this->data_ = NULL;
	    throw NoMemory("Matrix_TriangPacked::Reallocate(int, int)",
			   "Unable to reallocate memory for data_.");
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    throw NoMemory("Matrix_TriangPacked::Reallocate(int, int)",
			   "Unable to reallocate memory for data_.");
	  }
#endif

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
  inline void Matrix_TriangPacked<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  typename Matrix_TriangPacked<T, Prop, Storage, Allocator>
	  ::pointer data)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

    this->data_ = data;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Nullify()
  {
    this->data_ = NULL;
    this->m_ = 0;
    this->n_ = 0;
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
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>::value_type
  Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_TriangPacked::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_TriangPacked::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    if (Storage::UpLo())
      if (i > j)
	return 0;
      else
	return this->data_[Storage::GetFirst(i * this->n_
					     - (i * (i + 1)) / 2 + j,
					     (j * (j + 1)) / 2 + i)];
    else
      if (i < j)
	return 0;
      else
	return this->data_[Storage::GetFirst((i * (i + 1)) / 2 + j,
					     j * this->m_
					     - (j * (j + 1)) / 2 + i)];
  }


  //! Direct access method.
  /*!
    This method allows access to elements stored in memory, i.e. elements
    from the upper part. i <= j must be satisfied.
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>::reference
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_TriangPacked::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_TriangPacked::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (Storage::UpLo())
      {
	if (i > j)
	  throw WrongRow("Matrix_TriangPacked::Val(int, int)",
			 string("Attempted to access to element (")
			 + to_str(i) + ", " + to_str(j) + string(") but row")
			 + string(" index should not be strictly more")
			 + " than column index (upper triangular matrix).");
	return this->data_[Storage::GetFirst(i * this->n_
					     - (i * (i + 1)) / 2 + j,
					     (j * (j + 1)) / 2 + i)];
      }
    else
      {
	if (j > i)
	  throw WrongCol("Matrix_TriangPacked::Val(int, int)",
			 string("Attempted to access to element (")
			 + to_str(i) + ", " + to_str(j) + string(") but")
			 + string(" column index should not be strictly more")
			 + " than row index (lower triangular matrix).");
	return this->data_[Storage::GetFirst((i * (i + 1)) / 2 + j,
					     j * this->m_
					     - (j * (j + 1)) / 2 + i)];
      }
#endif

    if (Storage::UpLo())
      return this->data_[Storage::GetFirst(i * this->n_
					   - (i * (i + 1)) / 2 + j,
					   (j * (j + 1)) / 2 + i)];
    else
      return this->data_[Storage::GetFirst((i * (i + 1)) / 2 + j,
					   j * this->m_
					   - (j * (j + 1)) / 2 + i)];
  }


  //! Direct access method.
  /*!
    This method allows access to elements stored in memory, i.e. elements
    from the upper part. i <= j must be satisfied.
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_TriangPacked::Val(int, int) const",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_TriangPacked::Val(int, int)",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
    if (Storage::UpLo())
      {
	if (i > j)
	  throw WrongRow("Matrix_TriangPacked::Val(int, int) const",
			 string("Attempted to access to element (")
			 + to_str(i) + ", " + to_str(j) + string(") but row")
			 + string(" index should not be strictly more")
			 + " than column index (upper triangular matrix).");
	return this->data_[Storage::GetFirst(i * this->n_
					     - (i * (i + 1)) / 2 + j,
					     (j * (j + 1)) / 2 + i)];
      }
    else
      {
	if (j > i)
	  throw WrongCol("Matrix_TriangPacked::Val(int, int) const",
			 string("Attempted to access to element (")
			 + to_str(i) + ", " + to_str(j) + string(") but")
			 + string(" column index should not be strictly more")
			 + " than row index (lower triangular matrix).");
	return this->data_[Storage::GetFirst((i * (i + 1)) / 2 + j,
					     j * this->m_
					     - (j * (j + 1)) / 2 + i)];
      }
#endif

    if (Storage::UpLo())
      return this->data_[Storage::GetFirst(i * this->n_
					   - (i * (i + 1)) / 2 + j,
					   (j * (j + 1)) / 2 + i)];
    else
      return this->data_[Storage::GetFirst((i * (i + 1)) / 2 + j,
					   j * this->m_
					   - (j * (j + 1)) / 2 + i)];
  }


  //! Returns the element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>::reference
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    return this->Val(i, j);
  }


  //! Returns the element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return this->Val(i, j);
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>::reference
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_TriangPacked::operator[] (int)",
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
  inline typename Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->GetDataSize())
      throw WrongIndex("Matrix_TriangPacked::operator[] (int) const",
		       string("Index should be in [0, ")
		       + to_str(this->GetDataSize()-1) + "], but is equal to "
		       + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_TriangPacked<T, Prop, Storage, Allocator>&
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::
  operator= (const Matrix_TriangPacked<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    this->Val(i, j) = x;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_TriangPacked<T, Prop, Storage, Allocator>::
  Copy(const Matrix_TriangPacked<T, Prop, Storage, Allocator>& A)
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Zero()
  {
    this->allocator_.memoryset(this->data_, char(0),
			       this->GetDataSize() * sizeof(value_type));
  }


  //! Sets the matrix to the identity.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::SetIdentity()
  {
    this->Fill(T(0));

    T one(1);
    bool storage_col = (Storage::GetFirst(1,0) == 0);
    if ((Storage::UpLo() && storage_col)||(!Storage::UpLo() && !storage_col))
      {
	int index(-1);
	for (int i = 0; i < min(this->m_, this->n_); i++)
	  {
	    index += i+1;
	    this->data_[index] = one;
	  }
      }
    else if (Storage::UpLo() && !storage_col)
      {
	int index(0);
	for (int i = 0; i < min(this->m_, this->n_); i++)
	  {
	    this->data_[index] = one;
	    index += this->m_-i;
	  }
      }
    else if (!Storage::UpLo() && storage_col)
      {
	int index(0);
	for (int i = 0; i < min(this->m_, this->n_); i++)
	  {
	    this->data_[index] = one;
	    index += this->n_-i;
	  }
      }
  }


  //! Fills the matrix with 0, 1, 2, ...
  /*!
    On exit, the matrix is filled with 0, 1, 2, 3, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Fill()
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = i;
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x;
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_TriangPacked<T, Prop, Storage, Allocator>&
  Matrix_TriangPacked<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the matrix randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::FillRand()
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Print() const
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Print(int a, int b, int m, int n) const
  {
    for (int i = a; i < min(this->m_, a + m); i++)
      {
	for (int j = b; j < min(this->n_, b + n); j++)
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Print(int l) const
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_TriangPacked::Write(string FileName)",
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("Matrix_TriangPacked::Write(ofstream& FileStream)",
                    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->GetDataSize() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_TriangPacked::Write(ofstream& FileStream)",
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_TriangPacked::WriteText(string FileName)",
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_TriangPacked::WriteText(ofstream& FileStream)",
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
      throw IOError("Matrix_TriangPacked::WriteText(ofstream& FileStream)",
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.good())
      throw IOError("Matrix_TriangPacked::Read(string FileName)",
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_TriangPacked::Read(ifstream& FileStream)",
                    "Stream is not ready.");
#endif

    int new_m, new_n;
    FileStream.read(reinterpret_cast<char*>(&new_m), sizeof(int));
    FileStream.read(reinterpret_cast<char*>(&new_n), sizeof(int));
    this->Reallocate(new_m, new_n);

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    this->GetDataSize() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_TriangPacked::Read(ifstream& FileStream)",
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>::ReadText(string FileName)
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
  void Matrix_TriangPacked<T, Prop, Storage, Allocator>
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
    if (Storage::UpLo())
      for (int j = 0; j < n; j++)
	this->Val(0, j) = first_row(j);
    else
      this->Val(0, 0) = first_row(0);

    int nb = 0;
    if (Storage::UpLo())
      for (int i = 1; i < m; i++)
	{
	  for (int j = 0; j < i; j++)
	    nb++;

	  for (int j = i; j < n; j++)
	    this->Val(i, j) = other_rows(nb++);
	}
    else
      for (int i = 1; i < m; i++)
	{
	  for (int j = 0; j <= i; j++)
	    this->Val(i, j) = other_rows(nb++);

	  for (int j = i+1; j < n; j++)
	    nb++;
	}
  }



  ///////////////////////////////
  // MATRIX<COLUPTRIANGPACKED> //
  ///////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColUpTriangPacked, Allocator>::Matrix():
    Matrix_TriangPacked<T, Prop, ColUpTriangPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j column-major upper triangular matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColUpTriangPacked, Allocator>::Matrix(int i, int j):
    Matrix_TriangPacked<T, Prop, ColUpTriangPacked, Allocator>(i, i)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ColUpTriangPacked, Allocator>
  ::Resize(int i, int j)
  {

    // Storing the old values of the matrix.
    int nold = this->GetDataSize();
    Vector<T, VectFull, Allocator> xold(nold);
    for (int k = 0; k < nold; k++)
      xold(k) = this->data_[k];

    // Reallocation.
    this->Reallocate(i, j);

    // Filling the matrix with its old values.
    int nmin = min(nold, this->GetDataSize());
    for (int k = 0; k < nmin; k++)
      this->data_[k] = xold(k);
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, ColUpTriangPacked, Allocator>&
  Matrix<T, Prop, ColUpTriangPacked, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, ColUpTriangPacked, Allocator>&
  Matrix<T, Prop, ColUpTriangPacked, Allocator>::operator=
  (const Matrix<T, Prop, ColUpTriangPacked, Allocator>& A)
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
  Matrix<T, Prop, ColUpTriangPacked, Allocator>&
  Matrix<T, Prop, ColUpTriangPacked, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  ///////////////////////////////
  // MATRIX<COLLOTRIANGPACKED> //
  ///////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColLoTriangPacked, Allocator>::Matrix():
    Matrix_TriangPacked<T, Prop, ColLoTriangPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j column-major lower triangular matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColLoTriangPacked, Allocator>::Matrix(int i, int j):
    Matrix_TriangPacked<T, Prop, ColLoTriangPacked, Allocator>(i, i)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ColLoTriangPacked, Allocator>
  ::Resize(int i, int j)
  {

    // Storing the old values of the matrix.
    int nold = this->GetDataSize(), iold = this->m_;
    Vector<T, VectFull, Allocator> xold(nold);
    for (int k = 0; k < nold; k++)
      xold(k) = this->data_[k];

    // Reallocation.
    this->Reallocate(i, j);

    // Filling the matrix with its old values.
    int imin = min(iold, i);
    nold = 0;
    int n = 0;
    for (int k = 0; k < imin; k++)
      {
	for (int l = k; l < imin; l++)
	  this->data_[n+l-k] = xold(nold+l-k);

	n += i - k;
	nold += iold - k;
      }
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, ColLoTriangPacked, Allocator>&
  Matrix<T, Prop, ColLoTriangPacked, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, ColLoTriangPacked, Allocator>&
  Matrix<T, Prop, ColLoTriangPacked, Allocator>::operator=
  (const Matrix<T, Prop, ColLoTriangPacked, Allocator>& A)
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
  Matrix<T, Prop, ColLoTriangPacked, Allocator>&
  Matrix<T, Prop, ColLoTriangPacked, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  ///////////////////////////////
  // MATRIX<ROWUPTRIANGPACKED> //
  ///////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowUpTriangPacked, Allocator>::Matrix():
    Matrix_TriangPacked<T, Prop, RowUpTriangPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j row-major upper triangular matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowUpTriangPacked, Allocator>::Matrix(int i, int j):
    Matrix_TriangPacked<T, Prop, RowUpTriangPacked, Allocator>(i, i)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, RowUpTriangPacked, Allocator>
  ::Resize(int i, int j)
  {

    // Storing the old values of the matrix.
    int nold = this->GetDataSize(), iold = this->m_;
    Vector<T, VectFull, Allocator> xold(nold);
    for (int k = 0; k < nold; k++)
      xold(k) = this->data_[k];

    // Reallocation.
    this->Reallocate(i, j);

    // Filling the matrix with its old values.
    int imin = min(iold, i);
    nold = 0;
    int n = 0;
    for (int k = 0; k < imin; k++)
      {
	for (int l = k; l < imin; l++)
	  this->data_[n+l-k] = xold(nold+l-k);

	n += i - k;
	nold += iold - k;
      }
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, RowUpTriangPacked, Allocator>&
  Matrix<T, Prop, RowUpTriangPacked, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, RowUpTriangPacked, Allocator>&
  Matrix<T, Prop, RowUpTriangPacked, Allocator>::operator=
  (const Matrix<T, Prop, RowUpTriangPacked, Allocator>& A)
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
  Matrix<T, Prop, RowUpTriangPacked, Allocator>&
  Matrix<T, Prop, RowUpTriangPacked, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }



  ///////////////////////////////
  // MATRIX<ROWLOTRIANGPACKED> //
  ///////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowLoTriangPacked, Allocator>::Matrix():
    Matrix_TriangPacked<T, Prop, RowLoTriangPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j row-major lower triangular matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowLoTriangPacked, Allocator>::Matrix(int i, int j):
    Matrix_TriangPacked<T, Prop, RowLoTriangPacked, Allocator>(i, i)
  {
  }


  /*****************
   * OTHER METHODS *
   *****************/


  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, RowLoTriangPacked, Allocator>
  ::Resize(int i, int j)
  {

    // Storing the old values of the matrix.
    int nold = this->GetDataSize();
    Vector<T, VectFull, Allocator> xold(nold);
    for (int k = 0; k < nold; k++)
      xold(k) = this->data_[k];

    // Reallocation.
    this->Reallocate(i, j);

    // Filling the matrix with its old values.
    int nmin = min(nold, this->GetDataSize());
    for (int k = 0; k < nmin; k++)
      this->data_[k] = xold(k);
  }


  //! Fills the matrix with a given value.
  /*!
    \param x value to fill the matrix with.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  Matrix<T, Prop, RowLoTriangPacked, Allocator>&
  Matrix<T, Prop, RowLoTriangPacked, Allocator>::operator= (const T0& x)
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
  inline Matrix<T, Prop, RowLoTriangPacked, Allocator>&
  Matrix<T, Prop, RowLoTriangPacked, Allocator>::operator=
  (const Matrix<T, Prop, RowLoTriangPacked, Allocator>& A)
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
  Matrix<T, Prop, RowLoTriangPacked, Allocator>&
  Matrix<T, Prop, RowLoTriangPacked, Allocator>::operator*= (const T0& x)
  {
    for (int i = 0; i < this->GetDataSize();i++)
      this->data_[i] *= x;

    return *this;
  }

} // namespace Seldon.

#define SELDON_FILE_MATRIX_TRIANGPACKED_CXX
#endif
