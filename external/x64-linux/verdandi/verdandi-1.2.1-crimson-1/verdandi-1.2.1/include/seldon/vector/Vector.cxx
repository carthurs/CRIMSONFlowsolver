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


#ifndef SELDON_FILE_VECTOR_CXX

#include "Vector.hxx"

namespace Seldon
{


  ///////////////////
  // C_VECTOR_BASE //
  ///////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Nothing is allocated. Vector length is set to zero.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base()
  {
    m_ = 0;
    data_ = NULL;
  }


  //! Main constructor.
  /*!
    \param i length.
    \warning Nothing is allocated.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::Vector_Base(int i)
  {
    m_ = i;
    data_ = NULL;
  }


  //! Copy constructor.
  /*!
    \param A base vector to be copied.
    \warning Only the length is copied.
  */
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::
  Vector_Base(const Vector_Base<T, Allocator>& A)
  {
    m_ = A.GetM();
    data_ = NULL;
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  inline Vector_Base<T, Allocator>::~Vector_Base()
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->Clear();

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	m_ = 0;
	data_ = NULL;
      }
#endif

  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements.
  /*!
    \return The length of the vector.
  */
  template <class T, class Allocator>
  int Vector_Base<T, Allocator>::GetM() const
  {
    return m_;
  }


  //! Returns the number of elements.
  /*!
    \return The length of the vector.
  */
  template <class T, class Allocator>
  int Vector_Base<T, Allocator>::GetLength() const
  {
    return m_;
  }


  //! Returns the number of elements stored.
  /*!
    \return The length of the vector stored.
  */
  template <class T, class Allocator>
  int Vector_Base<T, Allocator>::GetSize() const
  {
    return m_;
  }


  //! Returns a pointer to data_ (stored data).
  /*!
    \return A pointer to the data_, i.e. the data array.
  */
  template <class T, class Allocator>
  typename Vector_Base<T, Allocator>::pointer
  Vector_Base<T, Allocator>::GetData() const
  {
    return data_;
  }


  //! Returns a const pointer to data_ (stored data).
  /*!
    \return A const pointer to the data_, i.e. the data array.
  */
  template <class T, class Allocator>
  typename Vector_Base<T, Allocator>::const_pointer
  Vector_Base<T, Allocator>::GetDataConst() const
  {
    return reinterpret_cast<typename Vector_Base<T,
      Allocator>::const_pointer>(data_);
  }


  //! Returns a pointer of type "void*" to the data array (data_).
  /*!
    \return A pointer of type "void*" to the data array.
  */
  template <class T, class Allocator>
  void* Vector_Base<T, Allocator>::GetDataVoid() const
  {
    return reinterpret_cast<void*>(data_);
  }


  //! Returns a pointer of type "const void*" to the data array (data_).
  /*!
    \return A pointer of type "const void*" to the data array.
  */
  template <class T, class Allocator>
  const void* Vector_Base<T, Allocator>::GetDataConstVoid() const
  {
    return reinterpret_cast<const void*>(data_);
  }


  ///////////////////////
  // VECTOR<VECT_FULL> //
  ///////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  Vector<T, VectFull, Allocator>::Vector():
    Vector_Base<T, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param i length of the vector.
  */
  template <class T, class Allocator>
  Vector<T, VectFull, Allocator>::Vector(int i):
    Vector_Base<T, Allocator>(i)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->vect_allocator_.allocate(i, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      this->m_ = 0;
    if (this->data_ == NULL && i != 0)
      throw NoMemory("Vector<VectFull>::Vector(int)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(i*sizeof(T)) + " bytes ("
		     + to_str(i) + " elements).");
#endif

  }


  //! Builds a vector using pre-existing data.
  /*!
    \param i length of the vector.
    \param data the data array. \a data contains the elements of the vector
    and must therefore contain \a i elements.
    \warning \a data has to be used carefully outside the object. Unless you
    use 'Nullify', \a data will be freed by the destructor, which means that
    \a data must have been allocated carefully. The vector allocator should be
    compatible.
    \note This constructor should only be used by advanced users.
  */
  template <class T, class Allocator>
  Vector<T, VectFull, Allocator>
  ::Vector(int i, typename Vector<T, VectFull, Allocator>::pointer data):
    Vector_Base<T, Allocator>()
  {
    SetData(i, data);
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, VectFull, Allocator>::
  Vector(const Vector<T, VectFull, Allocator>& V):
    Vector_Base<T, Allocator>(V)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = this->vect_allocator_.allocate(V.GetM(), this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->data_ = NULL;
      }
    if (this->data_ == NULL)
      this->m_ = 0;
    if (this->data_ == NULL && V.GetM() != 0)
      throw NoMemory("Vector<VectFull>::Vector(Vector<VectFull>&)",
		     string("Unable to allocate memory for a vector of size ")
		     + to_str(V.GetM()*sizeof(T)) + " bytes ("
		     + to_str(V.GetM()) + " elements).");
#endif

    this->vect_allocator_.memorycpy(this->data_, V.GetData(), V.GetM());

  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  Vector<T, VectFull, Allocator>::~Vector()
  {
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector. - not any more it doesn't - this was being called,
    destroying the vector, then further calls to the destroyed vector were 
    being used (in this file).
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  inline void Vector_Base<T, Allocator>::Clear()
  {
    // this->~Vector();
    if (data_ != NULL)
    {
      vect_allocator_.deallocate(data_, m_);
      m_ = 0;
      data_ = NULL;
    } 
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Reallocate(int i)
  {

    if (i != this->m_)
      {

	this->m_ = i;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    this->data_ =
	      reinterpret_cast<pointer>(this->vect_allocator_
					.reallocate(this->data_, i, this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->data_ = NULL;
	    return;
	  }
	if (this->data_ == NULL)
	  {
	    this->m_ = 0;
	    return;
	  }
#endif

      }
  }


  //! Changes the length of the vector, and keeps previous values.
  /*!
    Reallocates the vector to size i. Previous values are kept.
    \param n new length of the vector.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Resize(int n)
  {

    if (n == this->m_)
      return;

    Vector<T, VectFull, Allocator> X_new(n);
    for (int i = 0; i < min(this->m_, n); i++)
      X_new(i) = this->data_[i];

    SetData(n, X_new.GetData());
    X_new.Nullify();
  }


  //! Changes the length of the vector and sets its data array
  //! (low level method).
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param i new length of the vector.
    \param data the new data array. 'data' contains the new elements of the
    vector and must therefore contain 'i' elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The vector
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>
  ::SetData(int i, typename Vector<T, VectFull, Allocator>::pointer data)
  {
    this->Clear();

    this->m_ = i;

    this->data_ = data;
  }


  //! Lets the current vector point to the data of another vector.
  /*! Sets the length and the data of the current vector to that of \a V. On
    exit, the current vector shares it data array with \a V.
    \param V vector whose data should be shared with current instance.
    \warning On exit, \a V and the current instance share the same data array,
    and they are likely to free it. As a consequence, before the vectors are
    destructed, it is necessary to call 'Nullify' on either \a V or the
    current instance. In addition, if the current instance is to deallocate
    the data array, its allocator should be compatible with the allocator that
    allocated the data array (probably the allocator of \a V).
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  template <class Allocator0>
  inline void Vector<T, VectFull, Allocator>
  ::SetData(const Vector<T, VectFull, Allocator0>& V)
  {
    SetData(V.GetLength(), V.GetData());
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->data_ = NULL;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::reference
  Vector<T, VectFull, Allocator>::operator() (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectFull>::operator()",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access to element \a i.
  /*!
    \param i index.
    \return The value of the vector at \a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::reference
  Vector<T, VectFull, Allocator>::Get(int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectFull>::Get(int)",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access operator.
  /*!
    \param i index.
    \return The value of the vector at 'i'.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::const_reference
  Vector<T, VectFull, Allocator>::operator() (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectFull>::operator() const",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Access to element \a i.
  /*!
    \param i index.
    \return The value of the vector at \a i.
  */
  template <class T, class Allocator>
  inline typename Vector<T, VectFull, Allocator>::const_reference
  Vector<T, VectFull, Allocator>::Get(int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<VectFull>::Get(int) const",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif

    return this->data_[i];
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator= (const Vector<T, VectFull, Allocator>& X)
  {
    this->Copy(X);

    return *this;
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>
  ::Copy(const Vector<T, VectFull, Allocator>& X)
  {
    this->Reallocate(X.GetLength());

    this->vect_allocator_.memorycpy(this->data_, X.GetData(), this->m_);
  }


  //! Duplicates a vector.
  /*!
    \return A copy of the vector.
    \note Memory is duplicated: the returned vector is therefore independent
    from the current instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, VectFull, Allocator>
  Vector<T, VectFull, Allocator>::Copy() const
  {
    return Vector<T, VectFull, Allocator>(*this);
  }


  //! Multiplies a vector by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Allocator> template<class T0>
  inline Vector<T, VectFull, Allocator>& Vector<T, VectFull, Allocator>
  ::operator*= (const T0& alpha)
  {
    for (int i = 0; i < this->m_; i++)
      this->data_[i] *= alpha;

    return *this;
  }


  //! Appends an element to the vector.
  /*!
    \param x element to be appended.
    \warning This method will only work if the allocator preserves the
    elements while reallocating.
  */
  template <class T, class Allocator>
  inline void Vector<T, VectFull, Allocator>::Append(const T& x)
  {
    int i = this->GetLength();
    this->Reallocate(i + 1);
    this->data_[i] = x;
  }


  //! Appends an element at the end of the vector.
  /*!
    \param x element to be appended.
  */
  template <class T, class Allocator> template<class T0>
  inline void Vector<T, VectFull, Allocator>::PushBack(const T0& x)
  {
    Resize(this->m_+1);
    this->data_[this->m_-1] = x;
  }


  //! Appends a vector X at the end of the vector.
  /*!
    \param X vector to be appended.
  */
  template <class T, class Allocator> template<class Allocator0>
  inline void Vector<T, VectFull, Allocator>
  ::PushBack(const Vector<T, VectFull, Allocator0>& X)
  {
    int Nold = this->m_;
    Resize(this->m_ + X.GetM());
    for (int i = 0; i < X.GetM(); i++)
      this->data_[Nold+i] = X(i);
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  int Vector<T, VectFull, Allocator>::GetDataSize()
  {
    return this->m_;
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the vector stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::Zero()
  {
    this->vect_allocator_.memoryset(this->data_, char(0),
				    this->GetDataSize() * sizeof(value_type));
  }


  //! Fills the vector with 0, 1, 2, ...
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::Fill()
  {
    for (int i = 0; i < this->m_; i++)
      this->data_[i] = i;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  void Vector<T, VectFull, Allocator>::Fill(const T0& x)
  {
    for (int i = 0; i < this->m_; i++)
      this->data_[i] = x;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, VectFull, Allocator>&
  Vector<T, VectFull, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the vector randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::FillRand()
  {
    srand(time(NULL));
    for (int i = 0; i < this->m_; i++)
      this->data_[i] = rand();
  }


  //! Displays the vector.
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::Print() const
  {
    for (int i = 0; i < this->GetLength(); i++)
      cout << (*this)(i) << "\t";
    cout << endl;
  }


  /*********
   * NORMS *
   *********/


  //! Returns the infinite norm.
  /*!
    \return The infinite norm.
  */
  template <class T, class Allocator>
  typename Vector<T, VectFull, Allocator>::value_type
  Vector<T, VectFull, Allocator>::GetNormInf() const
  {
    value_type res = value_type(0);
    for (int i = 0; i < this->GetLength(); i++)
      {
	res = max(res, this->data_[i]);
	res = max(res, T(-this->data_[i]));
      }

    return res;
  }


  //! Returns the index of the highest absolute value.
  /*!
    \return The index of the element that has the highest absolute value.
  */
  template <class T, class Allocator>
  int Vector<T, VectFull, Allocator>::GetNormInfIndex() const
  {

#ifdef SELDON_CHECK_DIMENSIONS
    if (this->GetLength() == 0)
      throw WrongDim("Vector<VectFull>::GetNormInfIndex()",
		     "Vector is null.");
#endif

    value_type res = value_type(0), temp;
    int j = 0;
    for (int i = 0; i < this->GetLength(); i++)
      {
	temp = res;
	res = max(res, this->data_[i]);
	res = max(res, T(-this->data_[i]));
	if (temp != res) j = i;
      }

    return j;
  }


  /**************************
   * OUTPUT/INPUT FUNCTIONS *
   **************************/


  //! Writes the vector in a file.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::Write(string FileName, "
                    "bool with_size)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Write(ostream& FileStream, "
                    "bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		       sizeof(int));

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     this->m_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Write(ostream& FileStream, "
                    "bool with_size)",
                    "Output operation failed.");
#endif

  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the vector in a file stream.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    if (this->GetLength() != 0)
      FileStream << (*this)(0);

    for (int i = 1; i < this->GetLength(); i++)
      FileStream << "\t" << (*this)(i);

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Sets the vector from a file.
  /*!
    Sets the vector according to a binary file that stores the length of the
    vector (integer) and all elements.
    \param FileName file name.
    \param with_size if set to 'false', the length of the vector is not
    available in the file. In this case, the current size N of the vector is
    unchanged, and N elements are read in the file.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::Read(string FileName, bool with_size)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream, with_size);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets the vector according to a binary file stream that stores the length
    of the vector (integer) and all elements.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not
    available in the stream. In this case, the current size N of the vector is
    unchanged, and N elements are read in the stream.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Read(istream& FileStream, "
                    "bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int new_size;
        FileStream.read(reinterpret_cast<char*>(&new_size), sizeof(int));
        this->Reallocate(new_size);
      }

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    this->GetLength() * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::Read(istream& FileStream, "
                    "bool with_size)",
                    "Output operation failed.");
#endif

  }


  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<VectFull>::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, VectFull, Allocator>::ReadText(istream& FileStream)
  {
    // Previous values of the vector are cleared.
    this->Clear();

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Vector<VectFull>::ReadText(istream& FileStream)",
                    "The stream is not ready.");
#endif

    T entry;
    int number_element = 0;
    while (!FileStream.eof())
      {
	// Reads a new entry.
	FileStream >> entry;

	if (FileStream.fail())
	  break;
	else
	  {
	    number_element++;

	    // If needed, resizes the vector. Its size is already doubled so
	    // that the vector should be resized a limited number of times.
	    if (number_element > this->m_)
	      this->Resize(2 * number_element);

	    this->data_[number_element - 1] = entry;
	  }
      }

    // Resizes to the actual size.
    if (number_element > 0)
      this->Resize(number_element);
    else
      this->Clear();
  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Storage, Allocator>& V)
  {
    for (int i = 0; i < V.GetLength() - 1; i++)
      out << V(i) << '\t';
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);

    return out;
  }


} // namespace Seldon.

#define SELDON_FILE_VECTOR_CXX
#endif
