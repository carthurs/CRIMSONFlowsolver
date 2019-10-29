// Copyright (C) 2011-2012, INRIA
// Author(s): Marc Fragu
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


#ifndef SELDON_FILE_VECTOR_PETSCVECTOR_CXX


#include "PetscVector.hxx"


namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  PETScVector<T, Allocator>::PETScVector():
    Vector_Base<T, Allocator>()
  {
    this->m_ = 0;
    petsc_vector_deallocated_ = false;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  PETScVector<T, Allocator>::PETScVector(int i, MPI_Comm mpi_communicator)
    :Vector_Base<T, Allocator>(i)
  {
    this->m_ = i;
    petsc_vector_deallocated_ = false;
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] petsc_vector vector to be copied.
  */
  template <class T, class Allocator>
  PETScVector<T, Allocator>::
  PETScVector(Vec& petsc_vector)
  {
    petsc_vector_deallocated_ = true;
    Copy(petsc_vector);
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] V vector to be copied.
  */
  template <class T, class Allocator>
  PETScVector<T, Allocator>::
  PETScVector(const PETScVector<T, Allocator>& V)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  PETScVector<T, Allocator>::~PETScVector()
  {
    Clear();
  }


  //! Returns a reference on the inner petsc vector.
  /*!
    \return a reference on the inner petsc vector.
  */
  template <class T, class Allocator>
  Vec& PETScVector<T, Allocator>::GetPetscVector()
  {
    return petsc_vector_;
  }


  //! Returns a const reference on the inner petsc vector.
  /*!
    \return a const reference on the inner petsc vector.
  */
  template <class T, class Allocator>
  const Vec& PETScVector<T, Allocator>::GetPetscVector() const
  {
    return petsc_vector_;
  }


  //! Sets the MPI communicator.
  /*!
    \param[in] mpi_communicator the mpi communicator to be set.
  */
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::SetCommunicator(MPI_Comm mpi_communicator)
  {
    mpi_communicator_ = mpi_communicator;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Clears the vector.
  /*!
    Destructs the vector.
    \warning On exit, the vector is an empty vector.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Clear()
  {
    if (petsc_vector_deallocated_)
      return;
    int ierr;
    ierr = VecDestroy(&petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    petsc_vector_deallocated_ = true;
  }


  //! Changes the length of the vector, and keeps previous values.
  /*!
    Reallocates the vector to size i. Previous values are kept.
    \param[in] n new length of the vector.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Resize(int n)
  {
    throw Undefined("PETScVector<T, Allocator>::Resize(int n)");
  }


  /*! \brief Changes the length of the vector and sets its data array (low
    level method). */
  /*!
    Reallocates a vector and sets the new data array. It is useful to create
    a vector from pre-existing data.
    \param[in] i new length of the vector.
    \param[in] data the new data array. 'data' contains the new elements of
    the vector and must therefore contain 'i' elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The vector
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::SetData(int i, typename PETScVector<T, Allocator>::pointer data)
  {
    throw Undefined("PETScVector<T, Allocator>::SetData(int i, "
                    "typename PETScVector<T, Allocator>::pointer data)");
  }


  //! Clears the vector without releasing memory.
  /*!
    On exit, the vector is empty and the memory has not been released.
    It is useful for low level manipulations on a Vector instance.
    \warning Memory is not released.
  */
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::Nullify()
  {
    throw Undefined("PETScVector<T, Allocator>::Nullify()");
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
  inline typename PETScVector<T, Allocator>::value_type
  PETScVector<T, Allocator>::operator() (int i) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongIndex("Vector<PETSc>::operator()",
		       string("Index should be in [0, ") + to_str(this->m_-1)
		       + "], but is equal to " + to_str(i) + ".");
#endif
    int ierr;
    value_type ret[1];
    int index[1];
    index[0] = i;
    ierr = VecGetValues(petsc_vector_, 1, index, ret);
    CHKERRABORT(mpi_communicator_, ierr);
    return ret[0];
  }


  //! Inserts or adds values into certain locations of a vector.
  /*! \warning These values may be cached, so 'Flush' must be called after
    all calls to SetBuffer() have been completed.
    \param[in] i index where to insert the value.
    \param[in] value the value to insert.
    \param[in] insert_mode either INSERT_VALUES or ADD_VALUES, where
    ADD_VALUES adds the value to the entry, and INSERT_VALUES replaces
    existing entry with new value.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::SetBuffer(int i, T value, InsertMode insert_mode)
  {
    int ierr;
    int ix[1] = {i};
    T data[1] = {value};
    ierr = VecSetValues(petsc_vector_, 1, ix,
                        data, INSERT_VALUES);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Assembles the PETSc vector.
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::Flush()
  {
    int ierr;
    ierr = VecAssemblyBegin(petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = VecAssemblyEnd(petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Returns the range of indices owned by this processor.
  /*! The vectors are laid out with the first \f$n_1\f$ elements on the first
    processor, next \f$n_2\f$ elements on the second, etc. If the current
    processor is \f$k\f$, this method returns \f$n_k\f$ in \a i and
    \f$n_{k+1}\f$ in \a j. If \a i is set to PETSC_NULL on entry, it is not
    modified by this function. Same is true for \a j.
    \param[in,out] i the index of the first local element.
    \param[in,out] j the index of the last local element, plus 1.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::GetProcessorRange(int& i, int& j) const
  {
    int ierr;
    ierr = VecGetOwnershipRange(petsc_vector_, &i, &j);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::Copy(const PETScVector<T, Allocator>& X)
  {
    Copy(X.GetPetscVector());
  }


  //! Duplicates a vector.
  /*!
    \param petsc_vector vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>
  ::Copy(const Vec& petsc_vector)
  {
    Clear();

    int ierr;
    ierr = VecGetSize(petsc_vector, &this->m_);
    CHKERRABORT(mpi_communicator_, ierr);
    PetscObjectGetComm(reinterpret_cast<PetscObject>(petsc_vector),
                       &mpi_communicator_);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = VecDuplicate(petsc_vector, &petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    ierr = VecCopy(petsc_vector, petsc_vector_);
    CHKERRABORT(mpi_communicator_, ierr);
    petsc_vector_deallocated_ = false;
  }


  //! Appends an element to the vector.
  /*!
    \param x element to be appended.
    \warning This method will only work if the allocator preserves the
    elements while reallocating.
  */
  template <class T, class Allocator>
  inline void PETScVector<T, Allocator>::Append(const T& x)
  {
    throw Undefined("PETScVector<T, Allocator>::Append(const T& x)");
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  int PETScVector<T, Allocator>::GetDataSize()
  {
    return this->m_;
  }


  //! Returns the number of elements stored.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Allocator>
  int PETScVector<T, Allocator>::GetLocalM()
  {
    int size;
    VecGetLocalSize(petsc_vector_, &size);
    return size;
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
  void PETScVector<T, Allocator>::Zero()
  {
    int ierr;
    ierr = VecSet(petsc_vector_, 0.);
    CHKERRABORT(mpi_communicator_, ierr);
  }


  //! Fills the vector with 0, 1, 2, ...
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::Fill()
  {
    int ierr;
    Vector<int> index(this->m_);
    Vector<double> data(this->m_);
    index.Fill();
    data.Fill();
    ierr = VecSetValues(petsc_vector_, this->m_, index.GetData(),
                        data.GetData(), INSERT_VALUES);
    CHKERRABORT(mpi_communicator_, ierr);
    Flush();
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  void PETScVector<T, Allocator>::Fill(const T0& x)
  {
    int ierr;
    ierr = VecSet(petsc_vector_, double(x));
    CHKERRABORT(mpi_communicator_, ierr);
    Flush();
  }


  //! Fills the vector randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Allocator>
  void PETScVector<T, Allocator>::FillRand()
  {
    int ierr;
    Vector<int> index(this->m_);
    Vector<double> data(this->m_);
    index.Fill();
    data.FillRand();
    ierr = VecSetValues(petsc_vector_, this->m_, index.GetData(),
                        data.GetData(), INSERT_VALUES);
    CHKERRABORT(mpi_communicator_, ierr);
    Flush();
  }


  /*********
   * NORMS *
   *********/


  //! Returns the infinite norm.
  /*!
    \return The infinite norm.
  */
  template <class T, class Allocator>
  typename PETScVector<T, Allocator>::value_type
  PETScVector<T, Allocator>::GetNormInf() const
  {
    int ierr;
    value_type res;
    ierr = VecNorm(petsc_vector_, NORM_INFINITY, &res);
    CHKERRABORT(mpi_communicator_, ierr);
    return res;
  }


  //! Returns the index of the highest absolute value.
  /*!
    \return The index of the element that has the highest absolute value.
  */
  template <class T, class Allocator>
  int PETScVector<T, Allocator>::GetNormInfIndex() const
  {
    throw Undefined("PETScVector<T, Allocator>::GetNormInfIndex()");
  }


  /////////////////////////////
  // SEQUENTIAL PETSC VECTOR //
  /////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  Vector<T, PETScSeq, Allocator>::Vector():
    PETScVector<T, Allocator>()
  {
    this->mpi_communicator_ = MPI_COMM_WORLD;
    this->m_ = 0;
    int ierr;
    ierr = VecCreateSeq(PETSC_COMM_SELF, 0, &this->petsc_vector_);
    this->petsc_vector_deallocated_ = false;
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  Vector<T, PETScSeq, Allocator>::Vector(int i, MPI_Comm mpi_communicator)
    :PETScVector<T, Allocator>(i)
  {
    int ierr;
    this->mpi_communicator_ = mpi_communicator;
    this->m_ = i;
    ierr = VecCreateSeq(PETSC_COMM_SELF, i, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = VecSet(this->petsc_vector_, 0.);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Flush();
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param|in] petsc_vector vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, PETScSeq, Allocator>::
  Vector(Vec& petsc_vector): PETScVector<T, Allocator>(petsc_vector)
  {
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param V vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, PETScSeq, Allocator>::
  Vector(const Vector<T, PETScSeq, Allocator>& V)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  Vector<T, PETScSeq, Allocator>::~Vector()
  {
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Copy(const Vector<T, PETScSeq, Allocator>& X)
  {
    Copy(X.GetPetscVector());
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Copy(const Vec& X)
  {
    PETScVector<T, Allocator>::Copy(X);
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScSeq, Allocator>
  ::Reallocate(int i)
  {
    this->Clear();
    int ierr;
    this->m_ = i;
    ierr = VecCreateSeq(PETSC_COMM_SELF, i, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Fill(T(0));
    this->Flush();
    this->petsc_vector_deallocated_ = false;
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScSeq, Allocator>& Vector<T, PETScSeq, Allocator>
  ::operator= (const Vector<T, PETScSeq, Allocator>& X)
  {
    this->Copy(X);
    return *this;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, PETScSeq, Allocator>&
  Vector<T, PETScSeq, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
    return *this;
  }


  //! Multiplies a vector by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Allocator>
  template<class T0>
  inline Vector<T, PETScSeq, Allocator>& Vector<T, PETScSeq, Allocator>
  ::operator*= (const T0& alpha)
  {
    int ierr;
    ierr = VecScale(this->petsc_vector_, alpha);
    CHKERRABORT(this->mpi_communicator_, ierr);
    return *this;
  }


  //! Displays the vector.
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::Print() const
  {
    int ierr;
    ierr = VecView(this->petsc_vector_, PETSC_VIEWER_STDOUT_SELF);
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
  void Vector<T, PETScSeq, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::Write(string FileName, bool with_size) const");
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::Write(string FileName, bool with_size) const");
  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScSeq, Allocator>"
                    "::WriteText(string FileName)",
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
  void Vector<T, PETScSeq, Allocator>::WriteText(ostream& FileStream) const
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::WriteText(ostream& FileStream) const");
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
  void Vector<T, PETScSeq, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScSeq, Allocator>::Read(string FileName)",
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
  void Vector<T, PETScSeq, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::Read(istream& FileStream, bool with_size)");
  }


  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::ReadText(string FileName)
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::ReadText(string FileName)");
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScSeq, Allocator>::ReadText(istream& FileStream)
  {
    throw Undefined("Vector<T, PETScSeq, Allocator>"
                    "::ReadText(istream& FileStream)");
  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
                        const Vector<T, PETScSeq, Allocator>& V)
  {
    out << '\t';
    for (int i = 0; i < V.GetLength() - 1; i++)
      out << V(i) << '\t';
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);
    return out;
  }


  ///////////////////////////
  // PARALLEL PETSC VECTOR //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the vector is empty.
  */
  template <class T, class Allocator>
  Vector<T, PETScPar, Allocator>::Vector():
    PETScVector<T, Allocator>()
  {
    this->mpi_communicator_ = MPI_COMM_WORLD;
    int ierr;
    ierr = VecCreateMPI(this->mpi_communicator_,
                        PETSC_DECIDE, 0, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  Vector<T, PETScPar, Allocator>::Vector(int i, MPI_Comm mpi_communicator):
    PETScVector<T, Allocator>(i)
  {
    int ierr;
    this->mpi_communicator_ = mpi_communicator;
    ierr = VecCreateMPI(this->mpi_communicator_, PETSC_DECIDE, i,
                        &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = VecSet(this->petsc_vector_, 0.);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Flush();
  }


  //! Main constructor.
  /*! Builds a vector of a given size.
    \param[in] i length of the vector.
    \param[in] Nlocal size of local vector.
    \param[in] mpi_communicator MPI communicator to use.
  */
  template <class T, class Allocator>
  Vector<T, PETScPar, Allocator>::Vector(int i, int Nlocal, MPI_Comm
                                         mpi_communicator):
    PETScVector<T, Allocator>(i)
  {
    int ierr;
    this->mpi_communicator_ = mpi_communicator;
    ierr = VecCreateMPI(this->mpi_communicator_, Nlocal,
                        i, &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    ierr = VecSet(this->petsc_vector_, 0.);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Flush();
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] petsc_vector vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, PETScPar, Allocator>::
  Vector(Vec& petsc_vector): PETScVector<T, Allocator>(petsc_vector)
  {
  }


  //! Copy constructor.
  /*! Builds a copy of a vector.
    \param[in] V vector to be copied.
  */
  template <class T, class Allocator>
  Vector<T, PETScPar, Allocator>::
  Vector(const Vector<T, PETScPar, Allocator>& V)
  {
    Copy(V);
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Allocator>
  Vector<T, PETScPar, Allocator>::~Vector()
  {
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Copy(const Vector<T, PETScPar, Allocator>& X)
  {
    Copy(X.GetPetscVector());
  }


  //! Duplicates a vector.
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Copy(const Vec& X)
  {
    PETScVector<T, Allocator>::Copy(X);
  }


  //! Duplicates a vector (assignment operator).
  /*!
    \param X vector to be copied.
    \note Memory is duplicated: 'X' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Allocator>
  inline Vector<T, PETScPar, Allocator>& Vector<T, PETScPar, Allocator>
  ::operator= (const Vector<T, PETScPar, Allocator>& X)
  {
    this->Copy(X);
    return *this;
  }


  //! Fills the vector with a given value.
  /*!
    \param x value to fill the vector with.
  */
  template <class T, class Allocator>
  template <class T0>
  Vector<T, PETScPar, Allocator>&
  Vector<T, PETScPar, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
    return *this;
  }


  //! Multiplies a vector by a scalar.
  /*!
    \param alpha scalar.
  */
  template <class T, class Allocator>
  template<class T0>
  inline Vector<T, PETScPar, Allocator>& Vector<T, PETScPar, Allocator>
  ::operator*= (const T0& alpha)
  {
    int ierr;
    ierr = VecScale(this->petsc_vector_, alpha);
    CHKERRABORT(this->mpi_communicator_, ierr);
    return *this;
  }


  //! Displays the vector.
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::Print() const
  {
    int ierr;
    ierr = VecView(this->petsc_vector_, PETSC_VIEWER_STDOUT_WORLD);
  }


  //! Vector reallocation.
  /*!
    The vector is resized.
    \param i new length of the vector.
    \warning Depending on your allocator, initial elements of the vector may
    be lost.
  */
  template <class T, class Allocator>
  inline void Vector<T, PETScPar, Allocator>
  ::Reallocate(int i, int local_size)
  {
    this->Clear();
    int ierr;
    this->m_ = i;
    ierr = VecCreateMPI(this->mpi_communicator_, local_size, i,
                        &this->petsc_vector_);
    CHKERRABORT(this->mpi_communicator_, ierr);
    this->Fill(T(0));
    this->Flush();
    this->petsc_vector_deallocated_ = false;
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
  void Vector<T, PETScPar, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    throw Undefined("Vector<T, PETScPar, Allocator>"
                    "::Write(string FileName, bool with_size) const");
  }


  //! Writes the vector in a file stream.
  /*!
    The length of the vector (integer) and all elements of the vector are
    stored in binary format.
    \param FileStream file stream.
    \param with_size if set to 'false', the length of the vector is not saved.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {
    int local_n;
    VecGetLocalSize(this->petsc_vector_, &local_n);

    Vector<int> index(local_n);
    index.Fill();
    int i_start, i_end;
    this->GetProcessorRange(i_start, i_end);
    for (int i = 0; i < local_n; i++)
      index(i) += i_start;
    Vector<T> data(local_n);
    VecGetValues(this->petsc_vector_, local_n, index.GetData(),
                 data.GetData());
    data.Write(FileStream, with_size);
  }


  //! Writes the vector in a file.
  /*!
    All elements of the vector are stored in text format. The length is not
    stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScPar, Allocator>"
                    "::WriteText(string FileName)",
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
  void Vector<T, PETScPar, Allocator>::WriteText(ostream& FileStream) const
  {
    int local_n;
    VecGetLocalSize(this->petsc_vector_, &local_n);

    Vector<int> index(local_n);
    index.Fill();
    int i_start, i_end;
    this->GetProcessorRange(i_start, i_end);
    for (int i = 0; i < local_n; i++)
      index(i) += i_start;
    Vector<T> data(local_n);
    VecGetValues(this->petsc_vector_, local_n, index.GetData(),
                 data.GetData());
    data.WriteText(FileStream);
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
  void Vector<T, PETScPar, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Vector<T, PETScPar, Allocator>::Read(string FileName)",
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
  void Vector<T, PETScPar, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {
    throw Undefined("PETScVector<T, PETScPar, Allocator>"
                    "::Read(istream& FileStream, bool with_size)");
  }


  //! Sets the vector from a file.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileName file name.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::ReadText(string FileName)
  {
    throw Undefined("PETScVector<T, PETScPar, Allocator>"
                    "::ReadText(string FileName)");
  }


  //! Sets the vector from a file stream.
  /*!
    Sets all elements of the vector according to a text format. The length is
    not stored.
    \param FileStream file stream.
  */
  template <class T, class Allocator>
  void Vector<T, PETScPar, Allocator>::ReadText(istream& FileStream)
  {
    throw Undefined("PETScVector<T, PETScPar, Allocator>"
                    "::ReadText(istream& FileStream)");
  }


  //! operator<< overloaded for vectors.
  /*!
    \param out output stream.
    \param V vector to be put in the stream.
    \return The updated stream.
  */
  template <class T, class Allocator>
  ostream& operator << (ostream& out,
                        const Vector<T, PETScPar, Allocator>& V)
  {
    out << '\t';
    for (int i = 0; i < V.GetLength() - 1; i++)
      out << V(i) << '\t';
    if (V.GetLength() != 0)
      out << V(V.GetLength() - 1);
    return out;
  }


} // namespace Seldon.


#define SELDON_FILE_PETSCVECTOR_CXX
#endif
