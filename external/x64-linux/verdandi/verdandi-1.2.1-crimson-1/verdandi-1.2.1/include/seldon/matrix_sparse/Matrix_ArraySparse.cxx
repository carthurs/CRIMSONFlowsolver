// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_ARRAY_SPARSE_CXX

#include "Matrix_ArraySparse.hxx"

namespace Seldon
{

  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>::Matrix_ArraySparse()
    : val_()
  {
    this->m_ = 0;
    this->n_ = 0;
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Matrix_ArraySparse(int i, int j) :
    val_(Storage::GetFirst(i, j))
  {
    this->m_ = i;
    this->n_ = j;
  }


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocat>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocat>::~Matrix_ArraySparse()
  {
    this->m_ = 0;
    this->n_ = 0;
  }


  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix is
    empty (0 by 0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_ArraySparse();
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Data is lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Reallocate(int i, int j)
  {
    // Clears previous entries.
    Clear();

    this->m_ = i;
    this->n_ = j;

    int n = Storage::GetFirst(i, j);
    val_.Reallocate(n);
  }


  //! Reallocates additional memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    Data is kept
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    int n = Storage::GetFirst(this->m_, n_);
    int new_n = Storage::GetFirst(i, j);
    if (n != new_n)
      {
	Vector<Vector<T, VectSparse, Allocator>, VectFull,
	  NewAlloc<Vector<T, VectSparse, Allocator> > > new_val;

	new_val.Reallocate(new_n);

	for (int k = 0 ; k < min(n, new_n) ; k++)
	  Swap(new_val(k), this->val_(k));

	val_.SetData(new_n, new_val.GetData());
	new_val.Nullify();

      }

    this->m_ = i;
    this->n_ = j;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of rows.
  /*!
    \return the number of rows.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetM() const
  {
    return m_;
  }


  //! Returns the number of columns.
  /*!
    \return the number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetN() const
  {
    return n_;
  }


  //! Returns the number of rows of the matrix possibly transposed.
  /*!
    \param status assumed status about the transposition of the matrix.
    \return The number of rows of the possibly-transposed matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::GetM(const SeldonTranspose& status) const
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
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::GetN(const SeldonTranspose& status) const
  {
    if (status.NoTrans())
      return n_;
    else
      return m_;
  }


  //! Returns the number of non-zero entries.
  /*!
    \return The number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetNonZeros()
    const
  {
    int nnz = 0;
    for (int i = 0; i < this->val_.GetM(); i++)
      nnz += this->val_(i).GetM();

    return nnz;
  }


  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the number of non-zero entries.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetDataSize()
    const
  {
    return GetNonZeros();
  }


  //! Returns (row or column) indices of non-zero entries in row
  /*!
    \param[in] i row (or column) number.
    \return The array of column (or row) indices of non-zero entries
    of row (or column) i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int* Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetIndex(int i)
    const
  {
    return val_(i).GetIndex();
  }


  //! Returns values of non-zero entries of a row/column.
  /*!
    \param[in] i row (or column) number.
    \return The array of values of non-zero entries of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T*
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetData(int i) const
  {
    return val_(i).GetData();
  }


  //! Returns values of non-zero entries.
  /*!
    \return Array of sparse rows
    There is a different array for each row/column.
  */
  template <class T, class Prop, class Storage, class Allocat>
  inline Vector<T, VectSparse, Allocat>*
  Matrix_ArraySparse<T, Prop, Storage, Allocat>::GetData() const
  {
    return val_.GetData();
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param[in] i row index.
    \param[in] j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "],but is equal to " + to_str(i) + ".");

    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif

    return this->val_(Storage::GetFirst(i, j))(Storage::GetSecond(i, j));
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "],but is equal to " + to_str(i) + ".");

    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif

    return this->val_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "],but is equal to " + to_str(i) + ".");

    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif

    return this->val_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j));
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "], but is equal to " + to_str(i) + ".");

    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif

    return
      this->val_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->m_-1) +
		     "], but is equal to " + to_str(i) + ".");

    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ArraySparse::operator()",
		     "Index should be in [0, " + to_str(this->n_-1) +
		     "], but is equal to " + to_str(j) + ".");
#endif

    return
      this->val_(Storage::GetFirst(i, j)).Val(Storage::GetSecond(i, j));
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    this->Get(i, j) = x;
  }


  //! Returns j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return j-th non-zero entry of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const T& Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Value (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Storage::GetFirst(this->m_, this->n_))
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(Storage::GetFirst(this->m_, this->n_)-1)
                     + "], but is equal to " + to_str(i) + ".");

    if ((j < 0)||(j >= this->val_(i).GetM()))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->val_(i).GetM()-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return val_(i).Value(j);
  }


  //! Returns j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return j-th non-zero entry of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline T&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::Value (int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Storage::GetFirst(this->m_, this->n_))
      throw WrongRow("Matrix_ArraySparse::value", "Index should be in [0, "
		     + to_str(Storage::GetFirst(this->m_, this->n_)-1)
                     + "], but is equal to " + to_str(i) + ".");

    if ((j < 0)||(j >= this->val_(i).GetM()))
      throw WrongCol("Matrix_ArraySparse::value", "Index should be in [0, " +
		     to_str(this->val_(i).GetM()-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return val_(i).Value(j);
  }


  //! Returns column/row number of j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return Column/row number of j-th non-zero value of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  int Matrix_ArraySparse<T, Prop, Storage, Allocator>::Index(int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Storage::GetFirst(this->m_, this->n_))
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(Storage::GetFirst(this->m_, this->n_)-1)
                     + "], but is equal to " + to_str(i) + ".");

    if ((j < 0)||(j >= this->val_(i).GetM()))
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->val_(i).GetM()-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return val_(i).Index(j);
  }


  //! Returns column/row number of j-th non-zero value of row/column i.
  /*!
    \param[in] i row/column number.
    \param[in] j local number.
    \return Column/row number of j-th non-zero value of row/column i.
  */
  template <class T, class Prop, class Storage, class Allocator> inline
  int& Matrix_ArraySparse<T, Prop, Storage, Allocator>::Index(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= Storage::GetFirst(this->m_, this->n_))
      throw WrongRow("Matrix_ArraySparse::index", "Index should be in [0, "
		     + to_str(Storage::GetFirst(this->m_, this->n_)-1)
                     + "], but is equal to " + to_str(i) + ".");

    if (j < 0 || j >= this->val_(i).GetM())
      throw WrongCol("Matrix_ArraySparse::index", "Index should be in [0, " +
		     to_str(this->val_(i).GetM()-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return val_(i).Index(j);
  }


  //! Redefines a row/column of the matrix
  /*!
    \param[in] i row/col number
    \param[in] n number of non-zero entries in the row
    \param[in] val values
    \param[in] ind column numbers
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  SetData(int i, int n, T* val, int* ind)
  {
    val_(i).SetData(n, val, ind);
  }


  //!  Clears a row without releasing memory.
  /*!
    On exit, the row is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Nullify(int i)
  {
    val_(i).Nullify();
  }


  //! Redefines the matrix.
  /*!
    \param[in] m new number of rows.
    \param[in] n new number of columns.
    \param[in] val array of sparse rows/columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  SetData(int m, int n, Vector<T, VectSparse, Allocator>* val)
  {
    m_ = m;
    n_ = n;
    val_.SetData(Storage::GetFirst(m, n), val);
  }


  //!  Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Nullify()
  {
    m_ = 0;
    n_ = 0;
    val_.Nullify();
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Print() const
  {
    if (Storage::GetFirst(1, 0) == 1)
      for (int i = 0; i < this->m_; i++)
	{
	  for (int j = 0; j < this->val_(i).GetM(); j++)
	    cout << (i+1) << " " << this->val_(i).Index(j)+1
		 << " " << this->val_(i).Value(j) << endl;
	}
    else
      for (int i = 0; i < this->n_; i++)
	{
	  for (int j = 0; j < this->val_(i).GetM(); j++)
	    cout << this->val_(i).Index(j)+1 << " " << i+1
		 << " " << this->val_(i).Value(j) << endl;
	}
  }


  //! Assembles the matrix
  /*!
    All the column/row numbers are sorted.
    If same column/row numbers exist, values are added.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Assemble()
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Assemble();
  }


  //! Removes small coefficients from entries.
  /*!
    \param[in] epsilon entries whose values are below epsilon are removed.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  RemoveSmallEntry(const T0& epsilon)
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).RemoveSmallEntry(epsilon);
  }


  //! Matrix is initialized to the identity matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::SetIdentity()
  {
    this->n_ = this->m_;
    for (int i = 0; i < this->m_; i++)
      {
	val_(i).Reallocate(1);
	val_(i).Index(0) = i;
	val_(i).Value(0) = T(1);
      }
  }


  //! Non-zero entries are set to 0 (but not removed).
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Zero()
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Zero();
  }


  //! Non-zero entries are filled with values 0, 1, 2, 3 ...
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Fill()
  {
    int value = 0;
    for (int i = 0; i < val_.GetM(); i++)
      for (int j = 0; j < val_(i).GetM(); j++)
	val_(i).Value(j) = value++;
  }


  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allo> template<class T0>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allo>::Fill(const T0& x)
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Fill(x);
  }


  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  inline Matrix_ArraySparse<T, Prop, Storage, Allocator>&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
  }


  //! Non-zero entries take a random value.
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_ArraySparse<T, Prop, Storage, Allocator>::FillRand()
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).FillRand();
  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Write(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));

    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Write(FileStream);
  }


  //! Writes the matrix in a file.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row column value). 1-index convention is used.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  WriteText(string FileName) const
  {
    ofstream FileStream; FileStream.precision(14);
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row column value). 1-index convention is used.
    \param FileStream output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Write(ofstream& FileStream)",
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

    // If the last element a_{m,n} does not exist, we add a zero.
    int m = Storage::GetFirst(this->m_, this->n_);
    int n = Storage::GetSecond(this->m_, this->n_);
    if (m > 0 && n > 0)
      {
	bool presence_last_elt = false;
	if (this->val_(m-1).GetM() > 0)
	  {
	    int p = this->val_(m-1).GetM();
	    if (this->val_(m-1).Index(p-1) == n-1)
	      presence_last_elt = true;
	  }

	if (!presence_last_elt)
	  {
	    T zero;
	    SetComplexZero(zero);
	    FileStream << this->m_ << " " << this->n_ << " " << zero << '\n';
	  }
      }
  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Read(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.read(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		    sizeof(int));
    FileStream.read(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		    sizeof(int));

    val_.Reallocate(Storage::GetFirst(this->m_, this->n_));
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Read(FileStream);

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Read(istream& FileStream)",
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
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::ReadText(string FileName)",
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
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ReadText(istream& FileStream)
  {
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);

    T zero; int index = 1;
    ReadCoordinateMatrix(leaf_class, FileStream, zero, index);
  }


  /////////////////////////////
  // MATRIX<ARRAY_COLSPARSE> //
  /////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayColSparse, Allocator>(i, j)
  {
  }


  //! Clears column i.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::ClearColumn(int i)
  {
    this->val_(i).Clear();
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Alloc> inline
  void Matrix<T, Prop, ArrayColSparse, Alloc>::ReallocateColumn(int i,int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::ResizeColumn(int i,int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::SwapColumn(int i,int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  ReplaceIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSparse, Allocator>::
  GetColumnSize(int i) const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::PrintColumn(int i) const
  {
    this->val_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSparse, Allocator>::AssembleColumn(int i)
  {
    this->val_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    this->val_(j).AddInteraction(i, val);
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T, VectFull, Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      this->val_(col(j)).AddInteraction(i, val(j));
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T, VectFull, Alloc1>& val)
  {
    this->val_(i).AddInteractionRow(nb, row, val);
  }


  /////////////////////////////
  // MATRIX<ARRAY_ROWSPARSE> //
  /////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayRowSparse, Allocator>(i, j)
  {
  }


  //! Clears a row
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::ClearRow(int i)
  {
    this->val_(i).Clear();
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  ReallocateRow(int i, int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Changes the size of a row.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries of the row.
    \note Data is kept.
  */
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::ResizeRow(int i, int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two rows
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::SwapRow(int i,int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i column number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  ReplaceIndexRow(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator> inline
  int Matrix<T, Prop, ArrayRowSparse, Allocator>::GetRowSize(int i) const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a row.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::PrintRow(int i) const
  {
    this->val_(i).Print();
  }


  //! Assembles a row.
  /*!
    \param[in] i row number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::AssembleRow(int i)
  {
    this->val_(i).Assemble();
  }

  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    this->val_(i).AddInteraction(j, val);
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T, VectFull, Alloc1>& val)
  {
    this->val_(i).AddInteractionRow(nb, col, val);
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T, VectFull, Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      this->val_(row(j)).AddInteraction(i, val(j));
  }


  ////////////////////////////////
  // MATRIX<ARRAY_COLSYMSPARSE> //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayColSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayColSymSparse, Allocator>(i, j)
  {
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
  template <class T, class Prop, class Allocator>
  inline T
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    if (i < j)
      return this->val_(j)(i);

    return this->val_(i)(j);
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    if (i < j)
      return this->val_(j).Get(i);

    return this->val_(i).Get(j);
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    if (i < j)
      return this->val_(j).Get(i);

    return this->val_(i).Get(j);
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
    if (i > j)
      throw WrongArgument("Matrix::Val()", string("With this function, you ")
                          + "can only access upper part of matrix.");
#endif

    return this->val_(j).Val(i);
  }


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    return this->val_(j).Val(i);
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Allocator>
  inline void
  Matrix<T, Prop, ArrayColSymSparse, Allocator>::Set(int i, int j, const T& x)
  {
    if (i < j)
      this->val_(j).Get(i) = x;
    else
      this->val_(i).Get(j) = x;
  }


  //! Clears a column.
  template <class T, class Prop, class Allocator> inline
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::ClearColumn(int i)
  {
    this->val_(i).Clear();
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ReallocateColumn(int i, int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Reallocates column i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the column.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ResizeColumn(int i, int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two columns.
  /*!
    \param[in] i first column number.
    \param[in] j second column number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  SwapColumn(int i, int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets row numbers of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \param[in] new_index new row numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  ReplaceIndexColumn(int i, IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a column.
  /*!
    \param[in] i column number.
    \return The number of non-zero entries of the column i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  GetColumnSize(int i) const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  PrintColumn(int i) const
  {
    this->val_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AssembleColumn(int i)
  {
    this->val_(i).Assemble();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify();
    val.Nullify();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      this->val_(j).AddInteraction(i, val);
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T, VectFull, Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      if (i <= col(j))
	this->val_(col(j)).AddInteraction(i, val(j));
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T, VectFull, Alloc1>& val)
  {
    IVect new_row(nb);
    Vector<T, VectFull, Alloc1> new_val(nb);
    nb = 0;
    for (int j = 0; j < new_row.GetM(); j++)
      if (row(j) <= i)
	{
	  new_row(nb) = row(j);
	  new_val(nb) = val(j); nb++;
	}

    this->val_(i).AddInteractionRow(nb, new_row, new_val);
  }


  ////////////////////////////////
  // MATRIX<ARRAY_ROWSYMSPARSE> //
  ////////////////////////////////


  //! Default constructor.
  /*!
    Builds an empty matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Matrix():
    Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>()
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix
    \param i number of rows.
    \param j number of columns.
    \note Matrix values are not initialized.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Matrix(int i, int j):
    Matrix_ArraySparse<T, Prop, ArrayRowSymSparse, Allocator>(i, j)
  {
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
  template <class T, class Prop, class Allocator>
  inline T
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::operator() (int i, int j)
    const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_ArraySparse::operator()",
                     "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_ArraySparse::operator()",
                     "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    if (i < j)
      return this->val_(i)(j);

    return this->val_(j)(i);
  }


  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Get(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    if (i < j)
      return this->val_(i).Get(j);

    return this->val_(j).Get(i);
  }


  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Get(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
#endif

    if (i < j)
      return this->val_(i).Get(j);

    return this->val_(j).Get(i);
  }


  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
    if (i > j)
      throw WrongArgument("Matrix::Val()", string("With this function, you ")
                          + "can only access upper part of matrix.");
#endif

    return this->val_(i).Val(j);
  }


  //! Access to element (i, j)
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Allocator>
  inline const T&
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix::operator()", "Index should be in [0, "
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(j) + ".");
    if (i > j)
      throw WrongArgument("Matrix::Val()", string("With this function, you ")
                          + "can only access upper part of matrix.");
#endif

    return this->val_(i).Val(j);
  }


  //! Sets element (i, j) of the matrix
  /*!
    \param i row index.
    \param j column index.
    \param x A(i, j) = x
  */
  template <class T, class Prop, class Allocator>
  inline void
  Matrix<T, Prop, ArrayRowSymSparse, Allocator>::Set(int i, int j, const T& x)
  {
    if (i < j)
      this->val_(i).Get(j) = x;
    else
      this->val_(j).Get(i) = x;
  }


  //! Clears a row.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::ClearRow(int i)
  {
    this->val_(i).Clear();
  }


  //! Reallocates row i.
  /*!
    \param[in] i row number.
    \param[in] j new number of non-zero entries in the row.
    \warning Data may be lost.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ReallocateRow(int i,int j)
  {
    this->val_(i).Reallocate(j);
  }


  //! Reallocates row i.
  /*!
    \param[in] i column number.
    \param[in] j new number of non-zero entries in the row.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ResizeRow(int i,int j)
  {
    this->val_(i).Resize(j);
  }


  //! Swaps two rows.
  /*!
    \param[in] i first row number.
    \param[in] j second row number.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  SwapRow(int i,int j)
  {
    Swap(this->val_(i), this->val_(j));
  }


  //! Sets column numbers of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \param[in] new_index new column numbers.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  ReplaceIndexRow(int i,IVect& new_index)
  {
    for (int j = 0; j < this->val_(i).GetM(); j++)
      this->val_(i).Index(j) = new_index(j);
  }


  //! Returns the number of non-zero entries of a row.
  /*!
    \param[in] i row number.
    \return The number of non-zero entries of the row i.
  */
  template <class T, class Prop, class Allocator>
  inline int Matrix<T, Prop, ArrayRowSymSparse, Allocator>::GetRowSize(int i)
    const
  {
    return this->val_(i).GetSize();
  }


  //! Displays non-zero values of a column.
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::PrintRow(int i)
    const
  {
    this->val_(i).Print();
  }


  //! Assembles a column.
  /*!
    \param[in] i column number.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>
  ::AssembleRow(int i)
  {
    this->val_(i).Assemble();
  }


  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteraction(int i, int j, const T& val)
  {
    if (i <= j)
      this->val_(i).AddInteraction(j, val);
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<T, VectFull, Alloc1>& val)
  {
    IVect new_col(nb);
    Vector<T, VectFull, Alloc1> new_val(nb);
    nb = 0;
    for (int j = 0; j < new_col.GetM(); j++)
      if (i <= col(j))
	{
	  new_col(nb) = col(j);
	  new_val(nb) = val(j); nb++;
	}

    this->val_(i).AddInteractionRow(nb, new_col, new_val);
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator> template <class Alloc1>
  inline void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<T,VectFull,Alloc1>& val)
  {
    for (int j = 0; j < nb; j++)
      if (row(j) <= i)
        this->val_(row(j)).AddInteraction(i, val(j));
  }


  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayRowSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayColSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


  template <class T, class Prop, class Allocator>
  ostream& operator <<(ostream& out,
		       const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A)
  {
    A.WriteText(out);

    return out;
  }


} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_SPARSE_CXX
#endif
