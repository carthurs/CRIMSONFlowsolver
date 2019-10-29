// Copyright (C) 2010, INRIA
// Author(s): Marc Fragu, Vivien Mallet
// Copyright (C) 2011, Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_VECTOR2_CXX


#include "Vector2.hxx"


namespace Seldon
{


  /////////////
  // VECTOR2 //
  /////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::Vector2()
  {
  }


  //! Constructor.
  /*! The vector of vectors is allocated with \a length empty vectors.
    \param[in] length the length of the vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::Vector2(int length)
  {
    data_.Reallocate(length);
  }


  //! Constructor.
  /*! The vector of vectors and the inner vectors are allocated.
    \param[in] length the lengths of the inner vectors. The vector of vectors
    will obviously have as many elements as \a length has.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::Vector2(const Vector<int>& length)
  {
    data_.Clear();
    int m = length.GetSize();
    data_.Reallocate(m);
    for(int i = 0; i < m; i++)
      data_(i).Reallocate(length(i));
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor. The vector of vectors and the inner vectors are deallocated.
  template <class T, class Allocator0, class Allocator1>
  Vector2<T, Allocator0, Allocator1>::~Vector2()
  {
  }


  /*****************************
   * MANAGEMENT OF THE VECTORS *
   *****************************/


  //! Checks whether no elements are contained in the inner vectors.
  /*!
    \return True is no inner vector contains an element, false otherwise.
  */
  template <class T, class Allocator0, class Allocator1>
  bool Vector2<T, Allocator0, Allocator1>::IsEmpty() const
  {
    for (int i = 0; i < GetLength(); i++)
      if (GetLength(i) > 0)
        return false;
    return true;
  }


  //! Returns the size along dimension 1.
  /*!
    \return The size along dimension 1.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetSize() const
  {
    return data_.GetSize();
  }


  //! Returns the size along dimension 1.
  /*!
    \return The size along dimension 1.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetLength() const
  {
    return data_.GetLength();
  }


  //! Returns the size of the inner vector #\a i.
  /*!
    \param[in] i index.
    \return The size of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetSize(int i) const
  {
    return data_(i).GetSize();
  }


  //! Returns the size of the inner vector #\a i.
  /*!
    \param[in] i index.
    \return The size of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetLength(int i) const
  {
    return data_(i).GetLength();
  }


  //! Returns the total number of elements in the inner vectors.
  /*!
    \return The sum of the lengths of the inner vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetNelement() const
  {
    int total = 0;
    for (int i = 0; i < GetLength(); i++)
      total += GetLength(i);
    return total;
  }


  //! Returns the total number of elements in a range of inner vectors.
  /*! Returns the total number of elements in the range [\a beg, \a end[ of
    inner vectors.
    \param[in] beg inclusive lower-bound for the indexes.
    \param[in] end exclusive upper-bound for the indexes.
    \return The sum of the lengths of the inner vectors with index \a beg to
    \a end-1.
  */
  template <class T, class Allocator0, class Allocator1>
  int Vector2<T, Allocator0, Allocator1>::GetNelement(int beg, int end) const
  {
    if (beg > end)
      throw WrongArgument("Vector2::GetNelement(int beg, int end)",
                          "The lower bound of the range of inner vectors, ["
                          + to_str(beg) + ", " + to_str(end)
                          + "[, is strictly greater than its upper bound.");
    if (beg < 0 || end > GetLength())
      throw WrongArgument("Vector2::GetNelement(int beg, int end)",
        		  "The inner-vector indexes should be in [0,"
        		  + to_str(GetLength()) + "] but [" + to_str(beg)
                          + ", " + to_str(end) + "[ was provided.");

    int total = 0;
    for (int i = beg; i < end; i++)
      total += GetLength(i);
    return total;
  }


  //! Returns the shape.
  /*!
    \return A vector with the lengths of the inner vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<int> Vector2<T, Allocator0, Allocator1>::GetShape() const
  {
    Vector<int> shape(GetLength());
    for (int i = 0; i < GetLength(); i++)
      shape(i) = GetLength(i);
    return shape;
  }


  //! Returns the shape.
  /*!
    \param[out] shape the lengths of the inner vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::GetShape(Vector<int>& shape) const
  {
    shape.Reallocate(GetLength());
    for (int i = 0; i < GetLength(); i++)
      shape(i) = GetLength(i);
  }


  //! Reallocates the vector of vector.
  /*!
    \param[in] M the new size of the vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Reallocate(int M)
  {
    data_.Reallocate(M);
  }


  //! Reallocates the inner vector #\a i.
  /*!
    \param[in] i index of the inner vector to be reallocated.
    \param[in] N the new size of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Reallocate(int i, int N)
  {
    data_(i).Reallocate(N);
  }


  //! Reallocates the whole structure.
  /*!
    \param[in] length the new lengths of the inner vectors. The vector of
    vectors will obviously have as many elements as \a length has.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Reallocate(const Vector<int>& length)
  {
    int m = length.GetSize();
    data_.Reallocate(m);
    for(int i = 0; i < m; i++)
      data_(i).Reallocate(length(i));
  }


  //! Selects a range of inner vectors.
  /*! Only the inner vectors with index in [\a beg, \a end[ are kept. The
    other vectors are destroyed.
    \param[in] beg inclusive lower-bound for the indexes.
    \param[in] end exclusive upper-bound for the indexes.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Select(int beg, int end)
  {
    if (beg > end)
      throw WrongArgument("Vector2::SelectInnerVector(int beg, int end)",
                          "The lower bound of the range of inner vectors, ["
                          + to_str(beg) + ", " + to_str(end)
                          + "[, is strictly greater than its upper bound.");
    if (beg < 0 || end > GetLength())
      throw WrongArgument("Vector2::SelectInnerVector(int beg, int end)",
        		  "The inner-vector indexes should be in [0,"
        		  + to_str(GetLength()) + "] but [" + to_str(beg)
                          + ", " + to_str(end) + "[ was provided.");

    if (beg > 0)
      for (int i = 0; i < end - beg; i++)
        data_(i) = data_(beg + i);
    data_.Reallocate(end - beg);
  }


  //! Returns all values in a vector.
  /*! The output vector contains all inner vectors concatenated in the same
    order as they appear in the current Vector2 instance.
    \return All values from the current Vector2 instance.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<T, VectFull, Allocator0>
  Vector2<T, Allocator0, Allocator1>::Flatten() const
  {
    Vector<T, VectFull, Allocator0> data(GetNelement());
    int i, j, n(0);
    for (i = 0; i < GetLength(); i++)
      for (j = 0; j < GetLength(i); j++)
        data(n++) = data_(i)(j);
    return data;
  }


  //! Returns all values in a vector.
  /*! The output vector \a data contains all inner vectors concatenated in the
    same order as they appear in the current Vector2 instance.
    \param[out] data all values from the current Vector2 instance.
  */
  template <class T, class Allocator0, class Allocator1>
  template <class Td, class Allocatord>
  void Vector2<T, Allocator0, Allocator1>
  ::Flatten(Vector<Td, VectFull, Allocatord>& data) const
  {
    data.Reallocate(GetNelement());
    int i, j, n(0);
    for (i = 0; i < GetLength(); i++)
      for (j = 0; j < GetLength(i); j++)
        data(n++) = data_(i)(j);
  }


  //! Returns in a vector all values from a range of inner vectors.
  /*! The output vector \a data contains all inner vectors, in the index range
    [\a beg, \a end[, concatenated in the same order as they appear in the
    current Vector2 instance.
    \param[in] beg inclusive lower-bound for the indexes.
    \param[in] end exclusive upper-bound for the indexes.
    \param[out] data the values contained in the inner vectors [\a beg, \a
    end[.
  */
  template <class T, class Allocator0, class Allocator1>
  template <class Td, class Allocatord>
  void Vector2<T, Allocator0, Allocator1>
  ::Flatten(int beg, int end, Vector<Td, VectFull, Allocatord>& data) const
  {
    if (beg > end)
      throw WrongArgument("Vector2::Flatten(int beg, int end, Vector& data)",
                          "The lower bound of the range of inner vectors, ["
                          + to_str(beg) + ", " + to_str(end)
                          + "[, is strictly greater than its upper bound.");
    if (beg < 0 || end > GetLength())
      throw WrongArgument("Vector2::Flatten(int beg, int end, Vector& data)",
        		  "The inner-vector indexes should be in [0,"
        		  + to_str(GetLength()) + "] but [" + to_str(beg)
                          + ", " + to_str(end) + "[ was provided.");

    data.Reallocate(GetNelement(beg, end));
    int i, j, n(0);
    for (i = beg; i < end; i++)
      for (j = 0; j < GetLength(i); j++)
        data(n++) = data_(i)(j);
  }


  //! Appends an element at the end of the inner vector #\a i.
  /*!
    \param[in] i index of the inner vector to which \a x should be appended.
    \param[in] x element to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::PushBack(int i, const T& x)
  {
    data_(i).PushBack(x);
  }


  //! Appends an inner vector at the end of the vector.
  /*!
    \param[in] X vector to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::PushBack(const Vector<T, VectFull, Allocator0>& X)
  {
    data_.PushBack(X);
  }


  //! Appends a vector of vectors.
  /*! The inner vectors of \a V are appended to the current instance, in the
    same order as they appear in \a V.
    \param[in] V vector of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::PushBack(const Vector<Vector<T, VectFull, Allocator0>,
	     VectFull, Allocator1>& V)
  {
    for (int i = 0; i < V.GetLength(); i++)
      data_.PushBack(V(i));
  }


  //! Appends a vector of vectors.
  /*! The inner vectors of \a V are appended to the current instance, in the
    same order as they appear in \a V.
    \param[in] V vector of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::PushBack(const Vector2<T, Allocator0, Allocator1>& V)
  {
    PushBack(V.GetVector());
  }


  //! Clears the vector.
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Clear()
  {
    data_.Clear();
  }


  //! Clears a given vector.
  /*!
    \param[in] i index of the vector to be cleared.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Clear(int i)
  {
    data_(i).Clear();
  }


  //! Fills the vector with a given value.
  /*!
    \param[in] x value to fill the vector with.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Fill(const T& x)
  {
    for (int i = 0; i < data_.GetLength(); i++)
      data_(i).Fill(x);
  }


  //! Returns the vector of vectors.
  /*!
    \return The vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<Vector<T, VectFull, Allocator0>, VectFull, Allocator1>&
  Vector2<T, Allocator0, Allocator1>::GetVector()
  {
    return data_;
  }


  //! Returns the vector of vectors.
  /*!
    \return The vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1>
  const Vector<Vector<T, VectFull, Allocator0>, VectFull, Allocator1>
  Vector2<T, Allocator0, Allocator1>::GetVector() const
  {
    return data_;
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::GetVector(int i)
  {
    return data_(i);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  const Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::GetVector(int i) const
  {
    return data_(i);
  }


  //! Copies a Vector2 instance.
  /*!
    \param[in] V Vector2 instance to be copied.
    \note The current instance and \a V do not share memory on exit: \a V is
    duplicated in memory.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Copy(const Vector2<T, Allocator0, Allocator1>& V)
  {
    Clear();
    Reallocate(V.GetLength());
    for (int i = 0; i < V.GetLength(); i++)
      data_(i) = V(i);
  }


  /*********************************
   * ELEMENT ACCESS AND ASSIGNMENT *
   *********************************/


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  const Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::operator() (int i) const
  {
    return data_(i);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \return The inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  Vector<T, VectFull, Allocator0>&
  Vector2<T, Allocator0, Allocator1>::operator() (int i)
  {
    return data_(i);
  }


  //! Returns an element of a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the element in the inner vector #\a i.
    \return The element #\a j of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  typename Vector2<T, Allocator0, Allocator1>::const_reference
  Vector2<T, Allocator0, Allocator1>::operator() (int i, int j) const
  {
    return data_(i)(j);
  }


  //! Returns an element of a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the element in the inner vector #\a i.
    \return The element #\a j of the inner vector #\a i.
  */
  template <class T, class Allocator0, class Allocator1>
  typename Vector2<T, Allocator0, Allocator1>::reference
  Vector2<T, Allocator0, Allocator1>::operator() (int i, int j)
  {
    return data_(i)(j);
  }


  /**********************
   * CONVENIENT METHODS *
   **********************/


  //! Checks whether another Vector2 instance has the same shape.
  /*! Checks whether another Vector2 instance has the same shape as the
    current instance. The shapes are the same if both instances have the same
    number of inner vectors, and if the inner vectors have the same lengths.
    \param[in] V Vector2 instance whose shape is compared to that of the
    current instance.
    \return True if the current instance as the same shape as \a V, false
    otherwise.
  */
  template <class T, class Allocator0, class Allocator1>
  template <class V2>
  bool Vector2<T, Allocator0, Allocator1>::HasSameShape(const V2& V) const
  {
    if (V.GetLength() != GetLength())
      return false;
    for (int i = 0; i < GetLength(); i++)
      if (V.GetLength(i) != GetLength(i))
	return false;
    return true;
  }


  //! Displays the vector.
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>::Print() const
  {
    for (int i = 0; i < data_.GetLength(); i++)
      {
        cout << "Vector " << i << ": ";
        data_(i).Print();
      }
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the instance in a binary file.
  /*!
    The number of inner vectors (integer) is written first. Then for each
    vector, the length of the vector (integer) and all elements of the vector
    are written.
    \param[in] file_name file name.
    \param[in] with_size if set to 'false', the number of vectors and the
    lengths of the inner vectors are not saved.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Write(string file_name, bool with_size) const
  {
    ofstream file_stream;
    file_stream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!file_stream.is_open())
      throw IOError("Vector2::Write(string file_name, bool with_size)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->Write(file_stream, with_size);

    file_stream.close();
  }


  //! Writes the instance in a stream in a binary format.
  /*!
    The number of inner vectors (integer) is written first. Then for each
    vector, the length of the vector (integer) and all elements of the vector
    are written.
    \param[in,out] stream output stream.
    \param[in] with_size if set to 'false', the number of vectors and the
    lengths of the inner vectors are not saved.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Write(ostream& stream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector2::Write(ostream& stream, bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int m = GetLength();
        stream.write(reinterpret_cast<char*>(const_cast<int*>(&m)),
                     sizeof(int));
      }

    for (int i = 0; i < GetLength(); i++)
      data_(i).Write(stream, with_size);

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!stream.good())
      throw IOError("Vector2::Write(ostream& stream, bool with_size)",
                    "Output operation failed.");
#endif

  }


  //! Reads the Vector2 from a file.
  /*!
    Sets the current Vector2 instance according to a binary file that stores
    the total number of inner vectors, and, for each inner vector, the length
    of the vector (integer) and its elements.
    \param[in] file_name file name.
    \param[in] with_size if set to 'false', the total number of inner vectors
    and the lengths of the vectors are not available in the file. In this
    case, the shape of the current instance is unchanged and the values of the
    elements are directly read in the file.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Read(string file_name, bool with_size)
  {
    ifstream file_stream;
    file_stream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!file_stream.is_open())
      throw IOError("Vector2::Read(string file_name, bool with_size)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->Read(file_stream, with_size);

    file_stream.close();
  }


  //! Reads the Vector2 from a file.
  /*!
    Sets the current Vector2 instance according to a binary stream that stores
    the total number of inner vectors, and, for each inner vector, the length
    of the vector (integer) and its elements.
    \param[in,out] stream input stream.
    \param[in] with_size if set to 'false', the total number of inner vectors
    and the lengths of the vectors are not available in the stream. In this
    case, the shape of the current instance is unchanged and the values of the
    elements are directly read in the file.
  */
  template <class T, class Allocator0, class Allocator1>
  void Vector2<T, Allocator0, Allocator1>
  ::Read(istream& stream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector2::Read(istream& stream, bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int new_size;
        stream.read(reinterpret_cast<char*>(&new_size), sizeof(int));
        this->Reallocate(new_size);
      }

    for (int i = 0; i < GetLength(); i++)
      data_(i).Read(stream, with_size);

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!stream.good())
      throw IOError("Vector2::Read(istream& stream, bool with_size)",
                    "Output operation failed.");
#endif

  }


} // namespace Seldon.


#define SELDON_FILE_VECTOR_VECTOR2_CXX
#endif
