// Copyright (C) 2010, INRIA
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


#ifndef SELDON_FILE_VECTOR_VECTOR_3_CXX


#include "Vector3.hxx"


namespace Seldon
{


  /////////////
  // VECTOR3 //
  /////////////


  /***************
   * CONSTRUCTOR *
   ***************/


  //! Default constructor.
  /*!
    Nothing is allocated.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector3<T, Allocator0, Allocator1, Allocator2>::Vector3()
  {
  }


  //! Constructor.
  /*! The vector of vectors of vectors is allocated with \a length empty
    vectors of vectors.
    \param[in] length the length of the vector of vectors of vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector3<T, Allocator0, Allocator1, Allocator2>::Vector3(int length)
  {
    data_.Reallocate(length);
  }


  //! Constructor.
  /*! The vector of vectors of vectors and the inner vectors of vectors are
    allocated.
    \param[in] length the lengths of the inner vectors of vectors. The vector
    of vectors of vectors will obviously have as many elements as \a length
    has.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Vector3(Vector<int>& length)
  {
    data_.Clear();
    int n, m = length.GetSize();
    data_.Reallocate(m);
    for (int i = 0; i < m; i++)
      {
	n = length(i);
	data_(i).Reallocate(n);
      }
  }


  //! Constructor.
  /*! The vector of vectors of vectors, the inner vectors of vectors and the
    inner vectors are allocated.
    \param[in] length the lengths of the inner vectors of vectors and
    the inner vectors. The vector of vectors of vectors will obviously have as
    many elements as \a length has.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  template <class Allocator>
  Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Vector3(Vector<Vector<int>, Vect_Full, Allocator>& length)
  {
    data_.Clear();
    int n, m = length.GetSize();
    data_.Reallocate(m);
    for (int i = 0; i < m; i++)
      {
	n = length(i).GetSize();
	data_(i).Reallocate(n);
	for (int j = 0; j < n; j++)
	  data_(i)(j).Reallocate(length(i)(j));
      }
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  /*! The vector of vectors of vectors, the inner vectors of vectors and the
    inner vectors are deallocated.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector3<T, Allocator0, Allocator1, Allocator2>::~Vector3()
  {
  }


  /**********************
   * VECTORS MANAGEMENT *
   **********************/


  //! Returns size along dimension 1.
  /*!
    \return The size along dimension 1.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>::GetLength() const
  {
    return data_.GetLength();
  }


  //! Returns size along dimension 1.
  /*!
    \return The size along dimension 1.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>::GetSize() const
  {
    return data_.GetSize();
  }


  //! Returns the size of the inner vector of vectors #\a i.
  /*!
    \param[in] i index.
    \return The size of the inner vector of vectors #\a i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>::GetLength(int i) const
  {
    return data_(i).GetLength();
  }


  //! Returns the size of the inner vector of vectors #\a i.
  /*!
    \param[in] i index.
    \return The size of the inner vector of vectors #\a i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>::GetSize(int i) const
  {
    return data_(i).GetSize();
  }


  //! Returns the size of the inner vector #\a i #\a j.
  /*!
    \param[in] i index.
    \param[in] j index.
    \return The  size of the inner vector #\a i #\a j.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>
  ::GetLength(int i, int j) const
  {
    return data_(i)(j).GetLength();
  }


  //! Returns the size of the inner vector #\a i #\a j.
  /*!
    \param[in] i index.
    \param[in] j index.
    \return The  size of the inner vector #\a i #\a j.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>
  ::GetSize(int i, int j) const
  {
    return data_(i)(j).GetSize();
  }


  //! Returns the total number of elements in the inner vectors.
  /*!
    \return The sum of the lengths of the inner vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>::GetNelement() const
  {
    int total = 0;
    for (int i = 0; i < GetLength(); i++)
      for (int j = 0; j < GetLength(i); j++ )
	total += GetLength(i, j);
    return total;
  }


  //! Returns the total number of elements in a range of inner vectors.
  /*! Returns the total number of elements in the range [\a beg, \a end[ of
    inner vectors of vectors.
    \param[in] beg inclusive lower-bound for the indexes.
    \param[in] end exclusive upper-bound for the indexes.
    \return The sum of the lengths of the inner vectors of vectors with
    index \a beg to \a end-1.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  int Vector3<T, Allocator0, Allocator1, Allocator2>
  ::GetNelement(int beg, int end) const
  {
    if (beg > end)
      throw WrongArgument("Vector3::GetNelement(int beg, int end)",
                          "The lower bound of the range of inner vectors "
			  "of vectors, [" + to_str(beg) + ", " + to_str(end)
                          + "[, is strictly greater than its upper bound.");
    if (beg < 0 || end > GetLength())
      throw WrongArgument("Vector3::GetNelement(int beg, int end)",
        		  "The inner-vector of vectors indexes should be in "
			  "[0," + to_str(GetLength()) + "] but ["
			  + to_str(beg)
                          + ", " + to_str(end) + "[ was provided.");

    int total = 0;
    for (int i = beg; i < end; i++)
      for (int j = 0; j < GetLength(i); j++ )
	total += GetLength(i, j);
    return total;
  }


  //! Returns the shape of a vector of vectors.
  /*!
    \param[in] i index of the vector of vectors.
    \return A vector with the lengths of the inner vectors of the \a i th
    vector of vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector<int> Vector3<T, Allocator0, Allocator1, Allocator2>
  ::GetShape(int i) const
  {
    Vector<int> shape(GetLength(i));
    for (int j = 0; j < GetLength(i); j++)
      shape(j) = GetLength(i, j);
    return shape;
  }


  //! Returns the shape of a vector of vectors.
  /*!
    \param[in] i index of the vector of vectors.
    \param[out] shape the lengths of the inner vectors of the \a i th vector
    of vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::GetShape(int i, Vector<int>& shape) const
  {
    shape.Reallocate(GetLength(i));
    for (int j = 0; j < GetLength(i); j++)
      shape(j) = GetLength(i, j);
  }


  //! Reallocates the vector of vectors of vectors.
  /*!
    \param[in] N the new size of the vector of vectors of vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>::Reallocate(int N)
  {
    data_.Reallocate(N);
  }


  //! Reallocates the inner vector of vectors #\a i.
  /*!
    \param[in] i index of the inner vector of vectors to be reallocated.
    \param[in] N the new size of the inner vector of vectors #\a i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Reallocate(int i, int N)
  {
    data_(i).Reallocate(N);
  }


  //! Reallocates the inner vector #\a i #\a j.
  /*!
    \param[in] i index of the inner vector of vectors.
    \param[in] j index of the inner vector to be reallocated.
    \param[in] N the new vector size.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Reallocate(int i, int j, int N)
  {
    data_(i)(j).Reallocate(N);
  }


  //! Returns all values in a vector.
  /*! The output vector \a data contains all inner vectors concatenated in the
    same order as they appear in the current Vector2 instance.
    \param[out] data all values from the current Vector2 instance.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  template <class Td, class Allocatord>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Flatten(Vector<Td, VectFull, Allocatord>& data) const
  {
    data.Reallocate(GetNelement());
    int i, j, k, n(0);
    for (i = 0; i < GetLength(); i++)
      for (j = 0; j < GetLength(i); j++)
        for (k = 0; k < GetLength(i, j); k++)
	  data(n++) = data_(i)(j)(k);
  }


  //! Returns in a vector all values from a range of inner vectors of vectors.
  /*! The output vector \a data contains all inner vectors, in the index range
    [\a beg, \a end[, concatenated in the same order as they appear in the
    current Vector3 instance.
    \param[in] beg inclusive lower-bound for the indexes.
    \param[in] end exclusive upper-bound for the indexes.
    \param[out] data the values contained in the inner vectors [\a beg, \a
    end[.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  template <class Td, class Allocatord>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Flatten(int beg, int end, Vector<Td, VectFull, Allocatord>& data) const
  {
    if (beg > end)
      throw WrongArgument("Vector3:::Flatten(int beg, int end, Vector& data)",
                          "The lower bound of the range of inner vectors "
			  "of vectors, [" + to_str(beg) + ", " + to_str(end)
                          + "[, is strictly greater than its upper bound.");
    if (beg < 0 || end > GetLength())
      throw WrongArgument("Vector3:::Flatten(int beg, int end, Vector& data)",
        		  "The inner-vector of vectors indexes should be in "
			  "[0," + to_str(GetLength()) + "] but ["
			  + to_str(beg)
                          + ", " + to_str(end) + "[ was provided.");

    data.Reallocate(GetNelement());
    int i, j, k, n(0);
    for (i = beg; i < end; i++)
      for (j = 0; j < GetLength(i); j++)
        for (k = 0; k < GetLength(i, j); k++)
	  data(n++) = data_(i)(j)(k);
  }


  //! Appends an element at the end of the inner vector #\a i #\a j.
  /*!
    \param[in] i index of the inner vector of vectors.
    \param[in] j index of the inner vector to which \a x should be appended.
    \param[in] x element to be appended.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0,
	       Allocator1, Allocator2>::PushBack(int i, int j, const T& x)
  {
    data_(i)(j).PushBack(x);
  }


  //! Appends an inner vector at the end of the inner vector of vectors #\a i.
  /*!
    \param[in] i index of the inner vector of vectors.
    \param[in] X inner vector to be appended.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>::
  PushBack(int i, const Vector<T, Vect_Full, Allocator0>& X)
  {
    data_(i).PushBack(X);
  }


  //! Appends an inner vector of vectors at the end of the vector.
  /*!
    \param[in] X inner vector of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::PushBack(const Vector<Vector<T, Vect_Full, Allocator0>,
	     Vect_Full, Allocator1>& X)
  {
    data_.PushBack(X);
  }


  //! Appends a vector of vectors of vectors.
  /*! The inner vectors of vectors of \a X are appended to the current
    instance, in the same order as they appear in \a X.
    \param[in]  X vector of vectors of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::PushBack(const Vector<Vector<Vector<T, Vect_Full, Allocator0>,
             Vect_Full, Allocator1>, Vect_Full, Allocator2>& X)
  {
    for (int i = 0; i < X.GetLength(); i++)
      data_.PushBack(X(i));
  }


  //! Appends a vector of vectors of vectors.
  /*! The inner vectors of vectors of \a X are appended to the current
    instance, in the same order as they appear in \a X.
    \param[in]  X vector of vectors of vectors to be appended.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::PushBack(const Vector3<T, Allocator0, Allocator1, Allocator2>& X)
  {
    for (int i = 0; i < X.GetLength(); i++)
      data_.PushBack(X.GetVector());
  }


  //! Clears the vector.
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>::Clear()
  {
    data_.Clear();
  }


  //! Clears a given inner vector of vectors.
  /*!
    \param[in] i index of the vector of vectors to be cleared.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0,
	       Allocator1, Allocator2>::Clear(int i)
  {
    data_(i).Clear();
  }


  //! Clears a given inner vector.
  /*!
    \param[in] i index of the vector to be cleared.
    \param[in] j index of the vector to be cleared.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>::Clear(int i, int j)
  {
    data_(i)(j).Clear();
  }


  //! Fills the vector with a given value.
  /*!
    \param[in] x value to fill the vector with.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>::Fill(const T& x)
  {
    for (int i = 0; i < data_.GetSize(); i++)
      for (int j = 0; j < data_(i).GetSize(); j++)
	data_(i)(j).Fill(x);
  }


  //! Returns the vector of vectors of vectors.
  /*!
    \return The vector of vectors of vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector<Vector<Vector<T, Vect_Full, Allocator0>, Vect_Full, Allocator1>,
	 Vect_Full, Allocator2>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::GetVector()
  {
    return data_;
  }


  //! Returns the vector of vectors of vectors.
  /*!
    \return The vector of vectors of vectors.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  const Vector<Vector<Vector<T, Vect_Full, Allocator0>,
		      Vect_Full, Allocator1>, Vect_Full, Allocator2>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::GetVector() const
  {
    return data_;
  }


  //! Returns a given inner vector of vectors.
  /*!
    \param[in] i index of the inner vector of vectors.
    \return The inner vector of vectors i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::GetVector(int i)
  {
    return data_(i);
  }


  //! Returns a given inner vector of vectors.
  /*!
    \param[in] i index of the inner vector of vectors.
    \return The inner vector of vectors i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  const Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::GetVector(int i) const
  {
    return data_(i);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the inner vector.
    \return The inner vector #\a i #\a j.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector<T, Vect_Full, Allocator0>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::GetVector(int i, int j)
  {
    return data_(i)(j);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the inner vector.
    \return The inner vector #\a i #\a j.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  const Vector<T, Vect_Full, Allocator0>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::GetVector(int i, int j)
    const
  {
    return data_(i)(j);
  }


  /*********************************
   * ELEMENT ACCESS AND ASSIGNMENT *
   *********************************/


  //! Returns a given inner vector of vectors.
  /*!
    \param[in] i index of the inner vector of vectors.
    \return The inner vector of vectors i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  const Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
  Vector3<T, Allocator0,
	  Allocator1, Allocator2>::operator() (int i) const
  {
    return data_(i);
  }


  //! Returns a given inner vector of vectors.
  /*!
    \param[in] i index of the inner vector of vectors.
    \return The inner vector of vectors i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::operator() (int i)
  {
    return data_(i);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the inner vector.
    \return The inner vector #\a i #\a j.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  const Vector<T, Vect_Full, Allocator0>&
  Vector3<T, Allocator0,
	  Allocator1, Allocator2>::operator() (int i, int j) const
  {
    return data_(i)(j);
  }


  //! Returns a given inner vector.
  /*!
    \param[in] i index of the inner vector.
    \param[in] j index of the inner vector.
    \return The inner vector #\a i #\a j.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  Vector<T, Vect_Full, Allocator0>&
  Vector3<T, Allocator0, Allocator1, Allocator2>::operator() (int i, int j)
  {
    return data_(i)(j);
  }


  //! Returns an element of a given inner vector of a given vector of vectors.
  /*!
    \param[in] i index of the inner vector of vectors.
    \param[in] j index of the inner vector in the vector of vectors #\a i.
    \param[in] k index of the element in the inner vector #\a j.
    \return The element #\a k of vector #\a j of vector of vectors #\a i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  typename Vector3<T, Allocator0, Allocator1, Allocator2>::const_reference
  Vector3<T, Allocator0,
          Allocator1, Allocator2>::operator() (int i, int j, int k) const
  {
    return data_(i)(j)(k);
  }


  //! Returns an element of a given inner vector of a given vector of vectors.
  /*!
    \param[in] i index of the inner vector of vectors.
    \param[in] j index of the inner vector in the vector of vectors #\a i.
    \param[in] k index of the element in the inner vector #\a j.
    \return The element #\a k of vector #\a j of vector of vectors #\a i.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  typename Vector3<T, Allocator0, Allocator1, Allocator2>::reference
  Vector3<T, Allocator0, Allocator1, Allocator2>::operator()
    (int i, int j, int k)
  {
    return data_(i)(j)(k);
  }


  /**********************
   * CONVENIENT METHODS *
   *********************/


  //! Displays the vector.
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>::Print() const
  {
    for (int i = 0; i < data_.GetSize(); i++)
      for(int j = 0; j < data_(i).GetSize(); j++)
	{
	  cout << "Vector " << i << ", " << j << ": ";
	  data_(i)(j).Print();
	}
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/


  //! Writes the instance in a binary file.
  /*!
    \param[in] file_name file name.
    \param[in] with_size if set to 'false', the sizes are not saved so that
    the shape of the instance is lost.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Write(string file_name, bool with_size) const
  {
    ofstream file_stream;
    file_stream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!file_stream.is_open())
      throw IOError("Vector3::Write(string file_name, bool with_size)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->Write(file_stream, with_size);

    file_stream.close();
  }


  //! Writes the instance in a stream in a binary format.
  /*!
    \param[in,out] stream output stream.
    \param[in] with_size if set to 'false', the sizes are not saved so that
    the shape of the instance is lost.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Write(ostream& stream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector3::Write(ostream& stream, bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int m = GetLength();
        stream.write(reinterpret_cast<char*>(const_cast<int*>(&m)),
                     sizeof(int));
      }

    for (int i = 0; i < GetLength(); i++)
      {
        if (with_size)
          {
            int m = GetLength(i);
            stream.write(reinterpret_cast<char*>(const_cast<int*>(&m)),
                         sizeof(int));
          }
        for (int j = 0; j < GetLength(i); j++)
          data_(i)(j).Write(stream, with_size);
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!stream.good())
      throw IOError("Vector3::Write(ostream& stream, bool with_size)",
                    "Output operation failed.");
#endif

  }


  //! Reads the Vector3 from a file.
  /*!
    Sets the current Vector3 instance according to a binary file.
    \param[in] file_name file name.
    \param[in] with_size if set to 'false', the shape is not available in the
    file, the shape of the current instance is thus unchanged and the values
    of the elements are directly read in the file.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Read(string file_name, bool with_size)
  {
    ifstream file_stream;
    file_stream.open(file_name.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!file_stream.is_open())
      throw IOError("Vector3::Read(string file_name, bool with_size)",
		    string("Unable to open file \"") + file_name + "\".");
#endif

    this->Read(file_stream, with_size);

    file_stream.close();
  }


  //! Reads the Vector3 from a stream.
  /*!
    Sets the current Vector3 instance according to a binary stream.
    \param[in,out] stream input stream.
    \param[in] with_size if set to 'false', the shape is not available in the
    stream, the shape of the current instance is thus unchanged and the values
    of the elements are directly read in the stream.
  */
  template <class T, class Allocator0, class Allocator1, class Allocator2>
  void Vector3<T, Allocator0, Allocator1, Allocator2>
  ::Read(istream& stream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!stream.good())
      throw IOError("Vector3::Read(istream& stream, bool with_size)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int new_size;
        stream.read(reinterpret_cast<char*>(&new_size), sizeof(int));
        this->Reallocate(new_size);
      }

    for (int i = 0; i < GetLength(); i++)
      {
        if (with_size)
          {
            int new_size;
            stream.read(reinterpret_cast<char*>(&new_size), sizeof(int));
            this->Reallocate(i, new_size);
          }
        for (int j = 0; j < GetLength(i); j++)
          data_(i)(j).Read(stream, with_size);
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!stream.good())
      throw IOError("Vector3::Read(istream& stream, bool with_size)",
                    "Output operation failed.");
#endif

  }


} // namespace Seldon.


#define SELDON_FILE_VECTOR_VECTOR_3_CXX
#endif
