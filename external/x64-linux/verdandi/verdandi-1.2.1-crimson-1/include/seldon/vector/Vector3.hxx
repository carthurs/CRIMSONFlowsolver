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


#ifndef SELDON_FILE_VECTOR_VECTOR_3_HXX


#ifndef SELDON_VECTOR3_DEFAULT_ALLOCATOR_0
/*! \def SELDON_VECTOR3_DEFAULT_ALLOCATOR_0
  Default allocator for the inner vectors in Vector3.
*/
#define SELDON_VECTOR3_DEFAULT_ALLOCATOR_0 SELDON_DEFAULT_ALLOCATOR
#endif

#ifndef SELDON_VECTOR3_DEFAULT_ALLOCATOR_1
/*! \def SELDON_VECTOR3_DEFAULT_ALLOCATOR_1
  Default allocator for the inner vector of vectors in Vector3.
*/
#define SELDON_VECTOR3_DEFAULT_ALLOCATOR_1 MallocObject
#endif

#ifndef SELDON_VECTOR3_DEFAULT_ALLOCATOR_2
#define SELDON_VECTOR3_DEFAULT_ALLOCATOR_2 MallocObject
/*! \def SELDON_VECTOR3_DEFAULT_ALLOCATOR_2
  Default allocator for the vector of vectors of vectors in Vector3.
*/
#endif


namespace Seldon
{

  //! %Vector of vectors of vectors.
  /*! Vector3 is a structure that acts like a vector of vectors of
    vectors. Both inner vectors and inner vectors of vectors can be of any
    dimension, so that this structure is more flexible than an Array3D.
    \tparam T numerical type of the inner vectors.
    \tparam Allocator0 allocator for the inner vectors. The default allocator
    is SELDON_DEFAULT_ALLOCATOR.
    \tparam Allocator1 allocator for the vector of vectors. It is recommended
    to choose NewAlloc or, for more efficient in reallocations, MallocObject
    (default allocator here): these allocators can manage an array of inner
    vectors.
    \tparam Allocator2 allocator for the vector of vectors of vectors. It is
    recommended to choose NewAlloc or, for more efficient in reallocations,
    MallocObject (default allocator here).
  */
  template <class T,
	    class Allocator0 = SELDON_VECTOR3_DEFAULT_ALLOCATOR_0<T>,
	    class Allocator1 = SELDON_VECTOR3_DEFAULT_ALLOCATOR_1<
	      Vector<T, Vect_Full, Allocator0> >,
	    class Allocator2 =
	    SELDON_VECTOR3_DEFAULT_ALLOCATOR_2<
	      Vector<Vector<T, Vect_Full, Allocator0>,
		     Vect_Full, Allocator1> > >
  class Vector3
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

  protected:
    Vector<Vector<Vector<T, Vect_Full, Allocator0>, Vect_Full, Allocator1>,
	   Vect_Full, Allocator2> data_;

  public:

    /*** Constructors and destructor ***/

    Vector3();
    Vector3(int);
    Vector3(Vector<int>& length);
    template <class Allocator>
    Vector3(Vector<Vector<int>, Vect_Full, Allocator>& length);
    ~Vector3();

    /*** Management of the vectors ***/

    int GetLength() const;
    int GetSize() const;
    int GetLength(int i) const;
    int GetSize(int i) const;
    int GetLength(int i, int j) const;
    int GetSize(int i, int j) const;
    int GetNelement() const;
    int GetNelement(int beg, int end) const;
    Vector<int> GetShape(int i) const;
    void GetShape(int i, Vector<int>& shape) const;
    void Reallocate(int N);
    void Reallocate(int i, int N);
    void Reallocate(int i, int j, int N);

    template <class Td, class Allocatord>
    void Flatten(Vector<Td, VectFull, Allocatord>& data) const;
    template <class Td, class Allocatord>
    void Flatten(int beg, int end, Vector<Td, VectFull, Allocatord>& data)
      const;

    void PushBack(int i, int j, const T& x);
    void PushBack(int i, const Vector<T, Vect_Full, Allocator0>& X);
    void PushBack(const Vector<Vector<T, Vect_Full, Allocator0>,
		  Vect_Full, Allocator1>& X);
    void PushBack(const Vector<Vector<Vector<T, Vect_Full, Allocator0>,
                  Vect_Full, Allocator1>, Vect_Full, Allocator2>& X);
    void PushBack(const Vector3<T, Allocator0, Allocator1, Allocator2>& X);


    void Clear();
    void Clear(int i);
    void Clear(int i, int j);

    void Fill(const T& x);

    Vector<Vector<Vector<T, Vect_Full, Allocator0>, Vect_Full, Allocator1>,
	   Vect_Full, Allocator2>&
    GetVector();
    const Vector<Vector<Vector<T, Vect_Full, Allocator0>,
                        Vect_Full, Allocator1>, Vect_Full, Allocator2>&
    GetVector() const;

    Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    GetVector(int i);
    const Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    GetVector(int i) const;

    Vector<T, Vect_Full, Allocator0>& GetVector(int i, int j);
    const Vector<T, Vect_Full, Allocator0>& GetVector(int i, int j) const;

    /*** Element access and assignment ***/
    const
    Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    operator() (int i) const;
    Vector<Vector<T, Vect_Full, Allocator0>, VectFull, Allocator1>&
    operator() (int i);

    const Vector<T, Vect_Full, Allocator0>& operator() (int i, int j)
      const;
    Vector<T, Vect_Full, Allocator0>& operator() (int i, int j);

    const_reference operator() (int i, int j, int k) const;
    reference operator() (int i, int j, int k);

    /*** Convenient method ***/

    void Print() const;

    /*** Input/output functions ***/

    void Write(string file_name, bool with_size = true) const;
    void Write(ostream& file_stream, bool with_size = true) const;
    void Read(string file_name, bool with_size = true);
    void Read(istream& file_stream, bool with_size = true);
  };

}


#define SELDON_FILE_VECTOR_VECTOR_3_HXX
#endif
