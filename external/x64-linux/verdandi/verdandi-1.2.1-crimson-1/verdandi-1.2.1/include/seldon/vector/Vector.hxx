// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_VECTOR_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Base structure for all vectors.
  /*!
    It stores data and the vector size. It defines basic
    methods as well.
  */
  template <class T, class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Vector_Base
  {
    // typdef declarations.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    // Static attributes.
  protected:
    static Allocator vect_allocator_;

    // Attributes.
  protected:
    // Number of elements.
    int m_;
    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Vector_Base();
    explicit Vector_Base(int i);
    Vector_Base(const Vector_Base<T, Allocator>& A);

    // Destructor.
    ~Vector_Base();

    void Clear();

    // Basic methods.
    int GetM() const;
    int GetLength() const;
    int GetSize() const;
    pointer GetData() const;
    const_pointer GetDataConst() const;
    void* GetDataVoid() const;
    const void* GetDataConstVoid() const;

  };


  // Vector allocator.
  template <class T, class Allocator>
  Allocator Vector_Base<T, Allocator>::vect_allocator_;


  //! Full vector class.
  /*!
    Basic vector class (i.e. not sparse).
  */
  template <class T, class Allocator>
  class Vector<T, VectFull, Allocator>: public Vector_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    typedef VectFull storage;

    // Attributes.
  private:

    // Methods.
  public:
    // Constructor.
    explicit Vector();
    explicit Vector(int i);
    Vector(int i, pointer data);
    Vector(const Vector<T, VectFull, Allocator>& A);

    // Destructor.
    ~Vector();

    // Memory management.
    void Reallocate(int i);
    void Resize(int i);
    void SetData(int i, pointer data);
    template <class Allocator0>
    void SetData(const Vector<T, VectFull, Allocator0>& V);
    void Nullify();

    // Element access and affectation.
    reference operator() (int i);
    reference Get(int i);
#ifndef SWIG
    const_reference operator() (int i) const;
    const_reference Get(int i) const;
    Vector<T, VectFull, Allocator>& operator= (const Vector<T, VectFull,
					       Allocator>& X);
#endif
    void Copy(const Vector<T, VectFull, Allocator>& X);
    Vector<T, VectFull, Allocator> Copy() const;
    void Append(const T& x);
    template<class T0>
    void PushBack(const T0& x);
    template<class Allocator0>
    void PushBack(const Vector<T, VectFull, Allocator0>& X);

    // Basic functions.
    int GetDataSize();

    // Convenient functions.
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
#ifndef SWIG
    template <class T0>
    Vector<T, VectFull, Allocator>& operator= (const T0& X);
#endif
    template <class T0>
    Vector<T, VectFull, Allocator>& operator*= (const T0& X);
    void FillRand();
    void Print() const;

    // Norms.
    value_type GetNormInf() const;
    int GetNormInfIndex() const;

    // Input/output functions.
    void Write(string FileName, bool with_size = true) const;
    void Write(ostream& FileStream, bool with_size = true) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName, bool with_size = true);
    void Read(istream& FileStream, bool with_size = true);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);
  };

#ifndef SWIG
  template <class T, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Storage, Allocator>& V);
#endif


} // namespace Seldon.

#define SELDON_FILE_VECTOR_HXX
#endif
