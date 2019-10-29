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


// To be included by Seldon.hxx

#ifndef SELDON_FILE_MATRIX_TRIANGPACKED_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Triangular packed matrix class.
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_TriangPacked: public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::reference access_type;
    typedef typename Allocator::const_reference const_access_type;

    // Attributes.
  protected:

    // Methods.
  public:
    // Constructor.
    Matrix_TriangPacked();
    Matrix_TriangPacked(int i, int j = 0);
    Matrix_TriangPacked(const Matrix_TriangPacked<T, Prop, Storage,
			Allocator>& A);

    // Destructor.
    ~Matrix_TriangPacked();
    void Clear();

    // Basic methods.
    int GetDataSize() const;

    // Memory management.
    void Reallocate(int i, int j);
    void SetData(int i, int j, pointer data);
    void Nullify();

    // Element access and affectation.
    value_type operator() (int i, int j) const;
    reference Val(int i, int j);
    const_reference Val(int i, int j) const;
    reference Get(int i, int j);
    const_reference Get(int i, int j) const;
    reference operator[] (int i);
    const_reference operator[] (int i) const;
    Matrix_TriangPacked<T, Prop, Storage, Allocator>&
    operator= (const Matrix_TriangPacked<T, Prop, Storage, Allocator>& A);
    void Set(int i, int j, const T& x);
    void Copy(const Matrix_TriangPacked<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
    Matrix_TriangPacked<T, Prop, Storage, Allocator>& operator= (const T0& x);
    void FillRand();
    void Print() const;
    void Print(int a, int b, int m, int n) const;
    void Print(int l) const;

    // Input/output functions.
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);

  };


  //! Column-major upper-triangular packed matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColUpTriangPacked, Allocator>:
    public Matrix_TriangPacked<T, Prop, ColUpTriangPacked, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColUpTriangPacked storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j = 0);
    void Resize(int i, int j);

    template <class T0>
    Matrix<T, Prop, ColUpTriangPacked, Allocator>& operator= (const T0& x);
    Matrix<T, Prop, ColUpTriangPacked, Allocator>&
    operator= (const Matrix<T, Prop, ColUpTriangPacked, Allocator>& A);
    template <class T0>
    Matrix<T, Prop, ColUpTriangPacked, Allocator>& operator*= (const T0& x);
  };


  //! Column-major lower-triangular packed matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColLoTriangPacked, Allocator>:
    public Matrix_TriangPacked<T, Prop, ColLoTriangPacked, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColLoTriangPacked storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j = 0);
    void Resize(int i, int j);

    template <class T0>
    Matrix<T, Prop, ColLoTriangPacked, Allocator>& operator= (const T0& x);
    Matrix<T, Prop, ColLoTriangPacked, Allocator>&
    operator= (const Matrix<T, Prop, ColLoTriangPacked, Allocator>& A);
    template<class T0>
    Matrix<T, Prop, ColLoTriangPacked, Allocator>& operator*= (const T0& x);
  };


  //! Row-major upper-triangular packed matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowUpTriangPacked, Allocator>:
    public Matrix_TriangPacked<T, Prop, RowUpTriangPacked, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowUpTriangPacked storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j = 0);
    void Resize(int i, int j);

    template <class T0>
    Matrix<T, Prop, RowUpTriangPacked, Allocator>& operator= (const T0& x);
    Matrix<T, Prop, RowUpTriangPacked, Allocator>&
    operator= (const Matrix<T, Prop, RowUpTriangPacked, Allocator>& A);
    template<class T0>
    Matrix<T, Prop, RowUpTriangPacked, Allocator>& operator*= (const T0& x);
  };


  //! Row-major lower-triangular packed matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowLoTriangPacked, Allocator>:
    public Matrix_TriangPacked<T, Prop, RowLoTriangPacked, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowLoTriangPacked storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j = 0);
    void Resize(int i, int j);

    template <class T0>
    Matrix<T, Prop, RowLoTriangPacked, Allocator>& operator= (const T0& x);
    Matrix<T, Prop, RowLoTriangPacked, Allocator>&
    operator= (const Matrix<T, Prop, RowLoTriangPacked, Allocator>& A);
    template<class T0>
    Matrix<T, Prop, RowLoTriangPacked, Allocator>& operator*= (const T0& x);
  };


} // namespace Seldon.

#define SELDON_FILE_MATRIX_TRIANGPACKED_HXX
#endif
