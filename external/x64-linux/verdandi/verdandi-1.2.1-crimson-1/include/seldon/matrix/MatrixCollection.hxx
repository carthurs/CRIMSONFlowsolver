// Copyright (C) 2010 INRIA
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


#ifndef SELDON_FILE_MATRIX_COLLECTION_HXX


#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Matrix class made of a collection of matrices.
  template <class T, class Prop, class Storage,
            class Allocator = NewAlloc<T> >
  class MatrixCollection: public Matrix_Base<T, Allocator>
  {
    // Typedef declaration.
  public:
    typedef typename T::value_type value_type;
    typedef typename T::pointer pointer;
    typedef typename T::const_pointer const_pointer;
    typedef typename T::reference reference;
    typedef typename T::const_reference const_reference;

    typedef T matrix_type;
    typedef matrix_type* matrix_pointer;
    typedef const matrix_type* const_matrix_pointer;
    typedef matrix_type& matrix_reference;
    typedef const matrix_type& const_matrix_reference;

    typedef Matrix<matrix_type, Prop, Storage, NewAlloc<T> > collection_type;
    typedef const collection_type const_collection_type;
    typedef collection_type& collection_reference;
    typedef const collection_type& const_collection_reference;


    // Attributes.
  protected:
    //! Number of non-zero elements.
    int nz_;
    //! Number of rows of matrices.
    int Mmatrix_;
    //! Number of columns of matrices.
    int Nmatrix_;
    //! Number of rows in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Mlocal_;
    //! Cumulative number of rows in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Mlocal_sum_;
    //! Number of columns in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Nlocal_;
    //! Cumulative number of columns in the underlying matrices.
    Vector<int, VectFull, CallocAlloc<int> > Nlocal_sum_;

    //! Pointers of the underlying matrices.
    collection_type matrix_;


    // Methods.
  public:
    // Constructor.
    MatrixCollection();
    MatrixCollection(int i, int j);
    MatrixCollection(const MatrixCollection<T, Prop, Storage, Allocator>& A);

    // Destructor.
    ~MatrixCollection();
    void Clear();
    void Nullify();
    void Nullify(int i, int j);

    void Deallocate();

    // Basic methods.
    int GetM() const;
    int GetMmatrix() const;
    int GetM(int i) const;
    int GetN() const;
    int GetNmatrix() const;
    int GetN(int j) const;
    int GetSize() const;
    int GetDataSize() const;

    // Memory management.
    void Reallocate(int i, int j);

    // Management of the matrices.
    template <class T0, class Prop0, class Storage0, class Allocator0>
    void SetMatrix(int m, int n,
                   const Matrix<T0, Prop0, Storage0, Allocator0>&);
    template <class T0, class Prop0, class Allocator0>
    void SetMatrix(int m, int n,
                   const Matrix<T0, Prop0, RowSparse, Allocator0>&);

    // Element access and affectation.
    matrix_reference GetMatrix(int i, int j);
    const_matrix_reference GetMatrix(int i, int j) const;

    value_type operator() (int i, int j) const;

    MatrixCollection<T, Prop, Storage, Allocator>&
    operator= (const MatrixCollection<T, Prop, Storage, Allocator>& A);
    void Copy(const MatrixCollection<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Print() const;
    void Print(int m, int n) const;

    // Input/output functions.
    void Write(string FileName, bool with_size) const;
    void Write(ostream& FileStream, bool with_size) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;

    void Read(string FileName);
    void Read(istream& FileStream);
  };


  //! Column-major matrix collection class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColMajorCollection, Allocator>:
    public MatrixCollection<T, Prop, ColMajor, Allocator>
  {
    // typedef declaration.
  public:
    typedef T value_type;
    typedef Prop property;
    typedef ColMajorCollection storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
  };


  //! Row-major matrix collection class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowMajorCollection, Allocator>:
    public MatrixCollection<T, Prop, RowMajor, Allocator>
  {
    // typedef declaration.
  public:
    typedef T value_type;
    typedef Prop property;
    typedef RowMajorCollection storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
  };


  //! Symetric column-major matrix collection class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymPackedCollection, Allocator>:
    public MatrixCollection<T, Prop, ColSymPacked, Allocator>
  {
    // typedef declaration.
  public:
    typedef T value_type;
    typedef Prop property;
    typedef ColSymPackedCollection storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
  };


  //! Symetric row-major matrix collection class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymPackedCollection, Allocator>:
    public MatrixCollection<T, Prop, RowSymPacked, Allocator>
  {
    // typedef declaration.
  public:
    typedef T value_type;
    typedef Prop property;
    typedef RowSymPackedCollection storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
  };



} // namespace Seldon.


#define SELDON_FILE_MATRIX_COLLECTION_HXX
#endif
