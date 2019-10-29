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


// To be included by Seldon.hxx

#ifndef SELDON_FILE_MATRIX_SYMSPARSE_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Symmetric sparse-matrix class.
  /*!
    Symmetric sparse matrices are defined by: (1) the number of rows
    and columns; (2) the number of non-zero entries; (3) an array 'ptr_'
    of start indices (i.e. indices of the first element of each row or column,
    depending on the storage); (4) an array 'ind_' of column or row indices
    of each non-zero entry; (5) values of non-zero entries.\par
    Only values of the upper part are stored.
  */
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_SymSparse: public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef typename Allocator::value_type entry_type;
    typedef typename Allocator::value_type access_type;
    typedef typename Allocator::value_type const_access_type;

    // Attributes.
  protected:
    // Number of non-zero (stored) elements.
    int nz_;
    // Index (in data_) of first element stored for each row or column.
    int* ptr_;
    // Column or row index (in the matrix) each element.
    int* ind_;

    // Methods.
  public:
    // Constructors.
    Matrix_SymSparse();
    Matrix_SymSparse(int i, int j);
    Matrix_SymSparse(int i, int j, int nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix_SymSparse(int i, int j, Vector<T, Storage0, Allocator0>& values,
		     Vector<int, Storage1, Allocator1>& ptr,
		     Vector<int, Storage2, Allocator2>& ind);
    Matrix_SymSparse(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A);

    // Destructor.
    ~Matrix_SymSparse();
    void Clear();

    // Memory management.
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    void SetData(int i, int j,
		 Vector<T, Storage0, Allocator0>& values,
		 Vector<int, Storage1, Allocator1>& ptr,
		 Vector<int, Storage2, Allocator2>& ind);
    void SetData(int i, int j, int nz, pointer values, int* ptr, int* ind);
    void Nullify();
    void Reallocate(int i, int j);
    void Reallocate(int i, int j, int nz);
    void Resize(int i, int j);
    void Resize(int i, int j, int nz);
    void Copy(const Matrix_SymSparse<T, Prop, Storage, Allocator>& A);

    // Basic methods.
    int GetNonZeros() const;
    int GetDataSize() const;
    int* GetPtr() const;
    int* GetInd() const;
    int GetPtrSize() const;
    int GetIndSize() const;

    // Element acess and affectation.
    value_type operator() (int i, int j) const;
    value_type& Val(int i, int j);
    value_type& Get(int i, int j);
    const value_type& Val(int i, int j) const;
    const value_type& Get(int i, int j) const;
    void Set(int i, int j, const T& x);
    void AddInteraction(int i, int j, const T& x);
    Matrix_SymSparse<T, Prop, Storage, Allocator>&
    operator= (const Matrix_SymSparse<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    void FillRand();

    void Print() const;
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName) const;
    void WriteText(ostream& FileStream) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName);
    void ReadText(istream& FileStream);
  };


  //! Column-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymSparse, Allocator>:
    public Matrix_SymSparse<T, Prop, ColSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(int i, int j, int nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<int, Storage1, Allocator1>& ptr,
	   Vector<int, Storage2, Allocator2>& ind);
  };


  //! Row-major sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymSparse, Allocator>:
    public Matrix_SymSparse<T, Prop, RowSymSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowSymSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(int i, int j, int nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<int, Storage1, Allocator1>& ptr,
	   Vector<int, Storage2, Allocator2>& ind);
  };


} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMSPARSE_HXX
#endif
