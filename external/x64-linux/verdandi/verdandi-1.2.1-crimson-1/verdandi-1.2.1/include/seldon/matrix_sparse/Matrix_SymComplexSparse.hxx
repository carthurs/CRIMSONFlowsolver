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

#ifndef SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Symmetric complex sparse-matrix class.
  /*!
    Symmetric sparse matrices are defined by: (1) the number of rows
    and columns; (2) the number of non-zero entries; (3) an array 'ptr_' of
    start indices (i.e. indices of the first element of each row or column,
    depending on the storage); (4) an array 'ind_' of column or row indices
    of each non-zero entry; (5) values of non-zero entries.\par
    Complex sparse matrices are defined in the same way except that real
    and imaginary parts are splitted. It is as if two matrices were stored.
    There are therefore 6 arrays: 'real_ptr_', 'real_ind_', 'real_data_',
    'imag_ptr_', 'imag_ind_' and 'imag_data_'.\par
    Finally, since the matrix is symmetric, only its upper part is stored.
  */
  template <class T, class Prop, class Storage,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_SymComplexSparse: public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;
    typedef complex<value_type> entry_type;
    typedef complex<value_type> access_type;
    typedef complex<value_type> const_access_type;

    // Attributes.
  protected:
    // Number of non-zero (stored) elements.
    int real_nz_;
    int imag_nz_;
    // Index (in data_) of first element stored for each row or column.
    int* real_ptr_;
    int* imag_ptr_;
    // Column or row index (in the matrix) each element.
    int* real_ind_;
    int* imag_ind_;

    // Data.
    T* real_data_;
    T* imag_data_;

    // Methods.
  public:
    // Constructors.
    Matrix_SymComplexSparse();
    Matrix_SymComplexSparse(int i, int j);
    Matrix_SymComplexSparse(int i, int j, int real_nz, int imag_nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix_SymComplexSparse(int i, int j,
			    Vector<T, Storage0, Allocator0>& real_values,
			    Vector<int, Storage1, Allocator1>& real_ptr,
			    Vector<int, Storage2, Allocator2>& real_ind,
			    Vector<T, Storage0, Allocator0>& imag_values,
			    Vector<int, Storage1, Allocator1>& imag_ptr,
			    Vector<int, Storage2, Allocator2>& imag_ind);
    Matrix_SymComplexSparse(const Matrix_SymComplexSparse<T, Prop, Storage,
			    Allocator>& A);

    // Destructor.
    ~Matrix_SymComplexSparse();
    void Clear();

    // Memory management.
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    void SetData(int i, int j,
		 Vector<T, Storage0, Allocator0>& real_values,
		 Vector<int, Storage1, Allocator1>& real_ptr,
		 Vector<int, Storage2, Allocator2>& real_ind,
		 Vector<T, Storage0, Allocator0>& imag_values,
		 Vector<int, Storage1, Allocator1>& imag_ptr,
		 Vector<int, Storage2, Allocator2>& imag_ind);
    void SetData(int i, int j,
		 int real_nz, pointer real_values, int* real_ptr,
		 int* real_ind,
		 int imag_nz, pointer imag_values, int* imag_ptr,
		 int* imag_ind);
    void Nullify();
    void Reallocate(int i, int j);
    void Reallocate(int i, int j, int real_nz, int imag_nz);
    void Resize(int i, int j);
    void Resize(int i, int j, int real_nz, int imag_nz);
    void Copy(const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A);

    // Basic methods.
    int GetDataSize() const;
    int* GetRealPtr() const;
    int* GetImagPtr() const;
    int* GetRealInd() const;
    int* GetImagInd() const;
    int GetRealPtrSize() const;
    int GetImagPtrSize() const;
    int GetRealIndSize() const;
    int GetImagIndSize() const;
    int GetRealDataSize() const;
    int GetImagDataSize() const;
    T* GetRealData() const;
    T* GetImagData() const;

    // Element acess and affectation.
    complex<value_type> operator() (int i, int j) const;
    value_type& ValReal(int i, int j);
    const value_type& ValReal(int i, int j) const;
    value_type& ValImag(int i, int j);
    const value_type& ValImag(int i, int j) const;
    value_type& GetReal(int i, int j);
    const value_type& GetReal(int i, int j) const;
    value_type& GetImag(int i, int j);
    const value_type& GetImag(int i, int j) const;
    void Set(int i, int j, const complex<T>& x);
    void AddInteraction(int i, int j, const complex<T>& x);
    Matrix_SymComplexSparse<T, Prop, Storage, Allocator>&
    operator= (const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A);

    // Convenient functions.
    void Zero();
    void SetIdentity();
    void Fill();
    void Fill(const complex<T>& x);
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


  //! Column-major complex sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymComplexSparse, Allocator>:
    public Matrix_SymComplexSparse<T, Prop, ColSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef ColSymComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(int i, int j, int real_nz, int imag_nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<T, Storage0, Allocator0>& real_values,
	   Vector<int, Storage1, Allocator1>& real_ptr,
	   Vector<int, Storage2, Allocator2>& real_ind,
	   Vector<T, Storage0, Allocator0>& imag_values,
	   Vector<int, Storage1, Allocator1>& imag_ptr,
	   Vector<int, Storage2, Allocator2>& imag_ind);
  };


  //! Row-major complex sparse-matrix class.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymComplexSparse, Allocator>:
    public Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef Prop property;
    typedef RowSymComplexSparse storage;
    typedef Allocator allocator;

  public:
    Matrix();
    Matrix(int i, int j);
    Matrix(int i, int j, int real_nz, int imag_nz);
    template <class Storage0, class Allocator0,
	      class Storage1, class Allocator1,
	      class Storage2, class Allocator2>
    Matrix(int i, int j,
	   Vector<T, Storage0, Allocator0>& values,
	   Vector<int, Storage1, Allocator1>& ptr,
	   Vector<int, Storage2, Allocator2>& ind,
	   Vector<T, Storage0, Allocator0>& imag_values,
	   Vector<int, Storage1, Allocator1>& imag_ptr,
	   Vector<int, Storage2, Allocator2>& imag_ind);
  };


} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_HXX
#endif
