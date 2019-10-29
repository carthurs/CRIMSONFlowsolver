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

#ifndef SELDON_FILE_MATRIX_BASE_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Base class for all matrices.
  /*!
    It stores some data and matrix dimensions. It defines basic
    methods as well.
  */
  template <class T, class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix_Base
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
    static Allocator allocator_;

    // Attributes.
  protected:
    // Number of rows.
    int m_;
    // Number of columns.
    int n_;
    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Matrix_Base();
    explicit Matrix_Base(int i, int j);
    Matrix_Base(const Matrix_Base<T, Allocator>& A);

    // Destructor.
    ~Matrix_Base();

    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetM(const Seldon::SeldonTranspose& status) const;
    int GetN(const Seldon::SeldonTranspose& status) const;
#ifdef SELDON_WITH_BLAS
    int GetM(const CBLAS_TRANSPOSE& status) const;
    int GetN(const CBLAS_TRANSPOSE& status) const;
#endif
    int GetSize() const;
    pointer GetData() const;
    const_pointer GetDataConst() const;
    void* GetDataVoid() const;
    const void* GetDataConstVoid() const;

    Allocator& GetAllocator();

  };


  // Matrix allocator.
  template <class T, class Allocator>
  Allocator Matrix_Base<T, Allocator>::allocator_;


  template <class T, class Prop, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Matrix<T, Prop, Storage, Allocator>& A);


} // namespace Seldon.

#define SELDON_FILE_MATRIX_BASE_HXX
#endif
