// Copyright (C) 2009 Vivien Mallet
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


#ifndef SELDON_FILE_SUBMATRIX_BASE_HXX


namespace Seldon
{


  //! Sub-matrix base class.
  template <class T, class Prop, class M, class Allocator>
  class SubMatrix_Base: public Matrix_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename M::value_type value_type;
    typedef typename M::pointer pointer;
    typedef typename M::const_pointer const_pointer;
    typedef typename M::reference reference;
    typedef typename M::const_reference const_reference;
    typedef typename M::entry_type entry_type;
    typedef typename M::access_type access_type;
    typedef typename M::const_access_type const_access_type;

    // Attributes.
  protected:
    //! Pointer to the base matrix.
    M* matrix_;
    //! List of rows.
    Vector<int> row_list_;
    //! List of columns.
    Vector<int> column_list_;

    // Methods.
  public:
    // Constructor.
    SubMatrix_Base(M& A, Vector<int>& row_list, Vector<int>& column_list);

    // Destructor.
    ~SubMatrix_Base();

    // Element access and affectation.
    access_type operator() (int i, int j);
#ifndef SWIG
    const_access_type operator() (int i, int j) const;
#endif
    entry_type& Val(int i, int j);
#ifndef SWIG
    const entry_type& Val(int i, int j) const;
#endif

    // Basic methods.
    int GetM() const;
    int GetN() const;
    int GetM(const SeldonTranspose& status) const;
    int GetN(const SeldonTranspose& status) const;
    void Print() const;
  };


} // namespace Seldon.


#define SELDON_FILE_SUBMATRIX_BASE_HXX
#endif
