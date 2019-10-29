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


#ifndef SELDON_FILE_SUBMATRIX_HXX


namespace Seldon
{


  ////////////////////////
  // MATRIX<SUBSTORAGE> //
  ////////////////////////


  //! Sub-matrix class.
  template <class T, class Prop, class M, class Allocator>
  class Matrix<T, Prop, SubStorage<M>, Allocator>:
    public SubMatrix_Base<T, Prop, M, Allocator>
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

    // Methods.
  public:
    // Constructor.
    Matrix(M& A, Vector<int>& row_list, Vector<int>& column_list);

    // Destructor.
    ~Matrix();
  };


  ///////////////
  // SUBMATRIX //
  ///////////////


  //! Sub-matrix class.
  /*!
    \extends SubMatrix_Base<T, Prop, M, Allocator>
  */
  template <class M>
  class SubMatrix:
    public Matrix<typename M::value_type, typename M::property,
                  SubStorage<M>, typename M::allocator>
  {
  public:
    SubMatrix(M& A, Vector<int>& row_list, Vector<int>& column_list);
  };


} // namespace Seldon.


#define SELDON_FILE_SUBMATRIX_HXX
#endif
