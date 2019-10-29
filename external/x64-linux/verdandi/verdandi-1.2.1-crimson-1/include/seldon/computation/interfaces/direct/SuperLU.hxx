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


#ifndef SELDON_FILE_SUPERLU_HXX

extern "C"
{
#include "superlu_interface.h"
}


namespace Seldon
{

  //! class interfacing SuperLU functions
  template<class T>
  class MatrixSuperLU_Base
  {
  protected :
    //! objects of SuperLU
    SuperMatrix L, U, B;
    SCformat *Lstore;  //!< object of SuperLU
    NCformat *Ustore;  //!< object of SuperLU
    SuperLUStat_t stat; //!< statistics
    superlu_options_t options; //!< options
    //! permutation array
    Vector<int> perm_r, perm_c;

    colperm_t permc_spec; //!< ordering scheme
    int n; //!< number of rows
    bool display_info; //!< display information about factorization ?
    //! Error code returned by SuperLU.
    int info_facto;

  public :
    MatrixSuperLU_Base();
    ~MatrixSuperLU_Base();

    template<class Prop, class Allocator>
    void GetLU(Matrix<double, Prop, ColSparse, Allocator>& Lmat,
               Matrix<double, Prop, ColSparse, Allocator>& Umat,
               bool permuted = true);

    template<class Prop, class Allocator>
    void GetLU(Matrix<double, Prop, RowSparse, Allocator>& Lmat,
               Matrix<double, Prop, RowSparse, Allocator>& Umat,
               bool permuted = true);

    const Vector<int>& GetRowPermutation() const;
    const Vector<int>& GetColPermutation() const;

    void SelectOrdering(colperm_t type);
    void SetPermutation(const IVect&);

    void Clear();
    void HideMessages();
    void ShowMessages();

    int GetInfoFactorization() const;
  };


  //! empty matrix
  template<class T>
  class MatrixSuperLU : public MatrixSuperLU_Base<T>
  {
  };


  //! class interfacing SuperLU functions in double precision
  template<>
  class MatrixSuperLU<double> : public MatrixSuperLU_Base<double>
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<double>() {}

    template<class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<double, Prop, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<double, VectFull, Allocator2>& x);

    template<class TransStatus, class Allocator2>
    void Solve(const TransStatus& TransA,
               Vector<double, VectFull, Allocator2>& x);
  };


  //! class interfacing SuperLU functions in complex double precision
  template<>
  class MatrixSuperLU<complex<double> >
    : public MatrixSuperLU_Base<complex<double> >
  {
  public:
    MatrixSuperLU() : MatrixSuperLU_Base<complex<double> >() {}

    template<class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<complex<double>, Prop,
			 Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<complex<double>, VectFull, Allocator2>& x);

    template<class TransStatus, class Allocator2>
    void Solve(const TransStatus& TransA,
               Vector<complex<double>, VectFull, Allocator2>& x);

  };

}

#define SELDON_FILE_SUPERLU_HXX
#endif
