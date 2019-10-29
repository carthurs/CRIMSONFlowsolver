// Copyright (C) 2010 Marc Duruflé
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


#ifndef SELDON_FILE_ILUT_PRECONDITIONING_HXX

namespace Seldon
{

  template<class real, class cplx,
           class Allocator = SELDON_DEFAULT_ALLOCATOR<cplx> >
  class IlutPreconditioning
  {
  protected :
    //! Verbosity level.
    int print_level;
    //! True if symmetric matrix is constructed.
    bool symmetric_algorithm;
    //! Type of incomplete factorization.
    int type_ilu;
    /*! \brief Maximum number of elements on a row of L or U.
      For Ilu(k), fill_level = k
    */
    int fill_level;
    /*! \brief Additional number of elements for each row.
      This number is only used for ILUT(k)
    */
    int additional_fill;
    //! Size of block where the pivot is searched.
    int mbloc;
    //! Diagonal compensation parameter (alpha = 0 -> ILU, alpha = 1 -> MILU).
    real alpha;
    //! Threshold used for dropping small terms.
    real droptol;
    //! Threshold for pivoting.
    real permtol;
    //! Permutation arrays.
    IVect permutation_row, permutation_col;
    //! Symmetric matrix.
    Matrix<cplx, Symmetric, ArrayRowSymSparse, Allocator> mat_sym;
    //! Unsymmetric matrix.
    Matrix<cplx, General, ArrayRowSparse, Allocator> mat_unsym;
    //! Temporary vector.
    Vector<cplx, VectFull, Allocator> xtmp;

  public :

    //! Available types of incomplete factorization.
    enum {ILUT, ILU_D, ILUT_K, ILU_0, MILU_0, ILU_K};

    IlutPreconditioning();

    void Clear();

    int GetFactorisationType() const;
    int GetFillLevel() const;
    int GetAdditionalFillNumber() const;
    int GetPrintLevel() const;
    int GetPivotBlockInteger() const;

    void SetFactorisationType(int);
    void SetFillLevel(int);
    void SetAdditionalFillNumber(int);
    void SetPrintLevel(int);
    void SetPivotBlockInteger(int);
    void SetSymmetricAlgorithm();
    void SetUnsymmetricAlgorithm();

    real GetDroppingThreshold() const;
    real GetDiagonalCoefficient() const;
    real GetPivotThreshold() const;

    void SetDroppingThreshold(real);
    void SetDiagonalCoefficient(real);
    void SetPivotThreshold(real);

    template<class MatrixSparse>
    void FactorizeSymMatrix(const IVect& perm,
                            MatrixSparse& mat, bool keep_matrix = false);

    template<class MatrixSparse>
    void FactorizeUnsymMatrix(const IVect& perm,
                              MatrixSparse& mat, bool keep_matrix = false);

    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, General, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class T0, class Storage0, class Allocator0>
    void FactorizeMatrix(const IVect& perm,
                         Matrix<T0, Symmetric, Storage0, Allocator0>& mat,
                         bool keep_matrix = false);

    template<class Matrix1, class Vector1>
    void TransSolve(const Matrix1& A, const Vector1& r, Vector1& z);

    template<class Matrix1, class Vector1>
    void Solve(const Matrix1& A, const Vector1& r, Vector1& z);

    template<class Vector1>
    void TransSolve(Vector1& z);

    template<class Vector1>
    void Solve(Vector1& z);

    template<class TransStatus, class Vector1>
    void Solve(const TransStatus& transA, Vector1& z);

  };

}

#define SELDON_FILE_ILUT_PRECONDITIONING_HXX
#endif
