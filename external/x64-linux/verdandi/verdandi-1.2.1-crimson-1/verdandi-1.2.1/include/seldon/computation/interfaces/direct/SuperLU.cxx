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


#ifndef SELDON_FILE_SUPERLU_CXX

#include "SuperLU.hxx"

// The function comes from the Matlab interface to SuperLU. It is part of
// SuperLU package. Its copyright is held by University of California
// Berkeley, Xerox Palo Alto Research Center and Lawrence Berkeley National
// Lab. It is released under a license compatible with the GNU LGPL.
void LUextract(SuperMatrix *L, SuperMatrix *U, double *Lval, int *Lrow,
               int *Lcol, double *Uval, int *Urow, int *Ucol, int *snnzL,
               int *snnzU)
{
  int         i, j, k;
  int         upper;
  int         fsupc, istart, nsupr;
  int         lastl = 0, lastu = 0;
  SCformat    *Lstore;
  NCformat    *Ustore;
  double      *SNptr;

  Lstore = static_cast<SCformat*>(L->Store);
  Ustore = static_cast<NCformat*>(U->Store);
  Lcol[0] = 0;
  Ucol[0] = 0;

  /* for each supernode */
  for (k = 0; k <= Lstore->nsuper; ++k) {

    fsupc = L_FST_SUPC(k);
    istart = L_SUB_START(fsupc);
    nsupr = L_SUB_START(fsupc+1) - istart;
    upper = 1;

    /* for each column in the supernode */
    for (j = fsupc; j < L_FST_SUPC(k+1); ++j) {
      SNptr = &(static_cast<double*>(Lstore->nzval))[L_NZ_START(j)];

      /* Extract U */
      for (i = U_NZ_START(j); i < U_NZ_START(j+1); ++i) {
        Uval[lastu] = (static_cast<double*>(Ustore->nzval))[i];
        Urow[lastu++] = U_SUB(i);
      }
      for (i = 0; i < upper; ++i) { /* upper triangle in the supernode */
        Uval[lastu] = SNptr[i];
        Urow[lastu++] = L_SUB(istart+i);
      }
      Ucol[j+1] = lastu;

      /* Extract L */
      Lval[lastl] = 1.0; /* unit diagonal */
      Lrow[lastl++] = L_SUB(istart + upper - 1);
      for (i = upper; i < nsupr; ++i) {
        Lval[lastl] = SNptr[i];
        Lrow[lastl++] = L_SUB(istart+i);
      }
      Lcol[j+1] = lastl;

      ++upper;

    } /* for j ... */

  } /* for k ... */

  *snnzL = lastl;
  *snnzU = lastu;
}


namespace Seldon
{
  //! default constructor
  template<class T>
  MatrixSuperLU_Base<T>::MatrixSuperLU_Base()
  {
    n = 0;
    permc_spec = COLAMD;
    Lstore = NULL;
    Ustore = NULL;
    StatInit(&stat);
    set_default_options(&options);
    ShowMessages();
    display_info = false;
    info_facto = 0;
  }


  //! destructor
  template<class T>
  MatrixSuperLU_Base<T>::~MatrixSuperLU_Base()
  {
    Clear();
  }


  //! Returns the LU factorization.
  /*!
    \param[out] Lmat matrix L in the LU factorization.
    \param[out] Umat matrix U in the LU factorization.
    \param[in] permuted should the permuted matrices be provided? SuperLU
    permutes the rows and columns of the factorized matrix. If \a permuted is
    set to true, L and U are returned as SuperLU computed them, hence with
    permuted rows and columns. If \a permuted is set to false, the matrices L
    and U are "unpermuted" so that L times U is equal to the initial matrix.
  */
  template<class T>
  template<class Prop, class Allocator>
  void MatrixSuperLU_Base<T>
  ::GetLU(Matrix<double, Prop, ColSparse, Allocator>& Lmat,
          Matrix<double, Prop, ColSparse, Allocator>& Umat,
          bool permuted)
  {
    Lstore = static_cast<SCformat*>(L.Store);
    Ustore = static_cast<NCformat*>(U.Store);

    int Lnnz = Lstore->nnz;
    int Unnz = Ustore->nnz;

    int m = U.nrow;
    int n = U.ncol;

    Vector<double, VectFull, Allocator> Lval(Lnnz);
    Vector<int, VectFull, CallocAlloc<int> > Lrow(Lnnz);
    Vector<int, VectFull, CallocAlloc<int> > Lcol(n + 1);

    Vector<double, VectFull, Allocator> Uval(Unnz);
    Vector<int, VectFull, CallocAlloc<int> > Urow(Unnz);
    Vector<int, VectFull, CallocAlloc<int> > Ucol(n + 1);

    int Lsnnz;
    int Usnnz;
    LUextract(&L, &U, Lval.GetData(), Lrow.GetData(), Lcol.GetData(),
              Uval.GetData(), Urow.GetData(), Ucol.GetData(), &Lsnnz, &Usnnz);

    Lmat.SetData(m, n, Lval, Lcol, Lrow);
    Umat.SetData(m, n, Uval, Ucol, Urow);

    if (!permuted)
      {
        Vector<int> row_perm_orig = perm_r;
        Vector<int> col_perm_orig = perm_c;

        Vector<int> row_perm(n);
        Vector<int> col_perm(n);
        row_perm.Fill();
        col_perm.Fill();

        Sort(row_perm_orig, row_perm);
        Sort(col_perm_orig, col_perm);

        ApplyInversePermutation(Lmat, row_perm, col_perm);
        ApplyInversePermutation(Umat, row_perm, col_perm);
      }
  }


  //! Returns the LU factorization.
  /*!
    \param[out] Lmat matrix L in the LU factorization.
    \param[out] Umat matrix U in the LU factorization.
    \param[in] permuted should the permuted matrices be provided? SuperLU
    permutes the rows and columns of the factorized matrix. If \a permuted is
    set to true, L and U are returned as SuperLU computed them, hence with
    permuted rows and columns. If \a permuted is set to false, the matrices L
    and U are "unpermuted" so that L times U is equal to the initial matrix.
    \note This method will first retrieve the L and U matrices in 'ColSparse'
    format and then convert them into 'RowSparse'.
  */
  template<class T>
  template<class Prop, class Allocator>
  void MatrixSuperLU_Base<T>
  ::GetLU(Matrix<double, Prop, RowSparse, Allocator>& Lmat,
          Matrix<double, Prop, RowSparse, Allocator>& Umat,
          bool permuted)
  {
    Lmat.Clear();
    Umat.Clear();

    Matrix<double, Prop, ColSparse, Allocator> Lmat_col;
    Matrix<double, Prop, ColSparse, Allocator> Umat_col;
    GetLU(Lmat_col, Umat_col, permuted);

    Copy(Lmat_col, Lmat);
    Lmat_col.Clear();
    Copy(Umat_col, Umat);
    Umat_col.Clear();
  }


  //! Returns the permutation of rows.
  /*! In order to retain the sparsity as much as possible, SuperLU permutes
    rows and columns before the factorization. This method returns the row
    permutation that was employed in the factorization. This method is
    obviously to be called after the factorization has been performed.
    \return The permutation of the rows.
  */
  template<class T>
  const Vector<int>& MatrixSuperLU_Base<T>::GetRowPermutation() const
  {
    return perm_r;
  }


  //! Returns the permutation of columns.
  /*! In order to retain the sparsity as much as possible, SuperLU permutes
    rows and columns before the factorization. This method returns the column
    permutation that was employed in the factorization. This method is
    obviously to be called after the factorization has been performed.
    \return The permutation of the columns.
  */
  template<class T>
  const Vector<int>& MatrixSuperLU_Base<T>::GetColPermutation() const
  {
    return perm_c;
  }


  template<class T>
  void MatrixSuperLU_Base<T>::SelectOrdering(colperm_t type)
  {
    permc_spec = type;
  }


  template<class T>
  void MatrixSuperLU_Base<T>::SetPermutation(const IVect& permut)
  {
    permc_spec = MY_PERMC;
    perm_c = permut;
    perm_r = permut;
  }


  //! same effect as a call to the destructor
  template<class T>
  void MatrixSuperLU_Base<T>::Clear()
  {
    if (n > 0)
      {
	// SuperLU objects are cleared
        Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	perm_r.Clear();
        perm_c.Clear();
	n = 0;
      }
  }


  //! no message from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::HideMessages()
  {
    display_info = false;
  }


  //! allows messages from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::ShowMessages()
  {
    display_info = true;
  }


  template<class T>
  int MatrixSuperLU_Base<T>::GetInfoFactorization() const
  {
    return info_facto;
  }


  //! factorization of matrix in double precision using SuperLU
  template<class Prop, class Storage, class Allocator>
  void MatrixSuperLU<double>::
  FactorizeMatrix(Matrix<double, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    // clearing previous factorization
    Clear();

    // conversion in CSC format
    n = mat.GetN();
    Matrix<double, General, ColSparse> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();

    SuperMatrix A, AA;
    int nnz = Acsr.GetDataSize();
    dCreate_CompCol_Matrix(&AA, n, n, nnz, Acsr.GetData(), Acsr.GetInd(),
			   Acsr.GetPtr(), SLU_NC, SLU_D, SLU_GE);

    // we get renumbering vectors perm_r and perm_c
    options.ColPerm = permc_spec;
    if (permc_spec != MY_PERMC)
      {
        perm_r.Reallocate(n);
        perm_c.Reallocate(n);
        perm_r.Fill();
        perm_c.Fill();

        get_perm_c(permc_spec, &AA, perm_c.GetData());
      }

    // Original matrix AA is permuted to obtain matrix A.
    Vector<int> etree(n);
    sp_preorder(&options, &AA, perm_c.GetData(), etree.GetData(), &A);

    int panel_size = sp_ienv(1);
    int relax = sp_ienv(2);
    int lwork = 0;

    // Then calling factorization on permuted matrix.
    dgstrf(&options, &A, relax, panel_size, etree.GetData(),
           NULL, lwork, perm_c.GetData(), perm_r.GetData(), &L, &U, &stat,
           &info_facto);

    // Clearing matrices.
    Destroy_CompCol_Permuted(&A);
    Destroy_CompCol_Matrix(&AA);

    if (info_facto == 0 && display_info)
      {
	mem_usage_t mem_usage;
	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	cout << "Number of nonzeros in factor L = " << Lstore->nnz << endl;
	cout << "Number of nonzeros in factor U = " << Ustore->nnz << endl;
	cout << "Number of nonzeros in L+U     = "
             << Lstore->nnz + Ustore->nnz << endl;
	dQuerySpace(&L, &U, &mem_usage);
	cout << "Memory used for factorization in MB: "
             << mem_usage.total_needed / (1024. * 1024.) << endl;
      }

    Acsr.Nullify();
  }


  //! resolution of linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Vector<double, VectFull, Allocator2>& x)
  {
    trans_t trans = NOTRANS;
    int nb_rhs = 1, info;
    // Putting right hand side on SuperLU structure.
    dCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 x.GetData(), x.GetM(), SLU_DN, SLU_D, SLU_GE);

    SuperLUStat_t stat;
    StatInit(&stat);
    // Solving A x = b.
    dgstrs(trans, &L, &U, perm_r.GetData(),
	   perm_c.GetData(), &B, &stat, &info);
  }


  //! resolution of linear system A x = b
  template<class TransStatus, class Allocator2>
  void MatrixSuperLU<double>::Solve(const TransStatus& TransA,
                                    Vector<double, VectFull, Allocator2>& x)
  {
    if (TransA.NoTrans())
      {
        Solve(x);
        return;
      }

    trans_t trans = TRANS;
    int nb_rhs = 1, info;
    // Putting right hand side on SuperLU structure.
    dCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 x.GetData(), x.GetM(), SLU_DN, SLU_D, SLU_GE);

    SuperLUStat_t stat;
    StatInit(&stat);
    // Solving A^T x = b.
    dgstrs(trans, &L, &U, perm_r.GetData(),
	   perm_c.GetData(), &B, &stat, &info);
  }


  //! factorization of matrix in complex double precision using SuperLU
  template<class Prop, class Storage, class Allocator>
  void MatrixSuperLU<complex<double> >::
  FactorizeMatrix(Matrix<complex<double>, Prop, Storage, Allocator> & mat,
		  bool keep_matrix)
  {
    // clearing previous factorization
    Clear();

    // conversion in CSR format
    n = mat.GetN();
    Matrix<complex<double>, General, ColSparse> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();

    SuperMatrix AA, A;
    int nnz = Acsr.GetDataSize();
    zCreate_CompCol_Matrix(&AA, n, n, nnz,
			   reinterpret_cast<doublecomplex*>(Acsr.GetData()),
			   Acsr.GetInd(), Acsr.GetPtr(),
			   SLU_NC, SLU_Z, SLU_GE);

    // We get renumbering vectors perm_r and perm_c.
    options.ColPerm = permc_spec;
    if (permc_spec != MY_PERMC)
      {
        perm_r.Reallocate(n);
        perm_c.Reallocate(n);
        perm_r.Fill();
        perm_c.Fill();

        get_perm_c(permc_spec, &AA, perm_c.GetData());
      }

    // Permuting matrix.
    Vector<int> etree(n);
    sp_preorder(&options, &AA, perm_c.GetData(), etree.GetData(), &A);

    int panel_size = sp_ienv(1);
    int relax = sp_ienv(2);
    int lwork = 0;

    // Factorization.
    zgstrf(&options, &A, relax, panel_size, etree.GetData(),
           NULL, lwork, perm_c.GetData(), perm_r.GetData(), &L, &U, &stat,
           &info_facto);

    // Clearing matrices.
    Destroy_CompCol_Permuted(&A);
    Destroy_CompCol_Matrix(&AA);

    if (info_facto == 0 && display_info)
      {
	mem_usage_t mem_usage;
	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	cout << "Number of nonzeros in factor L = " << Lstore->nnz<<endl;
	cout << "Number of nonzeros in factor U = " << Ustore->nnz<<endl;
	cout << "Number of nonzeros in L+U     = "
             << Lstore->nnz + Ustore->nnz<<endl;
	zQuerySpace(&L, &U, &mem_usage);
	cout << "Memory used for factorization in MB: "
	     << mem_usage.total_needed / (1024. * 1024.) << endl;
      }

    Acsr.Nullify();
  }


  //! resolution of linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(Vector<complex<double>, VectFull, Allocator2>& x)
  {
    trans_t trans = NOTRANS;
    int nb_rhs = 1, info;
    zCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 reinterpret_cast<doublecomplex*>(x.GetData()),
			 x.GetM(), SLU_DN, SLU_Z, SLU_GE);

    zgstrs(trans, &L, &U, perm_c.GetData(),
           perm_r.GetData(), &B, &stat, &info);
  }


  //! resolution of linear system A x = b
  template<class TransStatus, class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(const TransStatus& TransA,
        Vector<complex<double>, VectFull, Allocator2>& x)
  {
    if (TransA.NoTrans())
      {
        Solve(x);
        return;
      }

    trans_t trans = TRANS;
    int nb_rhs = 1, info;
    zCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 reinterpret_cast<doublecomplex*>(x.GetData()),
			 x.GetM(), SLU_DN, SLU_Z, SLU_GE);

    zgstrs(trans, &L, &U, perm_c.GetData(),
	   perm_r.GetData(), &B, &stat, &info);
  }


  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T, Prop, Storage, Allocator>& A, MatrixSuperLU<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  template<class T, class Allocator>
  void SolveLU(const SeldonTranspose& TransA,
	       MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }

}

#define SELDON_FILE_SUPERLU_CXX
#endif
