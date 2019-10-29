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


#ifndef SELDON_FILE_MUMPS_CXX

#include "Mumps.hxx"

namespace Seldon
{

  //! Mumps is called in double precision
  template<>
  inline void MatrixMumps<double>::CallMumps()
  {
    dmumps_c(&struct_mumps);
  }


  //! Mumps is called in complex double precision
  template<>
  inline void MatrixMumps<complex<double> >::CallMumps()
  {
    zmumps_c(&struct_mumps);
  }


  //! initialization
  template<class T>
  inline MatrixMumps<T>::MatrixMumps()
  {
    struct_mumps.comm_fortran = -987654;

    // parameters for mumps
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 0; // 0 -> unsymmetric matrix

    // other parameters
    struct_mumps.n = 0;
    type_ordering = 7; // default : we let Mumps choose the ordering
    print_level = -1;
    out_of_core = false;
  }


  //! initialization of the computation
  template<class T> template<class MatrixSparse>
  inline void MatrixMumps<T>
  ::InitMatrix(const MatrixSparse& A, bool distributed)
  {
    // we clear previous factorization
    Clear();

#ifdef SELDON_WITH_MPI
    // MPI initialization for parallel version
    if (distributed)
      {
	// for distributed matrix, every processor is assumed to be involved
	struct_mumps.comm_fortran = MPI_Comm_c2f(MPI::COMM_WORLD);
      }
    else
      {
	// centralized matrix => a linear system per processor
        struct_mumps.comm_fortran = MPI_Comm_c2f(MPI::COMM_SELF);
      }
#endif

    // symmetry is specified during the initialization stage
    struct_mumps.job = -1;
    if (IsSymmetricMatrix(A))
      struct_mumps.sym = 2; // general symmetric matrix
    else
      struct_mumps.sym = 0; // unsymmetric matrix

    // mumps is called
    CallMumps();

    struct_mumps.icntl[13] = 20;
    struct_mumps.icntl[6] = type_ordering;
    if (type_ordering == 1)
      struct_mumps.perm_in = perm.GetData();

    // setting out of core parameters
    if (out_of_core)
      struct_mumps.icntl[21] = 1;
    else
      struct_mumps.icntl[21] = 0;

    struct_mumps.icntl[17] = 0;

    // the print level is set in mumps
    if (print_level >= 0)
      {
	struct_mumps.icntl[0] = 6;
	struct_mumps.icntl[1] = 0;
	struct_mumps.icntl[2] = 6;
	struct_mumps.icntl[3] = 2;
      }
    else
      {
	struct_mumps.icntl[0] = -1;
	struct_mumps.icntl[1] = -1;
	struct_mumps.icntl[2] = -1;
	struct_mumps.icntl[3] = 0;
      }
  }


  //! selects another ordering scheme
  template<class T>
  inline void MatrixMumps<T>::SelectOrdering(int num_ordering)
  {
    type_ordering = num_ordering;
  }


  template<class T>
  inline void MatrixMumps<T>::SetPermutation(const IVect& permut)
  {
    type_ordering = 1;
    perm.Reallocate(permut.GetM());
    for (int i = 0; i < perm.GetM(); i++)
      perm(i) = permut(i) + 1;
  }


  //! clears factorization
  template<class T>
  MatrixMumps<T>::~MatrixMumps()
  {
    Clear();
  }


  //! clears factorization
  template<class T>
  inline void MatrixMumps<T>::Clear()
  {
    if (struct_mumps.n > 0)
      {
	struct_mumps.job = -2;
	// Mumps variables are deleted.
        CallMumps();

	num_row_glob.Clear();
	num_col_glob.Clear();
	struct_mumps.n = 0;
      }
  }


  //! no display from Mumps
  template<class T>
  inline void MatrixMumps<T>::HideMessages()
  {
    print_level = -1;

    struct_mumps.icntl[0] = -1;
    struct_mumps.icntl[1] = -1;
    struct_mumps.icntl[2] = -1;
    struct_mumps.icntl[3] = 0;

  }


  //! standard display
  template<class T>
  inline void MatrixMumps<T>::ShowMessages()
  {
    print_level = 0;

    struct_mumps.icntl[0] = 6;
    struct_mumps.icntl[1] = 0;
    struct_mumps.icntl[2] = 6;
    struct_mumps.icntl[3] = 2;

  }


  //! Enables writing on the disk (out of core).
  template<class T>
  inline void MatrixMumps<T>::EnableOutOfCore()
  {
    out_of_core = true;
  }


  //! Disables writing on the disk (incore).
  template<class T>
  inline void MatrixMumps<T>::DisableOutOfCore()
  {
    out_of_core = false;
  }


  //! Computes an ordering for matrix renumbering.
  /*!
    \param[in,out] mat matrix whose we want to find the ordering
    \param[out] numbers new row numbers
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop,class Storage,class Allocator>
  void MatrixMumps<T>::FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
				    IVect& numbers, bool keep_matrix)
  {
    InitMatrix(mat);

    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format
    IVect num_row, num_col; Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);

    // no values needed to renumber
    values.Clear();
    if (!keep_matrix)
      mat.Clear();

    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.irn = num_row.GetData();
    struct_mumps.jcn = num_col.GetData();

    // Call the MUMPS package.
    struct_mumps.job = 1; // we analyse the system
    CallMumps();

    numbers.Reallocate(n);
    for (int i = 0; i < n; i++)
      numbers(i) = struct_mumps.sym_perm[i]-1;
  }


  //! factorization of a given matrix
  /*!
    \param[in,out] mat matrix to factorize
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>::FactorizeMatrix(Matrix<T,Prop,Storage,Allocator> & mat,
				       bool keep_matrix)
  {
    InitMatrix(mat);

    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format with fortran convention (1-index)
    IVect num_row, num_col; Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    if (!keep_matrix)
      mat.Clear();

    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.irn = num_row.GetData();
    struct_mumps.jcn = num_col.GetData();
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
  }


  //! Symbolic factorization
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>
  ::PerformAnalysis(Matrix<T, Prop, Storage, Allocator> & mat)
  {
    InitMatrix(mat);

    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format with fortran convention (1-index)
    Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row_glob, num_col_glob, values, 1);

    struct_mumps.n = n; struct_mumps.nz = nnz;
    struct_mumps.irn = num_row_glob.GetData();
    struct_mumps.jcn = num_col_glob.GetData();
    struct_mumps.a = reinterpret_cast<pointer>(values.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 1; // we analyse the system
    CallMumps();
  }


  //! Numerical factorization
  /*!
    Be careful, because no conversion is performed in the method,
    so you have to choose RowSparse/ColSparse for unsymmetric matrices
    and RowSymSparse/ColSymSparse for symmetric matrices.
    The other formats should not work
  */
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>
  ::PerformFactorization(Matrix<T, Prop, Storage, Allocator> & mat)
  {
    // we consider that the values are corresponding
    // to the row/column numbers given for the analysis
    struct_mumps.a = reinterpret_cast<pointer>(mat.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 2; // we factorize the system
    CallMumps();
  }


  //! returns information about factorization performed
  template<class T>
  int MatrixMumps<T>::GetInfoFactorization() const
  {
    return struct_mumps.info[0];
  }


  //! Computation of Schur complement.
  /*!
    \param[in,out] mat initial matrix.
    \param[in] num numbers to keep in Schur complement.
    \param[out] mat_schur Schur matrix.
    \param[in] keep_matrix if false, \a mat is cleared.
  */
  template<class T> template<class Prop1, class Storage1, class Allocator,
			     class Prop2, class Storage2, class Allocator2>
  void MatrixMumps<T>::
  GetSchurMatrix(Matrix<T, Prop1, Storage1, Allocator>& mat, const IVect& num,
		 Matrix<T, Prop2, Storage2, Allocator2> & mat_schur,
		 bool keep_matrix)
  {
    InitMatrix(mat);

    int n_schur = num.GetM(), n = mat.GetM();
    // Subscripts are changed to respect fortran convention
    IVect index_schur(n_schur);
    for (int i = 0; i < n_schur; i++)
      index_schur(i) = num(i)+1;

    // array that will contain values of Schur matrix
    Vector<T, VectFull, Allocator2> vec_schur(n_schur*n_schur);

    struct_mumps.icntl[18] = n_schur;
    struct_mumps.size_schur = n_schur;
    struct_mumps.listvar_schur = index_schur.GetData();
    struct_mumps.schur = reinterpret_cast<pointer>(vec_schur.GetData());

    // factorization of the matrix
    FactorizeMatrix(mat, keep_matrix);

    // resetting parameters related to Schur complement
    struct_mumps.icntl[18] = 0;
    struct_mumps.size_schur = 0;
    struct_mumps.listvar_schur = NULL;
    struct_mumps.schur = NULL;

    // schur complement stored by rows
    int nb = 0;
    mat_schur.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n ;j++)
	mat_schur(i,j) = vec_schur(nb++);

    vec_schur.Clear(); index_schur.Clear();
  }


  //! resolution of a linear system using the computed factorization
  /*!
    \param[in,out] x right-hand-side on input, solution on output
    It is assumed that a call to FactorizeMatrix has been done before
  */
  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixMumps<T>::Solve(const Transpose_status& TransA,
			     Vector<T, VectFull, Allocator2>& x)
  {
#ifdef SELDON_CHECK_DIMENSIONS
    if (x.GetM() != struct_mumps.n)
      throw WrongDim("Mumps::Solve(TransA, c)",
                     string("The length of x is equal to ")
                     + to_str(x.GetM())
                     + " while the size of the matrix is equal to "
                     + to_str(struct_mumps.n) + ".");
#endif

    if (TransA.Trans())
      struct_mumps.icntl[8] = 0;
    else
      struct_mumps.icntl[8] = 1;

    struct_mumps.rhs = reinterpret_cast<pointer>(x.GetData());
    struct_mumps.job = 3; // we solve system
    CallMumps();
  }



  template<class T> template<class Allocator2>
  void MatrixMumps<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }


  //! Resolution of a linear system using a computed factorization.
  /*!
    \param[in,out] x on entry, the right-hand-side; on exit, the solution.
    It is assumed that 'FactorizeMatrix' has already been called.
  */
  template<class T>
  template<class Allocator2, class Transpose_status, class Prop>
  void MatrixMumps<T>::Solve(const Transpose_status& TransA,
			     Matrix<T, Prop, ColMajor, Allocator2>& x)
  {

#ifdef SELDON_CHECK_BOUNDS
    if (x.GetM() != struct_mumps.n)
      throw WrongDim("Mumps::Solve", string("Row size of x is equal to ")
                     + to_str(x.GetM()) + " while size of matrix is equal to "
                     + to_str(struct_mumps.n));
#endif

    if (TransA.Trans())
      struct_mumps.icntl[8] = 0;
    else
      struct_mumps.icntl[8] = 1;

    struct_mumps.nrhs = x.GetN();
    struct_mumps.lrhs = x.GetM();
    struct_mumps.rhs = reinterpret_cast<pointer>(x.GetData());
    struct_mumps.job = 3; // we solve system
    CallMumps();
  }


#ifdef SELDON_WITH_MPI
  //! factorization of a given matrix in distributed form (parallel execution)
  /*!
    \param[in,out] comm_facto MPI communicator
    \param[in,out] Ptr start indices for each column
    \param[in,out] IndRow row indices
    \param[in,out] Val data
    \param[in] sym if true, the matrix is assumed to be symmetric (upper part
    is provided)
    \param[in] glob_number row numbers (in the global matrix)
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T>
  template<class Alloc1, class Alloc2, class Alloc3, class Tint>
  void MatrixMumps<T>::
  FactorizeDistributedMatrix(MPI::Comm& comm_facto,
                             Vector<Tint, VectFull, Alloc1>& Ptr,
                             Vector<Tint, VectFull, Alloc2>& IndRow,
                             Vector<T, VectFull, Alloc3>& Val,
                             const Vector<Tint>& glob_number,
                             bool sym, bool keep_matrix)
  {
    // Initialization depending on symmetry of the matrix.
    if (sym)
      {
        Matrix<T, Symmetric, RowSymSparse, Alloc3> Atest;
        InitMatrix(Atest, true);
      }
    else
      {
        Matrix<T, General, RowSparse, Alloc3> Atest;
        InitMatrix(Atest, true);
      }

    // Fortran communicator
    struct_mumps.comm_fortran = MPI_Comm_c2f(comm_facto);

    // distributed matrix
    struct_mumps.icntl[17] = 3;

    // finding the size of the overall system
    Tint nmax = 0, N = 0;
    for (int i = 0; i < glob_number.GetM(); i++)
      nmax = max(glob_number(i)+1, nmax);

    comm_facto.Allreduce(&nmax, &N, 1, MPI::INTEGER, MPI::MAX);

    // number of non-zero entries on this processor
    int nnz = IndRow.GetM();

    // conversion in coordinate format
    Vector<Tint, VectFull, Alloc2> IndCol(nnz);
    for (int i = 0; i < IndRow.GetM(); i++)
      IndRow(i)++;

    for (int i = 0; i < Ptr.GetM()-1; i++)
      for (int j = Ptr(i); j < Ptr(i+1); j++)
        IndCol(j) = glob_number(i) + 1;

    if (!keep_matrix)
      Ptr.Clear();

    // Define the problem on the host
    struct_mumps.n = N;
    struct_mumps.nz_loc = nnz;
    struct_mumps.irn_loc = IndRow.GetData();
    struct_mumps.jcn_loc = IndCol.GetData();
    struct_mumps.a_loc = reinterpret_cast<pointer>(Val.GetData());

    // Call the MUMPS package.
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();

    if ((comm_facto.Get_rank() == 0) && (print_level >= 0))
      cout<<"Factorization completed"<<endl;

  }


  //! solves linear system with parallel execution
  /*!
    \param[in] TransA we solve A x = b or A^T x = b
    \param[in,out] x right-hand-side then solution
    \param[in,out] glob_num global row numbers
  */
  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixMumps<T>::SolveDistributed(MPI::Comm& comm_facto,
                                        const Transpose_status& TransA,
					Vector<T, VectFull, Allocator2>& x,
					const IVect& glob_num)
  {
    Vector<T, VectFull, Allocator2> rhs;
    int cplx = sizeof(T)/8;
    // allocating the global right hand side
    rhs.Reallocate(struct_mumps.n); rhs.Zero();

    if (comm_facto.Get_rank() == 0)
      {
	// on the host, we retrieve datas of all the other processors
	int nb_procs = comm_facto.Get_size();
        MPI::Status status;
	if (nb_procs > 1)
	  {
	    // assembling the right hand side
	    Vector<T, VectFull, Allocator2> xp;
	    IVect nump;
	    for (int i = 0; i < nb_procs; i++)
	      {

		if (i != 0)
		  {
                    // On the host processor receiving components of right
                    // hand side.
		    int nb_dof;
		    comm_facto.Recv(&nb_dof, 1, MPI::INTEGER, i, 34, status);

                    xp.Reallocate(nb_dof);
		    nump.Reallocate(nb_dof);

                    comm_facto.Recv(xp.GetDataVoid(), cplx*nb_dof,
                                    MPI::DOUBLE, i, 35, status);

		    comm_facto.Recv(nump.GetData(), nb_dof, MPI::INTEGER,
                                    i, 36, status);
		  }
		else
		  {
		    xp = x; nump = glob_num;
		  }

		for (int j = 0; j < nump.GetM(); j++)
		  rhs(nump(j)) = xp(j);
	      }
	  }
	else
	  Copy(x, rhs);

	struct_mumps.rhs = reinterpret_cast<pointer>(rhs.GetData());
      }
    else
      {
	// On other processors, we send right hand side.
	int nb = x.GetM();
	comm_facto.Send(&nb, 1, MPI::INTEGER, 0, 34);
	comm_facto.Send(x.GetDataVoid(), cplx*nb, MPI::DOUBLE, 0, 35);
	comm_facto.Send(glob_num.GetData(), nb, MPI::INTEGER, 0, 36);
      }

    // we solve system
    if (TransA.Trans())
      struct_mumps.icntl[8] = 0;
    else
      struct_mumps.icntl[8] = 1;

    struct_mumps.job = 3;
    CallMumps();

    // we distribute solution on all the processors
    comm_facto.Bcast(rhs.GetDataVoid(), cplx*rhs.GetM(), MPI::DOUBLE, 0);

    // and we extract the solution on provided numbers
    for (int i = 0; i < x.GetM(); i++)
      x(i) = rhs(glob_num(i));
  }


  template<class T> template<class Allocator2, class Tint>
  void MatrixMumps<T>::SolveDistributed(MPI::Comm& comm_facto,
                                        Vector<T, VectFull, Allocator2>& x,
					const Vector<Tint>& glob_num)
  {
    SolveDistributed(comm_facto, SeldonNoTrans, x, glob_num);
  }
#endif


  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T,Symmetric,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T,General,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }


  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T, Symmetric, Storage, Allocator>& A,
                      MatrixMumps<T>& mat_lu, const IVect& num,
                      MatrixFull& schur_matrix, bool keep_matrix = false)
  {
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }


  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T, General, Storage, Allocator>& A,
                      MatrixMumps<T>& mat_lu, const IVect& num,
                      MatrixFull& schur_matrix, bool keep_matrix = false)
  {
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }


  template<class T, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }


  template<class T, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }


  template<class T, class Prop, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu,
               Matrix<T, Prop, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(SeldonNoTrans, x);
  }


  template<class T, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       MatrixMumps<T>& mat_lu, Matrix<T, ColMajor, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }

}

#define SELDON_FILE_MUMPS_CXX
#endif
