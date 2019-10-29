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


#ifndef SELDON_FILE_MUMPS_HXX

// including Mumps headers
extern "C"
{
#include "dmumps_c.h"
#include "zmumps_c.h"

  // including mpi from sequential version of Mumps if the
  // compilation is not made on a parallel machine
#ifndef SELDON_WITH_MPI
#include "mpi.h"
#endif

}

namespace Seldon
{
  template<class T>
  class TypeMumps
  {
  };


  //! class containing MUMPS data structure
  template<>
  class TypeMumps<double>
  {
  public :
    typedef DMUMPS_STRUC_C data;
    typedef double* pointer;
  };


  //! class containing MUMPS data structure
  template<>
  class TypeMumps<complex<double> >
  {
  public :
    typedef ZMUMPS_STRUC_C data;
    typedef mumps_double_complex* pointer;
  };


  //! object used to solve linear system by calling mumps subroutines
  template<class T>
  class MatrixMumps
  {
  protected :
    int type_ordering; //!< ordering scheme (AMD, Metis, etc)
    //! object containing Mumps data structure
    typename TypeMumps<T>::data struct_mumps;
    //! double* or complex<double>*
    typedef typename TypeMumps<T>::pointer pointer;
    int print_level;
    bool out_of_core;
    IVect num_row_glob, num_col_glob;
    IVect perm;

    // internal methods
    void CallMumps();
    template<class MatrixSparse>
    void InitMatrix(const MatrixSparse&, bool dist = false);

  public :
    MatrixMumps();
    ~MatrixMumps();

    void Clear();

    void SelectOrdering(int num_ordering);
    void SetPermutation(const IVect& permut);

    void HideMessages();
    void ShowMessages();

    void EnableOutOfCore();
    void DisableOutOfCore();

    int GetInfoFactorization() const;

    template<class Prop, class Storage, class Allocator>
    void FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
		      IVect& numbers, bool keep_matrix = false);

    template<class Prop, class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, Prop, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Prop, class Storage, class Allocator>
    void PerformAnalysis(Matrix<T, Prop, Storage, Allocator> & mat);

    template<class Prop, class Storage, class Allocator>
    void PerformFactorization(Matrix<T, Prop, Storage, Allocator> & mat);

    template<class Prop1, class Storage1, class Allocator1,
	     class Prop2, class Storage2, class Allocator2>
    void GetSchurMatrix(Matrix<T, Prop1, Storage1, Allocator1>& mat,
			const IVect& num,
			Matrix<T, Prop2, Storage2, Allocator2> & mat_schur,
			bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2, class Transpose_status>
    void Solve(const Transpose_status& TransA,
	       Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2, class Transpose_status, class Prop>
    void Solve(const Transpose_status& TransA,
	       Matrix<T, Prop, ColMajor, Allocator2>& x);

#ifdef SELDON_WITH_MPI
    template<class Alloc1, class Alloc2, class Alloc3, class Tint>
    void FactorizeDistributedMatrix(MPI::Comm& comm_facto,
                                    Vector<Tint, VectFull, Alloc1>&,
                                    Vector<Tint, VectFull, Alloc2>&,
                                    Vector<T, VectFull, Alloc3>&,
                                    const Vector<Tint>& glob_number,
				    bool sym, bool keep_matrix = false);

    template<class Allocator2, class Tint>
    void SolveDistributed(MPI::Comm& comm_facto,
                          Vector<T, Vect_Full, Allocator2>& x,
                          const Vector<Tint>& glob_num);

    template<class Allocator2, class Transpose_status>
    void SolveDistributed(MPI::Comm& comm_facto,
                          const Transpose_status& TransA,
			  Vector<T, VectFull, Allocator2>& x,
			  const IVect& glob_num);

#endif

  };

}

#define SELDON_FILE_MUMPS_HXX
#endif



