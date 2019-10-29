// Copyright (C) 2001-2010 Marc Duruflé
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

#ifndef SELDON_FILE_PASTIX_HXX

// including Pastix headers
extern "C"
{
#define _COMPLEX_H

#include "pastix.h"
}

namespace Seldon
{

  template<class T>
  class MatrixPastix
  {
  protected :
    //! pastix structure
    pastix_data_t* pastix_data;
    //! options (integers)
    pastix_int_t iparm[64];
    //! options (floats)
    double dparm[64];
    //! number of columns
    pastix_int_t n;
    //! permutation arrays
    Vector<pastix_int_t> perm, invp;
    //! local to global
    Vector<pastix_int_t> col_num;
    //! if true, resolution on several nodes
    bool distributed;
    //! level of display
    int print_level;
    //! if true, solution is refined
    bool refine_solution;

  public :

    MatrixPastix();
    ~MatrixPastix();

    void Clear();

    void CallPastix(const MPI_Comm&, pastix_int_t* colptr, pastix_int_t* row,
                    T* val, T* b, pastix_int_t nrhs);

    void HideMessages();
    void ShowMessages();
    void ShowFullHistory();

    void SelectOrdering(int type);
    void SetPermutation(const IVect& permut);

    void RefineSolution();
    void DoNotRefineSolution();

    template<class T0, class Prop, class Storage, class Allocator, class Tint>
    void FindOrdering(Matrix<T0, Prop, Storage, Allocator> & mat,
		      Vector<Tint>& numbers, bool keep_matrix = false);

    template<class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, General, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Storage, class Allocator>
    void FactorizeMatrix(Matrix<T, Symmetric, Storage, Allocator> & mat,
			 bool keep_matrix = false);

    template<class Allocator2>
    void Solve(Vector<T, VectFull, Allocator2>& x);

    template<class Allocator2, class Transpose_status>
    void Solve(const Transpose_status& TransA,
	       Vector<T, VectFull, Allocator2>& x);

    void SetNumberThreadPerNode(int num_thread);

    template<class Alloc1, class Alloc2, class Alloc3, class Tint>
    void FactorizeDistributedMatrix(MPI::Comm& comm_facto,
                                    Vector<pastix_int_t, VectFull, Alloc1>&,
                                    Vector<pastix_int_t, VectFull, Alloc2>&,
                                    Vector<T, VectFull, Alloc3>&,
                                    const Vector<Tint>& glob_number,
				    bool sym, bool keep_matrix = false);

    template<class Allocator2, class Tint>
    void SolveDistributed(MPI::Comm& comm_facto,
                          Vector<T, Vect_Full, Allocator2>& x,
                          const Vector<Tint>& glob_num);

    template<class Allocator2, class Transpose_status, class Tint>
    void SolveDistributed(MPI::Comm& comm_facto,
                          const Transpose_status& TransA,
			  Vector<T, Vect_Full, Allocator2>& x,
                          const Vector<Tint>& glob_num);

  };

}

#define SELDON_FILE_PASTIX_HXX
#endif

