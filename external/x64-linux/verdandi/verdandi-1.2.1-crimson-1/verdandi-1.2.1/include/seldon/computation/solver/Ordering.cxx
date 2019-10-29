// Copyright (C) 2011 Marc Duruflé
// Copyright (C) 2011 Vivien Mallet
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


#ifndef SELDON_FILE_SOLVER_ORDERING_CXX
#define SELDON_FILE_SOLVER_ORDERING_CXX


#include "Ordering.hxx"


namespace Seldon
{


  //! Constructs reverse Cuthill-McKee ordering from a given matrix.
  template<class T, class Prop, class Storage, class Allocator,
	   class Tint, class Alloc>
  void FindReverseCuthillMcKeeOrdering(const Matrix<T, Prop,
				       Storage, Allocator>& A,
				       Vector<Tint, VectFull, Alloc>& num)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	num.Clear();
	return;
      }

    // Pattern of A + A' is retrieved in CSC format.
    Vector<Tint, VectFull, Alloc> Ptr, Ind;
    Vector<T, VectFull, Allocator> Value;
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Value, true);
    Value.Clear();

    // Allocation of arrays.
    Vector<Tint, VectFull, Alloc> vertex(n), neighbor(n);
    Vector<bool> RowUsed(n);
    IVect nb_neighbor(n);

    num.Reallocate(n);
    num.Fill(0);
    RowUsed.Fill(false);
    vertex.Fill(0);
    neighbor.Fill(0);
    nb_neighbor.Fill(0);

    // Starting with the first row.
    int nb_row = 1;
    int nb_active_row = 1;
    vertex(0) = 0;
    num(0) = 0;
    RowUsed(0) = true;
    int first_free_row = 1;

    // Loop until all rows are found.
    while (nb_row < n)
      {
	// Searching all neighboring rows.
	int nb = 0;
	for (int i = 0; i < nb_active_row; i++)
	  {
	    int irow = vertex(i);
	    for (int j = Ptr(irow); j < Ptr(irow+1); j++)
	      {
		int icol = Ind(j);
		if (!RowUsed(icol))
		  {
		    RowUsed(icol) = true;
		    neighbor(nb) = icol;
		    nb_neighbor(nb) = Ptr(icol+1) - Ptr(icol);
		    nb++;
		  }
	      }
	  }

	if (nb == 0)
	  {
	    // No row has been found, therefore there are independent blocks
	    // in the matrix; searching next free row.
	    while (RowUsed(first_free_row))
	      first_free_row++;

	    RowUsed(first_free_row) = true;
	    num(nb_row) = first_free_row;
	    nb_active_row = 1;
	    vertex(0) = first_free_row;
	    nb_row++;
	  }
	else
	  {
	    // Sorting neighboring vertices by the number of adjacent
	    // vertices.
	    Sort(nb, nb_neighbor, neighbor);

	    // Then adding those rows.
	    nb_active_row = nb;
	    for (int i = 0; i < nb; i++)
	      {
		vertex(i) = neighbor(i);
		num(nb_row) = neighbor(i);
		nb_row++;
	      }
	  }
      }

    // Reverting final result.
    for (int i = 0; i < n; i++)
      vertex(i) = num(i);

    for (int i = 0; i < n; i++)
      num(n-1-i) = vertex(i);
  }


  //! Constructs an ordering for the factorization of a sparse matrix.
  template<class T, class Prop, class Storage, class Allocator,
	   class Tint, class Alloc>
  void FindSparseOrdering(Matrix<T, Prop, Storage, Allocator>& A,
			  Vector<Tint, VectFull, Alloc>& num, int type)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	num.Clear();
	return;
      }

    switch (type)
      {
      case SparseMatrixOrdering::IDENTITY :
	{
	  // No permutation.
	  num.Reallocate(n);
	  num.Fill();
	}
	break;

      case SparseMatrixOrdering::REVERSE_CUTHILL_MCKEE :
	{
	  // Reverse Cuthill-McKee.
	  FindReverseCuthillMcKeeOrdering(A, num);
	}
	break;

      case SparseMatrixOrdering::PORD :
	{
	  // Pord package provided by Mumps.
#ifdef SELDON_WITH_MUMPS
	  MatrixMumps<T> mat_lu;
	  mat_lu.SelectOrdering(4);
	  mat_lu.FindOrdering(A, num, true);
#else
          throw Error("FindSparseOrdering(Matrix&, Vector&, int)",
                      "PORD is supported when Mumps is available.");
#endif
	}
	break;

      case SparseMatrixOrdering::SCOTCH :
	{
	  // Scotch interface provided by Pastix.
#ifdef SELDON_WITH_PASTIX
	  MatrixPastix<T> mat_lu;
	  mat_lu.FindOrdering(A, num, true);
#else
          throw Error("FindSparseOrdering(Matrix&, Vector&, int)",
                      "SCOTCH is supported when Mumps is available.");
#endif
	}
	break;

      case SparseMatrixOrdering::METIS :
	{
	  // Metis ordering provided by Mumps.
#ifdef SELDON_WITH_MUMPS
	  MatrixMumps<T> mat_lu;
	  mat_lu.SelectOrdering(5);
	  mat_lu.FindOrdering(A, num, true);
#else
          throw Error("FindSparseOrdering(Matrix&, Vector&, int)",
                      "METIS is supported when Mumps is available.");
#endif
	}
	break;

      case SparseMatrixOrdering::AMD :
	{
	  // Camd package in SuiteSparse (containing UmfPack).
#ifdef SELDON_WITH_UMFPACK

	  num.Reallocate(n);

	  // pattern of A+A' is retrieved in CSC format
	  Vector<Tint, VectFull, Alloc> Ptr, Ind;
	  Vector<T, VectFull, Allocator> Value;
	  General sym;
	  ConvertToCSC(A, sym, Ptr, Ind, Value, true);
	  Value.Clear();

	  // then calling camd
	  Vector<Tint, VectFull, Alloc> C(n);
	  C.Fill(0);
	  double Control[CAMD_CONTROL], Info[CAMD_INFO];
	  camd_defaults(Control);
	  camd_order(n, Ptr.GetData(), Ind.GetData(),
		     num.GetData(), Control, Info, C.GetData());

#else
          throw Error("FindSparseOrdering(Matrix&, Vector&, int)",
                      "AMD is supported when UmfPack is available.");
#endif
	}
	break;

      case SparseMatrixOrdering::COLAMD :
	{
	  // Colamd package in SuiteSparse (containing UmfPack).
#ifdef SELDON_WITH_UMFPACK

	  num.Reallocate(n);

	  // Pattern of A + A' is retrieved in CSC format.
	  Vector<Tint, VectFull, Alloc> Ptr, Ind;
	  Vector<T, VectFull, Allocator> Value;
	  General sym;
	  ConvertToCSC(A, sym, Ptr, Ind, Value, true);
	  Value.Clear();

	  // then calling colamd
	  int Stats[COLAMD_STATS];
	  colamd(n, n, Ind.GetM(), Ind.GetData(), Ptr.GetData(), NULL, Stats);
#else
          throw Error("FindSparseOrdering(Matrix&, Vector&, int)",
                      "COLAMD is supported when UmfPack is available.");
#endif
	}
	break;

      case SparseMatrixOrdering::QAMD :
	{
	  // Qamd ordering provided by Mumps
#ifdef SELDON_WITH_MUMPS
	  MatrixMumps<T> mat_lu;
	  mat_lu.SelectOrdering(6);
	  mat_lu.FindOrdering(A, num, true);
#else
          throw Error("FindSparseOrdering(Matrix&, Vector&, int)",
                      "QAMD is supported when Mumps is available.");
#endif
	}
	break;

      case SparseMatrixOrdering::USER :
	// nothing to do
	break;
      }
  }

} // namespace Seldon.


#endif
