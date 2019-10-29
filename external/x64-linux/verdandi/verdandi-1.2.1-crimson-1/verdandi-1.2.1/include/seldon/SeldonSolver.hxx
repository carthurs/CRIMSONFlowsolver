// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_SELDON_SOLVER_HXX

// additional classes and functions for sparse matrices
#include "matrix_sparse/Matrix_Conversions.cxx"
#include "matrix_sparse/Matrix_ArraySparse.cxx"
#include "matrix_sparse/Matrix_ArrayComplexSparse.cxx"
#include "matrix_sparse/Permutation_ScalingMatrix.cxx"
#include "matrix_sparse/Relaxation_MatVect.cxx"
#include "matrix_sparse/Functions_MatrixArray.cxx"


// interfaces with direct solvers
#ifdef SELDON_WITH_MUMPS
#include "computation/interfaces/direct/Mumps.cxx"
#endif

#ifdef SELDON_WITH_UMFPACK
#include "computation/interfaces/direct/UmfPack.cxx"
#endif

#ifdef SELDON_WITH_SUPERLU
#include "computation/interfaces/direct/SuperLU.cxx"
#endif

#ifdef SELDON_WITH_PASTIX
#include "computation/interfaces/direct/Pastix.cxx"
#endif

#ifdef SELDON_WITH_PRECONDITIONING
#include "SeldonPreconditioner.hxx"
#endif

#include "computation/solver/SparseSolver.cxx"

// iterative solvers and preconditioning
#include "computation/solver/iterative/Iterative.cxx"
#include "computation/solver/preconditioner/Precond_Ssor.cxx"

// Cholesky Solver
#ifdef SELDON_WITH_CHOLMOD
#include "computation/interfaces/direct/Cholmod.cxx"
#endif

#include "computation/solver/SparseCholeskyFactorisation.cxx"

#define SELDON_FILE_SELDON_SOLVER_HXX
#endif
