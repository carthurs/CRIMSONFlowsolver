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


#ifndef SELDON_FILE_SELDON_HXX


#include "SeldonHeader.hxx"

namespace Seldon
{


  // Vector class - specialized for each used type.
  template <class T, class Storage, class Allocator>
  class Vector
  {
    // Nothing in it: no default vector is supplied so as to avoid suprises!
  };


  // Matrix class - specialized for each used type.
  template <class T, class Prop, class Storage, class Allocator>
  class Matrix
  {
    // Nothing in it: no default matrix is supplied so as to avoid suprises!
  };


} // namespace Seldon.


#include "share/Common.cxx"

// Memory management.
#include "share/Allocator.cxx"

// Storage type.
#include "share/Storage.cxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/MatrixFlag.cxx"
#include "share/Errors.cxx"
#endif

#include "array3d/Array3D.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#include "matrix/Matrix_Triangular.cxx"
#include "matrix/Matrix_Symmetric.cxx"
#include "matrix/Matrix_Hermitian.cxx"
#include "matrix_sparse/Matrix_Sparse.cxx"
#include "matrix_sparse/Matrix_ComplexSparse.cxx"
#include "matrix_sparse/Matrix_SymSparse.cxx"
#include "matrix_sparse/Matrix_SymComplexSparse.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_HermPacked.cxx"
#include "matrix/Matrix_TriangPacked.cxx"
#include "vector/Vector.cxx"
#include "vector/VectorCollection.cxx"
#include "vector/Functions_Arrays.cxx"
#include "vector/SparseVector.cxx"
#include "matrix/Functions.cxx"
#include "matrix_sparse/Matrix_Conversions.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"

#include "matrix/SubMatrix_Base.cxx"
#include "matrix/SubMatrix.cxx"

// Blas interface.
#ifdef SELDON_WITH_BLAS
#include "computation/interfaces/Blas_1.cxx"
#include "computation/interfaces/Blas_2.cxx"
#include "computation/interfaces/Blas_3.cxx"
#endif

// Lapack interface.
#ifdef SELDON_WITH_LAPACK
#include "computation/interfaces/Lapack_LinearEquations.cxx"
#include "computation/interfaces/Lapack_LeastSquares.cxx"
#include "computation/interfaces/Lapack_Eigenvalues.cxx"
#endif // SELDON_WITH_LAPACK.

#ifdef SELDON_WITH_COMPILED_LIBRARY
#include "lib/Common.cpp"
#include "lib/Vector.cpp"
#include "lib/MatrixPointers.cpp"
#endif


#define SELDON_FILE_SELDON_HXX
#endif
