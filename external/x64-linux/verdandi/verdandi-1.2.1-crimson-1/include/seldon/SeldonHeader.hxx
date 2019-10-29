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


#ifndef SELDON_FILE_SELDONHEADER_HXX

#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <stdexcept>

// For the compiled library.
#ifdef SELDON_WITH_COMPILED_LIBRARY
#define SELDON_EXTERN extern
#else
#define SELDON_EXTERN
#endif

// For backward compatibility.
#ifdef SELDON_WITH_CBLAS
#define SELDON_WITH_BLAS
#endif

#ifdef SELDON_WITH_BLAS
extern "C"
{
#include "computation/interfaces/cblas.h"
}
#endif


//////////////////
// DEBUG LEVELS //
//////////////////

#ifdef SELDON_DEBUG_LEVEL_4
#ifndef SELDON_DEBUG_LEVEL_3
#define SELDON_DEBUG_LEVEL_3
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_3
#ifndef SELDON_CHECK_BOUNDS
#define SELDON_CHECK_BOUNDS
#endif
#ifndef SELDON_DEBUG_LEVEL_2
#define SELDON_DEBUG_LEVEL_2
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_2
#ifndef SELDON_CHECK_DIMENSIONS
#define SELDON_CHECK_DIMENSIONS
#endif
#ifndef SELDON_DEBUG_LEVEL_1
#define SELDON_DEBUG_LEVEL_1
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_1
#ifndef SELDON_LAPACK_CHECK_INFO
#define SELDON_LAPACK_CHECK_INFO
#endif
#ifndef SELDON_CHECK_MEMORY
#define SELDON_CHECK_MEMORY
#endif
#ifndef SELDON_CHECK_IO
#define SELDON_CHECK_IO
#endif
#ifndef SELDON_DEBUG_LEVEL_0
#define SELDON_DEBUG_LEVEL_0
#endif
#endif

#ifdef SELDON_DEBUG_LEVEL_0
#ifndef SELDON_DEBUG_LEVEL_1
#ifndef SELDON_WITHOUT_THROW
#define SELDON_WITHOUT_THROW
#endif
#endif
#endif

// Convenient macros to catch exceptions.
#ifndef TRY
#define TRY try {
#endif
#ifndef END
#define END                                                             \
  }                                                                     \
    catch(Seldon::Error& Err)                                           \
      {                                                                 \
        Err.CoutWhat();                                                 \
        return 1;                                                       \
      }                                                                 \
    catch (std::exception& Err)                                         \
      {                                                                 \
        std::cout << "C++ exception: " << Err.what() << std::endl;      \
        return 1;                                                       \
      }                                                                 \
    catch (std::string& str)                                            \
      {                                                                 \
        std::cout << str << std::endl;                                  \
        return 1;                                                       \
      }                                                                 \
    catch (const char* str)                                             \
      {                                                                 \
        std::cout << str << std::endl;                                  \
        return 1;                                                       \
      }                                                                 \
    catch(...)                                                          \
      {                                                                 \
        std::cout << "Unknown exception..." << std::endl;               \
        return 1;                                                       \
      }
#endif

//! To display a message... call Hermes!
#ifndef ERR
#define ERR(x) std::cout << "Hermes - " #x << std::endl
#endif
//! To display a variable (with its name); same as DISPLAY.
#ifndef DISP
#define DISP(x) std::cout << #x ": " << x << std::endl
#endif
//! To display a variable (with its name); same as DISP.
#ifndef DISPLAY
#define DISPLAY(x) std::cout << #x ": " << x << std::endl
#endif

// For backward compatibility. These lines should be removed one day.
#define Vect_Full VectFull
#define Vect_Sparse VectSparse

//! Seldon namespace.
namespace Seldon
{
  using namespace std;
}

// Exceptions and useful functions.
#include "share/Errors.hxx"
#include "share/Common.hxx"

// Default allocator.
#ifndef SELDON_DEFAULT_ALLOCATOR
#define SELDON_DEFAULT_ALLOCATOR MallocAlloc
#endif
// Memory management.
#include "share/Allocator.hxx"

// Storage type.
#include "share/Storage.hxx"

#include "share/MatrixFlag.hxx"

// Properties.
#include "share/Properties.hxx"

namespace Seldon
{


  // Base structure for all vectors.
  template <class T, class Allocator>
  class Vector_Base;

  // Vector class - specialized for each used type.
  template <class T, class Storage = VectFull,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Vector;

  // Full vector.
  template <class T, class Allocator>
  class Vector<T, VectFull, Allocator>;

  // Sparse vector.
  template <class T, class Allocator>
  class Vector<T, VectSparse, Allocator>;

  // Vector collection.
  template <class T, class Allocator>
  class Vector<T, Collection, Allocator>;

  // Matrix class - specialized for each used type.
  template <class T, class Prop = General,
	    class Storage = RowMajor,
	    class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Matrix;

  // column-major matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColMajor, Allocator>;

  // row-major matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowMajor, Allocator>;

  // column-major symmetric packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymPacked, Allocator>;

  // row-major symmetric packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymPacked, Allocator>;

  // column-major upper-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColUpTriangPacked, Allocator>;

  // column-major lower-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColLoTriangPacked, Allocator>;

  // row-major upper-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowUpTriangPacked, Allocator>;

  // row-major lower-triangular packed matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowLoTriangPacked, Allocator>;

  // column-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSparse, Allocator>;

  // row-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSparse, Allocator>;

  // column-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymSparse, Allocator>;

  // row-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymSparse, Allocator>;

  // column-major complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColComplexSparse, Allocator>;

  // row-major complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowComplexSparse, Allocator>;

  // column-major symmetric complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColSymComplexSparse, Allocator>;

  // row-major symmetric complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, RowSymComplexSparse, Allocator>;

  // column-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSparse, Allocator>;

  // row-major sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSparse, Allocator>;

  // column-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymSparse, Allocator>;

  // row-major symmetric sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymSparse, Allocator>;

  // column-major complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ColComplexSparse, Allocator>;

  // row-major complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowComplexSparse, Allocator>;

  // column-major symmetric complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>;

  // row-major symmetric complex sparse matrix.
  template <class T, class Prop, class Allocator>
  class Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>;

  // 3D array.
  template <class T, class Allocator>
  class Array3D;


} // namespace Seldon.


#include "array3d/Array3D.hxx"
#include "matrix/Matrix_Base.hxx"
#include "matrix/Matrix_Pointers.hxx"
#include "matrix/Matrix_Triangular.hxx"
#include "matrix/Matrix_Symmetric.hxx"
#include "matrix/Matrix_Hermitian.hxx"
#include "matrix_sparse/Matrix_Sparse.hxx"
#include "matrix_sparse/Matrix_ComplexSparse.hxx"
#include "matrix_sparse/Matrix_SymSparse.hxx"
#include "matrix_sparse/Matrix_SymComplexSparse.hxx"
#include "matrix/Matrix_SymPacked.hxx"
#include "matrix/Matrix_HermPacked.hxx"
#include "matrix/Matrix_TriangPacked.hxx"
#include "vector/Vector.hxx"
#include "vector/SparseVector.hxx"
#include "matrix/Functions.hxx"
#include "matrix_sparse/Matrix_Conversions.hxx"
#include "computation/basic_functions/Functions_Matrix.hxx"
#include "computation/basic_functions/Functions_Vector.hxx"
#include "computation/basic_functions/Functions_MatVect.hxx"

#include "matrix/SubMatrix_Base.hxx"
#include "matrix/SubMatrix.hxx"

// Lapack interface.
#ifdef SELDON_WITH_LAPACK
#undef LAPACK_INTEGER
#define LAPACK_INTEGER int
#undef LAPACK_REAL
#define LAPACK_REAL float
#undef LAPACK_DOUBLEREAL
#define LAPACK_DOUBLEREAL double
#undef LAPACK_COMPLEX
#define LAPACK_COMPLEX void
#undef LAPACK_DOUBLECOMPLEX
#define LAPACK_DOUBLECOMPLEX void
#undef LAPACK_LOGICAL
#define LAPACK_LOGICAL int
#undef LAPACK_L_FP
#define LAPACK_L_FP int*
#undef LAPACK_FTNLEN
#define LAPACK_FTNLEN int*
extern "C"
{
#include "computation/interfaces/clapack.h"
}
#ifdef SELDON_LAPACK_CHECK_INFO
#ifndef SELDON_CHECK_INFO
#define SELDON_CHECK_INFO(f, lf) info.Check(f, lf)
#endif
#else
#ifndef SELDON_CHECK_INFO
#define SELDON_CHECK_INFO(f, lf)
#endif
#endif
#endif // SELDON_WITH_LAPACK.

namespace Seldon
{


  typedef Vector<int, VectFull, SELDON_DEFAULT_ALLOCATOR<int> > IVect;
  typedef Vector<float, VectFull, SELDON_DEFAULT_ALLOCATOR<float> > SVect;
  typedef Vector<double, VectFull, SELDON_DEFAULT_ALLOCATOR<double> > DVect;
  typedef Vector<complex<float>, VectFull,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CVect;
  typedef Vector<complex<double>, VectFull,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZVect;

  typedef Matrix<int, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGCMat;
  typedef Matrix<float, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGCMat;
  typedef Matrix<double, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGCMat;
  typedef Matrix<complex<float>, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGCMat;
  typedef Matrix<complex<double>, General, ColMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGCMat;

  typedef Matrix<int, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGRMat;
  typedef Matrix<float, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGRMat;
  typedef Matrix<double, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGRMat;
  typedef Matrix<complex<float>, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGRMat;
  typedef Matrix<complex<double>, General, RowMajor,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGRMat;

  typedef Matrix<int, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGRSMat;
  typedef Matrix<float, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGRSMat;
  typedef Matrix<double, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGRSMat;
  typedef Matrix<complex<float>, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGRSMat;
  typedef Matrix<complex<double>, General, RowSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGRSMat;

  typedef Matrix<int, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<int> > IGCSMat;
  typedef Matrix<float, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<float> > SGCSMat;
  typedef Matrix<double, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<double> > DGCSMat;
  typedef Matrix<complex<float>, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<float> > > CGCSMat;
  typedef Matrix<complex<double>, General, ColSparse,
		 SELDON_DEFAULT_ALLOCATOR<complex<double> > > ZGCSMat;


} // namespace Seldon.

#define SELDON_FILE_SELDONHEADER_HXX
#endif
