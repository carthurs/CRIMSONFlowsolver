// Copyright (C) 2010 Vivien Mallet
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


#include "SeldonHeader.hxx"


#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/Storage.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_HermPacked.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_TriangPacked.cxx"
#endif


#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;


namespace Seldon
{

  SELDON_EXTERN template class Matrix_HermPacked<@complex, General, @storage_blasHP, MallocAlloc<@complex> >;
  SELDON_EXTERN template class Matrix<@complex, General, @storage_blasHP, MallocAlloc<@complex> >;
  SELDON_EXTERN template class Matrix_SymPacked<@real_complex, General, @storage_blasSP, MallocAlloc<@real_complex> >;
  SELDON_EXTERN template class Matrix<@real_complex, General, @storage_blasSP, MallocAlloc<@real_complex> >;
  SELDON_EXTERN template class Matrix_TriangPacked<@real_complex, General, @storage_blasTP, MallocAlloc<@real_complex> >;
  SELDON_EXTERN template class Matrix<@real_complex, General, @storage_blasTP, MallocAlloc<@real_complex> >;

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<@complex, General, @storage_blasHP, MallocAlloc<@complex> >&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<@real_complex, General, @storage_blasSP, MallocAlloc<@real_complex> >&);
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<@real_complex, General, @storage_blasTP, MallocAlloc<@real_complex> >&);
#endif

}
