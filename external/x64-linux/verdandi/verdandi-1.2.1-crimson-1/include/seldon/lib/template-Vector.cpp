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
#include "share/Errors.hxx"

namespace Seldon
{
  const int& max(const int& a, const int& b)
  {
    return std::max(a, b);
  }
  const float& max(const float& a, const float& b)
  {
    return std::max(a, b);
  }
  const double& max(const double& a, const double& b)
  {
    return std::max(a, b);
  }
  const complex<float>& max(const complex<float>& a, const complex<float>& b)
  {
    throw Undefined("max(complex<float>, complex<float>)");
  }
  const complex<double>& max(const complex<double>& a, const complex<double>& b)
  {
    throw Undefined("max(complex<float>, complex<float>)");
  }
}

#include "share/Allocator.cxx"
#include "vector/Vector.cxx"
#endif


#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;


namespace Seldon
{

  SELDON_EXTERN template class MallocAlloc<@scalar>;
  SELDON_EXTERN template class Vector_Base<@scalar, MallocAlloc<@scalar> >;
  SELDON_EXTERN template class Vector<@scalar, VectFull, MallocAlloc<@scalar> >;

  // Function templates.
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::SetData(const Vector<@scalar, VectFull, MallocAlloc<@scalar> >&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::PushBack(const @scalar&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::PushBack(const Vector<@scalar, VectFull, MallocAlloc<@scalar> >&);
  SELDON_EXTERN template void Vector<@scalar, VectFull, MallocAlloc<@scalar> >::Fill(const @scalar&);
#ifndef SWIG
  SELDON_EXTERN template Vector<@scalar, VectFull, MallocAlloc<@scalar> >& Vector<@scalar, VectFull, MallocAlloc<@scalar> >::operator= (const @scalar&);
#endif
  SELDON_EXTERN template Vector<@scalar, VectFull, MallocAlloc<@scalar> >& Vector<@scalar, VectFull, MallocAlloc<@scalar> >::operator*= (const @scalar&);

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream& out, const Vector<@scalar, VectFull, MallocAlloc<@scalar> >& V);
#endif

}
