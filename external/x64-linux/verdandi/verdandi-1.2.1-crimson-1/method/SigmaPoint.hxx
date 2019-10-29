// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_METHOD_SIGMA_POINT_HXX

#include "seldon/vector/VectorCollection.hxx"

namespace Verdandi
{


    template <class SigmaPoint>
    void ComputeCanonicalSigmaPoint(int p,
                                    Vector<SigmaPoint, Collection>&
                                    sigma_point, SigmaPoint& alpha,
                                    bool& alpha_constant);
    template <class T,  template <class U> class Allocator>
    void ComputeCanonicalSigmaPoint(int p,
                                    Matrix<T, General, RowMajor,
                                    Allocator<T> >& sigma_point,
                                    Vector<T, VectFull, Allocator<T> >&
                                    alpha, bool& alpha_constant);
    template <class T,  template <class U> class Allocator>
    void ComputeCanonicalSigmaPoint(int p,
                                    Matrix<T, General, RowSparse,
                                    Allocator<T> >& sigma_point,
                                    Vector<T, VectFull, Allocator<T> >&
                                    alpha, bool& alpha_constant);

    template <class SigmaPoint>
    void ComputeStarSigmaPoint(int p,
                               Vector<SigmaPoint, Collection>& sigma_point,
                               SigmaPoint& alpha, bool& alpha_constant);
    template <class T,  template <class U> class Allocator>
    void ComputeStarSigmaPoint(int p,
                               Matrix<T, General, RowMajor, Allocator<T> >&
                               sigma_point,
                               Vector<T, VectFull, Allocator<T> >& alpha,
                               bool& alpha_constant);
    template <class T,  template <class U> class Allocator>
    void ComputeStarSigmaPoint(int p,
                               Matrix<T, General, RowSparse, Allocator<T> >&
                               sigma_point,
                               Vector<T, VectFull, Allocator<T> >& alpha,
                               bool& alpha_constant);
    template <class SigmaPoint>
    void ComputeSimplexSigmaPoint(int p,
                                  Vector<SigmaPoint, Collection>& sigma_point,
                                  SigmaPoint& alpha, bool& alpha_constant);
    template <class T,  template <class U> class Allocator>
    void ComputeSimplexSigmaPoint(int p,
                                  Matrix<T, General, RowMajor, Allocator<T> >&
                                  sigma_point,
                                  Vector<T, VectFull, Allocator<T> >& alpha,
                                  bool& alpha_constant);
    template <class T,  template <class U> class Allocator>
    void ComputeSimplexSigmaPoint(int p,
                                  Matrix<T, General, RowSparse, Allocator<T> >
                                  & sigma_point,
                                  Vector<T, VectFull, Allocator<T> >& alpha,
                                  bool& alpha_constant);

} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_SIGMA_POINT_HXX
#endif
