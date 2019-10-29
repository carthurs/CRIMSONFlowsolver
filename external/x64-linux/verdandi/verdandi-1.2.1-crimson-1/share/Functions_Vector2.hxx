// Copyright (C) 2010-2012 INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_SHARE_FUNCTIONS_VECTOR2_HXX


namespace Verdandi
{


    template <class T, class Allocator0, class Allocator1>
    T Minimum(const Vector2<T, Allocator0, Allocator1>& V);


    template <class T, class Allocator0, class Allocator1>
    T Maximum(const Vector2<T, Allocator0, Allocator1>& V);


    template <class T, class Allocator0, class Allocator1>
    T Mean(const Vector2<T, Allocator0, Allocator1>& V);


    template <class T,
              class TV, class Allocator0, class Allocator1>
    void RemoveData(T value, Vector2<TV, Allocator0, Allocator1>& V);


    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W>
    void RemoveData(T value, Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W);


    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W,
              class TZ, class Allocator0Z, class Allocator1Z>
    void RemoveData(T value, Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W,
                    Vector2<TZ, Allocator0Z, Allocator1Z>& Z);


    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W>
    void RemoveData(T value_min, T value_max,
                    Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W);


    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W,
              class TZ, class Allocator0Z, class Allocator1Z>
    void RemoveData(T value_min, T value_max,
                    Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W,
                    Vector2<TZ, Allocator0Z, Allocator1Z>& Z);


    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d>
    void
    RemoveEmptyLocation(Vector2<int, Allocator0l, Allocator1l>&
                        location_index,
                        Vector2<Td, Allocator0d, Allocator1d>& data,
                        double ratio = numeric_limits<double>::epsilon());


    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d>
    void
    RemoveEmptyLocation(Vector2<int, Allocator0l, Allocator1l>&
                        location_index,
                        Vector2<Td, Allocator0d, Allocator1d>& data,
                        Vector<bool>& keep,
                        double ratio = numeric_limits<double>::epsilon());


    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d>
    void RemoveLocation(const Vector<bool>& keep,
                        Vector2<int, Allocator0l, Allocator1l>&
                        location_index,
                        Vector2<Td, Allocator0d, Allocator1d>& data);


    template <class Allocator0r, class Allocator1r,
              class Allocator0i, class Allocator1i>
    void
    SelectLocation(const Vector2<int, Allocator0r, Allocator1r>& index_ref,
                   Vector2<int, Allocator0i, Allocator1i>& index);


    template <class Allocator0r, class Allocator1r,
              class Allocator0i, class Allocator1i,
              class Td, class Allocator0d, class Allocator1d>
    void
    SelectLocation(const Vector2<int, Allocator0r, Allocator1r>& index_ref,
                   Vector2<int, Allocator0i, Allocator1i>& index,
                   Vector2<Td, Allocator0d, Allocator1d>& data);


    template <class TV, class Allocator0V, class Allocator1V,
              class T, class Allocator>
    void Collect(const Vector2<TV, Allocator0V, Allocator1V>& V,
                 Vector<T, VectFull, Allocator>& data);


    template <class TV, class Allocator0V, class Allocator1V,
              class T, class Allocator>
    void CollectDate(const Vector2<TV, Allocator0V, Allocator1V>& V,
                     int beg, int end, Vector<T, VectFull, Allocator>& data);


    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d,
              class Allocators, class To, class Allocatoro>
    void
    CollectLocation(const Vector2<int, Allocator0l, Allocator1l>&
                    location_index,
                    const Vector2<Td, Allocator0d, Allocator1d>& data,
                    int location,
                    Vector<int, Allocators>& step,
                    Vector<To, Allocatoro>& output);


    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d,
              class Allocators, class To, class Allocatoro>
    void
    CollectLocation(const Vector2<int, Allocator0l, Allocator1l>&
                    location_index,
                    const Vector2<Td, Allocator0d, Allocator1d>& data,
                    int location, int beg, int end,
                    Vector<int, Allocators>& step,
                    Vector<To, Allocatoro>& output);


    template <class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W>
    void CheckShape(const Vector2<TV, Allocator0V, Allocator1V>& V,
                    const Vector2<TW, Allocator0W, Allocator1W>& W,
                    string function = "");


} // namespace Verdandi


#define VERDANDI_FILE_SHARE_FUNCTIONS_VECTOR2_HXX
#endif
