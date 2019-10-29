// Copyright (C) 2008, INRIA
// Author(s): Édouard Debry, Vivien Mallet
//
// This file is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This file is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.


#ifndef VERDANDI_FILE_SEQUENTIAL_AGGREGATION_AGGREGATE_HXX
#define VERDANDI_FILE_SEQUENTIAL_AGGREGATION_AGGREGATE_HXX


namespace Verdandi
{

    template <class Vector3, class T, class Vector2>
    void Aggregate(const Vector3& ensemble, const Vector<T>& weight,
                   Vector2& aggregated);

    template <class Vector3, class T, class Vector2>
    void Aggregate(const Vector3& ensemble, const Vector<T>& weight,
                   int begin, int end, Vector2& aggregated);

    template <class Vector3, class T, class Vector2>
    void Aggregate(const Vector3& ensemble, const Matrix<T>& weight,
                   Vector2& aggregated);

    template <class Vector2, class T>
    void Flatten(const Vector2& Vin, int begin, int end, Vector<T>& Vout);

    template <class Vector3, class T>
    void Flatten(const Vector3& Vin, int begin, int end, Matrix<T>& Vout);

}


#endif
