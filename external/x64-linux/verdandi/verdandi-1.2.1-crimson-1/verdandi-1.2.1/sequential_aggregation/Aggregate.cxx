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


#ifndef VERDANDI_FILE_SEQUENTIAL_AGGREGATION_AGGREGATE_CXX
#define VERDANDI_FILE_SEQUENTIAL_AGGREGATION_AGGREGATE_CXX


#include "Aggregate.hxx"


namespace Verdandi
{


    //! Aggregates ensemble data.
    /*! The aggregation is carried out with a given weight vector. Every
      member is given a single weight, independent of time and location.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] weight aggregation weights, indexed by member.
      \param[out] aggregated the aggregated output, indexed by time and
      location.
    */
    template <class Vector3, class T, class Vector2>
    void Aggregate(const Vector3& ensemble, const Vector<T>& weight,
                   Vector2& aggregated)
    {
        int Ns = ensemble.GetLength();

        if (Ns == 0)
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Vector<T>& weight, "
                                "Vector2& aggregated)",
                                "The ensemble contains no member.");
        if (Ns != weight.GetLength())
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Vector<T>& weight, "
                                "Vector2& aggregated)",
                                "The ensemble contains " + to_str(Ns)
                                + " member(s) while there are "
                                + to_str(weight.GetLength()) + " weight(s).");

        int Nt = ensemble.GetLength(0);
        aggregated.Reallocate(Nt);

        int Nl;
        int t, l, s;
        T aggregated_value;
        for (t = 0; t < Nt; t++)
        {
            Nl = ensemble.GetLength(0, t);
            aggregated(t).Reallocate(Nl);
            for (l = 0; l < Nl; l++)
            {
                aggregated_value = T(0);
                for (s = 0; s < Ns; s++)
                    aggregated_value += ensemble(s, t, l) * weight(s);
                aggregated(t, l) = aggregated_value;
            }
        }
    }


    //! Aggregates ensemble data.
    /*! The aggregation is carried out with a given weight vector. Every
      member is given a single weight, independent of time and location. The
      output aggregated values are kept only for the time indexes \a begin to
      \a end - 1.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] weight aggregation weights, indexed by member.
      \param[in] begin inclusive lower-bound for the time indexes.
      \param[in] end exclusive upper-bound for the time indexes.
      \param[out] aggregated the aggregated output, indexed by time and
      location, but only for the time period [\a beg, \a end[.
    */
    template <class Vector3, class T, class Vector2>
    void Aggregate(const Vector3& ensemble, const Vector<T>& weight,
                   int begin, int end, Vector2& aggregated)
    {
        int Ns = ensemble.GetLength();

        if (Ns == 0)
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Vector<T>& weight, "
                                "int begin, int end, Vector2& aggregated)",
                                "The ensemble contains no member.");
        if (Ns != weight.GetLength())
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Vector<T>& weight, "
                                "int begin, int end, Vector2& aggregated)",
                                "The ensemble contains " + to_str(Ns)
                                + " member(s) while there are "
                                + to_str(weight.GetLength()) + " weight(s).");
        if (begin < 0 || end > ensemble.GetLength(0) || begin > end)
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Vector<T>& weight, "
                                "int begin, int end, Vector2& aggregated)",
                                "The time period should be included in [0, "
                                + to_str(ensemble.GetLength(0))
                                + "[ but [" + to_str(begin) + ", "
                                + to_str(end) + "[ was provided.");

        Vector<int> shape(end - begin);

        for (int i = begin; i < end; i++)
            shape(i - begin) = ensemble(0, i).GetLength();
        aggregated.Reallocate(shape);

        int Nl;
        int t, l, s;
        T aggregated_value;
        for (t = begin; t < end; t++)
            for (l = 0; l < aggregated.GetLength(t - begin); l++)
            {
                aggregated_value = T(0);
                for (s = 0; s < Ns; s++)
                    aggregated_value += ensemble(s, t, l) * weight(s);
                aggregated(t - begin, l) = aggregated_value;
            }
    }


    //! Aggregates ensemble data.
    /*! The aggregation is carried out with a given weight vector. Every
      member is given a single weight per time step.
      \param[in] ensemble the ensemble data indexed by member, time and
      location.
      \param[in] weight aggregation weights, indexed by time and member.
      \param[out] aggregated the aggregated output, indexed by time and
      location.
    */
    template <class Vector3, class T, class Vector2>
    void Aggregate(const Vector3& ensemble, const Matrix<T>& weight,
                   Vector2& aggregated)
    {
        int Ns = ensemble.GetLength();

        int Nt = ensemble.GetLength(0);

        if (Ns == 0)
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Matrix<T>& weight, "
                                "Vector2& aggregated)",
                                "The ensemble contains no member.");
        if (Ns != weight.GetN())
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Matrix<T>& weight, "
                                "Vector2& aggregated)",
                                "The ensemble contains " + to_str(Ns)
                                + " member(s) while there are "
                                + to_str(weight.GetN())
                                + " weight(s) per time step.");
        if (Nt != weight.GetM())
            throw ErrorArgument("Aggregate(const Vector3& ensemble, "
                                "const Matrix<T>& weight, "
                                "Vector2& aggregated)",
                                "The ensemble contains data for "
                                + to_str(Nt)
                                + " time step(s) while there are weights for "
                                + to_str(weight.GetM()) + " time step(s).");

        aggregated.Reallocate(Nt);

        int Nl;
        int t, l, s;
        T aggregated_value;
        for (t = 0; t < Nt; t++)
        {
            Nl = ensemble.GetLength(0, t);
            aggregated(t).Reallocate(Nl);
            for (l = 0; l < Nl; l++)
            {
                aggregated_value = T(0);
                for (s = 0; s < Ns; s++)
                    aggregated_value += ensemble(s, t, l) * weight(t, s);
                aggregated(t, l) = aggregated_value;
            }
        }
    }


    /*! \brief Returns in a vector all values from a range of inner vectors of
      \a Vin. */
    /*! The output vector \a Vout contains all inner vectors of \a Vin, in the
      index range [\a beg, \a end[, concatenated in the same order as they
      appear in \a Vin.
      \param[out] Vin Vector2 instance from which the data is extracted.
      \param[in] begin inclusive lower-bound for the indexes.
      \param[in] end exclusive upper-bound for the indexes.
      \param[out] Vout the values contained in the inner vectors [\a beg, \a
      end[.
    */
    template <class Vector2, class T>
    void Flatten(const Vector2& Vin, int begin, int end, Vector<T>& Vout)
    {
        Vin.Flatten(begin, end, Vout);
    }


    /*! \brief Returns in a matrix all values from a range of vectors of \a
      Vin. */
    /*! Every row i the output matrix \a Mout contains all inner vectors of
      the i-th vector of vectors of \a Vin, in the index range [\a beg, \a
      end[, concatenated in the same order as they appear in the inner i-th
      vector of vectors of \a Vin. If \a Vin is a 3D array with shape (N0, N1,
      N2), the output matrix \a Mout is of shape (N0, (end - begin) ? N2).
      \param[in] Vin Vector3 instance from which the data is extracted. The
      total number of elements in every selected vector of vectors should be
      the same, since it is the number of columns in the output matrix.
      \param[in] begin inclusive lower-bound for the indexes.
      \param[in] end exclusive upper-bound for the indexes.
      \param[out] Mout the values contained in the inner vectors [\a beg, \a
      end[.
    */
    template <class Vector3, class T>
    void Flatten(const Vector3& Vin, int begin, int end, Matrix<T>& Mout)
    {
        int Ns = Vin.GetLength();
        Vector<T> Vout;
        for (int s = 0; s < Ns; s++)
        {
            // Number of elements to extract.
            int total = 0;
            for (int h = begin; h < end; h++)
                total += Vin.GetLength(s, h);
            Vout.Reallocate(total);
            int i = 0;
            for (int h = begin; h < end; h++)
                for (int l = 0; l < Vin.GetLength(s, h); l++)
                    Vout(i++) = Vin(s, h, l);

            if (s == 0)
                Mout.Reallocate(Ns, total);
            SetRow(Vout, s, Mout);
        }
    }


}


#endif
