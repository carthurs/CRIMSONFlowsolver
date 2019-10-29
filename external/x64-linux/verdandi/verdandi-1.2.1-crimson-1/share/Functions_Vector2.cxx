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


#ifndef VERDANDI_FILE_SHARE_FUNCTIONS_VECTOR2_CXX


#include "Functions_Vector2.hxx"


namespace Verdandi
{


    //! Returns the minimum value.
    /*!
      \param[in] V the Vector2 instance in which the minimum value is sought.
      \return The minimum scalar value found in the inner vectors of \a V.
    */
    template <class T, class Allocator0, class Allocator1>
    T Minimum(const Vector2<T, Allocator0, Allocator1>& V)
    {
        if (V.IsEmpty())
            throw WrongArgument("Minimum(Vector2& V)",
                                "The Vector2 instance is empty.");
        T minimum = numeric_limits<double>::max();
        int i, j;
        for (i = 0; i < V.GetLength(); i++)
            for (j = 0; j < V.GetLength(i); j++)
                minimum = min(minimum, V(i, j));
        return minimum;
    }


    //! Returns the maximum value.
    /*!
      \param[in] V the Vector2 instance in which the maximum value is sought.
      \return The maximum scalar value found in the inner vectors of \a V.
    */
    template <class T, class Allocator0, class Allocator1>
    T Maximum(const Vector2<T, Allocator0, Allocator1>& V)
    {
        if (V.IsEmpty())
            throw WrongArgument("Maximum(Vector2& V)",
                                "The Vector2 instance is empty.");
        T maximum = numeric_limits<double>::min();
        int i, j;
        for (i = 0; i < V.GetLength(); i++)
            for (j = 0; j < V.GetLength(i); j++)
                maximum = max(maximum, V(i, j));
        return maximum;
    }


    //! Returns the mean value.
    /*!
      \param[in] V the Vector2 instance whose mean is computed.
      \return The mean of the scalar values of the inner vectors of \a V.
    */
    template <class T, class Allocator0, class Allocator1>
    T Mean(const Vector2<T, Allocator0, Allocator1>& V)
    {
        if (V.IsEmpty())
            throw WrongArgument("Mean(Vector2& V)",
                                "The Vector2 instance is empty.");
        T total = 0;
        int i, j, N(0);
        for (i = 0; i < V.GetLength(); i++)
            for (j = 0; j < V.GetLength(i); j++, N++)
                total += V(i, j);
        return total / T(N);
    }


    //! Removes data equal to a given value.
    /*! Any data equal to \a value is removed from \a V.
      \param[in] value the value to be removed.
      \param[in,out] V instance of Vector2 in which the values equal to \a
      value are removed.
    */
    template <class T,
              class TV, class Allocator0, class Allocator1>
    void RemoveData(T value, Vector2<TV, Allocator0, Allocator1>& V)
    {
        int i, j, new_j;
        for (i = 0; i < V.GetLength(); i++)
        {
            for (new_j = j = 0; j < V.GetLength(i); j++)
                if (V(i, j) != value)
                    if (new_j == j)
                        new_j++;
                    else
                    {
                        V(i, new_j) = V(i, j);
                        new_j++;
                    }
            V.Reallocate(i, new_j);
        }
    }


    //! Removes data equal to a given value.
    /*! Any data equal to \a value is removed from \a V. The elements located
      at the same positions in \a W are removed as well (whatever their
      values).
      \param[in] value the value to be removed.
      \param[in,out] V instance of Vector2 in which the values equal to \a
      value are removed.
      \param[in,out] W instance of Vector2 which will undergo the same
      removals as in \a V.
    */
    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W>
    void RemoveData(T value, Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W)
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(V, W, "RemoveData(T value, Vector2& V, Vector2& W)");
#endif

        int i, j, new_j;
        for (i = 0; i < V.GetLength(); i++)
        {
            for (new_j = j = 0; j < V.GetLength(i); j++)
                if (V(i, j) != value)
                {
                    if (new_j == j)
                        new_j++;
                    else
                    {
                        V(i, new_j) = V(i, j);
                        W(i, new_j) = W(i, j);
                        new_j++;
                    }
                }
            V.Reallocate(i, new_j);
            W.Reallocate(i, new_j);
        }
    }


    //! Removes data equal to a given value.
    /*! Any data equal to \a value is removed from \a V. The elements located
      at the same positions in \a W and \a Z are removed as well (whatever
      their values).
      \param[in] value the value to be removed.
      \param[in,out] V instance of Vector2 in which the values equal to \a
      value are removed.
      \param[in,out] W instance of Vector2 which will undergo the same
      removals as in \a V.
      \param[in,out] Z instance of Vector2 which will undergo the same
      removals as in \a V.
    */
    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W,
              class TZ, class Allocator0Z, class Allocator1Z>
    void RemoveData(T value, Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W,
                    Vector2<TZ, Allocator0Z, Allocator1Z>& Z)
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(V, W,
                   "RemoveData(T value, Vector2& V, Vector2& W, Vector2& Z)");
        CheckShape(V, Z,
                   "RemoveData(T value, Vector2& V, Vector2& W, Vector2& Z)");
#endif

        int i, j, new_j;
        for (i = 0; i < V.GetLength(); i++)
        {
            for (new_j = j = 0; j < V.GetLength(i); j++)
                if (V(i, j) != value)
                {
                    if (new_j == j)
                        new_j++;
                    else
                    {
                        V(i, new_j) = V(i, j);
                        W(i, new_j) = W(i, j);
                        Z(i, new_j) = Z(i, j);
                        new_j++;
                    }
                }
            V.Reallocate(i, new_j);
            W.Reallocate(i, new_j);
            Z.Reallocate(i, new_j);
        }
    }


    //! Removes the values that are in a given interval.
    /*! Any data greater than \a value_min and strictly lower than \a
      value_max is removed from \a V. The elements located at the same
      positions in \a W are removed as well (whatever their values).
      \param[in] value_min the lower bound of the interval (inclusive).
      \param[in] value_max the upper bound of the interval (exclusive).
      \param[in,out] V instance of Vector2 in which the values greater than \a
      value_min and strictly lower than \a value_max are removed.
      \param[in,out] W instance of Vector2 which will undergo the same
      removals as in \a V.
    */
    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W>
    void RemoveData(T value_min, T value_max,
                    Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W)
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(V, W, "RemoveData(T value_min, T value_max, "
                   "Vector2& V, Vector2& W)");
#endif

        int i, j, new_j;
        for (i = 0; i < V.GetLength(); i++)
        {
            for (new_j = j = 0; j < V.GetLength(i); j++)
                if (V(i, j) < value_min || V(i, j) >= value_max)
                {
                    if (new_j == j)
                        new_j++;
                    else
                    {
                        V(i, new_j) = V(i, j);
                        W(i, new_j) = W(i, j);
                        new_j++;
                    }
                }
            V.Reallocate(i, new_j);
            W.Reallocate(i, new_j);
        }
    }


    //! Removes the values that are in a given interval..
    /*! Any data greater than \a value_min and strictly lower than \a
      value_max is removed from \a V. The elements located at the same
      positions in \a W and \a Z are removed as well (whatever their values).
      \param[in] value_min the lower bound of the interval (inclusive).
      \param[in] value_max the upper bound of the interval (exclusive).
      \param[in,out] V instance of Vector2 in which the values greater than \a
      value_min and strictly lower than \a value_max are removed.
      \param[in,out] W instance of Vector2 which will undergo the same
      removals as in \a V.
      \param[in,out] Z instance of Vector2 which will undergo the same
      removals as in \a V.
    */
    template <class T,
              class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W,
              class TZ, class Allocator0Z, class Allocator1Z>
    void RemoveData(T value_min, T value_max,
                    Vector2<TV, Allocator0V, Allocator1V>& V,
                    Vector2<TW, Allocator0W, Allocator1W>& W,
                    Vector2<TZ, Allocator0Z, Allocator1Z>& Z)
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(V, W, "RemoveData(T value_min, T value_max, "
                   "Vector2& V, Vector2& W, Vector2& Z)");
        CheckShape(V, Z, "RemoveData(T value_min, T value_max, "
                   "Vector2& V, Vector2& W, Vector2& Z)");
#endif

        int i, j, new_j;
        for (i = 0; i < V.GetLength(); i++)
        {
            for (new_j = j = 0; j < V.GetLength(i); j++)
                if (V(i, j) < value_min || V(i, j) >= value_max)
                {
                    if (new_j == j)
                        new_j++;
                    else
                    {
                        V(i, new_j) = V(i, j);
                        W(i, new_j) = W(i, j);
                        Z(i, new_j) = Z(i, j);
                        new_j++;
                    }
                }
            V.Reallocate(i, new_j);
            W.Reallocate(i, new_j);
            Z.Reallocate(i, new_j);
        }
    }


    //! Removes the locations without enough data.
    /*! This function applies to data in the form of (\a location_index, \a
      data) where \a location_index and \a data are Vector2 instances indexed
      by time and location. On exit, any location where the proportion of
      available data was strictly less than \a ratio has been removed from \a
      location_index and \a data. Also, the location indexes have been
      redefined so that there is no missing location indexes.
      \param[in,out] location_index the location indexes corresponding to \a
      data. Hence it must have the same shape on entry (and on exit) as \a
      data. On exit, the location indexes range from 0 to N-1 where N is the
      total number of remaining locations.
      \param[in,out] data Vector2 instance that contains some values at the
      locations.
      \param[in] ratio the ratio of availability (in time), strictly below
      which the locations are removed.
    */
    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d>
    void
    RemoveEmptyLocation(Vector2<int, Allocator0l, Allocator1l>&
                        location_index,
                        Vector2<Td, Allocator0d, Allocator1d>& data,
                        double ratio = numeric_limits<double>::epsilon())
    {
        int N = 0;
        if (!location_index.IsEmpty())
            N = Maximum(location_index) + 1;

        Vector<bool> keep(N);
        RemoveEmptyLocation(location_index, data, keep, ratio);
    }


    //! Removes the locations without enough data.

    /*! This function applies to data in the form of (\a location_index, \a
      data) where \a location_index and \a data are Vector2 instances indexed
      by time and location. On exit, any location where the proportion of
      available data was strictly less than \a ratio has been removed from \a
      location_index and \a data. Also, the location indexes have been
      redefined so that there is no missing location indexes.
      \param[in,out] location_index the location indexes corresponding to \a
      data. Hence it must have the same shape on entry (and on exit) as \a
      data. On exit, the location indexes range from 0 to N-1 where N is the
      total number of remaining locations.
      \param[in,out] data Vector2 instance that contains some values at the
      locations.
      \param[out] keep the vector, indexed by the locations (before any is
      removed), of the locations status: 'true' for a location that is kept,
      'false' for a location that is removed.
      \param[in] ratio the ratio of availability (in time), strictly below
      which the locations are removed.
    */
    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d>
    void
    RemoveEmptyLocation(Vector2<int, Allocator0l, Allocator1l>&
                        location_index,
                        Vector2<Td, Allocator0d, Allocator1d>& data,
                        Vector<bool>& keep,
                        double ratio = numeric_limits<double>::epsilon())
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(location_index, data, "RemoveEmptyLocation(Vector2& "
                   "location_index, Vector2& data, double ratio)");
#endif

        if (location_index.IsEmpty())
        {
            keep.Clear();
            return;
        }

        // Number of locations.
        int Nl = Maximum(location_index) + 1;
        // Number of available data per location.
        Vector<int> count(Nl);
        count.Zero();
        keep.Reallocate(Nl);
        keep.Fill(false);

        int t, l;
        int Nt = location_index.GetLength();
        for (t = 0; t < Nt; t++)
            for (l = 0; l < location_index.GetLength(t); l++)
                count(location_index(t, l))++;
        for (l = 0; l < Nl; l++)
            keep(l) = double(count(l)) / double(Nt) >= ratio;

        RemoveLocation(keep, location_index, data);
    }


    //! Removes given locations.
    /*! This function applies to data in the form of (\a location_index, \a
      data) where \a location_index and \a data are Vector2 instances indexed
      by time and location. On exit, the location indexes have been redefined
      so that there is no missing location indexes.
      \param[in] keep the vector, indexed by the locations (before any is
      removed), of the locations status: 'true' for a location that is kept,
      'false' for a location that is removed.
      \param[in,out] location_index the location indexes corresponding to \a
      data. Hence it must have the same shape on entry (and on exit) as \a
      data. On exit, the location indexes range from 0 to N-1 where N is the
      total number of remaining locations.
      \param[in,out] data Vector2 instance that contains some values at the
      locations.
    */
    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d>
    void RemoveLocation(const Vector<bool>& keep,
                        Vector2<int, Allocator0l, Allocator1l>&
                        location_index,
                        Vector2<Td, Allocator0d, Allocator1d>& data)
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(location_index, data,
                   "RemoveLocation(Vector<bool>, Vector2, Vector2)");
#endif

        int Nl = Maximum(location_index) + 1;

        if (keep.GetLength() != Nl)
            throw WrongArgument("RemoveLocation(Vector<bool>, "
                                "Vector2, Vector2)",
                                "There are " + to_str(Nl) + " locations "
                                + "but the filtering vector has "
                                + to_str(keep.GetLength()) + " elements.");

        int Nt = location_index.GetLength();

        int t, l;
        Nl = 0;
        for (l = 0; l < keep.GetLength(); l++)
            if (keep(l))
                Nl++;

        int new_l;
        for (t = 0; t < Nt; t++)
        {
            for (new_l = l = 0; l < location_index.GetLength(t); l++)
                if (keep(location_index(t, l)))
                {
                    if (new_l == l)
                        new_l++;
                    else
                    {
                        data(t, new_l) = data(t, l);
                        location_index(t, new_l) = location_index(t, l);
                        new_l++;
                    }
                }
            data.Reallocate(t, new_l);
            location_index.Reallocate(t, new_l);
        }

        // 'count(l)' is going to be the number of locations removed before
        // 'l'.
        Vector<int> count(keep.GetLength());
        int shift = 0;
        for (l = 0; l < keep.GetLength(); l++)
            if (!keep(l))
                shift++;
            else
                count(l) = l - shift;
        for (t = 0; t < Nt; t++)
            for (l = 0; l < location_index.GetLength(t); l++)
                location_index(t, l) = count(location_index(t, l));
    }


    //! Removes locations not available in \a index_ref.
    /*! It removes from \a index any location that is not in \a index_ref. The
      location indexes in \a index_ref and \a index should refer to the same
      locations.
      \param[in] index_ref the locations to be kept in \a index.
      \param[in,out] index the location indexes to be filtered.
    */
    template <class Allocator0r, class Allocator1r,
              class Allocator0i, class Allocator1i>
    void
    SelectLocation(const Vector2<int, Allocator0r, Allocator1r>& index_ref,
                   Vector2<int, Allocator0i, Allocator1i>& index)
    {
        int Nt = index.GetLength();

        if (Nt != index_ref.GetLength())
            throw WrongArgument("SelectLocation(const Vector2& index_ref, "
                                "Vector2& index)",
                                "The reference location indexes (index_ref) "
                                "have " + to_str(index_ref.GetLength())
                                + " inner vectors while the output location"
                                " indexes (index) have " + to_str(Nt)
                                + " inner vectors. There must be the same"
                                " number of inner vectors.");
        int l, ref_l;
        // 'old_l' is used to loop over the location indexes in 'index', and
        // 'new_l' will be the new index (in locations).
        int old_l, new_l, max_l;
        for (int t = 0; t < Nt; t++)
        {
            new_l = old_l = 0;
            for (l = 0; l < index_ref.GetLength(t); l++)
            {
                ref_l = index_ref(t, l);
                max_l = index.GetLength(t);
                // Searches for the location index in 'index'.
                while (old_l < max_l && index(t, old_l) < ref_l)
                    old_l++;
                if (old_l != max_l && index(t, old_l) == ref_l)
                    index(t, new_l++) = ref_l;
            }
            index.Reallocate(t, new_l);
        }
    }


    //! Removes locations not available in \a index_ref.
    /*! It removes from the pair (\a index, \a data) any location that is not
      in \a index_ref. The location indexes in \a index_ref and \a index
      should refer to the same locations.
      \param[in] index_ref the locations to be kept in the pair (\a index, \a
      data).
      \param[in,out] index the location index to be filtered.
      \param[in,out] data the data associated with \a index. It will be
      filtered consistently with \a index.
    */
    template <class Allocator0r, class Allocator1r,
              class Allocator0i, class Allocator1i,
              class Td, class Allocator0d, class Allocator1d>
    void
    SelectLocation(const Vector2<int, Allocator0r, Allocator1r>& index_ref,
                   Vector2<int, Allocator0i, Allocator1i>& index,
                   Vector2<Td, Allocator0d, Allocator1d>& data)
    {
#ifdef VERDANDI_CHECK_DIMENSIONS
        CheckShape(index, data, "RestrictLocation(const Vector2& index_ref, "
                   "Vector2& index, Vector2& data)");
#endif

        int Nt = index.GetLength();

        if (Nt != index_ref.GetLength())
            throw WrongArgument("RestrictLocation(const Vector2& index_ref, "
                                "Vector2& index, Vector2& data)",
                                "The reference location indexes (index_ref) "
                                "have " + to_str(index_ref.GetLength())
                                + " inner vectors while the output pair"
                                " (index, data) has " + to_str(Nt)
                                + " inner vectors. There must be the same"
                                " number of inner vectors.");
        int l, ref_l;
        // 'old_l' is used to loop over the location indexes in 'index', and
        // 'new_l' will be the new index (in locations) of the data for a
        // given inner vector.
        int old_l, new_l, max_l;
        for (int t = 0; t < Nt; t++)
        {
            new_l = old_l = 0;
            for (l = 0; l < index_ref.GetLength(t); l++)
            {
                ref_l = index_ref(t, l);
                max_l = index.GetLength(t);
                // Searches for the location index in 'index'.
                while (old_l < max_l && index(t, old_l) < ref_l)
                    old_l++;
                if (old_l != max_l && index(t, old_l) == ref_l)
                {
                    data(t, new_l) = data(t, old_l);
                    index(t, new_l) = ref_l;
                    new_l++;
                }
            }
            data.Reallocate(t, new_l);
            index.Reallocate(t, new_l);
        }
    }


    //! Returns in a vector all values from a Vector2 instance.
    /*! The output vector \a data contains all inner vectors concatenated in
      the same order as they appear in \a V.
      \param[in] V the Vector2 instance whose values are collected.
      \param[out] data all values from \a V.
    */
    template <class TV, class Allocator0V, class Allocator1V,
              class T, class Allocator>
    void Collect(const Vector2<TV, Allocator0V, Allocator1V>& V,
                 Vector<T, VectFull, Allocator>& data)
    {
        V.Flatten(data);
    }


    /*! \brief Returns in a vector the values of inner vectors from a Vector2
      instance. */
    /*! The output vector \a data contains the inner vectors, indexed in [\a
      beg, \a end[, concatenated in the same order as they appear in \a V.
      \param[in] V the Vector2 instance whose values are partially collected.
      \param[in] beg inclusive lower-bound for the indexes.
      \param[in] end exclusive upper-bound for the indexes.
      \param[out] data all values from the inner vectors [\a beg, \a end[ of
      \a V.
    */
    template <class TV, class Allocator0V, class Allocator1V,
              class T, class Allocator>
    void CollectDate(const Vector2<TV, Allocator0V, Allocator1V>& V,
                     int beg, int end, Vector<T, VectFull, Allocator>& data)
    {
        V.Flatten(beg, end, data);
    }


    //! Returns all available data at a given location.
    /*!
      \param[in] location_index the location indexes corresponding to \a
      data. It must have the same shape as \a data.
      \param[in] data Vector2 instance that contains some values at the
      locations.
      \param[in] location index of the selected location.
      \param[out] step steps at which the values are given: they are the
      indexes of the inner vectors from which the data was extracted.
      \param[out] output the available values for location \a location. \a
      step and \a data have the same length on exit.
    */
    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d,
              class Allocators, class To, class Allocatoro>
    void
    CollectLocation(const Vector2<int, Allocator0l, Allocator1l>&
                    location_index,
                    const Vector2<Td, Allocator0d, Allocator1d>& data,
                    int location,
                    Vector<int, Allocators>& step,
                    Vector<To, Allocatoro>& output)
    {
        int Nt = location_index.GetLength();

        step.Reallocate(Nt);
        output.Reallocate(Nt);

        int t, l, i(0);
        for (t = 0; t < Nt; t++)
        {
            l = 0;
            while (l < location_index.GetLength(t)
                   && location_index(t, l) < location)
                l++;
            if (l < location_index.GetLength(t)
                && location_index(t, l) == location)
            {
                step(i) = t;
                output(i) = data(t, l);
                i++;
            }
        }
        step.Reallocate(i);
        output.Reallocate(i);
    }


    /*! \brief Returns all available data at a given location and over a range
      of steps. */
    /*!
      \param[in] location_index the location indexes corresponding to \a
      data. It must have the same shape as \a data.
      \param[in] data Vector2 instance that contains some values at the
      locations.
      \param[in] location index of the selected location.
      \param[in] beg inclusive lower-bound for the selected steps.
      \param[in] end exclusive upper-bound for the selected steps.
      \param[out] step steps at which the values are given: they are the
      indexes of the inner vectors from which the data was extracted.
      \param[out] output the available values for location \a location, during
      the steps [\a beg, \a end[. \a step and \a data have the same length on
      exit.
    */
    template <class Allocator0l, class Allocator1l,
              class Td, class Allocator0d, class Allocator1d,
              class Allocators, class To, class Allocatoro>
    void
    CollectLocation(const Vector2<int, Allocator0l, Allocator1l>&
                    location_index,
                    const Vector2<Td, Allocator0d, Allocator1d>& data,
                    int location, int beg, int end,
                    Vector<int, Allocators>& step,
                    Vector<To, Allocatoro>& output)
    {
        if (beg > end)
            throw
                WrongArgument("CollectLocation(Vector2<int>& location_index, "
                              "Vector2& data, int location, int beg, int end,"
                              " Vector<int>& step, Vector& output)",
                              "The lower bound of the range of steps, ["
                              + to_str(beg) + ", " + to_str(end) + "[, is"
                              " strictly greater than its upper bound.");
        if (beg < 0 || end > location_index.GetLength())
            throw
                WrongArgument("CollectLocation(Vector2<int>& location_index, "
                              "Vector2& data, int location, int beg, int end,"
                              " Vector<int>& step, Vector& output)",
                              "The step indexes should be in [0,"
                              + to_str(location_index.GetLength()) + "] but ["
                              + to_str(beg) + ", " + to_str(end)
                              + "[ was provided.");

        step.Reallocate(end - beg);
        output.Reallocate(end - beg);

        int t, l, i(0);
        for (t = beg; t < end; t++)
        {
            l = 0;
            while (l < location_index.GetLength(t)
                   && location_index(t, l) < location)
                l++;
            if (l < location_index.GetLength(t)
                && location_index(t, l) == location)
            {
                step(i) = t;
                output(i) = data(t, l);
                i++;
            }
        }
        step.Reallocate(i);
        output.Reallocate(i);
    }


    /*! \brief Raises an exception in case two Vector2 instances do not have
      the same shape. */
    /*! This function is supposed to be called in another function that needs
      to check the consistency between two Vector2 instances before
      proceeding. It the shapes of \a V and \a W are not the same,
      \param[in] V instance of Vector2.
      \param[in] W instance of Vector2.
      \param[in] function the name of the function where the shapes are
      compared.
    */
    template <class TV, class Allocator0V, class Allocator1V,
              class TW, class Allocator0W, class Allocator1W>
    void CheckShape(const Vector2<TV, Allocator0V, Allocator1V>& V,
                    const Vector2<TW, Allocator0W, Allocator1W>& W,
                    string function)
    {
        if (V.HasSameShape(W))
            return;

        int i;
        ostringstream message;
        message << "Inconsistent shapes for two Vector2 instances.\n"
                << "The first Vector2 has " << V.GetLength()
                << " inner vector(s)";
        if (V.GetLength() == 0)
            message << ".\n";
        else
        {
            message << ", with size(s):\n";
            for (i = 0; i < V.GetLength() - 1; i++)
                message << V.GetLength(i) << " ";
            message << V.GetLength(V.GetLength() - 1) << "\n";
        }
        message << "The second Vector2 has " << W.GetLength()
                << " inner vector(s)";
        if (W.GetLength() == 0)
            message << ".\n";
        else
        {
            message << ", with size(s):\n";
            for (i = 0; i < W.GetLength() - 1; i++)
                message << W.GetLength(i) << " ";
            message << W.GetLength(W.GetLength() - 1) << "\n";
        }
        throw WrongDim(function, message.str());
    }


} // namespace Verdandi


#define VERDANDI_FILE_SHARE_FUNCTIONS_VECTOR2_CXX
#endif
