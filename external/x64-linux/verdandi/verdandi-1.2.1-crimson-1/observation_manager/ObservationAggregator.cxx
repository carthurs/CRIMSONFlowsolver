// Copyright (C) 2008-2009 INRIA
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_OSERVATIONAGGREGATOR_CXX


#include "ObservationAggregator.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    template <class T>
    ObservationAggregator<T>::ObservationAggregator()
    {
    }


    //! Constructor.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    ObservationAggregator<T>::ObservationAggregator(string configuration_file)
    {
    }


    //! Destructor.
    template <class T>
    ObservationAggregator<T>::~ObservationAggregator()
    {
    }


    //! Initializer.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ObservationAggregator<T>::Initialize(string configuration_file)
    {
        VerdandiOps configuration(configuration_file);

        configuration.SetPrefix("observation.aggregator.");

        string interpolation_type_str;
        configuration.Set("type",
                          "ops_in(v, {'step', 'triangle', 'interpolation'})",
                          interpolation_type_str);
        if (interpolation_type_str == "step")
            interpolation_type_ = type_step_;
        else if (interpolation_type_str == "triangle")
            interpolation_type_ = type_triangle_;
        else if (interpolation_type_str == "interpolation")
            interpolation_type_ = type_interpolation_;

        string width_property_str;
        configuration.Set("width_property", width_property_str);
        if (width_property_str == "constant")
            width_property_ = width_constant_;
        else
            width_property_ = width_per_observation_;
        configuration.Set("width_left_upper_bound",
                          width_left_upper_bound_);
        configuration.Set("width_right_upper_bound",
                          width_right_upper_bound_);

        configuration.Set("width_left", "", 0., width_left_);
        configuration.Set("width_right", "", 0., width_right_);
        configuration.Set("discard_observation", discard_observation_);

        active_track_index_ = CreateTrack();
    }


    //! Activates or deactivates the option 'discard_observation'.
    /*!
      \param[in] discard_observation if set to true, each observation will be
      used at most one time.
    */
    template <class T>
    void ObservationAggregator<T>
    ::DiscardObservation(bool discard_observation)
    {
        discard_observation_ = discard_observation;
    }


    //! Returns the contribution time interval corresponding to a given time.
    /*! This method returns the time interval into which observations have a
      non-zero contribution at time \a time. An integer is associated to this
      interval to indicate the observation selection policy.
      '0' indicates that all observations available in the given interval
      have to be considered.
      '-1' indicates that all observations available in the interval
      [\a time_inf ; \a time] have to be considered.
      '1' indicates that all observations available in the interval
      [\a time ; \a time_sup[ have to be considered.
      '2' indicates that only the closest left observation of the interval
      from time and the closest right observation are requested.
      '3' indicates that all observation in the given interval have to be
      considered, but, one should take into account non constant triangle
      widths.
      \param[in] time a given time.
      \param[out] time_inf lower bound of the time interval.
      \param[out] time_sup upper bound (excluded) of the time interval.
      \param[out] selection_policy interval selection policy.
    */
    template <class T>
    void ObservationAggregator<T>
    ::GetContributionInterval(double time, double& time_inf, double& time_sup,
                              int& selection_policy)
        const
    {
        if (interpolation_type_ == type_step_)
        {
            selection_policy = 0;

            // Always switch off observation.
            if (discard_observation_ && time_.GetSize() != 0)
            {
                int i_time_inf, i_time_sup;

                GetValueIndex(time_(active_track_index_), time, i_time_inf,
                              i_time_sup);

                if (i_time_inf == -1)
                {
                    time_inf = time - width_left_;
                    time_inf = time_inf > 0. ? time_inf : 0.;

                    time_sup = time + width_right_ <
                        time_(active_track_index_)(0)
                        - width_left_ ? time + width_right_ :
                        time_(active_track_index_)(0) - width_left_;

                    return;
                }

                if (i_time_sup == -1)
                {
                    time_inf = time_(active_track_index_)(i_time_inf);
                    time_inf = time - width_left_ > time_inf + width_right_ ?
                        time - width_left_ : time_inf + width_right_;
                    time_inf = time_inf > 0. ? time_inf : 0.;

                    time_sup = time + width_right_;
                    return;
                }

                time_inf = time_(active_track_index_)(i_time_inf);
                time_inf = time - width_left_ > time_inf + width_right_ ?
                    time - width_left_ : time_inf + width_right_;
                time_inf = time_inf > 0. ? time_inf : 0.;

                time_sup = time_(active_track_index_)(i_time_sup);
                time_sup = time + width_right_ < time_sup - width_left_ ?
                    time + width_right_ : time_sup - width_left_;

                return;
            }

            time_inf = time - width_left_;
            time_inf = time_inf > 0. ? time_inf : 0.;

            time_sup = time + width_right_;

            return;
        }

        if (interpolation_type_ == type_triangle_
            && width_property_ == width_constant_)
        {
            selection_policy = 0;

            // Always switch off observation.
            if (discard_observation_ && time_.GetSize() != 0)
            {
                int i_time_inf, i_time_sup;

                GetValueIndex(time_(active_track_index_), time, i_time_inf,
                              i_time_sup);

                if (i_time_inf == -1)
                {
                    time_inf = time - width_left_;
                    time_inf = time_inf > 0. ? time_inf : 0.;

                    time_sup = time + width_right_ <
                        time_(active_track_index_)(0)
                        - width_left_ ? time + width_right_ :
                        time_(active_track_index_)(0) - width_left_;

                    return;
                }

                if (i_time_sup == -1)
                {
                    time_inf = time_(active_track_index_)(i_time_inf);
                    time_inf = time - width_left_ > time_inf + width_right_ ?
                        time - width_left_ : time_inf + width_right_;
                    time_inf = time_inf > 0. ? time_inf : 0.;

                    time_sup = time + width_right_;

                    return;
                }

                time_inf = time_(active_track_index_)(i_time_inf);
                time_inf = time - width_left_ > time_inf + width_right_ ?
                    time - width_left_ : time_inf + width_right_;
                time_inf = time_inf > 0. ? time_inf : 0.;

                time_sup = time_(active_track_index_)(i_time_sup);
                time_sup = time + width_right_ < time_sup - width_left_ ?
                    time + width_right_ : time_sup - width_left_;

                return;
            }

            time_inf = time - width_left_;
            time_inf = time_inf > 0. ? time_inf : 0.;

            time_sup = time + width_right_;

            return;
        }

        if (interpolation_type_ == type_triangle_
            && width_property_ == width_per_observation_)
        {
            selection_policy = 3;

            time_inf = time - width_left_upper_bound_;
            time_inf = time_inf > 0. ? time_inf : 0.;

            time_sup = time + width_right_upper_bound_;

            return;
        }


        if (interpolation_type_ == type_interpolation_)
        {
            selection_policy = 2;

            // Always switch off observation.
            if (discard_observation_ && time_.GetSize() != 0)
            {
                int i_time_inf, i_time_sup;

                GetValueIndex(time_(active_track_index_), time, i_time_inf,
                              i_time_sup);

                if (i_time_inf == -1)
                {
                    time_inf = time - width_left_upper_bound_;
                    time_inf = time_inf > 0. ? time_inf : 0.;

                    time_sup = time + width_right_upper_bound_ <
                        time_(active_track_index_)(0)
                        - width_left_upper_bound_ ? time +
                        width_right_upper_bound_ :
                        time_(active_track_index_)(0) -
                        width_left_upper_bound_;

                    return;
                }

                if (i_time_sup == -1)
                {
                    time_inf = time_(active_track_index_)(i_time_inf);
                    time_inf = time - width_left_upper_bound_ > time_inf +
                        width_right_upper_bound_ ?
                        time - width_left_upper_bound_ : time_inf +
                        width_right_upper_bound_;
                    time_inf = time_inf > 0. ? time_inf : 0.;

                    time_sup = time + width_right_upper_bound_;

                    return;
                }

                time_inf = time_(active_track_index_)(i_time_inf);
                time_inf = time - width_left_upper_bound_ > time_inf +
                    width_right_upper_bound_ ?
                    time - width_left_upper_bound_ : time_inf +
                    width_right_upper_bound_;
                time_inf = time_inf > 0. ? time_inf : 0.;

                time_sup = time_(active_track_index_)(i_time_sup);
                time_sup = time + width_right_upper_bound_ < time_sup -
                    width_left_upper_bound_ ?
                    time + width_right_upper_bound_ : time_sup -
                    width_left_upper_bound_;

                return;
            }

            time_inf = time - width_left_upper_bound_;
            time_inf = time_inf > 0. ? time_inf : 0.;

            time_sup = time + width_right_upper_bound_;

            return;
        }
    }


    //! Computes an aggregated observation vector over a list of times.
    /*! The observations that have a non-zero contribution at time \a time are
      aggregated.
      \param[in] observation_time the times of \a observation.
      \param[in] contribution the contributions associated with observations.
      \param[in] observation the observations to be aggregated.
      \param[in] time the time at which the observations should be aggregated.
      \param[out] aggregated_observation the aggregated observation vector.
    */
    template <class T>
    template <class time_vector, class observation_vector2,
              class observation_vector>
    void ObservationAggregator<T>
    ::Aggregate(const time_vector& observation_time,
                const Vector<double>& contribution,
                const observation_vector2& observation,
                double time,
                observation_vector& aggregated_observation)
    {

        /*** Computes aggregated observations ***/

        // Assumes 'Nobservation' is constant.
        double sum = Norm1(contribution);
        int Nobservation = observation(0).GetSize();
        aggregated_observation.Reallocate(Nobservation);
        aggregated_observation.Fill(T(0.));
        int Nt = observation_time.GetSize();
        for (int h = 0; h < Nt; h++)
            Add(contribution(h), observation(h), aggregated_observation);
        Mlt(T(1. / sum), aggregated_observation);

        if (discard_observation_)
            PushTime(time);
    }


    //! Computes aggregated observation over a list of times and observations.
    /*! The observations that have a non-zero contribution at time \a time are
      aggregated. The variables associated with the new aggregated
      observations vector are stored in \a aggregated_variable.
      \param[in] observation_time the time of the given observations.
      \param[in] contribution the contributions associated with observations.
      \param[in] observation_variable variables associated with the
      observations.
      \param[in] observation the given observation observation.
      \param[in] time a given time.
      \param[out] aggregated_variable the variables associated with the
      aggregated observations.
      \param[out] aggregated_observation the aggregated observation.
    */
    template <class T>
    template <class time_vector, class variable_vector2,
              class observation_vector3,
              class variable_vector, class observation_vector2>
    void ObservationAggregator<T>
    ::Aggregate(const time_vector& observation_time,
                const Vector<double>& contribution,
                const variable_vector2& observation_variable,
                const observation_vector3& observation,
                double time,
                variable_vector& aggregated_variable,
                observation_vector2& aggregated_observation)
    {
        // Assumes 'Nobservation' is constant.
        int Nobservation = observation(0, 0).GetSize();

        /*** Builds variable vector ***/

        observation_variable.Flatten(aggregated_variable);
        RemoveDuplicate(aggregated_variable);
        int Nvariable = aggregated_variable.GetSize();

        map<int, int> variable_index_map;
        for (int v = 0; v < Nvariable; v++)
            variable_index_map[aggregated_variable(v)] = v;

        Vector<double> sum(Nvariable);
        sum.Fill(0.);

        /*** Builds aggregated observations ***/

        aggregated_observation.Reallocate(Nvariable);
        for (int v = 0; v < Nvariable; v++)
            aggregated_observation(v).Reallocate(Nobservation);
        aggregated_observation.Fill(T(0.));

        int Nt, variable_index;
        Nt = observation_time.GetSize();
        for (int h = 0; h < Nt; h++)
            for (int v = 0; v < observation(h).GetSize(); v++)
            {
                variable_index =
                    variable_index_map[observation_variable(h, v)];
                sum(variable_index) += contribution(h);
                Add(contribution(h), observation(h, v),
                    aggregated_observation(variable_index));
            }

        for (int v = 0; v < Nvariable; v++)
        {
            if (sum(v) == 0.)
                throw ErrorProcessing("ObservationAggregator::Aggregate()",
                                      "The sum of contributions for variable "
                                      + to_str(v) + " should not be zero.");
            Mlt(T(1. / sum(v)), aggregated_observation(v));
        }

        if (discard_observation_)
            PushTime(time);
    }


    //! Computes aggregated observation over a list of times and observations.
    /*! The observations that have a non-zero contribution at time \a time are
      aggregated. The variables associated with the new aggregated
      observations vector are stored in \a aggregated_variable. The indexes
      associated with the new aggregated observations vector are stored in
      \a aggregated_index.
      \param[in] observation_time the time of the given observations.
      \param[in] contribution the contributions associated with observations.
      \param[in] observation_variable variables associated with observations.
      \param[in] observation_index corresponding observation locations.
      \param[in] observation the given observation observation.
      \param[in] time a given time.
      \param[out] aggregated_variable variables associated with the
      aggregated observations.
      \param[out] aggregated_index the aggregated locations.
      \param[out] aggregated_observation the aggregated observation.
    */
    template <class T>
    template <class time_vector, class variable_vector2,
              class index_vector3, class observation_vector3,
              class variable_vector, class index_vector2,
              class observation_vector2>
    void ObservationAggregator<T>
    ::Aggregate(const time_vector& observation_time,
                const Vector<double>& contribution,
                const variable_vector2& observation_variable,
                const index_vector3& observation_index,
                const observation_vector3& observation,
                double time,
                variable_vector& aggregated_variable,
                index_vector2& aggregated_index,
                observation_vector2& aggregated_observation)
    {
        throw ErrorUndefined("ObservationAggregator::Aggregate"
                             "(const time_vector& observation_time, "
                             "const variable_vector2& observation_variable, "
                             "const index_vector3& observation_index, "
                             "const observation_vector3& observation, "
                             "double time, "
                             "variable_vector& aggregated_variable, "
                             "index_vector2& aggregated_index, "
                             "observation_vector2& aggregated_observation)");
    }


    //! Creates a new track.
    /*!
      \return The index of the new track.
    */
    template <class T>
    int ObservationAggregator<T>::CreateTrack()
    {
        Vector<double> track;
        track.PushBack(numeric_limits<double>::min());
        time_.PushBack(track);
        return time_.GetSize() - 1;
    }


    //! Sets the active track to a given track.
    /*!
      \param[in] track the given track.
    */
    template <class T>
    void ObservationAggregator<T>::SetTrack(int track)
    {
        if (track < 0 || track >= time_.GetSize())
            throw WrongIndex("ObservationAggregator<T>::SetTrack(int track)",
                             string("The track should be in [0, ") +
                             to_str(time_.GetSize() - 1) +
                             "], but is equal to " + to_str(track) + ".");
        active_track_index_ = track;
    }


    //! Returns the last time of the current track.
    /*!
      \return The last time of the current track.
    */
    template <class T>
    double ObservationAggregator<T>::LastTime() const
    {
        int Ntime = time_.GetSize(active_track_index_);
        return time_(active_track_index_, Ntime - 1);
    }


    //! Returns the last time of a given track.
    /*!
      \param[in] track a given track.
      \return The last time of the given track.
    */
    template <class T>
    double ObservationAggregator<T>::LastTime(int track) const
    {
        int Ntime = time_.GetSize(track);
        return time_(track, Ntime - 1);
    }


    //! Pushes a time to the current track.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    void ObservationAggregator<T>::PushTime(double time)
    {
        PushTime(time, active_track_index_);
    }



    //! Pushes a time to a given track.
    /*!
      \param[in] track index of a given track.
      \param[in] time a given time.
    */
    template <class T>
    void ObservationAggregator<T>::PushTime(double time, int track)
    {
        int Ntime = time_(track).GetSize();
        time_(track).Resize(Ntime + 1);
        for (int i = Ntime - 1; i > -1; i--)
        {
            time_(track)(i + 1) = time_(track)(i);
            if (time_(track)(i + 1) < time)
            {
                time_(track)(i + 1) = time;
                return;
            }
        }

        time_(0) = time;
    }


    //! Computes the contributions of given observations at a given time.
    /*!
      \param[in] time the given time.
      \param[in] observation_time the times associated with the given
      observations.
      \param[out] contribution the contributions computed.
    */
    template <class T>
    template <class time_vector>
    void ObservationAggregator<T>
    ::Contribution(double time, const time_vector& observation_time,
                   Vector<double>& contribution)
    {
        int Ntime = observation_time.GetSize();
        contribution.Reallocate(Ntime);

        if (interpolation_type_ == type_step_)
        {
            contribution.Fill(T(1.));
            return;
        }

        if (interpolation_type_ == type_triangle_
            && width_property_ == width_constant_)
        {
            for (int i = 0; i< Ntime; i++)
                contribution(i) = Contribution(time - observation_time(i));

            Logger::Log<3>(*this, "Contribution: {" +
                           to_str(contribution) + "}\n");
            return;
        }

        if (interpolation_type_ == type_interpolation_)
        {
            contribution.Fill(T(0.));

            if (Ntime == 2)
            {
                width_left_ = observation_time(1) - observation_time(0);
                width_right_ = observation_time(1) - observation_time(0);
                contribution(0) = Contribution(time - observation_time(0));
                contribution(1) = Contribution(observation_time(1) - time);
            }

            Logger::Log<3>(*this, "Contribution: {" +
                           to_str(contribution) + "}\n");
            return;
        }

        throw ErrorArgument("ObservationAggregator"
                            "::Contribution(double time, "
                            "const time_vector& observation_time,"
                            "Vector<double>& contribution) const");
    }


    //! Computes the contributions of given observations at a given time.
    /*!
      \param[in] time the given time.
      \param[in] observation_time the times associated with the given
      observations.
      \param[in] width_left the non constant triangle width left.
      \param[in] width_right the non constant triangle width right.
      \param[out] contribution the contributions computed.
    */
    template <class T>
    template <class time_vector>
    void ObservationAggregator<T>
    ::Contribution(double time, const time_vector& observation_time,
                   Vector<double>& width_left, Vector<double>& width_right,
                   Vector<double>& contribution)
    {
        int Ntime = observation_time.GetSize();
        contribution.Reallocate(Ntime);

        if (interpolation_type_ == type_triangle_
            && width_property_ == width_per_observation_)
        {
            for (int i = 0; i< Ntime; i++)
            {
                width_left_ = width_left(i);
                width_right_ = width_right(i);
                contribution(i) = Contribution(time - observation_time(i));
            }

            return;
        }

        throw ErrorArgument("ObservationAggregator"
                            "::Contribution(double time, "
                            "const time_vector& observation_time,"
                            "Vector<double>& contribution) const");
    }


    //! Computes the contribution of an observation given a time difference.
    /*!
      \param[in] delta_t the time difference.
    */
    template <class T>
    double ObservationAggregator<T>::Contribution(double delta_t) const
    {
        if (delta_t > 0.)
            return 1. - delta_t / width_left_;
        else
            return 1. - delta_t / width_right_;
    }


    //! Returns the least interval containing a given value.
    /*! This method returns the index of the closest value among the vector
      \a X elements that is lower than the given value \a value and the
      index of the closest value that is higher.
      \param[in] X a sorted vector.
      \param[in] value a given value.
      \param[out] index_inf index of the closest value lower than \a value.
      \param[out] index_sup index of the closest value higher than \a value.
    */
    template <class T>
    template <class time_vector>
    void ObservationAggregator<T>::GetValueIndex(time_vector& X, double value,
                                                 int& index_inf,
                                                 int& index_sup) const
    {
        int Nx = X.GetSize();
        index_inf = -1;
        index_sup = -1;
        for (int i = Nx - 1; i > -1; i--)
            if (X(i) < value)
            {
                index_inf = i;
                if (i < Nx - 1)
                    index_sup = i + 1;
                return;
            }
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string ObservationAggregator<T>::GetName() const
    {
        return "ObservationAggregator";
    }


}


#define VERDANDI_FILE_OBSERVATION_MANAGER_OSERVATIONAGGREGATOR_CXX
#endif
