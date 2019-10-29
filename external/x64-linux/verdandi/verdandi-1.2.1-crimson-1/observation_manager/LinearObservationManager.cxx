// Copyright (C) 2008-2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton, Marc Fragu
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_LINEAROBSERVATIONMANAGER_CXX

#include <cstdlib>
#include <limits>
#include "LinearObservationManager.hxx"
#include QUOTE(OBSERVATION_AGGREGATOR.cxx)


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    /*! It entirely defines the operator: no dimension or size is associated
      with this implementation.
    */
    template <class T>
    LinearObservationManager<T>::LinearObservationManager()
    {
    }


    //! Main constructor.
    /*!
      \param[in] model model.
      \param[in] configuration_file configuration_file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    LinearObservationManager<T>
    ::LinearObservationManager(const Model& model,
                               string configuration_file):
        observation_aggregator_(configuration_file)
    {
    }


    //! Destructor.
    template <class T>
    LinearObservationManager<T>
    ::~LinearObservationManager()
    {
    }


    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //! Initializes the observation manager.
    /*!
      \param[in] model model.
      \param[in] configuration_file configuration file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    void LinearObservationManager<T>
    ::Initialize(const Model& model, string configuration_file)
    {
        observation_aggregator_.Initialize(configuration_file);

        VerdandiOps configuration(configuration_file);

        Nstate_model_ = model.GetNstate();

        configuration.SetPrefix("observation.");
        configuration.Set("file", observation_file_);
        configuration.Set("type", "", "state", observation_type_);
        configuration.Set("Delta_t_constant", is_delta_t_constant_);
        if (!is_delta_t_constant_)
        {
            configuration.Set("observation_time_file",
                              observation_time_file_);
            observation_time_.Read(observation_time_file_);
        }
        else
            configuration.Set("Delta_t", "v > 0", Delta_t_);

        configuration.Set("Nskip", "v > 0", Nskip_);
        configuration.Set("initial_time", "", 0., initial_time_);
        configuration.Set("final_time", "", numeric_limits<double>::max(),
                          final_time_);

        time_ = numeric_limits<double>::min();

        configuration.Set("width_file", "", "", width_file_);
        configuration.Set("error.variance", "v > 0", error_variance_value_);

        /*** Building the matrices ***/

        configuration.Set("operator.scaled_identity",
                          operator_scaled_identity_);

        if (operator_scaled_identity_)
        {
            configuration.Set("operator.diagonal_value",
                              operator_diagonal_value_);

            Nobservation_ = Nstate_model_;
#ifdef VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
            build_diagonal_sparse_matrix(Nstate_model_,
                                         operator_diagonal_value_,
                                         tangent_operator_matrix_);
#else
            tangent_operator_matrix_.Reallocate(Nobservation_, Nstate_model_);
            tangent_operator_matrix_.SetIdentity();
            Mlt(operator_diagonal_value_, tangent_operator_matrix_);
#endif
        }
        else // Operator given in a Lua table or in a file.
        {
#ifdef VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
            if (configuration.Is<string>("operator.value"))
            {
                Matrix<T, General, ArrayRowSparse> tmp;
                tmp.Read(configuration.Get<string>("operator.value"));
                Copy(tmp, tangent_operator_matrix_);
            }
            else
            {
                Matrix<T> tmp;
                configuration.Set("operator.value", tmp);
                if (tmp.GetN() % model.GetNstate() != 0)
                    throw ErrorArgument("LinearObservationManager"
                                        "::Initialize()",
                                        "The total number of elements in the "
                                        "tangent operator matrix ("
                                        + to_str(tangent_operator_matrix_
                                                 .GetN())
                                        + ") is not a multiple of the"
                                        " dimension of the model state ("
                                        + to_str(model.GetNstate()) + ").");
                tmp.Resize(tmp.GetN() / model.GetNstate(), model.GetNstate());
                Matrix<T, General, ArrayRowSparse> tmp_array;
                ConvertDenseToArrayRowSparse(tmp, tmp_array);
                Copy(tmp_array, tangent_operator_matrix_);
            }
            if (tangent_operator_matrix_.GetN() != model.GetNstate())
                throw ErrorArgument("LinearObservationManager::Initialize()",
                                    "The number of columns of the tangent "
                                    "operator matrix ("
                                    + to_str(tangent_operator_matrix_.GetN())
                                    + ") is inconsistent with the"
                                    " dimension of the model state ("
                                    + to_str(model.GetNstate()) + ").");
            Nobservation_ = tangent_operator_matrix_.GetM();
#else
            if (configuration.Is<string>("operator.value"))
                tangent_operator_matrix_
                    .Read(configuration.Get<string>("operator.value"));
            else
            {
                tangent_operator_matrix_.Clear();
                configuration.Set("operator.value", tangent_operator_matrix_);
                if (tangent_operator_matrix_.GetN() % model.GetNstate() != 0)
                    throw ErrorArgument("LinearObservationManager"
                                        "::Initialize()",
                                        "The total number of elements in the "
                                        "tangent operator matrix ("
                                        + to_str(tangent_operator_matrix_
                                                 .GetN())
                                        + ") is not a multiple of the"
                                        " dimension of the model state ("
                                        + to_str(model.GetNstate()) + ").");
                tangent_operator_matrix_
                    .Resize(tangent_operator_matrix_.GetN()
                            / model.GetNstate(),
                            model.GetNstate());
            }
            if (tangent_operator_matrix_.GetN() != model.GetNstate())
                throw ErrorArgument("LinearObservationManager::Initialize()",
                                    "The number of columns of the tangent "
                                    "operator matrix ("
                                    + to_str(tangent_operator_matrix_.GetN())
                                    + ") is inconsistent with the"
                                    " dimension of the model state ("
                                    + to_str(model.GetNstate()) + ").");
            Nobservation_ = tangent_operator_matrix_.GetM();
#endif
        }

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        build_diagonal_sparse_matrix(Nobservation_, error_variance_value_,
                                     error_variance_);
        build_diagonal_sparse_matrix(Nobservation_,
                                     T(T(1) / error_variance_value_),
                                     error_variance_inverse_);
#else
        error_variance_.Reallocate(Nobservation_, Nobservation_);
        error_variance_.SetIdentity();
        Mlt(error_variance_value_, error_variance_);
        error_variance_inverse_.Reallocate(Nobservation_, Nobservation_);
        error_variance_inverse_.SetIdentity();
        Mlt(T(T(1)/ error_variance_value_), error_variance_inverse_);
#endif

        if (observation_type_ == "state")
            Nbyte_observation_ = Nstate_model_ * sizeof(T) + sizeof(int);

        if (observation_type_ == "observation")
            Nbyte_observation_ = Nobservation_ * sizeof(T) + sizeof(int);

        if (is_delta_t_constant_)
        {

            int expected_file_size;
            expected_file_size = Nbyte_observation_
                * int((final_time_ - initial_time_)
                      / (Delta_t_ * double(Nskip_)) + 1.);

            int file_size;
            ifstream file_stream;
            file_stream.open(observation_file_.c_str());

#ifdef VERDANDI_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw ErrorIO("LinearObservationManager"
                              "::Initialize(model, configuration_file)",
                              "Unable to open file \""
                              + observation_file_ + "\".");
#endif

            file_stream.seekg(0, ios_base::end);
            file_size = file_stream.tellg() ;
            file_stream.close();

            if (expected_file_size > file_size)
                throw IOError("LinearObservationManager"
                              "::Initialize(model, configuration_file)",
                              "Too few available observations, the size of \""
                              + observation_file_ + "\" must be greater than "
                              + to_str(expected_file_size) + "B.");
        }
        else
        {
            int file_size;
            ifstream file_stream;
            file_stream.open(observation_file_.c_str());

#ifdef VERDANDI_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw ErrorIO("LinearObservationManager"
                              "::Initialize(model, configuration_file)",
                              "Unable to open file \""
                              + observation_file_ + "\".");
#endif
            file_stream.close();
        }
    }


    /////////////////////////////
    // OBSERVATIONS MANAGEMENT //
    /////////////////////////////


    //! Activates or deactivates the option 'discard_observation'.
    /*!
      \param[in] discard_observation if set to true, each observation will be
      used at most one time.
    */
    template <class T>
    void LinearObservationManager<T>
    ::DiscardObservation(bool discard_observation)
    {
        observation_aggregator_.DiscardObservation(discard_observation);
    }


    //! Creates a new track.
    /*!
      \return The index of the new track.
    */
    template <class T>
    int LinearObservationManager<T>
    ::CreateTrack()
    {
        return observation_aggregator_.CreateTrack();
    }


    //! Sets the track to a given track.
    /*!
      \param[in] track the given track.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetTrack(int track)
    {
        observation_aggregator_.SetTrack(track);
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] model the model.
      \param[in] time a given time.
    */
    template <class T>
    template <class Model>
    void LinearObservationManager<T>
    ::SetTime(const Model& model, double time)
    {
        SetTime(time);
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetTime(double time)
    {
        if (time_ == time)
            return;

        time_ = time;
        SetAvailableTime(time_, available_time_);
    }


    //! Sets the available observation times at a given time.
    /*!
      \param[in] time the given time.
      \param[out] available_time the available observation times.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetAvailableTime(double time, time_vector&
                       available_time)
    {
        double time_inf, time_sup;
        int selection_policy;
        observation_aggregator_.GetContributionInterval(time, time_inf,
                                                        time_sup,
                                                        selection_policy);
        SetAvailableTime(time, time_inf, time_sup, selection_policy,
                         available_time);

        Logger::Log<3>(*this, to_str(time) + ", [" + to_str(time_inf) + " " +
                       to_str(time_sup) + "], {" + to_str(available_time) +
                       "}\n");
    }


    //! Sets available observation times at a given time interval.
    /*!
      \param[in] time_inf lower bound of the given time interval.
      \param[in] time_sup upper bound (excluded) of the given time interval.
      \param[in] selection_policy interval selection policy.
      \param[out] available_time the available observation times.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetAvailableTime(double time_inf, double time_sup,
                       time_vector& available_time)
    {
        available_time.Clear();

        if (is_delta_t_constant_)
        {
            double period = Delta_t_ * Nskip_;
            double available_time_0
                = initial_time_
                + floor((time_inf - initial_time_) / period) * period;
            if (available_time_0 == time_inf)
                available_time.PushBack(available_time_0);
            available_time_0 += period;
            for (double t = available_time_0; t < time_sup; t += period)
                available_time.PushBack(t);
        }
        else
            for (int i = 0; i < observation_time_.GetM(); i++)
                if (observation_time_(i) >= time_inf &&
                    observation_time_(i) <= time_sup)
                    available_time.PushBack(observation_time_(i));
        return;
    }


    //! Sets available observation times at a given time interval.
    /*!
      \param[in] time_inf lower bound of the given time interval.
      \param[in] time_sup upper bound (excluded) of the given time interval.
      \param[in] selection_policy interval selection policy.
      \param[out] available_time the available observation times.
    */
    template <class T>
    void LinearObservationManager<T>
    ::SetAvailableTime(double time, double time_inf, double time_sup,
                       int selection_policy,
                       time_vector& available_time)
    {
        available_time.Clear();
        time_inf = time_inf > initial_time_ ? time_inf : initial_time_;
        time_sup = time_sup < final_time_ ? time_sup : final_time_;


        if (is_delta_t_constant_)
        {
            double period = Delta_t_ * Nskip_;

            // All observations available in the given interval are
            // considered.
            if (selection_policy == 0)
            {
                double available_time_0
                    = initial_time_
                    + floor((time_inf - initial_time_) / period) * period;
                if (available_time_0 == time_inf)
                    available_time.PushBack(available_time_0);
                available_time_0 += period;
                for (double t = available_time_0; t < time_sup; t += period)
                    available_time.PushBack(t);
                observation_aggregator_.Contribution(time_, available_time_,
                                                     contribution_);
                return;
            }

            // Only the closest left observation and the closest right
            // observation are requested.
            if (selection_policy == 2)
            {
                double t1, t2;

                if (is_multiple(time - initial_time_, period))
                    t1 = time;
                else
                    t1 = initial_time_
                        + floor((time - initial_time_) / period) * period;

                t2 = t1 + period;

                if (t1 <= final_time_)
                    available_time.PushBack(t1);
                if (t2 <= final_time_)
                    available_time.PushBack(t2);
                observation_aggregator_.Contribution(time_, available_time_,
                                                     contribution_);
                return;
            }

            // All observations available in the given interval are considered
            // taking into account non constant triangle widths associated
            // with observations.
            if (selection_policy == 3)
            {
                double available_time_0, available_time_1;
                int Nobservation;

                Vector<double> width_left, width_right, available_width_left,
                    available_width_right;
                ReadObservationTriangleWidth(time_inf, time_sup, width_left,
                                             width_right);

                available_time_0 = initial_time_
                    + floor((time_inf - initial_time_) / period) * period;
                if (!is_equal(available_time_0, time_inf))
                    available_time_0 += period;

                available_time_1 = initial_time_
                    + floor((time_sup - initial_time_) / period) * period;
                if (is_equal(available_time_1, time_sup))
                    available_time_1 -= period;

                Nobservation = floor((available_time_1 - available_time_0)
                                     / period) + 1;

                double t = available_time_0;
                for (int i = 0; i < Nobservation; i++, t += period)
                {
                    if (t < time)
                    {
                        if (t + width_right(i) > time)
                        {
                            available_time.PushBack(t);
                            available_width_left.PushBack(width_left(i));
                            available_width_right.PushBack(width_right(i));
                        }
                    }
                    else if (t > time)
                    {
                        if (t - width_left(i) < time)
                        {
                            available_time.PushBack(t);
                            available_width_left.PushBack(width_left(i));
                            available_width_right.PushBack(width_right(i));
                        }
                    }
                    else
                    {
                        available_time.PushBack(t);
                        available_width_left.PushBack(width_left(i));
                        available_width_right.PushBack(width_right(i));
                    }

                }
                observation_aggregator_.
                    Contribution(time_, available_time_,
                                 available_width_left, available_width_right,
                                 contribution_);
                return;
            }
        }
        else
        {
            // All observations available in the given interval are
            // considered.
            if (selection_policy == 0)
            {
                for (int i = 0; i < observation_time_.GetM(); i++)
                    if (observation_time_(i) >= time_inf &&
                        observation_time_(i) <= time_sup)
                        available_time.PushBack(observation_time_(i));
                observation_aggregator_.Contribution(time_, available_time_,
                                                     contribution_);
                return;
            }
            // Only the closest left observation and the closest right
            // observation are requested.
            if (selection_policy == 2)
            {
                if (observation_time_(0) > time)
                    return;
                for (int i = 1; i < observation_time_.GetM(); i++)
                    if (observation_time_(i) >= time)
                    {
                        available_time.PushBack(observation_time_(i - 1));
                        available_time.PushBack(observation_time_(i));
                        observation_aggregator_.
                            Contribution(time_, available_time_,
                                         contribution_);
                        return;
                    }
                return;
            }
        }

        throw ErrorArgument("void LinearObservationManager<T>"
                            "::SetAvailableTime(double time,"
                            " double time_inf, double time_sup,"
                            " int selection_policy,"
                            "LinearObservationManager<T>"
                            "::time_vector&"
                            "available_time) const");
    }


    ////////////////////////////
    // FLATTENED OBSERVATIONS //
    ////////////////////////////


    /*** Gets observations ***/


    //! Gets observations flattened over a list of times.
    /*! The observations available at time \a time are loaded and concatenated
      in a vector.
      \param[in] time the given time.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double time,
                              observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetFlattenedObservation(available_time, observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations in the interval [\a time_inf, \a time_sup[
      are loaded and concatenated in a vector.
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double time_inf, double time_sup,
                              observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetFlattenedObservation(available_time, observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations available are loaded and concatenated in a vector.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(observation_vector& observation)
    {
        GetFlattenedObservation(available_time_, observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations available at the times \a available_time are loaded
      and concatenated in a vector.
      \param[in] available_time the given observation time vector.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(const time_vector& available_time,
                              observation_vector& observation)
    {
        observation_vector2 observation2;
        GetRawObservation(available_time, observation2);
        observation2.Flatten(observation);
    }


    /*** Gets observations and associated variables ***/


    //! Gets observations flattened over a list of times.
    /*! The observations available at time \a time are loaded and concatenated
      in a vector.
      \param[in] time the given time.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double time,
                              variable_vector& observation_variable,
                              observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetFlattenedObservation(available_time, observation_variable,
                                observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations in the interval [\a time_inf, \a time_sup[
      are loaded and concatenated in a vector.
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double time_inf, double time_sup,
                              variable_vector& observation_variable,
                              observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetFlattenedObservation(available_time, observation_variable,
                                observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations available are loaded and concatenated in a vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(variable_vector& observation_variable,
                              observation_vector& observation)
    {
        GetFlattenedObservation(available_time_, observation_variable,
                                observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations available at the times \a available_time are loaded
      and concatenated in a vector.
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(const time_vector& available_time,
                              variable_vector& observation_variable,
                              observation_vector& observation)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        GetRawObservation(available_time, observation_variable2,
                          observation3);
        observation3.Flatten(observation);
        observation_variable2.Flatten(observation_variable);
    }


    /*** Gets observations, associated variables and associated index ***/


    //! Gets observations flattened over a list of times.
    /*! The observations available at time \a time are loaded and concatenated
      in a vector.
      \param[in] time the given time.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double time,
                              variable_vector& observation_variable,
                              index_vector& observation_index,
                              observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetFlattenedObservation(available_time, observation_variable,
                                observation_index, observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations in the interval [\a time_inf, \a time_sup[
      are loaded and concatenated in a vector.
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be flattened.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(double time_inf, double time_sup,
                              variable_vector& observation_variable,
                              index_vector& observation_index,
                              observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetFlattenedObservation(available_time, observation_variable,
                                observation_index, observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations available are loaded and concatenated in a vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(variable_vector& observation_variable,
                              index_vector& observation_index,
                              observation_vector& observation)
    {
        GetFlattenedObservation(available_time_, observation_variable,
                                observation_index, observation);
    }


    //! Gets observations flattened over a list of times.
    /*! The observations available at the times \a available_time are loaded
      and concatenated in a vector.
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index indexes associated with the observations.
      \param[out] observation the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetFlattenedObservation(
        const time_vector& available_time,
        variable_vector& observation_variable,
        index_vector& observation_index,
        observation_vector& observation)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        index_vector3 observation_index3;
        GetRawObservation(available_time, observation_variable2,
                          observation_index3, observation3);
        observation3.Flatten(observation);
        observation_variable2.Flatten(observation_variable);
        observation_index3.Flatten(observation_index);
    }


    /////////////////////////////
    // AGGREGATED OBSERVATIONS //
    /////////////////////////////


    /*** Gets observations ***/


    //! Gets observations aggregated over a list of times.
    /*! The observations available at the given time are loaded and
      aggregated.
      \param[in] time the given time.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double time,
        observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetAggregatedObservation(available_time, observation);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations in the interval [\a time_inf, \a time_sup[ are loaded
      and aggregated.
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double time_inf, double time_sup,
        observation_vector& observation)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetAggregatedObservation(available_time, observation);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations available are loaded and
      aggregated.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(observation_vector& observation)
    {
        GetAggregatedObservation(available_time_, observation);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations available at the times \a available_time are loaded
      and aggregated.
      \param[in] available_time the given observation time vector.
      \param[out] observation the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(const time_vector& available_time,
                               observation_vector& observation)
    {
        observation_vector2 observation2;
        GetRawObservation(available_time, observation2);
        observation_aggregator_.Aggregate(available_time, contribution_,
                                          observation2, time_, observation);
    }


    /*** Gets observations and associated variables ***/


    //! Gets observations aggregated over a list of times.
    /*! The observations available at the given time are loaded and
      aggregated.
      \param[in] time the given time.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double time,
        variable_vector& observation_variable,
        observation_vector2& observation2)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetAggregatedObservation(available_time, observation_variable,
                                 observation2);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations in the interval [\a time_inf, \a time_sup[ are loaded
      and aggregated.
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double time_inf, double time_sup,
        variable_vector& observation_variable,
        observation_vector2& observation2)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetAggregatedObservation(available_time, observation_variable,
                                 observation2);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations available are loaded and
      aggregated.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(variable_vector& observation_variable,
                               observation_vector2& observation2)
    {
        GetAggregatedObservation(available_time_, observation_variable,
                                 observation2);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations available at the times \a available_time are loaded
      and aggregated.
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(const time_vector& available_time,
                               variable_vector& observation_variable,
                               observation_vector2& observation2)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        GetRawObservation(available_time, observation_variable2,
                          observation3);
        observation_aggregator_.Aggregate(available_time,
                                          contribution_,
                                          observation_variable2,
                                          observation3,
                                          time_,
                                          observation_variable,
                                          observation2);
    }


    /*** Gets observations, associated variables and associated index ***/


    //! Gets observations aggregated over a list of times.
    /*! The observations available at the given time are loaded and
      aggregated.
      \param[in] time the given time.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index2 indexes associated with the observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double time,
        variable_vector& observation_variable,
        index_vector2& observation_index2,
        observation_vector2& observation2)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetAggregatedObservation(available_time, observation_variable,
                                 observation_index2, observation2);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations in the interval [\a time_inf, \a time_sup[ are loaded
      and aggregated.
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index2 indexes associated with the observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        double time_inf, double time_sup,
        variable_vector& observation_variable,
        index_vector2& observation_index2,
        observation_vector2& observation2)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetAggregatedObservation(available_time, observation_variable,
                                 observation_index2, observation2);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations available are loaded and
      aggregated.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index2 indexes associated with the observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        variable_vector& observation_variable,
        index_vector2& observation_index2,
        observation_vector2& observation2)
    {
        GetAggregatedObservation(available_time_, observation_variable,
                                 observation_index2, observation2);
    }


    //! Gets observations aggregated over a list of times.
    /*! The observations available at the times \a available_time are loaded
      and aggregated.
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable variables associated with the
      observations.
      \param[out] observation_index2 indexes associated with the observations.
      \param[out] observation2 the aggregated observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetAggregatedObservation(
        const time_vector& available_time,
        variable_vector& observation_variable,
        index_vector2& observation_index2,
        observation_vector2& observation2)
    {
        observation_vector3 observation3;
        variable_vector2 observation_variable2;
        index_vector3 observation_index3;
        GetRawObservation(available_time, observation_variable2,
                          observation_index3, observation3);
        observation_aggregator_.Aggregate(available_time,
                                          contribution_,
                                          observation_variable2,
                                          observation_index3,
                                          observation3,
                                          time_,
                                          observation_variable,
                                          observation_index2,
                                          observation2);
    }


    //////////////////
    // OBSERVATIONS //
    //////////////////


    /*** Gets observations ***/


    //! Gets available observations at a given time.
    /*!
      \param[in] time the given time.
      \param[out] observation2 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double time,
        observation_vector2& observation2)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetRawObservation(available_time, observation2);
    }


    //! Gets observations available in a given interval.
    /*!
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation2 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double time_inf, double time_sup,
        observation_vector2& observation2)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetRawObservation(available_time, observation2);
    }


    //! Gets observations at the current time.
    /*!
      \param[out] observation2 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(observation_vector2& observation2)
    {
        GetRawObservation(available_time_, observation2);
    }


    //! Gets observations of a list of times.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation2 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        const time_vector& available_time,
        observation_vector2& observation2)
    {
        ReadObservation(available_time, observation2);
    }


    /*** Gets observations and associated variables ***/


    //! Gets available observations at a given time.
    /*!
      \param[in] time the given time.
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double time,
        variable_vector2& observation_variable2,
        observation_vector3& observation3)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetRawObservation(available_time, observation_variable2,
                          observation3);
    }


    //! Gets observations available in a given interval.
    /*!
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double time_inf, double time_sup,
        variable_vector2& observation_variable2,
        observation_vector3& observation3)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetRawObservation(available_time, observation_variable2,
                          observation3);
    }



    //! Gets observations at the current time.
    /*!
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        variable_vector2& observation_variable2,
        observation_vector3& observation3)
    {
        GetRawObservation(available_time_, observation_variable2,
                          observation3);
    }


    //! Gets observations of a list of times.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        const time_vector& available_time,
        variable_vector2& observation_variable2,
        observation_vector3& observation3)
    {
        ReadObservationVariable(available_time, observation_variable2);
        ReadObservation(available_time, observation_variable2, observation3);
    }


    /*** Gets observations, associated variables and associated index ***/


    //! Gets available observations at a given time.
    /*!
      \param[in] time the given time.
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation_index3 indexes associated with the observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double time,
        variable_vector2& observation_variable2,
        index_vector3& observation_index3,
        observation_vector3& observation3)
    {
        time_vector available_time;
        SetAvailableTime(time, available_time);
        GetRawObservation(available_time, observation_variable2,
                          observation_index3, observation3);
    }


    //! Gets observations available in a given interval.
    /*!
      \param[in] time_inf lower bound of the given interval.
      \param[in] time_sup upper bound (excluded) of the given interval.
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation_index3 indexes associated with the observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        double time_inf, double time_sup,
        variable_vector2& observation_variable2,
        index_vector3& observation_index3,
        observation_vector3& observation3)
    {
        time_vector available_time;
        SetAvailableTime(time_inf, time_sup, available_time);
        GetRawObservation(available_time, observation_variable2,
                          observation_index3, observation3);
    }



    //! Gets observations at the current time.
    /*!
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation_index3 indexes associated with the observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        variable_vector2& observation_variable2,
        index_vector3& observation_index3,
        observation_vector3& observation3)
    {
        GetRawObservation(available_time_, observation_variable2,
                          observation_index3, observation3);
    }


    //! Gets observations of a list of times.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable2 variables associated with the
      observations.
      \param[out] observation_index3 indexes associated with the observations.
      \param[out] observation3 the observation to be loaded.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetRawObservation(
        const time_vector& available_time,
        variable_vector2& observation_variable2,
        index_vector3& observation_index3,
        observation_vector3& observation3)
    {
        ReadObservationVariable(available_time, observation_variable2);
        ReadObservation(available_time, observation_variable2, observation3);
        ReadObservationIndex(available_time, observation_variable2,
                             observation_index3);
    }


    ///////////////////////////
    // READ OBSERVATION DATA //
    ///////////////////////////


    //! Builds variables vector associated with given observations.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation_variable2 variables associated with the
      observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservationVariable(
        const time_vector& available_time,
        variable_vector2& observation_variable2)
        const
    {
        int Nt = available_time.GetSize();
        observation_variable2.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
            observation_variable2(h).PushBack(0);
    }


    //! Builds observations associated with given times and variables.
    /*!
      \param[in] available_time the given observation time vector.
      \param[in] observation_variable2 variables associated with the
      observations.
      \param[out] observation3 the observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservation(
        const time_vector& available_time,
        const variable_vector2&
        observation_variable2,
        observation_vector3& observation3) const
    {
        ifstream file_stream;
        file_stream.open(observation_file_.c_str());
#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("LinearObservationManager"
                          "::LoadObservation(model)",
                          "Unable to open file \""
                          + observation_file_ + "\".");
#endif

        int Nvariable, Nt;
        Nt = available_time.GetSize();
        observation3.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
        {
            Nvariable = observation_variable2(h).GetSize();
            observation3.Reallocate(h, Nvariable);
            for (int v = 0; v < Nvariable; v++)
                ReadObservation(file_stream,
                                available_time(h),
                                observation_variable2(h, v),
                                observation3.GetVector(h, v));
        }

        file_stream.close();
    }


    //! Builds observations associated with given times.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation2 the observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservation(
        const time_vector& available_time,
        observation_vector2& observation2) const
    {
        ifstream file_stream;
        file_stream.open(observation_file_.c_str());
#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("LinearObservationManager"
                          "::LoadObservation(model)",
                          "Unable to open file \""
                          + observation_file_ + "\".");
#endif

        int Nt = available_time.GetSize();
        observation2.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
            ReadObservation(file_stream, available_time(h), 0,
                            observation2(h));

        file_stream.close();
    }


    //! Reads observation from observation file given a time and a variable.
    /*!
      \param[in] file_stream the observation file stream.
      \param[in] time the time.
      \param[in] variable the variable.
      \param[out] observation the observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservation(ifstream& file_stream, double time, int variable,
                      observation_vector&
                      observation) const
    {
        observation.Reallocate(Nobservation_);
        observation_vector input_data;

        streampos position;

        if (is_delta_t_constant_)
            position = (floor((time - initial_time_)
                              / (Delta_t_ * Nskip_) + 0.5) + variable)
                * Nbyte_observation_;
        else
        {
            int index;
            bool observation_available = false;
            for (index = 0; index < observation_time_.GetM(); index++)
                if (observation_time_(index) == time)
                {
                    observation_available = true;
                    break;
                }
            if (!observation_available)
                throw ErrorIO("LinearObservationManager::ReadObservation"
                              "(ifstream& file_stream, double time, "
                              "int variable, LinearObservationManager"
                              "::observation_vector& observation) const",
                              "No observation available at time "
                              + to_str(time) + ".");
            position = index * Nbyte_observation_;
        }
        file_stream.seekg(position);
        input_data.Read(file_stream);

        if (observation_type_ == "state")
        {
            if (input_data.GetSize() != Nstate_model_)
                throw ErrorIO("LinearObservationManager::ReadObservation"
                              "(ifstream& file_stream, double time, "
                              "int variable, LinearObservationManager"
                              "::observation_vector& observation) const",
                              "The observation type is 'state', so the whole"
                              " model state is supposed to be stored, but "
                              "the size of the observation read at time "
                              + to_str(time_) + " (Nread = "
                              + to_str(input_data.GetSize()) +
                              ") mismatches with the expected size (Nstate = "
                              + to_str(Nstate_model_) + ").");
            ApplyOperator(input_data, observation);
        }
        else
        {
            if (input_data.GetSize() != Nobservation_)
                throw ErrorIO("LinearObservationManager::ReadObservation"
                              "(ifstream& file_stream, double time, "
                              "int variable, LinearObservationManager"
                              "::observation_vector& observation) const",
                              "The observation type is 'observation', so "
                              "only observations are stored in the file, but"
                              " the size of the observation read at time "
                              + to_str(time_) + " (Nread = "
                              + to_str(input_data.GetSize())
                              + ") mismatches with the "
                              "expected size (Nobservation = "
                              + to_str(Nobservation_ ) + ").");
            Copy(input_data, observation);
        }
    }


    //! Reads observations indexes.
    /*!
      \param[in] available_time the available time.
      \param[in] observation_variable2 variable associated with the
      observations.
      \param[out] observation_index3 the indexes associated with the
      observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservationIndex(
        const time_vector& available_time,
        const variable_vector2&
        observation_variable2,
        index_vector3& observation_index3) const
    {
        int Nvariable, Nt;
        Nt = available_time.GetSize();
        observation_index3.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
        {
            Nvariable = observation_variable2.GetSize(h);
            observation_index3.Reallocate(h, Nvariable);
            for (int v = 0; v < Nvariable; v++)
            {
                observation_index3.Reallocate(h, v, Nobservation_);
                observation_index3(h, v).Fill(0);
            }
        }
    }


    //! Reads triangle width associated with observations of a given interval.
    /*!
      \param[in] time_inf lower bound of a given interval.
      \param[in] time_sup upper bound (excluded) of a given interval.
      \param[out] width_left left widths associated with observations.
      \param[out] width_right right widths associated with observations.
    */
    template <class T>
    void LinearObservationManager<T>
    ::ReadObservationTriangleWidth(double time_inf, double time_sup,
                                   Vector<double>& width_left,
                                   Vector<double>& width_right) const
    {
        double period, available_time_0, available_time_1;
        int Nwidth, Nbyte_width;
        Vector<double> input_data;

        ifstream file_stream;
        file_stream.open(width_file_.c_str());
        streampos position;
#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("LinearObservationManager"
                          "::ReadObservationTriangleWidth()",
                          "Unable to open file \""
                          + width_file_ + "\".");
#endif
        period = Delta_t_ * Nskip_;
        available_time_0 = initial_time_
            + floor((time_inf - initial_time_) / period) * period;
        if (!is_equal(available_time_0, time_inf))
            available_time_0 += period;
        available_time_1 = initial_time_
            + floor((time_sup - initial_time_) / period) * period;
        if (is_equal(available_time_1, time_sup))
            available_time_1 -= period;

        Nwidth = floor((available_time_1 - available_time_0) / period) + 1;
        Nbyte_width = 2 * sizeof(double) + sizeof(int);
        width_left.Reallocate(Nwidth);
        width_right.Reallocate(Nwidth);

        position = floor((available_time_0 - initial_time_) / Delta_t_)
            * Nbyte_width;
        for (int i = 0; i < Nwidth; i++, position += Nskip_ * Nbyte_width)
        {
            file_stream.seekg(position);
            input_data.Read(file_stream);
            width_left(i) = input_data(0);
            width_right(i) = input_data(1);
            input_data.Clear();
        }

        file_stream.close();
    }


    /////////////////
    // OBSERVATION //
    /////////////////


    //! Gets observation.
    /*!
      \param[out] observation observation vector.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetObservation(observation& observation)
    {
        GetAggregatedObservation(observation);
    }


    ////////////////
    // INNOVATION //
    ////////////////


    //! Gets innovation.
    /*!
      \param[in] state state vector.
      \param[out] innovation innovation vector.
    */
    template <class T>
    template <class state>
    void LinearObservationManager<T>
    ::GetInnovation(const state& x,
                    observation& innovation)
    {
        innovation.Reallocate(Nobservation_);
        observation observation;
        GetObservation(observation);
        ApplyOperator(x, innovation);
        Mlt(T(-1), innovation);
        Add(T(1), observation, innovation);
    }


    ////////////
    // ACCESS //
    ////////////


    //! Indicates if some observations are available at a given time.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    bool LinearObservationManager<T>::HasObservation(double time)
    {
        SetTime(time);
        return available_time_.GetSize() != 0
            && !is_equal(Norm1(contribution_), 0.);
    }


    //! Indicates if some observations are available at current time.
    template <class T>
    bool LinearObservationManager<T>::HasObservation() const
    {
        return available_time_.GetSize() != 0
            && !is_equal(Norm1(contribution_), 0.);
    }


    //! Gets Nobservation_ value.
    /*!
      \return The total number of observation at current time.
    */
    template <class T>
    int LinearObservationManager<T>::GetNobservation() const
    {
        return Nobservation_;
    }



    /*! \brief Checks whether the observation operator is available in a
      sparse matrix. */
    /*!
      \return True if the observation operator is available in a sparse
      matrix, false otherwise.
    */
    template <class T>
    bool LinearObservationManager<T>::IsOperatorSparse() const
    {
#ifdef VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
        return true;
#else
        return false;
#endif
    }


    //! Checks whether the observation error covariance matrix is sparse.
    /*!
      \return True if the observation error covariance matrix is sparse, false
      otherwise.
    */
    template <class T>
    bool LinearObservationManager<T>::IsErrorSparse() const
    {
#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        return true;
#else
        return false;
#endif
    }


    /*! \brief Checks whether the observation error covariance is available in
      a matrix. */
    /*!
      \return True if the observation error covariance is available in a
      matrix, false otherwise.
    */
    template <class T>
    bool LinearObservationManager<T>::HasErrorMatrix() const
    {
#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        return true;
#else
        return false;
#endif
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    template <class state>
    void LinearObservationManager<T>
    ::ApplyOperator(const state& x,
                    observation& y) const
    {
        if (x.GetSize() == 0)
            return;

        if (operator_scaled_identity_)
        {
            Copy(x, y);
            Mlt(operator_diagonal_value_, y);
        }

        // Operator defined in a file.
        else
            Mlt(tangent_operator_matrix_, x, y);
    }


    //! Applies the tangent linear operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator at \a x.
    */
    template <class T>
    template <class state>
    void LinearObservationManager<T>
    ::ApplyTangentLinearOperator(const state& x,
                                 observation& y)
        const
    {
        ApplyOperator(x, y);
    }


    //! Linearized observation operator.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the linearized operator.
    */
    template <class T>
    T LinearObservationManager<T>
    ::GetTangentLinearOperator(int i, int j) const
    {
        if (operator_scaled_identity_)
            if (i == j)
                return operator_diagonal_value_;
            else
                return T(0);
        else
            return tangent_operator_matrix_(i, j);
    }


    //! Linearized observation operator.
    /*!
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the linearized
      operator.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetTangentLinearOperatorRow(int row,
                                  tangent_linear_operator_row&
                                  tangent_operator_row)
        const
    {
        if (operator_scaled_identity_)
        {
            tangent_operator_row.Reallocate(Nobservation_);
            tangent_operator_row.Zero();
            tangent_operator_row(row) = operator_diagonal_value_;
        }
        else
            GetRow(tangent_operator_matrix_, row, tangent_operator_row);
    }


    //! Linearized observation operator.
    /*!
      \return The matrix of the linearized operator.
    */
    template <class T>
    const typename LinearObservationManager<T>
    ::tangent_linear_operator& LinearObservationManager<T>
    ::GetTangentLinearOperator() const
    {
        return tangent_operator_matrix_;
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x.
    */
    template <class T>
    template <class state>
    void LinearObservationManager<T>
    ::ApplyAdjointOperator(const state& x,
                           observation& y) const
    {
        if (operator_scaled_identity_)
        {
            y = x;
            Mlt(operator_diagonal_value_, y);
        }
        else
            Mlt(1., SeldonTrans, tangent_operator_matrix_, x, 0., y);
    }


    //! Checks whether a BLUE correction is available.
    /*!
      \return True if a BLUE correction is available, false otherwise.
    */
    template <class T>
    bool LinearObservationManager<T>
    ::HasBLUECorrection() const
    {
        throw ErrorUndefined("LinearObservationManager"
                             "::HasBLUECorrection()");
    }


    //! Gets the BLUE correction.
    /*!
      \param[out] BLUE_correction BLUE correction vector.
    */
    template <class T>
    void LinearObservationManager<T>
    ::GetBLUECorrection(Vector<T>& BLUE_correction) const
    {
        throw ErrorUndefined("LinearObservationManager"
                             "::GetBLUECorrection(correction)");
    }


    //! Observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error covariance.
    */
    template <class T>
    T LinearObservationManager<T>::GetErrorVariance(int i, int j) const
    {
        if (i == j)
            return error_variance_value_;
        else
            return T(0);
    }


    //! Observation error covariance matrix.
    /*!
      \return The matrix of the observation error covariance.
    */
    template <class T>
    const typename LinearObservationManager<T>
    ::error_variance& LinearObservationManager<T>::GetErrorVariance() const
    {
        return error_variance_;
    }


    //! Inverse of the observation error covariance matrix.
    /*!
      \return Inverse of the matrix of the observation error covariance.
    */
    template <class T>
    const typename LinearObservationManager<T>
    ::error_variance& LinearObservationManager<T>::GetErrorVarianceInverse()
        const
    {
        return error_variance_inverse_;
    }



    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string LinearObservationManager<T>::GetName() const
    {
        return "LinearObservationManager";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void LinearObservationManager<T>::Message(string message)
    {
    }


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATION_MANAGER_LINEAROBSERVATIONMANAGER_CXX
#endif
