// Copyright (C) 2008-2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_GRIDTONETWORKOBSERVATIONMANAGER_CXX


#include "GridToNetworkObservationManager.hxx"
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
    GridToNetworkObservationManager<T>
    ::GridToNetworkObservationManager()
    {
    }


    //! Main constructor.
    /*! It defines the operator for a 2D regular grid.
      \param[in] model model.
      \param[in] configuration_file configuration_file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class T>
    template <class Model>
    GridToNetworkObservationManager<T>
    ::GridToNetworkObservationManager(const Model& model,
                                      string configuration_file)
    {
    }


    //! Destructor.
    template <class T>
    GridToNetworkObservationManager<T>
    ::~GridToNetworkObservationManager()
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
    void GridToNetworkObservationManager<T>
    ::Initialize(const Model& model, string configuration_file)
    {
        observation_aggregator_.Initialize(configuration_file);

        VerdandiOps configuration(configuration_file);

        MessageHandler::AddRecipient("grid_to_network_observation_manager",
                                     reinterpret_cast<void*>(this),
                                     GridToNetworkObservationManager
                                     ::StaticMessage);

        //! First abscissa.
        double x_min_model = model.GetXMin();
        //! First ordinate.
        double y_min_model = model.GetYMin();

        //! Space step along x.
        double Delta_x_model = model.GetDeltaX();
        //! Space step along y.
        double Delta_y_model = model.GetDeltaY();

        Nx_model_ = model.GetNx();
        Ny_model_ = model.GetNy();
        Nbyte_observation_ = Nx_model_ * Ny_model_ * sizeof(double);

        configuration.SetPrefix("observation.");
        configuration.Set("file", observation_file_);
        configuration.Set("Delta_t", "v > 0", Delta_t_);
        configuration.Set("Nskip", "v > 0", Nskip_);
        configuration.Set("final_time", "", numeric_limits<double>::max(),
                          final_time_);

        time_ = -numeric_limits<double>::max();

        configuration.Set("error.variance", "v > 0", error_variance_value_);

        configuration.SetPrefix("observation.location.");

        vector<string> observation_location_vector;
        configuration.Set("observation_location",
                          observation_location_vector);
        int value;
        for (int i = 0; i < int(observation_location_vector.size() - 1);
             i += 2)
        {
            to_num(observation_location_vector[i], value);
            location_x_.PushBack(value);
            to_num(observation_location_vector[i + 1], value);
            location_y_.PushBack(value);
        }

        Nobservation_ = int(location_x_.GetSize());

        for (int i = 0; i < Nobservation_; i++)
            if (location_x_(i) < 0 || location_x_(i) >= Nx_model_)
                throw "Index along x should be in [0, "
                    + to_str(Nx_model_ - 1)
                    + "], but " + to_str(location_x_(i)) + " was provided.";
        for (int i = 0; i < Nobservation_; i++)
            if (location_y_(i) < 0 || location_y_(i) >= Ny_model_)
                throw "Index along y should be in [0, "
                    + to_str(Ny_model_ - 1)
                    + "], but " + to_str(location_y_(i)) + " was provided.";

        interpolation_index_.Reallocate(Nobservation_, 4);
        interpolation_weight_.Reallocate(Nobservation_, 4);

        int i, j;
        T weight_x, weight_y;
        for (int location = 0; location < Nobservation_; location++)
        {
            i = int((location_x_(location) - x_min_model) / Delta_x_model);
            j = int((location_y_(location) - y_min_model) / Delta_y_model);
            interpolation_index_(location, 0) = i * Ny_model_ + j;
            interpolation_index_(location, 1) = i * Ny_model_ + j + 1;
            interpolation_index_(location, 2) = (i + 1) * Ny_model_ + j;
            interpolation_index_(location, 3) = (i + 1) * Ny_model_ + j + 1;

            weight_x = (location_x_(location) - x_min_model
                        - T(i) * Delta_x_model) / Delta_x_model;
            weight_y = (location_y_(location) - y_min_model
                        - T(j) * Delta_y_model) / Delta_y_model;
            interpolation_weight_(location, 0) = (T(1) - weight_x)
                * (T(1) - weight_y);
            interpolation_weight_(location, 1) = (T(1) - weight_x) * weight_y;
            interpolation_weight_(location, 2) = weight_x * (T(1) - weight_y);
            interpolation_weight_(location, 3) = weight_x * weight_y;
        }

        // Sets all active by default.
        SetAllActive();

        int expected_file_size;
        expected_file_size = Nbyte_observation_ *
            int(final_time_ / (Delta_t_ * Nskip_));

        int file_size;
        ifstream file_stream;
        file_stream.open(observation_file_.c_str());

#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("GridToNetworkObservationManager"
                          "::Initialize(model, configuration_file)",
                          "Unable to open file \""
                          + observation_file_ + "\".");
#endif

        file_stream.seekg( 0 , ios_base::end );
        file_size = file_stream.tellg() ;
        file_stream.close();

        if (expected_file_size > file_size)
            throw IOError("GridToNetworkObservationManager"
                          "::Initialize(model, configuration_file)",
                          "Too few available observations, the size of \""
                          + observation_file_ + "\" must be greater than "
                          + to_str(expected_file_size) + ".");
    }


    //! Sets all observation locations as active.
    template <class T>
    void GridToNetworkObservationManager<T>::SetAllActive()
    {
        active_interpolation_index_ = interpolation_index_;
        active_interpolation_weight_ = interpolation_weight_;
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
    void GridToNetworkObservationManager<T>
    ::DiscardObservation(bool discard_observation)
    {
        observation_aggregator_.DiscardObservation(discard_observation);
    }


    //! Creates a new track.
    /*!
      \return The index of the new track.
    */
    template <class T>
    int GridToNetworkObservationManager<T>
    ::CreateTrack()
    {
        return observation_aggregator_.CreateTrack();
    }


    //! Sets the track to a given track.
    /*!
      \param[in] track the given track.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
    ::SetTime(const Model& model, double time)
    {
        SetTime(time);
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::SetTime(double time)
    {
        if (time_ == time)
            return;

        time_ = time;
        SetAvailableTime(time_, available_time_);
        observation_aggregator_.
            Contribution(time_, available_time_, contribution_);
    }


    //! Sets the available observation times at a given time.
    /*!
      \param[in] time the given time.
      \param[out] available_time the available observation times.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::SetAvailableTime(double time,
                       time_vector&
                       available_time) const
    {
        double time_inf, time_sup;
        int selection_policy;
        observation_aggregator_.GetContributionInterval(time, time_inf,
                                                        time_sup,
                                                        selection_policy);
        SetAvailableTime(time_inf, time_sup, available_time);

        Logger::Log<3>(*this, to_str(time) + ", [" + to_str(time_inf) + " " +
                       to_str(time_sup) + "], {" + to_str(available_time) +
                       "}\n");
    }


    //! Sets available observation times at a given time interval.
    /*!
      \param[in] time_inf lower bound of the given time interval.
      \param[in] time_sup upper bound (excluded) of the given time interval.
      \param[out] available_time the available observation times.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::SetAvailableTime(double time_inf, double time_sup,
                       time_vector&
                       available_time) const
    {
        available_time.Clear();

        time_sup = time_sup < final_time_ ? time_sup : final_time_;

        double period = Delta_t_ * Nskip_;
        double available_time_0 = floor(time_inf / period) * period;
        if (available_time_0 == time_inf)
            available_time.PushBack(available_time_0);
        available_time_0 += period;
        for (double t = available_time_0; t < time_sup; t += period)
            available_time.PushBack(t);
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
    ::GetFlattenedObservation(const time_vector& available_time,
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(double time,
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(double time_inf, double time_sup,
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(double time,
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(double time_inf, double time_sup,
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(double time,
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(double time_inf, double time_sup,
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(variable_vector& observation_variable,
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
    void GridToNetworkObservationManager<T>
    ::GetAggregatedObservation(const time_vector& available_time,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(double time,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(double time_inf, double time_sup,
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
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(const time_vector& available_time,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(double time,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(double time_inf, double time_sup,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(variable_vector2& observation_variable2,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(const time_vector& available_time,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(double time,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(double time_inf, double time_sup,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(variable_vector2& observation_variable2,
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
    void GridToNetworkObservationManager<T>
    ::GetRawObservation(const time_vector& available_time,
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
    void GridToNetworkObservationManager<T>
    ::ReadObservationVariable(const time_vector& available_time,
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
    void GridToNetworkObservationManager<T>
    ::ReadObservation(const time_vector& available_time,
                      const variable_vector2& observation_variable2,
                      observation_vector3& observation3)
        const
    {
        int Nvariable, Nt;
        Nt = available_time.GetSize();
        observation3.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
        {
            Nvariable = observation_variable2(h).GetSize();
            observation3.Reallocate(h, Nvariable);
            for (int v = 0; v < Nvariable; v++)
                ReadObservation(available_time(h),
                                observation_variable2(h, v),
                                observation3.GetVector(h, v));
        }
    }


    //! Builds observations associated with given times.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation2 the observations.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::ReadObservation(const time_vector& available_time,
                      observation_vector2& observation2)
        const
    {
        int Nt = available_time.GetSize();
        observation2.Reallocate(Nt);
        for (int h = 0; h < Nt; h++)
            ReadObservation(available_time(h), 0, observation2(h));
    }


    //! Reads observation from observation file given a time and a variable.
    /*!
      \param[in] time the time.
      \param[in] variable the variable.
      \param[out] observation the observations.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::ReadObservation(double time, int variable,
                      observation_vector& observation)
        const
    {
        observation.Reallocate(Nobservation_);
        Matrix<T> input_data(Nx_model_, Ny_model_);
        int iteration;

        iteration = int(time / (Delta_t_ * Nskip_));

        ifstream file_stream;
        file_stream.open(observation_file_.c_str());

#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("GridToNetworkObservationManager" \
                          "::LoadObservation(model)",
                          string("Unable to open file \"")
                          + observation_file_ + "\".");
#endif

        streampos position = iteration * Nbyte_observation_;
        file_stream.seekg(position);
        input_data.Read(file_stream, false);
        file_stream.close();

        for (int i = 0; i < Nobservation_; i++)
            observation(i) = input_data(location_x_(i), location_y_(i));

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
    void GridToNetworkObservationManager<T>
    ::ReadObservationIndex(const time_vector& available_time,
                           const variable_vector2& observation_variable2,
                           index_vector3& observation_index3)
        const
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


    /////////////////
    // OBSERVATION //
    /////////////////


    //! Gets observation.
    /*!
      \param[out] observation observation vector.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
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
    void GridToNetworkObservationManager<T>
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
    bool GridToNetworkObservationManager<T>::HasObservation(double time)
    {
        SetTime(time);
        return available_time_.GetSize() != 0
            && !is_equal(Norm1(contribution_), 0.);
    }


    //! Indicates if some observations are available at current time.
    template <class T>
    bool GridToNetworkObservationManager<T>::HasObservation() const
    {
        return available_time_.GetSize() != 0
            && !is_equal(Norm1(contribution_), 0.);
    }


    //! Gets Nobservation_ value.
    /*!
      \return The total number of observation at current time.
    */
    template <class T>
    int GridToNetworkObservationManager<T>::GetNobservation() const
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
    bool GridToNetworkObservationManager<T>::IsOperatorSparse() const
    {
#ifdef VERDANDI_TANGENT_OPERATOR_SPARSE
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
    bool GridToNetworkObservationManager<T>::IsErrorSparse() const
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
    bool GridToNetworkObservationManager<T>::HasErrorMatrix() const
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
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class T>
    template <class state>
    void GridToNetworkObservationManager<T>
    ::ApplyOperator(const state& x, observation& y) const
    {
        y.Reallocate(active_interpolation_index_.GetM());
        y.Zero();
        for (int i = 0; i < y.GetLength(); i++)
            y(i) = active_interpolation_weight_(i, 0)
                * x(active_interpolation_index_(i, 0))
                + active_interpolation_weight_(i, 1)
                * x(active_interpolation_index_(i, 1))
                + active_interpolation_weight_(i, 2)
                * x(active_interpolation_index_(i, 2))
                + active_interpolation_weight_(i, 3)
                * x(active_interpolation_index_(i, 3));
    }


    //! Applies the tangent linear operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator at \a x. It is
      resized if needed.
    */
    template <class T>
    template <class state>
    void GridToNetworkObservationManager<T>
    ::ApplyTangentLinearOperator(const state& x, observation& y) const
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
    T GridToNetworkObservationManager<T>
    ::GetTangentLinearOperator(int i, int j) const
    {
        if (j == active_interpolation_index_(i, 0))
            return active_interpolation_weight_(i, 0);
        else if (j == active_interpolation_index_(i, 1))
            return active_interpolation_weight_(i, 1);
        else if (j == active_interpolation_index_(i, 2))
            return active_interpolation_weight_(i, 2);
        else if (j == active_interpolation_index_(i, 3))
            return active_interpolation_weight_(i, 3);
        else
            return T(0);
    }


    //! Linearized observation operator.
    /*!
      \param[in] row row index.
      \param[out] tangent_operator_row the row \a row of the linearized
      operator.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::GetTangentLinearOperatorRow(int row,
                                  tangent_linear_operator_row&
                                  tangent_operator_row)
        const
    {
        for (int i = 0; i < tangent_operator_row.GetLength(); i++)
            tangent_operator_row(i) = GetTangentLinearOperator(row, i);
    }


    //! Linearized observation operator.
    /*!
      \return The matrix of the linearized operator.
    */
    template <class T>
    const typename
    GridToNetworkObservationManager<T>::tangent_linear_operator&
    GridToNetworkObservationManager<T>
    ::GetTangentLinearOperator() const
    {
        throw ErrorUndefined(
            "GridToNetworkObservationManager::GetTangentLinearOperator()");
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class T>
    template <class state>
    void GridToNetworkObservationManager<T>
    ::ApplyAdjointOperator(const state& x, observation& y) const
    {
        y.Reallocate(Nx_model_ * Ny_model_);
        y.Zero();
        for (int i = 0; i < Nx_model_ * Ny_model_; i++)
            for (int j = 0; j < active_interpolation_index_.GetM(); j++)
                if (i == active_interpolation_index_(j, 0))
                    y(i) += active_interpolation_weight_(j, 0) * x(j);
                else if (i == active_interpolation_index_(j, 1))
                    y(i) += active_interpolation_weight_(j, 1) * x(j);
                else if (i == active_interpolation_index_(j, 2))
                    y(i) += active_interpolation_weight_(j, 2) * x(j);
                else if (i == active_interpolation_index_(j, 3))
                    y(i) += active_interpolation_weight_(j, 3) * x(j);
    }


    //! Checks whether a BLUE correction is available.
    /*!
      \return True if a BLUE correction is available, false otherwise.
    */
    template <class T>
    bool GridToNetworkObservationManager<T>
    ::HasBLUECorrection() const
    {
        throw ErrorUndefined(
            "GridToNetworkObservationManager::HasBLUECorrection()");
    }


    //! Gets the BLUE correction.
    /*!
      \param[out] BLUE_correction BLUE correction vector.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::GetBLUECorrection(Vector<T>& BLUE_correction) const
    {
        throw ErrorUndefined(
            "GridToNetworkObservationManager::GetBLUECorrection(correction)");
    }


    //! Observation error covariance.
    /*!
      \param[in] i row index.
      \param[in] j column index.
      \return The element (\a i, \a j) of the observation error covariance.
    */
    template <class T>
    T GridToNetworkObservationManager<T>
    ::GetErrorVariance(int i, int j) const
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
    const typename
    GridToNetworkObservationManager<T>::error_variance&
    GridToNetworkObservationManager<T>
    ::GetErrorVariance() const
    {
        throw ErrorUndefined(
            string("GridToNetworkObservationManager")
            +"::GetObservationErrorVariance()");
    }


    //! Returns the inverse of the observation error covariance matrix.
    /*!
      \return The inverse of the matrix of the observation error covariance.
    */
    template <class T>
    const typename GridToNetworkObservationManager<T>::error_variance&
    GridToNetworkObservationManager<T>::GetErrorVarianceInverse() const
    {
        throw ErrorUndefined(
            "const typename"
            "GridToNetworkObservationManager::error_variance& "
            "GridToNetworkObservationManager::GetErrorVarianceInverse()"
            " const");
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string GridToNetworkObservationManager<T>::GetName() const
    {
        return "GridToNetworkObservationManager";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void GridToNetworkObservationManager<T>
    ::Message(string message)
    {
    }


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATION_MANAGER_GRIDTONETWORKOBSERVATIONMANAGER_CXX
#endif
