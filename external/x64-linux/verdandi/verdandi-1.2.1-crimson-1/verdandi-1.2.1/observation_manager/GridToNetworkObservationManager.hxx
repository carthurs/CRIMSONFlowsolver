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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_GRIDTONETWORKOBSERVATIONMANAGER_HXX

#include <iostream>
#include <list>

#ifndef OBSERVATION_AGGREGATOR
#define OBSERVATION_AGGREGATOR ObservationAggregator
#endif

#define _QUOTE(x) #x
#define QUOTE(x) _QUOTE(x)
#include QUOTE(OBSERVATION_AGGREGATOR.hxx)


namespace Verdandi
{


    /////////////////////////////////////
    // GRIDTONETWORKOBSERVATIONMANAGER //
    /////////////////////////////////////


    //! Observation operator that maps from a grid to a network.
    /*! The operator linearly interpolates from a regular grid to a list of
      locations (i.e., a network).
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class GridToNetworkObservationManager: public VerdandiBase
    {
    public:
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, RowSparse> tangent_linear_operator;
        //! Type of the observation error covariance matrix.
        typedef Matrix<T, General, RowSparse> error_variance;
        //! Type of a row of the tangent linear operator.
        typedef Vector<T> tangent_linear_operator_row;

        //! Type of the observation vector.
        typedef Vector<T> observation;
        //! Type of the observation vector.
        typedef Vector<T> observation_vector;
        //! Type of the observation vector 2.
        typedef Vector2<T> observation_vector2;
        //! Type of the observation vector 3.
        typedef Vector3<T> observation_vector3;

        //! Type of the variable vector.
        typedef Vector<int> variable_vector;
        //! Type of the variable vector 2.
        typedef Vector2<int> variable_vector2;
        //! Type of the variable vector 3.
        typedef Vector3<int> variable_vector3;

        //! Type of the index vector.
        typedef Vector<int> index_vector;
        //! Type of the index vector 2.
        typedef Vector2<int> index_vector2;
        //! Type of the index vector 3.
        typedef Vector3<int> index_vector3;

        //! Type of the time vector.
        typedef Vector<double> time_vector;
        //! Type of the time vector 2.
        typedef Vector2<double> time_vector2;
        //! Type of the time vector 3.
        typedef Vector3<double> time_vector3;

    protected:

        /*** Observations file structure ***/

        //! File that stores the observations.
        string observation_file_;
        //! How are stored the observations.
        string observation_type_;
        //! Number total of observations at current time.
        int Nobservation_;
        //! Size in byte of an observations vector.
        size_t Nbyte_observation_;
        //! Period with which observations are available.
        double Delta_t_;
        //! Period with which available observations are actually loaded.
        int Nskip_;
        //! Duration during which observations are assimilated.
        double final_time_;

        /*** Observation times ***/

        //! Requested time.
        double time_;
        //! Available observation time of the time interval.
        time_vector available_time_;
        //! Contribution associated with available observations.
        Vector<double> contribution_;
        //! Observations aggregator.
        OBSERVATION_AGGREGATOR<T> observation_aggregator_;

        /*** GridToNetwork parameters ***/

        //! Index along x.
        Vector<T> location_x_;
        //! Index along y.
        Vector<T> location_y_;

        //! Interpolation indices for all locations.
        Matrix<int> interpolation_index_;

        //! Interpolation weights for all locations.
        Matrix<T> interpolation_weight_;

        //! Interpolation indices for active locations.
        Matrix<int> active_interpolation_index_;

        //! Interpolation weights for active locations.
        Matrix<T> active_interpolation_weight_;

        //! Observation error variance.
        T error_variance_value_;

        /*** Model domain ***/

        //! Number of points along x in the model.
        int Nx_model_;
        //! Number of points along y in the model.
        int Ny_model_;

    public:
        // Constructors and destructor.
        GridToNetworkObservationManager();
        template <class Model>
        GridToNetworkObservationManager(const Model& model,
                                        string configuration_file);
        ~GridToNetworkObservationManager();

        // Initialization.
        template <class Model>
        void Initialize(const Model& model, string configuration_file);

        void DiscardObservation(bool discard_observation);
        void SetAllActive();

        int CreateTrack();
        void SetTrack(int track);

        template <class Model>
        void SetTime(const Model& model, double time);
        void SetTime(double time);
        void SetAvailableTime(double time, time_vector& available_time) const;
        void SetAvailableTime(double time_inf, double time_sup,
                              time_vector& available_time) const;


        ////////////////////////////
        // FLATTENED OBSERVATIONS //
        ////////////////////////////


        /*** Gets observations ***/

        void GetFlattenedObservation(double time,
                                     observation_vector& observation);
        void GetFlattenedObservation(double time_inf, double time_sup,
                                     observation_vector& observation);
        void GetFlattenedObservation(observation_vector& observation);
        void GetFlattenedObservation(const time_vector& available_time,
                                     observation_vector& observation);

        /*** Gets observations and associated variables ***/

        void GetFlattenedObservation(double time,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(double time_inf, double time_sup,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(variable_vector& observation_variable,
                                     observation_vector& observation);
        void GetFlattenedObservation(const time_vector& available_time,
                                     variable_vector& observation_variable,
                                     observation_vector& observation);

        /*** Gets observations, associated variables and associated index ***/

        void GetFlattenedObservation(double time,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(double time_inf, double time_sup,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);
        void GetFlattenedObservation(const time_vector& available_time,
                                     variable_vector& observation_variable,
                                     index_vector& observation_index,
                                     observation_vector& observation);


        /////////////////////////////
        // AGGREGATED OBSERVATIONS //
        /////////////////////////////


        /*** Gets observations ***/

        void GetAggregatedObservation(double time,
                                      observation_vector& observation);
        void GetAggregatedObservation(double time_inf, double time_sup,
                                      observation_vector& observation);
        void GetAggregatedObservation(observation_vector& observation);
        void GetAggregatedObservation(const time_vector& available_time,
                                      observation_vector& observation);

        /*** Gets observations and associated variables ***/

        void GetAggregatedObservation(double time,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(double time_inf, double time_sup,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(variable_vector& observation_variable,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(const time_vector& available_time,
                                      variable_vector& observation_variable,
                                      observation_vector2& observation2);

        /*** Gets observations, associated variables and associated index ***/

        void GetAggregatedObservation(double time,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(double time_inf, double time_sup,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);
        void GetAggregatedObservation(const time_vector& available_time,
                                      variable_vector& observation_variable,
                                      index_vector2& observation_index2,
                                      observation_vector2& observation2);


        //////////////////////
        // RAW OBSERVATIONS //
        //////////////////////


        /*** Gets observations ***/

        void GetRawObservation(double time,
                               observation_vector2& observation2);
        void GetRawObservation(double time_inf, double time_sup,
                               observation_vector2& observation2);
        void GetRawObservation(observation_vector2& observation2);
        void GetRawObservation(const time_vector& available_time,
                               observation_vector2& observation2);

        /*** Gets observations and associated variables ***/

        void GetRawObservation(double time,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(double time_inf, double time_sup,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(variable_vector2& observation_variable2,
                               observation_vector3& observation3);
        void GetRawObservation(const time_vector& available_time,
                               variable_vector2& observation_variable2,
                               observation_vector3& observation3);

        /*** Gets observations, associated variables and associated index ***/

        void GetRawObservation(double time,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(double time_inf, double time_sup,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);
        void GetRawObservation(const time_vector& available_time,
                               variable_vector2& observation_variable2,
                               index_vector3& observation_index3,
                               observation_vector3& observation3);


        ///////////////////////////////
        // READ OBSERVATIONS METHODS //
        ///////////////////////////////


        void ReadObservationVariable(const time_vector& available_time,
                                     variable_vector2& observation_variable2)
            const;
        void ReadObservation(const time_vector& available_time,
                             const variable_vector2& observation_variable2,
                             observation_vector3& observation3) const;
        void ReadObservation(const time_vector& available_time,
                             observation_vector2& observation2) const;
        void ReadObservation(double time, int variable,
                             observation_vector& observation) const;
        void ReadObservationIndex(const time_vector& available_time, const
                                  variable_vector2& observation_variable2,
                                  index_vector3& observation_index3) const;


        /////////////////
        // OBSERVATION //
        ////////////////


        void GetObservation(observation& observation);


        ////////////////
        // INNOVATION //
        ////////////////


        template <class state>
        void GetInnovation(const state& x, observation& innovation);


        ////////////
        // ACCESS //
        ////////////


        bool HasObservation() const;
        bool HasObservation(double time);
        int GetNobservation() const;
        bool IsOperatorSparse() const;
        bool IsErrorSparse() const;
        bool HasErrorMatrix() const;


        ///////////////
        // OPERATORS //
        ///////////////


        template <class state>
        void ApplyOperator(const state& x, observation& y) const;

        template <class state>
        void ApplyTangentLinearOperator(const state& x, observation& y) const;
        T GetTangentLinearOperator(int i, int j) const;
        void GetTangentLinearOperatorRow(int row, tangent_linear_operator_row&
                                         tangent_operator_row) const;
        const tangent_linear_operator& GetTangentLinearOperator() const;

        template <class state>
        void ApplyAdjointOperator(const state& x, observation& y) const;

        bool HasBLUECorrection() const;
        void GetBLUECorrection(Vector<T>& BLUE_correction) const;

        T GetErrorVariance(int i, int j) const;
        const error_variance& GetErrorVariance() const;
        const error_variance& GetErrorVarianceInverse() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATION_MANAGER_GRIDTONETWORKOBSERVATIONMANAGER_HXX
#endif
