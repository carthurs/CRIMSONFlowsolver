// Copyright (C) 2008-2010 INRIA
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


#ifndef VERDANDI_FILE_METHOD_REDUCEDORDEREXTENDEDKALMANFILTER_HXX


namespace Verdandi
{


    //////////////////////////////////////
    // REDUCEDORDEREXTENDEDKALMANFILTER //
    //////////////////////////////////////


    //! This class implements a reduced order extended Kalman filter.
    template <class T, class Model, class ObservationManager>
    class ReducedOrderExtendedKalmanFilter: public VerdandiBase
    {

    public:
        //! Type of a row of the background error variance.
        typedef typename Model::state_error_variance_row
        model_state_error_variance_row;
        //! Type of the model state vector.
        typedef typename Model::state model_state;
        //! Type of the model/observation crossed matrix.
        typedef typename Model::matrix_state_observation
        matrix_state_observation;
        //! Type of the background error variance.
        typedef typename Model::state_error_variance
        model_state_error_variance;
        //! Type of the tangent linear model.
        typedef typename Model::tangent_linear_operator
        model_tangent_linear_operator;
        //! Type of the tangent linear observation operator.
        typedef typename ObservationManager
        ::tangent_linear_operator observation_tangent_linear_operator;
        //! Type of the observation error variance.
        typedef typename ObservationManager
        ::error_variance observation_error_variance;
        //! Type of a row of the tangent linear observation operator.
        typedef typename ObservationManager::tangent_linear_operator_row
        observation_tangent_linear_operator_row;
        //! Type of the observation vector.
        typedef typename ObservationManager::observation
        observation;
        //! Type of dense matrix.
        typedef Matrix<T, General, RowMajor> dense_matrix;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;
        //! Observation manager.
        ObservationManager observation_manager_;
        //! Matrix L in the P SVD decomposition.
        dense_matrix L_;
        //! Matrix U in the P SVD decomposition.
        dense_matrix U_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;
        //! Path to the configuration file for the observation manager.
        string observation_configuration_file_;

        //! Display options.
        map<string, bool> option_display_;
        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;
        //! Dimension of the filtered state.
        int Nreduced_;
        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;
        //! Indicates how R is stored (matrix, matrix_inverse, vector).
        string observation_error_variance_;
        //! Current time.
        double time_;

#if defined(VERDANDI_WITH_MPI)

        /*** Parallel data ***/

        //! Process rank.
        int rank_;
        //! Number of processes.
        int Nprocess_;
        //! Number of local column in L_.
        int Nlocal_reduced_;
        //! Local column of L_.
        Vector<int> local_reduced_column_;
        //! Local columns sum.
        Vector<int> Nlocal_reduced_column_sum_;
        //! Parameter displs relative to the first MPI Allgatherv call.
        int *displacement_gather_1_;
        //! Parameter recvcounts relative to the first MPI Allgatherv call.
        int *recvcount_gather_1_;
        //! Parameter displs relative to the second MPI Allgatherv call.
        int *displacement_gather_2_;
        //! Parameter recvcounts relative to the second Allgatherv call.
        int *recvcount_gather_2_;
        //! Parameter displs relative to the third MPI Allgatherv call.
        int *displacement_gather_3_;
        //! Parameter recvcounts relative to the third MPI Allgatherv call.
        int *recvcount_gather_3_;
#endif


        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        ReducedOrderExtendedKalmanFilter();
        ~ReducedOrderExtendedKalmanFilter();

        /*** Methods ***/

        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);
        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);

        void InitializeStep();

        void Forward();
        void Analyze();

        void FinalizeStep();
        void Finalize();

        void PropagateCovarianceMatrix();

        bool HasFinished();

        // Access methods.
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();
        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDORDEREXTENDEDKALMANFILTER_HXX
#endif
