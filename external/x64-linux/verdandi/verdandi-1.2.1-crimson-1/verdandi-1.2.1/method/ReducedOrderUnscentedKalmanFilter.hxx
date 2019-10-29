// Copyright (C) 2008-2010 INRIA
// Author(s): Marc Fragu, Philippe Moireau, Vivien Mallet
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


#ifndef VERDANDI_FILE_METHOD_REDUCEDORDERUNSCENTEDKALMANFILTER_HXX


#include "seldon/vector/VectorCollection.hxx"


namespace Verdandi
{


    ///////////////////////////////////////
    // REDUCEDORDERUNSCENTEDKALMANFILTER //
    ///////////////////////////////////////


    //! This class implements a reduced order unscented Kalman filter.
    template <class T, class Model, class ObservationManager>
    class ReducedOrderUnscentedKalmanFilter: public VerdandiBase
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
        //! Type of the sigma point vector.
        typedef Vector<T, VectFull, MallocAlloc<T> > sigma_point;
        //! Type of the state vector collection.
        typedef Vector<model_state, Collection> state_collection;
        //! Type of the observation vector collection.
        typedef Vector<observation, Collection> observation_collection;
        //! Type of the sigma point matrix.
        typedef Matrix<T, General, RowMajor, MallocAlloc<T> >
        sigma_point_matrix;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;
        //! Observation manager.
        ObservationManager observation_manager_;

        //! Matrix L in the P SVD decomposition.
        model_state_error_variance L_;
        //! Matrix U in the P SVD decomposition.
        sigma_point_matrix U_;
        //! Inverse of matrix U.
        sigma_point_matrix U_inv_;
        //! Innovation covariance.
        sigma_point_matrix C_;

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
        //! Dimension of the filtered state.
        int Nreduced_;
        //! Number of observations.
        int Nobservation_;
        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;
        //! Should resampling be perfomed.
        bool with_resampling_;
        //! Indicates how R is stored (matrix, matrix_inverse, vector).
        string observation_error_variance_;
        /*! \brief Should the computations be carried out on observations
          (false) or directly on innovations (true)? */
        bool innovation_computation_;
        bool with_innovation_;

        /*** Sigma-points ***/

        //! Choice of sigma-points.
        string sigma_point_type_;
        //! Matrix of sigma-points.
        sigma_point_matrix I_;
        //! Matrix of sigma-points (transposed).
        sigma_point_matrix I_trans_;
#if defined(VERDANDI_WITH_MPI)
        //! Transpose of [X_n^(*)].
        sigma_point_matrix X_i_trans_;
#else
        //! [X_n^(*)].
        model_state_error_variance X_i_;
        //! [X_n^(*)] (transposed).
        model_state_error_variance X_i_trans_;
#endif
        //! Coefficient vector associated with sigma-points.
        sigma_point D_alpha_;
        //! P_alpha^{V}.
        sigma_point_matrix D_v_;
        //! Boolean to indicate if the coefficients alpha are constants.
        bool alpha_constant_;
        //! alpha.
        T alpha_;
        //! Number of sigma-points.
        int Nsigma_point_;

        //! Process rank.
        int rank_;

#if defined(VERDANDI_WITH_MPI)

        /*** Parallel data ***/

        //! Process rank.
        int algorithm_;
        //! Number of processes.
        int Nprocess_;
        //! Number of local sigma-points.
        int Nlocal_sigma_point_;
        //! Local sigma-points sum.
        Vector<int> Nlocal_sigma_point_sum_;
        //! Local sigma-points.
        Vector<int> local_sigma_point_;
        //! Number of local columns of covariance matrix.
        int Nlocal_filtered_column_;
        //! Local sigma-points.
        Vector<int> local_filtered_column_;
        //! Local columns sum.
        Vector<int> Nlocal_filtered_column_sum_;
        //! Global matrix of sigma-points.
        sigma_point_matrix I_trans_global_;
        //! Contribution of process 0.
        double master_process_contribution_;
#endif

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        ReducedOrderUnscentedKalmanFilter();
        ~ReducedOrderUnscentedKalmanFilter();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);

        void InitializeStep();

        void Forward();
        void Analyze();

        void FinalizeStep();
        void Finalize();

        bool HasFinished();

        // Access methods.
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDORDERUNSCENTEDKALMANFILTER_HXX
#endif
