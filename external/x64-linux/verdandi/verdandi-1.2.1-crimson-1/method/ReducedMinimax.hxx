// Copyright (C) 2008-2010 INRIA
// Author(s): Vivien Mallet, Serhiy Zhuk
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


#ifndef VERDANDI_FILE_METHOD_REDUCEDMINIMAX_HXX


namespace Verdandi
{


    ////////////////////
    // REDUCEDMINIMAX //
    ////////////////////


    //! This class implements a reduced minimax filter.
    template <class T, class Model, class ObservationManager>
    class ReducedMinimax: public VerdandiBase
    {

    public:
        //! Type of the model state vector.
        typedef typename Model::state model_state;
        //! Type of the background error variance.
        typedef typename Model::state_error_variance
        model_state_error_variance;
        //! Type of the background error variance.
        typedef typename Model::error_variance
        model_error_variance;
        //! Type of the observation operator.
        typedef typename ObservationManager::tangent_linear_operator
        observation_tangent_linear_operator;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;

        //! Observation manager.
        ObservationManager observation_manager_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;

        //! Should the iterations be displayed?
        bool show_iteration_;
        //! Should the current time be displayed?
        bool show_time_;

        //! Should the filter be applied?
        bool with_filtering_;

        //! Model configuration file.
        string model_configuration_file_;
        //! Configuration file for the observation manager.
        string observation_configuration_file_;

        //! Dimension of the state.
        int Nstate_;

        //! Number of observations.
        int Nobservation_;

        /*** POD ***/

        /*! \brief Width of the simulation window (in number of time steps)
          for the snapshot recording. */
        int Nstep_snapshot_;
        /*! \brief Should \f$ H^T y \f$ and \f$ H^T H x \f$ be included in the
          snapshots? */
        bool with_Hty_HtHx_;
        //! Approximate number of observations per state snapshots.
        double Nobservation_step_per_snapshot_;

        //! Snapshots to be used for POD.
        Matrix<T, General, RowMajor> snapshot_;
        //! Current number of snapshots stored in matrix 'snapshot_'.
        int Nsnapshot_;

        //! Singular values.
        Vector<T> singular_value_;
        //! Left singular vectors.
        Matrix<T, General, RowMajor> left_singular_vector_;
        //! Right singular vectors.
        Matrix<T, General, RowMajor> right_singular_vector_;

        /*** Reduction ***/

        //! Reduction method: "none" or "pod".
        string reduction_method_;
        //! Maximum dimension of the reduced space.
        int Nprojection_max_;
        //! Acceptable relative quadratic error.
        T acceptable_error_;
        //! Dimension of the reduced space.
        int Nprojection_;
        //! Projection matrix.
        Matrix<T, General, RowMajor> projection_;
        //! Previous dimension of the reduced space.
        int Nprevious_projection_;
        //! Previous projection matrix.
        Matrix<T, General, RowMajor> previous_projection_;

        /*** Management of time iterations ***/

        //! Global iteration.
        int iteration_;

        //! Mode: 0 for snapshot recording for POD, 1 for error computation.
        int mode_;
        /*! \brief Number of iteration in a given sequence. After 'Nsnapshot_'
          iterations, one switches to either snapshot recording or error
          computation. */
        int inner_iteration_;
        //! Is the current sequence in the first window?
        bool first_sequence_;

        //! Starting time of a new sequence of error computation.
        double starting_time_;
        /*! \brief Snapshot of the full state vector used to restart before an
          error computation sequence. This full state vector is a snapshot at
          'starting_time_'. */
        model_state full_state_;

        /*** Filter ***/

        //! State estimator.
        model_state state_;

        //! Minimax gain.
        Matrix<T, General, RowMajor> G_;

        //! Is the model error variance composed of its diagonal part only?
        bool is_model_error_variance_diagonal_;
        //! Diagonal part of the model error.
        Vector<T> D_tilde_inv_;
        //! Number of modes in the square root $\widetilde Q^{\frac 12}$.
        int Nmode_Q_;
        /*! \brief Form of the model error matrix Q: 0 = no model error, 1 =
          identity times a scalar. */
        int model_error_description_;
        //! In case, the model error is the identity multiplied by a scalar.
        T model_error_scalar_;
        //! Conversion factor between standard deviation and bound.
        T bound_over_standard_deviation_;

        //! Systematic error in the initial condition.
        Vector<T> e0_;
        //! Systematic error.
        Vector<T> e_;

        //! Observation error.
        Matrix<T, General, RowMajor> R_inv_;
        //! Systematic error in the observations.
        T eta_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:

        /*** Constructor and destructor ***/

        ReducedMinimax();
        ~ReducedMinimax();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);

        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);

        void InitializeStep();

        void Forward();
        void FinalizeStep();
        void Finalize();

        void FilterInitialization();
        void Propagation();
        void ComputeTangentLinearModel(Matrix<T, General, RowMajor>& M);

        void SnapshotRecording();

        bool HasFinished();

        // Access methods.
        int GetMode() const;
        Matrix<T, General, RowMajor>& GetProjection();
        Matrix<T, General, RowMajor>& GetPreviousProjection();
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDMINIMAX_HXX
#endif
