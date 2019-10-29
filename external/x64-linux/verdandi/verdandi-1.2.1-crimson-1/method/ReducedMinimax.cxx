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


#ifndef VERDANDI_FILE_METHOD_REDUCEDMINIMAX_CXX


#include "ReducedMinimax.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    template <class T, class Model, class ObservationManager>
    ReducedMinimax<T, Model, ObservationManager>::ReducedMinimax()
    {

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_,
                                     Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ReducedMinimax::StaticMessage);
    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    ReducedMinimax<T, Model, ObservationManager>::~ReducedMinimax()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager);

    }


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("reduced_minimax.");


        /*********************************
         * Model and observation manager *
         *********************************/


        configuration.Set("model.configuration_file",
                          model_configuration_file_);
        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        Nstate_ = model_.GetNstate();

        configuration.Set("observation_manager.configuration_file",
                          observation_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);


        /***************************
         * Reads the configuration *
         ***************************/


        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.show_iteration", show_iteration_);
        // Should current time be displayed on screen?
        configuration.Set("display.show_time", show_time_);

        iteration_ = 0;

        /*** Reduction options ***/

        configuration.Set("reduction_method", "ops_in(v, {'none', 'pod'})",
                          reduction_method_);
        if (reduction_method_ == "pod")
        {
            configuration.Set("pod.Nprojection_max", "v > 0",
                              Nprojection_max_);
            configuration.Set("pod.acceptable_error",
                              "v >= 0 and v <= 1",
                              acceptable_error_);
            configuration.Set("pod.Nstep_snapshot", "v > 1", Nstep_snapshot_);
            configuration.Set("pod.with_Hty_HtHx", with_Hty_HtHx_);
            if (with_Hty_HtHx_)
                configuration.Set("pod.Nobservation_step_per_snapshot",
                                  "v >= 0.", Nobservation_step_per_snapshot_);
        }

        /*** Model options ***/

        configuration.Set("model.bound_over_standard_deviation",
                          bound_over_standard_deviation_);

        /*** Filter options ***/

        configuration.Set("with_filtering", with_filtering_);

        configuration.Set("model_error.diagonal",
                          is_model_error_variance_diagonal_);
        T diagonal_part;
        if (configuration.Is<T>("model_error.diagonal_part"))
        {
            configuration.Set("model_error.diagonal_part", diagonal_part);
            D_tilde_inv_.Reallocate(Nstate_);
            D_tilde_inv_.Fill(diagonal_part);
        }
        else
            configuration.Set("model_error.diagonal_part", D_tilde_inv_);
        // 'D_tilde_inv_' is provided as a variance.
        for (int i = 0; i < D_tilde_inv_.GetLength(); i++)
            D_tilde_inv_(i) = bound_over_standard_deviation_
                * sqrt(D_tilde_inv_(i));
        // 'D_tilde_inv_' is inverted.
        for (int i = 0; i < D_tilde_inv_.GetLength(); i++)
            D_tilde_inv_(i) = 1. / D_tilde_inv_(i);

        // Systematic error in the initial condition.
        T error;
        if (configuration.Is<T>("systematic_error.initial_condition"))
        {
            configuration.Set("systematic_error.initial_condition", error);
            e0_.Reallocate(Nstate_);
            e0_.Fill(error);
        }
        else
            configuration.Set("systematic_error.initial_condition", e0_);

        // Systematic error in the model error.
        e_.Reallocate(Nstate_);
        e_.Fill(0.);

        // Systematic error in the observation error.
        configuration.Set("systematic_error.observation_error", eta_);

        /*** Ouput saver ***/

        configuration.SetPrefix("reduced_minimax.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("state_forecast");
        output_saver_.Empty("state_analysis");
        output_saver_.Empty("reduced_state_analysis");
        output_saver_.Empty("projected_state");
        output_saver_.Empty("projection");
        output_saver_.Empty("minimax_gain");
        output_saver_.Empty("minimax_gain_model_part");
        output_saver_.Empty("snapshot");
        output_saver_.Empty("singular_value");
        output_saver_.Empty("left_singular_vector");
        output_saver_.Empty("right_singular_vector");
        output_saver_.Empty("model_error_variance_sqrt");
        output_saver_.Empty("reduced_tangent_linear_model");
        output_saver_.Empty("observation");
        output_saver_.Empty("reduced_tangent_linear_observation_operator");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("reduced_minimax.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }


        /*************
         * Algorithm *
         *************/

        if (reduction_method_ == "pod")
            // Starting with in a POD sequence.
            mode_ = 0;
        else
        {
            Nprojection_ = Nstate_;
            projection_.Reallocate(Nprojection_, Nprojection_);
            projection_.SetIdentity();
            previous_projection_ = projection_;
            Nprevious_projection_ = previous_projection_.GetM();
            mode_ = 1;
        }
        inner_iteration_ = 0;
        first_sequence_ = true;

        /*** Filter ***/

        if (show_time_)
            Logger::StdOut(*this,
                           "Starting time: " + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this,
                            "Starting time: " + to_str(model_.GetTime()));
        if (show_iteration_)
            Logger::StdOut(*this, "Initialization");
        else
            Logger::Log<-3>(*this, "Initialization");

        if (initialize_model)
            MessageHandler::Send(*this, "model", "initial condition");
    }


    //! Initializes a step.
    /*! Initializes a step for the model. */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        if (mode_ == 1 && inner_iteration_ == 0)
        {
            if (show_time_ || show_iteration_)
                Logger::StdOut(*this, "Starting filtering sequence");
            else
                Logger::Log<-3>(*this, "Starting filtering sequence");

            /*** Filter ***/

            if (first_sequence_)
            {
                FilterInitialization();

                MessageHandler::Send(*this, "model", "state_analysis");
                MessageHandler::Send(*this, "observation_manager",
                                     "state_analysis");
                MessageHandler::Send(*this, "driver", "state_analysis");
                MessageHandler::Send(*this, "driver",
                                     "reduced_state_analysis");
            }

            // In other windows, there is no such special case. The error
            // associated with the first step in the window is the error
            // of the last step in the previous window.
            first_sequence_ = false;
        }

        model_.InitializeStep();

        if (mode_ == 0 && inner_iteration_ == 0)
        {
            // Saves the initial state for the next sequence of error
            // computation.
            starting_time_ = model_.GetTime();
            model_.GetFullState(full_state_);

            model_state state;
            model_.GetState(state);

            // Allocates the predicted size for the snapshots.
            if (with_Hty_HtHx_)
                snapshot_.Reallocate(Nstate_, Nstep_snapshot_
                                     + int(Nobservation_step_per_snapshot_
                                           * double(2 * Nstep_snapshot_)
                                           + 0.5));
            else
                snapshot_.Reallocate(Nstate_, Nstep_snapshot_);
            SetCol(state, 0, snapshot_);
            Nsnapshot_ = 1;
            inner_iteration_++;

            if (show_time_ || show_iteration_)
                Logger::StdOut(*this, "Starting POD sequence");
            else
                Logger::Log<-3>(*this, "Starting POD sequence");
        }

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        if (show_time_)
            Logger::StdOut(*this, "Time: " + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this,
                            "Time: " + to_str(model_.GetTime()));
        if (show_iteration_)
            Logger::StdOut(*this, "Iteration " + to_str(iteration_) + " -> "
                           + to_str(iteration_ + 1));
        else
            Logger::Log<-3>(*this, "Iteration " + to_str(iteration_) + " -> "
                            + to_str(iteration_ + 1));

        if (mode_ == 0)
            SnapshotRecording();
        else
        {
            Propagation();

            MessageHandler::Send(*this, "model", "state_analysis");
            MessageHandler::Send(*this, "observation_manager",
                                 "state_analysis");
            MessageHandler::Send(*this, "driver", "state_analysis");
            MessageHandler::Send(*this, "driver", "reduced_state_analysis");
        }

        iteration_++;

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Finalizes a step for the model.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Initialization of the filter.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::FilterInitialization()
    {
        MessageHandler::Send(*this, "all", "::FilterInitialization begin");

        if (!with_filtering_)
        {
            // Without minimax filtering, just projects the state on the
            // reduced space.
            model_state state(model_.GetNstate());
            model_state reduced_state(Nprojection_);
            model_.GetState(state);
            MltAdd(T(1), projection_, state, T(0), reduced_state);
            MltAdd(T(1), SeldonTrans, projection_, reduced_state,
                   T(0), state);
            model_.SetState(state);

            output_saver_.Save(projection_,  model_.GetTime(), "projection");

            MessageHandler::Send(*this, "model", "state_forecast");
            MessageHandler::Send(*this, "driver", "state_forecast");
            MessageHandler::Send(*this, "all", "::FilterInitialization end");

            return;
        }

        // Model state.
        double time = model_.GetTime();
        model_state full_state;
        model_.GetFullState(full_state);

        // Observation operator and data.
        observation_tangent_linear_operator H_tilde;
        observation_tangent_linear_operator H;
        Vector<T> y;

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            observation_manager_.GetObservation(y);
            Nobservation_ = y.GetLength();
            for (int i = 0; i < Nobservation_; i++)
                y(i) -= eta_;
            H_tilde = observation_manager_.GetTangentLinearOperator();
            H.Reallocate(Nobservation_, Nprojection_);
            MltAdd(T(1), SeldonNoTrans, H_tilde, SeldonTrans, projection_,
                   T(0), H);
            R_inv_ = observation_manager_.GetErrorVariance();
            GetInverse(R_inv_);
            Mlt(1. / (bound_over_standard_deviation_
                      * bound_over_standard_deviation_), R_inv_);
        }
        else
        {
            // If there are no observations, we add a fictitious observation
            // and a null H.
            y.Reallocate(1);
            y.Zero();
            Nobservation_ = 1;
            H.Reallocate(1, Nprojection_);
            H.Zero();
            R_inv_.Reallocate(1, 1);
            R_inv_(0, 0) = T(1);
        }

        output_saver_.Save(y, model_.GetTime(), "observation");
        output_saver_.Save(H, model_.GetTime(),
                           "reduced_tangent_linear_observation_operator");

        // Temporary variables.
        Matrix<T, General, RowMajor> mtmp, mtmp_1;
        Vector<T> vtmp;

        /*** Computes the minimax gain 'G_' ***/

        // Computes $\widecheck F_0$.
        Matrix<T, General, RowMajor> F_check(Nstate_, Nprojection_);
        vtmp.Reallocate(Nprojection_);
        for (int i = 0; i < Nstate_; i++)
        {
            GetCol(projection_, i, vtmp);
            Mlt(D_tilde_inv_(i), vtmp);
            SetRow(vtmp, i, F_check);
        }

        // Computes $\widecheck Q^{\frac 12}$.
        model_state_error_variance Q_sqrt_check;
        if (is_model_error_variance_diagonal_)
        {
            Q_sqrt_check.Reallocate(Nstate_, 1);
            Q_sqrt_check.Zero();
        }
        else
            Q_sqrt_check = model_.GetStateErrorVarianceSqrt();
        Nmode_Q_= Q_sqrt_check.GetN();
        for (int i = 0; i < Nstate_; i++)
            for (int j = 0; j < Nmode_Q_; j++)
                Q_sqrt_check(i, j) *= bound_over_standard_deviation_
                    * D_tilde_inv_(i);

        // Puts $\widecheck F_0^T \widecheck F_0$ into 'G_'.
        G_.Reallocate(Nprojection_, Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, SeldonNoTrans, F_check, T(0), G_);

        // Computes $\widecheck F_0^T \widecheck Q_^{\frac 12}$.
        Matrix<T, General, RowMajor> FtQ(Nprojection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, F_check,
               SeldonNoTrans, Q_sqrt_check, T(0), FtQ);

        // Computes $(I_{q\times q} + \widecheck Q^{\frac T2} \widecheck
        // Q^{\frac 12})^{-1}$.
        Matrix<T, General, RowMajor> IQtQinv(Nmode_Q_, Nmode_Q_);
        IQtQinv.SetIdentity();
        MltAdd(T(1), SeldonTrans, Q_sqrt_check,
               SeldonNoTrans, Q_sqrt_check, T(1), IQtQinv);
        // Note that IQtQinv is used for estimate's propagation too.
        GetInverse(IQtQinv);

        // Subtracts to 'G_' this part: $\widecheck F_{0}^T\widecheck Q^{\frac
        // 12} (I_{q\times q} + \widecheck Q^{\frac T2} \widecheck Q^{\frac
        // 12})^{-1} (\widecheck F_{0}^T \widecheck Q^{\frac 12})^T$.
        Matrix<T, General, RowMajor> FtQ_IQtQinv(Nprojection_, Nmode_Q_);
        MltAdd(T(1), FtQ, IQtQinv, T(0), FtQ_IQtQinv);
        MltAdd(T(-1), SeldonNoTrans, FtQ_IQtQinv, SeldonTrans, FtQ, T(1), G_);

        // Adds $H_0^T R_0^{-1} H_0$ to 'G_'.
        Matrix<T, General, RowMajor> Ht_Rinv(Nprojection_, Nobservation_);
        MltAdd(T(1), SeldonTrans, H, SeldonNoTrans, R_inv_, T(0), Ht_Rinv);
        output_saver_.Save(G_, model_.GetTime(), "minimax_gain_model_part");
        MltAdd(T(1), Ht_Rinv, H, T(1), G_);

        /*** Computes the minimax estimator ***/

        Vector<T> e0;
        model_.GetState(e0);
        Add(T(1), e0_, e0);
        for (int j = 0; j < Nstate_; j++)
            e0(j) *= D_tilde_inv_(j);

        MessageHandler::Send(*this, "model", "state_forecast");
        MessageHandler::Send(*this, "driver", "state_forecast");

        // Computes $F_0^T \widetilde D^{-\frac 12}\overline e$.
        state_.Reallocate(Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, e0, T(0), state_);

        // Computes $F_0^T \widecheck Q^{\frac 12} (I_{q \times q} +
        // \widecheck Q^{\frac T2} \widecheck Q^{\frac 12})^{-1} \widecheck
        // Q^{\frac T2} \widetilde D^{-\frac 12 } \overline e$.
        vtmp.Reallocate(Nmode_Q_);
        MltAdd(T(1), SeldonTrans, Q_sqrt_check, e0, T(0), vtmp);
        MltAdd(T(-1), FtQ_IQtQinv, vtmp, T(1), state_);

        // Adds $H^T_0 R_0^{-1} (y_0 - \eta_0)$.
        MltAdd(T(1), Ht_Rinv, y, T(1), state_);

        // Finally computes the minimax state, with $G_0^{-1} z_0$.
        mtmp = G_;
        GetAndSolveLU(mtmp, state_);

        vtmp.Reallocate(Nstate_);
        MltAdd(T(1), SeldonTrans, projection_, state_, T(0), vtmp);
        model_.SetState(vtmp);

        output_saver_.Save(G_, model_.GetTime(), "minimax_gain");
        output_saver_.Save(projection_,  model_.GetTime(), "projection");

        MessageHandler::Send(*this, "all", "::FilterInitialization end");
    }


    //! Propagates the state and the minimax gain.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Propagation()
    {
        MessageHandler::Send(*this, "all", "::Propagation begin");

        if (!with_filtering_)
        {
            // Without minimax filtering, just applies the model and projects
            // the state on the reduced space.
            model_.Forward();

            model_state state(model_.GetNstate());
            model_state reduced_state(Nprojection_);
            model_.GetState(state);
            MltAdd(T(1), projection_, state, T(0), reduced_state);
            MltAdd(T(1), SeldonTrans, projection_, reduced_state,
                   T(0), state);
            model_.SetState(state);

            output_saver_.Save(projection_, model_.GetTime(), "projection");

            inner_iteration_++;
            if (reduction_method_ != "none"
                && inner_iteration_ == Nstep_snapshot_ - 1)
            {
                inner_iteration_ = 0;
                mode_ = 0;
            }

            MessageHandler::Send(*this, "model", "state_forecast");
            MessageHandler::Send(*this, "driver", "state_forecast");
            MessageHandler::Send(*this, "all", "::Propagation end");

            return;
        }

        // Model state.
        double time = model_.GetTime();
        model_state full_state;
        model_.GetFullState(full_state);

        // Tangent linear model.
        Matrix<T, General, RowMajor> M(Nstate_, Nprevious_projection_);
        ComputeTangentLinearModel(M);

        output_saver_.Save(M, time, "reduced_tangent_linear_model");

        // Temporary variables.
        Vector<T> vtmp, vtmp_1;
        Matrix<T, General, RowMajor> mtmp, mtmp_1;

        /*** Model-related variables ***/

        // Computes $\widecheck M_t$.
        Matrix<T, General, RowMajor> M_check = M;
        for (int i = 0; i < Nstate_; i++)
            for (int j = 0; j < Nprevious_projection_; j++)
                M_check(i, j) *= D_tilde_inv_(i);

        // Model error.
        model_error_variance Q_sqrt;
        if (is_model_error_variance_diagonal_)
        {
            Q_sqrt.Reallocate(Nstate_, 1);
            Q_sqrt.Zero();
        }
        else
        {
            Q_sqrt = model_.GetErrorVarianceSqrt();
            output_saver_.Save(Q_sqrt, time, "model_error_variance_sqrt");
        }
        Nmode_Q_= Q_sqrt.GetN();

        // Forecast step.
        model_state forecast_state;
        model_.GetState(forecast_state);
        model_.ApplyOperator(forecast_state, true, false);

        MessageHandler::Send(*this, "model", "state_forecast");
        MessageHandler::Send(*this, "driver", "state_forecast");

        // Computes $\widecheck Q_t^{\frac 12}$.
        Matrix<T, General, RowMajor> Q_sqrt_check = Q_sqrt;
        for (int i = 0; i < Nstate_; i++)
            for (int j = 0; j < Nmode_Q_; j++)
                Q_sqrt_check(i, j) *= bound_over_standard_deviation_
                    * D_tilde_inv_(i);
        Q_sqrt.Clear();

        // Computes $\widecheck F_{t+1}$.
        Matrix<T, General, RowMajor> F_check(Nstate_, Nprojection_);
        vtmp.Reallocate(Nprojection_);
        for (int i = 0; i < Nstate_; i++)
        {
            GetCol(projection_, i, vtmp);
            Mlt(D_tilde_inv_(i), vtmp);
            SetRow(vtmp, i, F_check);
        }

        // Computes $\widecheck F_{t+1}^T \widecheck Q_t^{\frac 12}$.
        Matrix<T, General, RowMajor> FtQ(Nprojection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, F_check,
               SeldonNoTrans, Q_sqrt_check, T(0), FtQ);

        // Computes $G_t + \widecheck M_t^T \widecheck M_t$.
        Matrix<T, General, RowMajor> G_MtM(Nprevious_projection_,
                                           Nprevious_projection_);
        MltAdd(T(1), SeldonTrans, M_check,
               SeldonNoTrans, M_check, T(0), G_MtM);
        Add(T(1), G_, G_MtM);

        // Computes $\widecheck U_t$.
        Matrix<T, General, RowMajor> U_check(Nstate_, Nprevious_projection_);
        mtmp = G_MtM;
        GetInverse(mtmp);
        GetCholesky(mtmp);
        MltAdd(T(1), M_check, mtmp, T(0), U_check);

        // Computes $\widecheck U_t^T \widecheck Q_t^{\frac 12}$.
        Matrix<T, General, RowMajor> UtQ(Nprevious_projection_, Nmode_Q_);
        MltAdd(T(1), SeldonTrans, U_check,
               SeldonNoTrans, Q_sqrt_check, T(0), UtQ);

        // Computes $V_t$.
        Matrix<T, General, RowMajor> IQtQinv;
        Matrix<T, General, RowMajor> Vinv(Nmode_Q_, Nmode_Q_);
        Vinv.SetIdentity();
        MltAdd(T(1), SeldonTrans, Q_sqrt_check,
               SeldonNoTrans, Q_sqrt_check, T(1), Vinv);
        IQtQinv = Vinv;
        // Note that IQtQinv is used for estimate's propagation too.
        GetInverse(IQtQinv);

        mtmp.Reallocate(Nmode_Q_, Nprevious_projection_);
        MltAdd(T(1), SeldonTrans,
               Q_sqrt_check, SeldonNoTrans, U_check, T(0), mtmp);
        MltAdd(T(-1), SeldonNoTrans, mtmp, SeldonTrans, mtmp, T(1), Vinv);
        GetInverse(Vinv);

        /*** Observation-related variables ***/

        // Observation operator and data.
        observation_tangent_linear_operator H_tilde;
        observation_tangent_linear_operator H;
        Vector<T> y;

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            observation_manager_.GetObservation(y);
            Nobservation_ = y.GetLength();
            for (int i = 0; i < Nobservation_; i++)
                y(i) -= eta_;
            H_tilde = observation_manager_.GetTangentLinearOperator();
            H.Reallocate(Nobservation_, Nprojection_);
            MltAdd(T(1), SeldonNoTrans, H_tilde, SeldonTrans, projection_,
                   T(0), H);
            R_inv_ = observation_manager_.GetErrorVariance();
            GetInverse(R_inv_);
            Mlt(1. / (bound_over_standard_deviation_
                      * bound_over_standard_deviation_), R_inv_);
        }
        else
        {
            // If there are no observations, we add a fictitious observation
            // and a null H.
            y.Reallocate(1);
            y.Zero();
            Nobservation_ = 1;
            H.Reallocate(1, Nprojection_);
            H.Zero();
            R_inv_.Reallocate(1, 1);
            R_inv_(0, 0) = T(1);
        }

        output_saver_.Save(y, model_.GetTime(), "observation");
        output_saver_.Save(H, model_.GetTime(),
                           "reduced_tangent_linear_observation_operator");

        // Computes $H_{t+1}^T R_{t+1}^{-1} H_{t+1}$.
        Matrix<T, General, RowMajor> Ht_Rinv(Nprojection_, Nobservation_);
        MltAdd(T(1), SeldonTrans, H, SeldonNoTrans, R_inv_, T(0), Ht_Rinv);

        /*** Computes the minimax gain $G_t$ ***/

        // Computes $G_t$.
        G_.Reallocate(Nprojection_, Nprojection_);
        MltAdd(T(1), SeldonTrans, F_check, SeldonNoTrans, F_check, T(0), G_);

        Matrix<T, General, RowMajor> FtU(Nprojection_, Nprevious_projection_);
        MltAdd(T(1), SeldonTrans, F_check, SeldonNoTrans, U_check, T(0), FtU);
        MltAdd(T(-1), SeldonNoTrans, FtU, SeldonTrans, FtU, T(1), G_);

        mtmp_1 = FtQ;
        MltAdd(T(-1), FtU, UtQ, T(1), mtmp_1);
        Matrix<T, General, RowMajor> FtQ_FtUUtQ_Vinv(Nprojection_, Nmode_Q_);
        MltAdd(T(1), mtmp_1, Vinv, T(0), FtQ_FtUUtQ_Vinv);
        MltAdd(T(-1), SeldonNoTrans, FtQ_FtUUtQ_Vinv, SeldonTrans, mtmp_1,
               T(1), G_);

        output_saver_.Save(G_, model_.GetTime(), "minimax_gain_model_part");
        MltAdd(T(1), Ht_Rinv, H, T(1), G_);

        /*** Computes minimax state estimation ***/

        // Computes $\widecheck M_t \widehat x_t$.
        for (int i = 0; i < Nstate_; i++)
            forecast_state(i) *= D_tilde_inv_(i);

        // $\widecheck F_{t+1}^T (\widecheck M_t \widehat x_t)$.
        state_.Reallocate(Nprojection_);
        Mlt(SeldonTrans, F_check, forecast_state, state_);

        // $(\widecheck F_{t+1}^T \widecheck U_t) \widecheck U_t^T
        //  (\widecheck M_t \widehat x_t)$.
        vtmp.Reallocate(Nprevious_projection_);
        MltAdd(T(1), SeldonTrans, U_check, forecast_state, T(0), vtmp);
        MltAdd(T(-1), FtU, vtmp, T(1), state_);

        // $\left[ \widecheck Q_t^{\frac T2} - (\widecheck Q_t^{\frac T2}
        //  \widecheck U_t) \widecheck U_t^T \right]
        //  (\widecheck M_t \widehat x_t)$.
        Matrix<T, General, RowMajor> Q_UUtQ = Q_sqrt_check;
        MltAdd(T(-1), U_check, UtQ, T(1), Q_UUtQ);
        vtmp.Reallocate(Nmode_Q_);
        Mlt(SeldonTrans, Q_UUtQ, forecast_state, vtmp);
        MltAdd(T(-1), FtQ_FtUUtQ_Vinv, vtmp, T(1), state_);

        // Adds $H^T_{t+1} R_{t+1}^{-1} (y_{t+1} - \eta_{t+1})$.
        MltAdd(T(1), Ht_Rinv, y, T(1), state_);

        // Now multiplies everything by $G_{t+1}^{-1}$.
        mtmp = G_;
        GetAndSolveLU(mtmp, state_);

        vtmp.Reallocate(Nstate_);
        MltAdd(T(1), SeldonTrans, projection_, state_, T(0), vtmp);
        model_.SetState(vtmp);

        if (inner_iteration_ == 0)
        {
            previous_projection_ = projection_;
            Nprevious_projection_ = previous_projection_.GetM();
        }

        inner_iteration_++;
        if (reduction_method_ != "none"
            && inner_iteration_ == Nstep_snapshot_ - 1)
        {
            inner_iteration_ = 0;
            mode_ = 0;
        }

        output_saver_.Save(G_, model_.GetTime(), "minimax_gain");
        output_saver_.Save(projection_, model_.GetTime(), "projection");

        MessageHandler::Send(*this, "all", "::Propagation end");
    }


    //! Computes the full matrix of the tangent linear model.
    /*!
      \param[out] M the tangent linear model
    */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>
    ::ComputeTangentLinearModel(Matrix<T, General, RowMajor>& M)
    {
        Matrix<T, General, RowMajor> I(Nprevious_projection_,
                                       Nprevious_projection_);
        Vector<T> pi(Nstate_);
        Vector<T> Ei(Nprevious_projection_);

        for (int i = 0; i < Nprevious_projection_; i++)
        {
            Ei.Zero();
            Ei(i) = T(1);
            // Takes the i-th vector from the basis, forming the projected
            // subspace these vectors are stored as rows of the matrix
            // 'projection_'.
            MltAdd(T(1), SeldonTrans, previous_projection_, Ei, T(0), pi);
            // Applies tangent linear model to the i-th vector from the basis.
            model_.ApplyTangentLinearOperator(pi);
            SetCol(pi, i, M);
        }
    }


    //! Performs a step forward with snapshot recording.
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::SnapshotRecording()
    {

        /*** Model ***/

        model_.Forward();
        model_state state;
        model_.GetState(state);

        // Observation operator and data.
        observation_tangent_linear_operator H_tilde;
        Vector<T> y;

        observation_manager_.SetTime(model_, model_.GetTime());

        if (with_Hty_HtHx_ && observation_manager_.HasObservation())
        {
            observation_manager_.GetObservation(y);
            Nobservation_ = y.GetLength();
            for (int i = 0; i < Nobservation_; i++)
                y(i) -= eta_;
            H_tilde = observation_manager_.GetTangentLinearOperator();

            model_state Hty(Nstate_);
            MltAdd(T(1), SeldonTrans, H_tilde, y, T(0), Hty);

            if (Nsnapshot_ + 3 > snapshot_.GetN())
                snapshot_.Resize(Nstate_, Nsnapshot_ + 3);

            SetCol(Hty, Nsnapshot_++, snapshot_);

            Mlt(H_tilde, state, y);
            MltAdd(T(1), SeldonTrans, H_tilde, y, T(0), Hty);
            SetCol(Hty, Nsnapshot_++, snapshot_);
        }
        else if (Nsnapshot_ + 1 > snapshot_.GetN())
            snapshot_.Resize(Nstate_, Nsnapshot_ + 1);

        SetCol(state, Nsnapshot_++, snapshot_);

        inner_iteration_++;

        /*** POD ***/

        if (inner_iteration_ == Nstep_snapshot_)
        {
            snapshot_.Resize(snapshot_.GetM(), Nsnapshot_);

            output_saver_.Save(snapshot_, "snapshot");

            // SVD decomposition.
            int m = snapshot_.GetM();
            int n = snapshot_.GetN();
            char jobl('N'), jobr('S');
            int lwork = max(3 * min(m, n) + max(m, n), 5 * min(m, n));
            Vector<double> work(lwork);
            singular_value_.Reallocate(min(m, n));
            singular_value_.Zero();
            left_singular_vector_.Reallocate(m, n);
            left_singular_vector_.Zero();
            right_singular_vector_.Reallocate(0, 0);
            int one = 1;
            dgesvd_(&jobl, &jobr, &n, &m, snapshot_.GetData(),
                    &n, singular_value_.GetData(),
                    right_singular_vector_.GetData(), &one,
                    left_singular_vector_.GetData(), &n, work.GetData(),
                    &lwork, &lapack_info.GetInfoRef());

            output_saver_.Save(singular_value_, "singular_value");
            output_saver_.Save(left_singular_vector_, "left_singular_vector");
            output_saver_.Save(right_singular_vector_,
                               "right_singular_vector");

            // Projection matrix.
            T total_energy = Norm1(singular_value_);
            T energy = 0;
            // Collects the main modes for the projection.
            int i;
            for (i = 0;
                 i < min(Nprojection_max_, min(snapshot_.GetM(), Nsnapshot_))
                     && energy <= (1. - acceptable_error_) * total_energy;
                 i++)
                energy += singular_value_(i);
            Nprojection_ = i;
            // Saving the previous projection because it is still useful for
            // the first step of reduced minimax.
            if (!first_sequence_)
            {
                previous_projection_ = projection_;
                Nprevious_projection_ = previous_projection_.GetM();
            }
            projection_.Reallocate(Nprojection_, Nstate_);
            for (i = 0; i < Nprojection_; i++)
                for (int j = 0; j < Nstate_; j++)
                    projection_(i, j) = left_singular_vector_(j, i);
            if (first_sequence_)
            {
                previous_projection_ = projection_;
                Nprevious_projection_ = previous_projection_.GetM();
            }

            Logger::Log(GetName(),
                        "Total energy and percentage retained with "
                        + to_str(Nprojection_) + " modes: "
                        + to_str(energy) + ", "
                        + to_str(energy / total_energy));

            inner_iteration_ = 0;
            mode_ = 1;
            iteration_ -= Nstep_snapshot_ - 1;

            if (reduction_method_ != "none")
            {
                model_.SetFullState(full_state_);
                model_.SetTime(starting_time_);
            }
        }
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool ReducedMinimax<T, Model, ObservationManager>::HasFinished()
    {
        // The condition 'mode_ == 0 && inner_iteration_ != Nsnapshot_' means
        // that a sequence of error computation has just ended. The simulation
        // should therefore end with a complete sequence of error computation.
        return (reduction_method_ == "none"
                || (mode_ == 0 && inner_iteration_ == 0))
            && model_.HasFinished();
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the current mode (snapshot recording or minimax calculations).
    /*!
      \return The mode: 0 for snapshot recording and 1 for minimax
      calculations.
    */
    template <class T, class Model, class ObservationManager>
    int ReducedMinimax<T, Model, ObservationManager>::GetMode() const
    {
        return mode_;
    }


    //! Returns the current projection matrix.
    /*!
      \return The current projection matrix.
    */
    template <class T, class Model, class ObservationManager>
    Matrix<T, General, RowMajor>&
    ReducedMinimax<T, Model, ObservationManager>::GetProjection()
    {
        return projection_;
    }


    //! Returns the previous projection matrix.
    /*!
      \return The previous projection matrix.
    */
    template <class T, class Model, class ObservationManager>
    Matrix<T, General, RowMajor>&
    ReducedMinimax<T, Model, ObservationManager>::GetPreviousProjection()
    {
        return previous_projection_;
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class Model, class ObservationManager>
    Model& ReducedMinimax<T, Model, ObservationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class T, class Model, class ObservationManager>
    ObservationManager&
    ReducedMinimax<T, Model, ObservationManager>::GetObservationManager()
    {
        return observation_manager_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class T, class Model, class ObservationManager>
    OutputSaver&
    ReducedMinimax<T, Model, ObservationManager>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string ReducedMinimax<T, Model, ObservationManager>::GetName() const
    {
        return "ReducedMinimax";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void ReducedMinimax<T, Model, ObservationManager>::Message(string message)
    {
        model_state state;
        if (message.find("state_analysis") != string::npos
            && message.find("reduced_state_analysis") == string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, model_.GetTime(), "state_analysis");
        }
        else if (message.find("state_forecast") != string::npos)
        {
            model_.GetState(state);
            output_saver_.Save(state, model_.GetTime(), "state_forecast");
        }
        else if (message.find("reduced_state_analysis") != string::npos)
            output_saver_.Save(state_, model_.GetTime(),
                               "reduced_state_analysis");
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDMINIMAX_CXX
#endif
