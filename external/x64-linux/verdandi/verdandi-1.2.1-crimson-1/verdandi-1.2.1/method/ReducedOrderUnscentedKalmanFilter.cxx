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


#ifndef VERDANDI_FILE_METHOD_REDUCEDORDERUNSCENTEDKALMANFILTER_CXX

#include "ReducedOrderUnscentedKalmanFilter.hxx"

#include "seldon/vector/VectorCollection.cxx"

#include "SigmaPoint.cxx"

#include "seldon/computation/solver/SparseCholeskyFactorisation.cxx"


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
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::ReducedOrderUnscentedKalmanFilter()
    {

        /*** Initializations ***/

#if defined(VERDANDI_WITH_MPI)
        MPI::Init();
        rank_ = MPI::COMM_WORLD.Get_rank();
        Nprocess_ = MPI::COMM_WORLD.Get_size();
        if (rank_ == 0)
        {
#endif
            MessageHandler::AddRecipient("model", model_,
                                         Model::StaticMessage);
            MessageHandler::AddRecipient("observation_manager",
                                         observation_manager_,
                                         ObservationManager::StaticMessage);
            MessageHandler::AddRecipient("driver", *this,
                                         ReducedOrderUnscentedKalmanFilter
                                         ::StaticMessage);
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::~ReducedOrderUnscentedKalmanFilter()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration, initialize_model,
                   initialize_observation_manager);
    }

    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Initialize begin");


        /***************************
         * Reads the configuration *
         ***************************/


        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("reduced_order_unscented_kalman_filter.");

        /*** Model ***/

        configuration.Set("model.configuration_file", "", configuration_file_,
                          model_configuration_file_);

        /*** Observation manager ***/

        configuration.Set("observation_manager.configuration_file", "",
                          configuration_file_,
                          observation_configuration_file_);

        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.show_iteration",
                          option_display_["show_iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("display.show_time", option_display_["show_time"]);

        /*** Assimilation options ***/

        configuration.Set("data_assimilation.analyze_first_step",
                          analyze_first_step_);
        configuration.Set("data_assimilation.with_resampling",
                          with_resampling_);
        configuration.Set("data_assimilation.observation_error_variance",
                          "ops_in(v, {'matrix', 'matrix_inverse'})",
                          observation_error_variance_);
        configuration.Set("data_assimilation.with_innovation_computation",
                          innovation_computation_);

        /*** Sigma-points ***/

        configuration.Set("sigma_point.type",
                          "ops_in(v, {'canonical', 'star', 'simplex'})",
                          sigma_point_type_);

#if defined(VERDANDI_WITH_MPI)

        /*** MPI ***/

        configuration.Set("mpi.algorithm", "ops_in(v, {0, 1, 2})",
                          algorithm_);
        configuration.Set("mpi.master_process_contribution",
                          master_process_contribution_);
        if (master_process_contribution_ < 0. ||
            master_process_contribution_ > 1.)
            throw "Contribution of process 0 should be in [0, 1] "
                "], but " + to_str(master_process_contribution_) +
                " was provided.";

#endif

        /*** Ouput saver ***/

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            configuration.
                SetPrefix("reduced_order_unscented_kalman_filter"
                          ".output_saver.");
            output_saver_.Initialize(configuration);
            output_saver_.Empty("state_forecast");
            output_saver_.Empty("state_analysis");

            /*** Logger and read configuration ***/

            configuration.SetPrefix("reduced_order_unscented_kalman_filter.");

            if (configuration.Exists("output.configuration"))
            {
                string output_configuration;
                configuration.Set("output.configuration",
                                  output_configuration);
                configuration.WriteLuaDefinition(output_configuration);
            }
#if defined(VERDANDI_WITH_MPI)
        }
#endif

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        /*** Initializations ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        if (initialize_observation_manager)
        {
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);
            if (innovation_computation_)
                observation_manager_.DiscardObservation(false);
        }
        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

        sigma_point_matrix U;
        model_.GetStateErrorVarianceSqrt(L_, U);
        Copy(U, U_);
        U_inv_.Copy(U_);

        GetInverse(U_inv_);

        Nreduced_ = U_.GetN();

        /*** Sigma-points ***/

        sigma_point_matrix V_trans;
        if (sigma_point_type_ == "canonical")
            ComputeCanonicalSigmaPoint(Nreduced_, V_trans, D_alpha_,
                                       alpha_constant_);
        else if (sigma_point_type_ == "star")
            ComputeStarSigmaPoint(Nreduced_, V_trans, D_alpha_,
                                  alpha_constant_);
        else if (sigma_point_type_ == "simplex")
            ComputeSimplexSigmaPoint(Nreduced_, V_trans, D_alpha_,
                                     alpha_constant_);
        if (alpha_constant_)
            alpha_ = D_alpha_(0);

        Nsigma_point_ = V_trans.GetM();

        // Initializes transpose of I.
        sigma_point_matrix P_alpha_v(Nreduced_, Nreduced_);
        I_trans_.Reallocate(Nsigma_point_, Nreduced_);

        if (alpha_constant_)
        {
            MltAdd(T(alpha_), SeldonTrans, V_trans, SeldonNoTrans, V_trans,
                   T(0), P_alpha_v);
            GetInverse(P_alpha_v);
            GetCholesky(P_alpha_v);
            MltAdd(T(1), SeldonNoTrans, V_trans, SeldonTrans, P_alpha_v,
                   T(0), I_trans_);
        }
        else
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Initialize()", "Calculation not "
                                 "implemented for no constant alpha_i.");
        I_.Copy(I_trans_);
        Transpose(I_);

#if defined(VERDANDI_WITH_PETSC)
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank_);
#endif

        // Initializes D_v.
        D_v_.Reallocate(Nsigma_point_, Nsigma_point_);
        if (alpha_constant_)
            MltAdd(T(alpha_ * alpha_), SeldonNoTrans, I_trans_, SeldonTrans,
                   I_trans_, T(0), D_v_);
        else
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Initialize()", "Calculation not "
                                 "implemented for no constant alpha_i.");

#if defined(VERDANDI_WITH_MPI)

        /*** Local sigma-points ***/

        Nlocal_sigma_point_ = int(Nsigma_point_ / Nprocess_);
        int r = Nsigma_point_ % Nprocess_;
        if (Nprocess_ - rank_ - 1 < r)
            Nlocal_sigma_point_++;

        Nlocal_sigma_point_sum_.Reallocate(Nprocess_ + 1);
        Nlocal_sigma_point_sum_(0) = 0;
        for (int i = 0; i < Nprocess_; i++)
        {
            Nlocal_sigma_point_sum_(i + 1) = int(Nsigma_point_ / Nprocess_);
            if (Nprocess_ - i - 1 < r)
                Nlocal_sigma_point_sum_(i + 1) += 1;
            Nlocal_sigma_point_sum_(i + 1) += Nlocal_sigma_point_sum_(i);
        }

        int Nsigma_point_0;
        Nsigma_point_0 =  int(master_process_contribution_ *
                              Nlocal_sigma_point_sum_(1));
        int q;
        q = int((Nlocal_sigma_point_sum_(1) - Nsigma_point_0) /
                (Nprocess_ - 1));
        r = (Nlocal_sigma_point_sum_(1) - Nsigma_point_0) % (Nprocess_ - 1);
        for (int i = 1; i < Nprocess_; i++)
        {
            Nlocal_sigma_point_sum_(i + 1) -= Nlocal_sigma_point_sum_(1) -
                Nsigma_point_0;
            Nlocal_sigma_point_sum_(i + 1) += q * i;
            if (i <= r)
                Nlocal_sigma_point_sum_(i + 1) += i;
            else
                Nlocal_sigma_point_sum_(i + 1) += r;

        }
        Nlocal_sigma_point_sum_(1) = Nsigma_point_0;
        for (int i = Nlocal_sigma_point_sum_(rank_);
             i < Nlocal_sigma_point_sum_(rank_ + 1); i++)
            local_sigma_point_.PushBack(i);
        Nlocal_sigma_point_ = local_sigma_point_.GetM();

        I_trans_global_.Copy(I_trans_);
        I_trans_.Reallocate(Nlocal_sigma_point_, Nreduced_);
        sigma_point I_col;
        for (int i = 0; i < Nlocal_sigma_point_; i++)
        {
            GetRow(I_trans_global_, local_sigma_point_(i), I_col);
            SetRow(I_col, i, I_trans_);
        }

        /*** Local columns of covariance matrix ***/

        Nlocal_filtered_column_ = int(Nreduced_ / Nprocess_);
        r = Nreduced_ % Nprocess_;
        if (Nprocess_ - rank_ - 1 < r)
            Nlocal_filtered_column_++;

        Nlocal_filtered_column_sum_.Reallocate(Nprocess_ + 1);
        Nlocal_filtered_column_sum_(0) = 0;
        for (int i = 0; i < Nprocess_; i++)
        {
            Nlocal_filtered_column_sum_(i + 1) = int(Nreduced_ / Nprocess_);
            if (Nprocess_ - i - 1 < r)
                Nlocal_filtered_column_sum_(i + 1) += 1;
            Nlocal_filtered_column_sum_(i + 1)
                += Nlocal_filtered_column_sum_(i);
        }

        for (int i = Nlocal_filtered_column_sum_(rank_);
             i < Nlocal_filtered_column_sum_(rank_ + 1); i++)
            local_filtered_column_.PushBack(i);

#endif

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            if (initialize_model)
            {
                MessageHandler::Send(*this, "model", "initial condition");
                MessageHandler::Send(*this, "driver", "initial condition");
            }

            MessageHandler::Send(*this, "all", "::Initialize end");

#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Initializes a step for the unscented Kalman filter.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::InitializeStep()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::InitializeStep begin");

        model_.InitializeStep();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Forward()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Forward begin");

        if (sigma_point_type_ == "simplex")
        {
#if defined(VERDANDI_WITH_MPI)
            model_state& x = model_.GetState();

            if (algorithm_ == 0)
            {
                sigma_point_matrix X_i_trans_global;
                sigma_point x_col;
                X_i_trans_.Reallocate(Nlocal_sigma_point_, Nstate_);
                if (rank_ == 0)
                {
                    sigma_point_matrix tmp;
                    GetCholesky(U_inv_);
                    Copy(L_, tmp);
                    MltAdd(T(1), tmp, U_inv_, T(0), L_);

                    // Computes X_n^{(i)+}.
                    X_i_trans_global.Reallocate(Nsigma_point_, Nstate_);
                    for (int i = 0; i < Nsigma_point_; i++)
                        SetRow(x, i, X_i_trans_global);
                    MltAdd(T(1), SeldonNoTrans, I_trans_global_, SeldonTrans,
                           L_, T(1), X_i_trans_global);
                }

                int displacement[Nprocess_],  sendcount[Nprocess_];
                for (int i = 0; i < Nprocess_; i++)
                {
                    sendcount[i] = (Nlocal_sigma_point_sum_(i + 1)
                                    - Nlocal_sigma_point_sum_(i))
                        * Nstate_;
                    displacement[i] = Nlocal_sigma_point_sum_(i)
                        * Nstate_;
                }

                MPI::COMM_WORLD.Scatterv(X_i_trans_global.GetData(),
                                         sendcount, displacement, MPI::DOUBLE,
                                         X_i_trans_.GetData(),
                                         Nlocal_sigma_point_ * Nstate_,
                                         MPI::DOUBLE, 0);

                /*** Prediction ***/

                // Computes X_{n + 1}^-.
                x.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRow(X_i_trans_, i, x_col);
                    model_.ApplyOperator(x_col, i + 1 == Nlocal_sigma_point_,
                                         true);
                    Add(T(alpha_), x_col, x);
                    if (rank_ == 0)
                        SetRow(x_col, i, X_i_trans_global);
                    SetRow(x_col, i, X_i_trans_);
                }

                model_state working_vector;

                if (rank_ == 0)
                {
                    working_vector.Copy(x);
                    working_vector.Fill(T(0));
                }

                MPI::COMM_WORLD.Reduce(x.GetData(), working_vector.GetData(),
                                       x.GetM(), MPI::DOUBLE, MPI::SUM, 0);


                if (rank_ == 0)
                    Copy(working_vector, x);

                model_.StateUpdated();

                MPI::COMM_WORLD.
                    Gatherv(X_i_trans_.GetData(),
                            Nlocal_sigma_point_ * Nstate_, MPI::DOUBLE,
                            X_i_trans_global.GetData(),  sendcount,
                            displacement, MPI::DOUBLE, 0);

                if (rank_ == 0)
                    MltAdd(T(alpha_), SeldonTrans, X_i_trans_global,
                           SeldonNoTrans, I_trans_global_, T(0), L_);
            }
            else if (algorithm_ == 1 || algorithm_ == 2)
            {

                /*** Sampling ***/

                sigma_point_matrix tmp;
                GetCholesky(U_inv_);
                Copy(L_, tmp);
                MltAdd(T(1), tmp, U_inv_, T(0), L_);

                // Computes X_n^{(i)+}.
                X_i_trans_.Reallocate(Nlocal_sigma_point_, Nstate_);
                sigma_point x_col;
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                    SetRow(x, i, X_i_trans_);
                MltAdd(T(1), SeldonNoTrans, I_trans_, SeldonTrans, L_, T(1),
                       X_i_trans_);

                /*** Prediction ***/

                // Computes X_{n + 1}^-.
                x.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRow(X_i_trans_, i, x_col);
                    model_.ApplyOperator(x_col, i + 1 == Nlocal_sigma_point_,
                                         true);
                    Add(T(alpha_), x_col, x);
                    SetRow(x_col, i, X_i_trans_);
                }

                model_state working_vector;
                Copy(x, working_vector);
                MPI::COMM_WORLD.Allreduce(x.GetData(),
                                          working_vector.GetData(), x.GetM(),
                                          MPI::DOUBLE, MPI::SUM);

                int displacement[Nprocess_],  recvcount[Nprocess_];
                for (int i = 0; i < Nprocess_; i++)
                {
                    recvcount[i] = (Nlocal_sigma_point_sum_(i + 1)
                                    - Nlocal_sigma_point_sum_(i))
                        * Nstate_;
                    displacement[i] = Nlocal_sigma_point_sum_(i)
                        * Nstate_;
                }

                // Computes L_{n + 1}.
                if (algorithm_ == 1)
                {
                    sigma_point_matrix
                        X_i_trans_global(Nsigma_point_, Nstate_);

                    MPI::COMM_WORLD.
                        Allgatherv(X_i_trans_.GetData(),
                                   Nlocal_sigma_point_ * Nstate_, MPI::DOUBLE,
                                   X_i_trans_global.GetData(),  recvcount,
                                   displacement, MPI::DOUBLE);

                    MltAdd(T(alpha_), SeldonTrans, X_i_trans_global,
                           SeldonNoTrans, I_trans_global_, T(0), L_);
                }
                else
                {
                    sigma_point_matrix L_local(Nstate_, Nreduced_);
                    MltAdd(T(alpha_), SeldonTrans, X_i_trans_, SeldonNoTrans,
                           I_trans_, T(0), L_local);

                    L_.Fill(T(0));
                    MPI::COMM_WORLD.
                        Allreduce(L_local.GetData(), L_.GetData(),
                                  L_.GetSize(), MPI::DOUBLE, MPI::SUM);
                }

                model_.StateUpdated();
            }
#else
            model_state x;
            x.Copy(model_.GetState()); ///////////////////////////////////////////////////////////////////////////////////////////////////

            /*** Sampling ***/

            model_state_error_variance tmp;
            GetCholesky(U_inv_);

            Copy(L_, tmp);
            MltAdd(T(1), tmp, U_inv_, T(0), L_);

            // Computes X_n^{(i)+}.
            Reallocate(X_i_, Nstate_, Nsigma_point_, model_);
            for (int i = 0; i < Nsigma_point_; i++)
                SetCol(x, i, X_i_);

            MltAdd(T(1), L_, I_, T(1), X_i_);

            /*** Prediction ***/

            // Computes X_{n + 1}^-.
            x.Fill(T(0));
            model_state x_col;
            Reallocate(x_col, x.GetM(), model_);
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetCol(X_i_, i, x_col);
                model_.ApplyOperator(x_col, i + 1 == Nsigma_point_, false); ///////////////////////////////////////////////////////////////////////////////////////////////////
                Add(T(alpha_), x_col, x);
                SetCol(x_col, i, X_i_);
            }

            /*** Resampling ***/

            if (with_resampling_)
            {
#if defined(VERDANDI_WITH_PETSC)
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Forward()", "'resampling 'option "
                                     "not support yet.");
#else
                MltAdd(T(alpha_), X_i_, I_trans_, T(0), L_);
                for(int i = 0; i < Nsigma_point_; i++)
                    SetCol(x, i, X_i_);
                MltAdd(T(1), L_, I_, T(1), X_i_);
#endif
            }

            // Computes L_{n + 1}.
            MltAdd(T(alpha_), X_i_, I_trans_, T(0), L_);

            model_.GetState().Copy(x); ///////////////////////////////////////////////////////////////////////////////////////////////////
            model_.StateUpdated(); ///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
        }
        else
        {
#if defined(VERDANDI_WITH_MPI) || defined(VERDANDI_WITH_PETSC)
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Forward()", "Parallel algorithm not "
                                 "implemented yet for the 'no"
                                 " simplex' cases.");
#else
            model_state x =  model_.GetState();
            x.Copy(model_.GetState());

            /*** Sampling ***/

            sigma_point_matrix tmp;
            GetCholesky(U_inv_);
            Copy(L_, tmp);
            MltAdd(T(1), tmp, U_inv_, T(0), L_);

            // Computes X_n^{(i)+}.
            X_i_trans_.Reallocate(Nsigma_point_, Nstate_);
            sigma_point x_col;
            for (int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, X_i_trans_);

            MltAdd(T(1), SeldonNoTrans, I_trans_, SeldonTrans, L_, T(1),
                   X_i_trans_);

            /*** Prediction ***/

            // Computes X_{n + 1}^-.
            x.Fill(T(0));
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetRow(X_i_trans_, i, x_col);
                model_.ApplyOperator(x_col, i + 1 == Nsigma_point_, false);
                Add(T(alpha_), x_col, x);
                SetRow(x_col, i, X_i_trans_);
            }

            model_.GetState().Copy(x);
            model_.StateUpdated();

            /*** Resampling with SVD ***/

            sigma_point_matrix M_trans(Nsigma_point_, Nstate_);
            for (int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, M_trans);
            Mlt(T(-1), M_trans);
            Add(T(1), X_i_trans_, M_trans);

            if (alpha_constant_)
                Mlt(sqrt(alpha_), M_trans);
            else
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Forward()", "Calculation not "
                                     "implemented for no constant alpha_i.");

            sigma_point_matrix G(Nsigma_point_, Nsigma_point_);
            MltAdd(T(1), SeldonNoTrans, M_trans, SeldonTrans, M_trans,
                   T(0), G);

            Vector<T> lambda;
            Matrix<T> U, V;
            GetSVD(G, lambda, U, V);
            U.Resize(Nsigma_point_, Nreduced_);

            sigma_point_matrix
                working_matrix_rr(Nsigma_point_, Nsigma_point_),
                working_matrix_rN(X_i_trans_);

            MltAdd(T(sqrt(alpha_)), SeldonNoTrans, U, SeldonTrans, I_trans_,
                   T(0), working_matrix_rr);

            for(int i = 0; i < Nsigma_point_; i++)
                SetRow(x, i, X_i_trans_);
            Add(T(-1), X_i_trans_, working_matrix_rN);

            MltAdd(T(1), SeldonTrans, working_matrix_rr, SeldonNoTrans,
                   working_matrix_rN, T(1), X_i_trans_);

            // Computes L_{n + 1}.
            MltAdd(T(alpha_), SeldonTrans, X_i_trans_, SeldonNoTrans,
                   I_trans_, T(0), L_);
#endif
        }

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            MessageHandler::Send(*this, "model", "forecast");
            MessageHandler::Send(*this, "observation_manager", "forecast");
            MessageHandler::Send(*this, "driver", "forecast");

            MessageHandler::Send(*this, "all", "::Forward end");
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Analyze()
    {

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (!observation_manager_.HasObservation())
        {
#if defined(VERDANDI_WITH_MPI)
            if (rank_ == 0)
#endif
                MessageHandler::Send(*this, "all", "::Analyze end");
            return;
        }

        Nobservation_  = observation_manager_.GetNobservation();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            if (option_display_["show_time"] && rank_ == 0)
                cout << "Performing Reduced Order UKF at time step ["
                     << model_.GetTime() << "]..." << endl;
#if defined(VERDANDI_WITH_MPI)
        }
#endif

        if (sigma_point_type_ == "simplex")
        {
#if defined(VERDANDI_WITH_MPI)
            if (algorithm_ == 0)
            {
                // Computes [HX_{n+1}^{*}].
                sigma_point_matrix Z_i_trans,
                    Z_i_trans_local(Nlocal_sigma_point_, Nobservation_);
                sigma_point x_col;
                observation z_col(Nobservation_), z_local(Nobservation_), z;
                z_local.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    SetRow(z_col, i, Z_i_trans_local);
                    Add(T(alpha_), z_col, z_local);
                    x_col.Nullify();
                }

                if (rank_ == 0)
                {
                    z.Reallocate(Nobservation_);
                    z.Fill(T(0));
                    Z_i_trans.Reallocate(Nsigma_point_, Nobservation_);
                }

                MPI::COMM_WORLD.
                    Reduce(z_local.GetData(), z.GetData(),
                           z_local.GetM(), MPI::DOUBLE, MPI::SUM, 0);

                int displacement[Nprocess_],  sendcount[Nprocess_];
                for (int i = 0; i < Nprocess_; i++)
                {
                    sendcount[i] = (Nlocal_sigma_point_sum_(i + 1)
                                    - Nlocal_sigma_point_sum_(i))
                        * Nobservation_;
                    displacement[i] = Nlocal_sigma_point_sum_(i)
                        * Nobservation_;
                }

                MPI::COMM_WORLD.
                    Gatherv(Z_i_trans_local.GetData(),
                            Nlocal_sigma_point_ * Nobservation_, MPI::DOUBLE,
                            Z_i_trans.GetData(),  sendcount, displacement,
                            MPI::DOUBLE, 0);

                if (rank_ == 0)
                {
                    sigma_point_matrix HL_trans(Nreduced_, Nobservation_);
                    MltAdd(T(alpha_), SeldonTrans, I_trans_global_,
                           SeldonNoTrans, Z_i_trans, T(0), HL_trans);
                    sigma_point_matrix
                        working_matrix_po(Nreduced_, Nobservation_), tmp;

                    if (observation_error_variance_ == "matrix_inverse")
                        Mlt(HL_trans,
                            observation_manager_.GetErrorVarianceInverse(),
                            working_matrix_po);
                    else
                    {
                        observation_error_variance R_inv;
                        Copy(observation_manager_.GetErrorVariance(),
                             R_inv);
                        GetInverse(R_inv);
                        Mlt(HL_trans, R_inv, working_matrix_po);
                    }

                    U_inv_.SetIdentity();
                    MltAdd(T(1), SeldonNoTrans, working_matrix_po,
                           SeldonTrans, HL_trans, T(1), U_inv_);
                    GetInverse(U_inv_);

                    // Computes K.
                    sigma_point_matrix K(Nstate_, Nobservation_);
                    Copy(working_matrix_po, tmp);
                    MltAdd(T(1), U_inv_, working_matrix_po, T(0), tmp);
                    MltAdd(T(1), L_, tmp, T(0), K);

                    // Computes innovation.
                    observation innovation;
                    observation_manager_.GetObservation(innovation);
                    Add(T(-1), z, innovation);

                    // Updates.
                    model_state& x =  model_.GetState();
                    MltAdd(T(1), K, innovation, T(1), x);
                    model_.StateUpdated();
                }
            }
            else if (algorithm_ == 1 || algorithm_ == 2)
            {
                observation innovation(Nobservation_);
                MPI::Request send_request[Nprocess_ - 1], recv_request;
                if (rank_ == 0)
                {
                    observation_manager_.GetObservation(innovation);
                    for (int i = 1; i < Nprocess_; i++)
                        send_request[i] =
                            MPI::COMM_WORLD.
                            Isend(innovation.GetData(),
                                  Nobservation_, MPI::DOUBLE, i, 0);
                }
                else
                    recv_request =
                        MPI::COMM_WORLD.Irecv(innovation.GetData(),
                                              Nobservation_, MPI::DOUBLE, 0,
                                              MPI::ANY_TAG);
                // Computes [HX_{n+1}^{*}].
                sigma_point_matrix
                    Z_i_trans(Nlocal_sigma_point_, Nobservation_);
                sigma_point x_col;
                observation z_col, z(Nobservation_), z_local(Nobservation_);
                z_local.Fill(T(0));
                for (int i = 0; i < Nlocal_sigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    GetRowPointer(Z_i_trans, i, z_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    Add(T(alpha_), z_col, z_local);
                    x_col.Nullify();
                    z_col.Nullify();
                }

                sigma_point_matrix HL_local_trans(Nreduced_, Nobservation_),
                    HL_trans(Nreduced_, Nobservation_);
                MltAdd(T(alpha_), SeldonTrans, I_trans_, SeldonNoTrans,
                       Z_i_trans, T(0), HL_local_trans);

                MPI::COMM_WORLD.Allreduce(HL_local_trans.GetData(),
                                          HL_trans.GetData(),
                                          HL_trans.GetSize(), MPI::DOUBLE,
                                          MPI::SUM);

                sigma_point_matrix
                    working_matrix_po(Nreduced_, Nobservation_), tmp;
                sigma_point_matrix K;
                if (algorithm_ == 1)
                {
                    if (observation_error_variance_ == "matrix_inverse")
                        Mlt(HL_trans,
                            observation_manager_.GetErrorVarianceInverse(),
                            working_matrix_po);
                    else
                    {
                        observation_error_variance R_inv;
                        Copy(observation_manager_.GetErrorVariance(),
                             R_inv);
                        GetInverse(R_inv);
                        Mlt(HL_trans, R_inv, working_matrix_po);
                    }

                    U_inv_.SetIdentity();
                    MltAdd(T(1), SeldonNoTrans, working_matrix_po,
                           SeldonTrans, HL_trans, T(1), U_inv_);
                    GetInverse(U_inv_);

                    // Computes K.
                    K.Reallocate(Nstate_, Nobservation_);
                    Copy(working_matrix_po, tmp);
                    MltAdd(T(1), U_inv_, working_matrix_po, T(0), tmp);
                    MltAdd(T(1), L_, tmp, T(0), K);
                }
                else
                {
                    sigma_point_matrix working_matrix_po_local(
                        Nlocal_filtered_column_, Nobservation_);
                    sigma_point_matrix
                        U_local(Nlocal_filtered_column_, Nreduced_);

                    HL_local_trans.
                        Reallocate(Nlocal_filtered_column_, Nobservation_);
                    for (int i = 0; i < Nlocal_filtered_column_; i++)
                    {
                        GetRow(HL_trans, local_filtered_column_(i), z_col);
                        SetRow(z_col, i, HL_local_trans);
                    }

                    if (observation_error_variance_ == "matrix_inverse")
                        Mlt(HL_local_trans,
                            observation_manager_.GetErrorVarianceInverse(),
                            working_matrix_po_local);
                    else
                    {
                        observation_error_variance R_inv;
                        Copy(observation_manager_.GetErrorVariance(),
                             R_inv);
                        GetInverse(R_inv);
                        Mlt(HL_local_trans, R_inv, working_matrix_po_local);
                    }

                    U_local.Fill(T(0));
                    for (int i = 0; i < Nlocal_filtered_column_; i++)
                        U_local(i, local_filtered_column_(i)) = T(1);

                    MltAdd(T(1), SeldonNoTrans, working_matrix_po_local,
                           SeldonTrans, HL_trans, T(1), U_local);

                    int displacement[Nprocess_],  recvcount[Nprocess_];
                    for (int i = 0; i < Nprocess_; i++)
                    {
                        recvcount[i] = (Nlocal_filtered_column_sum_(i + 1)
                                        - Nlocal_filtered_column_sum_(i))
                            * Nreduced_;
                        displacement[i] = Nlocal_filtered_column_sum_(i)
                            * Nreduced_;
                    }

                    MPI::COMM_WORLD.
                        Allgatherv(U_local.GetData(),
                                   Nlocal_filtered_column_ * Nreduced_,
                                   MPI::DOUBLE, U_inv_.GetData(),  recvcount,
                                   displacement, MPI::DOUBLE);

                    for (int i = 0; i < Nprocess_; i++)
                    {
                        recvcount[i] = (Nlocal_filtered_column_sum_(i + 1)
                                        - Nlocal_filtered_column_sum_(i))
                            * Nobservation_;
                        displacement[i] = Nlocal_filtered_column_sum_(i)
                            * Nobservation_;
                    }

                    MPI::COMM_WORLD.
                        Allgatherv(working_matrix_po_local.GetData(),
                                   Nlocal_filtered_column_ * Nobservation_,
                                   MPI::DOUBLE, working_matrix_po.GetData(),
                                   recvcount, displacement, MPI::DOUBLE);

                    GetInverse(U_inv_);

                    // Computes K.
                    K.Reallocate(Nstate_, Nobservation_);
                    for (int i = 0; i < Nlocal_filtered_column_; i++)
                    {
                        GetRow(U_inv_, local_filtered_column_(i), z_col);
                        SetRow(z_col, i, U_local);
                    }

                    MltAdd(T(1), U_local, working_matrix_po, T(0),
                           working_matrix_po_local);

                    MPI::COMM_WORLD.
                        Allgatherv(working_matrix_po_local.GetData(),
                                   Nlocal_filtered_column_ * Nobservation_,
                                   MPI::DOUBLE, working_matrix_po.GetData(),
                                   recvcount, displacement, MPI::DOUBLE);

                    MltAdd(T(1), L_, working_matrix_po, T(0), K);
                }

                // Computes innovation.
                if (rank_ == 0)
                    MPI::Request::Waitall(Nprocess_ - 1, send_request);
                else
                    recv_request.Wait();

                MPI::COMM_WORLD.Allreduce(z_local.GetData(), z.GetData(),
                                          Nobservation_, MPI::DOUBLE,
                                          MPI::SUM);
                Add(T(-1), z, innovation);

                // Updates.
                model_state& x =  model_.GetState();
                MltAdd(T(1), K, innovation, T(1), x);
                model_.StateUpdated();
            }
#else
            // Computes [HX_{n+1}^{*}].
            sigma_point_matrix Z_i_trans(Nsigma_point_, Nobservation_);
            model_state x_col;
            Reallocate(x_col, Nstate_, model_);
            model_.StateUpdated();    ///////////////////////////////////////////////////////////////////////////////////////////////////

            observation z_col, z(Nobservation_);
            z.Fill(T(0));
            if (!innovation_computation_)
                for (int i = 0; i < Nsigma_point_; i++)
                {
                    GetCol(X_i_, i, x_col);
                    GetRowPointer(Z_i_trans, i, z_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    Add(T(alpha_), z_col, z);
                    z_col.Nullify();
                }
            else
                for (int i = 0; i < Nsigma_point_; i++)
                {
                    GetCol(X_i_, i, x_col);
                    observation_manager_.GetInnovation(x_col, z_col); ///////////////////////////////////////////////////////////////////////////////////////////////////
                    Add(T(alpha_), z_col, z);
                    SetRow(z_col, i, Z_i_trans);
                }

            sigma_point_matrix HL_trans(Nreduced_, Nobservation_);
            MltAdd(T(alpha_), SeldonTrans, I_trans_, SeldonNoTrans, Z_i_trans,
                   T(0), HL_trans);

            observation_error_variance R_inv;
            sigma_point_matrix working_matrix_po(Nreduced_, Nobservation_),
                tmp;

            if (observation_error_variance_ == "matrix_inverse")
                Mlt(HL_trans, observation_manager_.GetErrorVarianceInverse(),
                    working_matrix_po);
            else
            {
                observation_error_variance R_inv;
                Copy(observation_manager_.GetErrorVariance(), R_inv);
                GetInverse(R_inv);
                Mlt(HL_trans, R_inv, working_matrix_po);
            }

            U_inv_.SetIdentity();
            MltAdd(T(1), SeldonNoTrans, working_matrix_po,
                   SeldonTrans, HL_trans, T(1), U_inv_);
            GetInverse(U_inv_);

            tmp.Reallocate(Nreduced_, Nobservation_);
            tmp.Fill(T(0));

            observation reduced_innovation(Nreduced_);
            MltAdd(T(1), U_inv_, working_matrix_po, T(0), tmp);
            if (!innovation_computation_)
            {
                // Computes innovation.
                observation innovation;
                observation_manager_.GetObservation(innovation);
                Add(T(-1), z, innovation);
                MltAdd(T(1), tmp, innovation, T(0), reduced_innovation);
            }
            else
                MltAdd(T(-1), tmp, z, T(0), reduced_innovation);

            // Updates.
            model_state& x =  model_.GetState();
            MltAdd(T(1), L_, reduced_innovation, T(1), x);
            model_.StateUpdated(); ///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
        }
        else
        {
#if defined(VERDANDI_WITH_MPI) || defined(VERDANDI_WITH_PETSC)
            throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                 "Analyse()", "Parallel algorithm not"
                                 " implemented yet for the 'no"
                                 " simplex' cases.");
#else
            // Computes [HX_{n+1}^{*}].
            sigma_point_matrix Z_i_trans(Nsigma_point_, Nobservation_);
            sigma_point x_col;
            observation z_col, z(Nobservation_);
            z.Fill(T(0));
            if (!alpha_constant_)
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Analyse()", "Calculation not "
                                     "implemented for no constant alpha_i.");

            if (!innovation_computation_)
                for (int i = 0; i < Nsigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    GetRowPointer(Z_i_trans, i, z_col);
                    observation_manager_.ApplyOperator(x_col, z_col);
                    Add(T(alpha_), z_col, z);
                    x_col.Nullify();
                    z_col.Nullify();
                }
            else
                for (int i = 0; i < Nsigma_point_; i++)
                {
                    GetRowPointer(X_i_trans_, i, x_col);
                    observation_manager_.GetInnovation(x_col, z_col);
                    Add(T(alpha_), z_col, z);
                    SetRow(z_col, i, Z_i_trans);
                    x_col.Nullify();
                }

            // Computes [Z] = [HX_{n+1}^{*} - E(HX_{n+1}^{*})].
            for (int i = 0; i < Nsigma_point_; i++)
            {
                GetRowPointer(Z_i_trans, i, z_col);
                Add(T(-1), z, z_col);
                z_col.Nullify();
            }

            sigma_point_matrix
                working_matrix_ro(Nsigma_point_, Nobservation_),
                D_m(Nsigma_point_, Nsigma_point_);
            sigma_point_matrix HL_trans;

            observation_error_variance R_inv;
            if (observation_error_variance_ == "matrix_inverse")
                Mlt(Z_i_trans, observation_manager_.GetErrorVarianceInverse(),
                    working_matrix_ro);
            else
            {
                Copy(observation_manager_.GetErrorVariance(), R_inv);
                GetInverse(R_inv);
                Mlt(Z_i_trans, R_inv, working_matrix_ro);
            }

            // Computes D_m.
            MltAdd(T(1), SeldonNoTrans, working_matrix_ro, SeldonTrans,
                   Z_i_trans, T(0), D_m);

            // Computes U_{n+1}.
            sigma_point_matrix
                working_matrix_rp(Nsigma_point_, Nreduced_),
                working_matrix_rr(Nsigma_point_, Nsigma_point_),
                working_matrix_rr2(Nsigma_point_, Nsigma_point_),
                working_matrix_rr3(Nsigma_point_, Nsigma_point_);

            Copy(D_v_, working_matrix_rr);
            Mlt(T(-1), working_matrix_rr);
            if (alpha_constant_)
            {
                for(int i = 0; i < Nsigma_point_; i++ )
                    working_matrix_rr(i, i) += alpha_;
                MltAdd(T(1), D_m, working_matrix_rr, T(0),
                       working_matrix_rr2);
                for(int i = 0; i < Nsigma_point_; i++ )
                    working_matrix_rr2(i, i) += 1;
                GetInverse(working_matrix_rr2);
                MltAdd(T(1), working_matrix_rr2, D_m, T(0),
                       working_matrix_rr);
                MltAdd(T(alpha_), working_matrix_rr, I_trans_, T(0),
                       working_matrix_rp);
                U_.SetIdentity();
                MltAdd(T(alpha_), SeldonTrans, I_trans_, SeldonNoTrans,
                       working_matrix_rp, T(1), U_);

                Copy(U_, U_inv_);
                GetInverse(U_inv_);

                // Computes {HL}_{n+1}.
                HL_trans.Reallocate(Nreduced_, Nobservation_);
                working_matrix_rr2.SetIdentity();
                MltAdd(T(1), D_v_, working_matrix_rr, T(1),
                       working_matrix_rr2);

                working_matrix_rr.SetIdentity();
                Add(T(alpha_), D_m, working_matrix_rr);
                GetInverse(working_matrix_rr);

                Mlt(working_matrix_rr, working_matrix_rr2,
                    working_matrix_rr3);

                MltAdd(T(alpha_), working_matrix_rr3, I_trans_, T(0),
                       working_matrix_rp);

                MltAdd(T(1), SeldonTrans, working_matrix_rp, SeldonNoTrans,
                       Z_i_trans, T(0), HL_trans);
            }
            else
                throw ErrorUndefined("ReducedOrderUnscentedKalmanFilter::"
                                     "Analyse()", "Calculation not "
                                     "implemented for no constant alpha_i.");

            // Computes K.
            sigma_point_matrix K(Nstate_, Nobservation_),
                working_matrix_po(Nreduced_, Nobservation_),
                working_matrix_po2(Nreduced_, Nobservation_);

            if (observation_error_variance_ == "matrix_inverse")
                Mlt(HL_trans, observation_manager_.GetErrorVarianceInverse(),
                    working_matrix_po);
            else
                Mlt(HL_trans, R_inv, working_matrix_po);
            Mlt(U_inv_, working_matrix_po, working_matrix_po2);
            Mlt(L_, working_matrix_po2, K);

            if (!innovation_computation_)
            {
                // Computes innovation.
                observation innovation;
                observation_manager_.GetObservation(innovation);
                Add(T(-1), z, innovation);

                // Updates.
                model_state& x =  model_.GetState();
                MltAdd(T(1), K, innovation, T(1), x);
                model_.StateUpdated();
            }
            else
            {
                // Updates.
                model_state& x =  model_.GetState();
                MltAdd(T(-1), K, z, T(1), x);
                model_.StateUpdated();
            }
#endif
        }

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            if (option_display_["show_time"] & rank_ == 0)
                cout << " done." << endl;

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");

            MessageHandler::Send(*this, "all", "::Analyze end");

#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Finalizes a step for the model.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::FinalizeStep()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Finalize()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Finalize end");
#if defined(VERDANDI_WITH_MPI)
        MPI::Finalize();
#endif
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::HasFinished()
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class Model, class ObservationManager>
    Model&
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager..
    */
    template <class T, class Model, class ObservationManager>
    ObservationManager&
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::GetObservationManager()
    {
        return observation_manager_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class T, class Model, class ObservationManager>
    OutputSaver&
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string
    ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::GetName() const
    {
        return "ReducedOrderUnscentedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderUnscentedKalmanFilter<T, Model, ObservationManager>
    ::Message(string message)
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            if (message.find("initial condition") != string::npos)
                output_saver_.Save(model_.GetState(), double(model_.GetTime())
                                   , "state_forecast");
            if (message.find("forecast") != string::npos)
                output_saver_.Save(model_.GetState(), double(model_.GetTime())
                                   , "state_forecast");
            if (message.find("analysis") != string::npos)
                output_saver_.Save(model_.GetState(), double(model_.GetTime())
                                   , "state_analysis");
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_REDUCEDORDERUNSCENTEDKALMANFILTER_CXX
#endif
