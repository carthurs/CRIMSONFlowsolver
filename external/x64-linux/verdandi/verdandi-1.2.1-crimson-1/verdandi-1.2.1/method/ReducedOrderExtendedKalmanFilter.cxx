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


#ifndef VERDANDI_FILE_METHOD_REDUCEDORDEREXTENDEDKALMANFILTER_CXX

#include "ReducedOrderExtendedKalmanFilter.hxx"

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
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::ReducedOrderExtendedKalmanFilter()
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
                                         ReducedOrderExtendedKalmanFilter
                                         ::StaticMessage);
#if defined(VERDANDI_WITH_MPI)
        }
#endif
    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::~ReducedOrderExtendedKalmanFilter()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager);
    }


    //! Initializes the driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
            MessageHandler::Send(*this, "all", "::Initialize begin");

        /***************************
         * Reads the configuration *
         ***************************/

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

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


        /*** Ouput saver ***/

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
        {
#endif
            configuration.
                SetPrefix("reduced_order_extended_kalman_filter"
                          ".output_saver.");
            output_saver_.Initialize(configuration);
            output_saver_.Empty("state_forecast");
            output_saver_.Empty("state_analysis");

            /*** Logger and read configuration ***/

            configuration.SetPrefix("reduced_order_extended_kalman_filter.");

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
        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        /*** Initializations ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);
        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

        model_state_error_variance L, U;
        model_.GetStateErrorVarianceSqrt(L, U);
        Copy(L, L_);
        Copy(U, U_);
        Nreduced_ = U.GetN();

        if (Nprocess_ > Nreduced_)
            throw ErrorProcessing("ReducedOrderExtendedKalmanFilter"
                                  "::Initialize(bool initialize_model,"
                                  " bool initialize_observation_manager)",
                                  "The number of processes ("
                                  + to_str(Nprocess_) + ") has to be in [2, "
                                  +  to_str(Nreduced_) + "].");

        /*** Local column indexes of L and U ***/

        Nlocal_reduced_ = int(Nreduced_ / Nprocess_);
        int r = Nreduced_ % Nprocess_;
        if (Nprocess_ - rank_ - 1 < r)
            Nlocal_reduced_++;

        Nlocal_reduced_column_sum_.Reallocate(Nprocess_ + 1);
        Nlocal_reduced_column_sum_(0) = 0;
        for (int i = 0; i < Nprocess_; i++)
        {
            Nlocal_reduced_column_sum_(i + 1) = int(Nreduced_ / Nprocess_);
            if (Nprocess_ - i - 1 < r)
                Nlocal_reduced_column_sum_(i + 1) += 1;
            Nlocal_reduced_column_sum_(i + 1)
                += Nlocal_reduced_column_sum_(i);
        }

        for (int i = Nlocal_reduced_column_sum_(rank_);
             i < Nlocal_reduced_column_sum_(rank_ + 1); i++)
            local_reduced_column_.PushBack(i);

        U_.Reallocate(Nreduced_, Nlocal_reduced_);
        Vector<T> col(Nreduced_);
        for (int i = 0; i < Nlocal_reduced_; i++)
        {
            GetCol(U, local_reduced_column_(i), col);
            SetCol(col, i, U_);
        }

        L_.Reallocate(Nstate_, Nlocal_reduced_);
        col.Reallocate(Nstate_);
        for (int i = 0; i < Nlocal_reduced_; i++)
        {
            GetCol(L, local_reduced_column_(i), col);
            SetCol(col, i, L_);
        }

        displacement_gather_1_ = new int[Nprocess_];
        recvcount_gather_1_ = new int[Nprocess_];
        for (int i = 0; i < Nprocess_; i++)
        {
            recvcount_gather_1_[i] = (Nlocal_reduced_column_sum_(i + 1)
                                      - Nlocal_reduced_column_sum_(i))
                * Nobservation_;
            displacement_gather_1_[i] = Nlocal_reduced_column_sum_(i)
                * Nobservation_;
        }

        displacement_gather_2_ = new int[Nprocess_];
        recvcount_gather_2_ = new int[Nprocess_];
        for (int i = 0; i < Nprocess_; i++)
        {
            recvcount_gather_2_[i] = (Nlocal_reduced_column_sum_(i + 1)
                                      - Nlocal_reduced_column_sum_(i))
                * Nreduced_;
            displacement_gather_2_[i] = Nlocal_reduced_column_sum_(i)
                * Nreduced_;
        }

        displacement_gather_3_ = new int[Nprocess_];
        recvcount_gather_3_ = new int[Nprocess_];
        for (int i = 0; i < Nprocess_; i++)
        {
            recvcount_gather_3_[i] = (Nlocal_reduced_column_sum_(i + 1)
                                      - Nlocal_reduced_column_sum_(i));
            displacement_gather_3_[i] = Nlocal_reduced_column_sum_(i);
        }

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        if (rank_ == 0)
        {
            if (initialize_model)
            {
                MessageHandler::Send(*this, "model", "initial condition");
                MessageHandler::Send(*this, "driver", "initial condition");
            }

            MessageHandler::Send(*this, "all", "::Initialize end");
        }
#else
        MessageHandler::Send(*this, "all", "::Initialize begin");

        /***************************
         * Reads the configuration *
         ***************************/

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

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

        configuration.
            SetPrefix("reduced_order_extended_kalman_filter"
                      ".output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("state_forecast");
        output_saver_.Empty("state_analysis");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration",
                              output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        configuration.SetPrefix("reduced_order_extended_kalman_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        /*** Initializations ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);
        Nstate_ = model_.GetNstate();
        Nobservation_  = observation_manager_.GetNobservation();

        model_state_error_variance L, U;
        model_.GetStateErrorVarianceSqrt(L, U);
        Copy(L, L_);
        Copy(U, U_);
        Nreduced_ = U_.GetN();

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        if (initialize_model)
        {
            MessageHandler::Send(*this, "model", "initial condition");
            MessageHandler::Send(*this, "driver", "initial condition");
        }

        MessageHandler::Send(*this, "all", "::Initialize end");

#endif
    }


    //! Initializes a step for the extended Kalman filter.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::InitializeStep()
    {
        model_.InitializeStep();
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Forward()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            MessageHandler::Send(*this, "all", "::Forward begin");

        time_ = model_.GetTime();

        model_.Forward();

        PropagateCovarianceMatrix();

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
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::Analyze()
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
            MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            if (rank_ == 0)
                if (option_display_["show_time"])
                    cout << "Performing Reduced Order EKF at time step ["
                         << model_.GetTime() << "]..." << endl;

            model_state& x = model_.GetState();
            Nstate_ = model_.GetNstate();

            /*** Updated matrix U ***/

            dense_matrix HL_local_trans(Nobservation_, Nlocal_reduced_),
                working_matrix_or(Nobservation_, Nlocal_reduced_),
                HL_global_trans(Nreduced_, Nobservation_),
                U_global(Nreduced_, Nreduced_);


            Mlt(observation_manager_.GetTangentLinearOperator(), L_,
                HL_local_trans);
            Transpose(HL_local_trans);
            MPI::COMM_WORLD.Allgatherv(HL_local_trans.GetData(),
                                       Nlocal_reduced_ * Nobservation_,
                                       MPI::DOUBLE, HL_global_trans.GetData(),
                                       recvcount_gather_1_,
                                       displacement_gather_1_, MPI::DOUBLE);
            Transpose(HL_local_trans);
            Mlt(observation_manager_.GetErrorVarianceInverse(),
                HL_local_trans, working_matrix_or);
            MltAdd(T(1), HL_global_trans, working_matrix_or, T(1), U_);

            Transpose(U_);
            MPI::COMM_WORLD.Allgatherv(U_.GetData(),
                                       Nlocal_reduced_ * Nreduced_,
                                       MPI::DOUBLE, U_global.GetData(),
                                       recvcount_gather_2_,
                                       displacement_gather_2_, MPI::DOUBLE);
            Transpose(U_);

            observation y;
            observation_manager_.GetInnovation(x, y);
            Nobservation_ = y.GetSize();

            model_state_error_variance_row state_innovation(Nreduced_),
                state_innovation_local(Nlocal_reduced_),
                x_local(Nreduced_);
            MltAdd(T(1), SeldonTrans, working_matrix_or, y,
                   T(0), state_innovation_local);

            MPI::COMM_WORLD.Allgatherv(state_innovation_local.GetData(),
                                       Nlocal_reduced_,
                                       MPI::DOUBLE,
                                       state_innovation.GetData(),
                                       recvcount_gather_3_,
                                       displacement_gather_3_, MPI::DOUBLE);

            GetInverse(U_global);
            MltAdd(T(1), U_global, state_innovation, T(0), x_local);

            model_state_error_variance_row x_local_local(Nlocal_reduced_);
            for (int i = 0; i < Nlocal_reduced_; i++)
                x_local_local(i) =
                    x_local(Nlocal_reduced_column_sum_(rank_) + i);

            model_state_error_variance_row dx_local(Nstate_), dx(Nstate_);
            MltAdd(T(1), L_, x_local_local, T(0), dx_local);
            dx.Fill(T(0));
            MPI::COMM_WORLD.Allreduce(dx_local.GetData(), dx.GetData(),
                                      Nstate_, MPI::DOUBLE,
                                      MPI::SUM);
            Add(T(1), dx, x);

            model_.StateUpdated();


            if (rank_ == 0)
            {
                if (option_display_["show_time"])
                    cout << " done." << endl;

                MessageHandler::Send(*this, "model", "analysis");
                MessageHandler::Send(*this, "observation_manager",
                                     "analysis");
                MessageHandler::Send(*this, "driver", "analysis");
            }
        }
        if (rank_ == 0)
            MessageHandler::Send(*this, "all", "::Analyze end");
#else
        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            if (option_display_["show_time"])
                cout << "Performing Reduced Order EKF at time step ["
                     << model_.GetTime() << "]..." << endl;

            model_state& x = model_.GetState();
            Nstate_ = model_.GetNstate();

            observation y;
            observation_manager_.GetInnovation(x, y);
            Nobservation_ = y.GetSize();

            /*** Updated matrix U ***/

            dense_matrix HL(Nobservation_, Nreduced_),
                working_matrix_or(Nobservation_, Nreduced_);
            Mlt(observation_manager_.GetTangentLinearOperator(), L_, HL);
            Mlt(observation_manager_.GetErrorVarianceInverse(), HL,
                working_matrix_or);
            MltAdd(T(1), SeldonTrans, HL, SeldonNoTrans,
                   working_matrix_or, T(1), U_);

            /*** Updated K ***/

            dense_matrix U_inv(U_),
                working_matrix_ro2(Nreduced_, Nobservation_);

            Mlt(observation_manager_.GetErrorVarianceInverse(), HL,
                working_matrix_or);
            model_state_error_variance_row state_innovation(Nreduced_);
            MltAdd(T(1), SeldonTrans, working_matrix_or, y, T(0),
                   state_innovation);

            model_state_error_variance_row correction(Nreduced_);
            GetInverse(U_inv);
            MltAdd(T(1), U_inv, state_innovation, T(0), correction);
            MltAdd(T(1), L_, correction, T(1), x);

            model_.StateUpdated();

            if (option_display_["show_time"])
                cout << " done." << endl;

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
#endif
    }


    //! Finalizes a step for the model.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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


    //! Computes Covariance.
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::PropagateCovarianceMatrix()
    {
#ifdef VERDANDI_WITH_MPI
        double saved_time = model_.GetTime();
        model_.SetTime(time_);

        // One column of L.
        model_state_error_variance_row L_col(Nstate_);
        for (int j = 0; j < Nlocal_reduced_; j++)
        {
            GetCol(L_, j, L_col);
            model_.ApplyTangentLinearOperator(L_col);
            SetCol(L_col, j, L_);
        }

        model_.SetTime(saved_time);
#else
        double saved_time = model_.GetTime();
        model_.SetTime(time_);

        // One column of L.
        model_state_error_variance_row L_col(Nstate_);
        for (int j = 0; j < Nreduced_; j++)
        {
            GetCol(L_, j, L_col);
            model_.ApplyTangentLinearOperator(L_col);
            SetCol(L_col, j, L_);
        }

        model_.SetTime(saved_time);
#endif
    }



    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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
    ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
    ::GetName() const
    {
        return "ReducedOrderExtendedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void ReducedOrderExtendedKalmanFilter<T, Model, ObservationManager>
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


#define VERDANDI_FILE_METHOD_REDUCEDORDEREXTENDEDKALMANFILTER_CXX
#endif
