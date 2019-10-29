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


#ifndef VERDANDI_FILE_METHOD_OPTIMALINTERPOLATION_CXX


#include "OptimalInterpolation.hxx"
#include "BLUE.cxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    template <class T, class Model, class ObservationManager>
    OptimalInterpolation<T, Model, ObservationManager>
    ::OptimalInterpolation()
    {

        /*** Initializations ***/

#if defined(VERDANDI_WITH_MPI)
        MPI::Init();
	rank_ = MPI::COMM_WORLD.Get_rank();
#endif
        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     OptimalInterpolation::StaticMessage);
    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    OptimalInterpolation<T, Model, ObservationManager>
    ::~OptimalInterpolation()
    {
#if defined(VERDANDI_WITH_MPI)
        MPI::Finalize();
#endif
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the optimal interpolation driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager);
    }


    //! Initializes the optimal interpolation driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");
        configuration_file_ = configuration.GetFilePath();


        /***************************
         * Reads the configuration *
         ***************************/


        configuration.SetPrefix("optimal_interpolation.");

        /*** Model ***/

        configuration.Set("model.configuration_file", "",
                          configuration_file_,
                          model_configuration_file_);

        /*** Observation manager ***/

        configuration.Set("observation_manager.configuration_file", "",
                          configuration_file_,
                          observation_configuration_file_);

        /*** Display options ***/

        configuration.SetPrefix("optimal_interpolation.display.");
        // Should iterations be displayed on screen?
        configuration.Set("show_iteration",
                          option_display_["show_iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("show_time", option_display_["show_time"]);

        /*** Assimilation options ***/

        configuration.
            SetPrefix("optimal_interpolation.data_assimilation.");
        configuration.Set("analyze_first_step", analyze_first_step_);

        configuration.SetPrefix("optimal_interpolation.");
        configuration.Set("BLUE_computation",
#ifdef VERDANDI_WITH_DIRECT_SOLVER
                          "ops_in(v, {'vector', 'matrix'})",
#else
                          "ops_in(v, {'vector'})",
#endif
                          blue_computation_);

        /*** Ouput saver ***/

        configuration.SetPrefix("optimal_interpolation.output_saver.");
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
	{
#endif
            output_saver_.Initialize(configuration);
            output_saver_.Empty("state_forecast");
            output_saver_.Empty("state_analysis");

#if defined(VERDANDI_WITH_MPI)
	}
#endif
        /*** Logger and read configuration ***/

        configuration.SetPrefix("optimal_interpolation.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        /*** Initializations ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        if (initialize_model && initialize_observation_manager)
        {
            MessageHandler::Send(*this, "model", "initial condition");
            MessageHandler::Send(*this, "driver", "initial condition");
        }

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the optimal interpolation.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        if (option_display_["show_time"])
            cout << "Current step: "
                 << model_.GetTime() << endl;
        model_.InitializeStep();

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward, with optimal interpolation at the end.
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        model_.Forward();

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");


        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE.
     */
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>::Analyze()
    {
        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            if (option_display_["show_time"])
                cout << "Performing optimal interpolation at time step ["
                     << model_.GetTime() << "]..." << endl;

            model_state& state = model_.GetState();
            Nstate_ = model_.GetNstate();

            observation innovation;
            observation_manager_.GetInnovation(state, innovation);
            Nobservation_ = innovation.GetSize();

            matrix_state_observation matrix_state_observation_tmp;
            if (blue_computation_ == "vector")
                ComputeBLUE_vector(model_, observation_manager_,
                                   innovation, state);
#ifdef VERDANDI_WITH_DIRECT_SOLVER
            else
                ComputeBLUE_matrix(model_.GetStateErrorVariance(),
                                   observation_manager_
                                   .GetTangentLinearOperator(),
                                   matrix_state_observation_tmp,
                                   innovation,
                                   observation_manager_
                                   .GetErrorVariance(),
                                   state, true, false);
#endif

            model_.StateUpdated();

            if (option_display_["show_time"])
                cout << " done." << endl;

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Finalizes a step for the model.
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool OptimalInterpolation<T, Model, ObservationManager>::HasFinished()
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class Model, class ObservationManager>
    Model&
    OptimalInterpolation<T, Model, ObservationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class T, class Model, class ObservationManager>
    ObservationManager&
    OptimalInterpolation<T, Model, ObservationManager>
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
    OptimalInterpolation<T, Model, ObservationManager>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string OptimalInterpolation<T, Model, ObservationManager>::GetName() const
    {
        return "OptimalInterpolation";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void OptimalInterpolation<T, Model, ObservationManager>
    ::Message(string message)
    {
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
	{
#endif
            if (message.find("initial condition") != string::npos)
                output_saver_.Save(model_.GetState(), double(model_.GetTime()),
                                   "state_forecast");
            if (message.find("forecast") != string::npos)
                output_saver_.Save(model_.GetState(), model_.GetTime(),
                                   "state_forecast");
            if (message.find("analysis") != string::npos)
                output_saver_.Save(model_.GetState(), model_.GetTime(),
                                   "state_analysis");
#if defined(VERDANDI_WITH_MPI)
	}
#endif
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_OPTIMALINTERPOLATION_CXX
#endif
