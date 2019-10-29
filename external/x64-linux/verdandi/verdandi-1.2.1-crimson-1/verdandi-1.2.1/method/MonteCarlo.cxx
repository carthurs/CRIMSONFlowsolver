// Copyright (C) 2009 INRIA
// Author(s): Anne Tilloy, Vivien Mallet
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


#ifndef VERDANDI_FILE_METHOD_MONTECARLO_CXX


#include "MonteCarlo.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration_file configuration file.
    */
    template <class T, class ClassModel, class PerturbationManager>
    MonteCarlo<T, ClassModel, PerturbationManager>::MonteCarlo():
        iteration_(-1)
    {
        MessageHandler::AddRecipient("model", model_,
                                     ClassModel::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     MonteCarlo::StaticMessage);
    }


    //! Destructor.
    template <class T, class ClassModel, class PerturbationManager>
    MonteCarlo<T, ClassModel, PerturbationManager>::~MonteCarlo()
    {
        typename vector<uncertain_parameter>::iterator i;
        for (i = perturbation_.begin(); i != perturbation_.end(); i++)
            Clear(*i);
    }


    //! Clears a vector.
    /*!
      \param[in,out] V vector to be cleared.
    */
    template <class T, class ClassModel, class PerturbationManager>
    template <class T0, class Storage0, class Allocator0>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::Clear(Vector<T0, Storage0, Allocator0>& V)
    {
        V.Clear();
    }


    //! Clears a vector collection.
    /*! Each inner vector of the collection is cleared.
      \param[in,out] V vector to be cleared.
    */
    template <class T, class ClassModel, class PerturbationManager>
    template <class T0, class Allocator0>
    void
    MonteCarlo<T, ClassModel, PerturbationManager>
    ::Clear(Vector<T0, Collection, Allocator0>& V)
    {
        for (int i = 0; i < V.GetNvector(); i++)
            V.GetVector(i).Clear();
    }


    //////////////////
    // MAIN METHODS //
    //////////////////


    //! Allocates an output vector to the dimension of the input vector.
    /*!
      \param[in] in input vector.
      \param[out] out output vector.
    */
    template <class T, class ClassModel, class PerturbationManager>
    template <class T0, class Storage0, class Allocator0>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::SetDimension(Vector<T0, Storage0, Allocator0>& in,
                   Vector<T0, Storage0, Allocator0>& out)
    {
        out.Reallocate(in.GetLength());
    }


    /*! \brief Allocates an output vector collection to the dimension of the
      input vector collection. */
    /*!
      \param[in] in input collection vector.
      \param[out] out output collection vector.
    */
    template <class T, class ClassModel, class PerturbationManager>
    template <class T0, class Allocator0>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::SetDimension(Vector<T0, Collection, Allocator0>& in,
                   Vector<T0, Collection, Allocator0>& out)
    {
        T0 suboutput;
        for (int i = 0; i < in.GetNvector(); i++)
        {
            suboutput.Reallocate(in.GetVector(i).GetLength());
            out.AddVector(suboutput);
            suboutput.Nullify();
        }
    }


    /*! \brief Fills an input vector collection according to its probability
      distribution. */
    /*!
      \param[in,out] in input collection vector.
      \param[in] pdf probability density function: Normal, NormalHomogeneous,
      LogNormal or LogNormalHomogeneous.
    */
    template <class T, class ClassModel, class PerturbationManager>
    template <class T0, class Allocator0>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::Fill(Vector<T0, Collection, Allocator0>& in,
           string pdf)
    {
        if (pdf == "Normal" || pdf == "NormalHomogeneous")
            for (int i = 0; i < in.GetNvector(); i++)
                in.GetVector(i).Fill(typename T0::value_type(0));
        else if (pdf == "LogNormal" || pdf == "LogNormalHomogeneous")
            for (int i = 0; i < in.GetNvector(); i++)
                in.GetVector(i).Fill(typename T0::value_type(1));
    }


    /*! \brief Fills an input vector according to its probability
      distribution. */
    /*!
      \param[in,out] in input vector.
      \param[in] pdf probability density function: Normal, NormalHomogeneous,
      LogNormal or LogNormalHomogeneous.
    */
    template <class T, class ClassModel, class PerturbationManager>
    template <class T0, class Storage0, class Allocator0>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::Fill(Vector<T0, Storage0, Allocator0>& in, string pdf)
    {
        if (pdf == "Normal" || pdf == "NormalHomogeneous")
            in.Fill(T0(0));
        else if (pdf == "LogNormal" || pdf == "LogNormalHomogeneous")
            in.Fill(T0(1));
    }


    //! Initializes the simulation.
    /*! It reads the configuration and initializes the model and the
      perturbation manager.
      \param[in] configuration_file configuration file for the method.
      \param[in] initialize_model should the model be initialized with a call
      to Model::Initialize(string)?
      \param[in] initialize_perturbation_manager should the perturbation
      manager be initialized with a call to
      PerturbationManager::Initialize(string)?
      \warning If \a initialize_model is set to false, the model should be
      initialized before calling this function.
    */
    template <class T, class ClassModel, class PerturbationManager>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model,
                 bool initialize_perturbation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration);
    }


    //! Initializes the simulation.
    /*! It reads the configuration and initializes the model and the
      perturbation manager.
      \param[in] configuration_file configuration file for the method.
      \param[in] initialize_model should the model be initialized with a call
      to Model::Initialize(string)?
      \param[in] initialize_perturbation_manager should the perturbation
      manager be initialized with a call to
      PerturbationManager::Initialize(string)?
      \warning If \a initialize_model is set to false, the model should be
      initialized before calling this function.
    */
    template <class T, class ClassModel, class PerturbationManager>
    void MonteCarlo<T, ClassModel, PerturbationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model,
                 bool initialize_perturbation_manager)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("monte_carlo.");


        /***************************
         * Reads the configuration *
         ***************************/


        /*** Model ***/

        configuration.Set("model.configuration_file", "",
                          configuration_file_,
                          model_configuration_file_);

        /*** Perturbation managager ***/

        configuration.Set("perturbation_manager.configuration_file", "",
                          configuration_file_,
                          perturbation_manager_configuration_file_);

        /*** Display options ***/

        // Should the iteration be displayed on screen?
        configuration.Set("display.show_iteration", show_iteration_);
        // Should the time be displayed on screen?
        configuration.Set("display.show_time", show_time_);

        /*** Ouput saver ***/

        configuration.SetPrefix("monte_carlo.output_saver.");

        output_saver_.Initialize(configuration);
        output_saver_.Empty("perturbation");
        output_saver_.Empty("state_forecast");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("monte_carlo.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }


        if (show_time_)
            Logger::StdOut(*this, "Time: " + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this, "Time: " + to_str(model_.GetTime()));
        if (show_iteration_)
            Logger::StdOut(*this, "Initialization");
        else
            Logger::Log<-3>(*this, "Initialization");

        if (initialize_model)
            model_.Initialize(model_configuration_file_);

        if (initialize_perturbation_manager)
            perturbation_manager_
                .Initialize(perturbation_manager_configuration_file_);

        iteration_ = 0;

        MessageHandler::Send(*this, "model", "initial condition");
        for (int i = 0; i < model_.GetNparameter(); i++)
        {
            uncertain_parameter output;
            bool allocate;
            uncertain_parameter parameter = model_.GetParameter(i);

            if (model_.GetParameterPDF(i) == "Normal"
                || model_.GetParameterPDF(i) == "LogNormal"
                || model_.GetParameterPDF(i) == "BlockNormal"
                || model_.GetParameterPDF(i) == "BlockLogNormal")
            {
                SetDimension(parameter, output);
                Fill(output, model_.GetParameterPDF(i));
                perturbation_manager_
                    .Sample(model_.GetParameterPDF(i),
                            model_.GetParameterVariance(i),
                            model_.GetParameterParameter(i),
                            model_.GetParameterCorrelation(i),
                            output);
            }
            else if (model_.GetParameterPDF(i) == "NormalHomogeneous"
                     || model_.GetParameterPDF(i) == "LogNormalHomogeneous"
                     || model_.GetParameterPDF(i) == "BlockNormalHomogeneous"
                     || model_.GetParameterPDF(i)
                     == "BlockLogNormalHomogeneous")
            {
                SetDimension(parameter, output);
                Fill(output, model_.GetParameterPDF(i));
                perturbation_manager_
                    .Sample(model_.GetParameterPDF(i),
                            model_.GetParameterVariance(i)(0, 0),
                            model_.GetParameterParameter(i),
                            model_.GetParameterCorrelation(i),
                            output);
            }
            else
                throw ErrorConfiguration("MonteCarlo::Initialize(string)",
                                         "The probability distribution \""
                                         + model_.GetParameterPDF(i)
                                         + "\" is not supported.");
            perturbation_.push_back(output);

            if (model_.GetParameterOption(i) == "init_step")
                // Applies the perturbations.
                if (model_.GetParameterPDF(i) == "Normal"
                    || model_.GetParameterPDF(i) == "BlockNormal"
                    || model_.GetParameterPDF(i) == "NormalHomogeneous"
                    || model_.GetParameterPDF(i) == "BlockNormalHomogeneous")
                    Add(1., perturbation_[i], parameter);
                else if (model_.GetParameterPDF(i) == "LogNormal"
                         || model_.GetParameterPDF(i) == "BlockLogNormal"
                         || model_.GetParameterPDF(i)
                         == "LogNormalHomogeneous"
                         || model_.GetParameterPDF(i)
                         == "BlockLogNormalHomogeneous")
                    for (int k = 0; k < perturbation_[i].GetM(); k++)
                        parameter(k) *= perturbation_[i](k);

            model_.SetParameter(i, parameter);

            output_saver_.Save(perturbation_[i], "perturbation");
        }

        perturbation_manager_.Finalize();

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes the model before a time step.
    template <class T, class ClassModel, class PerturbationManager>
    void MonteCarlo<T, ClassModel, PerturbationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");
        model_.InitializeStep();
        for (int i = 0; i < model_.GetNparameter(); i++)
        {
            uncertain_parameter parameter = model_.GetParameter(i);
            if (model_.GetParameterOption(i) == "every_step")
                if (model_.GetParameterPDF(i) == "Normal"
                    || model_.GetParameterPDF(i) == "BlockNormal"
                    || model_.GetParameterPDF(i) == "NormalHomogeneous"
                    || model_.GetParameterPDF(i) == "BlockNormalHomogeneous")
                    Add(1., perturbation_[i], parameter);
                else if (model_.GetParameterPDF(i) == "LogNormal"
                         || model_.GetParameterPDF(i) == "BlockLogNormal"
                         || model_.GetParameterPDF(i)
                         == "LogNormalHomogeneous"
                         || model_.GetParameterPDF(i)
                         == "BlockLogNormalHomogeneous")
                    for (int k = 0; k < perturbation_[i].GetM(); k++)
                        parameter(k) *= perturbation_[i](k);

            model_.SetParameter(i, parameter);
        }

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward without optimal interpolation.
    template <class T, class ClassModel, class PerturbationManager>
    void MonteCarlo<T, ClassModel, PerturbationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        time_.PushBack(model_.GetTime());

        if (show_time_)
            Logger::StdOut(*this, "Time: " + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this, "Time: " + to_str(model_.GetTime()));
        if (show_iteration_)
            Logger::StdOut(*this, "Iteration " + to_str(iteration_) + " -> "
                           + to_str(iteration_ + 1));
        else
            Logger::Log<-3>(*this, "Iteration " + to_str(iteration_) + " -> "
                            + to_str(iteration_ + 1));

        model_.Forward();
        iteration_++;
        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Finalizes a step for the model.
    template <class T, class ClassModel, class PerturbationManager>
    void MonteCarlo<T, ClassModel, PerturbationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class T, class ClassModel, class PerturbationManager>
    void MonteCarlo<T, ClassModel, PerturbationManager>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if the simulation is finished, false otherwise.
    */
    template <class T, class ClassModel, class PerturbationManager>
    bool MonteCarlo<T, ClassModel, PerturbationManager>::HasFinished()
    {
        return model_.HasFinished();
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class ClassModel, class PerturbationManager>
    ClassModel& MonteCarlo<T, ClassModel, PerturbationManager>
    ::GetModel()
    {
        return model_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class T, class ClassModel, class PerturbationManager>
    OutputSaver& MonteCarlo<T, ClassModel, PerturbationManager>
    ::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class ClassModel, class PerturbationManager>
    string MonteCarlo<T, ClassModel, PerturbationManager>::GetName() const
    {
        return "MonteCarlo";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class ClassModel, class PerturbationManager>
    void  MonteCarlo<T, ClassModel, PerturbationManager>
    ::Message(string message)
    {
        if (message.find("forecast") != string::npos)
        {
            output_saver_.Save(model_.GetState(), model_.GetTime(), "state_forecast");
        }
    }


} // namespace Verdandi


#define VERDANDI_FILE_METHOD_MONTECARLO_CXX
#endif
