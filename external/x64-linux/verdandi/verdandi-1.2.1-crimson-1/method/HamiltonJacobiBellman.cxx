// Copyright (C) 2009-2010 INRIA
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
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_HAMILTONJACOBIBELLMAN_CXX


#include "HamiltonJacobiBellman.hxx"


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
    HamiltonJacobiBellman<T, Model, ObservationManager>
    ::HamiltonJacobiBellman():
        time_step_(0)
    {

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     HamiltonJacobiBellman::StaticMessage);

    }


    //! Destructor.
    template <class T, class Model, class ObservationManager>
    HamiltonJacobiBellman<T, Model, ObservationManager>
    ::~HamiltonJacobiBellman()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the solver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step.
      \param[in] configuration_file configuration file.
    */
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>
    ::Initialize(string configuration_file)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration);
    }


    //! Initializes the solver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step.
      \param[in] configuration_file configuration file.
    */
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration)
    {
        string configuration_file = configuration.GetFilePath();
        configuration.SetPrefix("hjb.");
        MessageHandler::Send(*this, "all", "::Initialize begin");

        if (option_display_["show_time"])
            Logger::StdOut(*this, "Time: "
                           + to_str(T(time_step_) * Delta_t_));
        else
            Logger::Log<-3>(*this, "Time: "
                            + to_str(T(time_step_) * Delta_t_));
        if (option_display_["show_iteration"])
            Logger::StdOut(*this, "Iteration " + to_str(time_step_) + " -> "
                           + to_str(time_step_ + 1));
        else
            Logger::Log<-3>(*this, "Iteration " + to_str(time_step_) + " -> "
                            + to_str(time_step_ + 1));


        /***************************
         * Reads the configuration *
         ***************************/


        /*** Display options ***/

        configuration.SetPrefix("hjb.display.");
        // Should iterations be displayed on screen?
        configuration.Set("show_iteration",
                          option_display_["show_iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("show_time", option_display_["show_time"]);

        /*** Domain definition ***/

        configuration.SetPrefix("hjb.domain.");

        // Discretization along all dimensions. A regular grid is assumed. In
        // every dimension, the format is x_min, delta_x, Nx for every
        // dimension.
        vector<string> discretization_vector;
        configuration.Set("discretization", discretization_vector);
        if (discretization_vector.size() % 3 != 0)
            throw ErrorConfiguration("HamiltonJacobiBellman::"
                                     "HamiltonJacobiBellman",
                                     "The entry \"discretization\" should be "
                                     "in format \"x_min delta_x Nx\" for "
                                     "every dimension.");
        Ndimension_ = discretization_vector.size() / 3;
        x_min_.Reallocate(Ndimension_);
        Delta_x_.Reallocate(Ndimension_);
        Nx_.Reallocate(Ndimension_);
        Npoint_ = 1;
        for (int i = 0; i < Ndimension_; i++)
        {
            to_num(discretization_vector[3 * i], x_min_(i));
            to_num(discretization_vector[3 * i + 1], Delta_x_(i));
            to_num(discretization_vector[3 * i + 2], Nx_(i));
            Npoint_ *= Nx_(i);
        }

        // Checks consistency of 'Ndimension_' with the model state.
        if (Ndimension_ != model_.GetNstate())
            throw ErrorConfiguration("HamiltonJacobiBellman::"
                                     "HamiltonJacobiBellman",
                                     "The dimension of the model ("
                                     + to_str(model_.GetNstate())
                                     + ") is incompatible with that of "
                                     " the HJB solver (" + to_str(Ndimension_)
                                     + ").");

        configuration.Set("initial_time", initial_time_);
        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("Nt", Nt_);

        /*** Equation coefficients ***/

        configuration.SetPrefix("hjb.equation_coefficient.");

        configuration.Set("with_quadratic_term", with_quadratic_term_);
        configuration.Set("with_advection_term", with_advection_term_);
        configuration.Set("with_source_term", with_source_term_);

        if (with_advection_term_)
            configuration.Set("model_time_dependent", model_time_dependent_);

        Q_0_.Reallocate(Ndimension_, Ndimension_);
        configuration.Set("Q_0", Q_0_);
        if (Q_0_.GetM() != Ndimension_ || Q_0_.GetN() != Ndimension_)
            throw ErrorConfiguration("HamiltonJacobiBellman::"
                                     "HamiltonJacobiBellman",
                                     "The entry \"Q_0\" should be "
                                     "a " + to_str(Ndimension_) + " x "
                                     + to_str(Ndimension_) + " matrix, "
                                     " but a " + to_str(Q_0_.GetM()) + " x "
                                     + to_str(Q_0_.GetN())
                                     + " matrix was provided.");

        x_0_.Reallocate(Ndimension_);
        configuration.Set("x_0", x_0_);
        if (x_0_.GetLength() != Ndimension_)
            throw ErrorConfiguration("HamiltonJacobiBellman::"
                                     "HamiltonJacobiBellman",
                                     "The entry \"x_0\" should contain "
                                     + to_str(Ndimension_) + " elements, but "
                                     + to_str(x_0_.GetLength())
                                     + " elements were provided.");

        if (with_quadratic_term_)
        {
            Q_inv_.Reallocate(Ndimension_, Ndimension_);
            configuration.Set("Q", Q_inv_);
            if (Q_inv_.GetM() != Ndimension_ || Q_inv_.GetN() != Ndimension_)
                throw ErrorConfiguration("HamiltonJacobiBellman::"
                                         "HamiltonJacobiBellman",
                                         "The entry \"Q\" should be "
                                         "a " + to_str(Ndimension_) + " x "
                                         + to_str(Ndimension_) + " matrix, "
                                         " but a " + to_str(Q_inv_.GetM())
                                         + " x " + to_str(Q_inv_.GetN())
                                         + " matrix was provided.");
            // Q must be diagonal, and its inverse is therefore trivial to
            // compute.
            for (int i = 0; i <  Ndimension_; i++)
                for (int j = 0; j <  Ndimension_; j++)
                    if (i != j && Q_inv_(i, j) != T(0))
                        throw ErrorConfiguration("HamiltonJacobiBellman::"
                                                 "HamiltonJacobiBellman",
                                                 "The entry \"Q\" should be "
                                                 "a diagonal matrix. Non-"
                                                 "diagonal matrix is not"
                                                 " supported.");
            for (int i = 0; i <  Ndimension_; i++)
                if (Q_inv_(i, i) <= T(0))
                    throw ErrorConfiguration("HamiltonJacobiBellman::"
                                             "HamiltonJacobiBellman",
                                             "The entry \"Q\" should be "
                                             "positive definite, but its "
                                             + to_str(i) + "-th diagonal "
                                             "element is "
                                             + to_str(Q_inv_(i, i)) + ".");
                else
                    Q_inv_(i, i) = Delta_t_
                        / (Q_inv_(i, i) * Delta_x_(i) * Delta_x_(i));
        }
        else
        {
            // Even when the quadratic term is not taken into account,
            // 'Q_inv_' is still needed. Deactivation of the quadratic term is
            // achieved with a diagonal of zeros in 'Q_inv_'.
            Q_inv_.Reallocate(Ndimension_, Ndimension_);
            Q_inv_.Zero();
        }

        if (with_source_term_)
        {
            R_.Reallocate(Ndimension_, Ndimension_);
            configuration.Set("R", R_);
            if (R_.GetM() != Ndimension_ || R_.GetN() != Ndimension_)
                throw ErrorConfiguration("HamiltonJacobiBellman::"
                                         "HamiltonJacobiBellman",
                                         "The entry \"R\" should be "
                                         "a " + to_str(Ndimension_) + " x "
                                         + to_str(Ndimension_) + " matrix, "
                                         " but a " + to_str(R_.GetM())
                                         + " x " + to_str(R_.GetN())
                                         + " matrix was provided.");
        }

        /*** Solver ***/

        configuration.SetPrefix("hjb.solver.");

        configuration.Set("scheme",
                          "ops_in(v, {'LxF', 'BrysonLevy', 'Godunov'})",
                          scheme_);
        if (with_quadratic_term_ && scheme_ != "Godunov")
            throw ErrorConfiguration("HamiltonJacobiBellman::"
                                     "HamiltonJacobiBellman",
                                     "The quadratic term can only be taken "
                                     "into account in Godunov scheme.");

        /*** Boundary condition ***/

        configuration.SetPrefix("hjb.boundary_condition.");

        configuration.Set("type",
                          "ops_in(v, {'Dirichlet', "
                          "'Extrapolation', 'Periodic'})",
                          boundary_condition_type_);
        if (boundary_condition_type_ == "Dirichlet")
            boundary_condition_index_ = 0;
        else if (boundary_condition_type_ == "Extrapolation")
            boundary_condition_index_ = 1;
        else
            boundary_condition_index_ = 2;
        configuration.Set("value", "v >= 0", boundary_condition_);

        /*** Lax-Friedrichs scheme ***/

        if (scheme_ == "LxF")
        {
            configuration.SetPrefix("hjb.lax_friedrichs.");

            configuration.Set("Upper_bound_model", upper_bound_model_);
            if (upper_bound_model_.GetLength() != Ndimension_)
                throw ErrorConfiguration("HamiltonJacobiBellman::"
                                         "HamiltonJacobiBellman",
                                         "The entry \"upper_bound_model\" "
                                         "should contain "
                                         + to_str(Ndimension_)
                                         + " elements, but "
                                         + to_str(upper_bound_model_
                                                  .GetLength())
                                         + " elements were provided.");
        }

        /*** Ouput saver ***/

        configuration.SetPrefix("hjb.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("value_function");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("hjb.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        /*** Initializations ***/

        if (with_advection_term_)
            model_.Initialize(configuration_file);
        if (with_source_term_)
            observation_manager_.Initialize(model_, configuration_file);

        /*** Initial value function ***/

        V_.Reallocate(Npoint_);

        // Initial value function: V(0, x) = <Q_0 x, x>.
        Vector<T> x(Ndimension_), Qx(Ndimension_);
        for (int i = 0; i < Npoint_; i++)
        {
            get_coordinate(i, x_min_, Delta_x_, Nx_, x);
            Add(T(-1.), x_0_, x);
            Mlt(Q_0_, x, Qx);
            V_(i) = DotProd(Qx, x);
        }

        /*** Model ***/

        Vector<T> Mx;
        if (with_advection_term_)
            Mx_.Reallocate(Npoint_, Ndimension_);
        if (with_advection_term_ && !model_time_dependent_)
        {
            courant_number_ = 0.;
            double time, time_step;
            for (int i_cell = 0; i_cell < Npoint_; i_cell++)
            {
                get_coordinate(i_cell, x_min_, Delta_x_, Nx_, x);

                model_.SetTime(initial_time_);
                model_.SetState(x);
                model_.Forward();
                model_.GetState(Mx);
                time = model_.GetTime();
                time_step = time - initial_time_;

                Add(-1., x, Mx);
                for (int d = 0; d < Ndimension_; d++)
                {
                    Mx(d) *= Delta_t_ / (Delta_x_(d) * time_step);
                    courant_number_ = max(courant_number_, abs(Mx(d)));
                }

                SetRow(Mx, i_cell, Mx_);
            }

            Logger::Log(*this, "Courant number: " + to_str(courant_number_));
        }

        /*** Evolution points ***/

        if (with_advection_term_ && scheme_ == "BrysonLevy")
        {
            // Location of the evolution points: (a, a, a, ..., a).
            T a = 0;
            for (int d = 0; d < Ndimension_; d++)
                a += 1. / (Delta_x_(d) * Delta_x_(d));
            a = sqrt(a);
            for (int d = 0; d < Ndimension_; d++)
                a += 1. / Delta_x_(d);
            a = 1. / a;

            // 'a' over Delta_x.
            a_Delta_x_.Reallocate(Ndimension_);
            for (int d = 0; d < Ndimension_; d++)
                a_Delta_x_ = a / Delta_x_(d);
        }

        MessageHandler::Send(*this, "all", "initial value");

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the time integration of HJB equation.
    /*! Initializes a step for the model.
     */
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        if (option_display_["show_time"])
            Logger::StdOut(*this, "Time: "
                           + to_str(T(time_step_) * Delta_t_));
        else
            Logger::Log<-3>(*this, "Time: "
                            + to_str(T(time_step_) * Delta_t_));
        if (option_display_["show_iteration"])
            Logger::StdOut(*this, "Iteration " + to_str(time_step_) + " -> "
                           + to_str(time_step_ + 1));
        else
            Logger::Log<-3>(*this, "Iteration " + to_str(time_step_) + " -> "
                            + to_str(time_step_ + 1));

        if (with_advection_term_)
            model_.InitializeStep();

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward.
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        if (with_advection_term_ || with_quadratic_term_)
            if (scheme_ == "LxF")
                AdvectionLxFForward();
            else if (scheme_ == "BrysonLevy")
                AdvectionBrysonLevyForward();
            else
                AdvectionGodunov();

        /*** Source term (norm of the innovation) ***/

        double time = initial_time_ + T(time_step_) * Delta_t_;
        if (with_source_term_)
        {
            model_.SetTime(time);
            observation_manager_.SetTime(model_, time);
        }
        if (with_source_term_ && observation_manager_.HasObservation())
        {
            Vector<T> x(Ndimension_);

            int Nobservation = observation_manager_.GetNobservation();
            Vector<T> y(Nobservation), innovation(Nobservation),
                Rinnovation(Nobservation);

            observation_manager_.GetObservation(y);

            x = x_min_;
            Vector<int> position(Ndimension_);
            position.Fill(0);
            for (int i_cell = 0; i_cell < Npoint_; i_cell++)
            {
                observation_manager_.ApplyOperator(x, innovation);
                Mlt(T(-1), innovation);
                Add(T(1), y, innovation);
                Mlt(R_, innovation, Rinnovation);
                V_(i_cell) += Delta_t_ * DotProd(Rinnovation, innovation);

                int d = Ndimension_ - 1;
                while (d != 0 && position(d) == Nx_(d) - 1)
                {
                    position(d) = 0;
                    x(d) = x_min_(d);
                    d--;
                }
                position(d)++;
                x(d) = x_min_(d) + T(position(d)) * Delta_x_(d);
            }
        }

        time_step_++;

        MessageHandler::Send(*this, "all", "forecast value");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Performs a step forward, using a first-order Lax-Friedrichs scheme.
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>
    ::AdvectionLxFForward()
    {
        MessageHandler::Send(*this, "all", "::AdvectionLxFForward begin");

        // Current value of V.
        Vector<T> V_cur = V_;

        // Position.
        Vector<int> position, position_left, position_right, position_cell;
        // Coordinates.
        Vector<T> x;
        // M(x).
        Vector<T> Mx;

        int i_left, i_right, i_cell;

        /*** Advection term ***/

        Vector<T> time_length_upper_bound(Ndimension_);
        for (int d = 0; d < Ndimension_; d++)
            time_length_upper_bound(d) = Delta_t_ / Delta_x_(d)
                * upper_bound_model_(d);

        // Model times.
        double initial_time = initial_time_ + T(time_step_) * Delta_t_;
        double time, time_step;

        T time_delta = 0.;
        while (time_delta != Delta_t_)
        {
            if (model_time_dependent_)
            {
                courant_number_ = 0.;
                for (int i_cell = 0; i_cell < Npoint_; i_cell++)
                {
                    get_coordinate(i_cell, x_min_, Delta_x_, Nx_, x);

                    model_.SetTime(initial_time);
                    model_.SetState(x);
                    model_.Forward();
                    model_.GetState(Mx);
                    time = model_.GetTime();
                    time_step = time - initial_time;

                    Add(-1., x, Mx);
                    for (int d = 0; d < Ndimension_; d++)
                    {
                        Mx(d) *= Delta_t_ / (Delta_x_(d) * time_step);
                        courant_number_ = max(courant_number_, abs(Mx(d)));
                    }

                    SetRow(Mx, i_cell, Mx_);
                }
            }

            for (int i_cell = 0; i_cell < Npoint_; i_cell++)
            {
                get_position(i_cell, Nx_, position);

                GetRow(Mx_, i_cell, Mx);

                for (int d = 0; d < Ndimension_; d++)
                {
                    if (position(d) == Nx_(d) - 1)
                    {
                        position_left = position;
                        position_left(d)--;
                        i_left = get_position(Nx_, position_left);
                        if (boundary_condition_index_ == 1)
                            boundary_condition_ = 2. * V_cur(i_cell)
                                - V_cur(i_left);
                        else if (boundary_condition_index_ == 2)
                        {
                            position_right = position;
                            position_right(d) = 0;
                            i_right = get_position(Nx_, position_right);
                            boundary_condition_ = V_cur(i_right);
                        }
                        V_(i_cell) += -Mx(d)
                            * .5 * (boundary_condition_ - V_cur(i_left))
                            + time_length_upper_bound(d)
                            * (boundary_condition_ + V_cur(i_left)
                               - 2. * V_cur(i_cell));
                    }
                    else if (position(d) == 0)
                    {
                        position_right = position;
                        position_right(d)++;
                        i_right = get_position(Nx_, position_right);
                        if (boundary_condition_index_ == 1)
                            boundary_condition_ = 2. * V_cur(i_cell)
                                - V_cur(i_right);
                        else if (boundary_condition_index_ == 2)
                        {
                            position_left = position;
                            position_left(d) = Nx_(d) - 1;
                            i_left = get_position(Nx_, position_left);
                            boundary_condition_ = V_cur(i_left);
                        }
                        V_(i_cell) += -Mx(d)
                            * .5 * (V_cur(i_right) - boundary_condition_)
                            + time_length_upper_bound(d)
                            * (V_cur(i_right) + boundary_condition_
                               - 2. * V_cur(i_cell));
                    }
                    else
                    {
                        position_left = position;
                        position_left(d)--;
                        i_left = get_position(Nx_, position_left);
                        position_right = position;
                        position_right(d)++;
                        i_right = get_position(Nx_, position_right);
                        V_(i_cell) += -Mx(d)
                            * .5 * (V_cur(i_right) - V_cur(i_left))
                            + time_length_upper_bound(d)
                            * (V_cur(i_right) + V_cur(i_left)
                               - 2. * V_cur(i_cell));
                    }
                }
            }

            T limit = 0.5;
            if (courant_number_ > limit)
            {
                Logger::Log(*this, "Courant number: "
                            + to_str(courant_number_));
                T division = int(courant_number_ / limit) + 1;
                T local_step = Delta_t_ / division;
                // Checks that Delta_t_ is not reached before.
                if (time_delta + local_step > Delta_t_)
                {
                    local_step = Delta_t_ - time_delta;
                    time_delta = Delta_t_;
                }
                else
                    time_delta += local_step;
                local_step /= Delta_t_;
                for (int i = 0; i < Npoint_; i++)
                    V_(i) = V_cur(i) + (V_(i) - V_cur(i)) * local_step;
                Logger::Log(*this, "Local time step: " + to_str(local_step));
            }
            else
                time_delta = Delta_t_;
        }

        MessageHandler::Send(*this, "all", "::AdvectionLxFForward end");
    }


    /*! \brief Performs a step forward, using a first-order central scheme
      introduced in Bryson and Levy (SIAM J. Sci. Comput., 2003).
    */
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>
    ::AdvectionBrysonLevyForward()
    {
        MessageHandler::Send(*this, "all",
                             "::AdvectionBrysonLevyForward begin");

        // Direction derivatives of V.
        Matrix<T> V_x_m(Npoint_, Ndimension_), V_x_p(Npoint_, Ndimension_);

        // Position.
        Vector<int> position, position_left, position_right, position_cell;
        // Coordinates.
        Vector<T> x;
        // M(x).
        Vector<T> Mx;

        int i_left, i_right, i_cell;

        /*** Advection term ***/

        // Model times.
        double initial_time = initial_time_ + T(time_step_) * Delta_t_;
        double time, time_step;

        if (model_time_dependent_)
        {
            courant_number_ = 0.;
            for (int i_cell = 0; i_cell < Npoint_; i_cell++)
            {
                get_coordinate(i_cell, x_min_, Delta_x_, Nx_, x);

                model_.SetTime(initial_time);
                model_.SetState(x);
                model_.Forward();
                model_.GetState(Mx);
                time = model_.GetTime();
                time_step = time - initial_time;

                Add(-1., x, Mx);
                for (int d = 0; d < Ndimension_; d++)
                {
                    Mx(d) *= Delta_t_ / (Delta_x_(d) * time_step);
                    courant_number_ = max(courant_number_, abs(Mx(d)));
                }

                SetRow(Mx, i_cell, Mx_);
            }
        }

        /*** Computing the directional derivatives of V ***/

        for (int i_cell = 0; i_cell < Npoint_; i_cell++)
        {
            get_position(i_cell, Nx_, position);

            for (int d = 0; d < Ndimension_; d++)
            {
                if (position(d) == Nx_(d) - 1)
                {
                    position_left = position;
                    position_left(d)--;
                    i_left = get_position(Nx_, position_left);
                    if (boundary_condition_index_ == 1)
                        boundary_condition_ = 2. * V_(i_cell) - V_(i_left);
                    else if (boundary_condition_index_ == 2)
                    {
                        position_right = position;
                        position_right(d) = 0;
                        i_right = get_position(Nx_, position_right);
                        boundary_condition_ = V_(i_right);
                    }
                    V_x_m(i_cell, d) = V_(i_cell) - V_(i_left);
                    V_x_p(i_cell, d) = boundary_condition_ - V_(i_cell);
                }
                else if (position(d) == 0)
                {
                    position_right = position;
                    position_right(d)++;
                    i_right = get_position(Nx_, position_right);
                    if (boundary_condition_index_ == 1)
                        boundary_condition_ = 2. * V_(i_cell) - V_(i_right);
                    else if (boundary_condition_index_ == 2)
                    {
                        position_left = position;
                        position_left(d) = Nx_(d) - 1;
                        i_left = get_position(Nx_, position_left);
                        boundary_condition_ = V_(i_left);
                    }
                    V_x_m(i_cell, d) = V_(i_cell) - boundary_condition_;
                    V_x_p(i_cell, d) = V_(i_right) - V_(i_cell);
                }
                else
                {
                    position_left = position;
                    position_left(d)--;
                    i_left = get_position(Nx_, position_left);
                    position_right = position;
                    position_right(d)++;
                    i_right = get_position(Nx_, position_right);
                    V_x_m(i_cell, d) = V_(i_cell) - V_(i_left);
                    V_x_p(i_cell, d) = V_(i_right) - V_(i_cell);
                }
            }
        }

        /*** Evolving the central values ***/

        int i_prev, i_next;
        for (int i_cell = 0; i_cell < Npoint_; i_cell++)
            for (int d = 0; d < Ndimension_; d++)
                V_(i_cell) += .25 * a_Delta_x_(d)
                    * (V_x_p(i_cell, d) - V_x_m(i_cell, d))
                    - .5 * Mx_(i_cell, d)
                    * (V_x_m(i_cell, d) + V_x_p(i_cell, d));

        MessageHandler::Send(*this, "all",
                             "::AdvectionBrysonLevyForward end");
    }


    //! Performs a step forward, using a first-order Godunov scheme.
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>
    ::AdvectionGodunov()
    {
        MessageHandler::Send(*this, "all", "::AdvectionGodunov begin");

        // Current value of V.
        Vector<T> V_cur = V_;

        // Position.
        Vector<int> position, position_left, position_right, position_cell;
        // Coordinates.
        Vector<T> x;
        // M(x).
        Vector<T> Mx(Ndimension_);
        Mx.Zero();

        int i_left, i_right, i_cell;

        /*** Advection term ***/

        // Model times.
        double initial_time = initial_time_ + T(time_step_) * Delta_t_;
        double time, time_step;

        T time_delta = 0.;
        if (with_advection_term_ && model_time_dependent_)
        {
            courant_number_ = 0.;
            for (int i_cell = 0; i_cell < Npoint_; i_cell++)
            {
                get_coordinate(i_cell, x_min_, Delta_x_, Nx_, x);

                model_.SetTime(initial_time);
                model_.SetState(x);
                model_.Forward();
                model_.GetState(Mx);
                time = model_.GetTime();
                time_step = time - initial_time;

                Add(-1., x, Mx);
                for (int d = 0; d < Ndimension_; d++)
                {
                    Mx(d) *= Delta_t_ / (Delta_x_(d) * time_step);
                    courant_number_ = max(courant_number_, abs(Mx(d)));
                }

                SetRow(Mx, i_cell, Mx_);
            }
        }

        for (int i_cell = 0; i_cell < Npoint_; i_cell++)
        {
            get_position(i_cell, Nx_, position);

            if (with_advection_term_)
                GetRow(Mx_, i_cell, Mx);

            for (int d = 0; d < Ndimension_; d++)
            {
                if (position(d) == Nx_(d) - 1)
                {
                    position_left = position;
                    position_left(d)--;
                    i_left = get_position(Nx_, position_left);
                    if (boundary_condition_index_ == 1)
                        boundary_condition_ = 2. * V_cur(i_cell)
                            - V_cur(i_left);
                    else if (boundary_condition_index_ == 2)
                    {
                        position_right = position;
                        position_right(d) = 0;
                        i_right = get_position(Nx_, position_right);
                        boundary_condition_ = V_cur(i_right);
                    }
                    V_(i_cell) -= GodunovFlux(Q_inv_(d, d), Mx(d),
                                              V_cur(i_left), V_cur(i_cell),
                                              boundary_condition_);
                }
                else if (position(d) == 0)
                {
                    position_right = position;
                    position_right(d)++;
                    i_right = get_position(Nx_, position_right);
                    if (boundary_condition_index_ == 1)
                        boundary_condition_ = 2. * V_cur(i_cell)
                            - V_cur(i_right);
                    else if (boundary_condition_index_ == 2)
                    {
                        position_left = position;
                        position_left(d) = Nx_(d) - 1;
                        i_left = get_position(Nx_, position_left);
                        boundary_condition_ = V_cur(i_left);
                    }
                    V_(i_cell) -= GodunovFlux(Q_inv_(d, d), Mx(d),
                                              boundary_condition_,
                                              V_cur(i_cell), V_cur(i_right));
                }
                else
                {
                    position_left = position;
                    position_left(d)--;
                    i_left = get_position(Nx_, position_left);
                    position_right = position;
                    position_right(d)++;
                    i_right = get_position(Nx_, position_right);
                    V_(i_cell) -= GodunovFlux(Q_inv_(d, d), Mx(d),
                                              V_cur(i_left), V_cur(i_cell),
                                              V_cur(i_right));
                }
            }
        }

        MessageHandler::Send(*this, "all", "::AdvectionGodunov end");
    }


    //! Finalizes a step for the model.
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Computes the Godunov flux along a given dimension.
    /*!
      \param[in] q value of Q along the given dimension.
      \param[in] M value of the model M(x) along the given dimension.
      \param[in] v_l value in the left cell.
      \param[in] v value in the center cell.
      \param[in] v_r value in the right cell.
      \return The Godunov flux.
    */
    template <class T, class Model, class ObservationManager>
    inline T HamiltonJacobiBellman<T, Model, ObservationManager>
    ::GodunovFlux(T q, T M, T v_l, T v, T v_r) const
    {
        if (q == T(0))
            if (M < 0.)
                return M * (v_r - v);
            else
                return M * (v - v_l);
        else
        {
            T v_p = v_r - v;
            T v_m = v - v_l;
            if (v_m < v_p)
                if (.5 * q * v_m > -M)
                    return (.25 * q * v_m + M) * v_m;
                else if (.5 * q * v_p < -M)
                    return (.25 * q * v_p + M) * v_p;
                else
                    return -M * M / q;
            else
                return max((.25 * q * v_m + M) * v_m,
                           (.25 * q * v_p + M) * v_p);
        }
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T, class Model, class ObservationManager>
    bool HamiltonJacobiBellman<T, Model, ObservationManager>::HasFinished()
    {
        return time_step_ == Nt_;
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class T, class Model, class ObservationManager>
    Model&
    HamiltonJacobiBellman<T, Model, ObservationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class T, class Model, class ObservationManager>
    ObservationManager&
    HamiltonJacobiBellman<T, Model, ObservationManager>
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
    HamiltonJacobiBellman<T, Model, ObservationManager>
    ::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T, class Model, class ObservationManager>
    string HamiltonJacobiBellman<T, Model, ObservationManager>::GetName()
        const
    {
        return "HamiltonJacobiBellman";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T, class Model, class ObservationManager>
    void HamiltonJacobiBellman<T, Model, ObservationManager>
    ::Message(string message)
    {
        if (message.find("initial value") != string::npos
            || message.find("forecast value") != string::npos)
            output_saver_.Save(V_, time_step_, "value_function");
    }


} // namespace Verdandi.


#define VERDANDI_FILE_HAMILTONJACOBIBELLMAN_CXX
#endif
