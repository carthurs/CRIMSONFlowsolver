// Copyright (C) 2008-2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton, Kevin Charpentier
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


#ifndef VERDANDI_FILE_MODEL_SHALLOWWATER_CXX


#include "ShallowWater.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Constructor.
    template <class T>
    ShallowWater<T>::ShallowWater():
        time_(0.), g_(9.81), urng_(0), current_row_(-1),
        current_column_(-1)
    {
    }


    //! Constructor.
    /*! It builds allocates the state vectors.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    ShallowWater<T>::ShallowWater(string configuration_file):
        time_(0.), g_(9.81), urng_(0), current_row_(-1),
        current_column_(-1)
    {
        Initialize(configuration_file);
    }


    //! Destructor.
    template <class T>
    ShallowWater<T>::~ShallowWater()
    {
        if (urng_ != 0)
        {
            delete urng_;
            urng_ = 0;
        }

        for (int i = 0; i < Nparameter_; i++)
            parameter_[i].Nullify();
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ShallowWater<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/
        VerdandiOps configuration(configuration_file);

        configuration.SetPrefix("shallow_water.domain.");

        configuration.Set("x_min", x_min_);
        configuration.Set("y_min", y_min_);
        configuration.Set("Delta_x", Delta_x_);
        configuration.Set("Delta_y", Delta_y_);
        configuration.Set("Nx", Nx_);
        configuration.Set("Ny", Ny_);

        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("final_time", final_time_);


        // Departure from the uniform initial condition.
        configuration.SetPrefix("shallow_water.initial_condition.");
        configuration.Set("value", value_);

        // Perturbations.
        configuration.SetPrefix("shallow_water.error.");
        configuration.Set("standard_deviation_bc", "v >= 0",
                          model_error_std_bc_);
        configuration.Set("standard_deviation_ic", "v >= 0",
                          model_error_std_ic_);

        if (value_ - 2. * model_error_std_ic_ <= T(0))
            throw "The model standard-deviation of "
                + to_str(model_error_std_ic_)
                + " for initial conditions is too high to avoid negative "
                + "water heights.";

        configuration.Set("random_seed", seed_);

        if (is_num(seed_) && urng_ == 0)
        {
            double seed_number;
            to_num(seed_, seed_number);
            if (seed_number <= 0. || seed_number >= 1.)
                throw "Error: the seed number must be in ]0, 1[.";
            urng_ = new NEWRAN::MotherOfAll(seed_number);
            NEWRAN::Random::Set(*urng_);
        }
        else if (seed_ == "current_time" && urng_ == 0)
        {
            srand(static_cast<unsigned long>(time(0)));
            double seed_number = rand() / double(RAND_MAX);
            urng_ = new NEWRAN::MotherOfAll(seed_number);
            NEWRAN::Random::Set(*urng_);
        }
        else if (urng_ == 0)
        {
            throw "In this version, the seed cannot be read from disk.";
            NEWRAN::Random::SetDirectory(seed_.c_str());
            urng_ = new NEWRAN::MotherOfAll;
            NEWRAN::Random::Set(*urng_);
            NEWRAN::Random::CopySeedFromDisk(true);
        }

        // Error statistics.
        configuration.SetPrefix("shallow_water.state_error.");

        configuration.Set("variance", "v >= 0", state_error_variance_value_);
        configuration.Set("scale", "v > 0", Balgovind_scale_background_);
#ifdef VERDANDI_STATE_ERROR_SPARSE
        build_diagonal_sparse_matrix(Nx_ * Ny_,
                                     state_error_variance_value_,
                                     state_error_variance_);
#endif

        // Description of boundary conditions.
        ReadConfigurationBoundaryCondition("left", configuration,
                                           boundary_condition_left_,
                                           value_left_, amplitude_left_,
                                           frequency_left_);
        if (boundary_condition_left_ == 2)
            // The opposite of the in-flow rate is used when computing the
            // values in the ghost cell on the left.
        {
            value_left_ *= -1;
            amplitude_left_ *= -1;
        }
        ReadConfigurationBoundaryCondition("right", configuration,
                                           boundary_condition_right_,
                                           value_right_, amplitude_right_,
                                           frequency_right_);
        ReadConfigurationBoundaryCondition("bottom", configuration,
                                           boundary_condition_bottom_,
                                           value_bottom_, amplitude_bottom_,
                                           frequency_bottom_);
        if (boundary_condition_bottom_ == 2)
            // The opposite of the in-flow rate is used when computing the
            // values in the ghost cell at the bottom.
        {
            value_bottom_ *= -1;
            amplitude_bottom_ *= -1;
        }
        ReadConfigurationBoundaryCondition("top", configuration,
                                           boundary_condition_top_,
                                           value_top_, amplitude_top_,
                                           frequency_top_);

        configuration.SetPrefix("shallow_water.data_assimilation.");

        configuration.Set("with_positivity_requirement",
                          with_positivity_requirement_);


        /*** Allocations ***/

        h_.Reallocate(Nx_, Ny_);
        u_.Reallocate(Nx_, Ny_);
        v_.Reallocate(Nx_, Ny_);

        // Initial conditions.
        h_.Fill(1.);
        NEWRAN::Random::Set(*urng_);
        value_ += max(-2., min(2., normal_.Next())) * model_error_std_ic_;
        configuration.SetPrefix("shallow_water.initial_condition.");
        configuration.Set("center", source_center_);

        if (source_center_)
        {
            int center_x = (Nx_ - 1) / 2;
            int center_y = (Ny_ - 1) / 2;
            h_(center_x, center_y) = value_;
            for (int i = -5; i < 6; i++)
                for (int j = -5; j < 6; j++)
                    if (center_x + i >= 0 && center_x + i < Nx_
                        && center_y + j >= 0 && center_y + j < Ny_)
                        h_(center_x + i, center_y + j) = value_;
        }
        u_.Fill(0.);
        v_.Fill(0.);

        configuration.Set("left", source_left_);
        if (source_left_)
        {
            int position_x = Nx_ / 10;
            int center_y = (Ny_ - 1) / 2;
            for (int i = -2; i < 3; i++)
                for (int j = -2; j < 3; j++)
                    if (position_x + i >= 0 && position_x + i < Nx_
                        && center_y + j >= 0 && center_y + j < Ny_)
                        h_(position_x + i, center_y + j) = value_;
        }

        hf_x_.Reallocate(Nx_ + 1, Ny_);
        hf_y_.Reallocate(Nx_, Ny_ + 1);
        uf_x_.Reallocate(Nx_ + 1, Ny_);
        uf_y_.Reallocate(Nx_, Ny_ + 1);
        vf_x_.Reallocate(Nx_ + 1, Ny_);
        vf_y_.Reallocate(Nx_, Ny_ + 1);

        state_error_variance_row_.Reallocate(Nx_ * Ny_);

        // Perturbations of step height and bc.
        configuration.ClearPrefix();
        Nparameter_ = 0;

        if (configuration.Exists("shallow_water.uncertainty"))
        {
            configuration.SetPrefix("shallow_water.uncertainty.");

            configuration.Set("uncertain_parameter_list",
                              "ops_in(v, {'step_height', 'bc'})",
                              uncertain_parameter_vector_);

            Nparameter_ = int(uncertain_parameter_vector_.size());
            if (Nparameter_ == 2)
            {
                // "step_height" must be before "bc" for later calls.
                uncertain_parameter_vector_[0] = "step_height";
                uncertain_parameter_vector_[1] = "bc";
            }
            parameter_.resize(Nparameter_);

            if (find(uncertain_parameter_vector_.begin(),
                     uncertain_parameter_vector_.end(),
                     "step_height") != uncertain_parameter_vector_.end())
            {
                configuration.Set("step_height.mean", step_height_mean_);
                value_ += step_height_mean_;

                step_height_variance_.Reallocate(1, 1);
                configuration.Set("step_height.variance",
                                  step_height_variance_(0, 0));

                configuration.Set("step_height.distribution",
                                  "ops_in(v, {'Normal', 'LogNormal', "
                                  "'NormalHomogeneous', "
                                  "'LogNormalHomogeneous'})",
                                  step_height_pdf_);

                configuration.Set("step_height.parameter",
                                  step_height_parameter_);

                parameter_[0].SetData(1, &value_);
            }

            if (find(uncertain_parameter_vector_.begin(),
                     uncertain_parameter_vector_.end(),
                     "bc") != uncertain_parameter_vector_.end())
            {
                configuration.Set("bc.mean", bc_mean_);
                value_left_ += bc_mean_;

                bc_variance_.Reallocate(1, 1);
                configuration.Set("bc.variance", bc_variance_(0, 0));

                configuration.Set("bc.distribution",
                                  "ops_in(v, {'Normal', 'LogNormal', "
                                  "'NormalHomogeneous', "
                                  "'LogNormalHomogeneous'})",
                                  bc_pdf_);

                configuration.Set("bc.parameter", bc_parameter_);

                if (Nparameter_ == 1)
                    parameter_[0].SetData(1, &value_left_);
                else
                    parameter_[1].SetData(1, &value_left_);
            }
        }

        /*** Output saver ***/

        output_saver_.
            Initialize(configuration_file, "shallow_water.output_saver.");
        output_saver_.Empty("u");
        output_saver_.Empty("v");
        output_saver_.Empty("h");
    }


    //! Initializes the current time step for the model.
    template <class T>
    void ShallowWater<T>::InitializeStep()
    {

    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    template <class T>
    void ShallowWater<T>::Forward()
    {
        int i, j;

        if (time_ == 0. && source_center_)
        {
            int center_x = (Nx_ - 1) / 2;
            int center_y = (Ny_ - 1) / 2;
            h_(center_x, center_y) = value_;
            for (int i = -5; i < 6; i++)
                for (int j = -5; j < 6; j++)
                    if (center_x + i >= 0 && center_x + i < Nx_
                        && center_y + j >= 0 && center_y + j < Ny_)
                        h_(center_x + i, center_y + j) = value_;
        }
        if (time_ == 0. && source_left_)
        {
            int position_x = Nx_ / 10;
            int center_y = (Ny_ - 1) / 2;
            for (int i = -2; i < 3; i++)
                for (int j = -2; j < 3; j++)
                    if (position_x + i >= 0 && position_x + i < Nx_
                        && center_y + j >= 0 && center_y + j < Ny_)
                        h_(position_x + i, center_y + j) = value_;
        }

        // Left and right values around the interface where a flux is
        // computed. Heights.
        T h_ghost;
        // Velocities.
        T u_ghost, v_ghost;

        /*** Computing the fluxes ***/

        // Model error related to boundary conditions.
        double model_error = 0;

        // Clipped Gaussian.
        if (model_error_std_bc_ != 0.)
        {
            NEWRAN::Random::Set(*urng_);
            model_error
                = max(-2., min(2., normal_.Next())) * model_error_std_bc_;

        }

        // Fluxes along x, inside the domain.
        for (i = 0; i < Nx_ - 1; i++)
            for (j = 0; j < Ny_; j++)
            {
                ComputeFluxHLL(h_(i, j), h_(i + 1, j),
                               u_(i, j), u_(i + 1, j),
                               v_(i, j), v_(i + 1, j),
                               hf_x_(i + 1, j), uf_x_(i + 1, j),
                               vf_x_(i + 1, j));
            }

        // Fluxes along x, boundary conditions.
        for (j = 0; j < Ny_; j++)
        {
            ComputeGhostCellValue(boundary_condition_left_,
                                  value_left_ + model_error,
                                  amplitude_left_, frequency_left_,
                                  h_(0, j), -u_(0, j), -v_(0, j),
                                  h_ghost, u_ghost, v_ghost);
            ComputeFluxHLL(h_ghost, h_(0, j),
                           -u_ghost, u_(0, j),
                           -v_ghost, v_(0, j),
                           hf_x_(0, j), uf_x_(0, j), vf_x_(0, j));
            ComputeGhostCellValue(boundary_condition_right_,
                                  value_right_ + model_error,
                                  amplitude_right_, frequency_right_,
                                  h_(Nx_ - 1, j), u_(Nx_ - 1, j),
                                  v_(Nx_ - 1, j),
                                  h_ghost, u_ghost, v_ghost);
            ComputeFluxHLL(h_(Nx_ - 1, j), h_ghost,
                           u_(Nx_ - 1, j), u_ghost,
                           v_(Nx_ - 1, j), v_ghost,
                           hf_x_(Nx_, j), uf_x_(Nx_, j), vf_x_(Nx_, j));
        }

        // Fluxes along y, inside the domain.
        for (i = 0; i < Nx_; i++)
            for (j = 0; j < Ny_ - 1; j++)
                ComputeFluxHLL(h_(i, j), h_(i, j + 1),
                               v_(i, j), v_(i, j + 1),
                               -u_(i, j), -u_(i, j + 1),
                               hf_y_(i, j + 1), vf_y_(i, j + 1),
                               uf_y_(i, j + 1));

        // Fluxes along y, boundary conditions.
        for (i = 0; i < Nx_; i++)
        {
            ComputeGhostCellValue(boundary_condition_bottom_,
                                  value_bottom_ + model_error,
                                  amplitude_bottom_, frequency_bottom_,
                                  h_(i, 0), -v_(i, 0), u_(i, 0),
                                  h_ghost, v_ghost, u_ghost);
            ComputeFluxHLL(h_ghost, h_(i, 0),
                           -v_ghost, v_(i, 0),
                           u_ghost, u_(i, 0),
                           hf_y_(i, 0), vf_y_(i, 0), uf_y_(i, 0));
            ComputeGhostCellValue(boundary_condition_top_,
                                  value_top_ + model_error,
                                  amplitude_top_, frequency_top_,
                                  h_(i, Ny_ - 1), v_(i, Ny_ - 1),
                                  -u_(i, Ny_ - 1),
                                  h_ghost, v_ghost, u_ghost);
            ComputeFluxHLL(h_(i, Ny_ - 1), h_ghost,
                           v_(i, Ny_ - 1), v_ghost,
                           -u_(i, Ny_ - 1), u_ghost,
                           hf_y_(i, Ny_), vf_y_(i, Ny_), uf_y_(i, Ny_));
        }

        /*** Updating the state ***/

        double factor = Delta_x_ * Delta_y_ * Delta_t_;
        for (i = 0; i < Nx_; i++)
            for (j = 0; j < Ny_; j++)
            {
                h_(i, j) += factor * (hf_x_(i, j) - hf_x_(i + 1, j)
                                      + hf_y_(i, j) - hf_y_(i, j + 1));
                u_(i, j) += factor * (uf_x_(i, j) - uf_x_(i + 1, j)
                                      - uf_y_(i, j) + uf_y_(i, j + 1));
                v_(i, j) += factor * (vf_x_(i, j) - vf_x_(i + 1, j)
                                      + vf_y_(i, j) - vf_y_(i, j + 1));
            }

        // Checking the CFL.
        T velocity = 0.;
        T tmp;
        for (i = 0; i < Nx_; i++)
            for (j = 0; j < Ny_; j++)
            {
                tmp = sqrt(u_(i, j) * u_(i, j) + v_(i, j) * v_(i, j))
                    + sqrt(g_ * h_(i, j));
                velocity = max(velocity, tmp);
            }
        if (2. * Delta_t_ * velocity / Delta_x_ > 1.)
        {
            cout << "Error! The Courant-Friedrichs-Lewy condition is not met."
                 << endl << "  Acceptable velocity: "
                 << 2. * Delta_t_ * velocity / Delta_x_ << endl
                 << "  Actual velocity:     " <<  velocity << endl;
            abort();
        }

        time_ += Delta_t_;
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool ShallowWater<T>::HasFinished() const
    {
        return time_ >= final_time_;
    }


    //! Saves the simulated data.
    /*! It saves the velocities 'u' and 'v', and the height 'h'.
     */
    template <class T>
    void ShallowWater<T>::Save()
    {
        output_saver_.Save(u_, time_, "u");
        output_saver_.Save(v_, time_, "v");
        output_saver_.Save(h_, time_, "h");
    }


    //! Finalizes the current time step for the model.
    template <class T>
    void ShallowWater<T>::FinalizeStep()
    {
    }


    //! Finalizes the model.
    template <class T>
    void ShallowWater<T>::Finalize()
    {
    }


    ///////////////////
    // ACCESS METHOD //
    ///////////////////


    //! Returns the current time.
    /*!
      \return The current time.
    */
    template <class T>
    double ShallowWater<T>::GetTime() const
    {
        return time_;
    }


    //! Sets the current time.
    /*!
      \param[in] time the current time.
    */
    template <class T>
    void ShallowWater<T>::SetTime(double time)
    {
        time_ = time;
    }


    //! Returns the number of points along x (in the grid for height).
    /*!
      \return The number of points along x (in the grid for height).
    */
    template <class T>
    int ShallowWater<T>::GetNx() const
    {
        return Nx_;
    }


    //! Returns the number of points along y (in the grid for height).
    /*!
      \return The number of points along y (in the grid for height).
    */
    template <class T>
    int ShallowWater<T>::GetNy() const
    {
        return Ny_;
    }


    //! Returns the number of points along z.
    /*!
      \return The number of points along z.
    */
    template <class T>
    int ShallowWater<T>::GetNz() const
    {
        return 1;
    }


    //! Returns the first abscissa.
    /*!
      \return The first abscissa.
    */
    template <class T>
    int ShallowWater<T>::GetXMin() const
    {
        return x_min_;
    }


    //! Returns the first ordinate.
    /*!
      \return The first ordinate.
    */
    template <class T>
    int ShallowWater<T>::GetYMin() const
    {
        return y_min_;
    }


    //! Returns the space step along x.
    /*!
      \return The space step along x.
    */
    template <class T>
    int ShallowWater<T>::GetDeltaX() const
    {
        return Delta_x_;
    }


    //! Returns the space step along y.
    /*!
      \return The space step along y.
    */
    template <class T>
    int ShallowWater<T>::GetDeltaY() const
    {
        return Delta_y_;
    }


    //! Returns the number of points in the grid.
    /*!
      \return The number of points in the grid.
    */
    template <class T>
    int ShallowWater<T>::GetNstate() const
    {
        return Nx_ * Ny_;
    }


    //! Returns the full state vector size.
    /*!
      \return The full state vector size.
    */
    template <class T>
    int ShallowWater<T>::GetNfull_state() const
    {
        return 2 * Nx_ * Ny_ + 2 * Nx_ + Ny_;
    }


    //! Provides the state vector.
    /*!
      \return A reference to the state vector
    */
    template <class T>
    typename ShallowWater<T>::state& ShallowWater<T>
    ::GetState()
    {
        state_.Reallocate(Nx_ * Ny_);
        int position = 0;
        for (int i = 0; i < Nx_; i++)
            for (int j = 0; j < Ny_; j++)
                state_(position++) = h_(i, j);
        return state_;
    }


    //! Performs some calculations when the update of the model state is done.
    template <class T>
    void ShallowWater<T>::StateUpdated()
    {
        // Positivity requirement.
        if (with_positivity_requirement_)
            for (int r = 0; r < Nx_ * Ny_; r++)
                if (state_(r) < T(0))
                    state_(r) = T(0);

        int position = 0;
        for (int i = 0; i < Nx_; i++)
            for (int j = 0; j < Ny_; j++)
                h_(i, j) = state_(position++);
    }


    //! Provides the full state vector.
    /*!
      \param[out] state the full state vector.
    */
    template <class T>
    typename ShallowWater<T>::state& ShallowWater<T>
    ::GetFullState()
    {
        state_.Reallocate(2 * Nx_ * Ny_ + 2 * Nx_ + Ny_);
        for (int i = 0; i < Nx_; i++)
            for (int j = 0; j < Ny_; j++)
            {
                state_(i * Ny_ + j) = h_(i, j);
                state_(Nx_ * Ny_ + i * Ny_ + j) = u_(i, j);
                state_(2 * Nx_ * Ny_ + i * Ny_ + j) = v_(i, j);
            }
        return state_;
    }


    //! Performs some calculations when the update of the model state is done.
    template <class T>
    void ShallowWater<T>
    ::FullStateUpdated()
    {
        for (int i = 0; i < Nx_; i++)
            for (int j = 0; j < Ny_; j++)
            {
                h_(i, j) = state_(i * Ny_ + j);
                u_(i, j) = state_(Nx_ * Ny_ + i * Ny_ + j);
                v_(i, j) = state_(2 * Nx_ * Ny_ + i * Ny_ + j);
            }
    }


    //! Computes a row of the background error covariance matrix B.
    /*!
      \param[in] row row index.
      \param[out] state_error_variance_row the value of row number \a row.
    */
    template <class T>
    void ShallowWater<T>
    ::GetStateErrorVarianceRow(int row, state_error_variance_row&
                               state_error_variance_row)
    {
#ifdef VERDANDI_STATE_ERROR_SPARSE
        {
            state_error_variance_row.Reallocate(Nx_ * Ny_);
            state_error_variance_row.Zero();
            state_error_variance_row(row) = state_error_variance_value_;
        }
# else
        {
#ifdef VERDANDI_STATE_ERROR_DENSE
            {
                state_error_variance_row.Reallocate(Nx_ * Ny_);
                state_error_variance_row.Zero();
                state_error_variance_row(row)
                    = state_error_variance_value_;
            }
#else
            {
                // The row has already been computed.
                if (row == current_row_)
                    state_error_variance_row = state_error_variance_row_;
                else
                {
                    int i, j;
                    current_row_ = row;
                    current_column_ = -1;

                    // Positions related to 'row'.
                    int i_row = row / Ny_;
                    int j_row = row - i_row * Ny_;

                    T distance_x, distance_y;
                    T distance;
                    int position = 0;
                    for (i = 0; i < Nx_; i++)
                        for (j = 0; j < Ny_; j++)
                        {
                            distance_x = Delta_x_ * T(i - i_row);
                            distance_y = Delta_y_ * T(j - j_row);
                            distance = sqrt(distance_x * distance_x
                                            + distance_y * distance_y)
                                / Balgovind_scale_background_;
                            state_error_variance_row_(position++)
                                = state_error_variance_value_
                                * (1. + distance) * exp(-distance);
                        }
                    state_error_variance_row = state_error_variance_row_;
                }
            }
#endif
        }
#endif
    }


    //! Returns the background error covariance matrix (B) if available.
    /*! Returns the background error covariance matrix (B) if available,
      raises an exception otherwise.
      \return The matrix of the background error covariance.
    */
    template <class T>
    const typename ShallowWater<T>::state_error_variance& ShallowWater<T>
    ::GetStateErrorVariance() const
    {
#ifdef VERDANDI_STATE_ERROR_SPARSE
        return state_error_variance_;
#else
        throw ErrorUndefined(
            "ShallowWater::GetStateErrorVariance()",
            "The state error covariance matrix is not available!");
#endif
    }


    //! Checks if the error covariance matrix is sparse.
    /*!
      \return True if there is a sparse error matrix, false otherwise.
    */
    template <class T>
    bool ShallowWater<T>::IsErrorSparse() const
    {
#ifdef VERDANDI_STATE_ERROR_SPARSE
        return true;
#else
        return false;
#endif
    }


    //! Returns the number of parameters to be perturbed.
    /*!
      \return The number of parameters to be perturbed.
    */
    template <class T>
    int ShallowWater<T>::GetNparameter()
    {
        return Nparameter_;
    }


    //! Gets the i-th uncertain parameter.
    /*!
      \param[in] i index of the parameter.
      \return The vector associated with the i-th parameter.
    */
    template<class T>
    typename ShallowWater<T>::uncertain_parameter&
    ShallowWater<T>::GetParameter(int k)
    {
        return parameter_[k];
    }


    //! Sets the i-th parameter.
    /*!
      \param[in] i index of the parameter.
      \param[in] parameter the parameter to assign.
    */
    template<class T>
    void ShallowWater<T>::SetParameter(int i, uncertain_parameter parameter)
    {
        if (i == 0 && uncertain_parameter_vector_[0] == "step_height")
            value_ = parameter(0);
        else
            value_left_ = parameter(0);
    }


    //! Returns the correlation between the uncertain parameters.
    /*! Since there is only one parameter, an empty vector is
      returned.
      \param[in] i parameter index.
      \return An empty vector.
    */
    template<class T>
    Vector<T>& ShallowWater<T>::GetParameterCorrelation(int i)
    {
        if (i == 0 && uncertain_parameter_vector_[0] == "step_height")
            return step_height_correlation_;
        else
            return bc_correlation_;
    }


    //! Returns the PDF of the i-th parameter.
    /*!
      \param[in] i uncertain-variable index.
      \return The PDF of the i-th parameter.
    */
    template<class T>
    string ShallowWater<T>::GetParameterPDF(int i)
    {
        if (i == 0 && uncertain_parameter_vector_[0] == "step_height")
            return step_height_pdf_;
        else
            return bc_pdf_;
    }


    /*! \brief Returns the covariance matrix associated with the i-th
      parameter. */
    /*!
      \param[in] i parameter index.
      \return The covariance matrix associated with the i-th parameter.
    */
    template<class T>
    Matrix<T, Symmetric, RowSymPacked>&
    ShallowWater<T>::GetParameterVariance(int i)
    {
        if (i == 0 && uncertain_parameter_vector_[0] == "step_height")
            return step_height_variance_;
        else
            return bc_variance_;
    }



    //! Returns parameters associated with the PDF of some model parameter.
    /*! In case of normal or log-normal distribution, the parameters are
      clipping parameters.
      \param[in] i model parameter index.
      \return The parameters associated with the i-th parameter.
    */
    template<class T>
    Vector<T>& ShallowWater<T>::GetParameterParameter(int i)
    {
        if (i == 0 && uncertain_parameter_vector_[0] == "step_height")
            return step_height_parameter_;
        else
            return bc_parameter_;
    }



    //! Returns the perturbation option of the i-th parameter.
    /*!
      \param[in] i parameter index.
      \return The perturbation option of the i-th parameter.
    */
    template<class T>
    string ShallowWater<T>::GetParameterOption(int i)
    {
        if (i == 0 && uncertain_parameter_vector_[0] == "step_height")
            return "init_step";
        else
            return "every_step";
    }


    //! Returns the model itself.
    /*!
      \return The model.
    */
    template <class T>
    const ShallowWater<T>& ShallowWater<T>::GetModel() const
    {
        return *this;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string ShallowWater<T>::GetName() const
    {
        return "ShallowWater";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void ShallowWater<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


    ////////////////////////////
    // PRIVATE USEFUL METHODS //
    ////////////////////////////


    /*** Configuration ***/

    //! Reads the configuration of a boundary condition.
    /*! It accepts the following boundary conditions: (0) "free" for
      homogeneous Neumann on all height and rate, (1) "wall" for impermeable
      conditions, (2) "flow" for in-flow or out-flow conditions, (3) "height"
      for Dirichlet conditions on the height.
      \param[in] side side for the boundary condition (left, right, bottom or
      top).
      \param[in] configuration VerdandiOps instance.
      \param[out] type type of the boundary condition (0, 1, 2 or 3).
      \param[out] inflow in-flow rate.
      \param[out] amplitude amplitude of the variations.
      \param[out] frequency angular frequency of the variations.
    */
    template <class T>
    void ShallowWater<T>
    ::ReadConfigurationBoundaryCondition(string side,
                                         VerdandiOps& configuration,
                                         int& type, T& value,
                                         T& amplitude, T& frequency)
    {
        configuration.SetPrefix("shallow_water.boundary_condition.");
        string description;
        configuration.Set(side, description);
        if (description == "free")
        {
            type = 0;
            return;
        }
        if (description == "wall")
        {
            type = 1;
            return;
        }

        vector<string> description_vector = split(description);
        if ((description_vector.size() != 2 && description_vector.size() != 4)
            || (description_vector[0] != "flow"
                && description_vector[0] != "height"))
            throw "Error in the description of " + side
                + " boundary conditions. Got:\n" + description + "\n"
                + "Six options are available:\n"
                + "free\nwall\nflow value\nflow value amplitude frequency\n"
                + "height value\nheight value amplitude frequency";

        if (description_vector[0] == "flow")
            type = 2;
        else // "height"
            type = 3;
        to_num(description_vector[1], value);
        if (type == 3 && value <= T(0))
            throw
                string("Any fixed height given for boundary condition ")
                + "should be strictly greater than 0, but " + to_str(value)
                + " was given.";
        if (description_vector.size() == 4)
        {
            to_num(description_vector[2], amplitude);
            // Checks that the height cannot be negative.
            if (type == 3 && value - fabs(amplitude) <= T(0))
                throw
                    string("Any height given for boundary condition should ")
                    + "be strictly greater than 0, but a mean value of "
                    + to_str(value)
                    + " and an amplitude of " + to_str(amplitude)
                    + " were given.";
            // Checks that the model error cannot be larger than the height.
            if (type == 3
                && value - fabs(amplitude) - 2. * model_error_std_bc_ <= T(0))
                throw
                    string("Any height given for boundary condition should ")
                    + "be strictly greater than 0, but a mean value of "
                    + to_str(value)
                    + ", an amplitude of " + to_str(amplitude)
                    + " and a model standard-deviation of "
                    + to_str(model_error_std_bc_)
                    + " were given.";
            to_num(description_vector[3], frequency);
            // Converts the frequency to an angular frequency.
            frequency *= 6.2831853071795862;
        }
        else
        {
            amplitude = T(0);
            frequency = T(0);
        }
    }

    /*** Processing ***/

    /*! \brief Determines the surrounding fictitious values to manage the
      boundary conditions. */
    template <class T>
    void ShallowWater<T>
    ::ComputeGhostCellValue(int type, T value, T amplitude, T frequency,
                            T h_l, T u_l, T v_l, T& h_r, T& u_r, T& v_r)
    {
        if (type == 0) // free
        {
            h_r = h_l;
            u_r = u_l;
            v_r = v_l;
        }
        if (type == 1) // wall
        {
            h_r = h_l;
            u_r = -u_l;
            v_r = v_l;
        }
        if (type == 2) // flow
        {
            h_r = h_l;
            u_r = value / h_l;
            v_r = v_l;
        }
        if (type == 3) // height
        {
            h_r = value;
            u_r = u_l + 2. * sqrt(g_) * (sqrt(h_l) - sqrt(h_r));
            v_r = v_l;
        }
    }


    //! Computes the HHL flux.
    template <class T>
    void ShallowWater<T>
    ::ComputeFluxHLL(T h_l, T h_r, T u_l, T u_r, T v_l, T v_r,
                     T& flux_h, T& flux_u, T& flux_v)
    {
        T tmp;

        // Wave velocity.
        T s_l, s_r;
        // Phase speed of the wave.
        T c_l, c_r, c_lr;
        // Other coefficients.
        T p_l, p_r;

        // Height at the interface.
        T h_interface;
        T h_tmp;

        T u_lr;

        c_l = sqrt(g_ * h_l);
        c_r = sqrt(g_ * h_r);
        c_lr = c_l + c_r;

        u_lr = u_l - u_r;

        h_tmp = .5 * (h_l + h_r) * (1. + .5 * u_lr / c_lr);

        if (h_tmp <= min(h_l, h_r))
        {
            tmp = .5 * c_lr + .25 * u_lr;
            h_interface = tmp * tmp / g_;
        }
        else if (h_tmp >= max(h_l, h_r))
        {
            T g_l, g_r;
            g_l = sqrt(.5 * g_ * (h_tmp + h_l) / (h_tmp * h_l));
            g_r = sqrt(.5 * g_ * (h_tmp + h_r) / (h_tmp * h_r));
            h_interface = (h_l * g_l + h_r * g_r + u_lr) / (g_l + g_r);
        }
        else
            h_interface = h_tmp;

        p_l = h_interface > h_l ?
            sqrt(.5 * h_interface * (h_interface + h_l)) / h_l
            : T(1.);
        s_l = u_l - c_l * p_l;

        if (s_l >= 0.)
        {
            ComputeFlux(h_l, u_l, v_l, flux_h, flux_u, flux_v);
            return;
        }

        p_r = h_interface > h_r ?
            sqrt(.5 * h_interface * (h_interface + h_r)) / h_r
            : T(1.);
        s_r = u_r + c_r * p_r;

        if (s_r <= 0.)
        {
            ComputeFlux(h_r, u_r, v_r, flux_h, flux_u, flux_v);
            return;
        }

        // Below, s_l < 0. and s_r > 0..
        T flux_h_l, flux_u_l, flux_v_l, flux_h_r, flux_u_r, flux_v_r;
        ComputeFlux(h_l, u_l, v_l, flux_h_l, flux_u_l, flux_v_l);
        ComputeFlux(h_r, u_r, v_r, flux_h_r, flux_u_r, flux_v_r);
        tmp = 1. / (s_r - s_l);
        flux_h
            = (s_r * flux_h_l - s_l * flux_h_r + s_r * s_l * (h_r - h_l))
            * tmp;
        flux_u
            = (s_r * flux_u_l - s_l * flux_u_r - s_r * s_l * u_lr) * tmp;
        flux_v
            = (s_r * flux_v_l - s_l * flux_v_r + s_r * s_l * (v_r - v_l))
            * tmp;
    }


    //! Computes the flux.
    template <class T>
    void ShallowWater<T>
    ::ComputeFlux(T h, T u, T v, T& flux_h, T& flux_u, T& flux_v)
    {
        T tmp = h * u;
        flux_h = tmp;
        flux_u = tmp * u + .5 * g_ * h * h;
        flux_v = tmp * v;
    }


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_SHALLOWWATER_CXX
#endif
