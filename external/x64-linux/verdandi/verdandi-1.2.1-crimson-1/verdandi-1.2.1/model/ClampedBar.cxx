// Copyright (C) 2008-2009 INRIA
// Author(s): Dominique Chapelle, Philippe Moireau, Marc Fragu,
//            Akos Matszangosz
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


#ifndef VERDANDI_FILE_MODEL_CLAMPEDBAR_CXX


#include "ClampedBar.hxx"

#include "OutputSaver.cxx"

#include "seldon/vector/VectorCollection.cxx"

namespace Verdandi
{


    ///////////////////
    // STATIC FIELDS //
    ///////////////////


    template <class T>
    const double ClampedBar<T>::Pi_ = 3.141592653589793238462;


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Constructor.
    template <class T>
    ClampedBar<T>::ClampedBar(): time_(0.)
    {
    }


    //! Constructor.
    /*! It builds allocates the state vectors.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    ClampedBar<T>::ClampedBar(string configuration_file):
        time_(0.)
    {
    }


    //! Destructor.
    template <class T>
    ClampedBar<T>::~ClampedBar()
    {
        x_.Nullify();
        x_full_.Nullify();
        // Deallocate q.
        if (is_adjoint_initialized_)
        {
            state working_vector;
            working_vector.SetData(q_.GetM(), q_.GetVector(0).GetData());
            q_.Nullify();
        }
    }


    ////////////////
    // INITIALIZE //
    ////////////////


    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void ClampedBar<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        VerdandiOps configuration(configuration_file);

        configuration.SetPrefix("clamped_bar.domain.");
        configuration.Set("bar_length", bar_length_);
        configuration.Set("Nx", Nx_);
        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("final_time", final_time_);
        configuration.SetPrefix("clamped_bar.error_statistics.");
        configuration.Set("state_error_variance", "v >= 0",
                          state_error_variance_value_);
        configuration.SetPrefix("clamped_bar.physics.");
        configuration.Set("Young_modulus", Young_modulus_);
        configuration.Set("mass_density", mass_density_);

        configuration.Set("theta_force", theta_force_);
        Ntheta_force_ = theta_force_.GetSize();
        Ndof_ = Nx_ + 1;
        BuildRegionIndex(Ndof_, Ntheta_force_, theta_force_index_);
        configuration.Set("theta_stiffness", theta_stiffness_);
        Ntheta_stiffness_ = theta_stiffness_.GetSize();
        BuildRegionIndex(Ndof_, Ntheta_stiffness_, theta_stiffness_index_);
        configuration.Set("theta_mass", theta_mass_);
        Ntheta_mass_ = theta_mass_.GetSize();
        BuildRegionIndex(Ndof_, Ntheta_mass_, theta_mass_index_);
        configuration.Set("theta_damp", theta_damp_);
        Ntheta_damp_ = theta_damp_.GetSize();
        BuildRegionIndex(Ndof_, Ntheta_damp_, theta_damp_index_);

        configuration.Set("alpha", alpha_);
        configuration.Set("beta", beta_);
        configuration.SetPrefix("clamped_bar.");
        vector<string> stable;
        configuration.Set("state","ops_in(v, {'displacement', 'velocity', "
                          "'theta_force', 'theta_stiffness', 'theta_mass', "
                          "'theta_damp'})",
                          stable);
        if (stable.size() <= 0)
            throw ErrorArgument("void ClampedBar<T>::Initialize"
                                "(string configuration_file)", "The name of"
                                " the underlying state vector (state) are "
                                "not defined.");
        for (unsigned int i = 0; i < stable.size(); i++)
            stable_.insert(stable[i]);

        configuration.Set("reduced_state", reduced_);

        /*** Ouput saver ***/

        output_saver_.Initialize(configuration_file,
                                 "clamped_bar.output_saver.");
        output_saver_.Empty("disp_0");
        output_saver_.Empty("velo_0");

        /*** Allocation ***/

        Delta_x_ = bar_length_ / Nx_;
        AllocateSparseMatrix();

        /*** Elementary mass matrix construction ***/

        mass_FEM_matrix_.Reallocate(2, 2);
        T mass_lin = T(mass_density_ * Delta_x_ / 3);
        mass_FEM_matrix_(0, 0) = mass_lin;
        mass_FEM_matrix_(1, 1) = mass_lin;
        mass_FEM_matrix_(0, 1) = mass_lin / 2;

        /*** Elementary stifness matrix construction ***/

        stiffness_FEM_matrix_.Reallocate(2, 2);
        T stiff_lin = T(Young_modulus_ / Delta_x_);
        stiffness_FEM_matrix_(0, 0) = stiff_lin;
        stiffness_FEM_matrix_(1, 1) = stiff_lin;
        stiffness_FEM_matrix_(0, 1) = -stiff_lin;

        /*** Elementary damping matrix construction ***/

        damp_FEM_matrix_.Reallocate(2, 2);
        damp_FEM_matrix_(0, 0) = alpha_ * mass_FEM_matrix_(0, 0) +
            beta_ * stiffness_FEM_matrix_(0, 0);
        damp_FEM_matrix_(1, 1) = alpha_ * mass_FEM_matrix_(1, 1) +
            beta_ * stiffness_FEM_matrix_(1, 1);
        damp_FEM_matrix_(0, 1) = alpha_ * mass_FEM_matrix_(0, 1) +
            beta_ * stiffness_FEM_matrix_(0, 1);

        /*** State initialization ***/

        disp_0_.Reallocate(Ndof_);
        velo_0_.Reallocate(Ndof_);
        force_.Reallocate(Ndof_);
        disp_0_.Fill(T(0));
        velo_0_.Fill(T(0));
        force_.Fill(T(1));

        for (int i = 0; i < Ndof_; i++)
            disp_0_(i) = T(0);

        // Initialize full vector collection.
        state working_vector;
        working_vector.SetData(Ndof_ - 1, disp_0_.GetData() + 1);
        x_full_.AddVector(working_vector, "displacement");
        working_vector.Nullify();
        working_vector.SetData(Ndof_ - 1, velo_0_.GetData() + 1);
        x_full_.AddVector(working_vector, "velocity");
        working_vector.Nullify();
        x_full_.AddVector(theta_force_, "theta_force");
        x_full_.AddVector(theta_stiffness_, "theta_stiffness");
        x_full_.AddVector(theta_mass_, "theta_mass");
        x_full_.AddVector(theta_damp_, "theta_damp");

        // Initialize vector collection.
        if (stable_.find("displacement") != stable_.end())
            x_.AddVector(x_full_.GetVector("displacement"), "displacement");
        if (stable_.find("velocity") != stable_.end())
            x_.AddVector(x_full_.GetVector("velocity"), "velocity");
        if (stable_.find("theta_force") != stable_.end())
            x_.AddVector(theta_force_, "theta_force");
        if (stable_.find("theta_stiffness") != stable_.end())
            x_.AddVector(theta_stiffness_, "theta_stiffness");
        if (stable_.find("theta_mass") != stable_.end())
            x_.AddVector(theta_mass_, "theta_mass");
        if (stable_.find("theta_damp") != stable_.end())
            x_.AddVector(theta_damp_, "theta_damp");

        Nstate_ = x_.GetM();

#ifdef VERDANDI_STATE_ERROR_SPARSE
        build_diagonal_sparse_matrix(GetNstate(),
                                     state_error_variance_value_,
                                     state_error_variance_);
        build_diagonal_sparse_matrix(GetNstate(),
                                     T(1) / state_error_variance_value_,
                                     state_error_variance_inverse_);
#else
        state_error_variance_.Reallocate(GetNstate(), GetNstate());
        state_error_variance_.SetIdentity();
        Mlt(T(state_error_variance_value_), state_error_variance_);
        state_error_variance_inverse_.Reallocate(GetNstate(), GetNstate());
        state_error_variance_inverse_.SetIdentity();
        Mlt(T(state_error_variance_value_), state_error_variance_inverse_);
#endif
        is_adjoint_initialized_ = false;
    }


    //! Initializes the adjoint model.
    template <class T>
    void ClampedBar<T>::InitializeAdjoint()
    {

        /*** Allocate adjoint variables ***/

        state working_vector(Nstate_);
        working_vector.Fill(T(0));
        SetShape(working_vector, q_);
        working_vector.Nullify();

        is_adjoint_initialized_ = true;
    }


    //! Initializes the current time step for the model.
    template <class T>
    void ClampedBar<T>::InitializeStep()
    {
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    /*
      \param[in] update_force Boolean to indicate if the force has to be
      updated or not.
    */
    template <class T>
    void ClampedBar<T>::Forward(bool update_force)
    {
        /*** Update time ***/

        time_ += Delta_t_;

        /*** Right hand side ***/

        if (update_force)
        {
            AssembleMassMatrix(theta_force_, theta_force_index_);
            state ones(Ndof_);
            ones.Fill(T(1));
            MltAdd(T(sin(Pi_ * (time_ + 0.5 * Delta_t_) / final_time_)),
                   mass_matrix_, ones, T(0), force_);
            force_(0) = T(0);
        }

        /*** Assembling process ***/

        AssembleMassMatrix(theta_mass_, theta_mass_index_);

        AssembleNewMarkMatrix0();
        AssembleNewMarkMatrix1();


        state disp_1(Ndof_);
        state velo_1(Ndof_);
        velo_1.Fill(T(0));

        MltAdd(T(2) / Delta_t_, mass_matrix_, velo_0_, T(1), force_);
        MltAdd(T(1), Newmark_matrix_0_, disp_0_, T(1), force_);

        // Dirichlet conditions .
        force_(0) = T(0);

#ifdef VERDANDI_WITH_DIRECT_SOLVER
        Vector<T> tmp(Ndof_);
        tmp = force_;
        Matrix<T, Symmetric, RowSymSparse> temporary_matrix;
        Copy(Newmark_matrix_1_, temporary_matrix);
        SparseSolve(temporary_matrix, tmp);
        disp_1 = tmp;
#else
        // Initialization of the Gmres parameters.
        int nb_max_iter = 1000;
        double tolerance = 1e-6;
        Iteration<double> iter(nb_max_iter, tolerance);
        Preconditioner_Base precond;
        // No preconditioning.
        iter.SetRestart(5);
        iter.HideMessages();
        Gmres(Newmark_matrix_1_, disp_1, force_, precond, iter);
#endif

        velo_1.Fill(T(0));
        Add(T(-1), velo_0_, velo_1);
        Add(T(2) / Delta_t_, disp_1, velo_1);
        Add(T(-2) / Delta_t_, disp_0_, velo_1);

        // Update.
        disp_0_ = disp_1;
        velo_0_ = velo_1;
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool ClampedBar<T>::HasFinished() const
    {
        return time_ >= final_time_;
    }


    //! Finalizes the current time step for the model.
    template <class T>
    void ClampedBar<T>::FinalizeStep()
    {
    }


    //! Finalizes the model.
    template <class T>
    void ClampedBar<T>::Finalize()
    {
    }


    //! Saves the simulated data.
    /*! It saves the displacement 'disp_0_' and  the velocity 'velo_0_'.
     */
    template <class T>
    void ClampedBar<T>::Save()
    {
        output_saver_.Save(disp_0_, time_, "disp_0");
        output_saver_.Save(velo_0_, time_, "velo_0");
    }


    //! Performs one step backward in adjoint model.
    /*!
      \param[in] observation_term \f$ H^T R^{-1}(y - Hx) \f$.
    */
    template <class T>
    void ClampedBar<T>::BackwardAdjoint(state& observation_term)
    {
        if (!is_adjoint_initialized_)
            InitializeAdjoint();

        state_collection source_term;
        SetShape(observation_term, source_term);

        state force_active;
        force_active.SetData(Ndof_ - 1, force_.GetData() + 1);

        /*** Computes p_f_n ***/

        if (stable_.find("theta_force") != stable_.end())
        {
            state e_i(theta_force_.GetM());
            for (int i = 0; i < theta_force_.GetM(); i++)
            {
                e_i.Fill(T(0));
                e_i(i) = T(1);
                AssembleMassMatrix(e_i, theta_force_index_);
                state ones(Ndof_);
                ones.Fill(T(1));
                MltAdd(T(sin(Pi_ * (time_ + 0.5 * Delta_t_) / final_time_)),
                       mass_matrix_, ones, T(0), force_);
                force_(0) = T(0);
                q_.GetVector("theta_force")(i)
                    += DotProd(q_.GetVector("velocity"), force_active);
            }
        }

        /*** Computes p_s_n and p_m_n ***/

        if (stable_.find("theta_stiffness") != stable_.end() ||
            stable_.find("theta_mass") != stable_.end())
        {

            // Computes x_{n+1} + x_n  and v_{n+1} - v_n.
            state disp(Ndof_), velo(Ndof_);
            Copy(disp_0_, disp);
            Copy(velo_0_, velo);

            Forward(true);

            time_ -= Delta_t_;

            Add(T(1), disp_0_, disp);
            Mlt(T(-1), velo);
            Add(T(1), velo_0_, velo);

            if (stable_.find("theta_stiffness") != stable_.end())
            {
                state e_i(Ntheta_stiffness_);
                for (int i = 0; i < Ntheta_stiffness_; i++)
                {
                    e_i.Fill(T(0));
                    e_i(i) = T(1);
                    AssembleStiffnessMatrix(e_i, theta_stiffness_index_);
                    MltAdd(T(0.5), stiffness_matrix_, disp, T(0), force_);
                    q_.GetVector("theta_stiffness")(i)
                        -= DotProd(q_.GetVector("velocity"), force_active);
                }
            }
            if (stable_.find("theta_mass") != stable_.end())
            {
                state e_i(Ntheta_mass_);
                for (int i = 0; i < Ntheta_mass_; i++)
                {
                    e_i.Fill(T(0));
                    e_i(i) = T(1);
                    AssembleMassMatrix(e_i, theta_mass_index_);
                    MltAdd(T(1) / Delta_t_, mass_matrix_, velo, T(0), force_);
                    q_.GetVector("theta_mass")(i)
                        -= DotProd(q_.GetVector("velocity"), force_active);
                }
            }
        }

        if (stable_.find("theta_damp") != stable_.end())
            throw ErrorUndefined("BackwardAdjoint(state& observation_term)",
                                 "theta_damp");

        state rhs(Ndof_), rhs_active;
        rhs.Fill(T(0));
        rhs_active.SetData(Ndof_ - 1, rhs.GetData() + 1);


        if (stable_.find("theta_force") != stable_.end())
            Add(T(1), source_term.GetVector("theta_force"),
                q_.GetVector("theta_force"));
        if (stable_.find("theta_stiffness") != stable_.end())
            Add(T(1), source_term.GetVector("theta_stiffness"),
                q_.GetVector("theta_stiffness"));
        if (stable_.find("theta_mass") != stable_.end())
            Add(T(1), source_term.GetVector("theta_mass"),
                q_.GetVector("theta_mass"));
        if (stable_.find("theta_damp") != stable_.end())
            Add(T(1), source_term.GetVector("theta_damp"),
                q_.GetVector("theta_damp"));

        /*** Assembling process ***/

        AssembleMassMatrix(theta_mass_, theta_mass_index_);
        AssembleDampMatrix(theta_damp_, theta_damp_index_);
        AssembleStiffnessMatrix(theta_stiffness_, theta_stiffness_index_);
        AssembleNewMarkMatrix0();
        AssembleNewMarkMatrix1();

        /*** Computes K_inv * source_term ***/

        state K_inv_Sd(Ndof_), K_inv_Sd_active;
        K_inv_Sd(0) = 0;
        K_inv_Sd_active.SetData(Ndof_ - 1, K_inv_Sd.GetData() + 1);
        Copy(source_term.GetVector("displacement"), K_inv_Sd_active);
#ifdef VERDANDI_WITH_DIRECT_SOLVER
        Matrix<T, Symmetric, RowSymSparse> temporary_matrix;
        Copy(stiffness_matrix_, temporary_matrix);
        SparseSolve(temporary_matrix, K_inv_Sd);
#else
        state b(K_inv_Sd);
        int nb_max_iter = 1000;
        double tolerance = 1e-6;
        Iteration<double> iter(nb_max_iter, tolerance);
        Preconditioner_Base precond;
        iter.SetRestart(5);
        iter.HideMessages();
        K_inv_Sd.Fill(T(0));
        Gmres(stiffness_matrix_, K_inv_Sd, b, precond, iter);
#endif

        state_collection q_disp_1, q_velo_1;
        state q_disp_1_0(1), q_velo_1_0(1);
        q_disp_1_0(0) = q_velo_1_0(0) = T(0);

        q_disp_1.AddVector(q_disp_1_0, "inactive");
        q_disp_1.AddVector(q_.GetVector("displacement"), "active");
        q_velo_1.AddVector(q_velo_1_0, "inactive");
        q_velo_1.AddVector(q_.GetVector("velocity"), "active");

        /*** Computes right-hand side ***/

        rhs.Fill(T(0));
        MltAdd(T(-2) / Delta_t_, mass_matrix_, q_velo_1, T(0), rhs);
        MltAdd(T(2) / Delta_t_, mass_matrix_, K_inv_Sd, T(1), rhs);
        MltAdd(T(1), damp_matrix_, K_inv_Sd, T(1), rhs);
        Add(T(-1), source_term.GetVector("velocity"), rhs_active);
        MltAdd(T(1), Newmark_matrix_0_, q_disp_1, T(1), rhs);
        rhs(0) = T(0);

        /*** Computes q_disp_n ***/

        state q_disp(Ndof_), q_disp_active;
        q_disp(0) = T(0);
        q_disp_active.SetData(Ndof_ - 1, q_disp.GetData() + 1);
#ifdef VERDANDI_WITH_DIRECT_SOLVER
        Copy(rhs_active, q_disp_active);
        Copy(Newmark_matrix_1_, temporary_matrix);
        GetAndSolveLU(temporary_matrix, q_disp);
#else
        Iteration<double> iter2(nb_max_iter, tolerance);
        Preconditioner_Base precond2;
        iter2.SetRestart(5);
        iter2.HideMessages();
        q_disp.Fill(T(0));
        Gmres(Newmark_matrix_1_, q_disp, rhs, precond2, iter2);
#endif

        /*** Computes q_velo_n ***/

        Add(T(-1), q_disp_active, q_disp_1.GetVector("active"));
        Mlt(T(-1), q_velo_1.GetVector("active"));
        Add(T(2), K_inv_Sd_active, q_velo_1.GetVector("active"));
        Add(T(2) / Delta_t_, q_disp_1.GetVector("active"),
            q_velo_1.GetVector("active"));

        Copy(q_disp_active, q_disp_1.GetVector("active"));

        source_term.Nullify();
        q_disp_1.Nullify();
        q_velo_1.Nullify();
        force_active.Nullify();
        rhs_active.Nullify();
        K_inv_Sd_active.Nullify();
        q_disp_active.Nullify();
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the model to a given state vector.
    /*!
      \param[in,out] x on entry, the state vector to which the model is
      applied; on exit, the state vector after the model is applied.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.
    */
    template <class T>
    void ClampedBar<T>::ApplyOperator(state& x, bool forward,
                                      bool preserve_state,
                                      bool update_force)
    {
        double saved_time = 0;
        state saved_state;
        if (!forward)
            saved_time = GetTime();

        if (preserve_state)
            saved_state.SetData(duplicated_state_);

        duplicated_state_.Nullify();
        duplicated_state_.SetData(x);
        StateUpdated();

        Forward(update_force);

        GetStateCopy(duplicated_state_);

        duplicated_state_.Nullify();

        if (!forward)
            SetTime(saved_time);

        if (preserve_state)
        {
            duplicated_state_.SetData(saved_state);
            StateUpdated();
            saved_state.Nullify();
        }
    }


    //! Applies the tangent linear model to a given vector.
    /*!
      \param[in,out] increment the increment.
    */
    template <class T>
    void ClampedBar<T>
    ::ApplyTangentLinearOperator(state& increment_state)
    {
        double saved_time = 0;
        state saved_state;
        saved_time = GetTime();
        saved_state.Copy(GetFullState());

        /*** Computes x_{n+1} + x_n and v_{n+1} - v_n ***/

        state disp(Ndof_), velo(Ndof_);
        Copy(disp_0_, disp);
        Copy(velo_0_, velo);

        Forward(true);

        Add(T(1), disp_0_, disp);
        Mlt(T(-1), velo);
        Add(T(1), velo_0_, velo);

        /*** Right hand side ***/

        SetTime(saved_time);
        state_collection increment;
        SetShape(increment_state, increment);

        force_.Fill(0);
        if (stable_.find("theta_force") != stable_.end())
        {
            AssembleMassMatrix(increment.GetVector("theta_force"),
                               theta_force_index_);
            state ones(Ndof_);
            ones.Fill(T(1));
            MltAdd(T(sin(Pi_ * (time_ + 0.5 * Delta_t_) / final_time_)),
                   mass_matrix_, ones, T(0), force_);
        }
        // Assemble K.
        if (stable_.find("theta_stiffness") != stable_.end())
        {
            if (stiffness_matrix_.GetM() == 0)
                Copy(mass_matrix_, stiffness_matrix_);
            Fill(T(0), stiffness_matrix_);

            for (int i = 0; i < Nx_; i++)
                AddMatrixPosition(pow(T(2), theta_stiffness_(
                                          theta_stiffness_index_(i)))
                                  * increment.GetVector("theta_stiffness")
                                  (theta_stiffness_index_(i)), i, i,
                                  stiffness_FEM_matrix_, stiffness_matrix_);

            // Boundary condition by pseudo-elimination.
            stiffness_matrix_.Val(0, 0) = 1;
            stiffness_matrix_.Val(0, 1) = 0;
            stiffness_matrix_.Val(1, 0) = 0;

            MltAdd(T(-0.5), stiffness_matrix_, disp, T(1), force_);
        }
        // Assemble M.
        if (stable_.find("theta_mass") != stable_.end())
        {
            AssembleMassMatrix(increment.GetVector("theta_mass"),
                               theta_mass_index_);
            MltAdd(T(-1) / Delta_t_, mass_matrix_, velo, T(1), force_);
        }
        force_(0) = T(0);

        /*** Computes delta_y ***/

        Copy(increment.GetVector("displacement"),
             x_.GetVector("displacement"));
        Copy(increment.GetVector("velocity"), x_.GetVector("velocity"));

        Forward(false);

        Copy(x_.GetVector("displacement"),
             increment.GetVector("displacement"));
        Copy(x_.GetVector("velocity"), increment.GetVector("velocity"));

        SetTime(saved_time);

        GetFullState().Copy(saved_state);
        FullStateUpdated();
    }


    //! Gets the matrix of the tangent linear model.
    /*!
      \param[out] A the matrix of the tangent linear model.
    */
    template <class T>
    void ClampedBar<T>
    ::GetTangentLinearOperator(tangent_linear_operator& A) const
    {
        throw ErrorUndefined("ClampedBar::GetTangentLinearOperator"
                             "(tangent_linear_operator& A) const");
    }


    ////////////////////
    // ACCESS METHODS //
    ////////////////////


    //! Returns the current time.
    /*!
      \return The current time.
    */
    template <class T>
    double ClampedBar<T>::GetTime() const
    {
        return time_;
    }


    //! Sets the time of the model to a given time.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    void ClampedBar<T>::SetTime(double time)
    {
        time_= time;
    }


    //! Returns the state vector size.
    /*!
      \return The state vector size.
    */
    template <class T>
    int ClampedBar<T>::GetNstate() const
    {
        return Nstate_;
    }


    //! Provides the reduced state vector.
    /*!
      The state vector is duplicated.
      \param[out] state the reduced state vector.
    */
    template <class T>
    void ClampedBar<T>
    ::GetStateCopy(state& x)
    {
        x.Reallocate(x_.GetM());
        for (int i = 0; i < x_.GetM(); i++)
            x(i) = x_(i);
    }


    //! Sets the reduced state vector.
    /*!
      Before setting the reduced state vector, special requirements can be
      enforced; e.g. positivity requirement or inferior and superior limits.
      \param[in] state the reduced state vector.
    */
    template <class T>
    void ClampedBar<T>
    ::SetStateCopy(state& x)
    {
        if (x_.GetM() != x.GetM())
            throw ErrorProcessing("ClampedBar::SetState()",
                                  "Operation not permitted:\n x_ is a vector"
                                  " of length " + to_str(x_.GetM()) +
                                  ";\n x is a vector of length "
                                  + to_str(x.GetM()) + ".");
        for (int i = 0; i < x_.GetM(); i++)
            x_(i) = x(i);
    }


    //! Provides the reduced state vector.
    /*!
      \return state the reduced state vector.
    */
    template <class T>
    typename ClampedBar<T>::state&  ClampedBar<T>
    ::GetState()
    {
        duplicated_state_.Reallocate(x_.GetM());
        for (int i = 0; i < x_.GetM(); i++)
            duplicated_state_(i) = x_(i);
        return duplicated_state_;
    }


    //! Performs some calculations when the update of the model state is done.
    template <class T>
    void ClampedBar<T>
    ::StateUpdated()
    {
        for (int i = 0; i < x_.GetM(); i++)
            x_(i) = duplicated_state_(i);
    }


    //! Provides the state lower bound.
    /*!
      \return The state lower bound (componentwise).
    */
    template <class T>
    typename ClampedBar<T>::state& ClampedBar<T>
    ::GetStateLowerBound()
    {
        return lower_bound_;
    }


    //! Provides the state upper bound.
    /*!
      \return The state upper bound (componentwise).
    */
    template <class T>
    typename ClampedBar<T>::state& ClampedBar<T>
    ::GetStateUpperBound()
    {
        return upper_bound_;
    }


    //! Provides the full state vector.
    /*!
      \return state the full state vector.
    */
    template <class T>
    typename ClampedBar<T>::state&
    ClampedBar<T>::GetFullState()
    {
        return GetState();
    }


     //! Performs some calculations when the update of the model state is done.
    template <class T>
    void ClampedBar<T>
    ::FullStateUpdated()
    {
        StateUpdated();
    }


    //! Returns the adjoint state vector.
    /*!
      \param[in] state_adjoint the adjoint state vector.
    */
    template <class T>
    void ClampedBar<T>::GetAdjointState(state& state_adjoint)
    {
        if (!is_adjoint_initialized_)
            InitializeAdjoint();

        state_collection p;
        SetShape(state_adjoint, p);

        AssembleMassMatrix(theta_mass_, theta_mass_index_);
        AssembleDampMatrix(theta_damp_, theta_damp_index_);
        AssembleStiffnessMatrix(theta_stiffness_, theta_stiffness_index_);

        state_collection p_disp, p_velo;
        state p_disp_0(1), p_velo_0(1);
        p_disp_0(0) = p_velo_0(0) = T(0);
        p.GetVector("displacement").Fill(T(0));
        p.GetVector("velocity").Fill(T(0));

        p_disp.AddVector(p_disp_0, "inactive");
        p_disp.AddVector(p.GetVector("displacement"), "active");
        p_velo.AddVector(p_velo_0, "inactive");
        p_velo.AddVector(p.GetVector("velocity"), "active");

        /*** Computes p_disp_n ***/

        state q_velo(Ndof_), q_velo_active;
        q_velo(0) = T(0);
        q_velo_active.SetData(Ndof_ - 1, q_velo.GetData() + 1);
        Copy(q_.GetVector("velocity"), q_velo_active);
        Mlt(T(0.5), q_velo_active);
        Add(T(1) / Delta_t_, q_.GetVector("displacement"),
            q_velo_active);
        MltAdd(T(1), stiffness_matrix_, q_velo, T(0), p_disp);

        /*** Computes p_velo_n ***/

        state q_disp(Ndof_), q_disp_active;
        q_disp(0) = T(0);
        q_disp_active.SetData(Ndof_ - 1, q_disp.GetData() + 1);
        Copy(q_.GetVector("displacement"), q_disp_active);
        MltAdd(T(-0.5), stiffness_matrix_, q_disp, T(0), p_velo);

        q_velo(0) = T(0);
        Copy(q_.GetVector("velocity"), q_velo_active);
        MltAdd(T(0.5), damp_matrix_, q_velo, T(1), p_velo);
        MltAdd(T(1) / Delta_t_, mass_matrix_, q_velo, T(1), p_velo);

        if (stable_.find("theta_force") != stable_.end())
            Copy(q_.GetVector("theta_force"), p.GetVector("theta_force"));
        if (stable_.find("theta_stiffness") != stable_.end())
            Copy(q_.GetVector("theta_stiffness"),
                 p.GetVector("theta_stiffness"));
        if (stable_.find("theta_mass") != stable_.end())
            Copy(q_.GetVector("theta_mass"), p.GetVector("theta_mass"));
        if (stable_.find("theta_damp") != stable_.end())
            Copy(q_.GetVector("theta_damp"), p.GetVector("theta_damp"));

        p.Nullify();
        q_disp_active.Nullify();
        q_velo_active.Nullify();
    }


    //! Sets the adjoint state vector.
    /*!
      \param[out] state_adjoint the adjoint state vector.
    */
    template <class T>
    void ClampedBar<T>::SetAdjointState(const state& state_adjoint)
    {
        if (!is_adjoint_initialized_)
            InitializeAdjoint();

        if (q_.GetM() != state_adjoint.GetM())
            throw ErrorProcessing("ClampedBar::SetStateAdjoint()",
                                  "Operation not permitted:\n p_ is a vector"
                                  " of length " + to_str(q_.GetM()) +
                                  ";\n  state_adjoint is a vector of length "
                                  + to_str(state_adjoint.GetM()) + ".");

        for (int i = 0; i < q_.GetM(); i++)
            q_(i) = state_adjoint(i);
    }


    //! Computes a row of the background error covariance matrix B.
    /*!
      \param[in] row row index.
      \param[out] error_covariance_row the value of row number \a row.
    */
    template <class T>
    void ClampedBar<T>
    ::GetStateErrorVarianceRow(int row, state_error_variance_row&
                               state_error_variance_row)
    {
        GetRow(state_error_variance_, row, state_error_variance_row);
    }


    //! Returns the background error covariance matrix (\f$B\f$).
    /*! Returns the background error covariance matrix (\f$B\f$).
      \return The matrix of the background error covariance.
    */
    template <class T>
    typename ClampedBar<T>::state_error_variance&
    ClampedBar<T>::GetStateErrorVariance()
    {
        return state_error_variance_;
    }


    //! Returns the background error covariance matrix (\f$B\f$).
    /*! Returns the background error covariance matrix (\f$B\f$).
      \return The matrix of the background error covariance.
    */
    template <class T>
    const typename ClampedBar<T>::state_error_variance&
    ClampedBar<T>::GetStateErrorVariance() const
    {
        return state_error_variance_;
    }


    /*! Returns a decomposition of the state error covariance matrix (\f$B\f$)
      as a product \f$LUL^T\f$.
    */
    /*!
      \param[out] L the matrix \f$L\f$.
      \param[out] U the matrix \f$U\f$.
    */
    template <class T>
    void ClampedBar<T>::GetStateErrorVarianceSqrt(
        state_error_variance& L,
        state_error_variance& U)
    {
        int Nreduced = 0;
        for (unsigned int i = 0; i < reduced_.size(); i++)
            Nreduced += x_.GetVector(reduced_[i]).GetSize();
#ifndef VERDANDI_STATE_ERROR_SPARSE
        // Initializes L.
        L.Reallocate(Nstate_, Nreduced);
        L.Fill(T(0));
        for (unsigned int i = 0, l = 0; i < reduced_.size(); i++)
            for(int k = x_.GetIndex(reduced_[i]);
                k < x_.GetIndex(reduced_[i]) +
                    x_.GetVector(reduced_[i]).GetSize(); k++)
                L(k, l++) = 1;

        // Initializes U.
        U.Reallocate(Nreduced,  Nreduced);
        U.Fill(T(0));
        for (int i = 0; i < Nreduced; i++)
            U(i, i) = T(T(1) / state_error_variance_value_);
#else
        // Initializes L.
        Matrix<T, General, ArrayRowSparse> L_array(Nstate_, Nreduced);
        for (unsigned int i = 0, l = 0; i < reduced_.size(); i++)
            for(int k = x_.GetIndex(reduced_[i]);
                k < x_.GetIndex(reduced_[i]) +
                    x_.GetVector(reduced_[i]).GetSize(); k++)
                L_array(k, l++) = 1;
        Copy(L_array, L);

        // Initializes U.
        Matrix<T, General, ArrayRowSparse> U_array(Nreduced,  Nreduced);
        for (int i = 0; i < Nreduced; i++)
            U_array(i, i) = T(T(1) / state_error_variance_value_);
        Copy(U_array, U);
#endif
    }


    //! Returns the inverse of the background error variance (\f$B^{-1}\f$).
    /*!
      \return The inverse of the background error variance (\f$B^{-1}\f$).
    */
    template <class T>
    const typename ClampedBar<T>::state_error_variance&
    ClampedBar<T>::GetStateErrorVarianceInverse() const
    {
        return state_error_variance_inverse_;
    }


    //! Returns the name of the class.
    template <class T>
    string ClampedBar<T>::GetName() const
    {
        return "ClampedBar";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void ClampedBar<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


    /*! Build a vector indexed by the points of the bar to indicate
      which region each point belongs to.
    */
    /*
      \param[in] N the size of the bar.
      \param[in] Nregion the number of region.
      \param[out] index_vector the vector of indexes.
    */
    template <class T>
    void ClampedBar<T>
    ::BuildRegionIndex(int N, int Nregion, Vector<int>& index_vector)
    {
        index_vector.Reallocate(N);
        for(int i = 0; i < N; i++)
            index_vector(i) = i % Nregion;
        Sort(index_vector);
    }


    //! Assembles the NewMark matrices.
    template <class T>
    void ClampedBar<T>
    ::AssembleNewMarkMatrix0()
    {
        Fill(T(0), Newmark_matrix_0_);
        for (int i = 0; i < Nx_; i++)
        {
            // Newmark's matrix at time n.
            AddMatrixPosition(theta_mass_(theta_mass_index_(i)) * 2. /
                              (Delta_t_ * Delta_t_), i, i, mass_FEM_matrix_,
                              Newmark_matrix_0_);

            AddMatrixPosition(theta_damp_(theta_damp_index_(i)) /
                              Delta_t_, i, i, damp_FEM_matrix_,
                              Newmark_matrix_0_);

            AddMatrixPosition(-.5 * pow(T(2.), theta_stiffness_(
                                            theta_stiffness_index_(i))), i, i,
                              stiffness_FEM_matrix_, Newmark_matrix_0_);
        }
    }


    //! Assembles the NewMark matrices.
    template <class T>
    void ClampedBar<T>
    ::AssembleNewMarkMatrix1()
    {
        Fill(T(0), Newmark_matrix_1_);
        for (int i = 0; i < Nx_; i++)
        {
            // Newmark's matrix at time n+1.
            AddMatrixPosition(theta_mass_(theta_mass_index_(i)) * 2. /
                              (Delta_t_ * Delta_t_), i, i, mass_FEM_matrix_,
                              Newmark_matrix_1_);

            AddMatrixPosition(theta_damp_(theta_damp_index_(i)) /
                              Delta_t_, i, i, damp_FEM_matrix_,
                              Newmark_matrix_1_);

            AddMatrixPosition(.5 * pow(T(2.), theta_stiffness_(
                                           theta_stiffness_index_(i))), i, i,
                              stiffness_FEM_matrix_, Newmark_matrix_1_);
        }

        // Boundary condition by pseudo-elimination.
        Newmark_matrix_1_.Val(0, 0) = 1;
        Newmark_matrix_1_.Val(0, 1) = 0;
        Newmark_matrix_1_.Val(1, 0) = 0;

    }

    //! Assembles the Mass matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void ClampedBar<T>
    ::AssembleMassMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        Fill(T(0), mass_matrix_);
        for (int i = 0; i < Nx_; i++)
            AddMatrixPosition(theta(theta_index(i)) , i, i,
                              mass_FEM_matrix_, mass_matrix_);

        // Boundary condition by pseudo-elimination.
        mass_matrix_.Val(0, 0) = 1;
        mass_matrix_.Val(0, 1) = 0;
        mass_matrix_.Val(1, 0) = 0;
    }

    //! Assembles damp matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void ClampedBar<T>
    ::AssembleDampMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        Fill(T(0), damp_matrix_);
        for (int i = 0; i < Nx_; i++)
            AddMatrixPosition(theta(theta_index(i)) , i, i,
                              damp_FEM_matrix_, damp_matrix_);

        // Boundary condition by pseudo-elimination.
        damp_matrix_.Val(0, 0) = 1;
        damp_matrix_.Val(0, 1) = 0;
        damp_matrix_.Val(1, 0) = 0;
    }


    //! Assembles stiffness matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void ClampedBar<T>
    ::AssembleStiffnessMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        if (stiffness_matrix_.GetM() == 0)
            Copy(mass_matrix_, stiffness_matrix_);

        Fill(T(0), stiffness_matrix_);
        for (int i = 0; i < Nx_; i++)
            AddMatrixPosition(pow(T(2), theta(theta_index(i))), i, i,
                              stiffness_FEM_matrix_, stiffness_matrix_);

        // Boundary condition by pseudo-elimination.
        stiffness_matrix_.Val(0, 0) = 1;
        stiffness_matrix_.Val(0, 1) = 0;
        stiffness_matrix_.Val(1, 0) = 0;
    }


    //! Builds skeleton newmark, mass and damp matrices.
    template <class T>
    void ClampedBar<T>
    ::AllocateSparseMatrix()
    {
        Matrix<T, General, ArrayRowSymSparse> tridiagonal(Ndof_, Ndof_);

        // Skeleton.
        for (int i = 0; i < Ndof_ - 1; i++)
        {
            tridiagonal.Get(i, i) = 0;
            tridiagonal.Get(i, i + 1) = 0;
        }
        tridiagonal.Get(Ndof_ - 1, Ndof_ - 1) = 0;

        Copy(tridiagonal, Newmark_matrix_0_);
        Copy(Newmark_matrix_0_, Newmark_matrix_1_);
        Copy(Newmark_matrix_0_, mass_matrix_);
        Copy(Newmark_matrix_0_, damp_matrix_);

    }


    //! Builds a state shape collection over \a x.
    /*
      \param[in] x a given state.
      \param[out] x_collection collection encapsulating \a x.
    */
    template <class T>
    void ClampedBar<T>
    ::SetShape(state& x, state_collection& x_collection) const
    {
        x_collection.Clear();
        state working_vector;
        int position = 0;
        if (stable_.find("displacement") != stable_.end())
        {
            working_vector.SetData(Ndof_ - 1, x.GetData());
            x_collection.AddVector(working_vector, "displacement");
            working_vector.Nullify();
            position += Ndof_ - 1;
        }
        if (stable_.find("velocity") != stable_.end())
        {
            working_vector.SetData(Ndof_ - 1, x.GetData() + position);
            x_collection.AddVector(working_vector, "velocity");
            working_vector.Nullify();
            position += Ndof_ - 1;
        }
        if (stable_.find("theta_force") != stable_.end())
        {
            working_vector.SetData(Ntheta_force_, x.GetData() + position);
            x_collection.AddVector(working_vector, "theta_force");
            working_vector.Nullify();
            position += Ntheta_force_;
        }
        if (stable_.find("theta_stiffness") != stable_.end())
        {
            working_vector.SetData(Ntheta_stiffness_, x.GetData() + position);
            x_collection.AddVector(working_vector, "theta_stiffness");
            working_vector.Nullify();
            position += Ntheta_stiffness_;
        }
        if (stable_.find("theta_mass") != stable_.end())
        {
            working_vector.SetData(Ntheta_mass_, x.GetData() + position);
            x_collection.AddVector(working_vector, "theta_mass");
            working_vector.Nullify();
            position += Ntheta_mass_;
        }
        if (stable_.find("theta_damp") != stable_.end())
        {
            working_vector.SetData(Ntheta_damp_, x.GetData() + position);
            x_collection.AddVector(working_vector, "theta_damp");
            working_vector.Nullify();
            position += Ntheta_damp_;
        }
    }


}

#define VERDANDI_FILE_MODEL_CLAMPEDBAR_CXX
#endif
