// Copyright (C) 2011-2012 INRIA
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


#ifndef VERDANDI_FILE_MODEL_PETSCCLAMPEDBAR_CXX


#include "PetscClampedBar.hxx"

#include "OutputSaver.cxx"

#include "seldon/vector/VectorCollection.cxx"


namespace Verdandi
{

    ///////////////////
    // STATIC FIELDS //
    ///////////////////


    template <class T>
    const double PetscClampedBar<T>::Pi_ = 3.141592653589793238462;


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Constructor.
    template <class T>
    PetscClampedBar<T>::PetscClampedBar(): time_(0.)
    {
    }


    //! Constructor.
    /*! It builds allocates the state vectors.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    PetscClampedBar<T>::PetscClampedBar(string configuration_file):
    time_(0.)
    {
    }


    //! Destructor.
    template <class T>
    PetscClampedBar<T>::~PetscClampedBar()
    {
        parameter_.Nullify();
    }


    ////////////////
    // INITIALIZE //
    ////////////////


    //! Initializes the model.
    /*!
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void PetscClampedBar<T>::Initialize(string configuration_file)
    {
        mpi_communicator_ = PETSC_COMM_WORLD;
        int ierr;
        ierr = MPI_Comm_rank(mpi_communicator_, &rank_);
        CHKERRABORT(mpi_communicator_, ierr);
        ierr = MPI_Comm_size(mpi_communicator_, &Nprocess_);
        CHKERRABORT(mpi_communicator_, ierr);

        /*** Configuration ***/

        VerdandiOps configuration(configuration_file);

        configuration.SetPrefix("petsc_clamped_bar.domain.");
        configuration.Set("bar_length", bar_length_);
        configuration.Set("Nx", Nx_);
        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("final_time", final_time_);
        configuration.SetPrefix("petsc_clamped_bar.error_statistics.");
        configuration.Set("state_error_variance", "v >= 0",
                          state_error_variance_value_);
        configuration.SetPrefix("petsc_clamped_bar.physics.");
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
        configuration.SetPrefix("petsc_clamped_bar.");
        vector<string> stable;
        configuration.Set("state","ops_in(v, {'displacement', 'velocity', "
                          "'theta_force', 'theta_stiffness', 'theta_mass', "
                          "'theta_damp'})",
                          stable);
        if (stable.size() <= 0)
            throw ErrorArgument("void PetscClampedBar<T>::Initialize"
                                "(string configuration_file)", "The name of"
                                " the underlying state vector (state) are "
                                "not defined.");
        for (unsigned int i = 0; i < stable.size(); i++)
            stable_.insert(stable[i]);

        configuration.Set("reduced_state", reduced_);

        /*** Ouput saver ***/

        output_saver_.Initialize(configuration_file,
                                 "petsc_clamped_bar.output_saver.");
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
        mass_FEM_matrix_(1, 0) = mass_FEM_matrix_(0, 1);

        /*** Elementary stifness matrix construction ***/

        stiffness_FEM_matrix_.Reallocate(2, 2);
        T stiff_lin = T(Young_modulus_ / Delta_x_);
        stiffness_FEM_matrix_(0, 0) = stiff_lin;
        stiffness_FEM_matrix_(1, 1) = stiff_lin;
        stiffness_FEM_matrix_(0, 1) = -stiff_lin;
        stiffness_FEM_matrix_(1, 0) = stiffness_FEM_matrix_(0, 1);

        /*** Elementary damping matrix construction ***/

        damp_FEM_matrix_.Reallocate(2, 2);
        damp_FEM_matrix_(0, 0) = alpha_ * mass_FEM_matrix_(0, 0) +
            beta_ * stiffness_FEM_matrix_(0, 0);
        damp_FEM_matrix_(1, 1) = alpha_ * mass_FEM_matrix_(1, 1) +
            beta_ * stiffness_FEM_matrix_(1, 1);
        damp_FEM_matrix_(0, 1) = alpha_ * mass_FEM_matrix_(0, 1) +
            beta_ * stiffness_FEM_matrix_(0, 1);
        damp_FEM_matrix_(1, 0) = damp_FEM_matrix_(0, 1);

        /*** State initialization ***/

        displacement_0_.Reallocate(Ndof_);
        displacement_0_.Fill(T(0));
        velocity_0_.Reallocate(Ndof_);
        velocity_0_.Fill(T(0));
        rhs_.Reallocate(Ndof_);
        rhs_.Fill(T(1));

        // Initialize parameter collection.
        if (stable_.find("theta_force") != stable_.end())
            parameter_.AddVector(theta_force_, "theta_force");
        if (stable_.find("theta_stiffness") != stable_.end())
            parameter_.AddVector(theta_stiffness_, "theta_stiffness");
        if (stable_.find("theta_mass") != stable_.end())
            parameter_.AddVector(theta_mass_, "theta_mass");
        if (stable_.find("theta_damp") != stable_.end())
            parameter_.AddVector(theta_damp_, "theta_damp");

        Nreduced_ = 0;
        for (unsigned int i = 0; i < reduced_.size(); i++)
            Nreduced_ += parameter_.GetVector(reduced_[i]).GetSize();
        Nstate_ = 2 * Ndof_ + Nreduced_;

        Nstate_local_ = 2 * displacement_0_.GetLocalM();
        if (rank_ == Nprocess_ - 1)
            Nstate_local_ += Nreduced_;

        state_.Reallocate(Nstate_, Nstate_local_);

        newmark_0_assembled_ = false;
        newmark_1_assembled_ = false;
    }


    //! Finalizes the model.
    template <class T>
    void PetscClampedBar<T>::Finalize()
    {
    }


    //! Initializes the adjoint model.
    template <class T>
    void PetscClampedBar<T>::InitializeAdjoint()
    {
        throw ErrorUndefined("PetscClampedBar<T>::InitializeAdjoint()");
    }


    //! Initializes the first time step for the model.
    template <class T>
    void PetscClampedBar<T>::InitializeFirstStep()
    {
    }


    //! Initializes the current time step for the model.
    template <class T>
    void PetscClampedBar<T>::InitializeStep()
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
    void PetscClampedBar<T>::Forward(bool update_force)
    {
        //Timer timer;
        /*** Update time ***/

        time_ += Delta_t_;

        /*** Right hand side ***/

        if (update_force)
        {
            AssembleMassMatrix(theta_force_, theta_force_index_);
            state one(Ndof_);
            one.Fill(T(1));
            MltAdd(T(sin(Pi_ * (time_ + 0.5 * Delta_t_) / final_time_)),
                   mass_, one, T(0), rhs_);
            rhs_.SetBuffer(0, T(0));
            rhs_.Flush();
        }

        /*** Assembling process ***/

        AssembleMassMatrix(theta_mass_, theta_mass_index_);

        AssembleNewMarkMatrix0();
        AssembleNewMarkMatrix1();

        state displacement_1(Ndof_);
        state velocity_1(Ndof_);
        velocity_1.Fill(T(0));

        MltAdd(T(2) / Delta_t_, mass_, displacement_0_, T(1), rhs_);
        MltAdd(T(1), newmark_0_, displacement_0_, T(1), rhs_);

        rhs_.SetBuffer(0, T(0));
        rhs_.Flush();

        PetscGmres(newmark_1_, displacement_1, rhs_);

        velocity_1.Fill(T(0));
        Add(T(-1), velocity_0_, velocity_1);
        Add(T(2) / Delta_t_, displacement_1, velocity_1);
        Add(T(-2) / Delta_t_, displacement_0_, velocity_1);

        displacement_0_.Copy(displacement_1);
        velocity_0_.Copy(velocity_1);
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool PetscClampedBar<T>::HasFinished() const
    {
        return time_ >= final_time_;
    }


    //! Saves the simulated data.
    /*! It saves the displacement 'disp_0_' and  the velocity 'velo_0_'.
     */
    template <class T>
    void PetscClampedBar<T>::Save()
    {
        output_saver_.Save(displacement_0_, time_, "displacement");
        output_saver_.Save(velocity_0_, time_, "velocity");
    }


    //! Performs one step backward in adjoint model.
    /*!
      \param[in] observation_term \f$ H^T R^{-1}(y - Hx) \f$.
    */
    template <class T>
    void PetscClampedBar<T>::BackwardAdjoint(state& observation_term)
    {
        throw ErrorUndefined("PetscClampedBar<T>::BackwardAdjoint");
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
    void PetscClampedBar<T>::ApplyOperator(state& x, bool forward,
                                           bool preserve_state,
                                           bool update_force)
    {
        double saved_time = 0;
        state saved_state;
        if (!forward)
            saved_time = GetTime();

        if (preserve_state)
            GetStateCopy(saved_state);

        state_.Copy(x);
        StateUpdated();

        Forward(update_force);

        x.Copy(GetState());

        if (!forward)
            SetTime(saved_time);

        if (preserve_state)
        {
            state_.Copy(saved_state);
            StateUpdated();
        }
    }


    //! Applies the tangent linear model to a given vector.
    /*!
      \param[in,out] increment the increment.
    */
    template <class T>
    void PetscClampedBar<T>
    ::ApplyTangentLinearOperator(state& increment_state)
    {
        throw
            ErrorUndefined("PetscClampedBar<T>::ApplyTangentLinearOperator");
    }


    //! Gets the matrix of the tangent linear model.
    /*!
      \param[out] A the matrix of the tangent linear model.
    */
    template <class T>
    void PetscClampedBar<T>
    ::GetTangentLinearOperator(tangent_linear_operator& A) const
    {
        throw ErrorUndefined("PetscClampedBar<T>::GetTangentLinearOperator"
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
    double PetscClampedBar<T>::GetTime() const
    {
        return time_;
    }


    //! Sets the time of the model to a given time.
    /*!
      \param[in] time a given time.
    */
    template <class T>
    void PetscClampedBar<T>::SetTime(double time)
    {
        time_= time;
    }


    //! Returns the state vector size.
    /*!
      \return The state vector size.
    */
    template <class T>
    int PetscClampedBar<T>::GetNstate() const
    {
        return Nstate_;
    }


    //! Returns the state vector local size.
    /*!
      \return The state vector local size.
    */
    template <class T>
    int PetscClampedBar<T>::GetLocalNstate() const
    {
        return Nstate_local_;
    }


    //! Provides the reduced state vector.
    /*!
      \param[out] state the reduced state vector.
    */
    template <class T>
    void PetscClampedBar<T>
    ::GetStateCopy(state& x)
    {
        int disp_start, disp_end;
        displacement_0_.GetProcessorRange(disp_start, disp_end);
        int state_start, state_end;
        state_.GetProcessorRange(state_start, state_end);

        for (int i = disp_start; i < disp_end; i++)
            state_.SetBuffer(state_start++, displacement_0_(i));

        int velo_start, velo_end;
        velocity_0_.GetProcessorRange(velo_start, velo_end);
        for (int i = velo_start; i < velo_end; i++)
            state_.SetBuffer(state_start++, velocity_0_(i));

        if (rank_ == Nprocess_ - 1)
            for (unsigned int i = 0; i < reduced_.size(); i++)
                for (int j = 0; j < parameter_.GetVector(reduced_[i])
                         .GetSize(); j++)
                    state_.SetBuffer(state_start++,
                                     parameter_.GetVector(reduced_[i])(j));
        state_.Flush();
        x.Copy(state_);
    }


    //! Sets the reduced state vector.
    /*! Before setting the reduced state vector, special requirements can be
      enforced; e.g. positivity requirement or inferior and superior limits.
      \param[in] state the reduced state vector.
    */
    template <class T>
    void PetscClampedBar<T>
    ::SetStateCopy(state& x)
    {
        int disp_start, disp_end;
        displacement_0_.GetProcessorRange(disp_start, disp_end);
        int velo_start, velo_end;
        velocity_0_.GetProcessorRange(velo_start, velo_end);
        int state_start, state_end;
        state_.GetProcessorRange(state_start, state_end);

        for (int i = disp_start; i < disp_end; i++)
            displacement_0_.SetBuffer(i, x(state_start++));

        for (int i = velo_start; i < velo_end; i++)
            velocity_0_.SetBuffer(i, x(state_start++));

        if (rank_ == Nprocess_ - 1)
            for (unsigned int i = 0; i < reduced_.size(); i++)
                for (int j = 0; j < parameter_.GetVector(reduced_[i])
                         .GetSize(); j++)
                    parameter_.GetVector(reduced_[i])(j) = x(state_start++);

        displacement_0_.Flush();
        velocity_0_.Flush();

        for (unsigned int i = 0; i < reduced_.size(); i++)
            MPI_Bcast(parameter_.GetVector(reduced_[i]).GetData(),
                      parameter_.GetVector(reduced_[i]).GetM(),
                      MPI_DOUBLE, Nprocess_ - 1, mpi_communicator_);
    }


    //! Provides the reduced state vector.
    /*!
      \return state the reduced state vector.
    */
    template <class T>
    typename PetscClampedBar<T>::state& PetscClampedBar<T>
    ::GetState()
    {
        int disp_start, disp_end;
        displacement_0_.GetProcessorRange(disp_start, disp_end);
        int state_start, state_end;
        state_.GetProcessorRange(state_start, state_end);

        for (int i = disp_start; i < disp_end; i++)
            state_.SetBuffer(state_start++, displacement_0_(i));

        int velo_start, velo_end;
        velocity_0_.GetProcessorRange(velo_start, velo_end);
        for (int i = velo_start; i < velo_end; i++)
            state_.SetBuffer(state_start++, velocity_0_(i));

        if (rank_ == Nprocess_ - 1)
            for (unsigned int i = 0; i < reduced_.size(); i++)
                for (int j = 0; j < parameter_.GetVector(reduced_[i])
                         .GetSize(); j++)
                    state_.SetBuffer(state_start++,
                                     parameter_.GetVector(reduced_[i])(j));
        state_.Flush();
        return state_;
    }


    //! Performs some calculations when the update of the model state is done.
    template <class T>
    void PetscClampedBar<T>
    ::StateUpdated()
    {
        int disp_start, disp_end;
        displacement_0_.GetProcessorRange(disp_start, disp_end);
        int velo_start, velo_end;
        velocity_0_.GetProcessorRange(velo_start, velo_end);
        int state_start, state_end;
        state_.GetProcessorRange(state_start, state_end);

        for (int i = disp_start; i < disp_end; i++)
            displacement_0_.SetBuffer(i, state_(state_start++));

        for (int i = velo_start; i < velo_end; i++)
            velocity_0_.SetBuffer(i, state_(state_start++));

        if (rank_ == Nprocess_ - 1)
            for (unsigned int i = 0; i < reduced_.size(); i++)
                for (int j = 0; j < parameter_.GetVector(reduced_[i])
                         .GetSize(); j++)
                    parameter_.GetVector(reduced_[i])(j)
                        = state_(state_start++);

        displacement_0_.Flush();
        velocity_0_.Flush();

        for (unsigned int i = 0; i < reduced_.size(); i++)
            MPI_Bcast(parameter_.GetVector(reduced_[i]).GetData(),
                      parameter_.GetVector(reduced_[i]).GetM(),
                      MPI_DOUBLE, Nprocess_ - 1, mpi_communicator_);
    }


    //! Provides the state lower bound.
    /*!
      \return The state lower bound (componentwise).
    */
    template <class T>
    typename PetscClampedBar<T>::state& PetscClampedBar<T>
    ::GetStateLowerBound()
    {
        throw ErrorUndefined("PetscClampedBar<T>::GetStateLowerBound");
    }


    //! Provides the state upper bound.
    /*!
      \return The state upper bound (componentwise).
    */
    template <class T>
    typename PetscClampedBar<T>::state& PetscClampedBar<T>
    ::GetStateUpperBound()
    {
        throw ErrorUndefined("PetscClampedBar<T>::GetStateUpperBound");
    }


    //! Provides the full state vector.
    /*!
      \return state the full state vector.
    */
    template <class T>
    typename PetscClampedBar<T>::state&
    PetscClampedBar<T>::GetFullState()
    {
        return GetState();
    }


    //! Performs some calculations when the update of the model state is done.
    template <class T>
    void PetscClampedBar<T>
    ::FullStateUpdated()
    {
        StateUpdated();
    }


    //! Returns the adjoint state vector.
    /*!
      \param[in] state_adjoint the adjoint state vector.
    */
    template <class T>
    void PetscClampedBar<T>::GetAdjointState(state& state_adjoint)
    {
        throw ErrorUndefined("PetscClampedBar<T>::GetAdjointState");
    }


    //! Sets the adjoint state vector.
    /*!
      \param[out] state_adjoint the adjoint state vector.
    */
    template <class T>
    void PetscClampedBar<T>::SetAdjointState(const state& state_adjoint)
    {
        throw ErrorUndefined("PetscClampedBar<T>::SetAdjointState");
    }


    //! Computes a row of the background error covariance matrix B.
    /*!
      \param[in] row row index.
      \param[out] error_covariance_row the value of row number \a row.
    */
    template <class T>
    void PetscClampedBar<T>
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
    typename PetscClampedBar<T>::state_error_variance&
    PetscClampedBar<T>::GetStateErrorVariance()
    {
        return state_error_variance_;
    }


    //! Returns the background error covariance matrix (\f$B\f$).
    /*! Returns the background error covariance matrix (\f$B\f$).
      \return The matrix of the background error covariance.
    */
    template <class T>
    const typename PetscClampedBar<T>::state_error_variance&
    PetscClampedBar<T>::GetStateErrorVariance() const
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
    template <class L_matrix, class U_matrix>
    void PetscClampedBar<T>
    ::GetStateErrorVarianceSqrt(L_matrix& L, U_matrix& U)
    {
        int Nreduced = 0;
        for (unsigned int i = 0; i < reduced_.size(); i++)
            Nreduced += parameter_.GetVector(reduced_[i]).GetSize();
        L.Reallocate(Nstate_, Nreduced_, Nstate_local_);
        if (rank_ == Nprocess_ - 1)
            for (int i = 0; i < Nreduced_; i++)
                L.SetBuffer(Nstate_ - 1 - i, Nreduced_ - 1 - i, T(1));
        L.Flush();

        U.Reallocate(Nreduced_, Nreduced_);
        U.Fill(T(0));
        for (int i = 0; i < Nreduced_; i++)
            U(i, i) = T(T(1) / state_error_variance_value_);
    }


    //! Returns the inverse of the background error variance (\f$B^{-1}\f$).
    /*!
      \return The inverse of the background error variance (\f$B^{-1}\f$).
    */
    template <class T>
    const typename PetscClampedBar<T>::state_error_variance&
    PetscClampedBar<T>::GetStateErrorVarianceInverse() const
    {
        return state_error_variance_inverse_;
    }


    //! Returns the name of the class.
    template <class T>
    string PetscClampedBar<T>::GetName() const
    {
        return "PetscClampedBar";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void PetscClampedBar<T>::Message(string message)
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
    void PetscClampedBar<T>
    ::BuildRegionIndex(int N, int Nregion, Vector<int>& index_vector)
    {
        index_vector.Reallocate(N);
        for(int i = 0; i < N; i++)
            index_vector(i) = i % Nregion;
        Sort(index_vector);
    }


    //! Assembles the NewMark matrices.
    template <class T>
    void PetscClampedBar<T>
    ::AssembleNewMarkMatrix0()
    {
        if (newmark_0_assembled_)
            return;
        T value;
        int i, j, loc_start, loc_end;

        newmark_0_.Zero();
        newmark_0_.GetProcessorRowRange(loc_start, loc_end);
        if (loc_start > 0)
        {
            i = loc_start;
            j = loc_start;
            value = theta_mass_(theta_mass_index_(i - 1)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 1) +
                theta_damp_(theta_damp_index_(i)) /
                Delta_t_ * damp_FEM_matrix_(1, 1) +
                -.5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(i)))
                * stiffness_FEM_matrix_(1, 1);
            newmark_0_.SetBuffer(i, j, value, ADD_VALUES);

            j = i - 1;
            value = theta_mass_(theta_mass_index_(i - 1)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 0) +
                theta_damp_(theta_damp_index_(i)) /
                Delta_t_ * damp_FEM_matrix_(1, 0) +
                -.5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(i)))
                * stiffness_FEM_matrix_(1, 0);
            newmark_0_.SetBuffer(i, j, value, ADD_VALUES);
        }
        for (int k = loc_start; k < loc_end && k < Ndof_ - 1; k++)
        {
            value = theta_mass_(theta_mass_index_(k)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(0, 0) +
                theta_damp_(theta_damp_index_(k)) /
                Delta_t_ * damp_FEM_matrix_(0, 0) +
                -.5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(k)))
                * stiffness_FEM_matrix_(0, 0);
            i = k;
            j = k;
            newmark_0_.SetBuffer(i, j, value, ADD_VALUES);

            j = i + 1;
            value = theta_mass_(theta_mass_index_(k)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(0, 1) +
                theta_damp_(theta_damp_index_(k)) /
                Delta_t_ * damp_FEM_matrix_(0, 1) +
                -.5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(k)))
                * stiffness_FEM_matrix_(0, 1);
            newmark_0_.SetBuffer(i, j, value, ADD_VALUES);

            if (i < loc_end - 1)
            {
                i = k + 1;
                j = i;
                value = theta_mass_(theta_mass_index_(k)) * 2. /
                    (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 1) +
                    theta_damp_(theta_damp_index_(k)) /
                    Delta_t_ * damp_FEM_matrix_(1, 1) +
                    -.5 * pow(T(2.),
                              theta_stiffness_(theta_stiffness_index_(k))) *
                    stiffness_FEM_matrix_(1, 1);
                newmark_0_.SetBuffer(i, j, value, ADD_VALUES);

                j = k;
                value = theta_mass_(theta_mass_index_(k)) * 2. /
                    (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 0) +
                    theta_damp_(theta_damp_index_(k)) /
                    Delta_t_ * damp_FEM_matrix_(1, 0) +
                    -.5 * pow(T(2.),
                              theta_stiffness_(theta_stiffness_index_(k))) *
                    stiffness_FEM_matrix_(1, 0);
                newmark_0_.SetBuffer(i, j, value, ADD_VALUES);
            }
        }
        newmark_0_.Flush();

        // Boundary condition by pseudo-elimination.
        newmark_0_.SetBuffer(0, 0, T(1), INSERT_VALUES);
        newmark_0_.SetBuffer(1, 0, T(0), INSERT_VALUES);
        newmark_0_.SetBuffer(0, 1, T(0), INSERT_VALUES);
        newmark_0_.Flush();
        newmark_0_assembled_ = true;
    }


    //! Assembles the NewMark matrices.
    template <class T>
    void PetscClampedBar<T>
    ::AssembleNewMarkMatrix1()
    {
        if (newmark_1_assembled_)
            return;
        T value;
        int i, j, loc_start, loc_end;

        newmark_1_.Zero();
        newmark_1_.GetProcessorRowRange(loc_start, loc_end);
        if (loc_start > 0)
        {
            i = loc_start;
            j = loc_start;
            value = theta_mass_(theta_mass_index_(i - 1)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 1) +
                theta_damp_(theta_damp_index_(i)) /
                Delta_t_ * damp_FEM_matrix_(1, 1) +
                .5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(i)))
                * stiffness_FEM_matrix_(1, 1);
            newmark_1_.SetBuffer(i, j, value, ADD_VALUES);

            j = i - 1;
            value = theta_mass_(theta_mass_index_(i - 1)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 0) +
                theta_damp_(theta_damp_index_(i)) /
                Delta_t_ * damp_FEM_matrix_(1, 0) +
                .5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(i)))
                * stiffness_FEM_matrix_(1, 0);
            newmark_1_.SetBuffer(i, j, value, ADD_VALUES);
        }
        for (int k = loc_start; k < loc_end && k < Ndof_ - 1; k++)
        {
            value = theta_mass_(theta_mass_index_(k)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(0, 0) +
                theta_damp_(theta_damp_index_(k)) /
                Delta_t_ * damp_FEM_matrix_(0, 0) +
                .5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(k)))
                * stiffness_FEM_matrix_(0, 0);
            i = k;
            j = k;
            newmark_1_.SetBuffer(i, j, value, ADD_VALUES);

            j = i + 1;
            value = theta_mass_(theta_mass_index_(k)) * 2. /
                (Delta_t_ * Delta_t_) * mass_FEM_matrix_(0, 1) +
                theta_damp_(theta_damp_index_(k)) /
                Delta_t_ * damp_FEM_matrix_(0, 1) +
                .5 * pow(T(2.), theta_stiffness_(theta_stiffness_index_(k)))
                * stiffness_FEM_matrix_(0, 1);
            newmark_1_.SetBuffer(i, j, value, ADD_VALUES);

            if (i < loc_end - 1)
            {
                i = k + 1;
                j = i;
                value = theta_mass_(theta_mass_index_(k)) * 2. /
                    (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 1) +
                    theta_damp_(theta_damp_index_(k)) /
                    Delta_t_ * damp_FEM_matrix_(1, 1) +
                    .5 * pow(T(2.),
                             theta_stiffness_(theta_stiffness_index_(k))) *
                    stiffness_FEM_matrix_(1, 1);
                newmark_1_.SetBuffer(i, j, value, ADD_VALUES);

                j = k;
                value = theta_mass_(theta_mass_index_(k)) * 2. /
                    (Delta_t_ * Delta_t_) * mass_FEM_matrix_(1, 0) +
                    theta_damp_(theta_damp_index_(k)) /
                    Delta_t_ * damp_FEM_matrix_(1, 0) +
                    .5 * pow(T(2.),
                             theta_stiffness_(theta_stiffness_index_(k))) *
                    stiffness_FEM_matrix_(1, 0);
                newmark_1_.SetBuffer(i, j, value, ADD_VALUES);
            }
        }
        newmark_1_.Flush();

        // Boundary condition by pseudo-elimination.
        newmark_1_.SetBuffer(0, 0, T(1), INSERT_VALUES);
        newmark_1_.SetBuffer(1, 0, T(0), INSERT_VALUES);
        newmark_1_.SetBuffer(0, 1, T(0), INSERT_VALUES);
        newmark_1_.Flush();
        newmark_1_assembled_ = true;
    }


    //! Assembles the Mass matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void PetscClampedBar<T>
    ::AssembleMassMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        T value;
        int i, j, loc_start, loc_end;
        mass_.Zero();
        mass_.GetProcessorRowRange(loc_start, loc_end);
        if (loc_start > 0)
        {
            i = loc_start;
            j = loc_start;
            value = theta(theta_index(i - 1)) * mass_FEM_matrix_(1, 1);
            mass_.SetBuffer(i, j, value, ADD_VALUES);
            j = i - 1;
            value = theta(theta_index(i - 1)) * mass_FEM_matrix_(1, 0);
            mass_.SetBuffer(i, j, value, ADD_VALUES);
        }
        for (int k = loc_start; k < loc_end && k < Ndof_ - 1; k++)
        {
            value = theta(theta_index(k)) * mass_FEM_matrix_(0, 0);
            i = k;
            j = k;
            mass_.SetBuffer(i, j, value, ADD_VALUES);
            j = i + 1;
            value = theta(theta_index(k)) * mass_FEM_matrix_(0, 1);
            mass_.SetBuffer(i, j, value, ADD_VALUES);
            if (i < loc_end - 1)
            {
                i = k + 1;
                j = i;
                value = theta(theta_index(k)) * mass_FEM_matrix_(1, 1);
                mass_.SetBuffer(i, j, value, ADD_VALUES);
                j = k;
                value = theta(theta_index(k)) * mass_FEM_matrix_(1, 0);
                mass_.SetBuffer(i, j, value, ADD_VALUES);
            }
        }
        mass_.Flush();

        // Boundary condition by pseudo-elimination.
        mass_.SetBuffer(0, 0, T(1), INSERT_VALUES);
        mass_.SetBuffer(1, 0, T(0), INSERT_VALUES);
        mass_.SetBuffer(0, 1, T(0), INSERT_VALUES);
        mass_.Flush();
    }


    //! Assembles damp matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void PetscClampedBar<T>
    ::AssembleDampMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        throw ErrorUndefined("void PetscClampedBar<T>::AssembleDampMatrix");
    }


    //! Assembles stiffness matrix.
    /*
      \param[in] theta vector of 'theta' value.
      \param[in] theta_index vector that indicates for each element
      the 'theta' value index of the element.
    */
    template <class T>
    void PetscClampedBar<T>
    ::AssembleStiffnessMatrix(Vector<T>& theta, Vector<int>& theta_index)
    {
        throw ErrorUndefined("void PetscClampedBar<T>"
                             "::AssembleStiffnessMatrix");
    }


    //! Builds skeleton newmark, mass and damp matrices.
    template <class T>
    void PetscClampedBar<T>
    ::AllocateSparseMatrix()
    {
        Matrix<T, General, ArrayRowSparse> tridiagonal_array(Ndof_, Ndof_);
        for (int i = 1; i < Ndof_ - 1; i++)
        {
            tridiagonal_array.Get(i, i) = T(0);
            tridiagonal_array.Get(i, i + 1) = T(0);
            tridiagonal_array.Get(i, i - 1) = T(0);
        }
        tridiagonal_array.Get(0, 0) = T(0);
        tridiagonal_array.Get(0, 1) = T(0);
        tridiagonal_array.Get(Ndof_ - 1, Ndof_ - 2) = T(0);
        tridiagonal_array.Get(Ndof_ - 1, Ndof_ - 1) = T(0);

        Matrix<T, General, RowSparse> tridiagonal_rs(Ndof_, Ndof_);
        Copy(tridiagonal_array, tridiagonal_rs);

        mass_.Copy(tridiagonal_rs);
        newmark_0_.Copy(mass_);
        newmark_1_.Copy(mass_);
        damp_.Copy(mass_);
        stiffness_.Copy(mass_);
    }


    template <class T>
    template <class MatrixSparse, class Vector1>
    void PetscClampedBar<T>
    ::PetscGmres(MatrixSparse& A, Vector1& x, const Vector1& b)
    {
        int ierr;
        KSP ksp;
        PC pc;
        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        ierr = KSPSetOperators(ksp, A.GetPetscMatrix(), A.GetPetscMatrix(),
                               DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        ierr = KSPGetPC(ksp,&pc);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        ierr = KSPSetType(ksp, KSPGMRES);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        ierr = KSPSetTolerances(ksp, 1.e-6, PETSC_DEFAULT,
                                PETSC_DEFAULT, PETSC_DEFAULT);
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        ierr = KSPSolve(ksp, b.GetPetscVector(), x.GetPetscVector());
        CHKERRABORT(PETSC_COMM_WORLD, ierr);
        KSPDestroy(&ksp);
    }

}

#define VERDANDI_FILE_MODEL_PETSCCLAMPEDBAR_CXX
#endif
