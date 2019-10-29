// Copyright (C) 2009-2010 INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_MODEL_QUADRATICMODEL_CXX


#include "QuadraticModel.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Constructor.
    template <class T>
    QuadraticModel<T>::QuadraticModel(): Delta_t_(1.), time_(0.)
    {
    }


    //! Constructor.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file path to the configuration file.
    */
    template <class T>
    QuadraticModel<T>::QuadraticModel(string configuration_file)
    {
        Initialize(configuration_file);
    }


    //! Destructor.
    template <class T>
    QuadraticModel<T>::~QuadraticModel()
    {
        if (parameter_.GetData() != NULL)
            parameter_.Nullify();
    }


    /////////////////////
    // INITIALIZATIONS //
    /////////////////////


    //! Initializes the model.
    /*! It reads the initial condition and the time settings.
      \param[in] configuration_file configuration file.
    */
    template <class T>
    void QuadraticModel<T>::Initialize(string configuration_file)
    {

        /*** Configuration ***/

        VerdandiOps configuration(configuration_file);

        configuration.SetPrefix("quadratic_model.definition.");

        configuration.Set("initial_state", state_);
        Nstate_ = state_.GetLength();

        configuration.Set("with_quadratic_term", with_quadratic_term_);
        configuration.Set("with_linear_term", with_linear_term_);
        configuration.Set("with_constant_term", with_constant_term_);

        if (with_quadratic_term_)
        {
            S_state_.Reallocate(Nstate_);
            S_.resize(Nstate_);
            for (int i = 0; i < Nstate_; i++)
                S_[i].Reallocate(Nstate_, Nstate_);
            configuration.Set("quadratic_term", S_);
            if (int(S_.size()) != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"quadratic_term\" has "
                                         + to_str(int(S_.size()))
                                         + " elements.");
            for (int i = 0; i < Nstate_; i++)
                if (S_[i].GetM() != Nstate_ || S_[i].GetN() != Nstate_)
                    throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                             "The initial state has "
                                             + to_str(Nstate_) + " elements, "
                                             "but the matrix " + to_str(i)
                                             + " of \"quadratic_term\" has "
                                             + to_str(int(S_[i].GetM()))
                                             + " rows and "
                                             + to_str(int(S_[i].GetN()))
                                             + " columns.");
        }

        if (with_linear_term_)
        {
            L_.Reallocate(Nstate_, Nstate_);
            configuration.Set("linear_term", L_);
            if (L_.GetM() != Nstate_ || L_.GetN() != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"linear_term\" is a "
                                         + to_str(int(L_.GetM())) + " x "
                                         + to_str(int(L_.GetN()))
                                         + " matrix.");
        }

        if (with_constant_term_)
        {
            configuration.Set("constant", b_);
            if (b_.GetLength() != Nstate_)
                throw ErrorConfiguration("QuadraticModel::QuadraticModel",
                                         "The initial state has "
                                         + to_str(Nstate_) + " elements, but "
                                         "the entry \"constant\" has "
                                         + to_str(b_.GetLength())
                                         + " elements.");
        }

        configuration.Set("Delta_t", Delta_t_);
        configuration.Set("initial_time", time_);
        configuration.Set("final_time", final_time_);

        /*** Distribution for the constant term ***/

        configuration.ClearPrefix();
        Nparameter_ = 0;

        if (configuration.Exists("quadratic_model.uncertainty"))
        {
            configuration.SetPrefix("quadratic_model.uncertainty.");

            configuration.Set("uncertain_parameter_list",
                              "ops_in(v, {'quadratic_term', 'linear_term', "
                              "'constant'})",
                              uncertain_parameter_vector_);

            Nglob_parameter_ = int(uncertain_parameter_vector_.size());
            is_quadratic_perturbed_ = false;
            is_linear_perturbed_ = false;
            is_constant_perturbed_ = false;

            if (find (uncertain_parameter_vector_.begin(),
                      uncertain_parameter_vector_.end(),
                      "quadratic_term")
                != uncertain_parameter_vector_.end())
                if (with_quadratic_term_)
                    is_quadratic_perturbed_ = true;
                else
                    Nglob_parameter_--;
            if (find (uncertain_parameter_vector_.begin(),
                      uncertain_parameter_vector_.end(),
                      "linear_term")
                != uncertain_parameter_vector_.end())
                if (with_linear_term_)
                    is_linear_perturbed_ = true;
                else
                    Nglob_parameter_--;
            if (find (uncertain_parameter_vector_.begin(),
                      uncertain_parameter_vector_.end(),
                      "constant")
                != uncertain_parameter_vector_.end())
                if (with_constant_term_)
                    is_constant_perturbed_ = true;
                else
                    Nglob_parameter_--;

            Nparameter_ = Nglob_parameter_;

            if (is_quadratic_perturbed_)
                Nparameter_ += Nstate_ * Nstate_ - 1;
            if (is_linear_perturbed_)
                Nparameter_ += Nstate_ - 1;

            if (is_constant_perturbed_)
            {
                constant_mean_.Reallocate(b_.GetLength());
                configuration.Set("constant.mean", constant_mean_);
                Add(T(1), constant_mean_, b_);

                constant_variance_.Reallocate(Nstate_, Nstate_);
                configuration.Set("constant.variance",
                                  constant_variance_);

                configuration.Set("constant.distribution",
                                  "ops_in(v, {'Normal', 'LogNormal', "
                                  "'NormalHomogeneous',"
                                  "'LogNormalHomogeneous'})",
                                  constant_pdf_);

                configuration.Set("constant.parameter",
                                  constant_parameter_);
            }

            if (is_linear_perturbed_)
            {
                linear_mean_.Reallocate(Nstate_);
                configuration.Set("linear_term.mean", linear_mean_);

                linear_variance_.Reallocate(Nstate_, Nstate_);
                configuration.Set("linear_term.variance",
                                  linear_variance_);

                configuration.Set("linear_term.distribution",
                                  "ops_in(v, {'Normal', 'LogNormal', "
                                  "'NormalHomogeneous', "
                                  "'LogNormalHomogeneous'})",
                                  linear_pdf_);

                configuration.Set("linear_term.parameter",
                                  linear_parameter_);

                for (int i = 0; i < Nstate_; i++)
                {
                    Vector<double> row;
                    GetRow(L_, i, row);
                    Add(T(1), linear_mean_, row);
                    SetRow(row, i, L_);
                }

            }

            if (is_quadratic_perturbed_)
            {
                quadratic_mean_.Reallocate(Nstate_);
                configuration.Set("quadratic_term.mean", quadratic_mean_);

                quadratic_variance_.Reallocate(Nstate_, Nstate_);
                configuration.Set("quadratic_term.variance",
                                  quadratic_variance_);

                configuration.Set("quadratic_term.distribution",
                                  "ops_in(v, {'Normal', 'LogNormal', "
                                  "'NormalHomogeneous',"
                                  "'LogNormalHomogeneous'})",
                                  quadratic_pdf_);

                configuration.Set("quadratic_term.parameter",
                                  quadratic_parameter_);

                Vector<double> row;
                for (int i = 0; i < Nstate_; i++)
                    for (int j = 0; j < Nstate_; j++)
                    {
                        GetRow(S_[i], j, row);
                        Add(T(1), quadratic_mean_, row);
                        SetRow(row, j, S_[i]);
                    }
            }
        }

        /*** Errors ***/

        configuration.SetPrefix("quadratic_model.");

        if (configuration.Exists("error"))
        {
            Q_.Reallocate(Nstate_, Nstate_);
            if (configuration.Get<bool>("error.scaled_identity"))
            {
                Q_.SetIdentity();
                Mlt(configuration.Get<T>("error.diagonal_value", "v >= 0"),
                    Q_);
            }
            else
                configuration.Set("error.value", Q_);
        }

        if (configuration.Exists("error_sqrt"))
        {
            Q_sqrt_.Reallocate(Nstate_, 0);
            configuration.Set("error_sqrt.value", Q_sqrt_);
        }

        if (configuration.Exists("state_error"))
        {
            P_.Reallocate(Nstate_, Nstate_);
            if (configuration.Get<bool>("state_error.scaled_identity"))
            {
                P_.SetIdentity();
                Mlt(configuration.Get<T>("state_error.diagonal_value",
                                         "v >= 0"),
                    P_);
            }
            else
                configuration.Set("state_error.value", P_);
        }

        if (configuration.Exists("state_error_sqrt"))
        {
            P_sqrt_.Reallocate(Nstate_, 0);
            configuration.Set("state_error_sqrt.value", P_sqrt_);
        }

        /*** Output saver ***/

        configuration.SetPrefix("quadratic_model.output_saver.");
        output_saver_.Initialize(configuration);
        if (with_quadratic_term_)
        {
            output_saver_.Empty("S");
            for (int i = 0; i < Nstate_; i++)
                output_saver_.Save(S_[i], "S");
        }
        if (with_linear_term_)
        {
            output_saver_.Empty("L");
            output_saver_.Save(L_, "L");
        }
        if (with_constant_term_)
        {
            output_saver_.Empty("b");
            output_saver_.Save(b_, "b");
        }
        output_saver_.Empty("state");
    }


    //! Initializes the current time step for the model.
    template <class T>
    void QuadraticModel<T>::InitializeStep()
    {
    }


    ////////////////
    // PROCESSING //
    ////////////////


    //! Advances one step forward in time.
    template <class T>
    void QuadraticModel<T>::Forward()
    {
        if (with_quadratic_term_)
        {
            state current_state = state_;
            for (int i = 0; i < Nstate_; i++)
            {
                MltAdd(Delta_t_, S_[i], current_state, T(0), S_state_);
                state_(i) += DotProd(S_state_, current_state);
            }
            if (with_linear_term_)
                MltAdd(Delta_t_, L_, current_state, T(1), state_);
        }
        else if (with_linear_term_)
            MltAdd(Delta_t_, L_, state_, T(1), state_);
        if (with_constant_term_)
            Add(Delta_t_, b_, state_);

        time_ += Delta_t_;
    }


    //! Applies the model to a given state vector.
    /*!
      \param[in,out] c on entry, the state vector to which the model is
      applied; on exit, the state vector after the model is applied.
      \param[in] forward Boolean to indicate if the model has to go on to the
      next step.
      \param[in] preserve_state Boolean to indicate if the model state has to
      be preserved.
    */
    template <class T>
    void QuadraticModel<T>
    ::ApplyOperator(state& x,
                    bool forward, bool preserve_state)
    {
        state saved_state;
        if (preserve_state)
            saved_state.SetData(state_);

        state_.Nullify();
        state_.SetData(x);
        Forward();
        if (!forward)
            time_ -= Delta_t_;
        state_.Nullify();

        if (preserve_state)
        {
            state_.SetData(saved_state);
            saved_state.Nullify();
        }
    }


    //! Applies the tangent linear model to a given vector.
    /*!
      \param[in,out] x on entry, a vector to which the tangent linear model
      should be applied; on exit, the result.
    */
    template <class T>
    void QuadraticModel<T>
    ::ApplyTangentLinearOperator(state& x)
    {
        state input = x;
        if (with_quadratic_term_)
        {
            for (int i = 0; i < Nstate_; i++)
            {
                MltAdd(Delta_t_, S_[i], state_, T(0), S_state_);
                MltAdd(Delta_t_, SeldonTrans, S_[i], state_, T(1), S_state_);
                x(i) += DotProd(S_state_, input);
            }
            if (with_linear_term_)
                MltAdd(Delta_t_, L_, input, T(1), x);
        }
        else if (with_linear_term_)
            MltAdd(Delta_t_, L_, input, T(1), x);
    }


    //! Returns the tangent linear model.
    /*!
      \param[out] M the matrix of the tangent linear model.
    */
    template <class T>
    void QuadraticModel<T>
    ::GetTangentLinearOperator(tangent_linear_operator& M) const
    {
        M.Reallocate(Nstate_, Nstate_);
        if (with_quadratic_term_)
        {
            Vector<T> M_row(Nstate_);
            for (int i = 0; i < Nstate_; i++)
            {
                MltAdd(T(1), S_[i], state_, T(0), M_row);
                MltAdd(T(Delta_t_), SeldonTrans, S_[i], state_,
                       T(Delta_t_), M_row);
                SetRow(M_row, i, M);
                M(i, i) += T(1);
            }
            if (with_linear_term_)
                Add(T(Delta_t_), L_, M);
        }
        else if (with_linear_term_)
        {
            M.Copy(L_);
            Mlt(Delta_t_, M);
            for (int i = 0; i < Nstate_; i++)
                M(i, i) += T(1);
        }
        else
            M.SetIdentity();
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class T>
    bool QuadraticModel<T>::HasFinished() const
    {
        return time_ >= final_time_;
    }


    //! Saves the simulated data.
    /*! It saves the state.
     */
    template <class T>
    void QuadraticModel<T>::Save()
    {
        output_saver_.Save(state_, time_, "state");
    }


    //! Finalizes the current time step for the model.
    template <class T>
    void QuadraticModel<T>::FinalizeStep()
    {
    }


    //! Finalizes the model.
    template <class T>
    void QuadraticModel<T>::Finalize()
    {
    }


    ///////////////////
    // ACCESS METHOD //
    ///////////////////


    //! Returns the time step.
    /*!
      \return The time step.
    */
    template <class T>
    T QuadraticModel<T>::GetDelta_t() const
    {
        return Delta_t_;
    }


    //! Returns the current time.
    /*!
      \return The current time.
    */
    template <class T>
    double QuadraticModel<T>::GetTime() const
    {
        return time_;
    }


    //! Sets the current time.
    /*!
      \param[in] time the current time.
    */
    template <class T>
    void QuadraticModel<T>::SetTime(double time)
    {
        time_ = time;
    }


    //! Returns the dimension of the state.
    /*!
      \return The dimension of the state.
    */
    template <class T>
    int QuadraticModel<T>::GetNstate() const
    {
        return Nstate_;
    }


    //! Returns the dimension of the full state.
    /*!
      \return The dimension of the full state.
    */
    template <class T>
    int QuadraticModel<T>::GetNfull_state() const
    {
        return GetNstate();
    }


    //! Provides the controlled state vector.
    /*!
      \return state the controlled state vector.
    */
    template <class T>
    typename QuadraticModel<T>::state& QuadraticModel<T>
    ::GetState()
    {
        return state_;
    }


    //! Performs some calculations when the update of the model state is done.
    template <class T>
    void QuadraticModel<T>::StateUpdated()
    {
    }


    //! Provides the full state vector.
    /*!
      \return state the full state vector.
    */
    template <class T>
    typename QuadraticModel<T>::state& QuadraticModel<T>
    ::GetFullState()
    {
        return state_;
    }


    /*! \brief Performs some calculations when the update of the full model
      state is done.*/
    template <class T>
    void QuadraticModel<T>::FullStateUpdated()
    {
    }


    /*! \brief Returns the appropriate pair associated with the i-th
      uncertain parameter.*/
    /*!
      \param[in] i uncertain parameter index.
      \return The pair associated with the i-th uncertain parameter.
      First element of the pair is the term to perturb (0 : constant,
      1 : linear_term, 2 : quadratic_term), second element is the row index
      of the term to perturb.
    */
    template<class T>
    std::pair<int, int> QuadraticModel<T>::GetParameterIndex(int i)
    {
        if (is_constant_perturbed_ && i == 0) return make_pair(0, i);

        if (Nglob_parameter_ == 1)
        {
            if (is_linear_perturbed_) return make_pair(1, i);
            else return make_pair(2, i);
        }
        else if (Nglob_parameter_ == 2)
        {
            if (is_constant_perturbed_)
                if (is_linear_perturbed_) return make_pair(1, i - 1);
                else return make_pair(2, i - 1);
            else
                if (i < Nstate_) return make_pair(1, i);
                else return make_pair(2, i - Nstate_);
        }
        else
        {
            if (i < 1 + Nstate_) return make_pair(1, i - 1);
            else return make_pair(2, i - 1 - Nstate_);
        }
    }


    //! Returns the number of parameters to be perturbed.
    /*!
      \return The number of parameters to be perturbed.
    */
    template <class T>
    int QuadraticModel<T>::GetNparameter()
    {
        return Nparameter_;
    }


    //! Gets the i-th uncertain parameter.
    /*!
      \param[in] i index of the parameter.
      \return The vector associated with the i-th parameter.
    */
    template<class T>
    typename QuadraticModel<T>::uncertain_parameter&
    QuadraticModel<T>::GetParameter(int i)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_.GetData() != NULL)
            parameter_.Nullify();
        if (parameter_index.first == 0)
            parameter_.SetData(b_.GetLength(), b_.GetData());
        else if (parameter_index.first == 1)
        {
            Vector<T> row;
            GetRow(L_, parameter_index.second, row);
            parameter_.SetData(row.GetLength(),
                               L_.GetMe()[parameter_index.second]);
        }
        else
        {
            Vector<T> row;
            GetRow(S_[parameter_index.second / Nstate_],
                   parameter_index.second % Nstate_, row);
            parameter_.SetData(row.GetLength(),
                               S_[parameter_index.second / Nstate_].GetMe()
                               [parameter_index.second % Nstate_]);
        }
        return parameter_;
    }


    //! Sets the i-th parameter.
    /*!
      \param[in] i index of the parameter.
      \param[in] parameter the parameter to assign.
    */
    template<class T>
    void QuadraticModel<T>::SetParameter(int i,
                                         uncertain_parameter& parameter)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_index.first == 0)
            b_ = parameter;
        else if (parameter_index.first == 1)
            SetRow(parameter, parameter_index.second, L_);
        else SetRow(parameter, parameter_index.second / Nstate_,
                    S_[parameter_index.second % Nstate_]);
    }


    //! Returns the correlation between the uncertain parameters.
    /*! Since there is only one parameter, an empty vector is
      returned.
      \param[in] i parameter index.
      \return An empty vector.
    */
    template<class T>
    Vector<T>& QuadraticModel<T>::GetParameterCorrelation(int i)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_index.first == 0)
            return constant_correlation_;
        else if (parameter_index.first == 1)
            return linear_correlation_;
        else
            return quadratic_correlation_;
    }


    //! Returns the PDF of the i-th parameter.
    /*!
      \param[in] i uncertain-variable index.
      \return The PDF of the i-th parameter.
    */
    template<class T>
    string QuadraticModel<T>::GetParameterPDF(int i)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_index.first == 0) return constant_pdf_;
        else if (parameter_index.first == 1)
            return linear_pdf_;
        else return quadratic_pdf_;
    }


    /*! \brief Returns the covariance matrix associated with the i-th
      parameter. */
    /*!
      \param[in] i parameter index.
      \return The covariance matrix associated with the i-th parameter.
    */
    template<class T>
    Matrix<T, Symmetric, RowSymPacked>&
    QuadraticModel<T>::GetParameterVariance(int i)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_index.first == 0) return constant_variance_;
        else if (parameter_index.first == 1)
            return linear_variance_;
        else return quadratic_variance_;
    }


    //! Returns parameters associated with the PDF of some model parameter.
    /*! In case of normal or log-normal distribution, the parameters are
      clipping parameters.
      \param[in] i model parameter index.
      \return The parameters associated with the i-th parameter.
    */
    template<class T>
    Vector<T>& QuadraticModel<T>::GetParameterParameter(int i)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_index.first == 0) return constant_parameter_;
        else if (parameter_index.first == 1)
            return linear_parameter_;
        else return quadratic_parameter_;
    }


    //! Returns the perturbation option of the i-th parameter.
    /*!
      \param[in] i parameter index.
      \return The perturbation option of the i-th parameter.
    */
    template<class T>
    string QuadraticModel<T>::GetParameterOption(int i)
    {
        pair<int, int> parameter_index = GetParameterIndex(i);
        if (parameter_index.first == 0) return "init_step";
        else if (parameter_index.first == 1)
            return "init_step";
        else return "init_step";
    }


    ////////////
    // ERRORS //
    ////////////


    //! Returns the model error variance.
    /*!
      \return The model error variance.
    */
    template <class T>
    typename QuadraticModel<T>::error_variance&
    QuadraticModel<T>::GetErrorVariance()
    {
        return Q_;
    }


    //! Returns the model error variance.
    /*!
      \return The model error variance.
    */
    template <class T>
    const typename QuadraticModel<T>::error_variance&
    QuadraticModel<T>::GetErrorVariance() const
    {
        return Q_;
    }


    //! Returns the square root of the model error variance.
    /*!
      \return The square root of the model error variance.
    */
    template <class T>
    typename QuadraticModel<T>::error_variance&
    QuadraticModel<T>::GetErrorVarianceSqrt()
    {
        return Q_sqrt_;
    }


    //! Returns the square root of the model error variance.
    /*!
      \return The square root of the model error variance.
    */
    template <class T>
    const typename QuadraticModel<T>::error_variance&
    QuadraticModel<T>::GetErrorVarianceSqrt() const
    {
        return Q_sqrt_;
    }


    //! Returns the state error variance.
    /*!
      \return The state error variance.
    */
    template <class T>
    typename QuadraticModel<T>::state_error_variance&
    QuadraticModel<T>::GetStateErrorVariance()
    {
        return P_;
    }


    //! Returns the state error variance.
    /*!
      \return The state error variance.
    */
    template <class T>
    const typename QuadraticModel<T>::state_error_variance&
    QuadraticModel<T>::GetStateErrorVariance() const
    {
        return P_;
    }


    //! Returns a row of the state error variance.
    /*!
      \param[in] row row index.
      \param[out] P_row the row with index \a row in the state error variance.
    */
    template <class T>
    void
    QuadraticModel<T>
    ::GetStateErrorVarianceRow
    (int row, state_error_variance_row& P_row)
    {
        GetRow(P_, row, P_row);
    }


    //! Returns the square root of the state error variance.
    /*!
      \return The square root of the state error variance.
    */
    template <class T>
    typename QuadraticModel<T>::state_error_variance&
    QuadraticModel<T>::GetStateErrorVarianceSqrt()
    {
        return P_sqrt_;
    }


    //! Returns the square root of the state error variance.
    /*!
      \return The square root of the state error variance.
    */
    template <class T>
    const typename QuadraticModel<T>::state_error_variance&
    QuadraticModel<T>::GetStateErrorVarianceSqrt() const
    {
        return P_sqrt_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class T>
    string QuadraticModel<T>::GetName() const
    {
        return "QuadraticModel";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class T>
    void QuadraticModel<T>::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
            Save();
    }


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_QUADRATICMODEL_CXX
#endif
