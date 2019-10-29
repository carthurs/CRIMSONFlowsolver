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


#ifndef VERDANDI_FILE_MODEL_QUADRATICMODEL_HXX


namespace Verdandi
{


    /////////////////////
    // QUADRATIC MODEL //
    /////////////////////


    //! This class is a quadratic model.
    /*! The model is defined as \f$\frac{\mathrm{d}x_i}{\mathrm{d}t} = x^T Q_i
      x + L_i x + b_i\f$, where \f$Q_i\f$ is a matrix, \f$L_i\f$ is the
      \f$i\f$-th row of the matrix \f$L\f$ and \f$b\f$ a vector.
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class QuadraticModel: public VerdandiBase
    {
    public:
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;
        typedef Matrix<T> tangent_linear_operator;
        typedef Matrix<T> state_error_variance;
        typedef Vector<T> state_error_variance_row;
        typedef Vector<T> state;
        typedef Matrix<T> matrix_state_observation;
        typedef Matrix<T> error_variance;
        typedef Vector<T> uncertain_parameter;

    protected:

        //! Dimension of the state.
        int Nstate_;

        //! State vector.
        Vector<T> state_;

        //! Should the quadratic term be applied?
        bool with_quadratic_term_;
        //! Should the linear term be applied?
        bool with_linear_term_;
        //! Should the constant term be applied?
        bool with_constant_term_;

        //! Quadratic terms.
        vector<Matrix<T> > S_;

        //! Matrix that defines the linear part of the model.
        Matrix<T> L_;

        //! Vector that defines the constant part of the model.
        Vector<T> b_;

        //! Time step.
        double Delta_t_;

        //! Final time of the simulation.
        double final_time_;

        //! Current time.
        double time_;

        //! Temporary variable that stores S times the state vector.
        Vector<T> S_state_;

        /*** Uncertainty on parameters ***/

        //! Should the quadratic term be perturbed?
        bool is_quadratic_perturbed_;
        //! Should the linear term be perturbed?
        bool is_linear_perturbed_;
        //! Should the constant term be perturbed?
        bool is_constant_perturbed_;

        //! Parameters to be perturbed.
        uncertain_parameter parameter_;

        //! List of parameters to be perturbed.
        vector<string> uncertain_parameter_vector_;

        //! Number of "global" parameters to be perturbed.
        int Nglob_parameter_;

        //! Number of parameters to be perturbed.
        int Nparameter_;

        //! Correlations between the constant term and the other terms.
        Vector<T> constant_correlation_;

        //! Name of the probability distribution for the constant term.
        string constant_pdf_;

        //! Mean of the probability distribution for the constant term.
        Vector<T> constant_mean_;

        //! Covariance matrix for the constant term.
        Matrix<T, Symmetric, RowSymPacked> constant_variance_;

        //! PDF parameters for the constant term.
        Vector<T> constant_parameter_;

        //! Correlations between the parameters for the linear term.
        Vector<T> linear_correlation_;

        //! Name of the probability distribution for the linear term.
        string linear_pdf_;

        //! Mean of the probability distribution for the linear term.
        Vector<T> linear_mean_;

        //! Covariance matrix for the linear term.
        Matrix<T, Symmetric, RowSymPacked> linear_variance_;

        //! PDF parameters for the linear term.
        Vector<T> linear_parameter_;

        //! Correlations between the quadratic term and the other terms.
        Vector<T> quadratic_correlation_;

        //! Name of the probability distribution for the quadratic term
        string quadratic_pdf_;

        //! Mean of the probability distribution for the quadratic term.
        Vector<T> quadratic_mean_;

        //! Covariance matrix for the quadratic term
        Matrix<T, Symmetric, RowSymPacked> quadratic_variance_;

        //! PDF parameters for the quadratic term.
        Vector<T> quadratic_parameter_;


        /*** Errors ***/

        //! Variance of the model error.
        error_variance Q_;

        //! Variance of the model error in square root form.
        error_variance Q_sqrt_;

        //! Variance of the state error.
        error_variance P_;

        //! Variance of the state error in square root form.
        error_variance P_sqrt_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:
        // Constructors and destructor.
        QuadraticModel();
        QuadraticModel(string configuration_file);
        ~QuadraticModel();
        // Initializations.
        void Initialize(string configuration_file);
        void InitializeStep();

        // Processing.
        void Forward();
        void ApplyOperator(state& x,
                           bool forward = false, bool preserve_state = true);
        void ApplyTangentLinearOperator(state& x);
        void GetTangentLinearOperator(tangent_linear_operator&) const;
        bool HasFinished() const;
        void Save();

        void FinalizeStep();
        void Finalize();

        // Access methods.
        T GetDelta_t() const;
        double GetTime() const;
        void SetTime(double time);
        int GetNstate() const;
        int GetNfull_state() const;
        state& GetState();
        void StateUpdated();
        state& GetFullState();
        void FullStateUpdated();
        std::pair<int, int> GetParameterIndex(int i);
        int GetNparameter();
        uncertain_parameter& GetParameter(int i);
        void SetParameter(int i, uncertain_parameter& parameter);
        Vector<T>& GetParameterCorrelation(int i);
        string GetParameterPDF(int i);
        Matrix<T, Symmetric, Seldon::RowSymPacked>&
        GetParameterVariance(int i);
        Vector<T>& GetParameterParameter(int i);
        string GetParameterOption(int i);

        // Errors.
        error_variance& GetErrorVariance();
#ifndef SWIG
        const error_variance& GetErrorVariance() const;
#endif
        error_variance& GetErrorVarianceSqrt();
#ifndef SWIG
        const error_variance& GetErrorVarianceSqrt() const;
#endif
        state_error_variance& GetStateErrorVariance();
#ifndef SWIG
        const state_error_variance& GetStateErrorVariance() const;
#endif
        void GetStateErrorVarianceRow(int row, state_error_variance_row&
                                      state_error_variance_row);
        state_error_variance& GetStateErrorVarianceSqrt();
#ifndef SWIG
        const state_error_variance& GetStateErrorVarianceSqrt() const;
#endif

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_QUADRATICMODEL_HXX
#endif
