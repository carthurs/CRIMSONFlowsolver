// Copyright (C) 2010 INRIA
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


#ifndef VERDANDI_FILE_MODEL_MODELTEMPLATE_HXX


namespace Verdandi
{


    ////////////////////
    // MODEL TEMPLATE //
    ////////////////////


    //! This class is a model template.
    class ModelTemplate: public VerdandiBase
    {
    public:
        //! The numerical type (e.g., double).
        typedef double value_type;
#ifdef VERDANDI_SPARSE
        //! Type of the state error variance.
        typedef Matrix<double, General, RowSparse> state_error_variance;
        //! Type of a row of the state error variance.
        typedef Vector<double, VectSparse> state_error_variance_row;
        //! Type of the state/observation crossed matrix.
        typedef Matrix<double, General, RowSparse> matrix_state_observation;
        //! Type of the tangent linear model.
        typedef Matrix<double, General, RowSparse> tangent_linear_operator;
#else
        //! Type of the state error variance.
        typedef Matrix<double> state_error_variance;
        //! Type of a row of the state error variance.
        typedef Vector<double> state_error_variance_row;
        //! Type of the state/observation crossed matrix.
        typedef Matrix<double> matrix_state_observation;
        //! Type of the tangent linear model.
        typedef Matrix<double> tangent_linear_operator;
        //! Type of the model error variance.
        typedef Matrix<double> error_variance;
#endif
        //! Type of the state vector.
        typedef Vector<double> state;
        //! Type of an uncertain parameter.
        typedef Vector<double> uncertain_parameter;



    public:

        // Constructor and destructor.
        ModelTemplate();
        ~ModelTemplate();
        void Initialize(string configuration_file);
        void InitializeStep();

        // Processing.
        void Forward();
        void BackwardAdjoint(state& observation_term);
        bool HasFinished() const;

        void FinalizeStep();
        void Finalize();

        // Operators.
        void ApplyOperator(state& x,
                           bool forward = false, bool preserve_state = true);
        void ApplyTangentLinearOperator(state& x);
        void GetTangentLinearOperator(tangent_linear_operator&) const;

        // Access methods.
        double GetTime() const;
        void SetTime(double time);
        int GetNstate() const;
        int GetNfull_state() const;
        state& GetState();
        void StateUpdated();
        state& GetFullState();
        void FullStateUpdated();
        state& GetStateLowerBound();
        state& GetStateUpperBound();
        void GetAdjointState(state& state_adjoint);
        void SetAdjointState(const state& state_adjoint);

        // Uncertainty.
        int GetNparameter();
        uncertain_parameter& GetParameter(int i);
        void SetParameter(int i, uncertain_parameter parameter);
        Vector<double>& GetParameterCorrelation(int i);
        string GetParameterPDF(int i);
        Matrix<double, Symmetric, RowSymPacked>& GetParameterVariance(int i);
        Vector<double>& GetParameterParameter(int i);
        string GetParameterOption(int i);

        // Errors.
        void GetStateErrorVarianceRow(int row,
                                      state_error_variance_row& P_row);
        state_error_variance& GetStateErrorVariance();
        error_variance& GetErrorVarianceSqrt();
        void GetStateErrorVarianceSqrt(state_error_variance& L,
                                       state_error_variance& U);
        const state_error_variance& GetStateErrorVarianceInverse() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_MODELTEMPLATE_HXX
#endif
