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


#ifndef VERDANDI_FILE_MODEL_CLAMPEDBAR_HXX

#include "seldon/vector/VectorCollection.hxx"

#include "OutputSaver.hxx"

#include <iostream>

#include <set>


namespace Verdandi
{


    //////////////////////
    // CLAMPEDBAR MODEL //
    //////////////////////


    //! This class is a clamped-bar model.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class ClampedBar: public VerdandiBase
    {
    public:
        //! The numerical type (e.g., double).
        typedef T value_type;
        //! Pointer to the numerical type.
        typedef T* pointer;
        //! Const pointer to the numerical type.
        typedef const T* const_pointer;
        //! Reference to the numerical type.
        typedef T& reference;
        //! Const reference to the numerical type.
        typedef const T& const_reference;
#ifdef VERDANDI_SPARSE
        //! Type of the background error covariance matrix.
        typedef Matrix<T, General, RowSparse> state_error_variance;
        //! Type of a row of the background error variance.
        typedef Vector<T> state_error_variance_row;
        //! Type of the model/observation crossed matrix.
        typedef Matrix<T, General, RowSparse> matrix_state_observation;
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, RowSparse> tangent_linear_operator;
#else
        //! Type of the background error covariance matrix.
        typedef Matrix<T> state_error_variance;
        //! Type of a row of the background error variance.
        typedef Vector<T> state_error_variance_row;
        //! Type of the model/observation crossed matrix.
        typedef Matrix<T> matrix_state_observation;
        //! Type of the tangent linear operator.
        typedef Matrix<T> tangent_linear_operator;
#endif
        //! Type of the model state vector.
        typedef Vector<T> state;
        //! Collection of vector state.
        typedef Vector<state, Collection> state_collection;

    protected:
        //! Bar length.
        double bar_length_;
        //! Space step along x.
        double Delta_x_;
        //! Number of elements along x.
        int Nx_;
        //! Number of degrees of freedom (dofs).
        int Ndof_;
        //! Size of the state vector.
        int Nstate_;
        //! Time step.
        double Delta_t_;
        //! Current time.
        double time_;
        //! Simulation duration.
        double final_time_;
        //! Mass parameter.
        double mass_density_;
        //! Young's Modulus.
        double Young_modulus_;

        //! Force parameter.
        Vector<T> theta_force_;
        //! Number of force parameter regions.
        int Ntheta_force_;
        //! Force parameter region of elements.
        Vector<int> theta_force_index_;

        //! Stiffness parameter.
        Vector<T> theta_stiffness_;
        //! Number of stiffness parameter regions.
        int Ntheta_stiffness_;
        //! Stiffness parameter region of elements.
        Vector<int> theta_stiffness_index_;

        //! Damp parameter.
        Vector<T> theta_damp_;
        //! Number of damp parameter regions.
        int Ntheta_damp_;
        //! Damp parameter region of elements.
        Vector<int> theta_damp_index_;

        //! State collection.
        state_collection x_;
        //! Full state collection.
        state_collection x_full_;

        //! State lower bound.
        state lower_bound_;
        //! State upper bound.
        state upper_bound_;

        //! Mass parameter
        Vector<T> theta_mass_;
        //! Number of mass parameter regions.
        int Ntheta_mass_;
        //! Mass parameter region of elements.
        Vector<int> theta_mass_index_;

        //! FEM Vector (disp 0).
        state disp_0_;
        //! FEM Vector (velo 0).
        state velo_0_;
        //! FEM Vector (force).
        state force_;
        //! State.
        set<string> stable_;
        //! Reduced state.
        vector<string> reduced_;

        // Mass FEM matrix.
        Matrix<T, General, RowMajor> mass_FEM_matrix_;
        // Stiffness FEM matrix.
        Matrix<T, General, RowMajor> stiffness_FEM_matrix_;
        // Damp FEM matrix.
        Matrix<T, General, RowMajor> damp_FEM_matrix_;

        //! Newmark Global FEM matrix (mass matrix).
        Matrix<T, Symmetric, RowSymSparse> mass_matrix_;
        //! Newmark Global FEM matrix (Newmark matrix 0).
        Matrix<T, Symmetric, RowSymSparse> Newmark_matrix_0_;
        //! Newmark Global FEM matrix (Newmark matrix 1).
        Matrix<T, Symmetric, RowSymSparse> Newmark_matrix_1_;

        //! Damp matrix (C).
        Matrix<T, Symmetric, RowSymSparse> damp_matrix_;
        //! Damp alpha coefficient.
        double alpha_;
        //! Damp beta coefficient.
        double beta_;

        //! Stiffness matrix (K).
        Matrix<T, Symmetric, RowSymSparse> stiffness_matrix_;

        //! Background error variance.
        double state_error_variance_value_;

        //! Background error covariance matrix (B).
        state_error_variance state_error_variance_;
        //! Inverse of the background error covariance matrix (B^-1).
        state_error_variance state_error_variance_inverse_;

        //! Number of the row of B currently stored.
        int current_row_;
        //! Number of the column of Q currently stored.
        int current_column_;
        //! Value of the row of B currently stored.
        state_error_variance_row state_error_variance_row_;
        //! PI.
        static const double Pi_;

        /***  Adjoint model ***/

        //! To indicate if the adjoint variables are allocated or not.
        bool is_adjoint_initialized_;
        //! Adjoint variables;
        state_collection q_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

        //! Duplicated state.
        state duplicated_state_;

    public:
        // Constructor and destructor.
        ClampedBar();
        ClampedBar(string configuration_file);
        ~ClampedBar();
        void Initialize(string configuration_file);
        void InitializeStep();
        void InitializeAdjoint();

        // Processing.
        void Forward(bool update_force = true);
        bool HasFinished() const;
        void FinalizeStep();
        void Finalize();
        void Save();
        void BackwardAdjoint(state& state_innovation);

        // Operators.
        void ApplyOperator(state& x, bool forward = false,
                           bool preserve_state = true,
                           bool update_force = true);
        void ApplyTangentLinearOperator(state& x);
        void GetTangentLinearOperator(tangent_linear_operator&) const;

        // Access methods.
        double GetTime() const;
        void SetTime(double time);
        int GetNstate() const;
        void GetStateCopy(state& state);
        void SetStateCopy(state& state);
        state& GetState();
        void StateUpdated();
        state& GetStateLowerBound();
        state& GetStateUpperBound();
        state& GetFullState();
        void FullStateUpdated();
        void GetAdjointState(state& state_adjoint);
        void SetAdjointState(const state& state_adjoint);

        void GetStateErrorVarianceRow(int row, state_error_variance_row&
                                      state_error_variance_row);
        state_error_variance& GetStateErrorVariance();
#ifndef SWIG
        const state_error_variance& GetStateErrorVariance() const;
#endif
        void GetStateErrorVarianceSqrt(state_error_variance& L,
                                       state_error_variance& U);
        const state_error_variance& GetStateErrorVarianceInverse() const;

        string GetName() const;
        void Message(string message);

    private:
        void BuildRegionIndex(int N, int Nregion, Vector<int>& index_vector);
        void AssembleMassMatrix(Vector<T>& theta, Vector<int>& theta_index);
        void AssembleNewMarkMatrix0();
        void AssembleNewMarkMatrix1();
        void AssembleDampMatrix(Vector<T>& theta, Vector<int>& theta_index);

        void AssembleStiffnessMatrix(Vector<T>& theta,
                                     Vector<int>& theta_index);
        void AllocateSparseMatrix();
        void SetShape(state& x, state_collection& x_collection) const;
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_CLAMPEDBAR_HXX
#endif
