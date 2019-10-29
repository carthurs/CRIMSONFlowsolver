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


#ifndef VERDANDI_FILE_MODEL_PETSCCLAMPEDBAR_HXX

#include "seldon/vector/VectorCollection.hxx"
#include "OutputSaver.hxx"


namespace Verdandi
{


    ///////////////////////////
    // PETSCCLAMPEDBAR MODEL //
    ///////////////////////////


    //! This class is a clamped-bar model.
    /*!
      \tparam T the type of floating-point numbers.
    */
    template <class T>
    class PetscClampedBar: public VerdandiBase
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

        //! Type of the background error covariance matrix.
        typedef Matrix<T, General, PETScMPIDense> state_error_variance;
        //! Type of a row of the background error variance.
        typedef Vector<T, PETScPar> state_error_variance_row;
        //! Type of the model/observation crossed matrix.
        typedef Matrix<T, General, PETScMPIDense> matrix_state_observation;
        //! Type of the tangent linear operator.
        typedef Matrix<T, General, PETScMPIDense> tangent_linear_operator;

        //! Type of the model state vector.
        typedef Vector<T, PETScPar> state;
        //! Type of the model state vector.
        typedef Vector<T> parameter;
        //! Collection of vector state.
        typedef Vector<parameter, Collection> parameter_collection;

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
        //! Size of the reduced state vector.
        int Nreduced_;
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
        parameter theta_force_;
        //! Number of force parameter regions.
        int Ntheta_force_;
        //! Force parameter region of elements.
        Vector<int> theta_force_index_;
        //! Stiffness parameter.
        parameter theta_stiffness_;
        //! Number of stiffness parameter regions.
        int Ntheta_stiffness_;
        //! Stiffness parameter region of elements.
        Vector<int> theta_stiffness_index_;
        //! Damp parameter.
        parameter theta_damp_;
        //! Number of damp parameter regions.
        int Ntheta_damp_;
        //! Damp parameter region of elements.
        Vector<int> theta_damp_index_;
        //! Mass parameter
        parameter theta_mass_;
        //! Number of mass parameter regions.
        int Ntheta_mass_;
        //! Mass parameter region of elements.
        Vector<int> theta_mass_index_;

        //! Parameter collection.
        parameter_collection parameter_;
        //! Full parameter collection.
        parameter_collection all_parameter_;

        //! State.
        set<string> stable_;
        //! Reduced state.
        vector<string> reduced_;

        //! Displacement.
        state displacement_0_;
        //! Velocity.
        state velocity_0_;
        //! Force.
        state rhs_;
        //! Local size of state vector.
        int Nstate_local_;
        //! Petsc state.
        state state_;

        // Mass FEM matrix.
        Matrix<T, General, RowMajor> mass_FEM_matrix_;
        // Stiffness FEM matrix.
        Matrix<T, General, RowMajor> stiffness_FEM_matrix_;
        // Damp FEM matrix.
        Matrix<T, General, RowMajor> damp_FEM_matrix_;

        //! The MPI communicator to use.
        MPI_Comm mpi_communicator_;
        //! Process rank.
        int rank_;
        //! Process rank.
        int Nprocess_;

        //! Mass matrix.
        Matrix<T, General, PETScMPIAIJ> mass_;
        //! Newmark matrix 0.
        Matrix<T, General, PETScMPIAIJ> newmark_0_;
        bool newmark_0_assembled_;

        //! Newmark matrix 1.
        Matrix<T, General, PETScMPIAIJ> newmark_1_;
        bool newmark_1_assembled_;
        //! Damp matrix.
        Matrix<T, General, PETScMPIAIJ> damp_;
        //! Stiffness matrix.
        Matrix<T, General, PETScMPIAIJ> stiffness_;

        //! Damp alpha coefficient.
        double alpha_;
        //! Damp beta coefficient.
        double beta_;

        //! Background error variance.
        double state_error_variance_value_;

        //! Background error covariance matrix (B).
        state_error_variance state_error_variance_;
        //! Inverse of the background error covariance matrix (B^-1).
        state_error_variance state_error_variance_inverse_;

        static const double Pi_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

    public:
        // Constructor and destructor.
        PetscClampedBar();
        PetscClampedBar(string configuration_file);
        ~PetscClampedBar();
        void Initialize(string configuration_file);
        void Finalize();
        void InitializeFirstStep();
        void InitializeStep();
        void InitializeAdjoint();

        // Processing.
        void Forward(bool update_force = true);
        bool HasFinished() const;
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
        int GetLocalNstate() const;
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
        template <class L_matrix, class U_matrix>
        void GetStateErrorVarianceSqrt(L_matrix& L,
                                       U_matrix& U);
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
        template <class MatrixSparse, class Vector1>
        void PetscGmres(MatrixSparse& A, Vector1& x, const Vector1& b);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_MODEL_PETSCCLAMPEDBAR_HXX
#endif
