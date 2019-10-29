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


#ifndef VERDANDI_FILE_HAMILTONJACOBIBELLMAN_HXX


namespace Verdandi
{


    ///////////////////////////
    // HAMILTONJACOBIBELLMAN //
    ///////////////////////////


    //! This class is a solver for Hamilton-Jacobi-Bellman equation.
    template <class T, class Model, class ObservationManager>
    class HamiltonJacobiBellman: public VerdandiBase
    {

    public:
        typedef typename Model::state model_state;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;

        //! Observation manager.
        ObservationManager observation_manager_;

        /*** Configuration ***/

        //! Display options.
        map<string, bool> option_display_;

        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

        /*** Equation coefficients ***/

        //! Should the quadratic term be taken into account?
        bool with_quadratic_term_;
        //! Should the advection term be taken into account?
        bool with_advection_term_;
        //! Should the source term be taken into account?
        bool with_source_term_;

        //! Q_0.
        Matrix<T> Q_0_;
        //! Position of the minimum of the initial parabola.
        Vector<T> x_0_;
        //! Inverse of Q.
        Matrix<T> Q_inv_;
        //! R.
        Matrix<T> R_;

        /*** Value function ***/

        //! Number of dimensions.
        int Ndimension_;

        //! Total number of points in the domain.
        int Npoint_;

        //! First coordinate along each dimension.
        Vector<T> x_min_;
        //! Space step in each dimension.
        Vector<T> Delta_x_;
        //! Number of points in each dimension.
        Vector<int> Nx_;

        //! Number of time steps.
        int Nt_;
        //! Initial time.
        T initial_time_;
        //! Time step.
        T Delta_t_;
        //! Current iteration.
        int time_step_;

        //! Value function V(t, x).
        Vector<T> V_;

        /*** Numerical scheme ***/

        /*! Name of the numerical scheme: LxF (for first-order
          Lax-Friedrichs), BrysonLevy (first-order central scheme) or Godunov
          (first-order). */
        string scheme_;

        //! Location of the evolution points in Bryson-Levy scheme.
        Vector<T> a_Delta_x_;

        // Is the model time-dependent?
        bool model_time_dependent_;
        //! M(x) in every point in space.
        Matrix<T> Mx_;
        //! Courant number.
        T courant_number_;

        /*! Type of the boundary condition ('Dirichlet', 'Extrapolation' and
          'Period' are supported). */
        string boundary_condition_type_;
        /*! 0 for Dirichlet boundary conditions, 1 for the extrapolation, and
          2 for periodic. */
        int boundary_condition_index_;
        /*! Boundary condition (constant value imposed outside the domain) in
          case of Dirichlet boundary conditions. */
        T boundary_condition_;

        /*! Upper bound on the absolute value of the model operator, in every
          direction. Useful in Lax-Friedrichs scheme. */
        Vector<T> upper_bound_model_;

    public:

        /*** Constructor and destructor ***/

        HamiltonJacobiBellman();
        ~HamiltonJacobiBellman();

        /*** Methods ***/

        void Initialize(string configuration_file);
        void Initialize(VerdandiOps& configuration);

        void InitializeStep();

        void Forward();
        void AdvectionLxFForward();
        void AdvectionBrysonLevyForward();
        void AdvectionGodunov();

        void FinalizeStep();
        void Finalize();

        T GodunovFlux(T q, T M, T v_l, T v, T v_r) const;

        bool HasFinished();

        // Access methods.
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_HAMILTONJACOBIBELLMAN_HXX
#endif
