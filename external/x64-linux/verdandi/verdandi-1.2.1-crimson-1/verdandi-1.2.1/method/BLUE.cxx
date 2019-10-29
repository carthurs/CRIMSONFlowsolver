// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu, Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_METHOD_BLUE_CXX


#include "BLUE.hxx"
#include "seldon/computation/solver/SparseSolver.cxx"


namespace Verdandi
{


    //! Computes BLUE.
    /*! It computes the BLUE (best linear unbiased estimator).
      \param[in] model the model.
      \param[in] observation_manager the observation manager.
      \param[in] innovation the innovation vector.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
    */
    template <class Model, class ObservationManager,
              class Innovation, class State>
    void ComputeBLUE_vector(Model& model,
                            ObservationManager& observation_manager,
                            const Innovation& innovation, State& state)
    {
        typedef typename State::value_type T;

        int Nobservation, Nstate;
        Nobservation = observation_manager.GetNobservation();
        Nstate = model.GetNstate();
        int Nlocal_state = Nstate;

        if (Nobservation == 0) // No observations.
            return;

        int r, c;

        int global_state_number = 0;

#if defined(VERDANDI_WITH_MPI)
        int rank = MPI::COMM_WORLD.Get_rank();
        int Nprocess = MPI::COMM_WORLD.Get_size();
        Nlocal_state = Nstate / Nprocess;
        if (rank < Nstate % Nprocess)
            Nlocal_state ++;
        int div = Nstate % Nprocess;
        if (rank < div)
            global_state_number = rank * Nlocal_state;
        else
            global_state_number = div * (Nlocal_state + 1)
                + (rank - div) * Nlocal_state;
#endif

        // One row of background matrix B.
        typename Model::state_error_variance_row
            state_error_variance_row(Nstate);

        // One row of tangent operator matrix.
        typename ObservationManager::tangent_linear_operator_row
            tangent_operator_row(Nstate);

        // Temporary matrix and vector.
        // 'HBHR_inv' will eventually contain the matrix (HBH' + R)^(-1).
        Matrix<T> HBHR_inv(Nobservation, Nobservation);
        HBHR_inv.Fill(T(0));

        Vector<T> row(Nobservation);

        // Computes HBH'.
        T H_entry;
        for (int j = 0; j < Nlocal_state; j++)
        {
            model.GetStateErrorVarianceRow(j + global_state_number,
                                           state_error_variance_row);
            // Computes the j-th row of BH'.
            for (r = 0; r < Nobservation; r++)
            {
                observation_manager
                    .GetTangentLinearOperatorRow(r, tangent_operator_row);
                row(r) = DotProd(state_error_variance_row,
                                 tangent_operator_row);
            }

            // Keeps on building HBH'.
            for (r = 0; r < Nobservation; r++)
            {
                H_entry = observation_manager.
                    GetTangentLinearOperator(r, j + global_state_number);
                for (c = 0; c < Nobservation; c++)
                    HBHR_inv(r, c) += H_entry * row(c);
            }
        }

#if defined(VERDANDI_WITH_MPI)
	Matrix<T> HBHR_recv(Nobservation, Nobservation);
        MPI::COMM_WORLD.Allreduce(HBHR_inv.GetData(),
                                  HBHR_recv.GetData(),
                                  Nobservation * Nobservation,
                                  MPI::DOUBLE,
                                  MPI::SUM);
	HBHR_inv = HBHR_recv;
#endif

        // Computes (HBH' + R).
        for (r = 0; r < Nobservation; r++)
            for (c = 0; c < Nobservation; c++)
                HBHR_inv(r, c) += observation_manager.GetErrorVariance(r, c);

        // Computes (HBH' + R)^{-1}.
        GetInverse(HBHR_inv);

        // Computes HBHR_inv * innovation.
        Vector<T> HBHR_inv_innovation(Nobservation);
        MltAdd(T(1), HBHR_inv, innovation, T(0), HBHR_inv_innovation);

        // Computes new state.
        Vector<T> BHt_row(Nobservation);
#if defined(VERDANDI_WITH_MPI)
        typename Model::state state_update_send(Nstate);
        state_update_send.Fill(T(0.));
#endif
        BHt_row.Fill(T(0));
        for (r = 0; r < Nlocal_state; r++)
        {
            // Computes the r-th row of BH'.
            model.GetStateErrorVarianceRow(r + global_state_number,
                                           state_error_variance_row);
            for (c = 0; c < Nobservation; c++)
            {
                observation_manager
                    .GetTangentLinearOperatorRow(c, tangent_operator_row);
                BHt_row(c) = DotProd(state_error_variance_row,
                                     tangent_operator_row);
            }

#if defined(VERDANDI_WITH_MPI)
            state_update_send(r + global_state_number)
                += DotProd(BHt_row, HBHR_inv_innovation);
#else
            state(r) += DotProd(BHt_row, HBHR_inv_innovation);
#endif
        }

#if defined(VERDANDI_WITH_MPI)
        typename Model::state state_update_recv(Nstate);
        state_update_recv.Fill(T(0.));
        MPI::COMM_WORLD.Allreduce(state_update_send.GetDataVoid(),
                                  state_update_recv.GetDataVoid(),
                                  Nstate, MPI::DOUBLE, MPI::SUM);
        Add(T(1.), state_update_recv, state);
#endif
    }


    //! Computes BLUE using operations on matrices.
    /*! This method is mainly intended for cases where the covariance matrices
      are sparse matrices. Otherwise, the manipulation of the matrices may
      lead to unreasonable memory requirements and to high computational
      costs.
      \param[in] B error variance associated with \a state.
      \param[in] H observation operator.
      \param[in] y the vector of observations or innovations.
      \param[in] R error variance associated with \a observation.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
      \param[in] is_y_innovation Boolean to indicate if the parameter \a y is
      a vector of observations or innovations.
      \param[in] compute_variance Boolean to indicate if the covariance matrix
      has to be updated.
    */
    template <class StateErrorVariance, class ObservationOperator,
              class Observation, class ObservationErrorVariance,
              class State>
    void ComputeBLUE_matrix(StateErrorVariance& B,
                            const ObservationOperator& H,
                            const Observation& y,
                            const ObservationErrorVariance& R,
                            State& x,
                            bool is_y_innovation,
                            bool compute_variance)
    {
        ComputeBLUE_matrix(B, H, H, y, R, x, is_y_innovation,
                           compute_variance);
    }


    //! Computes BLUE using operations on matrices.
    /*! This method is mainly intended for cases where the covariance matrices
      are sparse matrices. Otherwise, the manipulation of the matrices may
      lead to unreasonable memory requirements and to high computational
      costs.
      \param[in] B error variance associated with \a state.
      \param[in] H observation operator.
      \param[in] cm this parameter is only used to determine the type of an
      intermediate matrix in the computations. \a cm is not modified nor
      read. Its type CrossedMatrix will be the type of the intermediate matrix
      BH', whose size is Nx times Ny, if Nx is the length of \a x and Ny is
      the length of \a y.
      \param[in] y the vector of observations or innovations.
      \param[in] R error variance associated with \a observation.
      \param[in,out] x on entry, the background vector; on exit, the analysis.
      \param[in] is_y_innovation Boolean to indicate if the parameter \a y is
      a vector of observations or innovations.
      \param[in] compute_variance Boolean to indicate if the covariance matrix
      has to be updated.
    */
    template <class StateErrorVariance,
              class ObservationOperator, class MatrixStateObservation,
              class Observation, class ObservationErrorVariance,
              class State>
    void ComputeBLUE_matrix(StateErrorVariance& B,
                            const ObservationOperator& H,
                            const MatrixStateObservation& cm,
                            const Observation& y,
                            const ObservationErrorVariance& R,
                            State& x,
                            bool is_y_innovation,
                            bool compute_variance)
    {
        typedef typename State::value_type T;

        int Ny = y.GetLength();
        int Nx = x.GetLength();

        if (B.GetM() != Nx || B.GetN() != Nx)
            throw ErrorArgument("ComputeBLUE_matrix",
                                "The state error variance has "
                                "dimensions " + to_str(B.GetM())
                                + " x " + to_str(B.GetN())
                                + " while the state vector has "
                                + to_str(Nx) + " elements.");

        if (H.GetN() != Nx)
            throw ErrorArgument("ComputeBLUE_matrix",
                                "The observation operator has "
                                + to_str(H.GetN()) + " columns "
                                "while the state vector has "
                                + to_str(Nx) + " elements.");

        if (H.GetM() != Ny)
            throw ErrorArgument("ComputeBLUE_matrix",
                                "The observation operator has "
                                + to_str(H.GetM()) + " rows "
                                "while the observation vector has "
                                + to_str(Ny) + " elements.");

        if (R.GetM() != Ny || R.GetN() != Ny)
            throw ErrorArgument("ComputeBLUE_matrix",
                                "The observation error variance has "
                                "dimensions " + to_str(R.GetM())
                                + " x " + to_str(R.GetN())
                                + " while the observation vector has "
                                + to_str(Ny) + " elements.");

        if (Ny == 0) // No observations.
            return;

        // Temporary matrices.
        MatrixStateObservation working_matrix_xy(Nx, Ny);

        ObservationErrorVariance working_matrix_yy(Ny, Ny);

        // Computes BH'.
        MltAdd(T(1), SeldonNoTrans, B, SeldonTrans, H, T(0),
               working_matrix_xy);

        // Computes HBH'.
        Mlt(H, working_matrix_xy, working_matrix_yy);

        // Computes (HBH' + R).
        Add(T(1), R, working_matrix_yy);

        // Innovation.
        Observation innovation;
        if (is_y_innovation)
            innovation.SetData(y);
        else
        {
            innovation = y;
            MltAdd(T(-1), H, x, T(1), innovation);
        }

        if (!compute_variance)
        {
            // Computes inc = (HBH' + R)^{-1} * innovation by solving the
            // linear system (HBH' + R) * inc = innovation.
            GetAndSolveLU(working_matrix_yy, innovation);
            MltAdd(T(1), working_matrix_xy, innovation, T(1), x);
        }
        else
        {
            // Kalman Gain.
            MatrixStateObservation K(Nx, Ny);
            K.Fill(T(0));

            GetInverse(working_matrix_yy);

            MltAdd(T(1), working_matrix_xy, working_matrix_yy, T(0), K);

            MltAdd(T(1), K, innovation, T(1), x);

            MltAdd(T(-1), SeldonNoTrans, K, SeldonTrans,
                   working_matrix_xy, T(1), B);
        }

        if (is_y_innovation)
            innovation.Nullify();
    }

} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BLUE_CXX
#endif
