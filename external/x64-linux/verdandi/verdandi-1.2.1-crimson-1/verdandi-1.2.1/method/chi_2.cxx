// Copyright (C) 2010 INRIA
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


#ifndef VERDANDI_FILE_METHOD_CHI_2_CXX


#include "chi_2.hxx"
#include "seldon/computation/solver/SparseSolver.cxx"


namespace Verdandi
{


    /*! \brief Computes \f$ \left(y - H x\right)^{\rm T} \left(R + H B H^{\rm
      T} \right)^{-1} \left(y - H x\right) \f$ */
    /*! This method can help to check the consistency between an innovation
      and its statistics.
      \param[in] x the state vector.
      \param[in] B error variance associated with \a x.
      \param[in] H observation operator.
      \param[in] y the vector of observations.
      \param[in] R error variance associated with \a y.
      \return \f$ \left(y - H x\right)^{\rm T} \left(R + H B H^{\rm T}
      \right)^{-1} \left(y - H x\right) \f$
    */
    template <class StateErrorVariance, class ObservationOperator,
              class ObservationVector, class ObservationErrorVariance,
              class StateVector>
    typename StateVector::value_type
    chi_2(const StateVector& x,
          const StateErrorVariance& B,
          const ObservationOperator& H,
          const ObservationVector& y,
          const ObservationErrorVariance& R)
    {
        typedef typename StateVector::value_type T;

        int Ny = y.GetLength();
        int Nx = x.GetLength();

        if (Ny == 0) // No observations.
            throw ErrorArgument("chi_2(x, B, H, y, R)",
                                "The observation vector 'y' is empty.");

        // Temporary matrix and vector.
        ObservationOperator working_matrix_xy(Nx, Ny);

        ObservationErrorVariance working_matrix_yy(Ny, Ny);

        // Computes BH'.
        MltAdd(T(1), SeldonNoTrans, B, SeldonTrans, H, T(0),
               working_matrix_xy);

        // Computes HBH'.
        Mlt(H, working_matrix_xy, working_matrix_yy);

        // Computes (HBH' + R).
        Add(T(1), R, working_matrix_yy);

        // Computes inc = (HBH' + R)^{-1} * innovation by solving the linear
        // system (HBH' + R) * inc = innovation.
        ObservationVector increment = y;
        MltAdd(T(-1), H, x, T(1), increment);
        ObservationVector innovation = increment;
        GetAndSolveLU(working_matrix_yy, increment);
        return DotProd(innovation, increment);
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_CHI_2_CXX
#endif
