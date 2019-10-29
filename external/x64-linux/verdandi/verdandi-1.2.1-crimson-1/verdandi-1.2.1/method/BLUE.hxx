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


#ifndef VERDANDI_FILE_METHOD_BLUE_HXX


namespace Verdandi
{


    template <class Model, class ObservationManager,
              class Innovation, class State>
    void ComputeBLUE_vector(Model& model,
                            ObservationManager& observation_manager,
                            const Innovation& innovation, State& state);


    template <class StateErrorVariance, class ObservationOperator,
              class Observation, class ObservationErrorVariance,
              class State>
    void ComputeBLUE_matrix(StateErrorVariance& B,
                            const ObservationOperator& H,
                            const Observation& y,
                            const ObservationErrorVariance& R,
                            State& x,
                            bool is_y_innovation = false,
                            bool compute_variance = false);


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
                            bool is_y_innovation = false,
                            bool compute_variance = false);


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_BLUE_HXX
#endif
