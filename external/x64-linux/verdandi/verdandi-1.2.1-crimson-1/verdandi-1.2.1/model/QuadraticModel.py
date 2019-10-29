# Copyright (C) 2011-2012, INRIA
# Author(s): Kevin Charpentier, Vivien Mallet
#
# This file is an example part of the data assimilation library Verdandi.
#
# Verdandi is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Verdandi. If not, see http://www.gnu.org/licenses/.
#
# For more information, visit the Verdandi web site:
#      http://verdandi.gforge.inria.fr/

from numpy import *
import ops

## This class is a quadratic model written in Python.
class QuadraticModel:
    ## Initializes the model.
    # @param[in] path The path to the python file
    def __init__(self, path):
        configuration = ops.Ops(path)

        configuration.SetPrefix("quadratic_model.")

        self.state_ = array(
            configuration.GetVectDouble("definition.initial_state"))

        self.with_quadratic_term_ \
            = configuration.GetBool("definition.with_quadratic_term")
        self.with_linear_term_ \
            = configuration.GetBool("definition.with_linear_term")
        self.with_constant_term_ \
            = configuration.GetBool("definition.with_constant_term")

        # Dimension of the state.
        self.Nstate_ = self.state_.size

        # Quadratic terms.
        self.S_ \
            = array(configuration.GetVectDouble("definition.quadratic_term"))
        self.S_ = self.S_.reshape(self.Nstate_, self.Nstate_, self.Nstate_)

        # Matrix that defines the linear part of the model.
        self.L_ = array(configuration.GetVectDouble("definition.linear_term"))
        self.L_ = self.L_.reshape(self.Nstate_, self.Nstate_)

        # Vector that defines the constant part of the model.
        self.b_ = array(configuration.GetVectDouble("definition.constant"))

        # Time step.
        self.Delta_t_ = configuration.GetDouble("definition.Delta_t")

        # Current time.
        self.time_ = configuration.GetDouble("definition.initial_time")

        # Final time of the simulation.
        self.final_time_ = configuration.GetDouble("definition.final_time")

        # Errors
        if configuration.Exists("error"):
            if configuration.GetBool("error.scaled_identity"):
                self.Q_ = configuration.GetDouble("error.diagonal_value",
                                                  "v >= 0")
                self.Q_ = self.Q_ * identity(self.Nstate_, float)
            else:
                self.Q_ = array(configuration.GetVectDouble("error.value"))
                self.Q_ = self.Q_.reshape(self.Nstate_, self.Nstate_)

        if configuration.Exists("error_sqrt"):
            self.Q_sqrt \
                = array(configuration.GetVectDouble("error_sqrt.value"))
            self.Q_sqrt = self.Q_sqrt.reshape(self.Nstate_, 0)

        if configuration.Exists("state_error"):
            if configuration.GetBool("state_error.scaled_identity"):
                self.P_ \
                    = configuration.GetDouble("state_error.diagonal_value",
                                              "v >= 0")
                self.P_ = self.P_ * identity(self.Nstate_, float)
            else:
                self.P_ \
                    = array(configuration.GetVectDouble("state_error.value"))
                self.P_ = self.P_.reshape(self.Nstate_, self.Nstate_)

        if configuration.Exists("state_error_sqrt"):
            self.P_sqrt \
                = array(configuration.GetVectDouble("state_error_sqrt.value"))
            self.P_sqrt = self.P_sqrt.reshape(self.Nstate_, 0)


    ## Initializes the current time step for the model.
    def InitializeStep(self):
        return None


    ### Processing ###


    ## Advances one step forward in time.
    # \f[x^f_{h+1} = \mathcal{M}_h(x^a_h, p_h)\,.\f]
    def Forward(self):
        if self.with_quadratic_term_:
            current_state = self.state_.copy()
            for i in range(0, self.Nstate_):
                S_state_ = dot(self.S_[i], current_state)
                self.state_[i] += self.Delta_t_ * dot(S_state_, current_state)

            if self.with_linear_term_:
                self.state_ += self.Delta_t_ * dot(self.L_, current_state)

        elif self.with_linear_term_:
            self.state_ += self.Delta_t_ * dot(self.L_, self.state_)

        if self.with_constant_term_:
            self.state_ += self.Delta_t_ * self.b_

        self.time_ += self.Delta_t_


    ## Applies the model to a given vector.
    # The current state of the model is modified.
    # @param[in, out] state a vector.
    # @param[in] forward Boolean to indicate if the model has to go on to the
    #  next step.
    # @param[in] preserve_state Boolean to indicate if the model state has to
    #  be preserved.
    def ApplyOperator(self, state, forward = False, preserve_state = True):
        if preserve_state:
            current_state = self.state_.copy()
        self.state_ = state
        self.Forward()
        if not forward:
            self.time_ -= self.Delta_t_
        state = self.state_
        if preserve_state:
            self.state_ = current_state


    ## Applies the tangent linear model to a given vector.
    # @param[in,out] x on entry, a vector to which the tangent linear model
    #  should be applied; on exit, the result.
    def ApplyTangentLinearOperator(self, x):
        input = x.copy()
        if self.with_quadratic_term_:
            for i in range(0, self.Nstate_):
                S_state_ = self.Delta_t_ * dot(self.S_[i], self.state_)
                S_state_ += self.Delta_t_ \
                    * dot(self.S_[i].transpose(), self.state_)
                x[i] += dot(S_state_, input)

        if self.with_linear_term_:
            x += self.Delta_t_ * dot(self.L_, input)

        elif self.with_linear_term_:
            x += self.Delta_t_ * dot(self.L_, input)


    ## Returns the tangent linear model.
    # @param[out] M the matrix of the tangent linear model.
    def GetTangentLinearOperator(self):
        M = identity(self.Nstate_, float)
        if self.with_quadratic_term_:
            for i in range(0, self.Nstate_):
                M_row = self.Delta_t_ * dot(self.S_[i], self.state_)
                M_row += self.Delta_t_ \
                    * dot(self.S_[i].transpose(), self.state_)
                M[i] = M_row
                M[i, i] += 1.

            if self.with_linear_term_:
                M += self.Delta_t_ * self.L_

        elif self.with_linear_term_:
            M = self.L_.copy()
            M *= self.Delta_t_
            for i in range(0, self.Nstate_):
                M[i, i] += 1.

        return M


    ## Checks whether the model has finished.
    # @return True if the simulation is done, false otherwise.
    def HasFinished(self):
        return self.time_ >= self.final_time_


    ## Saves the simulated data.
    def Save(self):
        return

    ## Finalizes the current time step for the model.
    def FinalizeStep(self):
        return

    ## Finalizes the model.
    def Finalize(self):
        return


    ### Access methods ###


    ## Returns the time step.
    # @return The time step.
    def GetDelta_t(self):
        return self.Delta_t_


    ## Returns the current time.
    # @return The current time.
    def GetTime(self):
        return self.time_


    ## Sets the current time.
    # @param[in] time the current time.
    def SetTime(self, time):
        self.time_ = time


    ## Returns the dimension of the state.
    # @return The dimension of the state.
    def GetNstate(self):
        return self.Nstate_


    ## Returns the dimension of the full state.
    # @return The dimension of the full state.
    def GetNFullstate(self):
        return self.Nstate_


    ## Provides the controlled state vector.
    # @return The controlled state vector.
    def GetState(self):
        return self.state_


    ## Performs some calculations when the update of the model state is done.
    def StateUpdated(self):
        return


    ## Provides the full state vector.
    # @return The full state vector.
    def GetFullState(self):
        return self.state_


    ## Performs some calculations when the update of the model state is done.
    def FullStateUpdated(self):
        return


    ### Errors ###


    ## Returns the model error variance.
    # @return The model error variance.
    def GetErrorVariance(self):
        return self.Q_


    ## Returns the square root of the model error variance.
    # @return The square root of the model error variance.
    def GetErrorVarianceSqrt(self):
        return self.Q_sqrt_


    ## Returns the state error variance.
    # @return The state error variance.
    def GetStateErrorVariance(self):
        return self.P_


    ## Returns a row of the state error variance.
    # @param[in] row row index.
    # @return The row with index \a row in the state error variance.
    def GetStateErrorVarianceRow(self, row):
        return self.P_[row]


    ## Returns the name of the class.
    # @return The name of the class.
    def GetName(self):
        return "QuadraticModel"


    ## Receives and handles a message.
    # @param[in] message the received message
    def Message(self, message):
        return
