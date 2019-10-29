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
import ops, sys, os, ctypes, array

## This class is a linear observation manager written in Python.
class LinearObservationManager:
    ## Initializes the model.
    # @param[in] path The path to the python file
    def __init__(self, path, Nstate_model):
        configuration = ops.Ops(path)

        self.Nstate_model_ = Nstate_model;

        configuration.SetPrefix("observation.")
        self.observation_file_ = configuration.GetString("file")
        self.observation_type_ = configuration.GetString("type", "", "state")
        self.Delta_t_ = configuration.GetDouble("Delta_t", "v > 0")
        self.Nskip_ = configuration.GetInt("Nskip", "v > 0")
        self.initial_time_ = configuration.GetDouble("initial_time", "", 0.)
        self.final_time_ = configuration.GetDouble("final_time", "",
                                                   sys.float_info.max)

        self.time_ = sys.float_info.min

        self.width_file_ = configuration.GetString("width_file", "", "")
        self.error_variance_value_ \
            = configuration.GetDouble("error.variance", "v > 0")

        # Building the matrices
        self.operator_scaled_identity_ \
            = configuration.GetBool("operator.scaled_identity")
        if self.operator_scaled_identity_:
            self.operator_diagonal_value_ \
                = configuration.GetDouble("operator.diagonal_value")
            self.Nobservation_ = self.Nstate_model_
            self.tangent_operator_matrix_ = self.operator_diagonal_value_ \
                * identity(self.Nobservation_, float)
        else:
            print "The operator can only be a scaled identity matrix."
            return

        self.error_variance_ = self.error_variance_value_ \
            * identity(self.Nobservation_, float)
        self.error_variance_inverse_ = 1. / self.error_variance_value_ \
            * identity(self.Nobservation_, float)

        if self.observation_type_ == "state":
            self.Nbyte_observation_ = Nstate_model \
                * ctypes.sizeof(ctypes.c_double) \
                + ctypes.sizeof(ctypes.c_int)
        if self.observation_type_ == "observation":
            self.Nbyte_observation_ = self.Nobservation_ \
                * ctypes.sizeof(ctypes.c_double) \
                + sizeof(ctypes.c_int)

        expected_file_size = self.Nbyte_observation_ \
            * (int((self.final_time_ - self.initial_time_)
                   / (self.Delta_t_ * self.Nskip_)) + 1)

        print self.observation_file_
        actual_size = os.path.getsize(self.observation_file_)

        if expected_file_size > actual_size:
            print "Too few available observations. The size of \"", \
                self.observation_file_, "\" must be greater than ", \
                expected_file_size, " B."


    ## Sets the current time.
    # @param[in] time the current time.
    def SetTime(self, time):
        if self.time_ == time:
            return
        self.time_ = time


    ## Indicates if some observations are available at a given time.
    # @param[in] time a given time.
    def HasObservation(self, time = None):
        if time is not None:
            self.SetTime(time)
        return self.time_ < self.final_time_


    ## Returns the number of available observations.
    # @return The total number of observation at current time.
    def GetNobservation(self):
        return self.Nobservation_


    ## Returns the observations.
    # This method is called after 'SetTime' set the time at which the
    # observations are requested.
    # @return observation observation vector.
    def GetObservation(self):
        return self.ReadObservation(self.time_)


    # Builds observations associated with given time.
    # @param[in] time the given time.
    # @return the observations.
    def ReadObservation(self, time):
        f = open(self.observation_file_, "rb")
        position = int(floor((time - self.initial_time_)
                             / (self.Delta_t_ * self.Nskip_) + 0.5)
                       * self.Nbyte_observation_)
        f.seek(position)
        a = array.array('i')
        a.read(f, 1)
        size = a[0]
        obs = array.array('d')
        obs.read(f, size)
        f.close()
        return obs


    ## Applies the operator to a given vector.
    # @param[in] x a vector.
    # @param[out] y the value of the operator at \a x.
    def ApplyOperator(self, state, observation = None):
        if state.size == 0:
            return
        if self.operator_scaled_identity_:
            observation = self.operator_diagonal_value_ * state
        return observation


    ## Applies the tangent linear operator to a given vector.
    # This method is called after 'SetTime' set the time at which the
    # observations are requested.
    # @param[in] x a vector.
    # @return the value of the tangent linear operator at \a x.
    def ApplyTangentLinearOperator(self, state):
        return self.ApplyOperator(state)


    ## Gets innovation.
    # @param[in] state state vector.
    # @return innovation innovation vector.
    def GetInnovation(self, state):
        return self.GetObservation() - self.ApplyOperator(state)


     ## Returns an element of the tangent linear operator, or the full tangent
     # linear operator if no row and column index is provided.
     # This method is called after 'SetTime' set the time at which the
     # operator is defined.
     # @param[in] i row index.
     # @param[in] j column index.
     # @return The element (\a i, \a j) of the tangent linear operator.
    def GetTangentLinearOperator(self, i = None, j = None):
        if i is None and j is None:
            return self.tangent_operator_matrix_
        else:
            if i == j:
                return self.operator_diagonal_value_
            else:
                return 0.


    ## Returns a row of the tangent linear operator.
    # This method is called after 'SetTime' set the time at which the
    # operator is defined.
    # @param[in] row row index.
    # @return The row \a row of the tangent linear operator.
    def GetTangentLinearOperatorRow(self, row):
        return self.tangent_operator_matrix_[row]


    ## Applies the adjoint operator to a given vector.
    # @param[in] x a vector
    # @return the value of the operator at \a x.
    def ApplyAdjointOperator(self, state):
        if self.operator_scaled_identity_:
            return self.operator_diagonal_value_ * state


    ## Returns an element of the observation error covariance matrix, or the
    # full matrix.
    # @param[in] i row index.
    # @param[in] j column index.
    # @return The element (\a i, \a j) of the observation error covariance
    # matrix, or the full matrix.
    def GetErrorVariance(self, i = None, j = None):
        if i is None and j is None:
            return self.error_variance_
        else:
            if i == j:
                return self.error_variance_value_
            else:
                return 0


    ## Returns the inverse of the observation error covariance matrix.
    # @return The inverse of the observation error covariance matrix.
    def GetErrorVarianceInverse(self):
        return self.error_variance_inverse_


    ## Returns the name of the class.
    # @return The name of the class.
    def GetName(self):
        return "LinearObservationManager"


   ## Receives and handles a message.
    # @param[in] message the received message
    def Message(self, message):
        return
