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


## This class is a template observation manager in Python.
class PythonObservationManagerTemplate:
    ## Initializes the observation manager.
    # @param[in] path The path to the python file
    # @param[in] Nstate_model The size of the state vector model
    def __init__(self, path, Nstate_model):
        return


    ## Activates or deactivates the option 'discard_observation'.
    # @param[in] discard_observation if set to true, each observation will be
    # used at most one time.
    def DiscardObservation(discard_observation):
        return


    ## Sets the time of observations to be loaded.
    # @param[in] time The current time.
    def SetTime(self, time):
        return


    ### Observation ###


    ## Returns the observations.
    # This method is called after 'SetTime' set the time at which the
    # observations are requested.
    # @return observation observation vector.
    def GetObservation(self):
        return observation


    ### Innovation ###


    ## Returns an innovation
    # This method is called after 'SetTime' set the time at which the
    # innovation is requested.
    #  @param[in] state state vector.
    #  @param[out] innovation innovation vector.
    def BackwardAdjoint(self, state):
        return innovation


    ### Access ###


    ## Indicates if some observations are available at a given time or at
    # current time.
    # @param[in] time a given time.
    def HasObservation(self, time = None):
        return hasobservation


    ## Returns the number of available observations.
    # @return The total number of observation at current time.
    def GetNobservation(self):
        return Nobservation


    ### Operators ###


    ## Applies the observation operator to a given vector.
    # This method is called after 'SetTime' set the time at which the
    # operator is defined.
    # @param[in] x a vector.
    # @param[out] y the value of the operator applied to x.
    def ApplyOperator(self, state, y = None):
        return y


    ## Applies the tangent linear operator to a given vector.
    # This method is called after 'SetTime' set the time at which the
    # operator is defined.
    # @param[in] x a vector to which the tangent linear operator
    #  should be applied
    # @return The value of the operator applied to x.
    def ApplyTangentLinearOperator(self, x):
        return


    ## Returns an element of the tangent linear operator, or the full tangent
    # linear operator if no row and column index is provided.
    # This method is called after 'SetTime' set the time at which the
    # operator is defined.
    # @param[in] i row index.
    # @param[in] j column index.
    # @return The element (\a i, \a j) of the tangent linear operator.
    def GetTangentLinearOperator(self, i = None, j = None)
        return


    ## Returns a row of the tangent linear operator.
    # This method is called after 'SetTime' set the time at which the
    # operator is defined.
    # @param[in] row row index.
    # @return The row \a row of the tangent linear operator.
    def GetTangentLinearOperatorRow(self, row)
        return


    ## Applies the adjoint operator to a given vector.
    # This method is called after 'SetTime' set the time at which the
    # operator is defined.
    # @param[in] x a vector.
    # @return The value of the operator applied to x.
    def ApplyAdjointOperator(self, x):
        return


    ### Errors ###


    ## Returns an element of the observation error covariance matrix, or the
    # full matrix.
    # @param[in] i row index.
    # @param[in] j column index.
    # @return The element (\a i, \a j) of the observation error covariance
    # matrix, or the full matrix.
    def GetErrorVariance(self, i = None, j = None):
        return


    ## Returns the inverse of the observation error covariance matrix.
    # @return The inverse of the observation error covariance matrix.
    def GetErrorVarianceInverse(self):
        return


    ## Returns the name of the class.
    # @return The name of the class.
    def GetName(self):
        return "PythonObservationManagerTemplate"


    ## Receives and handles a message.
    # @param[in] message the received message
    def Message(self, message):
        return
