# Copyright (C) 2010, INRIA
# Author(s): Claire Mouton, Vivien Mallet
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


# Run this example either with "python example.py" or
# from ipython and then "run -i example.py".


import seldon
import verdandi


#############################################################################


print "*** With the forward driver"

configuration_file = "configuration/truth.lua"

print "Configuration file:", configuration_file

forward_driver = verdandi.Method()

forward_driver.Initialize(configuration_file)
forward_driver.InitializeStep()

model = forward_driver.GetModel()
print "Initial time:", model.GetTime()

# 'Forward()' can be processed either by calling the driver method or directly
# by calling the model method.
forward_driver.Forward()
print "Time after one time step:", model.GetTime()

forward_driver.InitializeStep()
model.Forward()
print "Time after another time step:", model.GetTime()

# Time loop.
for i in range(18):
    forward_driver.InitializeStep()
    forward_driver.Forward()

# In order to perform the whole simulation.
while not forward_driver.HasFinished():
    forward_driver.InitializeStep()
    forward_driver.Forward()


#############################################################################


print
print "*** With the optimal interpolation method"

configuration_file = "configuration/assimilation.lua"

print "Configuration file", configuration_file

oi = verdandi.Method1()
model = oi.GetModel()

oi.Initialize(configuration_file)

# One step forward.
oi.InitializeStep()
oi.Forward()
oi.Analyze()

# Another step forward.
oi.InitializeStep()
oi.Forward()

# Computes the analysis and monitors the change in the state vector.
print "State vector before analysis:"
state_vector = seldon.VectorDouble()
state_vector = model.GetState()
state_vector.Print()

# The analysis is now performed (if observations are available).
oi.Analyze()

print "State vector after analysis:"
state_vector.Print()

# Time loop (similar to that of the forward driver).
for i in range(19):
    oi.InitializeStep()
    oi.Forward()
    oi.Analyze()

# Note that the state may be modified directly in Python.
state_vector[0] = 0 # just changes the first component.
model.StateUpdated() # The state has been updated.

# And the simulation can go on.
while not oi.HasFinished():
    oi.InitializeStep()
    oi.Forward()
    oi.Analyze()
