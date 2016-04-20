#!/bin/bash

ulimit -s unlimited

export LES_LICENSE_SERVER="blackbook2v2"
export PHASTA_CONFIG="/home/ADD/PATH/estimator/configs"
# add path to boost libraries to LD_LIBRARY_PATH
export LD_LIBRARY_PATH="/home/ADD/PATH/boost_1_57_0/install/lib":$LD_LIBRARY_PATH
export CRIMSON_FLOWSOLVER_HOME="/home/ADD/PATH/estimator"

# on FLUX/SGI you will need to remove mpirun and call that in the job script instead
mpirun -np $1 /home/ADD/PATH/estimator/bin/flowsolver $2 $3
