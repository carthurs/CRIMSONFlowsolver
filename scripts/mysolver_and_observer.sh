#!/bin/bash
ulimit -s unlimited

export LES_LICENSE_SERVER="sgiuv"
export PHASTA_CONFIG="/home/klau/dev/simvascular_flowsolver_estimator/configs"

/home/klau/dev/simvascular_flowsolver_estimator/bin/flowsolver_and_observer $1 $2 $3

