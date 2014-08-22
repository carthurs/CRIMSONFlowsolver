#!/bin/bash
ulimit -s unlimited

export LES_LICENSE_SERVER="blackbook2v2"
export PHASTA_CONFIG="/home/klau/dev/simvascular_flowsolver_estimator/configs"

/home/klau/dev/simvascular_flowsolver_estimator/bin/flowsolver $1 $2 $3

