script_raw = """
#!/bin/bash

ulimit -s unlimited

export PHASTA_CONFIG="%s/configs"
export CRIMSON_FLOWSOLVER_HOME="%s"

mpirun -np $1 $CRIMSON_FLOWSOLVER_HOME/bin/flowsolver $2 $3
"""