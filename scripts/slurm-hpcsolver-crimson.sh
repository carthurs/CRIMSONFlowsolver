#!/bin/bash
# Script to launch the CRIMSON Flowsolver on SLURM
# Chris Arthurs 2016/07/25

if [[ $# -gt 1 ]]; then
	printf "\n"
	echo "Cores: $1"
	echo "Wall time: $2"
	if [[ $# -gt 2 ]]; then
		jobname=$3
	else
		jobname="CRIMSON"
	fi
	echo "Job name: ${jobname}"
	
	run="#!/bin/bash\
	 \n# MOAB/Torque submission script for SciNet GPC \
	 \n#\
	 \n#SBATCH --ntasks=${1}\
	 \n#SBATCH --time=00-${2}\
	 \n#SBATCH --job-name=${jobname}\
	 \n#SBATCH --mem-per-cpu=15360  #27531168\
	 \n\
	 \n# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from\
	 \n# Note: SLURM defaults to running jobs in the directory\
	 \n# where they are submitted, no need for $PBS_O_WORKDIR\
	 \n \
	 \n ulimit -s unlimited\
	 \n \
	 \n unset LD_LIBRARY_PATH\
	 \n unset PYTHONPATH\
	 \n \
	 \n export LD_LIBRARY_PATH="/usr/local/openmpi-intel/1.8.8/lib:/usr/local/intel/composer_xe_2015.3.187/compiler/lib/intel64:/usr/local/python/2.7.5-gcc/lib:/usr/local/boost/1.58.0-gcc/lib:/usr/local/intel/mkl/lib/intel64:/usr/local/gcc/4.9.3/lib64:/usr/local/openblas/0.2.9-gcc/lib"\
	 \n
	 \n module purge\
	 \n module load slurm/latest\
	 \n module load modules\
	 \n module load showq/0.15\
	 \n module load gmp/5.1.2\
	 \n module load mpfr/3.1.2\
	 \n module load mpc/1.0.3\
	 \n module load scons/2.5.0\
	 \n module load openmpi-intel/1.8.8\
	 \n module unload intel\
	 \n module load lapack-intel/3.5.0\
	 \n module load boost-intel/1.58.0\
	 \n module load gcc/4.9.3\
	 \n module load intel/15.3.187\
	 \n module load python-gcc/2.7.5\
	 \n module unload openmpi\
	 \n module unload openblas-gcc\
	 \n module unload R-gcc\
	 \n module load openblas-gcc/0.2.9\
	 \n \
	 \n \
	 \n export PHASTA_CONFIG="/vlsci/VR0285/jmynard/CRIMSON-Flowsolver/simvascular_flowsolver_estimator/configs"\
	 \n export CRIMSON_FLOWSOLVER_HOME=/vlsci/VR0285/jmynard/CRIMSON-Flowsolver/simvascular_flowsolver_estimator\
	 \n \
	 \n \
	 \n mpirun -np ${1} /vlsci/VR0285/jmynard/CRIMSON-Flowsolver/simvascular_flowsolver_estimator/bin/flowsolver solver.inp"

	echo -e $run
	echo -e $run | sbatch
else
	printf "\n"
	echo "CRIMSON Flowsolver SLURM Submission Script"
	echo "Version 2016.7.25 - christopher.arthurs@kcl.ac.uk"
	printf "\n"
	echo "Usage: hpcsolver-crimson.sh <processors> <wall time> <job name>"
	printf "\n"
fi
