# Installation Instructions

These instructions explain how to build CRIMSON Flowsolver on openSUSE Leap 15.1 (Linux). They have been tested on a completely fresh install of the OS. You will need to adjust them accordingly for other flavours of Linux.
    
1. Install the dependencies using the system pacakge manager
 - `sudo zypper install kernel-default-devel git scons openmpi openmpi-devel gcc-c++ gcc-fortran python-devel cblas-devel blas-devel lapack-devel cmake libboost_headers1_66_0-devel libboost_thread1_66_0-devel libboost_filesystem1_66_0-devel python2-numpy python2-scipy`
2. Setup paths
 - Add the mpi compilers `bin` folder to your system path (in this case, `export PATH=$PATH:/usr/lib64/mpi/gcc/openmpi/bin` - but this will almost certainly differ on your system if it is not opensuse Leap 15.1) and their runtime libraries to the `LD_LIBRARY_PATH` (in the opensuse case: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/mpi/gcc/openmpi/lib64` - but check for your system where the mpi runtime libraries are located)
 - You probably want to put these exports in your `~/.bashrc` file and then open a fresh terminal in order to activate it.
3. Install PetSc
 - Get petsc-3.2-p7 from the [PetSc website](https://www.mcs.anl.gov/petsc/download/index.html). Untar it in some folder and compile it, but note that the values for the lib locations may vary on your system, and use a `PETSC_ARCH` according to your gcc version (although this is just for your own records): `./configure PETSC_ARCH=gcc_7.4.1-opt --with-memalign=64 --with-debugging=0 COPTFLAGS=-O2 -fp-model precise FOPTFLAGS=-O2 -fp-model source --with-x=0 CXXOPTFLAGS=-O2 --with-blas-lib=/usr/lib64/libblas.so --with-lapack-lib=/usr/lib64/liblapack.so --with-mpi-dir=/usr/lib64/mpi/gcc/openmpi/`
 - After configuring is complete, follow the instructions in the console for compiling PetSc
 - After compilation, follow the instructions in the console for testing PetSc. Note that the fortran compiler may issue some warnings, which the petsc test system will report as failures - but if the only actual error output are the warnings, it is most likely fine
4. Install VTK version 5.8.0
 - Get vtk version 5.8.0 - it is available with CRIMSON modifications [here](https://github.com/carthurs/VTK-5.8.0)
 - Run the following configuration command, but make sure you set `/some/install/path` to be where you want vtk to be installed - and remember it for later setting up of CRIMSON `cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DVTK_USE_RENDERING=OFF -DCMAKE_INSTALL_PREFIX=/some/install/path -DCMAKE_CXX_FLAGS=-Wno-error=narrowing`
 - make && make install
5. Configure CRIMSON Flowsolver
 - Copy `buildoptions.py.example` from the CRIMSON root to `buildoptions.py`.
 - Edit `buildoptions.py`, setting the `#setme` fields. Not all of these will be required, particularly if some of the libraries are in the system lib64 folder. As a minimum, you'll need to set PETSC_TOP (to point at the petsc root folder) PETSC_BUILD (to be the value you set in PETSC_ARCH during the petsc configure call), VTK_TOP (to be whatever CMAKE_INSTALL_PREFIX you set during the vtk cmake step). Be careful to note how the paths are assembled in this script: relative paths and special characters such as ~ are not supported.
6. Compile CRIMSON Flowsolver
 - In the CRIMSON root, run `scons -jN test=1 debug=0` to compile using N processes
 - After a successful build, go into the `testbin` subfolder, and run `./mytest N`, to test on N processes. Allow the test suite to finish.
7. Use the CRIMSON Flowsolver on your own simulations
 - To run CRIMSON Flowsolver, use the script in the bin folder: `mysolver`. Invoke this from the directory containing your simulation files, as `/path/to/mytest N solver.inp` to run on N processes