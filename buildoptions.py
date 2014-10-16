import sys

SHELL           = '/bin/sh'
CXX             = 'mpicxx'
CC              = 'mpicc'
F90             = 'mpif90'
F77             = 'mpif90'
AR              = 'xiar'# 'xiar cr' no longer necessary, as scons deals with this automatically
 
GLOBAL_DEFINES = ['-DUNIX']
 
#STATICEXT = a
#OBJECTEXT = o

MAKE_WITH_PRESOLVER = 1
MAKE_WITH_POSTSOLVER = 1
MAKE_WITH_FLOWSOLVER = 1
MAKE_WITH_ESTIMATOR = 1

#TARGET_FLOWSOLVER = flowsolver
#TARGET_ESTIMATOR = estimator
#TARGET_FORWARDANDOBSERVE = flowsolver_and_observer
#TARGET_POSTSOLVER = postsolver
#TARGET_PRESOLVER = presolver

# This should be set to zero on the sgi
EXTRA_CONSOLE_OUTPUT_ON = 0

# MAKE_OPTIMIZED - this is now a command line option: e.g. "scons -j3 debug" does a debug build, otherwise, it builds optimised.

# if MAKE_OPTIMIZED==1:
   # DEBUG_FLAGS     = ['']
   # DEBUG_FFLAGS    = ['']
OPT_FLAGS       = ['-O3', '-std=c++0x']
OPT_FFLAGS      = ['-O3','-align','array64byte','-fpp','-D\'EXTRA_CONSOLE_OUTPUT='+str(EXTRA_CONSOLE_OUTPUT_ON)+'\'']
   #LINK_EXE        = '$(F90) -nofor-main -cxxlib -o'
   ##LINK_EXE        = $(CXX) -o
# else:
DEBUG_FLAGS     = ['-O0','-g','-std=c++0x']
DEBUG_FFLAGS    = ['-O0','-g','-align','array64byte','-traceback','-fp-stack-check','-fpp','-DEXTRA_CONSOLE_OUTPUT='+str(EXTRA_CONSOLE_OUTPUT_ON)]
   # OPT_FLAGS       = ['']
   # OPT_FFLAGS      = ['']
   ##LINK_EXE        = $(CXX) -o
   #LINK_EXE        = '$(F90) -nofor-main -cxxlib -o'

 
BUILDFLAGS      = GLOBAL_DEFINES
# GLOBAL_CXXFLAGS = BUILDFLAGS+DEBUG_FLAGS+OPT_FLAGS
# GLOBAL_CFLAGS   = BUILDFLAGS+DEBUG_FLAGS+OPT_FLAGS
# GLOBAL_FFLAGS   = BUILDFLAGS+DEBUG_FFLAGS+OPT_FFLAGS
CXX_LIBS    = ''
CXX_LIBS    = '-lrt'
F90_LIBS    = ['mpi_f90']
 
HOME_DIR = '/home/carthurs/workspace/merged_simvascular/simvascular_flowsolver_estimator'
EXTERNAL_LIB_DIR = HOME_DIR+'/external/x64-linux'
 
SOLVERIO_INCDIR   = [HOME_DIR+'/solverio/src']
SOLVERIO_LIBSDIR  = [HOME_DIR+'/lib']
SOLVERIO_LIBSLIST = ['simvascular_solverio']

FLOWSOLVER_INCDIR = [HOME_DIR+'/flowsolver/src']
FLOWSOLVER_LIBSDIR = [HOME_DIR+'/lib']
FLOWSOLVER_LIBSLIST = ['simvascular_flowsolver']
 
# MPI_TOP        = '/opt/sgi/mpt/mpt-2.05'
# MPI_INCDIR     = ' -I '+MPI_TOP+'/include'
# MPI_LIBSDIR    = MPI_TOP+'lib'
# MPI_LIBSLIST   = 'mpi'

INTEL_TOP      = '/opt/intel/composerxe'
INTEL_INCDIR   = [INTEL_TOP+'/include']
INTEL_LIBSDIR  = [INTEL_TOP+'/lib/intel64']
INTEL_LIBSLIST = ['ifcore','ifport','imf','irc']


MPI_TOP        = '/usr/local/OpenMPI-intel'
MPI_INCDIR     = [MPI_TOP+'/include']
MPI_LIBSDIR    = [MPI_TOP+'/lib']
MPI_LIBSLIST   = ['mpi']

METIS_TOP      = EXTERNAL_LIB_DIR+'/metis-4.0'
METIS_INCDIR   = [METIS_TOP]
METIS_LIBSDIR  = [METIS_TOP]
METIS_LIBSLIST = ['metis']
 
LESLIB_DEFS    = ['-DACUSIM_LINUX','-DACUSIM_LESLIB_VER_1_4']
LESLIB_TOP     = EXTERNAL_LIB_DIR+'/leslib-1.4'
LESLIB_INCDIR  = [LESLIB_TOP]
LESLIB_LIBSDIR = [LESLIB_TOP,LESLIB_TOP+'/deps']
LESLIB_LIBSLIST= ['les','ifcore','ifport','imf','irc'] #-L $(LESLIB_TOP)/deps -lifcore -lifport -limf -lirc

#MEMLS_TOP        = $(HOME_DIR)/memLS
#MEMLS_INCDIR     = -I $(MEMLS_TOP)/include
#MEMLS_DEF        =
#MEMLS_LIBS       = -L $(MEMLS_TOP)/lib -lmemLS

NSPCG_TOP      = EXTERNAL_LIB_DIR+'/NSPCG'
NSPCG_INCDIR   = [NSPCG_TOP]
NSPCG_LIBSDIR  = [NSPCG_TOP]
NSPCG_LIBSLIST = ['nspcg']

SPARSE_TOP     = EXTERNAL_LIB_DIR+'/sparse-1.4'
SPARSE_INCDIR  = [SPARSE_TOP]
SPARSE_LIBSDIR = [SPARSE_TOP]
SPARSE_LIBSLIST= ['sparse']

ZLIB_TOP       = EXTERNAL_LIB_DIR+'/zlib-1.2.3'
ZLIB_INCDIR    = [ZLIB_TOP]
ZLIB_LIBSDIR   = [ZLIB_TOP]
ZLIB_LIBSLIST  = ['z']

VERDANDI_TOP  = HOME_DIR+'/../verdandi-1.2.1'
VERDANDI_INCDIR  = [VERDANDI_TOP, \
		 VERDANDI_TOP+'/container', \
		 VERDANDI_TOP+'/error', \
		 VERDANDI_TOP+'/method', \
		 VERDANDI_TOP+'/model', \
		 VERDANDI_TOP+'/observation_manager', \
		 VERDANDI_TOP+'/output_saver', \
		 VERDANDI_TOP+'/share', \
		 VERDANDI_TOP+'/include', \
		 VERDANDI_TOP+'/include/lua/src', \
		 VERDANDI_TOP+'/include/seldon']

SELDON_TOP  = VERDANDI_TOP
SELDON_INCDIR  = [SELDON_TOP,'/include/seldon']

LUA_LIBSDIR = [VERDANDI_TOP+'/include/lua/src/']
LUA_LIBSLIST = ['lua']
# LUA_LIBSLIST = 'liblua.a'

PETSC_TOP    = '/home/carthurs/workspace/merged_simvascular/petsc-3.2-p7'
#PETSC_BUILD = arch-linux-intel-refblaslapack-opt
PETSC_BUILD = 'arch-linux-intel-intelMKL-opt'
PETSC_INCDIR = [PETSC_TOP+'/include',PETSC_TOP+'/'+PETSC_BUILD+'/include']
PETSC_LIBSDIR   = [PETSC_TOP+'/'+PETSC_BUILD+'/lib']
PETSC_LIBSLIST = ['petsc']

#CGAL_TOP     = $(HOME_DIR)/../CGAL-install
#CGAL_INCDIR  = -I $(CGAL_TOP)/include
#CGAL_LIBS    = -L $(CGAL_TOP)/lib -lCGAL -lCGAL_Core -lCGAL_ImageIO

BOOSTCPP_TOP    = HOME_DIR+'/../boost_1_40_0'
BOOSTCPP_INCDIR = [BOOSTCPP_TOP]
BOOSTCPP_LIBS   = '-L '+BOOSTCPP_TOP+'/stage/lib -lboost_thread'

VTK_TOP    = '/home/carthurs/workspace/merged_simvascular/vtk-5.8.0-install'
VTK_INCDIR = [VTK_TOP+'/include/vtk-5.8']
VTK_LIBSDIR   = [VTK_TOP+'/lib/vtk-5.8/']
VTK_LIBSLIST = ['vtkGraphics','vtkFiltering','vtkGenericFiltering','vtkIO','vtkCommon','vtksys']

BLASLAPACK_TOP    = '/opt/intel/mkl'
BLASLAPACK_INCDIR = [BLASLAPACK_TOP+'/include']
BLASLAPACK_LIBSDIR   = [BLASLAPACK_TOP+'/lib/intel64']
BLASLAPACK_LIBSLIST =  ['mkl_intel_lp64', 'mkl_sequential', 'mkl_core', 'pthread']

#BLASLAPACK_TOP    = /home/nxiao/dev/refblaslapack
#BLASLAPACK_INCDIR = 
#BLASLAPACK_LIBS   = -L $(BLASLAPACK_TOP) -llapack -lcblas -lblas


# %.$(OBJECTEXT): %.cxx
# 	$(CXX) $(CXXFLAGS) -c $<

# %.$(OBJECTEXT): %.c
# 	$(CC) $(CFLAGS) -c $<

# %.$(OBJECTEXT): %.f90
# 	$(F90) $(FFLAGS) -c $<

# %.$(OBJECTEXT): %.f
# 	$(F77) $(FFLAGS) -c $<