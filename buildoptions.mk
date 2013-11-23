SHELL           = /bin/sh
CXX             = mpicxx
CC              = mpicc
F90             = mpif90
F77             = mpif77
AR              = xiar cr 
 
GLOBAL_DEFINES = -DUNIX
 
STATICEXT = a
OBJECTEXT = o

MAKE_WITH_PRESOLVER = 1
MAKE_WITH_POSTSOLVER = 1
MAKE_WITH_FLOWSOLVER = 1
MAKE_WITH_MAIN = 1
MAKE_WITH_ESTIMATOR = 1

TARGET_FLOWSOLVER = flowsolver
TARGET_ESTIMATOR = estimator
TARGET_FORWARDANDOBSERVE = flowsolver_and_observer
TARGET_POSTSOLVER = postsolver
TARGET_PRESOLVER = presolver
 
MAKE_OPTIMIZED = 1
 
ifeq ($(MAKE_OPTIMIZED),1)
   DEBUG_FLAGS     =
   DEBUG_FFLAGS    =
   OPT_FLAGS       = -O3 #-fp-model precise
   OPT_FFLAGS      = -O3 -align array64byte #-fp-model source
   LINK_EXE        = $(F90) -nofor-main -cxxlib -o 
   #LINK_EXE        = $(CXX) -o 
else
   DEBUG_FLAGS     = -O0 -g -fp-model precise
   DEBUG_FFLAGS    = -O0 -g -align array64byte -fp-model source -traceback -fp-stack-check
   OPT_FLAGS       =
   OPT_FFLAGS      =
   #LINK_EXE        = $(CXX) -o 
   LINK_EXE        = $(F90) -nofor-main -cxxlib -o 
endif
 
BUILDFLAGS      = $(GLOBAL_DEFINES)
GLOBAL_CXXFLAGS = $(BUILDFLAGS) $(DEBUG_FLAGS) $(OPT_FLAGS)
GLOBAL_CFLAGS   = $(BUILDFLAGS) $(DEBUG_FLAGS) $(OPT_FLAGS)
GLOBAL_FFLAGS   = $(BUILDFLAGS) $(DEBUG_FFLAGS) $(OPT_FFLAGS)
#CXX_LIBS    =
CXX_LIBS    = -lrt
#F90_LIBS    =
F90_LIBS    =
 
HOME_DIR = ../..
EXTERNAL_LIB_DIR = $(HOME_DIR)/external/x64-linux
 
SOLVERIO_INCDIR = -I $(HOME_DIR)/solverio/src
SOLVERIO_LIBS = -L $(HOME_DIR)/lib -lsimvascular_solverio

FLOWSOLVER_INCDIR = -I $(HOME_DIR)/flowsolver/src
FLOWSOLVER_LIBS = -L $(HOME_DIR)/lib -lsimvascular_flowsolver
 
MPI_TOP        = 
MPI_INCDIR     = #-I
MPI_LIBS       = #-lmpi_f90 -lmpi_f77 -lmpi_cxx -lmpi
 
METIS_TOP      = $(EXTERNAL_LIB_DIR)/metis-4.0
METIS_INCDIR   = -I $(METIS_TOP)
METIS_LIBS     = -L $(METIS_TOP) -lmetis
 
LESLIB_DEFS    = -DACUSIM_LINUX -DACUSIM_LESLIB_VER_1_4
LESLIB_TOP     = $(EXTERNAL_LIB_DIR)/leslib-1.4
LESLIB_INCDIR  = -I $(LESLIB_TOP)
LESLIB_LIBS    = -L $(LESLIB_TOP) -lles #-L $(LESLIB_TOP)/deps -lifcore -lifport -limf -lirc

NSPCG_TOP      = $(EXTERNAL_LIB_DIR)/NSPCG
NSPCG_INCDIR   = -I $(NSPCG_TOP)
NSPCG_LIBS     = -L $(NSPCG_TOP) -lnspcg

SPARSE_TOP     = $(EXTERNAL_LIB_DIR)/sparse-1.4
SPARSE_INCDIR  = -I $(SPARSE_TOP)
SPARSE_LIBS    = -L $(SPARSE_TOP) -lsparse

ZLIB_TOP       = $(EXTERNAL_LIB_DIR)/zlib-1.2.3
ZLIB_INCDIR    = -I $(ZLIB_TOP)
ZLIB_LIBS      = -L $(ZLIB_TOP) -lz

VERDANDI_TOP  = $(HOME_DIR)/../verdandi-1.2.1
VERDANDI_INCDIR  = -I $(VERDANDI_TOP) \
		 -I $(VERDANDI_TOP)/container \
		 -I $(VERDANDI_TOP)/error \
		 -I $(VERDANDI_TOP)/method \
		 -I $(VERDANDI_TOP)/model \
		 -I $(VERDANDI_TOP)/observation_manager \
		 -I $(VERDANDI_TOP)/output_saver \
		 -I $(VERDANDI_TOP)/share \
		 -I $(VERDANDI_TOP)/include \
		 -I $(VERDANDI_TOP)/include/lua/src \
		 -I $(VERDANDI_TOP)/include/seldon \

SELDON_TOP  = $(VERDANDI_TOP)
SELDON_INCDIR  = -I $(SELDON_TOP)/include/seldon

LUA_LIBS = $(VERDANDI_TOP)/include/lua/src/liblua.a

PETSC_TOP    = $(HOME_DIR)/../petsc-3.2-p7
#PETSC_BUILD = arch-linux-intel-refblaslapack-opt
PETSC_BUILD = arch-linux-intel-intelMKL-opt
PETSC_INCDIR = -I $(PETSC_TOP)/include -I $(PETSC_TOP)/$(PETSC_BUILD)/include
PETSC_LIBS   = -L $(PETSC_TOP)/$(PETSC_BUILD)/lib -lpetsc

CGAL_TOP     = $(HOME_DIR)/../CGAL-install
CGAL_INCDIR  = -I $(CGAL_TOP)/include
CGAL_LIBS    = -L $(CGAL_TOP)/lib -lCGAL -lCGAL_Core -lCGAL_ImageIO

BOOSTCPP_TOP    = $(HOME_DIR)/../boost_1_40_0
BOOSTCPP_INCDIR = -I $(BOOSTCPP_TOP)
BOOSTCPP_LIBS   = -L $(BOOSTCPP_TOP)/stage/lib -lboost_thread

VTK_TOP    = $(HOME_DIR)/../VTK-install
VTK_INCDIR = -I $(VTK_TOP)/include/vtk-5.8
VTK_LIBS   = -L $(VTK_TOP)/lib/vtk-5.8 -lvtkGraphics -lvtkFiltering -lvtkGenericFiltering -lvtkIO -lvtkCommon -lvtksys 

BLASLAPACK_TOP    = /opt/intel/mkl
BLASLAPACK_INCDIR = -I $(BLASLAPACK_TOP)/include
BLASLAPACK_LIBS   = -L $(BLASLAPACK_TOP)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread

#BLASLAPACK_TOP    = /home/nxiao/dev/refblaslapack
#BLASLAPACK_INCDIR = 
#BLASLAPACK_LIBS   = -L $(BLASLAPACK_TOP) -llapack -lcblas -lblas


%.$(OBJECTEXT): %.cxx
	$(CXX) $(CXXFLAGS) -c $<

%.$(OBJECTEXT): %.c
	$(CC) $(CFLAGS) -c $<

%.$(OBJECTEXT): %.f90
	$(F90) $(FFLAGS) -c $<

%.$(OBJECTEXT): %.f
	$(F77) $(FFLAGS) -c $<
