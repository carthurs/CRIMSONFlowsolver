SHELL           = /bin/sh
CXX             = mpicxx
CC              = mpicc
F90             = mpif90
AR              = ar -cru 
 
GLOBAL_DEFINES = -DUSE_NOTIMER -DUNIX -Dsgi
 
STATICEXT = a
OBJECTEXT = o

MAKE_WITH_PRESOLVER = 1
MAKE_WITH_POSTSOLVER = 1
MAKE_WITH_FLOWSOLVER = 1

TARGET_FLOWSOLVER = flowsolver
TARGET_POSTSOLVER = postsolver
TARGET_PRESOLVER = presolver
 
MAKE_OPTIMIZED = 1
 
ifeq ($(MAKE_OPTIMIZED),1)
   DEBUG_FLAGS     =
   DEBUG_FFLAGS    =
   OPT_FLAGS       = -O2
   OPT_FFLAGS      = -O2
   LINK_EXE        = $(F90) -nofor-main -cxxlib -o 
else
   DEBUG_FLAGS     = -O0 -debug -g2   
   DEBUG_FFLAGS    = -g -zero -fpstkchk
   OPT_FLAGS       =
   OPT_FFLAGS      =
   LINK_EXE        = $(F90) -nofor-main -cxxlib -o 
endif
 
BUILDFLAGS      = $(GLOBAL_DEFINES)
GLOBAL_CXXFLAGS = $(BUILDFLAGS) $(DEBUG_FLAGS) $(OPT_FLAGS) 
GLOBAL_CFLAGS   = $(BUILDFLAGS) $(DEBUG_FLAGS) $(OPT_FLAGS)
GLOBAL_FFLAGS   = $(BUILDFLAGS) $(DEBUG_FFLAGS) $(OPT_FFLAGS)
 
CXX_LIBS    = 
F90_LIBS    = 
 
HOME_DIR = ../..
EXTERNAL_LIB_DIR = $(HOME_DIR)/external/x64-linux
 
SOLVERIO_INCDIR = -I $(HOME_DIR)/solverio/src
 
MPI_TOP        = 
MPI_INCDIR     = 
MPI_LIBS       = 
 
METIS_TOP      = $(EXTERNAL_LIB_DIR)/metis-4.0
METIS_INCDIR   = -I $(METIS_TOP)
METIS_LIBS     = -L $(METIS_TOP) -lmetis
 
LESLIB_DEFS    = -DACUSIM_LINUX -DACUSIM_LESLIB_VER_1_4
LESLIB_TOP     = $(EXTERNAL_LIB_DIR)/leslib-1.4
LESLIB_INCDIR  = -I $(LESLIB_TOP)
LESLIB_LIBS    = -L $(LESLIB_TOP) -lles

NSPCG_TOP      = $(EXTERNAL_LIB_DIR)/NSPCG
NSPCG_INCDIR   = -I $(NSPCG_TOP)
NSPCG_LIBS     = -L $(NSPCG_TOP) -lnspcg

SPARSE_TOP     = $(EXTERNAL_LIB_DIR)/sparse-1.4
SPARSE_INCDIR  = -I $(SPARSE_TOP)
SPARSE_LIBS    = -L $(SPARSE_TOP) -lsparse

ZLIB_TOP       = $(EXTERNAL_LIB_DIR)/zlib-1.2.3
ZLIB_INCDIR    = -I $(ZLIB_TOP)
ZLIB_LIBS      = -L $(ZLIB_TOP) -lz

%.$(OBJECTEXT): %.cxx
	$(CXX) $(CXXFLAGS) -c $<

%.$(OBJECTEXT): %.c
	$(CC) $(CFLAGS) -c $<

%.$(OBJECTEXT): %.f
	$(F90) $(FFLAGS) -c $<
