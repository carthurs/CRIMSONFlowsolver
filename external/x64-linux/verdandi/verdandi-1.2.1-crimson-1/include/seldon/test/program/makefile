# C++ compiler and options
CC = mpicxx
FLAGS = -g -DSELDON_WITH_MPI -DSELDON_WITH_UMFPACK

# file to compile and executable
#MAIN_CC = direct_test.cpp
MAIN_CC = matrix_test.cpp
BIN = main.x

# directory where Seldon is 
INCLUDE := -I../../

# basic libraries : BLAS, LAPACK
LIB := -L$(HOME)/Lapack -llapack -L$(HOME)/Cblas -lcblas -L$(HOME)/Blas -lblas -lgfortran

# Mumps
SOLVE = $(HOME)/Solve
DIR_MUMPS = $(SOLVE)/MUMPS_4.8.4
INCLUDE := -I$(DIR_MUMPS)/include \
           -I$(DIR_METIS)/metis-4.0/Lib $(INCLUDE)
DIR_SCOTCH = $(SOLVE)/scotch_5.0_esmumps
LIB := -L$(DIR_MUMPS)/lib -ldmumps -lzmumps -lmumps_common -lpord \
  -L$(DIR_SCOTCH)/lib -lesmumps -lfax -lsymbol -ldof -lorder -lgraph -lscotch -lscotcherr -lcommon \
  -L$(HOME)/common/lib -lm -lpthread -L/cvos/shared/apps/scalapack/gcc/1.8.0/lib64 -lscalapack \
  -L$(HOME)/Blacs -lblacs \
    -L$(SOLVE)/metis-4.0 -lmetis $(LIB)

# UmfPack
DIR_UMF = $(SOLVE)/UMFPACKv5.02
INCLUDE := -I$(DIR_UMF)/AMD/Include -I$(DIR_UMF)/UFconfig \
           -I$(DIR_UMF)/UMFPACK/Include $(INCLUDE)
LIB := -L$(DIR_UMF)/UMFPACK/Lib -lumfpack -L$(DIR_UMF)/AMD/Lib -lamd $(LIB)

# SuperLU
DIR_SUPERLU = $(SOLVE)/SuperLU
INCLUDE := -I$(DIR_SUPERLU)/SRC $(INCLUDE)
LIB := -L$(DIR_SUPERLU)/lib -lsuperlu $(LIB)

# compilation directive
all :
	$(CC) $(FLAGS) $(MAIN_CC) -o test.x $(INCLUDE) $(LIB)
	mv test.x $(BIN)
