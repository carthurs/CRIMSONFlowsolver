TOP=../..

include $(TOP)/buildoptions.mk

# --------------
# Compiler flags
# --------------

CXXFLAGS = $(GLOBAL_CXXFLAGS) $(SOLVERIO_INCDIR) $(NSPCG_INCDIR) $(ZLIB_INCDIR) $(SELDON_INCDIR)
CFLAGS  = $(GLOBAL_CXXFLAGS) $(SPARSE_INCDIR) $(ZLIB_INCDIR)
FFLAGS   = $(GLOBAL_FFLAGS) $(NSPCG_INCDIR)

HDRS	= cmd.h

CSRCS   = directsolve.c

CXXSRCS = cmd.cxx displacements.cxx helpers.cxx supre-cmds.cxx supre.cxx

FSRCS   = stannspcg.f

SRCS	= $(CSRCS) $(CXXSRCS)

DEPS	= $(CSRCS:.c=.d) $(CXXSRCS:.cxx=.d)

OBJS	= $(CSRCS:.c=.$(OBJECTEXT)) $(CXXSRCS:.cxx=.$(OBJECTEXT)) $(FSRCS:.f=.$(OBJECTEXT))

$(TOP)/Bin/$(TARGET_PRESOLVER): $(OBJS)
	$(LINK_EXE)$(TARGET_PRESOLVER) $(GLOBAL_LFLAGS) \
	$(OBJS) $(F90_LIBS) $(CXX_LIBS) \
	$(SPARSE_LIBS) $(NSPCG_LIBS) $(ZLIB_LIBS) \
	$(TOP)/lib/libsimvascular_solverio.$(STATICEXT)
	mv $(TARGET_PRESOLVER)* $(TOP)/bin

clean:
	-$(RM) $(OBJS) $(TOP)/bin/$(TARGET_PRESOLVER)*
