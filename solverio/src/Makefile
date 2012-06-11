TOP=../../

include $(TOP)/buildoptions.mk

# --------------
# Compiler flags
# --------------

CXXFLAGS = $(GLOBAL_CXXFLAGS)

HDRS	= cvSolverIO.h
CXXSRCS = cvSolverIO.cxx

SRCS	= $(CSRCS) $(CXXSRCS)

OBJS	= $(CSRCS:.c=.$(OBJECTEXT)) $(CXXSRCS:.cxx=.$(OBJECTEXT))

TARGET_LIB = libsimvascular_solverio.$(STATICEXT)

$(TOP)/lib/$(TARGET_LIB) : $(OBJS)
	$(AR)$(TARGET_LIB) $(OBJS)
	mv $(TARGET_LIB) $(TOP)/lib
	
clean:
	-$(RM) $(OBJS) $(TOP)/lib/$(TARGET_LIB) *.pdb *.manifest
