TOP=../..

include $(TOP)/buildoptions.mk

# --------------
# Compiler flags
# --------------

CXXFLAGS = $(GLOBAL_CXXFLAGS) $(SOLVERIO_INCDIR)

CXXSRCS = Reduce.cxx wrtc.cxx

SRCS	= $(CXXSRCS)

OBJS	= $(CXXSRCS:.cxx=.$(OBJECTEXT))

$(TOP)/Bin/$(TARGET_POSTSOLVER): $(OBJS)
	$(LINK_EXE)$(TARGET_POSTSOLVER) $(GLOBAL_LFLAGS) $(OBJS) $(SOLVERIO_LIBS)
	mv $(TARGET_POSTSOLVER)* $(TOP)/bin

clean:
	-$(RM) $(OBJS) $(TOP)/bin/$(TARGET_POSTSOLVER)*
