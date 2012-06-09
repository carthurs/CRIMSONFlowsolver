# Hey emacs, this is a -*- makefile -*-

# Copyright (c) 2000-2007, Stanford University, 
#    Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
#    Charles A. Taylor (see SimVascular Acknowledgements file 
#    for additional contributors to the source code).
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions 
# are met:
#
# Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer. 
# Redistributions in binary form must reproduce the above copyright 
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution. 
# Neither the name of the Stanford University or Rensselaer Polytechnic
# Institute nor the names of its contributors may be used to endorse or
# promote products derived from this software without specific prior 
# written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.

TOP=../../

include $(TOP)/buildoptions

# --------------
# Compiler flags
# --------------

CXXFLAGS = $(GLOBAL_CXXFLAGS) $(SOLVERIO_INCDIR) $(NSPCG_INCDIR) $(ZLIB_INCDIR)
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
             $(TOP)/lib/lib_simvascular_solverio.$(STATICEXT)
	mv $(TARGET_PRESOLVER)* $(TOP)/bin

clean:
	-$(RM) $(OBJS) $(TOP)/bin/$(TARGET_PRESOLVER) *.pdb $(TOP)/bin/$(TARGET_PRESOLVER).manifest