TOP = ../..

include $(TOP)/buildoptions.mk

CFLAGS = $(GLOBAL_CFLAGS) $(MPI_INCDIR) $(LESLIB_DEFS) $(LESLIB_INCDIR) $(SOLVERIO_INCDIR)
CXXFLAGS  = $(GLOBAL_CXXFLAGS) $(MPI_INCDIR) $(METIS_INCDIR) $(LESLIB_DEFS) $(LESLIB_INCDIR) $(SOLVERIO_INCDIR) $(VTK_INCDIR)
FFLAGS   = $(GLOBAL_FFLAGS) $(MPI_INCDIR)

# required link flag because of leslib 1.4
#CXXFLAGS += /MT
#CFLAGS += /MT
#FFLAGS += /MT
#CXX_LIBS = libcmt.lib libcpmt.lib

HDRS  =   CInput.h \
          common_c.h globalArrays.h input.h input_fform.h itrdrv.h \
          partition.h phSolver.h proces.h \
          shapeFuncInternals.h solverFields.h topo_shapedefs.h\
          usr.h ValType.h

CXXSRCS = newshape.cxx Cinput.cxx input_fform.cxx \
          partition.cxx \
          ValType.cxx \
          multidom.cxx fileReaders.cxx fortranPointerManager.cxx  \

CSRCS   = blend.c blendDrv.c Bn.c BnDrv.c En.c EnDrv.c \
          Fn.c FnDrv.c getIntPnts.c lestools.c modeShape.c \
          modeShapeDrv.c new_interface.c parDrv.c shp6w.c shphex.c shppyr.c \
          shptet.c shptri.c symhex.c symline.c sympyr.c symquad.c symquadw.c \
          symtet.c symtri.c symtripyr.c symwdg.c tmrc.c uniformP.c usr.c write_hessian.c

FSRCS   = cpp_interface.f90 common.f90 pointer.f90 globalArrays.f90 shapeTable.f90 memLS.f90 \
		turbsa.f90 \
		perprep.f90 local_mass.f90 \
		pvsqbi.f90 \
		boundarydata.f90 \
		dtn.f90 \
		timedata.f90 \
		bctint.f90 readnblk.f90 \
		distmeas.f90 \
          vlmwsttri.f90 \
          stats.f90 \
		  general_rcrbc.f90 \
		  debuggingTools.f90 \
		  datatypes.f90 \
		  externalDataTools.f90 \
		  multidomain.f90 \
          itrdrv.f90 \
          error.f90 \
          asadj.f90 asbflx.f90 asbmfg.f90 asbnabi.f90 asigmr.f90 \
          asiq.f90 bc3diag.f90 bc3lhs.f90 bc3per.f90 \
          bc3res.f90 bctcardio.f90 bflux.f90 cname.f90 commu.f90 ctypes.f90 \
          e3.f90 e3b.f90 e3bvar.f90 e3dc.f90 e3ivar.f90 e3lhs.f90 e3metric.f90 e3q.f90 e3ql.f90 \
          e3qvar.f90 e3res.f90 e3stab.f90 e3sts.f90 \
          elm3keps.f90 elmdist.f90 elmgmr.f90 elmStats.f90 errsmooth.f90 f3lhs.f90 \
          fillsparse.f90 forces.f90 ftools.f90 genadj.f90 genbc.f90 genbc1.f90 genbkb.f90 \
          genblk.f90 gendat.f90 genibc.f90 genini.f90 genint.f90 genlmass.f90 gensav.f90 gensvb.f90 genshp.f90 \
          genshpb.f90 gensvbDef.f90 get_a_not_hex.f90 getdiff.f90 \
          gtnods.f90 hessian.f90 hierarchic.f90 \
          input.f90 itrbc.f90 itrPC.f90 lesSparse.f90 local.f90 localy.f90 lubksb.f90 ludcmp.f90 \
          mpitools.f90 mpoint.f90 proces.f90 qpbc.f90 ramg_driver.f90 restar.f90 \
          rotabc.f90 rstatic.f90 solfar.f90 \
          solvecon.f90 timeseries.f90 \
          estimation_helpers.f90

SRCS  = $(FSRCS) $(CSRCS) $(CXXSRCS)

OBJS  = $(FSRCS:.f90=.$(OBJECTEXT)) $(CSRCS:.c=.$(OBJECTEXT)) $(CXXSRCS:.cxx=.$(OBJECTEXT))

#$(TARGET_FLOWSOLVER): $(OBJS)
#	$(LINK_EXE)$(TARGET_FLOWSOLVER) $(GLOBAL_LFLAGS) \
#						 $(OBJS) $(METIS_LIBS) \
#						 $(LESLIB_LIBS) \
#						 $(MPI_LIBS) \
#						 $(F90_LIBS) \
#						 $(CXX_LIBS) \
#						 $(TOP)/lib/lib_simvascular_solverio.$(STATICEXT)
#	mv $(TARGET_FLOWSOLVER) $(TOP)/bin
	
TARGET_FLOWSOLVERLIB = libsimvascular_flowsolver.$(STATICEXT)
	
$(TOP)/lib/$(TARGET_FLOWSOLVERLIB) : $(OBJS)
	$(AR)$(TARGET_FLOWSOLVERLIB) $(OBJS)
	mv $(TARGET_FLOWSOLVERLIB) $(TOP)/lib	

.PHONY : clean
clean : 
	-$(RM) $(OBJS) $(TARGET_FLOWSOLVERLIB) *.mod *.exe *.pdb *.manifest
	-$(RM) $(TOP)/lib/$(TARGET_FLOWSOLVERLIB) 
