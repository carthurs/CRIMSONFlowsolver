TOP = ../..

include $(TOP)/buildoptions.mk

CFLAGS = $(GLOBAL_CFLAGS) $(MPI_INCDIR) $(LESLIB_DEFS) $(LESLIB_INCDIR) $(SOLVERIO_INCDIR)
CXXFLAGS  = $(GLOBAL_CXXFLAGS) $(MPI_INCDIR) $(METIS_INCDIR) $(LESLIB_DEFS) $(LESLIB_INCDIR) $(SOLVERIO_INCDIR)
FFLAGS   = $(GLOBAL_FFLAGS) $(MPI_INCDIR)

# required link flag because of leslib 1.4
#CXXFLAGS += /MT
#CFLAGS += /MT
#FFLAGS += /MT
#CXX_LIBS = libcmt.lib libcpmt.lib

HDRS  =   Input.h ValType.h auxmpi.h common.h \
          common_c.h topo_shapedefs.h usr.h \
          vec_func.h \
          shapeFuncInternals.h \

CXXSRCS = Cinput.cxx ValType.cxx input_fform.cxx \
          newshape.cxx partition.cxx phasta.cxx \
          TriShapeAndDrv.cxx \

CSRCS   = GaussLegendre1D.c GaussLegendreSimplex.c \
          flush.c getIntPnts.c getSol.c gltri.c \
          lestools.c main.c missing_win32.c new_interface.c \
          shp6w.c shphex.c shppyr.c shptet.c symhex.c \
          sympyr.c symquad.c symquadw.c symtet.c symtri.c \
          symtripyr.c symwdg.c tmrc.c usr.c write_hessian.c \
          Bn.c BnDrv.c En.c EnDrv.c Fn.c \
          FnDrv.c Lagrange.c blend.c \
          blendDrv.c modeShape.c \
          modeShapeDrv.c parDrv.c \
          uniformP.c \
          rand.c

FSRCS   = spebc.f pointer.f \
          turbsa.f \
          perprep.f local_mass.f \
          filtprep.f  pvsqbi.f \
          dtn.f bardmc.f \
          timedata.f \
          turbke.f \
          bctint.f readnblk.f \
          distmeas.f \
          vlmwsttri.f \
          asadj.f asbnabi.f asbwmod.f asithf.f aveprep.f \
          blkdta.f clear.f cmass.f cname.f commu.f ctypes.f \
          e3metric.f elem-search.f elm3keps.f eqn_plane.f error.f \
          f3lhs.f fillsparse.f genadj.f genbc.f genbc1.f \
          genbkb.f genblk.f gendat.f genibc.f genini.f genint.f \
          gensav.f genscale.f genshp.f genshpb.f gensvb.f gensvbDef.f \
          get_a_not_hex.f \
          get_coeff.f getdmc.f getstrl.f gettab.f getvel.f gtnods.f \
          hierarchic.f input.f local.f localy.f lubksb.f \
          ludcmp.f mpiset.f mpitools.f mpoint.f \
          proces.f qpbc.f renum.f rerun_check.f \
          restar.f rotabc.f rwvelb.f scatnu.f settauw.f shp10t.f \
          shp4t.f solvecon.f sonfath.f timeseries.f \
          stats.f \
          advLES.f asbflx.f asbmfg.f asigmr.f asiq.f bc3diag.f \
          bc3global.f bc3lhs.f bc3per.f bc3res.f bflux.f \
          e3.f e3b.f e3bvar.f e3dc.f e3ivar.f e3lhs.f e3q.f e3ql.f \
          e3qvar.f e3res.f e3source.f e3stab.f e3sts.f elmStats.f \
          elmgmr.f errsmooth.f filters.f forces.f ftools.f genlmass.f \
          getdiff.f getfld.f hessian.f iclear.f itrPC.f itrbc.f \
          itrdrv.f lesSparse.f rstatic.f soldir.f solfar.f \
          ramg_driver.f \
          elmdist.f 

SRCS  = $(FSRCS) $(CSRCS) $(CXXSRCS)

OBJS  = $(FSRCS:.f=.$(OBJECTEXT)) $(CSRCS:.c=.$(OBJECTEXT)) $(CXXSRCS:.cxx=.$(OBJECTEXT))

$(TARGET_FLOWSOLVER): $(OBJS)
	$(LINK_EXE)$(TARGET_FLOWSOLVER) $(GLOBAL_LFLAGS) \
						 $(OBJS) $(METIS_LIBS) \
						 $(LESLIB_LIBS) \
						 $(MPI_LIBS) \
						 $(F90_LIBS) \
						 $(CXX_LIBS) \
						 $(TOP)/lib/lib_simvascular_solverio.$(STATICEXT)
	mv $(TARGET_FLOWSOLVER) $(TOP)/bin

.PHONY : clean
clean : 
	-$(RM) $(OBJS) $(TARGET_FLOWSOLVER) *.mod *.exe *.pdb *.manifest
	-$(RM) $(TOP)/bin/$(TARGET_FLOWSOLVER) $(TOP)/bin/$(TARGET_FLOWSOLVER).manifest
