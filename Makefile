TOP=.

include $(TOP)/buildoptions.mk

SUBDIRS = solverio/src 

ifeq ($(MAKE_WITH_PRESOLVER),1)
   SUBDIRS += presolver/src
endif

ifeq ($(MAKE_WITH_POSTSOLVER),1)
   SUBDIRS += postsolver/src
endif

#ifeq ($(MAKE_WITH_MEMLS),1)
#   SUBDIRS += memLS/src
#endif

ifeq ($(MAKE_WITH_FLOWSOLVER),1)
   SUBDIRS += flowsolver/src
endif

ifeq ($(MAKE_WITH_ADAPTOR),1)
   SUBDIRS += lu/src meshSimShapeFunction/src solution/src adapt/src
endif

#ifeq ($(MAKE_WITH_MAIN),1)
#   SUBDIRS += forwardmain/src
#endif

ifeq ($(MAKE_WITH_ESTIMATOR),1)
   SUBDIRS += estimation/src
endif

all:
	@for i in ${SUBDIRS}; do ( \
	  cd $$i; \
	  $(MAKE)) ; done

clean:
	for i in ${SUBDIRS}; do ( \
	  cd $$i; \
	  $(MAKE) clean ) ; done
