#ifndef __SUPERLU_INTERFACE /* allow multiple inclusions */
#define __SUPERLU_INTERFACE

#include "slu_zdefs.h"

/*! \brief Supernodal LU factor related */
extern void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);

extern void    dgstrf (superlu_options_t*, SuperMatrix*,
                       int, int, int*, void *, int, int *, int *,
                       SuperMatrix *, SuperMatrix *, SuperLUStat_t*, int *);

extern void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                        SuperMatrix *, SuperLUStat_t*, int *);

extern int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

#endif /* __SUPERLU_INTERFACE */
