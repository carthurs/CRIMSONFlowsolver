#ifndef _SOLVER_FIELDS_H_
#define _SOLVER_FIELDS_H_
//this struct contains pointers that will point to solver solutions fields that are needed by the adaptation routine
typedef struct SolverFields
{
  int * nshg;
  int * ndof;
  int * nvar;
  int * numnp;
  int * numel;
  int * numelb;
  int * numpbc;
  int * nen;
  int * nelblk;
  int * nelblb;
  int * nlwork;

  int * ilwork;
  int * iper;
  int * nBC;
  int * iBCtmp;
  int * iBC;

  double * sol_time;
  double * sol;
  double * err_ind;
  double * xcoord;
  double * BCinp;
  double * BCtmp;
  double * BC;

  int * numflx;
  int * nfaces;
  int * numEBC;
  int * numNBC;
  int * ftype;

} SolverFields;
#endif
