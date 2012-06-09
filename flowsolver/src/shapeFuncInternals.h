/*********************************************************************

Copyright (c) 2000-2007, Stanford University, 
    Rensselaer Polytechnic Institute, Kenneth E. Jansen, 
    Charles A. Taylor (see SimVascular Acknowledgements file 
    for additional contributors to the source code).

All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions 
are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in the 
documentation and/or other materials provided with the distribution. 
Neither the name of the Stanford University or Rensselaer Polytechnic
Institute nor the names of its contributors may be used to endorse or
promote products derived from this software without specific prior 
written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

**********************************************************************/

#ifndef _H_shapeFuncInternals
#define _H_shapeFuncInternals

#ifndef _ELEM_TOP_TYPE
#define _ELEM_TOP_TYPE

#define Svertex 0
#define Sedge   1
#define Stri    2
#define Squad   3
#define Stet    4
#define Shex    5
#define Swedge  6
#define Spyr    7
#endif /*_ELEM_TOP_TYPE*/

#ifdef __cplusplus
extern "C" {
#endif

int F_parDrv(int i, int j, int k, int type, int (**drv)[2]);
int E_parDrv(int i, int j, int type, double drv[][2]);

double V_blendOnEntity(int vid, int etype, double *L);
double V_blendIndexed(int i, double *L);
double V_blendIndexedOnEdge(int i, double *L);
double E_blendOnFace(int eindex[], int etype, double *L);
double F_edgeBlendTri(int index[2], double *L);
double F_edgeBlendQuad(int *index, double *L);
double E_blendOnRegion(int eindex[], int etype, double *L);
double R_edgeBlendTet(int index[2], double *L);
double F_blendOnRegion(int index[], int etype, double *L);

int V_blendOnEntityDrv(int index, int etype, double *L, double mdrv[3]);
int E_blendOnFaceDrv(int index[], int etype, double *L, double bdrv[2]);
int F_edgeBlendTriDrv(int index[2], double *L, double drv[]);
int E_blendOnRegionDrv(int index[], int etype, double *L, double bdrv[3]);
int F_blendOnRegionDrv(int index[], int etype, double *L, double drv[3]);
int R_edgeBlendTetDrv(int *index, double *L, double drv[]);
double E_blendOnRegion(int eindex[], int etype, double *L);
int R_faceBlendTetDrv(int *index, double *L, double drv[]);
int F_modeShapeTriDrv(int p, int i, double *L, double mdrv[2]);
int R_modeShapeTetDrv(int p, int i, double *L, double mdrv[3]);

double E_modeShape(int p, double *L);
double F_modeShapeTri(int p, int i, double *L);
double F_modeShapeQuad(int p, int i, double *L);
double F_modeShapeQuad(int p, int i, double *L);
double R_modeShapeTet(int p, int i, double *L);
double R_modeShapeHex(int p, int i, double *L);


/* these are low level entity functions */
double En(int ip, double r, double s);
int EnDrv(int ip, double r, double s, double *drv);
double Fn(int i, int j, double r, double s);
int FnDrv(int i,int j, double r, double s, double drv[2]);
double Bn(int i, int j, int k, double r, double s, double t);
int BnDrv(int i,int j, int k, double r, double s, double t, double drv[3]);

#ifdef __cplusplus
}
#endif

#endif

