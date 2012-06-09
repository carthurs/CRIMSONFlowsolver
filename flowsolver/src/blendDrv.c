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

/*

This function evaluate derivatives of entity blend functions over other 
entities.

*/

#include <math.h>
#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

int V_blendOnEntityDrv(int index, int etype, double *L, double mdrv[3]) {
  /* derivative of blend a vertex mode on a mesh entity */
  if( etype == Sedge ) {
    if(index == 0) /* N = 0.5(1-zi) */
      mdrv[0] = -0.5;
    else           /* N = 0.5(1+zi) */
      mdrv[0] = +0.5;
    return 1 ;
  } else if( etype == Stri ) {
    mdrv[0] = mdrv[1] = 0.0;
    if( index == 2 )
      mdrv[0] = mdrv[1] = -1.0;
    else
      mdrv[index] = 1.0;
    return 1 ;
  } else if ( etype == Stet ) {
    mdrv[0] = mdrv[1] = mdrv[2] = 0.0;
    if( index == 3 ) 
      mdrv[0] = mdrv[1] = mdrv[2] = -1.0;
    else
      mdrv[index] = 1.0;
    return 1;
  } else
    return 0;  
}

int E_blendOnFaceDrv(int index[], int etype, double *L, double bdrv[2]) {
  if( etype == Stri ) 
    return F_edgeBlendTriDrv(index, L, bdrv);
  else
    return 0;
}

int F_edgeBlendTriDrv(int index[2], double *L, double drv[]) {
  double r,s;

  drv[0] = drv[1] = 0.0;
  r = L[0] ; s = L[1] ;

  /* figure out which edge we are dealing with */
   /* v0=0, V1=1 */
   if( (index[0]==0 && index[1]==1) || (index[0]==1 && index[1]==0) ) {
      drv[0] = -2.0*s;
      drv[1] = -2.0*r;
   /* v0=1, V1=2 */
   } else if( (index[0]==1 && index[1]==2) || (index[0]==2 && index[1]==1) ) {
      drv[0] = 2.0*s;
      drv[1] = -2.0+2.0*r+4.0*s;
   /* v0=2, V1=0 */
   } else if( (index[0]==2 && index[1]==0) || (index[0]==0 && index[1]==2) ) {
      drv[0] = 4.0*r-2.0+2.0*s;
      drv[1] = 2.0*r;
    } else
      return 0;

  return 1 ;
}

int E_blendOnRegionDrv(int index[], int etype, double *L, double bdrv[3]) {
  if( etype == 4 )
      return R_edgeBlendTetDrv(index, L, bdrv);
  else
    return 0;
}

int F_blendOnRegionDrv(int index[], int etype, double *L, double drv[3]) {
  /* blend a face mode on a region */
  drv[0] = drv[1] = drv[2] = 0.0;
  if(etype==Stet )
    return 1;
  else
    return 0;
}

int R_edgeBlendTetDrv(int *index, double *L, double drv[]) {
  double r=L[0],s=L[1],t=L[2];

  /* the blend is given by -2*L[i]*L[j] */
   /* v0=0, V1=1 */
   if( (index[0]==0 && index[1]==1) || (index[0]==1 && index[1]==0) ) {
      drv[0] = -2.0*s;
      drv[1] = -2.0*r;
      drv[2] = 0.0;
   /* v0=0, V1=3 */
   } else if( (index[0]==0 && index[1]==3) || (index[0]==3 && index[1]==0) ) {
      drv[0] = -2.0+4.0*r+2.0*s+2.0*t;
      drv[1] = 2.0*r;
      drv[2] = 2.0*r;
   /* v0=1, V1=2 */
   } else if( (index[0]==1 && index[1]==2) || (index[0]==2 && index[1]==1) ) {
      drv[0] = 0.0;
      drv[1] = -2.0*t;
      drv[2] = -2.0*s;
   /* v0=1, V1=3 */
   } else if( (index[0]==1 && index[1]==3) || (index[0]==3 && index[1]==1) ) {
      drv[0] = 2.0*s;
      drv[1] = -2.0+2.0*r+4.0*s+2.0*t;
      drv[2] = 2.0*s;
   /* v0=2, V1=0 */
   } else if( (index[0]==2 && index[1]==0) || (index[0]==0 && index[1]==2) ) {
      drv[0] = -2.0*t;
      drv[1] = 0.0;
      drv[2] = -2.0*r;
   /* v0=2, V1=3 */
   } else if( (index[0]==2 && index[1]==3) || (index[0]==3 && index[1]==2) ) {
      drv[0] = 2.0*t;
      drv[1] = 2.0*t;
      drv[2] = -2.0+2.0*r+2.0*s+4.0*t;
    } else
      return 0;

  return 1;
}

int R_faceBlendTetDrv(int *index, double *L, double drv[]) {
  int isum = index[0]+index[1]+index[2];

  drv[0]=drv[1]=drv[2]=0;
  if( isum == 3) /* r*s*t */ {
    drv[0] = L[1]*L[2] ;
    drv[1] = L[0]*L[2] ;
    drv[2] = L[1]*L[0] ;
  } else if( isum == 5) /* r*t*(1-r-s-t) */ {
    drv[0] = L[2]*(1.0-2.0*L[0]-L[1]-L[2]);
    drv[1] = -L[0]*L[2] ;
    drv[2] = L[0]*(1.0-L[0]-L[1]-2.0*L[2]);
  } else if( isum == 4) /* r*s*(1-r-s-t) */ {
    drv[0] = L[1]*(1.0-2.0*L[0]-L[1]-L[2]);
    drv[1] = L[0]*(1.0-L[0]-2.0*L[1]-L[2]);
    drv[2] = -L[0]*L[1] ;
  } else if( isum == 6) /* s*t*(1-r-s-t) */ {
    drv[0] = -L[1]*L[2] ;
    drv[1] = L[2]*(1.0-L[0]-2.0*L[1]-L[2]);
    drv[2] = L[1]*(1.0-L[0]-L[1]-2.0*L[2]);
  } else
    return 0;
  return 1;
}
#ifdef __cplusplus
}
#endif
