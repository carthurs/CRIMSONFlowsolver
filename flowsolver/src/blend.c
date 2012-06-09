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

This function evaluates the blend funtion for a given mesh entity over
another mesh entity

*/

#include <math.h>
#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

double V_blendOnEntity(int vid, int etype, double *L) {
  /* blend a vertex mode on a mesh entity */
  /* vid = local vertex index 
     etype = entity type
     L = par. coords
  */
  if( etype == Sedge )
    return V_blendIndexedOnEdge(vid,L);
  else if(etype == Stri || etype == Stet)
    return V_blendIndexed(vid,L);
  else
    return 0.0e0;
}

double V_blendIndexed(int i, double *L) {
  return L[i] ;
}

double V_blendIndexedOnEdge(int i, double *L) {
  if( i == 0 )
    return 0.5*(1.0-L[0]);
  else
    return 0.5*(1.0+L[0]);
}

double E_blendOnFace(int eindex[], int etype, double *L) {
  /* blend a vertex mode on a mesh face */
  /* vid = local vertex index 
     etype = entity type
     L = par. coords
  */
  if( etype == Stri)
    return F_edgeBlendTri(eindex, L);
  else if(etype == Squad)
    return F_edgeBlendQuad(eindex, L);
  else
    return 0.0e0;
}

double F_edgeBlendTri(int index[2], double *L) {
  return -2.0*L[index[0]]*L[index[1]];
}

double F_edgeBlendQuad(int *index, double *L) {
  return 0.0;
}

double E_blendOnRegion(int eindex[], int etype, double *L) {
  /* blend a mesh edge mode on a tetra. region */
 if( etype == Stet ) 
    return R_edgeBlendTet(eindex, L);   
 else
  return 0.0;
}

double R_edgeBlendTet(int index[2], double *L) {
  return -2.0*L[index[0]]*L[index[1]];
}

double F_blendOnRegion(int index[], int etype, double *L) {
  /* blend a face mode on a tet. region */
  if( etype == Stet ) {
    return L[index[0]]*L[index[1]]*L[index[2]] ;
  } else
    return 0.0;
}

#ifdef __cplusplus
}
#endif
