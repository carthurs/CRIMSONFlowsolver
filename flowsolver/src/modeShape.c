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

This function returns hierarchic mode shapes associated with a given mesh 
entity of a specified polynomial order.

*/

#include <math.h>
#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

double E_modeShape(int p, double *L) {
/* Return edge mode shape function evaluated along an edge at a point (Li,Lj) 
for polynomial order p, note returned mode is of polynomial order p-2 */
  if( p < 2 ) 
    return 0.0;
  else
    return En(p-2,L[0],L[1]);
}

double F_modeShapeTri(int p, int i, double *L) {
  int alpha,beta,found,count,nskip;
  /* return i-th triangular face mode of polynomial order p 
     note: there are p-2 modes with polynomial order p */

  if( p < 3 || i < 0 || i > p-3 )
    return 0.0 ;
  
  count = found = 0;
  for(alpha=0; alpha <= p-3; alpha++) {
    for(beta=0; beta <= p-3; beta++) {
      if( alpha+beta == p-3 ) {
        if( count == i )
          found=1;
        else
          count++;
      } 
      if( found )
        break;   
    }
    if( found )
      break;   
  }
  if( found ) 
    return Fn(alpha, beta, L[0], L[1]);
  else
    return 0.0;
}

double F_modeShapeQuad(int p, int i, double *L) {
  /* return ith p-order quad face mode evaluated at L */
  return 0.0;
}

double R_modeShapeTet(int p, int i, double *L) {
  int alpha,beta,gamma,count,found;

  /* return the i-th mode shape of polynomial order p for this region , there are
     (p-2)*(p-3)/2 mode shapes of polynomial order p */
  if( p < 4 || i < 0 || i > (((p-2)*(p-3)/2)-1) )
    return 0.0;

  count = 0;
  found = 0;
  for( alpha=0; alpha <= p-4; alpha++) {
    for( beta=0; beta <= p-4; beta++) {
      for( gamma=0; gamma <= p-4; gamma++) {
        if( alpha+beta+gamma == p-4 ) {
          if( count == i )
            found = 1;
          else
            count++;
	}
        if( found ) 
          break ;
      }
      if( found ) 
        break ;
    }
    if( found ) 
      break ;
  }
  if( found ) 
   return Bn(alpha,beta,gamma,L[0],L[1],L[2]);
  else
    return 0.0;
}

double R_modeShapeHex(int p, int i, double *L) {
  /* ith p-order mode shape for a hex element */
  return 0.0;
}

#ifdef __cplusplus
}
#endif






