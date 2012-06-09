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
  
This function returns the derivative of the parametric transformation
from element to entity coordinate systems.

*/

#include <math.h>
#include "shapeFuncInternals.h"

#ifdef __cplusplus
extern "C" {
#endif

static int SF_par_TABLE[][3][2] = {
   {{0,0},{0,0},{0,0}},
   {{1,0},{0,1},{0,0}},
   {{1,0},{0,0},{0,1}},
   {{1,-1},{0,-1},{0,-1}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,1},{1,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{1,0},{0,1}},
   {{0,-1},{1,-1},{0,-1}}, 
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,1},{0,0},{1,0}},
   {{0,0},{0,1},{1,0}},
   {{0,0},{0,0},{0,0}},
   {{0,-1},{0,-1},{1,-1}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{0,0},{0,0},{0,0}},
   {{-1,1},{-1,0},{-1,0}},
   {{-1,0},{-1,1},{-1,0}},
   {{-1,0},{-1,0},{-1,1}}
};

int E_parDrv(int i, int j, int type, double drv[][2]) {
  int ind[2],ii;
  /* find out type of element */
  if( type == Sedge )
    return 0;
  else if( type == Stri ) {
    drv[0][0] = drv[0][1] = drv[1][0] = drv[1][1] = 0.0;
    if( i == 0 && j == 1 ) { /* edge +0 , xi = s-r */
      drv[0][0] = 1.0;
      drv[1][1] = 1.0;
    } else if( i == 1 && j == 0 ) { /* edge -0 , xi = r-s */
      drv[0][1] = 1.0;
      drv[1][0] = 1.0;
    } else if( i == 1 && j == 2 ) { /* edge +1 , xi = t-s */
      drv[0][1] = -1.0;
      drv[1][0] = 1.0;
      drv[1][1] = -1.0;
    } else if( i == 2 && j == 1 ) { /* edge -1 , xi = s-t */
      drv[0][0] = -1.0;
      drv[1][0] = -1.0;
      drv[1][1] = 1.0;
    } else if( i == 2 && j == 0 ) { /* edge +2 , xi = r-t */
      drv[0][0] = -1.0;
      drv[0][1] = 1.0;
      drv[1][0] = -1.0;
    } else if( i == 0 && j == 2 ) { /* edge -2 , xi = t-r */
      drv[0][0] = 1.0;
      drv[0][1] = -1.0;
      drv[1][1] = -1.0;
    } else
      return 0;
    return 1;
  } else if( type == Stet ) {
    drv[0][0]=drv[0][1]=drv[1][0]=drv[1][1]=drv[2][0]=drv[2][1]=0.0;
    ind[0] = i; ind[1] = j;
    for(ii=0; ii<2; ii++) {
      if(ind[ii] == 0) /* r'= r */
	drv[0][ii] = 1.0;
      else if(ind[ii] == 1) /* r'= s */
	drv[1][ii] = 1.0;
      else if(ind[ii] == 2) /* r'= t */
	drv[2][ii] = 1.0;
      else if(ind[ii] == 3) /* r'= 1-r-s-t */
	drv[0][ii]=drv[1][ii]=drv[2][ii]=-1.0;
    }
    return 1;
  } else
    return 0;
}

int F_parDrv(int i, int j, int k, int type, int (**drv)[2]) {
  int index = 10*i + j;
  *drv = SF_par_TABLE[index] ;
  return 1;
}

#ifdef __cplusplus
}
#endif
