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

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void TetQuadLag(double *par, double *N,double dN[][3]);

int TetShapeAndDrvL(int p,double par[3],double N[],double dN[][3]) {
  /*
  static int TetEMAP[6][2]={{0,1},{1,2},{2,0},{0,3},{1,3},{2,3}};
  static int TetFMAP[4][3]={{0,1,2},{0,3,1},{1,3,2},{0,2,3}};
  */
  static int Nshp[3]={4,10,20};
  if( p == 2) {
    TetQuadLag(par,N,dN);
  } else {
    fprintf(stderr,"p != 2 Not implemented for lagrange\n");    
  }
  return Nshp[p];
}

void TetQuadLag(double *par, double *N,double dN[][3]) {
  double r=par[0],s=par[1],t=par[2],t8,t4;
  t8 = 1.0-r-s-t;
  N[0] = (2.0*r-1.0)*r;
  N[1] = (2.0*s-1.0)*s;
  N[2] = (2.0*t-1.0)*t;
  N[3] = (1.0-2.0*r-2.0*s-2.0*t)*t8;
  N[4] = 4.0*r*s;
  N[5] = 4.0*s*t;
  N[6] = 4.0*t*r;
  N[7] = 4.0*r*t8;
  N[8] = 4.0*s*t8;
  N[9] = 4.0*t*t8;

  t4 = -3.0+4.0*r+4.0*s+4.0*t;
  dN[0][0] = 4.0*r-1.0;
  dN[0][1] = 0.0;
  dN[0][2] = 0.0;
  dN[1][0] = 0.0;
  dN[1][1] = 4.0*s-1.0;
  dN[1][2] = 0.0;
  dN[2][0] = 0.0;
  dN[2][1] = 0.0;
  dN[2][2] = 4.0*t-1.0;
  dN[3][0] = t4;
  dN[3][1] = t4;
  dN[3][2] = t4;
  dN[4][0] = 4.0*s;
  dN[4][1] = 4.0*r;
  dN[4][2] = 0.0;
  dN[5][0] = 0.0;
  dN[5][1] = 4.0*t;
  dN[5][2] = 4.0*s;
  dN[6][0] = 4.0*t;
  dN[6][1] = 0.0;
  dN[6][2] = 4.0*r;
  dN[7][0] = 4.0-8.0*r-4.0*s-4.0*t;
  dN[7][1] = -4.0*r;
  dN[7][2] = -4.0*r;
  dN[8][0] = -4.0*s;
  dN[8][1] = 4.0-4.0*r-8.0*s-4.0*t;
  dN[8][2] = -4.0*s;
  dN[9][0] = -4.0*t;
  dN[9][1] = -4.0*t;
  dN[9][2] = 4.0-4.0*r-4.0*s-8.0*t;
}

#ifdef __cplusplus
}
#endif

