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
   calculate the shape functions and their derivatives for
   a triangular face
*/

int TriShapeAndDrv(int p,double par[2],double N[],double dN[][2]){
  int i,j,nshp=0;
  double L[3];
  
  if(p != 2)
    return nshp;
  L[0]=par[0];
  L[1]=par[1];
  L[2]=1.0e0-par[0]-par[1];
    
  /* define shape functions for a quadratic triangle */

  /* collect all vertex modes */
  for(i=0; i<3; i++) {
    N[i] = L[i];
    if(i==2)
      dN[i][0]=dN[i][1]=-1.0e0;
    else {
      for(j=0; j<2; j++) {
	if(i==j)
	  dN[i][j] = 1.0e0;
	else
	  dN[i][j] = 0.0e0;
      }
    }
  }
  nshp=3;
  if( p > 1 ){
    /* collect edge modes (only quadratic for now) */
    N[3] = -2.0e0*L[0]*L[1];
    N[4] = -2.0e0*L[1]*L[2];
    N[5] = -2.0e0*L[0]*L[2];
    dN[3][0] = -2.0e0*L[1];
    dN[3][1] = -2.0e0*L[0];
    dN[4][0] = 2.0e0*L[1];
    dN[4][1] = -2.0e0+2.0e0*L[0]+4.0e0*L[1];
    dN[5][0] = -2.0e0+4.0e0*L[0]+2.0e0*L[1];
    dN[5][1] = 2.0e0*L[0];
    nshp=6;
  }
  return nshp;
}
    
    
  
  
