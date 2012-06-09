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

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

void brickToTet(double xi,double eta, double zeta,
                double *r, double *s, double *t, double *J);

void quadToTri(double xi,double eta,
               double *r, double *s,double *J);
double quadToTriJac(double xi,double eta);
int GaussLegendre1D(int, double **,double **);

int GaussLegendreTet(int n1,int n2,int n3,double GLrstw[][4],double *GLwt) {
  /* get degenerate n1Xn2Xn3 Gauss-Legendre scheme to integrate over a tet */
  int i,j,k,index=0;
  double *pt1,*pt2,*pt3,*wt1,*wt2,*wt3,dJ;
  const double six=6.000000000000000;

  GaussLegendre1D(n1,&pt1,&wt1);
  GaussLegendre1D(n2,&pt2,&wt2);
  GaussLegendre1D(n3,&pt3,&wt3);
  for(i=0; i < n1; i++) {
    for(j=0; j < n2; j++) {
      for(k=0; k < n3; k++) {
        brickToTet(pt1[i],pt2[j],pt3[k],&GLrstw[index][0],
                   &GLrstw[index][1],&GLrstw[index][2],&dJ);
        GLrstw[index][3] = 1.0e0-GLrstw[index][0]-GLrstw[index][1]
                                -GLrstw[index][2];
        GLwt[index++] = dJ*wt1[i]*wt2[j]*wt3[k]*six;
      }
    }
  }
  return index;
}
int GaussLegendreTri(int n1,int n2,double GLr[][3],double *GLwt) {
  /* get degenerate n1Xn2 Gauss-Legendre scheme to integrate over a tet */
  int i,j,index=0;
  double *pt1,*pt2,*wt1,*wt2,dJ;
  const double two = 2.0000000000000000;

  GaussLegendre1D(n1,&pt1,&wt1);
  GaussLegendre1D(n2,&pt2,&wt2);
  for(i=0; i < n1; i++) {
    for(j=0; j < n2; j++) {
      quadToTri(pt1[i],pt2[j],&GLr[index][0],&GLr[index][1],&dJ);
      GLr[index][2] = 1.0e0-GLr[index][0]-GLr[index][1];
      GLwt[index++] = dJ*wt1[i]*wt2[j]*two;
    }
  }
  return index;
}

void brickToTet(double xi,double eta, double zeta,
                      double *r, double *s, double *t, double *J) {
  double r1,rs1;
  *r = 0.5e0*(1.0e0+xi);
  r1 = 1.0e0-(*r);
  *s = 0.5e0*(1.0e0+eta)*r1;
  rs1 = 1.0e0-(*r)-(*s);
  *t = 0.5e0*(1.0e0+zeta)*rs1;
  *J = 0.125e0*r1*rs1;
}

void quadToTri(double xi,double eta,double *r, double *s, double *J) {
  double r1;
  *r = 0.5e0*(1.0e0+xi);
  r1 = 1.0e0-(*r);
  *s = 0.5e0*(1.0e0+eta)*r1;
  *J = 0.25e0*r1;  
}

double quadToTriJac(double xi,double eta) {
  return 0.125e0*(1.0e0-eta);
}

#ifdef __cplusplus
}
#endif

