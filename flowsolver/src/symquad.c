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

#include <malloc.h>

typedef double DARR3[3];

int quadIntPnt(int,DARR3**,double **);

#ifdef sun4_5
symquad_(int *n1, double pt[][4], double wt[], int *err)
#elif ibm6000
symquad(int *n1, double pt[][4], double wt[], int *err)
#elif sgi
void symquad_(int *n1, double pt[][4], double wt[], int *err)
#elif decalp
void symquad_(int *n1, double pt[][4], double wt[], int *err)
#elif intel
void  SYMQUAD(int *n1, double pt[][4], double wt[], int *err)
#endif
{
  double *lwt;
  DARR3 *lpt;
  int i,j;
  *err = quadIntPnt(*n1, &lpt, &lwt);
  for(i=0; i < *n1; i++) {
    wt[i] = lwt[i];
    for(j=0; j < 3; j++) 
      pt[i][j]=lpt[i][j];
  }
}

/*$Id: symquad.c,v 1.3 2005/06/29 19:36:04 nwilson Exp $*/
#include <stdio.h>

#define QpNon -1.000000000000000
#define Qp21  -0.577350269189626
#define Qp22   0.577350269189626
#define Qw2    1.000000000000000

#define Qp31  -0.774596669241483
#define Qp32   0.000000000000000
#define Qp33   0.774596669241483
#define Qw31   0.555555555555556
#define Qw32   0.888888888888889
#define Qw33   0.555555555555556
#define Qw311  0.308641975308642
#define Qw321  0.493827160493831
#define Qw322  0.790123456790124
#define Qw312  0.493827160493831
#define Qw313  0.308641975308642
#define Qw323  0.493827160493831
#define Qw331  0.308641975308642
#define Qw332  0.493827160493831
#define Qw333  0.308641975308642


static double rstw4[][3] = {
  {Qp21, Qp21, QpNon},
  {Qp22, Qp21, QpNon},
  {Qp21, Qp22, QpNon},
  {Qp22, Qp22, QpNon}
};

static double twt4[] = {Qw2,Qw2,Qw2,Qw2};

static double rstw9[][3] = {
   { Qp31,  Qp31,  QpNon},  
   { Qp32,  Qp31,  QpNon},
   { Qp33,  Qp31,  QpNon},
   { Qp31,  Qp32,  QpNon},
   { Qp32,  Qp32,  QpNon},
   { Qp33,  Qp32,  QpNon},
   { Qp31,  Qp33,  QpNon},
   { Qp32,  Qp33,  QpNon},
   { Qp33,  Qp33,  QpNon}
};

static double twt9[] = { Qw311, Qw321, Qw331, Qw312, Qw322, Qw332, Qw313, 
			 Qw323, Qw333 };


#ifdef __cplusplus
extern "C" {
#endif

int quadIntPnt(int nint, DARR3 **bcord, double **wt)
{
  int retval = 1;
  int i,j,l;
  DARR3* rstw;
  double* twt;

  /* Rule 4 & 5*/
  /* rule 5 has not been tested */
  double bp4[]={-0.8611363115940526,-0.3399810435848563,0.339981043584856,0.8611363115940526};
  double bp5[]={-0.9061798459386640,-0.5384693101056831,0,0.5384693101056831,0.9061798459386640};
  double bw4[]={0.3478548451374539,0.6521451548625461,0.6521451548625461,0.3478548451374539};
  double bw5[]={0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891};
  
  if( nint == 4 ){*bcord = rstw4 ; *wt = twt4; }
  else if (nint == 9){*bcord = rstw9 ; *wt = twt9; }
  else if (nint == 16){
    
    rstw = (DARR3 *)malloc(16*sizeof(DARR3));
    twt  = (double *)malloc(16*sizeof(double));
    
    i=j=0;
    for(l=0;l<16;l++){
      rstw[l][0]=bp4[i];
      rstw[l][1]=bp4[j];
      rstw[l][2]=QpNon;
      twt[l]= bw4[i]*bw4[j];
      i++;
      if( i == 4) { i=0; j++;}
    }
    *bcord = rstw; *wt = twt;

  }else if(nint == 25){

    rstw = (DARR3 *)malloc(25*sizeof(DARR3));
    twt  = (double *)malloc(25*sizeof(double));
    
    i=j=0;
    for(l=0;l<25;l++){
      rstw[l][0]=bp5[i];
      rstw[l][1]=bp5[j];
      rstw[l][2]=QpNon;
      twt[l]= bw5[i]*bw5[j];
      i++;
      if( i == 5) { i=0; j++;}
    }
    *bcord = rstw; *wt = twt;

  } else {

    fprintf(stderr,"\n%d integration points unsupported in symquad.c; give {4,9,16,25}\n",nint);
    retval = 0;
  }
  return retval ;
}

#ifdef __cplusplus
}
#endif
