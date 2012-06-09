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
#include <stdio.h>

typedef double DARR4[4];

int hexIntPnt(int,DARR4**,double **);

#ifdef sun4_5
symhex_(int *n1, double pt[][4], double wt[], int *err)
#elif ibm6000
symhex(int *n1, double pt[][4], double wt[], int *err)
#elif sgi
void symhex_(int *n1, double pt[][4], double wt[], int *err)
#elif decalp
void symhex_(int *n1, double pt[][4], double wt[], int *err)
#elif intel
void   SYMHEX(int *n1, double pt[][4], double wt[], int *err)
#endif
{
  double *lwt;
  DARR4 *lpt;
  int i,j;
  *err = hexIntPnt(*n1, &lpt, &lwt);
  for(i=0; i < *n1; i++) {
    wt[i] = lwt[i];
    for(j=0; j <4; j++) 
      pt[i][j]=lpt[i][j];
  }
}

/* these are the rule 2 int points and weights */

#define Qp21 -0.577350269189626
#define Qp22  0.577350269189626
#define Qw2   1.000000000000000

/* these are the rule 3 int points and weights */

#define Qp31  -0.774596669241483
#define Qp32   0.000000000000000
#define Qp33   0.774596669241483
#define Qw31   0.555555555555556
#define Qw32   0.888888888888889
#define Qw33   0.555555555555556
#define Qw311  0.308641975308642
#define Qw321  0.493827160493831
#define Qw322  0.790123456790124


#define Qw3111 0.171467764060361
#define Qw3112 0.274348422496571
#define Qw3113 0.171467764060361
#define Qw3121 0.274348422496571
#define Qw3122 0.438957475994513
#define Qw3123 0.274348422496571
#define Qw3131 0.171467764060361
#define Qw3132 0.274348422496571
#define Qw3133 0.171467764060361

#define Qw3211 0.274348422496571
#define Qw3212 0.438957475994513
#define Qw3213 0.274348422496571
#define Qw3221 0.438957475994513
#define Qw3222 0.702331961591221
#define Qw3223 0.438957475994513
#define Qw3231 0.274348422496571
#define Qw3232 0.438957475994513
#define Qw3233 0.274348422496571

#define Qw3311 0.171467764060361
#define Qw3312 0.274348422496571
#define Qw3313 0.171467764060361
#define Qw3321 0.274348422496571
#define Qw3322 0.438957475994513
#define Qw3323 0.274348422496571
#define Qw3331 0.171467764060361
#define Qw3332 0.274348422496571
#define Qw3333 0.171467764060357

/* Rule 1*/

static double  rstw1[][4] = {
  { 0.0, 0.0, 0.0, 0.0 }
};

static double twt1[] = { 8.0000000 };


/* Rule 2*/


static double rstw8[][4] = {
  {Qp21,Qp21,Qp21,0.0},
  {Qp22,Qp21,Qp21,0.0},
  {Qp21,Qp22,Qp21,0.0},
  {Qp22,Qp22,Qp21,0.0},
  {Qp21,Qp21,Qp22,0.0},
  {Qp22,Qp21,Qp22,0.0},
  {Qp21,Qp22,Qp22,0.0},
  {Qp22,Qp22,Qp22,0.0}
};

static double twt8[] = {Qw2,Qw2,Qw2,Qw2,Qw2,Qw2,Qw2,Qw2};

/* Rule 3*/

static double rstw27[][4] = {
  { Qp31,  Qp31,  Qp31,  0.0 },
  { Qp32,  Qp31,  Qp31,  0.0 },
  { Qp33,  Qp31,  Qp31,  0.0 },
  { Qp31,  Qp32,  Qp31,  0.0 },
  { Qp32,  Qp32,  Qp31,  0.0 },
  { Qp33,  Qp32,  Qp31,  0.0 },
  { Qp31,  Qp33,  Qp31,  0.0 },
  { Qp32,  Qp33,  Qp31,  0.0 },
  { Qp33,  Qp33,  Qp31,  0.0 },
  { Qp31,  Qp31,  Qp32,  0.0 }, 
  { Qp32,  Qp31,  Qp32,  0.0 },
  { Qp33,  Qp31,  Qp32,  0.0 },
  { Qp31,  Qp32,  Qp32,  0.0 },
  { Qp32,  Qp32,  Qp32,  0.0 },
  { Qp33,  Qp32,  Qp32,  0.0 },
  { Qp31,  Qp33,  Qp32,  0.0 },
  { Qp32,  Qp33,  Qp32,  0.0 },
  { Qp33,  Qp33,  Qp32,  0.0 },
  { Qp31,  Qp31,  Qp33,  0.0 },
  { Qp32,  Qp31,  Qp33,  0.0 },
  { Qp33,  Qp31,  Qp33,  0.0 },
  { Qp31,  Qp32,  Qp33,  0.0 },
  { Qp32,  Qp32,  Qp33,  0.0 },
  { Qp33,  Qp32,  Qp33,  0.0 },
  { Qp31,  Qp33,  Qp33,  0.0 },
  { Qp32,  Qp33,  Qp33,  0.0 },
  { Qp33,  Qp33,  Qp33,  0.0 }
};


static double twt27[] = 
{ Qw3111, Qw3211, Qw3311, Qw3121, Qw3221, Qw3321, Qw3131, Qw3231, Qw3331,
  Qw3112, Qw3212, Qw3312, Qw3122, Qw3222, Qw3322, Qw3132, Qw3232, Qw3332,  
  Qw3113, Qw3213, Qw3313, Qw3123, Qw3223, Qw3323, Qw3133, Qw3233, Qw3333 };


               	  
#ifdef __cplusplus
extern "C" {
#endif

int hexIntPnt(int nint, DARR4 **bcord, double **wt)
{
  int retval = 1;
  int i,j,k,l;
  DARR4 *rstw;
  double *twt;

  /* Rule 4 & 5*/
  /* rule 5 has not been tested */ 
  double bp4[]={-0.8611363115940526,-0.3399810435848563,0.339981043584856,0.8611363115940526};
  double bp5[]={-0.9061798459386640,-0.5384693101056831,0,0.5384693101056831,0.9061798459386640};
  double bw4[]={0.3478548451374539,0.6521451548625461,0.6521451548625461,0.3478548451374539};
  double bw5[]={0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891};


  if( nint == 1){ *bcord = rstw1; *wt = twt1;}

  else if( nint == 8 ){*bcord = rstw8 ; *wt = twt8; }

  else if( nint == 27 ) {*bcord = rstw27 ; *wt = twt27; }

  else if( nint == 64 ) {

    rstw = (DARR4 *)malloc(64*sizeof(DARR4));
    twt   = (double *)malloc(64*sizeof(double));

    i=j=k=0;

    for(l=0;l<64;l++){
      rstw[l][0]=bp4[i];
      rstw[l][1]=bp4[j];
      rstw[l][2]=bp4[k];
      rstw[l][3]=0.0;
      twt[l]=bw4[i]*bw4[j]*bw4[k];
      i++;
      if(i == 4){
	i=0 ; j++ ;
	if(j == 4) { j=0; k++;}
      }
    }
    *bcord = rstw; *wt = twt; 

  } else if( nint == 125) {
    
    rstw = (DARR4 *)malloc(125*sizeof(DARR4));
    twt   = (double *)malloc(125*sizeof(double));
    
    i=j=k=0;
    
    for(l=0;l<125;l++){
      rstw[l][0]=bp5[i];
      rstw[l][1]=bp5[j];
      rstw[l][2]=bp5[k];
      rstw[l][3]=0.0;
      twt[l]=bw5[i]*bw5[j]*bw5[k];
      i++;
      if(i == 5){
	i=0 ; j++ ;
	if(j == 5) { j=0; k++;}
      }
    }
    *bcord = rstw; *wt = twt; 

  } else {

    fprintf(stderr,"\n%d integration points unsupported; give {8,27,64,125,.}\n",nint);
    retval = 0;
  }
  return retval ;
}

#ifdef __cplusplus
}
#endif
