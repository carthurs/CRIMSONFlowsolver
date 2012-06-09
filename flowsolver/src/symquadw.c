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

typedef double DARR3[3];

int quadwIntPnt(int,DARR3**,double **);

#ifdef sun4_5
symquadw_(int *n1, double pt[][4], double wt[], int *err)
#elif ibm6000
symquadw(int *n1, double pt[][4], double wt[], int *err)
#elif sgi
void symquadw_(int *n1, double pt[][4], double wt[], int *err)
#elif decalp
void symquadw_(int *n1, double pt[][4], double wt[], int *err)
#elif intel
void  SYMQUADW(int *n1, double pt[][4], double wt[], int *err)
#endif
{
  double *lwt;
  DARR3 *lpt;
  int i,j;
  *err = quadwIntPnt(*n1, &lpt, &lwt);
  for(i=0; i < *n1; i++) {
    wt[i] = lwt[i];
    for(j=0; j < 3; j++) 
      pt[i][j]=lpt[i][j];
  }
}

/*$Id: symquadw.c,v 1.3 2005/06/29 19:36:05 nwilson Exp $*/
#include <stdio.h>

/* Rule 2 constants */

#define Qp21  -0.577350269189626
#define Qp22   0.577350269189626
#define Pp21  0.211324865405187
#define Pp22  0.788675134594813
#define Qw2   1.000000000000000

/* Rule 3 constants */

#define Qp31  -0.774596669241483
#define Qp32   0.000000000000000
#define Qp33   0.774596669241483
#define Qw311  0.308641975308642   /* Qw31 * Qw31  note: Qw31=Qw33*/
#define Qw321  0.493827160493831   /* Qw32 * Qw31 */
#define Qw322  0.790123456790124   /* Qw32 * Qw32 */
#define Pp31   0.112701665379259
#define Pp32   0.500000000000000
#define Pp33   0.887298334620741

/* Rule 4 constants */

#define Qp41  -0.861136311594053
#define Qp42  -0.339981043584856
#define Qp43   0.339981043584856
#define Qp44   0.861136311594053
#define Qw41   0.347854845137454
#define Qw42   0.652145154862544
#define Qw43   0.652145154862544
#define Qw44   0.347854845137544
#define Pp41   0.873821971016996
#define Pp42   0.063089014491502
#define Pp43   0.501426509658179
#define Pp44   0.249286745170910
#define Qw4141 0.121002993285602 /* Qw41 * Qw41 */
#define Qw4142 0.226851851851851 /* Qw41 * Qw42 */
#define Qw4143 0.226851851851851 /* Qw41 * Qw43 */
#define Qw4144 0.121002993285602 /* Qw41 * Qw44 */
#define Qw4241 0.226851851851851 /* etc..       */
#define Qw4242 0.425293303010692
#define Qw4243 0.425293303010692
#define Qw4244 0.226851851851851
#define Qw4341 0.226851851851851
#define Qw4342 0.425293303010692
#define Qw4343 0.425293303010692
#define Qw4344 0.226851851851851
#define Qw4441 0.121002993285602
#define Qw4442 0.226851851851851
#define Qw4443 0.226851851851851
#define Qw4444 0.121002993285602


/* typedef double DARR3[3] ; */

/* Rule 2 */

static double rstw4[][3] = {
  {Pp21, 0.0 , Qp21},
  {Pp22, 0.0 , Qp21},
  {Pp21, 0.0 , Qp22},
  {Pp22, 0.0 , Qp22}
};

static double twt4[] = {Qw2,Qw2,Qw2,Qw2};

/* Rule 3 */

static double rstw9[][3] = {
  {Pp31,0.0,Qp31},
  {Pp32,0.0,Qp31},
  {Pp33,0.0,Qp31},
  {Pp31,0.0,Qp32},
  {Pp32,0.0,Qp32},
  {Pp33,0.0,Qp32},
  {Pp31,0.0,Qp33},
  {Pp32,0.0,Qp33},
  {Pp33,0.0,Qp33}
};

static double twt9[] = 
{Qw311, Qw321, Qw311, Qw321, Qw322, Qw321, Qw311, Qw321, Qw311}; 

/* Rule 4 */

static double rstw16[][3] = {
  {Pp41,0.0,Qp41},
  {Pp42,0.0,Qp41},
  {Pp43,0.0,Qp41},
  {Pp44,0.0,Qp41},
  {Pp41,0.0,Qp42},
  {Pp42,0.0,Qp42},
  {Pp43,0.0,Qp42},
  {Pp44,0.0,Qp42},
  {Pp41,0.0,Qp43},
  {Pp42,0.0,Qp43},
  {Pp43,0.0,Qp43},
  {Pp44,0.0,Qp43},
  {Pp41,0.0,Qp44},
  {Pp42,0.0,Qp44},
  {Pp43,0.0,Qp44},
  {Pp44,0.0,Qp44}
};

static double twt16[] =
{Qw4141, Qw4241, Qw4341, Qw4441, Qw4142, Qw4242, Qw4342, Qw4442,
 Qw4143, Qw4243, Qw4343, Qw4443, Qw4144, Qw4244, Qw4344, Qw4444};


#ifdef __cplusplus
extern "C" {
#endif

int quadwIntPnt(int nint, DARR3 **bcord, double **wt)
{
  int retval = 1;

  if( nint == 4 ){*bcord = rstw4 ; *wt = twt4; }
  else if( nint == 9 ){*bcord = rstw9 ; *wt = twt9; }
  else if( nint == 16){*bcord = rstw16 ; *wt = twt16;}
  else
  {
    fprintf(stderr,"\n%d integration points unsupported in symquadw.c; give {4,9,16}\n",nint);
    retval = 0;
  }
  return retval ;
}

#ifdef __cplusplus
}
#endif
