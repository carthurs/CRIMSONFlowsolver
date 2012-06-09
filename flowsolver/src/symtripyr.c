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

int triIntPntPyr(int, DARR3**,double**);

#ifdef sun4_5
symtripyr_(int *n1, double pt[][4], double wt[], int *err)
#endif
#ifdef ibm6000
symtripyr(int *n1, double pt[][4], double wt[], int *err)
#endif
#ifdef sgi
void symtripyr_(int *n1, double pt[][4], double wt[], int *err)
#endif
#ifdef decalp
void symtripyr_(int *n1, double pt[][4], double wt[], int *err)
#endif
#ifdef intel
void  SYMTRIPYR(int *n1, double pt[][4], double wt[], int *err)
#endif
{
  double *lwt;
  DARR3 *lpt;
  int i,j;
  *err = triIntPntPyr(*n1, &lpt, &lwt);
  for(i=0; i < *n1; i++){
    wt[i] = lwt[i];
    for(j=0; j < 3; j++)
      pt[i][j] = lpt[i][j];
  }
}

#include <stdio.h>
/* 1pt constants */
#define zero 0.0000000000000000
#define ppt500 0.5000000000000000
#define p2pt00 2.0000000000000000

/* 4pt constants */
#define ppt455 0.4553418012614798
#define ppt788 0.7886751345948130
#define ppt577 0.5773502691896260
#define ppt122 0.1220084679281462
#define ppt211 0.2113248654051870

/* Rule 3 */
#define mpt68 -0.6872983346207413
#define pt68   0.6872983346207413
#define mpt77 -0.7745966692414831
#define pt77   0.7745966692414831
#define mpt38 -0.3872983346207416
#define pt38   0.3872983346207416
#define mpt08 -0.08729833462074179
#define pt08   0.08729833462074179
#define mpt5  -0.5000000000000000

#define mpt11 -0.1127016653792584
#define mpt88 -0.8872983346207415

#define Qw311  0.308641975308642
#define Qw321  0.493827160493831
#define Qw322  0.790123456790124
#define Qw312  0.493827160493831
#define Qw313  0.308641975308642
#define Qw323  0.493827160493831
#define Qw331  0.308641975308642
#define Qw332  0.493827160493831
#define Qw333  0.308641975308642

/*  #define J 0.108103018168070 */
/*  #define K 0.445948490915965 */

/*  #define L 0.225000000000000 */
/*  #define M 0.125939180544827 */
/*  #define N 0.132394152788506 */

/*  #define O 0.797426985353087 */
/*  #define P 0.101286507323456 */
/*  #define Q 0.470142064105115 */
/*  #define R 0.059715871789770 */

/*  #define S 0.050844906370207 */
/*  #define T 0.116786275726379 */
/*  #define U 0.082851075618374 */

/*  #define V 0.873821971016996 */
/*  #define W 0.063089014491502 */
/*  #define X 0.501426509658179 */
/*  #define Y 0.249286745170910 */
/*  #define Z 0.636502499121399 */
/*  #define AA 0.310352451033785 */
/*  #define BB 0.053145049844816 */

/*  typedef double DARR3[3] ; */

static double rst1[][3] = {{zero,-1,zero}};
static double wt1[] = {2.2360679774997897};
/* only rule 2 is implemented correctly at this time */
static double rst4[][3] ={
{-.455341801261479547, -.788675134594812883 ,-.577350269189625763},
{ .455341801261479547, -.788675134594812883 ,-.577350269189625763},
{-.122008467928146216, -.211324865405187118 , .577350269189625764},
{ .122008467928146216, -.211324865405187118 , .577350269189625764}};
 
static double wt4[] = {0.394337567297406440, 0.394337567297406440, 
                       0.105662432702593559, 0.105662432702593559};

static double rst9[][3] = {
                           {mpt68, mpt88, mpt77},
                           { zero, mpt88, mpt77},
                           { pt68, mpt88, mpt77},
                           {mpt38,  mpt5,  zero},
                           { zero,  mpt5,  zero},
                           { pt38,  mpt5,  zero},
                           {mpt08, mpt11,  pt77},
                           { zero, mpt11,  pt77},
                           { pt08, mpt11,  pt77}
};
static double wt9[] = {Qw311, Qw321, Qw331, Qw312, Qw322, Qw332, Qw313, 
			 Qw323, Qw333};

/*  static double rst7[][3] = {{A,A,A}, */
/*                               {O,P,P},{P,O,P},{P,P,O}, */
/*                               {Q,Q,R},{Q,R,Q},{R,Q,Q}}; */
/*  static double wt7[] = {L,M,M,M,N,N,N}; */

/*  static double rst12[][3] = {{V,W,1.0-V-W},{W,V,1.0-V-W},{W,W,1.0-W-W}, */
/*  			    {X,Y,1.0-X-Y},{Y,X,1.0-X-Y},{Y,Y,1.0-Y-Y}, */
/*  			    {Z,AA,1.0-Z-AA},{AA,Z,1.0-Z-AA},{Z,BB,1.0-Z-BB}, */
/*  			    {AA,BB,1.0-AA-BB},{BB,AA,1.0-BB-AA},{BB,Z,1.0-BB-Z} */

/*  }; */

/*  static double wt12[] = {S,S,S,T,T,T,U,U,U,U,U,U}; */
                             
int triIntPntPyr(int nint, DARR3 **bcord, double **wt)
{
  int retval = 1 ;
  if( nint == 1 ){*bcord = rst1 ; *wt = wt1; }
  else if( nint == 4 ){*bcord = rst4 ; *wt = wt4; }
  else if( nint == 9 ){*bcord = rst9 ; *wt = wt9; }
    else
  {
/*      fprintf(stderr,"\n%d integration points unsupported in symtri.c; give {1,3,4,6,7,12}\n",nint); */
    fprintf(stderr,"\n%d integration points unsupported in symtri.c; give {1,4}\n",nint);
    retval = 0;
  }
  return retval ;
}



