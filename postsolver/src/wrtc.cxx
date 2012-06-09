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
#ifndef intel
#include <unistd.h>
#endif
#ifdef intel
#include <direct.h>
#define getcwd _getcwd
#endif
/* 
   The code calling this has glued together the solution and geometry
   files and has passed in the following:

   nv is the number of degrees of freedom per node
   numel is the total number of elements
   ien is the TRUE global connectivity  
       (index local element node number and element
        number to TRUE global node number)
   numnp is the total number of nodes                          
   x is the coordinates (3*numnp)
   q is the solution (nv*numnp)
   nen is the number of nodes in an element
*/

void wrtc_(int *nv,    int *numel, int ien[],
           int *numnp, float x[], float q[], int *nen)
{
  int i,j,count;
  char *eltype;			/* element type */
  FILE *f1 = fopen("connect.bin","wb");
  FILE *f2 = fopen("points.bin","wb");
  FILE *f3 = fopen("data.bin","wb");
  FILE *f4 = fopen("header.dx","w");

  /* element type */
  if (*nen == 4) eltype = "tetrahedra";
  else eltype = "cubes";
     
  /* write the positions array */
  count=0;
  for(i=0; i < *numnp; i++){
    for(j=0; j < 3; j++) {
      fwrite((void *)&(x[count]),sizeof(float),1,f2);
      count++;
    }
  }

  /* write the connections array */
  count =0;
  for(i=0; i < *numel; i++) {
    for(j=0; j < *nen; j++) {
      fwrite((void *)&(ien[count]),sizeof(int),1,f1);
      count++;
    }
  }

/* write the data array, rho,u1,u2,u3,T */
  count = 0;
  for(i=0; i< *numnp; i++){
    for(j=0; j < *nv; j++){
      fwrite((void *)&(q[count]),sizeof(float),1,f3);
      count++;
    }
  }
  
  /*  write the header file */
  fprintf(f4,"object 1 class array type float rank 1 shape 3 items %d msb binary\n",
	  *numnp);
  fprintf(f4," data file %s/points.bin,0 \n",getcwd(NULL,128));
  fprintf(f4," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f4,"object 2 class array type int rank 1 shape %d items %d msb binary\n",
	  *nen,*numel);
  fprintf(f4," data file %s/connect.bin,0 \n",getcwd(NULL,128));
  fprintf(f4," attribute \"element type\" string \"%s\" \n",eltype);
  fprintf(f4," attribute \"ref\" string \"positions\" \n\n");

  fprintf(f4,"object 3 class array type float rank 1 shape %d items %d msb binary\n",
	  *nv,*numnp);
  fprintf(f4," data file %s/data.bin,0 \n",getcwd(NULL,128));
  fprintf(f4," attribute \"dep\" string \"positions\" \n\n");

  fprintf(f4,"object 4 class array type int rank 0 items %d \n",
	  *numel);
  fprintf(f4," data file %s/partit.out \n",getcwd(NULL,128));
  fprintf(f4," attribute \"dep\" string \"connections\" \n\n");

  fprintf(f4,"object \"irregular positions irregular connections\" class field\n");
  fprintf(f4," component \"positions\" value 1\n");
  fprintf(f4," component \"connections\" value 2\n");
  fprintf(f4," component \"data\" value 3\n");
  fprintf(f4,"\n end \n");

  printf("\n Files successfully written... \n");
  printf("\n Number of elements = %d \n",*numel);
  printf(" Number of nodes = %d \n\n",*numnp);


  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
}
