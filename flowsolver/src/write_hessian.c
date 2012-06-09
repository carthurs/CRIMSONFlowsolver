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

#ifdef intel
#define write_hessian_ WRITE_HESSIAN
#endif

void
write_hessian_( double* hessian, double* gradian, int* nshg ) {

    FILE* idmap = fopen( "lihessmap.dat","r");
    int* map = ( int* )malloc( (*nshg)*sizeof( int ) );
    int i,j,k;
    int x;
   
    FILE* uhess = fopen( "uhessian.dat", "w" );
    FILE* vhess = fopen( "vhessian.dat", "w" );
    FILE* whess = fopen( "whessian.dat", "w" );


    /* FILE* ugrad = fopen( "ugradian.dat", "w" ); */
    /* FILE* vgrad = fopen( "vgradian.dat", "w" ); */
    /* FILE* wgrad = fopen( "wgradian.dat", "w" ); */

    /* int ug[3] = { 1, 4, 7 }; */
    /* int vg[3] = { 2, 5, 8 }; */
    /* int wg[3] = { 3, 6, 9 }; */
    int u[6] = { 1, 10, 19, 13, 22, 25 };
    int v[6] = { 2, 11, 20, 14, 23, 26 };
    int w[6] = { 3, 12, 21, 15, 24, 27 };

    for( x=0; x < *nshg; x++ )  fscanf(idmap,"%d\n", map+x );
    fclose( idmap );

    fprintf( uhess,"%d\n", *nshg );
    fprintf( vhess,"%d\n", *nshg );
    fprintf( whess,"%d\n", *nshg );
    /* fprintf( ugrad,"%d\n", *nshg ); */
    /* fprintf( vgrad,"%d\n", *nshg ); */
    /* fprintf( wgrad,"%d\n", *nshg ); */

    for( i=0; i< *nshg; i++ ) {

        k = map[ i ]-1;

       /* for( j=0; j<3; j++ ) { */
         /* fprintf( ugrad, "%f ", gradian[i+(ug[j]-1)*(*nshg)] ); */
         /* fprintf( vgrad, "%f ", gradian[i+(vg[j]-1)*(*nshg)] ); */
         /* fprintf( wgrad, "%f ", gradian[i+(wg[j]-1)*(*nshg)] ); */
       /* } */
         
       for( j=0; j<6; j++ ) {
         fprintf( uhess, "%f ", hessian[k+(u[j]-1)*(*nshg)] );
         fprintf( vhess, "%f ", hessian[k+(v[j]-1)*(*nshg)] );
         fprintf( whess, "%f ", hessian[k+(w[j]-1)*(*nshg)] );
        }
       fprintf( uhess, "\n" ) ;
       fprintf( vhess, "\n" ) ;
       fprintf( whess, "\n" ) ;
       /* fprintf( ugrad, "\n" ) ; */
       /* fprintf( vgrad, "\n" ) ; */
       /* fprintf( wgrad, "\n" ) ; */
    }          
  
    free( map );
    fclose( uhess );
    fclose( vhess );
    fclose( whess );
    /* fclose( ugrad ); */
    /* fclose( vgrad ); */
    /* fclose( wgrad ); */
}
