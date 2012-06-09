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

/* Functions to generate the Legandre Polynomials and their derivatives */

extern "C" double LP(int j, double x);
extern "C" double LPdrv(int j, double x);
extern "C" double phi(int p, double x);
extern "C" double phiDrv(int p,double x);
extern "C" int HexShapeAndDrv(int p, double par[3], double N[], double
			      dN[][3]);
extern "C" int WedgeShapeAndDrv(int p, double Inputpar[3], double N[], double 
				dN[][3]);
extern "C" int PyrShapeAndDrv (int p, double Inputpar[3], double N[], double 
				dN[][3]);

/* Blending functions */

extern "C" double Line_eB(double xi1);
extern "C" double dLEBdxi1(double xi1);
extern "C" double dLEBdxi2(double xi1);
extern "C" double dLEBdxi3(double xi1);

extern "C" double Quad_eB(double xi1, double xi2, int sign);
extern "C" double dQEBdxi1(double xi1, double xi2, int sign);
extern "C" double dQEBdxi2(double xi1, double xi2, int sign);
extern "C" double dQEBdxi3(double xi1, double xi2, int sign);

extern "C" double Quad_fB(double xi1, double xi2);
extern "C" double dQFBdxi1(double xi1, double xi2);
extern "C" double dQFBdxi2(double xi1, double xi2);
extern "C" double dQFBdxi3(double xi1, double xi2);

extern "C" double Hex_eB(double xi[3], int sign2, int sign3);
extern "C" double dHEBdxi1(double xi[3], int sign2, int sign3);
extern "C" double dHEBdxi2(double xi[3], int sign2, int sign3);
extern "C" double dHEBdxi3(double xi[3], int sign2, int sign3);

extern "C" double Hex_fB(double xi[3], int sign3);
extern "C" double dHFBdxi1(double xi[3], int sign3);
extern "C" double dHFBdxi2(double xi[3], int sign3);
extern "C" double dHFBdxi3(double xi[3], int sign3);

/* Entity Level functions */

extern "C" int mesh_edge(double xi1,int gOrd[3], int p, double* entfn,
			 double** edrv);
extern "C" int quad_face(double xi[3], int gOrd[3], int p, double*
			 entfn, double** edrv); 
extern "C" int hex_regn(double xi[3], int p, double*
			entfn, double** edrv);

