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

#include <iostream>
#include <stdio.h>
#include <string>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "cvSolverIO.h"
#ifndef intel
#include <unistd.h>
#include <strings.h>
#endif
#ifdef intel
#include <string.h>
void  bzero_old(void* ptr, size_t sz) {
	int i;
	char *cptr;
	cptr = (char*) ptr;
	for (i=0; i < sz; i++) {
		cptr[i]=0;
	}
	return;
}
#endif

using namespace std;

void wrtc_(int *nv,    int *numel, int ien[],
		int *numnp, float x[], float q[], int *nen);
extern FILE* fgeombc;
extern FILE* frest;

/* 

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



int main(int argc, char* argv[])
{
	int stepnumber,iorig,iv;
	FILE *nstart,*volcheck;
	FILE* fgeombc;
	FILE* frest;
	int **ncorp2d, *ient, *ien;
	int nshgl, numvar;
	int maxnshg, numprocs , nshgtot,i,j,k,l, lstep;
	int np, numel, nen, nelblk;
	int neltot, nelsofar, neltp, nenl, ipordl, nshl;
	int gnum, opnum, nendx, junique, nshlmin;
	int iarray[10];

	int nshgl_bf,numvar_bf,lstep_bf;
	int nshgl_disp,numvar_disp,lstep_disp;
	int nshgl_wss,numvar_wss,lstep_wss;
	int nshgl_ybar,numvar_ybar,lstep_ybar;
	int nshgl_res,numvar_res,lstep_res;
	int newstepnumber = 0;

	double *qglobal, *qlocal, *xlocal, *xglobal, *aglobal, *fglobal, *dglobal, *dglobal_ref, *distglobal, *wglobal;
	double *yglobal;
	double *resglobal, *relativeVelocityGlobal, *updatedMeshCoordinatesGlobal; // ALE variables KDL, MA 2016
	float *qdx, *xdx;
	int *iendx, nodes[8];
	char rfname[40];
	char gfname[40];
	char fname1[255], rfile[255];
	char iotype[80];
	int intfromfile[50];
	int ione=1, itwo=2, ithree=3,  iseven=7;
	int igeom;
	int irstin;
	int ierr;
	int ixsiz, nsd, iientsiz, iqsiz;

	double qmax[20],qmin[20];

	double e1[3],e2[3],e3[3];
	double vol;
	int vtxnum0,vtxnum1,vtxnum2,vtxnum3;
	/* variables used in processing the commandline */
	int iarg, arglength;
	string tmpstr;
	bool StepNumberAvailable;
	int indxu1, indxu2, indxx1, indxx2;

	char visfn[1024];
	char visfnmesh[1024];

	for(i=0; i< 20; i++){
		qmax[i]=-100000000;
		qmin[i]=100000000;
	}

	igeom=1;
	irstin=2;

	iotype[0]='\0';
	sprintf(iotype,"%s","binary");

	/* BEGIN PROCESSING COMMAND-LINE ARGUMENTS */
	/* Assume the command line is okay */
	bool BogusCmdLine = false;
	/* Assume no options specified at command line */
	bool RequestedHelp = false;
	bool RequestedStepNumber = false;
	bool RequestedAcceleration = false;
	bool RequestedVolCheck = false;
	bool RequestedDX = false;
	bool RequestedASCII = false;
	bool RequestedVIS = false;
	bool RequestedVISmesh = false;
	bool RequestedPhasta = false;
	bool RequestedBoundaryFluxes = false;
	bool RequestedDisplacements  = false;
	bool RequestedDistances  = false;
	bool RequestedWSS  = false;
	bool RequestedYbar = false;
	bool RequestedNewSn = false;
	bool aleOn = false;

	/* argc is the number of strings on the command-line */
	/*  starting with the program name */
	for(iarg=1; iarg<argc; iarg++){
		arglength = strlen(argv[iarg]);
		/* replace 0..arglength-1 with argv[iarg] */
		tmpstr.replace(0,arglength,argv[iarg],0,arglength);
		if(tmpstr=="-h"){
			RequestedHelp = true;
			cout << endl;
			cout << "usage:" <<endl;
			cout << "  Reduce [-sn stepnumber][-ph][-td][-vol][-dx][-ASCII]" << endl;
			cout << endl;
			cout << "COMMAND-LINE ARGUMENT SUMMARY" << endl;
			cout << "  -h                  : Display usage and command-line argument summary"<< endl;
			cout << "  -sn stepnumber      : Specify step number to reduce"<< endl;
			cout << "  -ph                 : Write phasta-format file restart.<stepnumber>.0"<<endl;
			cout << "  -td                 : Reduce time-derivative field"<< endl;
			cout << "  -vol                : Calculate element volumes and write volcheck.dat" <<endl;
			cout << "  -dx                 : Write files for visualization in Data Explorer" <<endl;
			cout << "  -ASCII              : Write ASCII-format restart and geometry files" <<endl;
			cout << "  -vis filename       : Write Vis-format results file" <<endl;
			cout << "  -vismesh filename   : Write Vis-format mesh file" <<endl;
			cout << "  -nonbinary          : Read/Write files in ASCII format" <<endl;
			cout << "  -bflux              : reduce Boundary Fluxes"<<endl;
			cout << "  -disp               : reduce displacements"<<endl;
			cout << "  -dist               : reduce distances"<<endl;
			cout << "  -wss                : reduce Wall Shear Stress"<<endl;
			cout << "  -ybar               : reduce ybar field"<<endl;
			cout << "  -newsn stepnumber   : override step number in file"<<endl;
			cout << "  -ale                : reduce ALE variables in file"<<endl;
			cout << "END COMMAND-LINE ARGUMENT SUMMARY" << endl;
			cout << endl;
		}
		else if(tmpstr=="-sn"){
			RequestedStepNumber = true;
			iarg++;
			stepnumber = atoi(argv[iarg]);
			StepNumberAvailable = true;
		}
		else if(tmpstr=="-newsn"){
			RequestedNewSn = true;
			iarg++;
			newstepnumber = atoi(argv[iarg]);
		}
		else if(tmpstr=="-td"){
			RequestedAcceleration = true;
		}
		else if(tmpstr=="-dx"){
			RequestedDX = true;
		}
		else if(tmpstr=="-ph"){
			RequestedPhasta = true;
		}
		else if(tmpstr=="-vol"){
			RequestedVolCheck = true;
		}
		else if(tmpstr=="-ASCII"){
			RequestedASCII = true;
		}
		else if(tmpstr=="-nonbinary"){
			iotype[0]='\0';
			sprintf(iotype,"%s","ascii");
		}
		else if(tmpstr=="-vis"){
			RequestedVIS = true;
			iarg++;
			visfn[0]='\0';
			sprintf(visfn,"%s",argv[iarg]);
		}
		else if(tmpstr=="-vismesh"){
			RequestedVISmesh = true;
			iarg++;
			visfnmesh[0]='\0';
			sprintf(visfnmesh,"%s",argv[iarg]);
		}
		else if(tmpstr=="-bflux"){
			RequestedBoundaryFluxes = true;
		}
		else if(tmpstr=="-disp"){
			RequestedDisplacements = true;
		}
		else if(tmpstr=="-dist"){
			RequestedDistances = true;
		}
		else if(tmpstr=="-wss"){
			RequestedWSS = true;
		}
		else if(tmpstr=="-ybar"){
			RequestedYbar = true;
		}
		else if(tmpstr=="-ALE"){
			aleOn = true;
		}
        else if(tmpstr=="-ale"){
			aleOn = true;
		}		
		else {
			BogusCmdLine = true;
		}
		/* reset tmpstr for next argument */
		tmpstr.erase(0,arglength);
	}
	if(RequestedStepNumber){
		cout << "Will use requested step number " << stepnumber << endl;
	} else {
		/* User did not make a request */
		/* Try to use the number in numstart.dat */
		nstart = fopen("numstart.dat","r");
		if(nstart){
			/* The file numstart.dat is present*/
			fscanf(nstart,"%d", &stepnumber);
			fclose(nstart);
			StepNumberAvailable=true;
			cout << "Will use numstart.dat step number " << stepnumber << endl;
		}
	}
	if(RequestedAcceleration){
		cout << "Will reduce time-derivative field as requested" << endl;
	}
	if(RequestedDX){
		cout << "Will write files for visualization in Data Explorer as requested" << endl;
	}
	if(RequestedVolCheck){
		cout << "Will calculate element volumes and write volcheck.dat as requested" << endl;
	}
	if(RequestedASCII){
		cout << "Will write ASCII restart and geometry as requested" << endl;
	}
	if(RequestedBoundaryFluxes){
		cout << "Will reduce boundary flux field as requested" << endl;
	}
	if(RequestedDisplacements){
		cout << "Will reduce displacement field as requested" << endl;
	}
	if(RequestedDistances){
		cout << "Will reduce distances field as requested" << endl;
	}
	if(RequestedWSS){
		cout << "Will reduce wall shear stress field as requested" << endl;
	}
    if(aleOn){
		cout << "Will reduce ALE fields as requested" << endl;
	}	
	/* Keep these last */
	if(!StepNumberAvailable){
		cout << "No step number or range of steps given, so exiting." << endl;
		return(0);
	}
	if(RequestedHelp){
		cout << endl;
		cout << "Exiting before performing any of the above operations due to -h flag";
		cout << endl;
		return(0);
	}
	/* END PROCESSING COMMAND-LINE ARGUMENTS */


	/* scanning restart.<stepnum>.1 for numvar */
	sprintf(rfname,"restart.%d.1",stepnumber);
	cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
	openfile_(rfname, "read", &irstin   );
	readheader_(&irstin,"solution",(void*)intfromfile,&ithree,"double",iotype);
	numvar=intfromfile[1];  /* pushing the second int into numvar */
	closefile_( &irstin, "read" ); //the check below is not implemented-OK?
	sprintf(gfname,"geombc.dat.%d",1); /* geometry and bc database */
	cout << "Opening " << gfname << " to scan for number of processors and global modes..."<<endl;
	openfile_(gfname, "read", &igeom);
	readheader_(&igeom,"number of processors",(void*)iarray,&ione,"integer",iotype);
	numprocs=iarray[0];
	readheader_(&igeom,"number of global modes",(void*)iarray,&ione,"integer",iotype);
	nshgtot = iarray[0];
	closefile_(&igeom, "read");

	/* scanning geom.dat.<procnum> to add up neltot and to determine maxnshg */
	neltot=0;
	maxnshg=0;
	for(i=0; i< numprocs; i++){
		bzero_old( (void*)gfname, 40);
		sprintf(gfname,"geombc.dat.%d",i+1); /* geometry and bc database */
		openfile_(gfname, "read", &igeom);
		readheader_(&igeom,"number of interior elements",(void*)iarray,&ione,"integer",iotype);
		numel=iarray[0];
		readheader_(&igeom,"maximum number of element nodes",(void*)iarray,&ione,"integer",iotype);
		nen=iarray[0];
		readheader_(&igeom,"number of modes",(void*)iarray,&ione,"integer",iotype);
		nshgl=iarray[0];
		if(nshgl>maxnshg) maxnshg=nshgl;
		closefile_(&igeom, "read" );
		neltot+=numel;
	}

	/* now that we have our size of array information we can allocate
     our  local and global geometry arrays (here local means on a
     given processor and global means the total array assembled 
     across all processors) */

	xglobal = (double *) malloc( 3*nshgtot * sizeof(double));
	xlocal = (double *) malloc( 3*maxnshg *sizeof(double));
	ncorp2d = (int ** )malloc(sizeof(int *)*numprocs);
	nendx=nen;
	if(nen>4) nendx=8;  /* DX thinks of things as ALL tets or hexes ONLY */
	ien = (int *) malloc( nendx*neltot * sizeof(int));
	nshlmin=nendx;
	nelsofar=0;

	if(RequestedVolCheck){
		volcheck = fopen("volcheck.dat","w");
		fprintf(volcheck,"volume  row of x's, then y's then z's \n");
	}

	/*Next we loop over the processors and read each processors
    geometry database.  Using the ncorp2d array, we reconstruct the
    global geometry structures (coordinates and connectivity) */

	for(i=0; i< numprocs; i++){

		/* open geom file  and read header*/
		sprintf(gfname,"geombc.dat.%d",i+1);
		printf("Reducing : %s \n", gfname);
		openfile_( gfname, "read", &igeom );
		readheader_(&igeom,"number of nodes",(void*)iarray,&ione,"integer", iotype);
		np=iarray[0];
		readheader_(&igeom,"number of modes",(void*)iarray,&ione,"integer", iotype);
		nshgl=iarray[0];
		readheader_(&igeom,"number of interior tpblocks",(void*)iarray,&ione,"integer",iotype);
		nelblk=iarray[0];

		/* read coordinates and fill into global array */
		readheader_(&igeom,"co-ordinates",(void*)intfromfile,&itwo,"double",iotype);
		np=intfromfile[0];
		nsd=intfromfile[1];
		ixsiz=np*nsd;
		readdatablock_(&igeom,"co-ordinates",(void*)xlocal,&ixsiz, "double",iotype);
		/* get the map from partion numbering to global numbering */
		if ( numprocs > 1 ) {
			readheader_(&igeom,"mode number map from partition to global",(void*)intfromfile,&ione,"integer",iotype);
			np=intfromfile[0];
			ncorp2d[i] = (int * )malloc(sizeof(int)*np);
			readdatablock_(&igeom,"mode number map from partition to global",(void*)ncorp2d[i],&np,"integer",iotype);
		} else {
			ncorp2d[i] = (int * )malloc(sizeof(int)*nshgl);
			for(j=0; j< nshgl ; j++)
				ncorp2d[i][j]=j+1;
		}

		/* map it to global numbering */
		for(k=0; k< 3; k++)
			for(j=0; j< nshgl ; j++)
				xglobal[k*nshgtot+ncorp2d[i][j]-1] = xlocal[k*nshgl+j];

		/*read connectivity data */
		for(k=0; k< nelblk; k++){
			/* keyphrase identifying interior connectivity element block */
			readheader_(&igeom,"connectivity interior",(void*)intfromfile,&iseven,"integer",iotype);
			neltp  =intfromfile[0];
			nenl   =intfromfile[1];
			ipordl =intfromfile[2];
			nshl   =intfromfile[3];

			if(nshl < nshlmin) nshlmin=nshl;
			/* allocate the array to the right size */
			ient = (int *) malloc( nshl*neltp * sizeof(int));
			iientsiz=neltp*nshl;
			/* now read the array */
			readdatablock_(&igeom,"connectivity interior",(void*)ient,&iientsiz,"integer", iotype);


			/* Now we need to bring ien to the global numbering */
			for(l=0; l< neltp; l++){
				for(j=0; j<nshl ; j++){
					gnum=nelsofar+j*neltot+l;
					opnum=ient[j*neltp+l];
					ien[gnum]=ncorp2d[i][opnum-1];
				}

				if(RequestedVolCheck){
					if(nshl==4){
						vtxnum3=ient[3*neltp+l]-1;
						vtxnum2=ient[2*neltp+l]-1;
						vtxnum1=ient[neltp+l]-1;
						vtxnum0=ient[l]-1;
						for(j=0; j<3; j++)
							e3[j]= xlocal[j*nshgl+vtxnum3]-xlocal[j*nshgl+vtxnum0];
						for(j=0; j<3; j++)
							e2[j]= xlocal[j*nshgl+vtxnum2]-xlocal[j*nshgl+vtxnum0];
						for(j=0; j<3; j++)
							e1[j]= xlocal[j*nshgl+vtxnum1]-xlocal[j*nshgl+vtxnum0];

						vol=-(e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])
								-e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])
								+e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]))/6;

						if(vol<1.0e-10){
							fprintf(volcheck,"%10.8e \n",vol);
							for (j=0; j<3; j++)
								fprintf(volcheck,"%10.8e %10.8e %10.8e %10.8e \n",
										xlocal[j*nshgl+vtxnum0],
										xlocal[j*nshgl+vtxnum1],
										xlocal[j*nshgl+vtxnum2],
										xlocal[j*nshgl+vtxnum3]);
						}
					}
				}
				/* pad to largest nendx since dx has to treat pyr and wedg as
	   degenerate hex */
				opnum=ien[gnum];  /* hijack opnum to keep last real node */
				for(j=nshl; j<nendx ; j++){
					gnum=nelsofar+j*neltot+l;
					ien[gnum]=opnum;
				}
			}
			nelsofar+=neltp;
			free(ient);
		}
		closefile_( &igeom, "read" );
	}
	free(xlocal);

	qglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
	qlocal = (double *) malloc(numvar*maxnshg *sizeof(double));

	if(RequestedAcceleration){
		aglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
	}
	if(RequestedBoundaryFluxes){
		fglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
	}
	if(RequestedDisplacements){
		dglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
		dglobal_ref = (double *) malloc( numvar*nshgtot * sizeof(double));
	}
	if(RequestedDistances){
		distglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
	}
	if(RequestedWSS){
		wglobal = (double *) malloc( numvar*nshgtot * sizeof(double));
	}

	if(aleOn){
#if DEBUG_ALE == 1	
		resglobal = (double *) malloc( 4*nshgtot * sizeof(double));
#endif
		relativeVelocityGlobal = (double *) malloc( 3*nshgtot * sizeof(double));
		updatedMeshCoordinatesGlobal = (double *) malloc( 3*nshgtot * sizeof(double));		
	}

	if (RequestedYbar) {
		/* scanning restart.<stepnum>.1 for numvar */
		sprintf(rfname,"restart.%d.1",stepnumber);
		cout << "Opening " << rfname << " to scan for number of variables..."<<endl;
		openfile_(rfname, "read", &irstin   );
		iarray[0] = -1;
		readheader_(&irstin,"ybar",(void*)iarray,&ithree,"double",iotype);
		if(iarray[0]==-1){
			cout << "No ybar in restart." << stepnumber <<".1 -- Exiting"<<endl;
			return 1;
		}
		numvar_ybar=iarray[1];  /* pushing the second int into numvar */
		closefile_( &irstin, "read" ); //the check below is not implemented-OK?

		yglobal = (double *) malloc( numvar_ybar*nshgtot * sizeof(double));
	}

	/*We loop over the processors and read each processors
      solution database.  Using the ncorp2d array, we reconstruct
      the global solution data */
	for(i=0; i<numprocs; i++){
		/* read in solution for current processor */
		sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
		printf("Reducing : %s \n", rfname);
		openfile_(rfname, "read", &irstin   );
		readheader_(&irstin,"solution",(void*)intfromfile,&ithree,"double",iotype);
		nshgl=intfromfile[0];
		numvar=intfromfile[1];
		lstep=intfromfile[2];
		iqsiz=nshgl*numvar;
		readdatablock_(&irstin,"solution",(void*)qlocal, &iqsiz, "double", iotype);
		closefile_( &irstin, "read" ); //the check below is not implemented-OK?

		/* map solution to global */
		for(k=0; k< numvar; k++){
			for(j=0; j< nshgl ; j++){
				qglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl+j];
				// I don't understand why they skip calculating the min/max
				// if you request boundary fluxes -nate
				//if(!(RequestedBoundaryFluxes)){
				if(qmax[k] < qlocal[k*nshgl+j]) qmax[k]=qlocal[k*nshgl+j];
				if(qmin[k] > qlocal[k*nshgl+j]) qmin[k]=qlocal[k*nshgl+j];
				//} // commented out by nate 11-10-03
			}
		}
		/* read in time derivative of solution for current processor */
		if(RequestedAcceleration){
			for(k=0;k<numvar*maxnshg;k++){ // assigning first zero acceleration
				qlocal[k]=0.0;
			}
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s \n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"time derivative of solution",(void*)intfromfile,&ithree,"double",iotype);
			nshgl=intfromfile[0];
			numvar=intfromfile[1];
			lstep=intfromfile[2];
			iqsiz=nshgl*numvar;
			readdatablock_(&irstin,"time derivative of solution",(void*)qlocal,&iqsiz, "double", iotype);
			closefile_(&irstin, "read");


			/* map solution to global */
			for(k=0; k< numvar; k++){
				for(j=0; j< nshgl ; j++){
					aglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl+j];
				}
			}
		}

		if(RequestedBoundaryFluxes){
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for boundary fluxes\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"boundary flux",(void*)intfromfile,&ithree,"double",iotype);
			nshgl_bf=intfromfile[0];
			numvar_bf=intfromfile[1];
			lstep_bf=intfromfile[2];
			iqsiz=nshgl_bf*numvar_bf;
			readdatablock_(&irstin,"boundary flux",(void*)qlocal,&iqsiz, "double", iotype);
			closefile_(&irstin, "read");
			/* map solution to global */
			for(k=0; k< numvar; k++){
				for(j=0; j< nshgl_bf ; j++){
					fglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_bf+j];
					//if(qmax[k] < qlocal[k*nshgl+j]) qmax[k]=qlocal[k*nshgl+j];
					//if(qmin[k] > qlocal[k*nshgl+j]) qmin[k]=qlocal[k*nshgl+j];
				}
			}
		}

		if(RequestedDisplacements){
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for displacements\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"displacement",(void*)intfromfile,&ithree,"double",iotype);
			nshgl_disp=intfromfile[0];
			numvar_disp=intfromfile[1];
			lstep_disp=intfromfile[2];
			iqsiz=nshgl_disp*numvar_disp;
			readdatablock_(&irstin,"displacement",(void*)qlocal,&iqsiz, "double", iotype);
			closefile_(&irstin, "read");
			/* map solution to global */
			for(k=0; k< numvar; k++){
				for(j=0; j< nshgl_disp ; j++){
					dglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_disp+j];
				}
			}

			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"displacement_ref",(void*)intfromfile,&ithree,"double",iotype);
			nshgl_disp=intfromfile[0];
			numvar_disp=intfromfile[1];
			lstep_disp=intfromfile[2];
			iqsiz=nshgl_disp*numvar_disp;
			readdatablock_(&irstin,"displacement_ref",(void*)qlocal,&iqsiz, "double", iotype);
			closefile_(&irstin, "read");
			/* map solution to global */
			for(k=0; k< numvar; k++){
				for(j=0; j< nshgl_disp ; j++){
					dglobal_ref[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_disp+j];
				}
			}

		}

		if(RequestedDistances){
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for distances\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"distances",(void*)intfromfile,&ithree,"double",iotype);
			nshgl_disp=intfromfile[0];
			numvar_disp=intfromfile[1];
			lstep_disp=intfromfile[2];
			iqsiz=nshgl_disp*numvar_disp;
			readdatablock_(&irstin,"distances",(void*)qlocal,&iqsiz, "double", iotype);
			closefile_(&irstin, "read");
			/* map solution to global */
			for(k=0; k< numvar; k++){
				for(j=0; j< nshgl_disp ; j++){
					distglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_disp+j];
				}
			}
		}

		if(RequestedWSS){
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for wall shear stresses\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"wall shear stresses",(void*)intfromfile,&ithree,"double",iotype);
			nshgl_wss=intfromfile[0];
			numvar_wss=intfromfile[1];
			lstep_wss=intfromfile[2];
			iqsiz=nshgl_wss*numvar_wss;
			readdatablock_(&irstin,"wall shear stresses",(void*)qlocal,&iqsiz, "double", iotype);
			closefile_(&irstin, "read");
			/* map solution to global */
			for(k=0; k< numvar; k++){
				for(j=0; j< nshgl_wss ; j++){
					wglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_wss+j];
				}
			}
		}

		/* read in ybar field for current processor */
		if(RequestedYbar){
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for ybar\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"ybar",(void*)iarray,&ithree,"double",iotype);
			nshgl_ybar=iarray[0];
			numvar_ybar=iarray[1];
			lstep_ybar=iarray[2];
			iqsiz=nshgl_ybar*numvar_ybar;
			readdatablock_(&irstin,"ybar",(void*)qlocal, &iqsiz, "double", iotype);
			closefile_( &irstin, "read" );
			/* map solution to global */
			for(k=0; k< numvar_ybar; k++){
				for(j=0; j< nshgl_ybar ; j++){
					yglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_ybar+j];
				}
			}
		}

		if(aleOn){

#if DEBUG_ALE == 1		
			// read in residual field for current processor 
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for residual\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"residual",(void*)iarray,&ithree,"double",iotype);
			nshgl_res=iarray[0];
			numvar_res=iarray[1];
			lstep_res=iarray[2];
			iqsiz=nshgl_res*numvar_res;
			readdatablock_(&irstin,"residual",(void*)qlocal, &iqsiz, "double", iotype);
			closefile_( &irstin, "read" );
			/* map solution to global */
			for(k=0; k< numvar_res; k++){
				for(j=0; j< nshgl_res ; j++){
					resglobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_res+j];
				}
			}
#endif

			// read in relative velocity field for current processor 
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for relative velocity\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"relative velocity",(void*)iarray,&ithree,"double",iotype);
			nshgl_res=iarray[0];
			numvar_res=iarray[1];
			lstep_res=iarray[2];
			iqsiz=nshgl_res*numvar_res;
			readdatablock_(&irstin,"relative velocity",(void*)qlocal, &iqsiz, "double", iotype);
			closefile_( &irstin, "read" );
			/* map solution to global */
			for(k=0; k< numvar_res; k++){
				for(j=0; j< nshgl_res ; j++){
					relativeVelocityGlobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_res+j];
				}
			}

			// read in updated mesh coordinate field for current processor 
			sprintf(rfname,"restart.%d.%d",stepnumber, i+1);
			printf("Reducing : %s for updated mesh coordinates\n", rfname);
			openfile_(rfname, "read", &irstin   );
			readheader_(&irstin,"updated mesh coordinates",(void*)iarray,&ithree,"double",iotype);
			nshgl_res=iarray[0];
			numvar_res=iarray[1];
			lstep_res=iarray[2];
			iqsiz=nshgl_res*numvar_res;
			readdatablock_(&irstin,"updated mesh coordinates",(void*)qlocal, &iqsiz, "double", iotype);
			closefile_( &irstin, "read" );
			/* map solution to global */
			for(k=0; k< numvar_res; k++){
				for(j=0; j< nshgl_res ; j++){
					updatedMeshCoordinatesGlobal[k*nshgtot+ncorp2d[i][j]-1] = qlocal[k*nshgl_res+j];
				}
			}

		}

	}
	for(k=0; k< numvar; k++)
		printf("var.=%d, min=%f, max=%f \n", k,qmin[k],qmax[k]);

	/* Write out the results */
	int magic_number = 362436;
	int* mptr = &magic_number;

	if(RequestedPhasta){
		/* write the total restart into restart.lstep.0  NOTE 0 is not
       used in our partitioned system which goes from 1..nproc */
		i=-1;

		// allow user to reset the step number
		if (RequestedNewSn) {
			lstep = newstepnumber;
		}

		sprintf(rfile,"restart.%d.%d",lstep,i+1);
		openfile_(rfile, "write" , &irstin);

		writestring_( &irstin,"# PHASTA Input File Version 2.0\n");
		writestring_( &irstin, "# Byte Order Magic Number : 362436 \n");

		bzero_old( (void*)rfile, 255 );
		sprintf(rfile,"# Output generated by phPost version 2.25:  \n");
		writestring_( &irstin, rfile );

		int one=1;
		int size = 1;
		int nitems = 1;
		iarray[ 0 ] = 1;
		writeheader_( &irstin, "byteorder magic number ",(void*)iarray, &nitems, &size, "integer", iotype );

		writedatablock_( &irstin, "byteorder magic number ",(void*)mptr, &nitems, "integer", iotype );

		bzero_old( (void*)rfile, 255 );
		sprintf(rfile,"number of modes : < 0 > %d\n", nshgtot);
		writestring_( &irstin, rfile );

		bzero_old( (void*)rfile, 255 );
		sprintf(rfile,"number of variables : < 0 > %d\n", numvar);
		writestring_( &irstin, rfile );

		size =  numvar*nshgtot;
		nitems = 3;
		iarray[ 0 ] = nshgtot;
		iarray[ 1 ] = numvar;
		iarray[ 2 ] = lstep;

		writeheader_( &irstin, "solution ",
				( void* )iarray, &nitems, &size,"double", iotype);

		nitems = numvar*nshgtot;
		writedatablock_( &irstin, "solution ",
				( void* )(qglobal), &nitems, "double", iotype );
		// if the wall shear stress is requested, write it to the file
		if(RequestedWSS) {
			size =  numvar_wss*nshgtot;
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = numvar_wss;
			iarray[ 2 ] = lstep;

			writeheader_( &irstin, "wall shear stresses ",
					( void* )iarray, &nitems, &size,"double", iotype);

			nitems = numvar_wss*nshgtot;
			writedatablock_( &irstin, "wall shear stresses ",
					( void* )(wglobal), &nitems, "double", iotype );
		}
		// if the traction is requested, write it to the file
		if(RequestedBoundaryFluxes) {
			size =  numvar_bf*nshgtot;
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = numvar_bf;
			iarray[ 2 ] = lstep;

			writeheader_( &irstin, "boundary fluxes ",
					( void* )iarray, &nitems, &size,"double", iotype);

			nitems = numvar_bf*nshgtot;
			writedatablock_( &irstin, "boundary fluxes ",
					( void* )(fglobal), &nitems, "double", iotype );
		}
		/* finished writing the solution in restart */

		// if the displacement is requested, write it to the file
		if(RequestedDisplacements) {
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = nsd;
			iarray[ 2 ] = lstep;
			size = nsd*nshgtot;
			writeheader_( &irstin, "displacement ",
					( void* )iarray, &nitems, &size,"double", iotype );
			nitems = size;
			writedatablock_( &irstin, "displacement ",
					( void* )(dglobal), &nitems, "double", iotype );

			nitems = 3;
			writeheader_( &irstin, "displacement_ref ",
					( void* )iarray, &nitems, &size,"double", iotype );
			nitems = size;
			writedatablock_( &irstin, "displacement_ref ",
					( void* )(dglobal_ref), &nitems, "double", iotype );
		}

		// if the distances are requested, write it to the file
		if(RequestedDistances) {
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = nsd;
			iarray[ 2 ] = lstep;
			size = nsd*nshgtot;
			writeheader_( &irstin, "distances ",
					( void* )iarray, &nitems, &size,"double", iotype );
			nitems = size;
			writedatablock_( &irstin, "distances ",
					( void* )(distglobal), &nitems, "double", iotype );
		}

        if (aleOn)
        {
        
			// write out global residuals
#if DEBUG_ALE == 1		
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = 4;
			iarray[ 2 ] = lstep;
			size = 4*nshgtot;
			writeheader_( &irstin, "residual ",
						( void* )iarray, &nitems, &size,"double", iotype );
			nitems = size;
			writedatablock_( &irstin, "residual ",
						( void* )(resglobal), &nitems, "double", iotype );
#endif

			// write out relative velocity
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = 3;
			iarray[ 2 ] = lstep;
			size = 3*nshgtot;
			writeheader_( &irstin, "relative velocity ",
						( void* )iarray, &nitems, &size,"double", iotype );
			nitems = size;
			writedatablock_( &irstin, "relative velocity ",
						( void* )(relativeVelocityGlobal), &nitems, "double", iotype );

			// write out updated mesh coordinates
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = 3;
			iarray[ 2 ] = lstep;
			size = 3*nshgtot;
			writeheader_( &irstin, "updated mesh coordinates ",
						( void* )iarray, &nitems, &size,"double", iotype );
			nitems = size;
			writedatablock_( &irstin, "updated mesh coordinates ",
						( void* )(updatedMeshCoordinatesGlobal), &nitems, "double", iotype );			
		}

		// if the acceleration is requested, before closing the file write the acceleration
		if(RequestedAcceleration){
			size =  numvar*nshgtot;
			nitems = 3;
			iarray[ 0 ] = nshgtot;
			iarray[ 1 ] = numvar;
			iarray[ 2 ] = lstep;

			writeheader_( &irstin, "time derivative of solution ",
					( void* )iarray, &nitems, &size,"double", iotype);

			nitems = numvar*nshgtot;
			writedatablock_( &irstin, "time derivative of solution ",
					( void* )(aglobal), &nitems, "double", iotype );

			closefile_( &irstin, "write" );

		} else {  // no acceleration is requested, then simply close the restart file
			closefile_( &irstin, "write" );
		}
	}
	free(qlocal);

	/* the balance of the file is specific to a particular
     postprocessor (IBM Data Explorer) format that we must convert
     our information into.  Other postprocessors have different
     requirements that will require changes */


	for(i=0; i< numprocs; i++) free(ncorp2d[i]);
	free(ncorp2d);

	if(RequestedDX){
		/* transpose and float data for  dx */
		xdx = (float *) malloc( 3*nshgtot * sizeof(float));
		qdx = (float *) malloc( numvar*nshgtot * sizeof(float));
		iendx = (int *) malloc( nendx*neltot * sizeof(int));


		for(i=0; i< nshgtot; i++)
			for(j=0; j < numvar; j++)
				qdx[i*numvar+j]=(float)qglobal[j*nshgtot+i];

		for(i=0; i< nshgtot; i++)
			for(j=0; j < 3; j++)
				xdx[i*3+j]=(float)xglobal[j*nshgtot+i];

		/* connectivity in dx is rather strange. */

		if (nendx==nshlmin) {
			for(i=0; i< neltot; i++)
				for(j=0; j < nendx; j++)
					iendx[i*nendx+j]=ien[j*neltot+i]-1;
		}
		/* multitopology */
		else {
			for(i=0; i< neltot; i++){

				junique=4;
				for(j=4; j < nendx; j++)
					if(ien[j*neltot+i] != ien[(j-1)*neltot+i]) junique=j+1;

				/* is this a tet */

				if(junique==4 && nendx==8) {
					iendx[i*nendx+7]=ien[3*neltot+i]-1;
					iendx[i*nendx+6]=ien[3*neltot+i]-1;
					iendx[i*nendx+5]=ien[3*neltot+i]-1;
					iendx[i*nendx+4]=ien[3*neltot+i]-1;
					iendx[i*nendx+3]=ien[1*neltot+i]-1;
					iendx[i*nendx+2]=ien[2*neltot+i]-1;
					iendx[i*nendx+1]=ien[2*neltot+i]-1;
					iendx[i*nendx+0]=ien[0*neltot+i]-1;
				}

				/* is this a wedge */

				else if(junique==6) {
					iendx[i*nendx+7]=ien[5*neltot+i]-1;
					iendx[i*nendx+6]=ien[5*neltot+i]-1;
					iendx[i*nendx+5]=ien[4*neltot+i]-1;
					iendx[i*nendx+4]=ien[3*neltot+i]-1;
					iendx[i*nendx+3]=ien[2*neltot+i]-1;
					iendx[i*nendx+2]=ien[2*neltot+i]-1;
					iendx[i*nendx+1]=ien[1*neltot+i]-1;
					iendx[i*nendx+0]=ien[0*neltot+i]-1;
				}

				/*        if(junique==8) {   this is a hex */

				else{
					for(j=0; j < nendx; j++)
						iendx[i*nendx+j]=ien[j*neltot+i]-1;
				}
			}
		}


		if(nendx==8) {
			for(i=0; i< neltot; i++){
				nodes[0]=iendx[i*nendx+4];
				nodes[1]=iendx[i*nendx+0];
				nodes[2]=iendx[i*nendx+7];
				nodes[3]=iendx[i*nendx+3];
				nodes[4]=iendx[i*nendx+5];
				nodes[5]=iendx[i*nendx+1];
				nodes[6]=iendx[i*nendx+6];
				nodes[7]=iendx[i*nendx+2];
				for(j=0; j < nendx; j++)
					iendx[i*nendx+j]=nodes[j];
			}
		}
		wrtc_(&numvar, &neltot, iendx,
				&nshgtot, xdx, qdx,&nendx);

		free(xdx);
		free(qdx);
		free(iendx);
	}
	if(RequestedASCII){
		//      if(nshgtot<1000) {
		/* echo out restart results */

		frest = fopen("restart.asc.out","w");
		fprintf(frest,"%d %d %d \n", nshgtot, numvar, lstep);

		for(i=0; i< nshgtot; i++){
			for(j=0; j < numvar; j++)
				fprintf(frest,"%f  ",qglobal[j*nshgtot+i]);
			fprintf(frest,"\n");
		}

		fclose(frest);

		/* echo out geom results */
		fgeombc = fopen(gfname,"r");


		fgeombc = fopen("geombc.asc.out","w");
		fprintf(fgeombc, "%d %d %d %d %d %d \n",
				nshgtot, nshgtot, nsd, neltot, nen, nelblk);

		for(i=0; i< nshgtot; i++){
			for(j=0; j < 3; j++)
				fprintf(fgeombc,"%22.7e  ",xglobal[j*nshgtot+i]);
			fprintf(fgeombc,"\n");
		}

		fprintf(fgeombc, "%d %d %d %d \n",
				neltot, nendx, ipordl, nshl);
		for(i=0; i< neltot; i++){
			for(j=0; j < nendx; j++)
				fprintf(fgeombc,"%d  ",ien[j*neltot+i]);
			fprintf(fgeombc,"\n");
		}
		fclose(fgeombc);
		//      }
	}

	// write out mesh in Spectrum Vis format
	// NOTE: this code assumes there is only 1 material region.
	// NOTE: assumes only linear tets!

	if (RequestedVISmesh) {

		frest = fopen(visfnmesh,"w");

		// output the nodes
		char s[80];
		fprintf(frest, "problem  \"%s\"  \n", "phasta mesh");
		fprintf(frest, "  time information \n");
		fprintf(frest, "    number of time steps %d \n", 1);
		fprintf(frest, "    time steps \n");
		fprintf(frest, "    %d  %g \n", 1, 1.0);
		fprintf(frest, "    end time steps \n");
		fprintf(frest, "  end time information \n\n");
		s[0]='\0';
		sprintf (s, "region_%d", 1);
		fprintf(frest, "  region  \"%s\" \n", s);
		fprintf(frest, "    nature \"%s\" \n", "solid");

		fprintf(frest, "    number of nodal coordinates %d \n", nshgtot);
		fprintf(frest, "    nodal coordinates \n");

		for(i=0; i< nshgtot; i++){
			fprintf(frest,"%i ",i+1);
			for(j=0; j < 3; j++) fprintf(frest,"%22.7le  ",xglobal[j*nshgtot+i]);
			fprintf(frest,"\n");
		}

		fprintf(frest, "    end node coordinates \n");

		// output the elements

		fprintf(frest, "    element set \"%s\"  \n", "eset_1");
		fprintf(frest, "      material id 0  \n");
		fprintf(frest, "      nodes per element %d \n", 4);
		fprintf(frest, "      topology \"%s\" \n", "tet");
		fprintf(frest, "      number of elements %d \n", neltot);
		fprintf(frest, "      connectivity \n" );

		for(i=0; i< neltot; i++){
			fprintf(frest,"%i ", i+1);
			for(j=0; j < nendx; j++) fprintf(frest,"%d  ",ien[j*neltot+i]);
			fprintf(frest,"\n");
		}


		fprintf(frest, "      end connectivity \n");
		fprintf(frest, "    end element set \n");
		fprintf(frest, "  end region \n\n");
		fprintf(frest, "end problem  \n");
		fclose(frest);

	}


	// write out Spectrum Vis format
	if(RequestedVIS) {
		int numNodes = nshgtot;
		frest = fopen(visfn,"w");

		// write out initial header
		fprintf(frest,"  region \"fluid_region\"\n");
		fprintf(frest,"    time step 1\n");
		fprintf(frest,"    time 1.0\n");
		fprintf(frest,"\n");

		// pressure
		fprintf(frest,"    analysis results \"pressure\"\n");
		fprintf(frest,"      number of data %i\n",numNodes);
		fprintf(frest, "      type \"nodal\"\n");
		fprintf(frest, "      order \"scalar\"\n");
		fprintf(frest, "      length 1\n");
		fprintf(frest, "      data\n");
		for(i=0; i< nshgtot; i++){
			fprintf(frest,"%25.15le \n",qglobal[i]);
		}
		fprintf(frest, "      end data\n");
		fprintf(frest, "    end analysis results\n");
		fprintf(frest, "\n");

		// velocity
		fprintf(frest, "    analysis results \"velocity\"\n");
		fprintf(frest, "      number of data %i\n",numNodes);
		fprintf(frest, "      type \"nodal\"\n");
		fprintf(frest, "      order \"vector\"\n");
		fprintf(frest, "      number of components 3\n");
		fprintf(frest, "      components\n");
		fprintf(frest, "      \"x\"\n");
		fprintf(frest, "      \"y\"\n");
		fprintf(frest, "      \"z\"\n");
		fprintf(frest, "      end components\n");
		fprintf(frest, "      length 3\n");
		fprintf(frest, "      data\n");
		for (i=0; i< nshgtot; i++) {
			for (j=1; j < 4; j++) {
				fprintf(frest,"%25.15le ",qglobal[j*nshgtot+i]);
			}
			fprintf(frest,"\n");
		}
		fprintf(frest, "      end data\n");
		fprintf(frest, "    end analysis results\n");
		fprintf(frest, "\n");

		if (numvar > 5) {
			// transport
			fprintf(frest,"    analysis results \"transport\"\n");
			fprintf(frest,"      number of data %i\n",numNodes);
			fprintf(frest, "      type \"nodal\"\n");
			fprintf(frest, "      order \"scalar\"\n");
			fprintf(frest, "      length 1\n");
			fprintf(frest, "      data\n");
			for(i=0; i< nshgtot; i++){
				fprintf(frest,"%25.15le \n",qglobal[5*nshgtot+i]);
			}
			fprintf(frest, "      end data\n");
			fprintf(frest, "    end analysis results\n");
			fprintf(frest, "\n");
		}

		// traction
		if(RequestedBoundaryFluxes) {
			fprintf(frest, "    analysis results \"traction\"\n");
			fprintf(frest, "      number of data %i\n",numNodes);
			fprintf(frest, "      type \"nodal\"\n");
			fprintf(frest, "      order \"vector\"\n");
			fprintf(frest, "      number of components 3\n");
			fprintf(frest, "      components\n");
			fprintf(frest, "      \"x\"\n");
			fprintf(frest, "      \"y\"\n");
			fprintf(frest, "      \"z\"\n");
			fprintf(frest, "      end components\n");
			fprintf(frest, "      length 3\n");
			fprintf(frest, "      data\n");
			for (i=0; i< nshgtot; i++) {
				for (j=1; j < 4; j++) {
					fprintf(frest,"%25.15le ",fglobal[j*nshgtot+i]);
				}
				fprintf(frest,"\n");
			}
			fprintf(frest, "      end data\n");
			fprintf(frest, "    end analysis results\n");
			fprintf(frest, "\n");
		}

		// displacements
		if(RequestedDisplacements) {
			fprintf(frest, "    analysis results \"displacement\"\n");
			fprintf(frest, "      number of data %i\n",numNodes);
			fprintf(frest, "      type \"nodal\"\n");
			fprintf(frest, "      order \"vector\"\n");
			fprintf(frest, "      number of components 3\n");
			fprintf(frest, "      components\n");
			fprintf(frest, "      \"x\"\n");
			fprintf(frest, "      \"y\"\n");
			fprintf(frest, "      \"z\"\n");
			fprintf(frest, "      end components\n");
			fprintf(frest, "      length 3\n");
			fprintf(frest, "      data\n");
			for (i=0; i< nshgtot; i++) {
				for (j=0; j < 3; j++) {
					fprintf(frest,"%25.15le ",dglobal[j*nshgtot+i]);
				}
				fprintf(frest,"\n");
			}
			fprintf(frest, "      end data\n");
			fprintf(frest, "    end analysis results\n");
			fprintf(frest, "\n");
			fprintf(frest, "    analysis results \"displacement_ref\"\n");
			fprintf(frest, "      number of data %i\n",numNodes);
			fprintf(frest, "      type \"nodal\"\n");
			fprintf(frest, "      order \"vector\"\n");
			fprintf(frest, "      number of components 3\n");
			fprintf(frest, "      components\n");
			fprintf(frest, "      \"x\"\n");
			fprintf(frest, "      \"y\"\n");
			fprintf(frest, "      \"z\"\n");
			fprintf(frest, "      end components\n");
			fprintf(frest, "      length 3\n");
			fprintf(frest, "      data\n");
			for (i=0; i< nshgtot; i++) {
				for (j=0; j < 3; j++) {
					fprintf(frest,"%25.15le ",dglobal_ref[j*nshgtot+i]);
				}
				fprintf(frest,"\n");
			}
			fprintf(frest, "      end data\n");
			fprintf(frest, "    end analysis results\n");
			fprintf(frest, "\n");
		}

		// distances
		if(RequestedDistances) {
			fprintf(frest, "    analysis results \"distances\"\n");
			fprintf(frest, "      number of data %i\n",numNodes);
			fprintf(frest, "      type \"nodal\"\n");
			fprintf(frest, "      order \"scalar\"\n");
			fprintf(frest, "      length 1\n");
			fprintf(frest, "      data\n");
			for (i=0; i< nshgtot; i++) {
				fprintf(frest,"%25.15le ",distglobal[nshgtot+i]);
				fprintf(frest,"\n");
			}
			fprintf(frest, "      end data\n");
			fprintf(frest, "    end analysis results\n");
			fprintf(frest, "\n");
		}

		// Wall Shear Stress
		if(RequestedWSS) {
			fprintf(frest, "    analysis results \"wall shear stress\"\n");
			fprintf(frest, "      number of data %i\n",numNodes);
			fprintf(frest, "      type \"nodal\"\n");
			fprintf(frest, "      order \"vector\"\n");
			fprintf(frest, "      number of components 3\n");
			fprintf(frest, "      components\n");
			fprintf(frest, "      \"x\"\n");
			fprintf(frest, "      \"y\"\n");
			fprintf(frest, "      \"z\"\n");
			fprintf(frest, "      end components\n");
			fprintf(frest, "      length 3\n");
			fprintf(frest, "      data\n");
			for (i=0; i< nshgtot; i++) {
				for (j=0; j < 3; j++) {
					fprintf(frest,"%25.15le ",wglobal[j*nshgtot+i]);
				}
				fprintf(frest,"\n");
			}
			fprintf(frest, "      end data\n");
			fprintf(frest, "    end analysis results\n");
			fprintf(frest, "\n");
		}

		fclose(frest);

	}

	// ybar
	if(RequestedYbar) {

		int size, nitems;

		sprintf(rfile,"%s.%d.0","ybar",lstep_ybar);
		openfile_(rfile, "write" , &irstin);

		writestring_( &irstin,"# PHASTA Input File Version 2.0\n");
		writestring_( &irstin, "# Byte Order Magic Number : 362436 \n");

		bzero_old( (void*)rfile, 255 );
		sprintf(rfile,"# Output generated by phPost version 2.7:  \n");
		writestring_( &irstin, rfile );

		size = 1;
		nitems = 1;
		iarray[0] = 1;
		writeheader_( &irstin, "byteorder magic number ",
				(void*)iarray, &nitems, &size, "integer", iotype );

		writedatablock_( &irstin, "byteorder magic number ",
				(void*)mptr, &nitems, "integer", iotype );

		bzero_old( (void*)rfile, 255 );
		sprintf(rfile,"number of modes : < 0 > %d\n", nshgtot);
		writestring_( &irstin, rfile );

		bzero_old( (void*)rfile, 255 );
		sprintf(rfile,"number of variables : < 0 > %d\n", numvar_ybar);
		writestring_( &irstin, rfile );

		size =  numvar_ybar*nshgtot;
		nitems = 3;
		iarray[0] = nshgtot;
		iarray[1] = numvar_ybar;
		iarray[2] = lstep_ybar;

		writeheader_( &irstin, "ybar",
				(void*)iarray, &nitems, &size,"double", iotype);

		nitems = numvar_ybar*nshgtot;
		writedatablock_( &irstin, "ybar",
				(void*)(yglobal), &nitems, "double", iotype );
		closefile_( &irstin, "write" );

	}


	free(qglobal);
	if(RequestedAcceleration){
		free(aglobal);
	}
	if(RequestedBoundaryFluxes){
		free(fglobal);
	}
	if(RequestedDisplacements){
		free(dglobal);
		free(dglobal_ref);
	}
	if(RequestedDistances){
		free(distglobal);
	}
	if(RequestedWSS){
		free(wglobal);
	}
	if(RequestedYbar){
		free(yglobal);
	}
	if (aleOn){
#if DEBUG_ALE == 1
		free(resglobal);
#endif	
		free(relativeVelocityGlobal);
		free(updatedMeshCoordinatesGlobal);
	}
	free(xglobal);
	free(ien);
	return 0;
}
