#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cstring> 
#include <map>
#include <utility>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

#include "cvSolverIO.h"
#include "common_c.h"

#ifdef intel 
#include <direct.h>
#define chdir _chdir
#define stat _stat
#define sonfath_ SONFATH

#else
#include <unistd.h>
#include <strings.h>
#endif

//using namespace PHSOLVER;

char keyphrase[100];
char filename[255];
char _directory_name[256];

extern "C" {
void
METIS_PartGraphKway(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*,
		int*);
void sonfath_(int*, int*, int*, int*, int*, int*, int*, int*);
}

void
Gather_Headers(int* fileDescriptor, std::vector<std::string>& headers);

using namespace std;

typedef vector<int> block;

struct lessKey {
	bool operator()(block b1, block b2) const {
		return (100 * b1[1] + 10 * b1[5] < 100 * b2[1] + 10 * b2[5]);
	}
};

typedef struct lessKey lessKey;
typedef map<block, vector<vector<int> >, lessKey> Block_Map;

struct CommuTask {
	int tag;
	int type;
	int other_pid;
	vector<int> shpFn;
};

typedef struct CommuTask CommuTask;

// since ParallelData[x] is a sorted data structure, the lowest pid 
// adjacency always  becomes the master.

int getMaster(int x, map<int, map<int, int> > &ParallelData) {
	map<int, int>::iterator iter = ParallelData[x].begin();
	return (*iter).first;

}
void generate_keyphrase(char* target, const char* prefix, block& tpblock) {

	strcpy(target, prefix); /* interior or boundary */

	switch (tpblock[1]) {
	case 1:
		strcat(target, "linear ");
		break;
	case 2:
		strcat(target, "quadratic ");
		break;
	case 3:
		strcat(target, "cubic ");
		break;
	case 4:
		strcat(target, "quartic ");
		break;
	}

	switch (tpblock[5]) {
	case 1:
		strcat(target, "tetrahedron ");
		break;
	case 2:
		strcat(target, "hexahedron ");
		break;
	case 3:
		if (!strcmp(prefix, "boundary "))
			strcat(target, "wedge triface ");
		else
			strcat(target, "wedge ");
		break;
	case 4:
		strcat(target, "wedge quadface ");
		break;
	case 5:
		if (!strcmp(prefix, "boundary "))
			strcat(target, "pyramid quadface ");
		else
			strcat(target, "pyramid ");
		break;
	case 6:
		strcat(target, "pyramid triface ");
		break;
	}
}

void Partition_Problem(int numProcs) {
	int rank = 0;
	int stepno;
	int igeombc; /* file handle for geombc */
	int irestart; /* file handle for restart */
	int iarray[10];
	int ione = 1;
	int itwo = 2;
	int ithree = 3;
	int iseven = 7;
	int ieight = 8;
	int isize;
	int ierr;
	int nitems;

	char iformat[80];
	char oformat[80];
	int SONFATH_VAR = turbvar.sonfathvar;

	map<int, map<int, int> > ParallelData;

	// ParallelData is the most important data structure of this code. It maintains
	// the global to local mapping of shapefuntion numbers.
	//
	// ParallelData[ global_number][ processor_id ] = processor_local_no
	//
	// all shape function numbers, global and local are stored in fortran indexing

	strcpy(iformat, outpar.iotype);
	strcpy(oformat, outpar.iotype);

	/* Lets us first load the dual and partition it */

	_directory_name[0] = '\0';

	sprintf(_directory_name, "%d-procs-case/", numProcs);

#ifdef intel
	struct _stat buf;
	char probdir[256];
	probdir[0]='\0';
	sprintf(probdir,"%d-procs-case",numProcs);
	int result = _stat( probdir, &buf );

	//    if ( (result == 0 ) &&  ( 01 & ( buf.st_mode >> _S_IFDIR ) ) ) 
	if ( (result == 0 ) && ( buf.st_mode & _S_IFDIR ) )
#else
	struct stat buf;
	if (!stat(_directory_name, &buf) && S_ISDIR(buf.st_mode))
#endif
			{
		cout << "Directory " << _directory_name << " already exists " << endl;
		cout << "using the existing inputfiles" << endl;
		//if ( !chdir( _directory_name ) ) {/* cd successful */
		//	cout <<"changing to the problem directory " << _directory_name << endl;
		return;
		//}else { 
		//	cout << "Please check the permissions for the problem";
		//	cout << " directory " << _directory_name << endl;
		//	exit(0);
		//}

	} else {
		cout << _directory_name << " does not exist or is unusable" << endl;
		cout << "creating a new one" << endl;
#ifdef intel
		_unlink( _directory_name );
#else
		unlink(_directory_name);
#endif
	}
	// else if the case does not already exist

#ifdef intel
	if ( _mkdir( _directory_name ) )
#else
	if (mkdir(_directory_name, 00755))
#endif
			{
		cerr << "cannot create directory " << _directory_name << endl;
		exit(0);
	}
	filename[0] = '\0';
	sprintf(filename, "%snumpe.in", _directory_name);
	ofstream npe(filename);
	npe << numProcs << endl;
	npe.close();

	// copy the bct.dat file if one exists

	char systemcmd[2048];
	systemcmd[0] = '\0';

	sprintf(systemcmd, "cp *.dat %s", _directory_name);

	system(systemcmd);

	if (numProcs < 2) {
		sprintf(systemcmd, "cp *.1 %s", _directory_name);
		system(systemcmd);
		return;
	}

	openfile_("geombc.dat.1", "read", &igeombc);
	vector < string > headers;
	Gather_Headers(&igeombc, headers);

	int numel, nIfaces, spebc;
	readheader_(&igeombc, "keyword xadj", (void*) iarray, &itwo, "integer",
			iformat);
	numel = iarray[0];
	spebc = iarray[1];
	isize = numel + 1;
	int* xadj = new int[isize];
	readdatablock_(&igeombc, "keyword xadj", (void*) xadj, &isize, "integer",
			iformat);

	readheader_(&igeombc, "keyword adjncy", (void*) iarray, &ione, "integer",
			iformat);
	nIfaces = iarray[0];
	isize = 2 * nIfaces;
	int* adjncy = new int[isize];
	readdatablock_(&igeombc, "keyword adjncy", (void*) adjncy, &isize,
			"integer", iformat);

	isize = numel;
	int* vwgt = new int[numel];
	readheader_(&igeombc, "keyword vwgt", (void*) iarray, &ione, "integer",
			iformat);
	readdatablock_(&igeombc, "keyword vwgt", (void*) vwgt, &isize, "integer",
			iformat);

	int wgtflag = 2;
	int numflag = 0;
	int edgecut;
	int options[10];
	options[0] = 0;
	int* epart = new int[numel];
	METIS_PartGraphKway(&numel, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag,
			&numProcs, options, &edgecut, epart);

	/* putting SPEBC elements together on the processor 0 */
	if (spebc > 0) {
		int swap;
		for (int i = 0; i < numel; i++)
			if (vwgt[i] > 1) {
				swap = epart[i];
				break;
			}
		for (int i = 0; i < numel; i++) {
			if (epart[i] == swap)
				epart[i] = 0;
			else if (epart[i] == 0)
				epart[i] = swap;
			if (vwgt[i] == 0)
				epart[i] = 0;
		}
	}

	delete[] xadj;
	delete[] adjncy;
	delete[] vwgt;

	/* now we have the epart array which is all the partition Info we need */

	int ndof = 5;
	ifstream numstart("numstart.dat");
	numstart >> stepno;
	numstart.close();

	readheader_(&igeombc, "number of global modes", (void*) iarray, &ione,
			"integer", iformat);
	int nshgTot = iarray[0];
	int nshg = iarray[0];

	readheader_(&igeombc, "number of variables", (void*) iarray, &ione,
			"integer", iformat);

	ndof = iarray[0];

	readheader_(&igeombc, "maximum number of element nodes", (void*) iarray,
			&ione, "integer", iformat);
	int nenmax = iarray[0];

	readheader_(&igeombc, "number of interior tpblocks", (void*) iarray, &ione,
			"integer", iformat);
	int nblock = iarray[0];

	readheader_(&igeombc, "number of nodes in the mesh", (void*) iarray, &ione,
			"integer", iformat);
	int numnp = iarray[0];

	readheader_(&igeombc, "number of edges in the mesh", (void*) iarray, &ione,
			"integer", iformat);
	int numedges = iarray[0];

	vector<Block_Map> ien(numProcs);

	int* vcount = new int[numProcs];
	int* ecount = new int[numProcs];
	int* fcount = new int[numProcs];

	for (int e = 0; e < numProcs; e++)
		vcount[e] = ecount[e] = fcount[e] = 1;

	for (int b = 0; b < nblock; b++) {

		readheader_(&igeombc, "connectivity interior?", (void*) iarray, &iseven,
				"integer", iformat);
		block CurrentBlock;
		for (int w = 1; w < iseven; w++)
			CurrentBlock.push_back(iarray[w]);
		isize = iarray[0] * iarray[3];
		int* ient = new int[isize];

		readdatablock_(&igeombc, "connectivity interior?", (void*) ient, &isize,
				"integer", iformat);

		readheader_(&igeombc, "ien to sms?", (void*) iarray, &ione, "integer",
				iformat);

		isize = iarray[0];
		int* ient_sms = new int[isize];
		readdatablock_(&igeombc, "ien to sms?", (void*) ient_sms, &isize,
				"integer", iformat);

		for (int c = 0; c < iarray[0]; c++) {
			vector<int> element;
			for (int d = 0; d < iarray[3]; d++) {
				element.push_back(ient[d * iarray[0] + c]);
			}
			int gid = ient_sms[c];
			int pid = epart[gid];

			ien[pid][CurrentBlock].push_back(element);

			for (vector<int>::iterator iter = element.begin();
					iter != element.end(); iter++) {
				if (!ParallelData[abs(*iter)][pid])
					ParallelData[abs(*iter)][pid] =
							(abs(*iter) > numnp) ?
									((abs(*iter) > numnp + numedges) ?
											fcount[pid]++ : ecount[pid]++) :
									vcount[pid]++;
			}
			element.clear();
		}
		delete[] ient;
		delete[] ient_sms;
	}

	/* At this point the values stored in vcount ecount and fcount are 1 more 
	 * than  the number of vertices edges and faces respectively , sp we need
	 * to decrement  all of them by one so that we can directly added them to
	 * the lower order counts  in the next  step.
	 */

	for (int e = 0; e < numProcs; e++) {
		--vcount[e];
		--ecount[e];
		--fcount[e];
	}

	/* our parallel data structure has indepenent numbering for vertices, faces
	 * and edges, now we need to correct it and add the correct offsets , we
	 * cannot write out ien before we make this correction
	 */

	for (map<int, map<int, int> >::iterator iter = ParallelData.begin();
			iter != ParallelData.end(); iter++) {
		if ((*iter).first > numnp) { // it is not a vertex
			if ((*iter).first > (numnp + numedges)) { // it is not a edge,
				// must be a face
				for (map<int, int>::iterator modeIter = (*iter).second.begin();
						modeIter != (*iter).second.end(); modeIter++)
					(*modeIter).second = (*modeIter).second
							+ vcount[(*modeIter).first]
							+ ecount[(*modeIter).first];
			} else { // is an edge

				for (map<int, int>::iterator modeIter = (*iter).second.begin();
						modeIter != (*iter).second.end(); modeIter++)
					(*modeIter).second = (*modeIter).second
							+ vcount[(*modeIter).first];
			}
		}
	}

	delete[] vcount;
	delete[] ecount;
	delete[] fcount;

	/* after correction all local mode numbers are accurate ( hopefully ) and
	 * are in fortran style ( starting 1 ) */

	/* at this point, ien still has global numbering */
	/* we now use the parallel data structure to correct ien */

	for (int p = 0; p < numProcs; p++) {
		for (Block_Map::iterator iter = ien[p].begin(); iter != ien[p].end();
				iter++) {
			for (vector<vector<int> >::iterator iter2 =
					((*iter).second).begin(); iter2 != ((*iter).second).end();
					iter2++) {
				for (vector<int>::iterator iter3 = (*iter2).begin();
						iter3 != (*iter2).end(); iter3++) {
					*iter3 =
							(*iter3 > 0) ?
									ParallelData[*iter3][p] :
									-1 * ParallelData[abs(*iter3)][p];
				}
			}
		}
	}

	/* now we can write out ien and delete that data structure. */
	/* ien is stored as ien[pid][block][iel][nshl] */
	/* we need to write this block by block for each pid */

	int magic_number = 362436;
	int* mptr = &magic_number;
	int fgeom;
	int frest;

	for (int p = 0; p < numProcs; p++) {
		int numel = 0;

		bzero_old((void*) filename, 255);
		sprintf(filename, "%sgeombc.dat.%d", _directory_name, p + 1);
		openfile_(filename, "write", &fgeom);

		for (vector<string>::iterator iter = headers.begin();
				iter != headers.end(); iter++) {
			writestring_(&fgeom, (*iter).c_str());
			writestring_(&fgeom, "\n");
		}

#if defined ( DEBUG )
		sprintf(filename,"%sascii_out.%d",_directory_name,p+1);
		ofstream fascii( filename );
		fascii <<"Interior Element Modal Connectivity:"<<endl;
		fascii <<"------------------------------------"<<endl;
#endif

		isize = 1;
		nitems = 1;
		iarray[0] = 1;
		writeheader_(&fgeom, "byteorder magic number ", (void*) iarray, &nitems,
				&isize, "integer", oformat);

		nitems = 1;
		writedatablock_(&fgeom, "byteorder magic number ", (void*) mptr,
				&nitems, "integer", oformat);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of interior tpblocks : < 0 > %d \n",
				(int) ien[p].size());
		writestring_(&fgeom, filename);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of global modes : < 0 > %d \n", nshgTot);
		writestring_(&fgeom, filename);

		for (Block_Map::iterator bIter = ien[p].begin(); bIter != ien[p].end();
				bIter++) {
			block CurrentBlock = (*bIter).first;
			vector < vector<int> > blockIEN = (*bIter).second;
			numel += blockIEN.size();
			isize = blockIEN.size() * CurrentBlock[2];
			iarray[0] = blockIEN.size();

#if defined ( DEBUG )
			fascii << endl;
			fascii << "********************"<< endl;
#endif

			int xct = 1;
			for (block::iterator ibtr = CurrentBlock.begin();
					ibtr != CurrentBlock.end(); ibtr++) {
				iarray[xct++] = *ibtr;

#if defined ( DEBUG )
				fascii << *ibtr <<" ";
#endif
			}

			generate_keyphrase(keyphrase, "connectivity interior ",
					CurrentBlock);
			writeheader_(&fgeom, keyphrase, (void*) iarray, &xct, &isize,
					"integer", oformat);

#if defined ( DEBUG )
			fascii << endl;
			fascii << "********************"<< endl;
#endif
			/* now we need to create the fortran style array for this block of
			 * ien, so need to invert it */

			int* ient = new int[blockIEN.size() * CurrentBlock[2]];
			for (int h = 0; h < blockIEN.size(); h++) {
				for (int y = 0; y < CurrentBlock[2]; y++) {
					ient[y * blockIEN.size() + h] = blockIEN[h][y];
#if defined ( DEBUG )
					fascii << ient[ y*blockIEN.size() + h ] <<" ";
#endif
				}
#if defined ( DEBUG )
				fascii << endl;
#endif
			}
			nitems = blockIEN.size() * CurrentBlock[2];
			writedatablock_(&fgeom, keyphrase, (void*) (ient), &nitems,
					"integer", oformat);

			delete[] ient;
			blockIEN.clear();
		} // loop over blocks with in a processor
		ien[p].clear();

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of interior elements : < 0 > %d \n", numel);
		writestring_(&fgeom, filename);

		closefile_(&fgeom, "write");
#if defined ( DEBUG )
		fascii.close();
#endif
	} // loop over procs
	  //ien.clear();

	/* co-ords */

	readheader_(&igeombc, "co-ordinates", (void*) iarray, &itwo, "double",
			iformat);
	isize = iarray[0] * iarray[1];
	double* xloc = new double[isize];
	readdatablock_(&igeombc, "co-ordinates", (void*) xloc, &isize, "double",
			iformat);
	vector < map<int, vector<double> > > Xpart(numProcs);

	for (int x = 1; x < iarray[0] + 1; x++) {
		for (map<int, int>::iterator iter = ParallelData[x].begin();
				iter != ParallelData[x].end(); iter++) {
			for (int s = 0; s < 3; s++) {
				Xpart[(*iter).first][(*iter).second].push_back(
						xloc[x - 1 + s * iarray[0]]);
				//cout << (*iter).first << " " << (*iter).second << " " << xloc[x - 1 + s * iarray[0]] << endl;
			}
		}
	}
	delete[] xloc;

	// if we need to calculate sonfath later we will need to have vnumnp[i]
	int* vnumnp;
	if (SONFATH_VAR > 0)
		vnumnp = new int[numProcs];

	for (int p = 0; p < numProcs; p++) {

		bzero_old((void*) filename, 255);
		sprintf(filename, "%sgeombc.dat.%d", _directory_name, p + 1);
		openfile_(filename, "append", &fgeom);

		if (SONFATH_VAR > 0)
			vnumnp[p] = Xpart[p].size();

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of nodes : < 0 > %d \n",
				(int) Xpart[p].size());
		writestring_(&fgeom, filename);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of processors : < 0 > %d \n", numProcs);
		writestring_(&fgeom, filename);

		double* xf = new double[Xpart[p].size() * 3];
		for (map<int, vector<double> >::iterator mIter = Xpart[p].begin();
				mIter != Xpart[p].end(); mIter++)
			for (int f = 0; f < 3; f++)
				xf[f * Xpart[p].size() + (*mIter).first - 1] =
						(*mIter).second[f];

		isize = Xpart[p].size() * 3;
		nitems = 2;
		iarray[0] = Xpart[p].size();
		iarray[1] = 3;
		writeheader_(&fgeom, "co-ordinates", (void*) iarray, &nitems, &isize,
				"double", oformat);

		nitems = Xpart[p].size() * 3;
		writedatablock_(&fgeom, "co-ordinates", (void*) (xf), &nitems, "double",
				oformat);

#if defined ( DEBUG )
		sprintf( filename, "%sascii_out.%d",_directory_name, p+1 );
		FILE* cfascii=fopen(filename,"a");
		fprintf(cfascii,"\n");
		fprintf(cfascii, "Nodal coordinates:\n");
		fprintf(cfascii, "------------------\n");
		int numnp_tmp = Xpart[p].size();
		for( int c=0; c< numnp_tmp; c++ )
		fprintf( cfascii,"%f %f %f \n",
				xf[ c ], xf[ c + numnp_tmp ],xf[c+2*numnp_tmp] );
		fclose( cfascii );
#endif          

		delete[] xf;
		Xpart[p].clear();
		closefile_(&fgeom, "append");
	}
	Xpart.clear();

	/* let us take care of Essential BCs.*/
	// BCs are in the indirect numbering of nBC so only nBC needs to be sorted 
	// according to nshg
	// BC abd iBC are numbered as numpbc which us given by nBC.
	int numEBC = ndof + 7;
	readheader_(&igeombc, "bc mapping array", (void*) iarray, &ione, "integer",
			iformat);
	nshg = iarray[0];
	int* nBC = new int[nshg];
	readdatablock_(&igeombc, "bc mapping array", (void*) nBC, &nshg, "integer",
			iformat);

	readheader_(&igeombc, "bc codes array", (void*) iarray, &ione, "integer",
			iformat);
	int numpbc = iarray[0];
	int* iBC = new int[numpbc];
	readdatablock_(&igeombc, "bc codes array", (void*) iBC, &numpbc, "integer",
			iformat);

	readheader_(&igeombc, "boundary condition array", (void*) iarray, &ione,
			"double", iformat);
	double* BC = new double[numpbc * numEBC];
	readdatablock_(&igeombc, "boundary condition array", (void*) BC, &iarray[0],
			"double", iformat);

	int* inumpbc = new int[numProcs];
	for (int c = 0; c < numProcs; c++)
		inumpbc[c] = 0;
	vector < map<int, int> > nBCpart(numProcs);
	vector < vector<int> > iBCpart(numProcs);
	vector < vector<double> > BCpart(numProcs);

	for (int x = 1; x < nshg + 1; x++) {
		if (nBC[x - 1] > 0) {
			for (map<int, int>::iterator iter = ParallelData[x].begin();
					iter != ParallelData[x].end(); iter++) {
				nBCpart[(*iter).first][(*iter).second] =
						inumpbc[(*iter).first]++;
				iBCpart[(*iter).first].push_back(iBC[nBC[x - 1] - 1]);
				for (int g = 0; g < numEBC; g++)
					BCpart[(*iter).first].push_back(
							BC[nBC[x - 1] - 1 + g * numpbc]);
			}
		} else {
			for (map<int, int>::iterator iter = ParallelData[x].begin();
					iter != ParallelData[x].end(); iter++)
				nBCpart[(*iter).first][(*iter).second] = -1;
		}
	}

	delete[] nBC;
	delete[] iBC;
	delete[] BC;



	// we will calculate a maxnshg here for later use with sonfath.

	int maxnshg = 0;

	for (int p = 0; p < numProcs; p++) {

		bzero_old((void*) filename, 255);
		sprintf(filename, "%sgeombc.dat.%d", _directory_name, p + 1);
		openfile_(filename, "append", &fgeom);

#if defined ( DEBUG )
		sprintf( filename, "%sascii_out.%d",_directory_name, p+1 );
		FILE* cfascii = fopen( filename, "a" );
		fprintf(cfascii,"\nBoundary Condition Mapping Array (nBC):\n" );
		fprintf(cfascii,"---------------------------------------\n" );
#endif
		if ((SONFATH_VAR > 0) && (maxnshg < nBCpart[p].size()))
			maxnshg = nBCpart[p].size();
		bzero_old((void*) filename, 255);
		sprintf(filename, "maximum number of element nodes : < 0 > %d\n",
				nenmax);
		writestring_(&fgeom, filename);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of modes : < 0 > %d\n",
				(int) nBCpart[p].size());
		writestring_(&fgeom, filename);

		bzero_old((void*) filename, 255);
		sprintf(filename,
				"number of shapefunctions saved on processor : < 0 > %d\n",
				(int) nBCpart[p].size());
		writestring_(&fgeom, filename);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of nodes with Dirichlet BCs : < 0 > %d\n",
				inumpbc[p]);
		writestring_(&fgeom, filename);

		isize = nBCpart[p].size();
		nitems = 1;
		iarray[0] = nBCpart[p].size();
		writeheader_(&fgeom, "bc mapping array", (void*) iarray, &nitems,
				&isize, "integer", oformat);

		int* inBC = new int[nBCpart[p].size()];
		int count = 0;
		for (map<int, int>::iterator iter = nBCpart[p].begin();
				iter != nBCpart[p].end(); iter++) {
			inBC[count++] = ((*iter).second + 1);

#if defined ( DEBUG )
			fprintf( cfascii, "%d\n", inBC[ count - 1 ] );
#endif

		}
		writedatablock_(&fgeom, "bc mapping array", (void*) (inBC), &count,
				"integer", oformat);
		delete[] inBC;
		nBCpart[p].clear();

#if defined ( DEBUG )
		fprintf(cfascii,"\nBoundary Condition Codes Array (iBC):\n" );
		fprintf(cfascii,"-------------------------------------\n" );
#endif
		isize = inumpbc[p];
		nitems = 1;
		iarray[0] = inumpbc[p];

		writeheader_(&fgeom, "bc codes array ", (void*) iarray, &nitems, &isize,
				"integer", oformat);

		int* iiBC = new int[inumpbc[p]];
		if (inumpbc[p] != iBCpart[p].size())
			cerr << "Error 1" << __LINE__ << endl;
		for (int i = 0; i < inumpbc[p]; i++) {
			iiBC[i] = iBCpart[p][i];
#if defined ( DEBUG )
			fprintf( cfascii, "%d\n", iiBC[ i ] );
#endif
		}
		iBCpart[p].clear();
		nitems = inumpbc[p];
		writedatablock_(&fgeom, "bc codes array ", (void*) (iiBC), &nitems,
				"integer", oformat);
		delete[] iiBC;

#if defined ( DEBUG )
		fprintf(cfascii,"\nBoundary Condition Array (BC):\n" );
		fprintf(cfascii,"--------------------------------\n" );
#endif
		isize = inumpbc[p] * numEBC;
		nitems = 1;
		iarray[0] = inumpbc[p] * numEBC;
		writeheader_(&fgeom, "boundary condition array ", (void*) iarray,
				&nitems, &isize, "double", oformat);

		double* BCf = new double[inumpbc[p] * numEBC];
		for (int a = 0; a < inumpbc[p]; a++) {
			for (int b = 0; b < numEBC; b++) {
				BCf[b * inumpbc[p] + a] = BCpart[p][a * numEBC + b];

#if defined ( DEBUG )
				fprintf( cfascii,"%.3f ",BCpart[p][a*numEBC+b]);
#endif

			}

#if defined ( DEBUG )
			fprintf( cfascii, "\n" );
#endif

		}
		BCpart[p].clear();
		nitems = inumpbc[p] * numEBC;
		writedatablock_(&fgeom, "boundary condition array ", (void*) (BCf),
				&nitems, "double", oformat);

		delete[] BCf;
		closefile_(&fgeom, "append");

#if defined ( DEBUG )
		fclose( cfascii );
#endif
	}

	iBCpart.clear();
	BCpart.clear();

	/* done writing ebcs */

	/* ienb and BCB and iBCB */

	int numNBC = ndof + 1;
	int nPropshere = nomodule.nProps;
	int nPropsTS = 2;
	int nPropsEst = 1;

	readheader_(&igeombc, "number of boundary tpblocks", (void*) iarray, &ione,
			"integer", iformat);
	nblock = iarray[0];

	//ienb -- Boundary Element Nodal Connectivity.
	//ibcb -- Natural Boundary Condition Codes
	//bcb  -- Natural Boundary Condition Values
	//swb  -- Vessel Wall Properties array
	//twb  -- Tissue Support Properties array

	vector < map<block, vector<int>, lessKey> > iBCBpart(numProcs);
	vector < map<block, vector<double>, lessKey> > BCBpart(numProcs);
	vector < map<block, vector<double>, lessKey> > SWBpart(numProcs);
	vector < map<block, vector<double>, lessKey> > TWBpart(numProcs);

	for (int b = 0; b < nblock; b++) {

		readheader_(&igeombc, "connectivity boundary?", (void*) iarray, &ieight,
				"integer", iformat);

		block CurrentBlock;
		for (int w = 1; w < ieight; w++)
			CurrentBlock.push_back(iarray[w]);

		isize = iarray[0] * iarray[3];
		int* ient = new int[isize];
		readdatablock_(&igeombc, "connectivity boundary?", (void*) ient, &isize,
				"integer", iformat);

		readheader_(&igeombc, "ienb to sms?", (void*) iarray, &ione, "integer",
				iformat);

		isize = iarray[0];
		int* ient_sms = new int[isize];
		readdatablock_(&igeombc, "ienb to sms?", (void*) ient_sms, &isize,
				"integer", iformat);

		readheader_(&igeombc, "nbc codes?", (void*) iarray, &ieight, "integer",
				iformat);

		isize = iarray[0] * 2;
		int* iBCB = new int[isize];
		readdatablock_(&igeombc, "nbc codes?", (void*) iBCB, &isize, "integer",
				iformat);

		readheader_(&igeombc, "nbc values?", (void*) iarray, &ieight, "double",
				iformat);

		isize = iarray[0] * numNBC;
		double* BCB = new double[isize];
		readdatablock_(&igeombc, "nbc values?", (void*) BCB, &isize, "double",
				iformat);

		double* SWB = NULL;
		double* TWB = NULL;

		if (nomodule.ideformwall != 0 && nomodule.iUseSWB != 0) {
			readheader_(&igeombc, "SWB array?", (void*) iarray, &itwo, "double",
					iformat);

			isize = iarray[0] * nPropshere;
			SWB = new double[isize];
			readdatablock_(&igeombc, "SWB array?", (void*) SWB, &isize,
					"double", iformat);
		}

		if (nomodule.ideformwall != 0 && nomodule.iUseTWB != 0) {
			// added by Nan 06/26/09
			readheader_(&igeombc, "TWB array?", (void*) iarray, &itwo, "double",
					iformat);

			isize = iarray[0] * nPropsTS;
			TWB = new double[isize];
			readdatablock_(&igeombc, "TWB array?", (void*) TWB, &isize,
					"double", iformat);

		}

		for (int c = 0; c < iarray[0]; c++) {
			vector<int> element;
			for (int d = 0; d < iarray[3]; d++)
				element.push_back(ient[d * iarray[0] + c]);

			int gid = ient_sms[c];
			int pid = epart[gid];

			for (vector<int>::iterator iter = element.begin();
					iter != element.end(); iter++) {
				*iter = (*iter > 0) ?
						ParallelData[*iter][pid] :
						-1 * ParallelData[abs(*iter)][pid];

			}
			ien[pid][CurrentBlock].push_back(element);
			element.clear();

			iBCBpart[pid][CurrentBlock].push_back(iBCB[c]);
			iBCBpart[pid][CurrentBlock].push_back(iBCB[c + iarray[0]]);

			for (int o = 0; o < numNBC; o++)
				BCBpart[pid][CurrentBlock].push_back(BCB[c + iarray[0] * o]);

			if (nomodule.ideformwall != 0 && nomodule.iUseSWB != 0) {
				for (int o = 0; o < nPropshere; o++)
					SWBpart[pid][CurrentBlock].push_back(
							SWB[c + iarray[0] * o]);
			}

			if (nomodule.ideformwall != 0 && nomodule.iUseTWB != 0) {
				for (int o = 0; o < nPropsTS; o++)
					TWBpart[pid][CurrentBlock].push_back(
							TWB[c + iarray[0] * o]);
			}

		}

		delete[] ient;
		delete[] ient_sms;
		delete[] iBCB;
		delete[] BCB;

		if (nomodule.ideformwall != 0 && nomodule.iUseSWB != 0) {
			delete[] SWB;
		}

		if (nomodule.ideformwall != 0 && nomodule.iUseTWB != 0) {
			delete[] TWB;
		}

	}

	delete[] epart;

	/* format and write out */

	for (int p = 0; p < numProcs; p++) {

		bzero_old((void*) filename, 255);
		sprintf(filename, "%sgeombc.dat.%d", _directory_name, p + 1);
		openfile_(filename, "append", &fgeom);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of boundary tpblocks : < 0 > %d\n",
				(int) ien[p].size());
		writestring_(&fgeom, filename);

#if defined ( DEBUG )
		sprintf(filename,"%sascii_out.%d",_directory_name,p+1);
		ofstream fascii( filename, ios_base::app );
		fascii << endl;
		fascii <<"Boundary Element Modal Connectivity:"<<endl;
		fascii <<"------------------------------------"<<endl;
#endif
		int numelb = 0;
		for (Block_Map::iterator bIter = ien[p].begin(); bIter != ien[p].end();
				bIter++) {
			block CurrentBlock = (*bIter).first;
			vector < vector<int> > blockIEN = (*bIter).second;
			numelb += blockIEN.size();

			isize = blockIEN.size() * CurrentBlock[2];
			iarray[0] = blockIEN.size();

#if defined ( DEBUG )
			fascii << endl;
			fascii <<"************************"<< endl;
#endif
			int xct = 1;
			for (block::iterator ibtr = CurrentBlock.begin();
					ibtr != CurrentBlock.end(); ibtr++) {
				iarray[xct++] = *ibtr;
#if defined ( DEBUG )
				fascii << *ibtr <<" ";
#endif
			}
			generate_keyphrase(keyphrase, "connectivity boundary ",
					CurrentBlock);
			writeheader_(&fgeom, keyphrase, (void*) iarray, &xct, &isize,
					"integer", oformat);

#if defined ( DEBUG )
			fascii << endl;
			fascii <<"************************"<< endl;
#endif

			/* now we need to create the fortran style array for this block of
			 * ien, so need to invert it */

			int* ient = new int[blockIEN.size() * CurrentBlock[2]];
			for (int h = 0; h < blockIEN.size(); h++) {
				for (int y = 0; y < CurrentBlock[2]; y++) {
					ient[y * blockIEN.size() + h] = blockIEN[h][y];
#if defined ( DEBUG )
					fascii << ient[ y*blockIEN.size() + h ] <<" ";
#endif
				}
#if defined ( DEBUG )
				fascii << endl;
#endif
			}

			nitems = blockIEN.size() * CurrentBlock[2];
			writedatablock_(&fgeom, keyphrase, (void*) (ient), &nitems,
					"integer", oformat);

			delete[] ient;

			/* need to also write IBCB and BCB now */
#if defined ( DEBUG )
			fascii << "------------------------" << endl;
			fascii << "iBCB for the above block" << endl;
			fascii << "------------------------" << endl;
#endif

			int* iBCBf = new int[blockIEN.size() * 2];
			for (int u = 0; u < blockIEN.size(); u++) {
				iBCBf[u] = iBCBpart[p][CurrentBlock][u * 2];
				iBCBf[u + blockIEN.size()] =
						iBCBpart[p][CurrentBlock][u * 2 + 1];

#if defined ( DEBUG )
				fascii << iBCBf[u] <<"  "<<iBCBf[ u + blockIEN.size() ]<< endl;
#endif

			}

			iBCBpart[p][CurrentBlock].clear();

			isize = blockIEN.size() * 2;
			nitems = xct;
			generate_keyphrase(keyphrase, "nbc codes ", CurrentBlock);
			writeheader_(&fgeom, keyphrase, (void*) iarray, &xct, &isize,
					"integer", oformat);

			nitems = blockIEN.size() * 2;
			writedatablock_(&fgeom, keyphrase, (void*) (iBCBf), &nitems,
					"integer", oformat);

			delete[] iBCBf;

#if defined ( DEBUG )
			fascii.precision( 8 );
			fascii << "------------------------" << endl;
			fascii << "BCB for the above block" << endl;
			fascii << "------------------------" << endl;
			fascii.precision( 8 );
#endif
			double* BCBf = new double[blockIEN.size() * numNBC];
			for (int u = 0; u < blockIEN.size(); u++) {
				for (int v = 0; v < numNBC; v++) {
					BCBf[v * blockIEN.size() + u] = BCBpart[p][CurrentBlock][u
							* numNBC + v];

#if defined ( DEBUG )
					fascii << BCBf[v* blockIEN.size()+u] << " ";
#endif

				}

#if defined ( DEBUG )
				fascii << endl;
#endif

			}

			BCBpart[p][CurrentBlock].clear();

			isize = blockIEN.size() * numNBC;
			nitems = xct;
			generate_keyphrase(keyphrase, "nbc values ", CurrentBlock);
			writeheader_(&fgeom, keyphrase, (void*) iarray, &xct, &isize,
					"double", oformat);

			nitems = blockIEN.size() * numNBC;
			writedatablock_(&fgeom, keyphrase, (void*) (BCBf), &nitems,
					"double", oformat);
			delete[] BCBf;

			if (nomodule.ideformwall != 0) {
				/* need to also write SWB if deformable wall exists now */

				if (nomodule.iUseSWB != 0) {

#if defined ( DEBUG )
					fascii.precision( 8 );
					fascii << "------------------------" << endl;
					fascii << "SWB for the above block" << endl;
					fascii << "------------------------" << endl;
					fascii.precision( 8 );
#endif
					double* SWBf = new double[blockIEN.size() * nPropshere];
					for (int u = 0; u < blockIEN.size(); u++) {
						for (int v = 0; v < nPropshere; v++) {
							SWBf[v * blockIEN.size() + u] =
									SWBpart[p][CurrentBlock][u * nPropshere + v];

#if defined ( DEBUG )   
							fascii << SWBf[v* blockIEN.size()+u] << " ";
#endif

						}

#if defined ( DEBUG )
						fascii << endl;
#endif
					}

					SWBpart[p][CurrentBlock].clear();

					isize = blockIEN.size() * nPropshere;
					xct = 2;
					iarray[1] = nPropshere;
					generate_keyphrase(keyphrase, "SWB array ", CurrentBlock);
					writeheader_(&fgeom, keyphrase, (void*) iarray, &xct,
							&isize, "double", oformat);

					nitems = blockIEN.size() * nPropshere;
					writedatablock_(&fgeom, keyphrase, (void*) (SWBf), &nitems,
							"double", oformat);
					delete[] SWBf;
				}

				if (nomodule.iUseTWB != 0) {

#if defined ( DEBUG )
					fascii.precision( 8 );
					fascii << "------------------------" << endl;
					fascii << "TWB for the above block" << endl;
					fascii << "------------------------" << endl;
					fascii.precision( 8 );
#endif
					double* TWBf = new double[blockIEN.size() * nPropsTS];
					for (int u = 0; u < blockIEN.size(); u++) {
						for (int v = 0; v < nPropsTS; v++) {
							TWBf[v * blockIEN.size() + u] =
									TWBpart[p][CurrentBlock][u * nPropsTS + v];

#if defined ( DEBUG )   
							fascii << TWBf[v* blockIEN.size()+u] << " ";
#endif

						}

#if defined ( DEBUG )
						fascii << endl;
#endif
					}

					TWBpart[p][CurrentBlock].clear();

					isize = blockIEN.size() * nPropsTS;
					xct = 2;
					iarray[1] = nPropsTS;
					generate_keyphrase(keyphrase, "TWB array ", CurrentBlock);
					writeheader_(&fgeom, keyphrase, (void*) iarray, &xct,
							&isize, "double", oformat);

					nitems = blockIEN.size() * nPropsTS;
					writedatablock_(&fgeom, keyphrase, (void*) (TWBf), &nitems,
							"double", oformat);
					delete[] TWBf;
				}

			} // end of ideformwall check

			blockIEN.clear();
		} // loop over blocks with in a processor
		ien[p].clear();

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of boundary elements : < 0 > %d \n", numelb);
		writestring_(&fgeom, filename);

		closefile_(&fgeom, "append");

#if defined ( DEBUG )
		fascii.close();
#endif
	} // loop over procs
	ien.clear();

	/* now starts commu and periodicity */
	/* We need to generate ilwork and ncorp and the iper array for each
	 * processor . To achieve this we first read in the complete iper array and
	 * then go through "ParallelData" iterating all the global shapefunctions
	 * and then generate a consistant communication trace for all
	 * for on processor periodiciy we just dump an array */

	readheader_(&igeombc, "periodic masters array?", (void*) iarray, &ione,
			"integer", iformat);
	int* periodic = new int[nshg];
	readdatablock_(&igeombc, "periodic masters array?", (void*) periodic, &nshg,
			"integer", iformat);

	vector < map<int, int> > PeriodicPart(numProcs);

	CommuTask*** rtask;
	CommuTask*** stask;

	rtask = new CommuTask**[numProcs];
	stask = new CommuTask**[numProcs];
	for (int b = 0; b < numProcs; b++) {
		rtask[b] = new CommuTask*[numProcs];
		stask[b] = new CommuTask*[numProcs];
		for (int c = 0; c < numProcs; c++) {
			rtask[b][c] = NULL;
			stask[b][c] = NULL;
		}
	}

	int newtag = 1;
	bool isPERIODIC = false;
	for (int x = 1; x < nshg + 1; x++) {

		// first we need to fill a 0 in all the iper entries of x and its images 
		// since this will be the case even if the node has no periodicity at all 

		isPERIODIC = false;

		for (map<int, int>::iterator iter = ParallelData[x].begin();
				iter != ParallelData[x].end(); iter++)
			PeriodicPart[(*iter).first][(*iter).second] = 0;

		int master_pid = getMaster(x,ParallelData); // assuming no periodicity
		int Global_Master = ParallelData[x][master_pid];

		if (periodic[x - 1] > 0) {

			// if there is periodicity, then the real master node and processor are 
			// not that of the current mode but of the master mode,  and need 
			// updating...also now the image which shares this master_pid should get 
			// an iper and not a commu.

			master_pid = getMaster(periodic[x - 1],ParallelData);

			// master_pid is the processor on which the master is real ( present ) 

			Global_Master = ParallelData[periodic[x - 1]][master_pid];

			// Global_Master is the processor local shape function number of the master
			// mode on the master processor

			isPERIODIC = true;
		}
		// periodicity has been taken care of 
		// now we do the simple commu to the master_pid 

		for (map<int, int>::iterator iter = ParallelData[x].begin();
				iter != ParallelData[x].end(); iter++)
			if ((*iter).first != master_pid) {

				int sender = (*iter).first;
				int receiver = master_pid;

				// if we ever add a send task from A to B then we also immediately add a 
				// receive task from B to A so its enuf if we check one of them 
				if (!stask[sender][receiver]) {
					stask[sender][receiver] = new CommuTask;
					rtask[receiver][sender] = new CommuTask;
					stask[sender][receiver]->tag = newtag;
					rtask[receiver][sender]->tag = newtag++;
					stask[sender][receiver]->type = 0;
					rtask[receiver][sender]->type = 1;
					stask[sender][receiver]->other_pid = receiver;
					rtask[receiver][sender]->other_pid = sender;
				}

				stask[sender][receiver]->shpFn.push_back((*iter).second);
				rtask[receiver][sender]->shpFn.push_back(Global_Master);

			} else {
				// we  have an image of the original slave which is on the master
				// processor of the master node for this image we need to fill iper
				// instead of commu. ( ofcourse only for periodic modes )

				if (isPERIODIC) {
					PeriodicPart[(*iter).first][(*iter).second] = Global_Master;
				}
			}
	}

	delete[] periodic;

	// now to generate ILwork and Periodicity

	vector<int> ilwork;

#if defined ( DEBUG)
	sprintf(filename,"%silwork.info",_directory_name );
	ofstream ilwf( filename );
#endif

	for (int a = 0; a < numProcs; a++) { //outer loop over procs
		ilwork.push_back(0);

#if defined ( DEBUG )
		ilwf << a+1 <<" receives  ";
#endif

		for (int b = 0; b < numProcs; b++) { // inner loop 1 over procs
			if (rtask[a][b]) {
				ilwork[0]++;
				ilwork.push_back(rtask[a][b]->tag);
				ilwork.push_back(rtask[a][b]->type);
				ilwork.push_back(rtask[a][b]->other_pid + 1);
				ilwork.push_back(rtask[a][b]->shpFn.size());

#if defined ( DEBUG )
				ilwf << " from " << rtask[a][b]->other_pid+1
				<< " tag " << rtask[a][b]->tag
				<< " numSeg " << rtask[a][b]->shpFn.size() << endl;
#endif

				for (vector<int>::iterator iter = rtask[a][b]->shpFn.begin();
						iter != rtask[a][b]->shpFn.end(); iter++) {
					ilwork.push_back(*iter);
					ilwork.push_back(1);

#if defined ( DEBUG )
					ilwf << *iter << " ,";
#endif

				}

#if defined ( DEBUG )
				ilwf << endl;
#endif

				delete rtask[a][b];

			}

		} // end inner loop 1 over procs

#if defined ( DEBUG )

		ilwf << endl;

		ilwf << a+1 << " sends  ";
#endif

		for (int b = 0; b < numProcs; b++) { // inner loop 2 over procs
			if (stask[a][b]) {
				ilwork[0]++;
				ilwork.push_back(stask[a][b]->tag);
				ilwork.push_back(stask[a][b]->type);
				ilwork.push_back(stask[a][b]->other_pid + 1);
				ilwork.push_back(stask[a][b]->shpFn.size());

#if defined ( DEBUG )
				ilwf << " to " << stask[a][b]->other_pid+1
				<< " tag " << stask[a][b]->tag
				<< " numSeg " << stask[a][b]->shpFn.size() << endl;
#endif

				for (vector<int>::iterator iter = stask[a][b]->shpFn.begin();
						iter != stask[a][b]->shpFn.end(); iter++) {
					ilwork.push_back(*iter);
					ilwork.push_back(1);

#if defined ( DEBUG )
					ilwf << *iter << " ,";
#endif

				}

#if defined ( DEBUG )
				ilwf << endl;
#endif

				delete stask[a][b];
			}
		} // end inner 2loop over procs

#if defined ( DEBUG )
		ilwf << endl;
		ilwf << "Size of ILwork for "<< a+1 << " : " << ilwork.size() << endl;
		if ( ilwork.size() > 0 ) {
			for( vector<int>::iterator iter = ilwork.begin();
					iter != ilwork.end();
					iter++ )
			ilwf << *iter <<" ";

			ilwf << endl;
			ilwf << endl;
		}
#endif

		delete[] rtask[a];
		delete[] stask[a];

		bzero_old((void*) filename, 255);
		sprintf(filename, "%sgeombc.dat.%d", _directory_name, a + 1);
		openfile_(filename, "append", &fgeom);

		bzero_old((void*) filename, 255);
		sprintf(filename, "size of ilwork array : < 0 > %d\n",
				(int) ilwork.size());
		writestring_(&fgeom, filename);

		isize = ilwork.size();
		nitems = 1;
		iarray[0] = ilwork.size();
		writeheader_(&fgeom, "ilwork ", (void*) iarray, &nitems, &isize,
				"integer", oformat);

		int* filwork = new int[ilwork.size()];
		for (int r = 0; r < ilwork.size(); r++)
			filwork[r] = ilwork[r];
		nitems = ilwork.size();
		writedatablock_(&fgeom, "ilwork ", (void*) (filwork), &nitems,
				"integer", oformat);
		ilwork.clear();
		delete[] filwork;

		isize = PeriodicPart[a].size();
		nitems = 1;
		iarray[0] = PeriodicPart[a].size();
		writeheader_(&fgeom, "periodic masters array ", (void*) iarray, &nitems,
				&isize, "integer", oformat);

		int* fper = new int[PeriodicPart[a].size()];
		int count = 0;
		for (map<int, int>::iterator iter = PeriodicPart[a].begin();
				iter != PeriodicPart[a].end(); iter++)
			fper[count++] = (*iter).second;

#if defined ( DEBUG )

		sprintf( filename, "%sascii_out.%d",_directory_name, a+1 );
		ofstream fascii( filename, ios_base::app );
		fascii <<"\n------------------------"<<endl;
		fascii <<"Periodic Partners Array:"<<endl;
		fascii <<"------------------------"<<endl;
		for( int shg=0; shg < count; shg++ )
		fascii << fper[ shg ] << endl;
		fascii.close();
#endif

		nitems = count;
		writedatablock_(&fgeom, "periodic masters array ", (void*) (fper),
				&nitems, "integer", oformat);

		delete[] fper;
		PeriodicPart[a].clear();
		closefile_(&fgeom, "append");

	} // end outer loop over procs

#if defined ( DEBUG )
	ilwf.close();
#endif

	PeriodicPart.clear();

	delete[] rtask;
	delete[] stask;

	// write ncorp
	// generating ncorp 
	// ncorp is our map between the "partition local" and the global (sms) 
	// numbering of modes  ncorp[ proc ][ local number ] = sms_number

	vector < map<int, int> > ncorp(numProcs);
	for (int x = 1; x < nshg + 1; x++) {
		for (map<int, int>::iterator iter = ParallelData[x].begin();
				iter != ParallelData[x].end(); iter++)
			ncorp[(*iter).first][(*iter).second] = x;
	}

	// generate an index of the local nodes with no duplication of nodes across
	// processors -- this is used when interfacing with external routines that require
	// passing a solution array that is uniquely partitioned among processors

	vector < map<int, int> > NodesUniquePart(numProcs);

	for (int x = 1; x < nshg + 1; x++) {
		map<int, int>::iterator pIter = ParallelData[x].begin();

		NodesUniquePart[(*pIter).first][(*pIter).second] = (*pIter).second;

		//cout << (*pIter).first << " " << (*pIter).second << endl;
	}

	//
	// generate the partitioned structure for the simple observation functions
	//
	int* linobs_soln;
	vector < map<int, vector<int> > > linobs_soln_Part(numProcs);

	readheader_(&igeombc, "observation function solution?", (void*) iarray, &itwo, "integer",
			iformat);

	isize = nshg * ndof;

	linobs_soln = new int[isize];

	readdatablock_(&igeombc, "observation function solution?", (void*) linobs_soln, &isize,
			"integer", iformat);

	for (int x = 1; x < nshg + 1; x++) {
		for (map<int, int>::iterator pIter = ParallelData[x].begin();
				pIter != ParallelData[x].end(); pIter++) {
			for (int v = 0; v < ndof; v++)
				linobs_soln_Part[(*pIter).first][(*pIter).second].push_back(
						linobs_soln[v * nshg + x - 1]);
		}
	}

	delete[] linobs_soln;

	int* linobs_acc;
	vector < map<int, vector<int> > > linobs_acc_Part(numProcs);

	readheader_(&igeombc, "observation function time derivative of solution?", (void*) iarray,
			&itwo, "integer", iformat);

	isize = nshg * ndof;

	linobs_acc = new int[isize];

	readdatablock_(&igeombc, "observation function time derivative of solution?",
			(void*) linobs_acc, &isize, "integer", iformat);

	for (int x = 1; x < nshg + 1; x++) {
		for (map<int, int>::iterator pIter = ParallelData[x].begin();
				pIter != ParallelData[x].end(); pIter++) {
			for (int v = 0; v < ndof; v++)
				linobs_acc_Part[(*pIter).first][(*pIter).second].push_back(
						linobs_acc[v * nshg + x - 1]);
		}
	}

	delete[] linobs_acc;

	vector < map<int, vector<int> > > linobs_disp_Part(numProcs);
	vector < map<int, vector<int> > > obs_dist_Part(numProcs);

	if (nomodule.ideformwall != 0) {
		int* linobs_disp;

		readheader_(&igeombc, "observation function displacement?",
				(void*) iarray, &itwo, "integer", iformat);

		isize = nshg * 3;

		linobs_disp = new int[isize];

		readdatablock_(&igeombc,
				"observation function displacement?",
				(void*) linobs_disp, &isize, "integer", iformat);

		for (int x = 1; x < nshg + 1; x++) {
			for (map<int, int>::iterator pIter = ParallelData[x].begin();
					pIter != ParallelData[x].end(); pIter++) {
				for (int v = 0; v < 3; v++) {
					linobs_disp_Part[(*pIter).first][(*pIter).second].push_back(
							linobs_disp[v * nshg + x - 1]);

					//cout << (*pIter).first << " " << (*pIter).second << " " << linobs_disp[v * nshg + x - 1] << endl;
				}
			}
		}

		delete[] linobs_disp;

		int* obs_dist;

		readheader_(&igeombc, "observation function distance?",
				(void*) iarray, &itwo, "integer", iformat);

		isize = nshg;

		obs_dist = new int[isize];

		readdatablock_(&igeombc,
				"observation function distance?",
				(void*) obs_dist, &isize, "integer", iformat);

		for (int x = 1; x < nshg + 1; x++) {
			for (map<int, int>::iterator pIter = ParallelData[x].begin();
					pIter != ParallelData[x].end(); pIter++) {

				obs_dist_Part[(*pIter).first][(*pIter).second].push_back(
						obs_dist[x - 1]);

			}
		}

		delete[] obs_dist;


	}


	// SONFATH stuff

	if (SONFATH_VAR > 0) {
		int* ncvec = new int[maxnshg * numProcs];
		for (int x = 0; x < numProcs * maxnshg; x++)
			ncvec[x] = 0;

		for (int p = 0; p < numProcs; p++)
			for (map<int, int>::iterator niter = ncorp[p].begin();
					niter != ncorp[p].end(); niter++)
				ncvec[numProcs * (*niter).first + p] = (*niter).second;

		int* nsons = new int[SONFATH_VAR];
		int* ifath = new int[maxnshg * numProcs];

		sonfath_(ncvec, vnumnp, ifath, nsons, &SONFATH_VAR, &numProcs, &nshgTot,
				&maxnshg);

		delete[] vnumnp;
		delete[] ncvec;

		// now we have ifath and nsons for the whole mesh, need to split it 
		// and write out.

		int* ifath1d;

		for (int a = 0; a < numProcs; a++) {

			bzero_old((void*) filename, 255);
			sprintf(filename, "%sgeombc.dat.%d", _directory_name, a + 1);
			openfile_(filename, "append", &fgeom);

			bzero_old((void*) filename, 255);
			sprintf(filename, "number of father-nodes : < 0 > %d\n",
					SONFATH_VAR);
			writestring_(&fgeom, filename);

			isize = SONFATH_VAR;
			nitems = 1;
			iarray[0] = SONFATH_VAR;
			writeheader_(&fgeom, "number of son-nodes for each father ",
					(void*) iarray, &nitems, &isize, "integer", oformat);
			nitems = SONFATH_VAR;
			writedatablock_(&fgeom, "number of son-nodes for each father ",
					(void*) (nsons), &nitems, "integer", oformat);

			ifath1d = new int[ncorp[a].size()];
			for (int j = 0; j < ncorp[a].size(); j++)
				ifath1d[j] = ifath[j * numProcs + a];

			isize = ncorp[a].size();
			nitems = 1;
			iarray[0] = ncorp[a].size();
			writeheader_(&fgeom, "keyword ifath ", (void*) iarray, &nitems,
					&isize, "integer", oformat);
			nitems = ncorp[a].size();
			writedatablock_(&fgeom, "keyword ifath ", (void*) (ifath1d),
					&nitems, "integer", oformat);

#if defined ( DEBUG )
			sprintf( filename, "%sascii_out.%d",_directory_name, a+1 );
			ofstream fascii( filename, ios_base::app );
			fascii <<"\n-------------------------------"<< endl;
			fascii <<"IFATH array "<<endl;
			fascii <<"-------------------------------"<< endl;
			for( int j=0; j < ncorp[a].size(); j++ )
			fascii << ifath1d[j] << endl;
			fascii <<"-------------------------------"<< endl;
			fascii.close();
#endif

			closefile_(&fgeom, "append");
			delete[] ifath1d;
		}

		delete[] nsons;
		delete[] ifath;
	}

	for (int a = 0; a < numProcs; a++) {

		bzero_old((void*) filename, 255);
		sprintf(filename, "%sgeombc.dat.%d", _directory_name, a + 1);
		openfile_(filename, "append", &fgeom);

		isize = ncorp[a].size();
		nitems = 1;
		iarray[0] = ncorp[a].size();
		writeheader_(&fgeom, " mode number map from partition to global",
				(void*) iarray, &nitems, &isize, "integer", oformat);

		int* fncorp = new int[ncorp[a].size()];
		int count = 0;
		for (map<int, int>::iterator iter = ncorp[a].begin();
				iter != ncorp[a].end(); iter++)
			fncorp[count++] = (*iter).second;

		nitems = count;
		writedatablock_(&fgeom, " mode number map from partition to global",
				(void*) (fncorp), &nitems, "integer", oformat);

		delete[] fncorp;

		isize = NodesUniquePart[a].size();
		nitems = 1;
		iarray[0] = NodesUniquePart[a].size();
		writeheader_(&fgeom, "local index of unique nodes",
				(void*) iarray, &nitems, &isize, "integer", oformat);

		int* fNodesUnique = new int[NodesUniquePart[a].size()];
		count = 0;
		for (map<int, int>::iterator iter = NodesUniquePart[a].begin();
				iter != NodesUniquePart[a].end(); iter++)
			fNodesUnique[count++] = (*iter).second;

		nitems = count;
		writedatablock_(&fgeom, "local index of unique nodes",
				(void*) (fNodesUnique), &nitems, "integer", oformat);

		delete[] fNodesUnique;


		int nshgLocal = linobs_soln_Part[a].size();
		int nsd = 3;

		//
		isize = nshgLocal * ndof;
		nitems = 2;
		iarray[0] = nshgLocal;
		iarray[1] = ndof;

		writeheader_(&fgeom, "observation function solution", (void*) iarray, &nitems, &isize,
				"integer", oformat);

		int* flinobs_soln = new int[isize];
		for (int w = 0; w < ndof; w++)
			for (int y = 1; y < nshgLocal + 1; y++) {
				flinobs_soln[w * nshgLocal + (y - 1)] =
						linobs_soln_Part[a][y][w];
			}

		nitems = isize;
		writedatablock_(&fgeom, "observation function solution", (void*) (flinobs_soln),
				&nitems, "integer", oformat);
		delete[] flinobs_soln;

        //
		isize = nshgLocal * ndof;
		nitems = 2;
		iarray[0] = nshgLocal;
		iarray[1] = ndof;

		writeheader_(&fgeom, "observation function time derivative of solution", (void*) iarray,
				&nitems, &isize, "integer", oformat);

		int* flinobs_acc = new int[isize];
		for (int w = 0; w < ndof; w++)
			for (int y = 1; y < nshgLocal + 1; y++) {
				flinobs_acc[w * nshgLocal + (y - 1)] =
						linobs_acc_Part[a][y][w];
			}

		nitems = isize;
		writedatablock_(&fgeom, "observation function time derivative of solution",
				(void*) (flinobs_acc), &nitems, "integer", oformat);
		delete[] flinobs_acc;

		//

		if (nomodule.ideformwall != 0) {

			isize = nshgLocal * nsd;
			nitems = 2;
			iarray[0] = nshgLocal;
			iarray[1] = nsd;

			writeheader_(&fgeom, "observation function displacement",
					(void*) iarray, &nitems, &isize, "integer", oformat);

			int* flinobs_disp = new int[isize];
			for (int w = 0; w < nsd; w++)
				for (int y = 1; y < nshgLocal + 1; y++) {
					flinobs_disp[w * nshgLocal + (y - 1)] = linobs_disp_Part[a][y][w];
				}

			nitems = isize;
			writedatablock_(&fgeom,
					"observation function displacement",
					(void*) (flinobs_disp), &nitems, "integer", oformat);
			delete[] flinobs_disp;

			isize = nshgLocal;
			nitems = 1;
			iarray[0] = nshgLocal;

			writeheader_(&fgeom, "observation function distance",
					(void*) iarray, &nitems, &isize, "integer", oformat);

			int* fobs_dist = new int[isize];

			for (int y = 1; y < nshgLocal + 1; y++) {
				fobs_dist[(y - 1)] = obs_dist_Part[a][y][0];
			}

			nitems = isize;
			writedatablock_(&fgeom,
					"observation function distance",
					(void*) (fobs_dist), &nitems, "integer", oformat);
			delete[] fobs_dist;

		}


		closefile_(&fgeom, "append");
		ncorp[a].clear();
		NodesUniquePart[a].clear();
		linobs_soln_Part[a].clear();
		linobs_acc_Part[a].clear();
		linobs_disp_Part[a].clear();
	}

	//write restart

	sprintf(filename, "restart.%d.1", stepno);
	openfile_(filename, "read", &irestart);
	readheader_(&irestart, "solution?", (void*) iarray, &ithree, "double",
			iformat);

	nshg = iarray[0];
	ndof = iarray[1];
	isize = nshg * ndof;

	double* solution = new double[isize];

	readdatablock_(&irestart, "solution?", (void*) solution, &isize, "double",
			iformat);

	// now we need to create the partitioned data structure.
	// keeping in mind that unlike BCs this has a direct corelation 
	// to each nshg so we need to use a 
	// sorted data structure and we choose a std::map for simplicity

	vector < map<int, vector<double> > > solPart(numProcs);

	for (int x = 1; x < nshg + 1; x++) {
		for (map<int, int>::iterator pIter = ParallelData[x].begin();
				pIter != ParallelData[x].end(); pIter++) {
			for (int v = 0; v < ndof; v++)
				solPart[(*pIter).first][(*pIter).second].push_back(
						solution[v * nshg + x - 1]);
		}
	}

	delete[] solution;

	readheader_(&irestart, "time derivative of solution?", (void*) iarray, &ithree, "double",
				iformat);

	int iexistAccel = 0;

	double* fAccel;
	vector < map<int, vector<double> > > accelPart(numProcs);

	if (iarray[0] != 0) {

		nshg = iarray[0];
		ndof = iarray[1];
		isize = nshg * ndof;

		double* accel = new double[isize];

		readdatablock_(&irestart, "time derivative of solution?", (void*) accel, &isize, "double",
				iformat);

		for (int x = 1; x < nshg + 1; x++) {
			for (map<int, int>::iterator pIter = ParallelData[x].begin();
					pIter != ParallelData[x].end(); pIter++) {
				for (int v = 0; v < ndof; v++)
					accelPart[(*pIter).first][(*pIter).second].push_back(
							accel[v * nshg + x - 1]);
			}
		}

		delete [] accel;

		iexistAccel = 1;
	}

	int nsd = 3;
	double* displacement;
	double* fDisplacement;
	int nshgLocalDisp;
	vector < map<int, vector<double> > > dispPart(numProcs);

	// check flag to see if we expect displacements to exist
	if (nomodule.ideformwall != 0) {

		//same deal for the displacements

		readheader_(&irestart, "displacement?", (void*) iarray, &ithree,
				"double", iformat);

		nshg = iarray[0];
		nsd = iarray[1];
		isize = nshg * nsd;

		displacement = new double[isize];

		readdatablock_(&irestart, "displacement?", (void*) displacement, &isize,
				"double", iformat);

		for (int x = 1; x < nshg + 1; x++) {
			for (map<int, int>::iterator pIter = ParallelData[x].begin();
					pIter != ParallelData[x].end(); pIter++) {
				for (int v = 0; v < nsd; v++)
					dispPart[(*pIter).first][(*pIter).second].push_back(
							displacement[v * nshg + x - 1]);
			}
		}

		delete[] displacement;

	}

	// now we have the partitioned data structure we need to write it out.

	for (int a = 0; a < numProcs; a++) {

		int nshgLocal = solPart[a].size();
		double* fSolution = new double[nshgLocal * ndof];

		if (iexistAccel) {
			fAccel = new double[nshgLocal * ndof];
		}

		if (nomodule.ideformwall != 0) {
			nshgLocalDisp = dispPart[a].size();
			fDisplacement = new double[nshgLocalDisp * nsd];
		}

#if defined ( DEBUG )
		sprintf( filename,"%srestart.asc.%d",_directory_name, a+1 );
		FILE* rascii = fopen( filename,"w");
		fprintf( rascii, "nshgLocal: %d\n", nshgLocal);
		fprintf( rascii, "numVars: %d\n", ndof );
		fprintf( rascii, "Step Number: %d\n", stepno);
		fprintf( rascii, "---------------------------------------------\n");
		for( int y=1; y < nshgLocal+1; y++ ) {
			fprintf( rascii, "%d :", y );
			for( int w=0; w< ndof; w++ ) {
				fprintf(rascii," %f ", solPart[a][y][w] );
			}
			fprintf( rascii, "\n");
		}
		if (nomodule.ideformwall != 0) {
			fprintf( rascii, "nshgLocalDisp: %d\n", nshgLocalDisp);
			fprintf( rascii, "numVars: %d\n", nsd );
			fprintf( rascii, "Step Number: %d\n", stepno);
			fprintf( rascii, "---------------------------------------------\n");
			for( int y=1; y < nshgLocalDisp+1; y++ ) {
				fprintf( rascii, "%d :", y );
				for( int w=0; w< nsd; w++ ) {
					fprintf(rascii," %f ", dispPart[a][y][w] );
				}
				fprintf( rascii, "\n");
			}
		}
		fclose( rascii );
#endif
		for (int w = 0; w < ndof; w++)
			for (int y = 1; y < nshgLocal + 1; y++) {
				fSolution[w * nshgLocal + (y - 1)] = solPart[a][y][w];
			}

		solPart[a].clear();

		if (iexistAccel) {
			for (int w = 0; w < ndof; w++)
				for (int y = 1; y < nshgLocal + 1; y++) {
					fAccel[w * nshgLocal + (y - 1)] = accelPart[a][y][w];
				}
			accelPart[a].clear();
		}

		if (nomodule.ideformwall != 0) {
			for (int w = 0; w < nsd; w++)
				for (int y = 1; y < nshgLocalDisp + 1; y++) {
					fDisplacement[w * nshgLocalDisp + (y - 1)] =
							dispPart[a][y][w];
				}

			dispPart[a].clear();
		}

		bzero_old((void*) filename, 255);
		sprintf(filename, "%srestart.%d.%d", _directory_name, stepno, a + 1);
		openfile_(filename, "write", &frest);

		for (vector<string>::iterator iter = headers.begin();
				iter != headers.end(); iter++) {
			writestring_(&frest, (*iter).c_str());
			writestring_(&frest, "\n");
		}
		headers.clear();

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of modes : < 0 > %d\n", nshgLocal);
		writestring_(&frest, filename);

		bzero_old((void*) filename, 255);
		sprintf(filename, "number of variables : < 0 > %d\n", ndof);
		writestring_(&frest, filename);

		isize = 1;
		nitems = 1;
		iarray[0] = 1;
		writeheader_(&frest, "byteorder magic number ", (void*) iarray, &nitems,
				&isize, "integer", oformat);

		nitems = 1;
		writedatablock_(&frest, "byteorder magic number ", (void*) mptr,
				&nitems, "integer", oformat);

		isize = nshgLocal * ndof;
		nitems = 3;
		iarray[0] = nshgLocal;
		iarray[1] = ndof;
		iarray[2] = stepno;
		writeheader_(&frest, "solution ", (void*) iarray, &nitems, &isize,
				"double", oformat);

		nitems = nshgLocal * ndof;
		writedatablock_(&frest, "solution ", (void*) (fSolution), &nitems,
				"double", oformat);

		if (iexistAccel) {
			isize = nshgLocal * ndof;
			nitems = 3;
			iarray[0] = nshgLocal;
			iarray[1] = ndof;
			iarray[2] = stepno;
			writeheader_(&frest, "time derivative of solution ", (void*) iarray, &nitems, &isize,
					"double", oformat);

			nitems = nshgLocal * ndof;
			writedatablock_(&frest, "time derivative of solution ", (void*) (fAccel), &nitems,
					"double", oformat);
		}

		if (nomodule.ideformwall != 0) {
			isize = nshgLocalDisp * nsd;
			nitems = 3;
			iarray[0] = nshgLocalDisp;
			iarray[1] = nsd;
			iarray[2] = stepno;
			writeheader_(&frest, "displacement ", (void*) iarray, &nitems,
					&isize, "double", oformat);

			nitems = nshgLocalDisp * nsd;
			writedatablock_(&frest, "displacement ", (void*) (fDisplacement),
					&nitems, "double", oformat);
		}

		closefile_(&frest, "write");

		delete[] fAccel;
		delete[] fSolution;

		if (nomodule.ideformwall != 0) {
			delete[] fDisplacement;
		}

	}

	filename[0] = '\0';
	sprintf(filename, "%snumstart.dat", _directory_name);
	ofstream nstart(filename);
	nstart << stepno;
	nstart.close();

	closefile_(&igeombc, "read");
	closefile_(&irestart, "read");
	//if ( !chdir( _directory_name ) ) /* cd successful */
	//	cout <<" changing to the problem directory " << _directory_name << endl;
	ParallelData.clear();
}
