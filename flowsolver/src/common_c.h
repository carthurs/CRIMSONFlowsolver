// Routine contains the structures for reading the user input through
// input_fform. The default values for all these variables are defined in
// input.config.
//
// Input variables that have been previously declared in common.h have to be
// re-declared here, in a consistant structure. 

#ifndef _COMMON_C_H_
#define _COMMON_C_H_

#include "mpi.h"

#define MAXBLK   50000
#define MAXSURF  199
#define MAXREGIONS 255
#define MAXTS   100
#define MAXTOP   5
#define MAXQPT   125
#define MAXSH    32
#define NSD      3
//#define machin   'RS/6000'
//#define machfl   4
//#define zero   0.0000000000000000000000000000000d0
//#define pt125   0.1250000000000000000000000000000d0
//#define pt25   0.2500000000000000000000000000000d0
//#define pt33   0.3333333333333333333333333333333d0
//#define pt39   0.3968502629920498686879264098181d0
//#define pt5   0.5000000000000000000000000000000d0
//#define pt57   0.5773502691896257645091487805020d0
//#define pt66   0.6666666666666666666666666666667d0
//#define pt75   0.7500000000000000000000000000000d0
//#define one   1.0000000000000000000000000000000d0
//#define sqt2   1.4142135623730950488016887242097d0
//#define onept5   1.5000000000000000000000000000000d0
//#define two   2.0000000000000000000000000000000d0
//#define three   3.0000000000000000000000000000000d0
//#define four   4.0000000000000000000000000000000d0
//#define five   5.0000000000000000000000000000000d0
//#define pi   3.1415926535897932384626433832795d0

//namespace PHSOLVER {

#ifdef __cplusplus
   extern "C" {
#endif

      extern struct { 
         int master;
         int numpe;
         int myrank;
      } workfc ;

      extern struct { 
         int maxfront;
         int nlwork;
      } fronts ;

      extern struct { 
         int numper;
         int nshgt;
         int nshg0;
      } newdim ;

      extern struct { 
         double birth;
         double death;
         double comtim;
      } timer4 ;

      extern struct { 
         double ttim[100];
      } extrat ;

//      extern struct {
//         double zoutsponge, radsponge, zinsponge, grthosponge, grthisponge;
//         double betamax;
//         int spongecontinuity, spongemomentum1, spongemomentum2;
//         int spongeenergy, spongemomentum3;
//      } spongevar ;

      extern struct {
//         double eles;
//         double ylimit[9][3]; /* 9 = 5 + 4 = puvwT + 4Scalars */
//         double rmutarget;
//         double pzero;
//         double wtavei;
//         double dtavei;
//         double dke;
//         double fwr1;
//         double flump;
         int ierrcalc;
         int ihessian;
//         int itwmod;
//         int ngaussf;
//         int idim;
//         int nlist;
//         int nintf[MAXTOP];
//         int sonfathvar;
      } turbvar ;

      extern struct {
//         int irans, iles, isubmod;
//         int ifproj;
//         int i2filt;
//         int modlstats;
//         int idis;
//         int nohomog;
         int ierrsmooth;
         /* wonder if we should include nintf(MAXTOP) and MAXTOP since its
         in common.h */

         /*      int itwmod; */
         /*      double rtavei; */
         /*      int ierrcalc; */
      } turbvari ;

//      extern struct {
//         int irscale;
//         int intpres;
//         double plandist;
//         double thetag;
//         double ds;
//         double tolerence;
//         double radcyl;
//         double rbltin;
//         double rvscal;
//      } spebcvr ;

      extern struct {
         double scdiff[5];
         double tdecay;
         int nsclr, isclr,nsolt, nosource;
         int consrv_sclr_conv_vel;
      } sclrs;

      extern struct { 
         double flxID[MAXSURF+1][10] ;
         double Force[3];
         double HFlux;
         double flxIDsclr[MAXSURF][4];
         int nsrflist[MAXSURF+1];
         int isrfIM;
      } aerfrc ;

//      extern struct {
//         double a[100000];
//      } astore ;

      extern struct { 
         int numnp;
         int numel;
         int numelb;
         int numpbc;
         int nen;
         int nwnp;
         int nfaces;
         int numflx;
         int ndof;
         int iALE;
         int icoord;
         int navier;
         int irs;
         int iexec;
         int necho;
         int ichem;
         int iRK;
         int nedof;
         int nshg;
         int nnz;
         int istop;
         int nflow;
         int nnz_tot;
         int idtn;
         int nshguniq;
         int icurrentblk;
      } conpar ;

      extern struct {
         double epsilon_ls;
         double epsilon_lsd;
         double dtlset;
         int iLSet;
         int ivconstraint;
         int iExpLSSclr1;
         int iExpLSSclr2;
      } levlset;

      extern struct { 
         int nshape;
         int nshapeb;
         int maxshb;
         int nshl;
         int nshlb;
         int nfath;
         int ntopsh;
         int nsonmax;
      } shpdat ;

      extern struct { 
         int mshp;
         int mshgl;
         int mwght;
         int mshpb;
         int mshglb;
         int mwghtb;
         int mmut;
         int mrhot;
         int mxst;
      } datpnt ;

      extern struct { 
         int lelCat;
         int lcsyst;
         int iorder;
         int nenb;
         int nelblk;
         int nelblb;
         int ndofl;
         int nsymdl;
         int nenl;
         int nfacel;
         int nenbl;
         int intind;
         int mattyp;
      } elmpar ;

      extern struct { 
         double E3nsd;
         int I3nsd;
         int nsymdf;
         int ndofBC;
         int ndiBCB;
         int ndBCB;
         int Jactyp;
         int jump;
         int ires;
         int iprec;
         int iprev;
         int ibound;
         int idiff;
         int lhs;
         int itau;
         int ipord;
         int ipred;
         int lstres;
         int iepstm;
         double dtsfct;
         double taucfct;
         int ibksiz;
         int iabc;
         int isurf;
         int idflx;
         double Bo;
         int EntropyPressure;
      } genpar ;

      extern struct { 
         double epstol[8];  /* 1+ max number of scalars  (beginning of the
                            end of time sequences) */
         double Delt[MAXTS];
         double CFLfl[MAXTS];
         double CFLsl[MAXTS];
         int nstep[MAXTS];
         int niter[MAXTS];
         int impl[MAXTS];
         double rhoinf[MAXTS];
         int LHSupd[6];
         int loctim[MAXTS];
         double deltol[2][MAXTS];
         int memLSFlag;
      } inpdat ;

      extern struct { 
         int iin;
         int igeom;
         int ipar;
         int ibndc;
         int imat;
         int iecho;
         int iout;
         int ichmou;
         int irstin;
         int irstou;
         int ihist;
         int iflux;
         int ierror;
         int itable;
         int iforce;
         int igraph;
         int itime;
      } mio ;

      extern struct { 
         double fin;
         double fgeom;
         double fpar;
         double fbndc;
         double fmat;
         double fecho;
         double frstin;
         double frstou;
         double fhist;
         double ferror;
         double ftable;
         double fforce;
         double fgraph;
         double ftime;
      } mioname ;

      extern struct { 
         double eGMRES;
         int lGMRES;
         int iKs;
         int ntotGM;
      } itrpar ;

      extern struct { 
         int mHBrg;
         int meBrg;
         int myBrg;
         int mRcos;
         int mRsin;
      } itrpnt ;

      extern struct { 
         double datmat[MAXTS][7][3];
         int matflg[MAXTS][6];
         int nummat;
         int mexist;
      } matdat ;

      extern struct {
         double pr, Planck, Stephan, Nh, Rh, Rgas;
         double gamma, gamma1, s0;
         //, const, xN2, xO2;
         //double yN2,    yO2,    Msh[5], cpsh[5],s0sh[5],h0sh[5];
         //double Rs[5],  cps[5], cvs[5], h0s[5], Trot[5],sigs[5];
         //double Tvib[5],g0s[5], dofs[5],ithm;
      } mmatpar ;

      extern struct { 
         double ro;
         double vel;
         double temper;
         double press;
         double entrop;
         int ntout;
         int ioform;
         int iowflux;
         int iofieldv;
         char iotype[80];
         int ioybar;
         /*  int iostats; */
         /*      int ipresref; */
      } outpar ;

      extern struct { 
         int mbeg;
         int mend;
         int mprec;
      } point ;

      extern struct { 
         double epsM;
         int iabres;
      } precis ;

      extern struct { 
         int npro;
      } propar ;

      extern struct { 
         double resfrt;
      } resdat ;

      extern struct { 
         int imap;
         int ivart;
         int iDC;
         int iPcond;
         int Kspace;
         int nGMRES;
         int iconvflow;
         int iconvsclr;
         int idcsclr[2];
      } solpar ;

      extern struct { 
         double time;
         double CFLfld;
         double CFLsld;
         double Dtgl;
         double Dtmax;
         double alpha;
         double etol;
         int lstep;
         int ifunc;
         int itseq;
         int istep;
         int iter;
         int nitr;
         double almi;
         double alfi;
         double gami;
         double flmpl;
         double flmpr;
         double dtol[2];
         int iCFLworst;
      } timdat ;

      extern struct { 
         int LCtime;
         int ntseq;
      } timpar ;

      extern struct { 
         int numeqns[100];
         int minIters;
         int maxIters;
         int iprjFlag;
         int nPrjs;
         int ipresPrjFlag;
         int nPresPrjs;
         double prestol;
         double statsflow[6];
         double statssclr[6];
         int iverbose;
      } incomp ;

      extern struct { 
         double ccode[13];
      } mtimer1 ;

      extern struct { 
         double flops;
         double gbytes;
         double sbytes;
         int iclock;
         int icd;
         int icode;
         int icode2;
         int icode3;
      } mtimer2 ;

      extern struct { 
         double cpu[11];
         double cpu0[11];
         int nacess[11];
      } timer3 ;

      extern struct { 
         double title;
         int ititle;
      } title ;

      extern struct {
         int intg[MAXTS][2];
      }intdat;

      extern struct {
         double bcttimescale;    
         double ValueListResist[MAXSURF+1];
         double rhovw;
         double thicknessvw;
         double evw;
         double rnuvw;
         double rshearconstantvw;
         double betai;
         double ValueListWallE[MAXREGIONS+1];
         double ValueListWallh[MAXREGIONS+1];
         double tissSuppStiffCoeff; 
         double tissSuppDampCoeff;
         double tissSuppRingStiffCoeff;
         double tissSuppRingDampCoeff;
         double stateFilterCoeff;
         double rescontrol;
         double ResCriteria;
         double heartparam[16];
         double stabflux_coeff;
         int icardio;
         int itvn;
         int ipvsq;
         int incp;
         int numINCPSrfs;
         int nsrflistINCP[MAXSURF+1];
         int incpfile;
         int numResistSrfs;
         int nsrflistResist[MAXSURF+1];
         int numImpSrfs;
         int nsrflistImp[MAXSURF+1];
         int impfile;
         int numRCRSrfs;
         int nsrflistRCR[MAXSURF+1];
         int ircrfile;
         int numTRCRSrfs;
         int nsrflistTRCR[MAXSURF+1];
         int itrcrfile;
         int numCORSrfs;
         int nsrflistCOR[MAXSURF+1];
         int icorfile;
         int numVisFluxSrfs;
         int nsrflistVisFlux[MAXSURF+1];
         int numCalcSrfs;
         int nsrflistCalc[MAXSURF+1];
         int numDirCalcSrfs;
         int nsrflistDirCalc[MAXSURF+1];
         int Lagrange; 
         int numLagrangeSrfs;
         int nsrflistLagrange[MAXSURF+1];
         int iLagfile;
         int MinNumIter;
         int ideformwall;
         int iwallmassfactor;
         int iwallstiffactor;
         int nProps;
         int iUseSWB;
         int iUseSWBthickonly;
         int numWallRegions;
         int nsrflistWallRegions[MAXREGIONS+1];
         int nWallETagID;
         int nWallhTagID;
         int iwalldamp;
         int iwallsupp;
         int indsurf;
         int iringdamp;
         int iringsupp;
         int imeasdist;
         int idistancenudge;
         int iinitialprestress;
         int iupdateprestress;
         int iuseBET;
         int numBETFields;
         int iestimator;
         int iheart;
         int numControlledCoronarySrfs;
         int indicesOfCoronarySurfaces[MAXSURF+1];
         int numNetlistLPNSrfs;
         int indicesOfNetlistSurfaces[MAXSURF+1];
         int numLoopClosingCircuits;
         int inputHRandSP;
         int geombcHasObservationFields;
         int geombcHasNodeTags;
         int pureZeroDSimulation;
      } nomodule;


      extern struct {
    	  int numGRCRSrfs;
    	  int nsrflistGRCR[MAXSURF+1];
    	  int igrcrfile;
      } grcrbccom;

      extern struct {
         int seqsize;
         int stepseq[100];
      } sequence;

      extern struct {
         MPI_Fint iNewComm;
      } newcom;

#ifdef __cplusplus
   }
#endif


#endif
//}
