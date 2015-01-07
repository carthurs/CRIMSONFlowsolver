/*===========================================================================
 *
 * "usr.c":  user's function
 *
 *===========================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include "usr.hxx"
 extern "C"
 {
    #include "les.h"
 }
#include "cvSolverIO.h"

/*===========================================================================
 *
 * "usrNew":  Put all the values in usrHd
 *
 * From FORTRAN
 *
 *	integer		usr(100)
 *	dimension	aperm(numnp,nperm)
 *	...
 *	call usrnew( usr, aperm, ..., numnp, ...)
 *	
 *
 *===========================================================================
 */
#include "mpi.h"
#include "common_c.h"

static int lmNum = 0; 
static LesHd lesArray[8];

#ifdef ACUSIM_LESLIB_VER_1_4

extern "C"
{
    void   usrNew(	UsrHd	  usrHd,
                            int*      eqnType,
                            double*	  aperm,
                            double*	  atemp,
                            double*   resf,
                            double*   solinc,
                            double*   flowDiag,
                            double*   sclrDiag,
                            double*   lesP,
                            double*   lesQ,
                            int*      iBC,
                            double*   BC,
                            int*      iper,
                            int*      ilwork,
                            int*      numpe,
                            int*      nNodes,
                            int*      nenl,
                            int*	  nPermDims,
                            int*	  nTmpDims,
                            int*	  rowp,
                            int*	  colm,
                            double*   lhsK,
                            double*   lhsP,
                            double*   lhsS,
                            int       nnz_tot
        )
    {
        char*	funcName = "usrNew" ;	/* function name		*/
        
    /*---------------------------------------------------------------------------
     * Stick the parameters
     *---------------------------------------------------------------------------
     */
        usrHd->eqnType      = *eqnType ;
        usrHd->aperm	= aperm ;
        usrHd->atemp	= atemp ;
        usrHd->resf         = resf ;
        usrHd->solinc       = solinc ;
        usrHd->flowDiag     = flowDiag ;
        usrHd->sclrDiag     = sclrDiag ;
        usrHd->lesP         = lesP ;
        usrHd->lesQ         = lesQ ;
        usrHd->iBC          = iBC  ;
        usrHd->BC           = BC   ;
        usrHd->iper         = iper ;
        usrHd->ilwork       = ilwork ;
        usrHd->numpe        = *numpe ;
        usrHd->nNodes	= *nNodes ;
        usrHd->nenl         = *nenl ;
        usrHd->nPermDims	= *nPermDims ;
        usrHd->nTmpDims	= *nTmpDims ;
        usrHd->rowp	        = rowp ;
        usrHd->colm	        = colm ;
        usrHd->lhsK	        = lhsK ;
        usrHd->lhsP	        = lhsP ;
        usrHd->lhsS         = lhsS ;
        usrHd->nnz_tot      = nnz_tot ;

    } /* end of usrNew() */
}

#else

extern "C"
{
    void   usrNew(	UsrHd	  usrHd,
                            int*      eqnType,
                            double*	  aperm,
                            double*	  atemp,
                            double*   resf,
                            double*   solinc,
                            double*   flowDiag,
                            double*   sclrDiag,
                            double*   lesP,
                            double*   lesQ,
                            int*      iBC,
                            double*   BC,
                            int*      iper,
                            int*      ilwork,
                            int*      numpe,
                            int*      nNodes,
                            int*      nenl,
                            int*	  nPermDims,
                            int*	  nTmpDims,
                            int*	  rowp,
                            int*	  colm,
                            double*   lhsK,
                            double*   lhsP,
                            double*   lhsS,
                            int       nnz_tot,
                            double*   CGsol
        )
    {
        char*	funcName = "usrNew" ;	/* function name		*/
        
    /*---------------------------------------------------------------------------
     * Stick the parameters
     *---------------------------------------------------------------------------
     */
        usrHd->eqnType      = *eqnType ;
        usrHd->aperm	= aperm ;
        usrHd->atemp	= atemp ;
        usrHd->resf         = resf ;
        usrHd->solinc       = solinc ;
        usrHd->flowDiag     = flowDiag ;
        usrHd->sclrDiag     = sclrDiag ;
        usrHd->lesP         = lesP ;
        usrHd->lesQ         = lesQ ;
        usrHd->iBC          = iBC  ;
        usrHd->BC           = BC   ;
        usrHd->iper         = iper ;
        usrHd->ilwork       = ilwork ;
        usrHd->numpe        = *numpe ;
        usrHd->nNodes	= *nNodes ;
        usrHd->nenl         = *nenl ;
        usrHd->nPermDims	= *nPermDims ;
        usrHd->nTmpDims	= *nTmpDims ;
        usrHd->rowp	        = rowp ;
        usrHd->colm	        = colm ;
        usrHd->lhsK	        = lhsK ;
        usrHd->lhsP	        = lhsP ;
        usrHd->lhsS         = lhsS ;
        usrHd->nnz_tot      = nnz_tot ;
        usrHd->CGsol        = CGsol;
    } /* end of usrNew() */
}

#endif

/*===========================================================================
 *
 * "usrPointer":  Get the pointer
 *
 *===========================================================================
 */
Real*
usrPointer(	UsrHd	usrHd,
            Integer	id,
            Integer	offset,
            Integer	nDims )
{
    char*	funcName = "usrPointer";/* function name		*/
    Real*	pnt ;			/* pointer			*/

/*---------------------------------------------------------------------------
 * Get the head of the memory
 *---------------------------------------------------------------------------
 */
    if ( id == LES_RES_PNT ) {
	
        pnt	= usrHd->resf ;
        id	= 0 ;

    } else if ( id == LES_SOL_PNT ) {

        pnt	= usrHd->solinc ;
        id	= 0 ;

    } else if ( id < 0 ) {

        pnt	= usrHd->aperm ;
        id	= id + usrHd->nPermDims ;

    } else {

        pnt	= usrHd->atemp ;
        id	= id ;

    }
/*---------------------------------------------------------------------------
 * Get the offset
 *---------------------------------------------------------------------------
 */
    pnt		= pnt + (id + offset) * usrHd->nNodes ;

/*---------------------------------------------------------------------------
 * Return the pointer
 *---------------------------------------------------------------------------
 */
    return( pnt ) ;

} /* end of usrPointer() */

#ifdef ACUSIM_LESLIB_VER_1_4

#ifdef intel
#define myfLesNew MYFLESNEW
#define myflessolve_ MYFLESSOLVE
#define savelesrestart_ SAVELESRESTART
#define readlesrestart_ READLESRESTART
#define solverlicenseserver_ SOLVERLICENSESERVER

void    myfLesNew(  Integer*   lesId,
                    Integer*   lmport,
                    Integer*   eqnType,
                    Integer*   nDofs,
                    Integer*   minIters,
                    Integer*   maxIters,
                    Integer*   nKvecs,
                    Integer*   prjFlag,
                    Integer*   nPrjs,
                    Integer*   presPrjFlag,
                    Integer*   nPresPrjs,
                    Real*      tol,
                    Real*      presTol,
                    Integer*   verbose,
                    Real*      stats,
                    Integer*   nPermDims,
                    Integer*   nTmpDims,
                    char*     	fileName,
                    Integer   	len		) {

        lesArray[ *lesId ] = lesNew( fileName, *lmport, &lmNum, *eqnType,
                                     *nDofs, *minIters, *maxIters, *nKvecs, 
                                     *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,
                                     *tol, *presTol, *verbose, stats, nPermDims,
                                     nTmpDims );
    return ;}
/* the following is a fake function that was required when we moved to
   a C++ main on in the MS Visual Studio environment.  It fails to
   link because it is looking for this function
*/
//void  _CrtDbgReport() {
//    return ;}

double __vcos_(double fg) { fflush(stdout); printf(" vcos got called \n"); fflush(stdout);}
double __vlog_(double fg)  { fflush(stdout); printf(" vlog got called \n"); fflush(stdout);}


#else /* we are in unix land... whew.  secretly we have equivalenced fileName and  */

/* #ifdef LINUX*/
/* void flush_(int* junk ){ return; }*/
/* #endif*/
extern "C"
{
    void    myflesnew_(	     Integer*	lesId,
                             Integer*	lmport,                             
                             Integer*	eqnType,
                             Integer*	nDofs,
                             Integer*	minIters,
                             Integer*	maxIters,
                             Integer*	nKvecs,
                             Integer*	prjFlag,
                             Integer*	nPrjs,
                             Integer*	presPrjFlag,
                             Integer*	nPresPrjs,
                             Real*	    tol,
                             Real*     	presTol,
                             Integer*	verbose,
                             Real*     	stats,
                             Integer*	nPermDims,
                             Integer*	nTmpDims,
                             char*      lmhost,
                             Integer*   len             ) {
        int procId;
        MPI_Comm iNewComm_C = MPI_Comm_f2c(newcom.iNewComm);
        MPI_Comm_rank( iNewComm_C, &procId ) ;
        if(lmNum==0){
            if(procId==0){
                lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                             *nDofs, *minIters, *maxIters, *nKvecs, 
                                             *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,
                                             *tol, *presTol, *verbose, stats, nPermDims,
                                             nTmpDims );
                MPI_Bcast( &lmNum, 1, MPI_INT, 0, iNewComm_C ) ;
            } else {
                MPI_Bcast( &lmNum, 1, MPI_INT, 0, iNewComm_C ) ;
                lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                             *nDofs, *minIters, *maxIters, *nKvecs, 
                                             *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,
                                             *tol, *presTol, *verbose, stats, nPermDims,
                                             nTmpDims );
            }
        } else {
            lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                         *nDofs, *minIters, *maxIters, *nKvecs, 
                                         *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs,
                                         *tol, *presTol, *verbose, stats, nPermDims,
                                         nTmpDims );
        }
        return ;
    }
}

#endif

#else

#ifdef intel
#define myfLesNew MYFLESNEW
#define myflessolve_ MYFLESSOLVE
#define savelesrestart_ SAVELESRESTART
#define readlesrestart_ READLESRESTART
#define solverlicenseserver_ SOLVERLICENSESERVER

extern "C" {
    void    myfLesNew(  Integer*   lesId,
                        Integer*   lmport,
                        Integer*   eqnType,
                        Integer*   nDofs,
                        Integer*   minIters,
                        Integer*   maxIters,
                        Integer*   nKvecs,
                        Integer*   prjFlag,
                        Integer*   nPrjs,
                        Integer*   presPrjFlag,
                        Integer*   nPresPrjs,
                        Real*      tol,
                        Real*      presTol,
                        Integer*   verbose,
                        Real*      stats,
                        Integer*   nPermDims,
                        Integer*   nTmpDims,
                        char*     	fileName        ) {

            int ppePreCond=0; /* =1 will invoke precondition, =0 will be as v1.4 */
            lesArray[ *lesId ] = lesNew( fileName, *lmport, &lmNum, *eqnType,
                                         *nDofs, *minIters, *maxIters, *nKvecs, 
                                         *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs, ppePreCond,
                                         *tol, *presTol, *verbose, stats, nPermDims,
                                         nTmpDims );
        return ;}
    /* the following is a fake function that was required when we moved to
       a C++ main on in the MS Visual Studio environment.  It fails to
       link because it is looking for this function
    */
    void  _CrtDbgReport() {
        return ;}

    double __vcos_(double fg) { fflush(stdout); printf(" vcos got called \n"); fflush(stdout);}
    double __vlog_(double fg)  { fflush(stdout); printf(" vlog got called \n"); fflush(stdout);}
}


#else /* we are in unix land... whew.  secretly we have equivalenced fileName and  */

/* #ifdef LINUX*/
/* void flush_(int* junk ){ return; }*/
/* #endif*/
extern "C"
{
    void    myflesnew_(	     Integer*	lesId,
    		MPI_Bcast   Integer*	lmport,
                             Integer*	eqnType,
                             Integer*	nDofs,
                             Integer*	minIters,
                             Integer*	maxIters,
                             Integer*	nKvecs,
                             Integer*	prjFlag,
                             Integer*	nPrjs,
                             Integer*	presPrjFlag,
                             Integer*	nPresPrjs,
                             Real*	    tol,
                             Real*     	presTol,
                             Integer*	verbose,
                             Real*     	stats,
                             Integer*	nPermDims,
                             Integer*	nTmpDims,
                             char*      lmhost          ) {
        int procId;
        int ppePreCond=0; /* =1 will invoke precondition, =0 will be as v1.4 */
        MPI_Comm iNewComm_C = MPI_Comm_f2c(newcom.iNewComm);
        MPI_Comm_rank( iNewComm_C, &procId ) ;
        if(lmNum==0){
            if(procId==0){
                lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                             *nDofs, *minIters, *maxIters, *nKvecs, 
                                             *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs, ppePreCond,
                                             *tol, *presTol, *verbose, stats, nPermDims,
                                             nTmpDims );
                MPI_Bcast( &lmNum, 1, MPI_INT, 0, iNewComm_C ) ;
            } else {
                MPI_Bcast( &lmNum, 1, MPI_INT, 0, iNewComm_C ) ;
                lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                             *nDofs, *minIters, *maxIters, *nKvecs, 
                                             *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs, ppePreCond,
                                             *tol, *presTol, *verbose, stats, nPermDims,
                                             nTmpDims );
            }
        } else {
            lesArray[ *lesId ] = lesNew( lmhost, *lmport, &lmNum, *eqnType,
                                         *nDofs, *minIters, *maxIters, *nKvecs, 
                                         *prjFlag, *nPrjs, *presPrjFlag, *nPresPrjs, ppePreCond,
                                         *tol, *presTol, *verbose, stats, nPermDims,
                                         nTmpDims );
        }
        return ;
    }
}

#endif

#endif

extern "C"
{
    void savelesrestart_( Integer* lesId,
                     Real*    aperm,
                     Integer* nshg,
                     Integer* myrank,
                     Integer* lstep,
                     Integer* nPermDims ) {

        int nPrjs, PrjSrcId;
        int nPresPrjs, PresPrjSrcId;
        char filename[255];
        int fileHandle=0;
        int iarray[3];
        int size, nitems;
        double* projVec;
        int i, j, count;

        sprintf( filename,"restart.%d.%d", *lstep, *myrank+1 );
        openfile_( filename, "append", &fileHandle );
        
        nPrjs = (Integer) lesGetPar( lesArray[ *lesId ], LES_ACT_PRJS );
        PrjSrcId = (Integer) lesGetPar( lesArray[ *lesId ], LES_PRJ_VEC_ID );

        if ( PrjSrcId < 0 ) PrjSrcId += *nPermDims;

        projVec = (double*)malloc( nPrjs * ( *nshg ) * sizeof( double ) );

        count = 0;
        for( i = PrjSrcId; i < PrjSrcId+nPrjs; i ++ ) {
            for( j = 0 ; j < *nshg; j++ ) {
                projVec[ count++ ] = aperm[ (*nshg) * i + j ];
            }
        }

        iarray[ 0 ] = *nshg;
        iarray[ 1 ] = nPrjs;
        nitems = 2;
        size = (*nshg)*nPrjs;

        writeheader_( &fileHandle, "projection vectors ", (void*)iarray, 
                      &nitems, &size, "double", outpar.iotype );
        nitems = size;
        writedatablock_( &fileHandle, "projection vectors ", (void*)projVec, 
                         &nitems, "double", outpar.iotype );
        free(projVec);

        /*************************************************************************/

        nPresPrjs = (Integer) lesGetPar( lesArray[ *lesId ], LES_ACT_PRES_PRJS );
        PresPrjSrcId =(Integer)lesGetPar( lesArray[ *lesId ], LES_PRES_PRJ_VEC_ID );
        if ( PresPrjSrcId < 0 ) PresPrjSrcId += *nPermDims;

        projVec = (double*)malloc( nPresPrjs * ( *nshg ) * sizeof( double ) );

        count = 0;
        for( i = PresPrjSrcId; i < (PresPrjSrcId + nPresPrjs) ; i ++ ) {
            for( j = 0 ; j < *nshg; j++ ) {
                projVec[ count++ ] = aperm[ (*nshg) * i + j ];
            }
        }

        iarray[ 0 ] = *nshg;
        iarray[ 1 ] = nPresPrjs;
        nitems = 2;
        size = (*nshg)*nPresPrjs;

        writeheader_( &fileHandle, "pressure projection vectors ", (void*)iarray, 
                      &nitems, &size, "double", outpar.iotype );
        nitems = size;

        writedatablock_( &fileHandle, "pressure projection vectors ", 
                         (void*)projVec, &nitems, "double", outpar.iotype );
        free( projVec);

        closefile_( &fileHandle, "append" );
    }

    void readlesrestart_( Integer* lesId,
                     Real*    aperm,
                     Integer* nshg,
                     Integer* myrank,
                     Integer* lstep ,
                     Integer* nPermDims ) {

        int nPrjs, PrjSrcId;
        int nPresPrjs, PresPrjSrcId;
        char filename[255];
        int fileHandle=0;
        int iarray[3]={-1,-1,-1};
        int size, nitems;
        int itwo=2;
        int lnshg;
        double* projVec;
        int i,j,count;

        sprintf( filename,"restart.%d.%d", *lstep, *myrank+1 );
        openfile_( filename, "read", &fileHandle );

        if ( fileHandle == 0 ) return;

        readheader_( &fileHandle, "projection vectors", (void*)iarray, 
                     &itwo, "integer", outpar.iotype );

        if ( iarray[0] != *nshg ) {
            closefile_( &fileHandle, "read" );
            printf("projection vectors are being initialized to zero (SAFE)\n");
            return;
        }

        lnshg = iarray[ 0 ] ;
        nPrjs = iarray[ 1 ] ;

        size = (*nshg)*nPrjs;
        projVec = (double*)malloc( size * sizeof( double ));

        readdatablock_( &fileHandle, "projection vectors", (void*)projVec, 
                        &size, "double", outpar.iotype );
        
        lesSetPar( lesArray[ *lesId ], LES_ACT_PRJS, (Real) nPrjs );
        PrjSrcId = (Integer) lesGetPar( lesArray[ *lesId ], LES_PRJ_VEC_ID );
        if ( PrjSrcId < 0 ) PrjSrcId += *nPermDims;

        count = 0;
        for( i = PrjSrcId; i < PrjSrcId+nPrjs; i ++ ) {
            for( j = 0 ; j < *nshg; j++ ) {
                aperm[ (*nshg) * i + j ] = projVec[ count++ ] ;
            }
        }

        free( projVec );

        /************************************************************************/


        readheader_( &fileHandle, "pressure projection vectors", (void*)iarray, 
                     &itwo, "integer", outpar.iotype );

        lnshg = iarray[ 0 ] ;
        nPresPrjs = iarray[ 1 ] ;

        if ( lnshg != *nshg )  { 
            closefile_( &fileHandle, "read" );
            printf("projection vectors are being initialized to zero (SAFE)\n");
            return;
        }

        size = (*nshg)*nPresPrjs;
        projVec = (double*)malloc( size * sizeof( double ));
        
        readdatablock_( &fileHandle, "pressure projection vectors", (void*)projVec, 
                        &size, "double", outpar.iotype );
        
        lesSetPar( lesArray[ *lesId ], LES_ACT_PRES_PRJS, (Real) nPresPrjs );
        PresPrjSrcId=(Integer)lesGetPar( lesArray[ *lesId ], LES_PRES_PRJ_VEC_ID );
        if ( PresPrjSrcId < 0 ) PresPrjSrcId += *nPermDims;
        
        count = 0;
        for( i = PresPrjSrcId; i < PresPrjSrcId+nPresPrjs; i ++ ) {
            for( j = 0 ; j < *nshg; j++ ) {
                aperm[ (*nshg) * i + j ] = projVec[ count++ ] ;
            }
        }

        free( projVec );

        closefile_( &fileHandle, "read" );
    }

    void  myflessolve_( Integer* lesId,
                        UsrHd    usrHd){
        lesSolve( lesArray[ *lesId ], usrHd );
    }

    void getSol ( UsrHd usrHd,
                  double* Dy  )
    {

         Dy = usrHd->solinc;

    }


    int solverlicenseserver_(char key[]){
        char* env_server_name;
        env_server_name = getenv("LES_LICENSE_SERVER");
        if(env_server_name) strcpy(key, env_server_name);
        else { 
            fprintf(stderr,
    				"environment variable LES_LICENSE_SERVER not defined \n");
            //fprintf(stderr,"using server1 as default \n");
            //strcpy(key, "server1.scorec.rpi.edu");
        }
        return 1;
    }
}
