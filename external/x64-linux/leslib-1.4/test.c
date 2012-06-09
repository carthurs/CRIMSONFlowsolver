/****************************************************************************
**  Copyright (c) 1994-2003 ACUSIM Software, Inc.
**  All rights reserved.
**  This source code is confidential and may not be disclosed.
****************************************************************************/

/*===========================================================================
**
** "main.c":  The main program of the demo
**
** Original: Farzin Shakib (Aug 94)
**===========================================================================
*/

/*===========================================================================
 *
 * Include "les.h"
 *
 * UsrHd is the handle to user supplied structure containing information
 * needed to carry out the lower level operations.
 *
 *===========================================================================
 */
typedef struct _Usr* UsrHd ;
#include "les.h"

/*===========================================================================
 *
 * "main":  Main program
 *
 *===========================================================================
 */
int main( int, char** ) ;
int
main(	int	argc,
	char*	argv[] )
{
    char*	funcName = "main" ;	/* function name		*/
    Integer	eqnType ;		/* equation type		*/
    Integer	lmNum ;			/* license number		*/
    Integer	maxIters ;		/* Max number of iterations	*/
    Integer	minIters ;		/* Min number of iterations	*/
    Integer	nDofs ;			/* No. dofs			*/
    Integer	nKvecs ;		/* No. Krylov vectors		*/
    Integer	nPermDims ;		/* No. permanent vecs of nNodes	*/
    Integer	nPresPrjs ;		/* No. pressure projections	*/
    Integer	nPrjs ;			/* No. projections		*/
    Integer	nTmpDims ;		/* No. temporary vecs of nNodes	*/
    Integer	port ;			/* port				*/
    Integer	presPrjFlag ;		/* pressure projection flag	*/
    Integer	prjFlag ;		/* projection flag		*/
    Integer	verbose ;		/* Verbose level		*/
    LesHd	lesHd ;			/* les driver			*/
    Real	presTol ;		/* pressure projection tol.	*/
    Real	stats[6] ;		/* solver stats			*/
    Real	tol ;			/* linear solver tol.		*/
    UsrHd	usrHd ;			/* user structure		*/
    char*	host = "pelican" ;	/* license host			*/

/*---------------------------------------------------------------------------
 * Call les for flow
 *---------------------------------------------------------------------------
 */
    usrHd	= NULL ;

    eqnType	= LES_EQN_FLOW ;
    nDofs	= 4 ;
    minIters	= 10 ;
    maxIters	= 400 ;
    nKvecs	= 10 ;
    prjFlag	= 0 ;
    nPrjs	= 5 ;
    presPrjFlag	= 1 ;
    nPresPrjs	= 10 ;
    tol		= 1.e-1 ;
    presTol	= 1.e-2 ;
    verbose	= 1 ;
    port	= 41994 ;
    lmNum	= 0 ;

    sysMess( "calling lesNew" ) ;
    lesHd	= lesNew(	host,
				port,
				&lmNum,
				eqnType,
				nDofs,
				minIters,
				maxIters,
				nKvecs,
				prjFlag,
				nPrjs,
				presPrjFlag,
				nPresPrjs,
				tol,
				presTol,
				verbose,
				stats,
				&nPermDims,
				&nTmpDims				) ;

    sysMess( "calling lesSolve" ) ;
    lesSolve(			lesHd,
				usrHd					) ;

    sysMess( "calling lesFree" ) ;
    lesFree(			lesHd					) ;

} /* end of main() */

