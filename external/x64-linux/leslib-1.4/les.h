/****************************************************************************
**  Copyright (c) 1994-2003 ACUSIM Software, Inc.
**  All rights reserved.
**  This source code is confidential and may not be disclosed.
****************************************************************************/

/*===========================================================================
**
** "les.h":  Linear Equation Solvers.
**
** Original: Farzin Shakib (Aug 94)
**===========================================================================
*/

#ifndef	__LES_H__
#define	__LES_H__

/*===========================================================================
 *
 * Get the needed include files
 *
 *===========================================================================
 */
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*===========================================================================
 *
 * Main data types
 *
 *===========================================================================
 */
typedef	int	Integer	;		/* integer type			*/
typedef	double	Real ;			/* real    type			*/
typedef	float	Float ;			/* float   type			*/
typedef	char*	String ;		/* string  type			*/
typedef	void*	Data ;			/* data    type			*/
typedef	void	Void ;			/* void    type			*/

#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

typedef enum { 
    FALSE = 0 , 
    TRUE  = 1 
} Bool ;				/* boolean type			*/

/*===========================================================================
 *
 * Nil
 *
 *===========================================================================
 */
#define	Nil(T)	((T)NULL)

/*===========================================================================
 *
 * "Les":  Acusim LES data structure
 *
 *===========================================================================
 */
typedef struct _Les* LesHd ;

/*===========================================================================
 *
 * Equation type
 *
 *===========================================================================
 */
#define	LES_EQN_FLOW	1			/* Solve flow vel/pres	*/
#define	LES_EQN_SCLR	2			/* Solve multi-scalar	*/

/*===========================================================================
 *
 * Parameter names
 *
 *===========================================================================
 */
#define	LES_ACT_PRJS		1	/* No. active projection vecs.	*/
#define	LES_ACT_PRES_PRJS	2	/* No. act. pres. prj. vecs.	*/
#define	LES_PRJ_VEC_ID		3	/* Prj vec srcId into permVec	*/
#define	LES_PRES_PRJ_VEC_ID	4	/* Pres. Prj vec srcId 		*/

/*===========================================================================
 *
 * Pointers
 *
 *===========================================================================
 */
#define	LES_RES_PNT		100000		/* residual pointer	*/
#define	LES_SOL_PNT		200000		/* solution pointer	*/

/*===========================================================================
 *
 * NT calling convention
 *
 *===========================================================================
 */
#if defined(_WIN32) || defined(WIN32) || defined(WINNT)
#define WIN_API __cdecl
#else
#define WIN_API
#endif

/*===========================================================================
 *
 * Acusim callable funtions in C
 *
 *===========================================================================
 */
LesHd WIN_API	lesNew(			char*		lmHost,
					Integer		lmPort,
					Integer*	lmNum,
					Integer		eqnType,
					Integer		nDofs,
					Integer		minIters,
					Integer		maxIters,
					Integer		nKvecs,
					Integer		prjFlag,
					Integer		nPrjs,
					Integer		presPrjFlag,
					Integer		nPresPrjs,
					Real		tol,
					Real		presTol,
					Integer		verbose,
					Real*		stats,
					Integer*	nPermDims,
					Integer*	nTmpDims	) ;
Void WIN_API	lesFree(		LesHd		lesHd		) ;
Void WIN_API	lesSolve(		LesHd		lesHd,
					UsrHd		usrHd		) ;
Real WIN_API	lesGetPar(		LesHd		lesHd,
					Integer		parName		) ;
Void WIN_API	lesSetPar(		LesHd		lesHd,
					Integer		parName,
					Real		parVal		) ;

/*===========================================================================
 *
 * Fortran calling convention
 *
 *===========================================================================
 */
#if	defined(ACUSIM_SGI)     ||	defined(ACUSIM_SGI64)	|| \
	defined(ACUSIM_HAL)	||	defined(ACUSIM_SUN)	|| \
	defined(ACUSIM_ALPHA)	||	defined(ACUSIM_LINUX)	|| \
	defined(ACUSIM_LINUXIPF)
#define	fLesNew		flesnew_
#define	fLesFree	flesfree_
#define	fLesSolve	flessolve_
#elif	defined(ACUSIM_HP)
#define	fLesNew		flesnew
#define	fLesFree	flesfree
#define	fLesSolve	flessolve
#elif	defined(ACUSIM_NT)
#define	fLesNew		FLESNEW
#define	fLesFree	FLESFREE
#define	fLesSolve	FLESSOLVE
#else /* dummy */
#define	fLesNew		fLesNewXX
#define	fLesFree	fLesFreeXX
#define	fLesSolve	fLesSolveXX
#endif

/*===========================================================================
 *
 * Acusim callable funtions in Fortran
 *
 *===========================================================================
 */

Void WIN_API	fLesNew(		Integer*	lesId,
					char*		lmhost,
					Integer*	lmport,
					Integer*	lmNum,
					Integer*	eqnType,
					Integer*	nDofs,
					Integer*	minIters,
					Integer*	maxIters,
					Integer*	nKvecs,
					Integer*	prjFlag,
					Integer*	nPrjs,
					Integer*	presPrjFlag,
					Integer*	nPresPrjs,
					Real*		tol,
					Real*		presTol,
					Integer*	verbose,
					Real*		stats,
					Integer*	nPermDims,
					Integer*	nTmpDims,
					Integer		len_lmHost	) ;
Void WIN_API	fLesFree(		Integer*	lesId		) ;
Void WIN_API	fLesSolve(		Integer*	lesId,
					UsrHd		usrHd		) ;

/*===========================================================================
 *
 * Functions to be provided by user
 *
 *===========================================================================
 */
Void WIN_API	lesPrepDiag(		UsrHd		usrHd		) ;
Void WIN_API	lesDiagScaleCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		diagOff,
					Integer		nDims		) ;

Void WIN_API	lesZero(		UsrHd		usrHd,
					Integer		dstId, 
					Integer		nDims		) ;
Void WIN_API	lesCp(			UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nDims		) ;
Void WIN_API	lesScale(		UsrHd		usrHd,
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void WIN_API	lesScaleCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void WIN_API	lesAdd(			UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Integer		nDims		) ;
Void WIN_API	lesSub(			UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Integer		nDims		) ;
Real WIN_API	lesDot1(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		nDims		) ;
Real WIN_API	lesDot2(		UsrHd		usrHd,
					Integer		src1Id, 
					Integer		src2Id, 
					Integer		nDims		) ;
Void WIN_API	lesDaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void WIN_API	lesDxpay(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void WIN_API	lesInv(			UsrHd		usrHd,
					Integer		dstId, 
					Integer		nDims		) ;
Void WIN_API	lesBlkDot2(		UsrHd		usrHd,
					Integer		src1Id, 
					Integer		src2Id, 
					Real*		values,
					Integer		mDims,
					Integer		nDims		) ;
Void WIN_API	lesBlkDaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) ;
Void WIN_API	lesBlkDyeax(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) ;
Void WIN_API	lesBlkDmaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) ;
Void WIN_API	lesVdimCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		nDims		) ;
Void WIN_API	lesVdimDot2(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id, 
					Real*		coef,
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDims		) ;
Void WIN_API	lesVdimDaxpy(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Real*		coef,
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		nDims		) ;

Void WIN_API	lesApG(			UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void WIN_API	lesApKG(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id,
					Integer		dstId, 
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void WIN_API	lesApNGt(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void WIN_API	lesApNGtC(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id,
					Integer		dstId, 
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void WIN_API	lesApFull(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void WIN_API	lesApSclr(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;

/*===========================================================================
 *
 * End of the file
 *
 *===========================================================================
 */

#endif	/* __LES_H__ */
