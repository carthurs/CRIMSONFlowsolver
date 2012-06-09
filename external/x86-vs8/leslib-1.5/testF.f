	program test
c
c.... declare
c
	integer	usr(100)
c
	real*8	stats(6)
c
	integer	lesId,		eqnType,	nDofs,		
     &		minIters,	maxIters,	nKvecs,		
     &		prjFlag,	nPrjs,		presPrjFlag,	
     &		nPresPrjs,	presPrecFlag,	verbose
	real*8	tol,		presTol
c
	integer	nTempDims,	nPermDims
c
	character*16	host
	integer		port,	lmNum
c
c.... initialize
c
	lesId		= 1
	host		= "pelican"
	port		= 41994
	lmNum		= 0
	eqnType		= 1
	nDofs		= 4
	minIters	= 10
	maxIters	= 400
	nKvecs		= 10
	prjFlag		= 0
	nPrjs		= 5
	presPrjFlag	= 1
	nPresPrjs	= 10
	presPrecFlag	= 0
	tol		= 1.d-1
	presTol		= 1.d-2
	verbose		= 1
c
	call fLesNew(	lesId,		host,		port,		
     1			lmNum,		eqnType,	nDofs,		
     2			minIters,	maxIters,	nKvecs, 	
     3			prjFlag,	nPrjs,		presPrjFlag,	
     4			nPresPrjs,	presPrecFlag,	tol,		
     5			presTol,	verbose,	stats,		
     6			nPermDims,	nTempDims	)
c
c.... call solver
c
	call fLesSolve(	lesId,		usr				)
c
c.... free up the memory
c
	call fLesFree(	lesId						)
c
c.... end
c
	end
