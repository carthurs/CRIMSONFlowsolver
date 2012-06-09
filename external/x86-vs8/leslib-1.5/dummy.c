typedef struct _Usr* UsrHd ;
#include "les.h"


Void		lesPrepDiag(		UsrHd		usrHd		) 
{ return ; }
Void		lesDiagScaleCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		diagOff,
					Integer		nDims		) 
{ return ; }

Void		lesZero(		UsrHd		usrHd,
					Integer		dstId, 
					Integer		nDims		) 
{ return ; }
Void		lesCp(			UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nDims		) 
{ return ; }
Void		lesScale(		UsrHd		usrHd,
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) 
{ return ; }
Void		lesScaleCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) 
{ return ; }
Void		lesAdd(			UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Integer		nDims		) 
{ return ; }
Void		lesSub(			UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Integer		nDims		) 
{ return ; }
Real		lesDot1(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		nDims		) 
{ return(1) ; }
Real		lesDot2(		UsrHd		usrHd,
					Integer		src1Id, 
					Integer		src2Id, 
					Integer		nDims		) 
{ return(1) ; }
Void		lesDaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) 
{ return ; }
Void		lesDxpay(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) 
{ return ; }
Void		lesInv(			UsrHd		usrHd,
					Integer		dstId, 
					Integer		nDims		) 
{ return ; }
Void		lesBlkDot2(		UsrHd		usrHd,
					Integer		src1Id, 
					Integer		src2Id, 
					Real*		values,
					Integer		mDims,
					Integer		nDims		) 
{ return ; }
Void		lesBlkDaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) 
{ return ; }
Void		lesBlkDyeax(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) 
{ return ; }
Void		lesBlkDmaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) 
{ return ; }
Void		lesVdimCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		nDims		) 
{ return ; }
Void		lesVdimDot2(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id, 
					Real*		coef,
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDims		) 
{ return ; }
Void		lesVdimDaxpy(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Real*		coef,
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		nDims		) 
{ return ; }

Void		lesApG(			UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }
Void		lesApKG(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id,
					Integer		dstId, 
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }
Void		lesApNGt(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }
Void		lesApNGtC(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id,
					Integer		dstId, 
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }
Void		lesApFull(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }
Void		lesApSclr(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }

Void		lesPrecPPE(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) 
{ return ; }

Void dgees_()  { return ; }
Void dgeev_()  { return ; }
Void dgesvd_() { return ; }
Void dggev_()  { return ; }
Void dlarnv_() { return ; }
Void dsyev_()  { return ; }
Void dsygv_()  { return ; }
Void dtrexc_() { return ; }
