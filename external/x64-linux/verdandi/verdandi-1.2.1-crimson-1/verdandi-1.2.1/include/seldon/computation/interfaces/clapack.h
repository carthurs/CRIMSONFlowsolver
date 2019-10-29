// Headers from CLAPACK, downloaded at http://www.netlib.org/clapack/.

// Modifications (by Vivien Mallet):
// Replacements:
//    integer       --> LAPACK_INTEGER
//    real          --> LAPACK_REAL
//    doublereal    --> LAPACK_DOUBLEREAL
//    complex       --> LAPACK_COMPLEX
//    doublecomplex --> LAPACK_DOUBLECOMPLEX
//    logical       --> LAPACK_LOGICAL
//    L_fp          --> LAPACK_L_FP
//    ftnlen        --> LAPACK_FTNLEN

#ifndef __CLAPACK_H
#define __CLAPACK_H
 
/* Subroutine */ int cbdsqr_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ncvt, LAPACK_INTEGER *
			     nru, LAPACK_INTEGER *ncc, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_COMPLEX *vt, LAPACK_INTEGER *ldvt, 
			     LAPACK_COMPLEX *u, LAPACK_INTEGER *ldu, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbbrd_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *ncc,
			     LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_COMPLEX *pt, LAPACK_INTEGER *ldpt, 
			     LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbcon_(char *norm, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL 
			     *colcnd, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *afb, LAPACK_INTEGER *
			     ldafb, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *
			     ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cgbsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_INTEGER *
			    nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *
			    ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl,
			     LAPACK_INTEGER *ku, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *afb,
			     LAPACK_INTEGER *ldafb, LAPACK_INTEGER *ipiv, char *equed, LAPACK_REAL *r__, LAPACK_REAL *c__, 
			     LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL 
			     *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbtf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbtrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgbtrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX 
			     *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgebak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_REAL *scale, LAPACK_INTEGER *m, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cgebal_(char *job, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgebd2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_COMPLEX *tauq, LAPACK_COMPLEX *taup, LAPACK_COMPLEX *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cgebrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_COMPLEX *tauq, LAPACK_COMPLEX *taup, LAPACK_COMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgecon_(char *norm, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL *colcnd, LAPACK_REAL *amax, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cgees_(char *jobvs, char *sort, LAPACK_L_FP select, LAPACK_INTEGER *n, 
			    LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, LAPACK_COMPLEX *w, LAPACK_COMPLEX *vs, 
			    LAPACK_INTEGER *ldvs, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_LOGICAL *
			    bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeesx_(char *jobvs, char *sort, LAPACK_L_FP select, char *
			     sense, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, LAPACK_COMPLEX *
			     w, LAPACK_COMPLEX *vs, LAPACK_INTEGER *ldvs, LAPACK_REAL *rconde, LAPACK_REAL *rcondv, LAPACK_COMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			    LAPACK_INTEGER *lda, LAPACK_COMPLEX *w, LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, 
			    LAPACK_INTEGER *ldvr, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *
			    info);
 
/* Subroutine */ int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *w, LAPACK_COMPLEX *vl, 
			     LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi,
			     LAPACK_REAL *scale, LAPACK_REAL *abnrm, LAPACK_REAL *rconde, LAPACK_REAL *rcondv, LAPACK_COMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgegs_(char *jobvsl, char *jobvsr, LAPACK_INTEGER *n, LAPACK_COMPLEX *
			    a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *
			    beta, LAPACK_COMPLEX *vsl, LAPACK_INTEGER *ldvsl, LAPACK_COMPLEX *vsr, LAPACK_INTEGER *ldvsr, 
			    LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgegv_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			    LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta,
			    LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_COMPLEX *
			    work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgehd2_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgehrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER 
			     *info);
 
/* Subroutine */ int cgelq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgelqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgels_(char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *
			    nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *
			    work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgelsx_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *jpvt, LAPACK_REAL *rcond,
			     LAPACK_INTEGER *rank, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgelsy_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *jpvt, LAPACK_REAL *rcond,
			     LAPACK_INTEGER *rank, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cgeql2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeqlf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeqp3_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *jpvt, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeqpf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *jpvt, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cgeqr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgeqrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgerfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *
			     b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgerq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgerqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgesc2_(LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *
			     rhs, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv, LAPACK_REAL *scale);
 
/* Subroutine */ int cgesv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *
			    lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgesvd_(char *jobu, char *jobvt, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_REAL *s, LAPACK_COMPLEX *u, LAPACK_INTEGER *
			     ldu, LAPACK_COMPLEX *vt, LAPACK_INTEGER *ldvt, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgesvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *
			     ipiv, char *equed, LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgetc2_(LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *
			     ipiv, LAPACK_INTEGER *jpiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgetf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgetrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgetri_(LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *
			     ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgetrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cggbak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_REAL *lscale, LAPACK_REAL *rscale, LAPACK_INTEGER *m, LAPACK_COMPLEX *v, 
			     LAPACK_INTEGER *ldv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggbal_(char *job, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *lscale, 
			     LAPACK_REAL *rscale, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgges_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			    selctg, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *
			    ldb, LAPACK_INTEGER *sdim, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *vsl, 
			    LAPACK_INTEGER *ldvsl, LAPACK_COMPLEX *vsr, LAPACK_INTEGER *ldvsr, LAPACK_COMPLEX *work, LAPACK_INTEGER *
			    lwork, LAPACK_REAL *rwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggesx_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			     selctg, char *sense, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b,
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *sdim, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *
			     vsl, LAPACK_INTEGER *ldvsl, LAPACK_COMPLEX *vsr, LAPACK_INTEGER *ldvsr, LAPACK_REAL *rconde, LAPACK_REAL 
			     *rcondv, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			    LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta,
			    LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_COMPLEX *
			    work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb,
			     LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *
			     vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *lscale, LAPACK_REAL *
			     rscale, LAPACK_REAL *abnrm, LAPACK_REAL *bbnrm, LAPACK_REAL *rconde, LAPACK_REAL *rcondv, LAPACK_COMPLEX 
			     *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, LAPACK_LOGICAL *bwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cggglm_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *x, 
			     LAPACK_COMPLEX *y, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgghrd_(char *compq, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb,
			     LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgglse_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *c__, LAPACK_COMPLEX *d__, 
			     LAPACK_COMPLEX *x, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggqrf_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *taua, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *taub, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *taua, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *taub, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggsvd_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_COMPLEX *a, LAPACK_INTEGER *
			     lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_REAL *alpha, LAPACK_REAL *beta, LAPACK_COMPLEX *u, 
			     LAPACK_INTEGER *ldu, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cggsvp_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER 
			     *ldb, LAPACK_REAL *tola, LAPACK_REAL *tolb, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_COMPLEX *u, 
			     LAPACK_INTEGER *ldu, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, 
			     LAPACK_INTEGER *iwork, LAPACK_REAL *rwork, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cgtcon_(char *norm, LAPACK_INTEGER *n, LAPACK_COMPLEX *dl, LAPACK_COMPLEX *
			     d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *
			     rcond, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgtrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     dl, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *dlf, LAPACK_COMPLEX *df, LAPACK_COMPLEX *
			     duf, LAPACK_COMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *
			     x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cgtsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *dl, LAPACK_COMPLEX *
			    d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgtsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *dl, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *dlf, LAPACK_COMPLEX *
			     df, LAPACK_COMPLEX *duf, LAPACK_COMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *
			     ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgttrf_(LAPACK_INTEGER *n, LAPACK_COMPLEX *dl, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *
			     du, LAPACK_COMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgttrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     dl, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *
			     b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cgtts2_(LAPACK_INTEGER *itrans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_COMPLEX *dl, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *du2, LAPACK_INTEGER *ipiv, 
			     LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int chbev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			    LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			    LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chbevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chbevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, 
			     LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *
			     m, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int chbgst_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			     LAPACK_INTEGER *kb, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *bb, LAPACK_INTEGER *ldbb, 
			     LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chbgv_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			    LAPACK_INTEGER *kb, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *bb, LAPACK_INTEGER *ldbb, 
			    LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int chbgvx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ka, LAPACK_INTEGER *kb, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *bb, 
			     LAPACK_INTEGER *ldbb, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *
			     il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_COMPLEX *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int chbtrd_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_COMPLEX *q, LAPACK_INTEGER *
			     ldq, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int checon_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cheev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *w, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int cheevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *w, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cheevr_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *
			     iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_INTEGER *isuppz, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *
			     lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cheevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *
			     iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int chegs2_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int chegst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int chegv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_REAL *w, 
			    LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chegvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_REAL *w, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chegvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *
			     m, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int cherfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *
			     b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chesv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *a,
			    LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *work,
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chesvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *
			     ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond,
			     LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int chetf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int chetrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int chetrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chetri_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int chetrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int chgeqz_(char *job, char *compq, char *compz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq,
			     LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, 
			    LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int chpevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, 
			     LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int chpevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *ap, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *
			     abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *
			     work, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpgst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *
			     ap, LAPACK_COMPLEX *bp, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpgv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *bp, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			    LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpgvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *bp, LAPACK_REAL *w, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpgvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *bp, LAPACK_REAL *vl, LAPACK_REAL *vu, 
			     LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_COMPLEX *
			     z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int chprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     ap, LAPACK_COMPLEX *afp, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x,
			     LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int chpsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			    ap, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int chpsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *afp, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *
			     ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int chptrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_COMPLEX *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int chptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int chptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int chptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     ap, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int chsein_(char *side, char *eigsrc, char *initv, LAPACK_LOGICAL *
			     select, LAPACK_INTEGER *n, LAPACK_COMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_COMPLEX *w, LAPACK_COMPLEX *
			     vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *
			     m, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *ifaill, LAPACK_INTEGER *ifailr, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int chseqr_(char *job, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo,
			     LAPACK_INTEGER *ihi, LAPACK_COMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_COMPLEX *w, LAPACK_COMPLEX *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int clabrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_COMPLEX *tauq, LAPACK_COMPLEX *taup, 
			     LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_COMPLEX *y, LAPACK_INTEGER *ldy);
 
/* Subroutine */ int clacgv_(LAPACK_INTEGER *n, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx);
 
/* Subroutine */ int clacon_(LAPACK_INTEGER *n, LAPACK_COMPLEX *v, LAPACK_COMPLEX *x, LAPACK_REAL *est, 
			     LAPACK_INTEGER *kase);
 
/* Subroutine */ int clacp2_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int clacpy_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int clacrm_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *rwork);
 
/* Subroutine */ int clacrt_(LAPACK_INTEGER *n, LAPACK_COMPLEX *cx, LAPACK_INTEGER *incx, LAPACK_COMPLEX *
			     cy, LAPACK_INTEGER *incy, LAPACK_COMPLEX *c__, LAPACK_COMPLEX *s);
 
/* Subroutine */ int claed0_(LAPACK_INTEGER *qsiz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_COMPLEX *qstore, LAPACK_INTEGER *ldqs, LAPACK_REAL *rwork,
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int claed7_(LAPACK_INTEGER *n, LAPACK_INTEGER *cutpnt, LAPACK_INTEGER *qsiz, 
			     LAPACK_INTEGER *tlvls, LAPACK_INTEGER *curlvl, LAPACK_INTEGER *curpbm, LAPACK_REAL *d__, LAPACK_COMPLEX *
			     q, LAPACK_INTEGER *ldq, LAPACK_REAL *rho, LAPACK_INTEGER *indxq, LAPACK_REAL *qstore, LAPACK_INTEGER *
			     qptr, LAPACK_INTEGER *prmptr, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *
			     givcol, LAPACK_REAL *givnum, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int claed8_(LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER *qsiz, LAPACK_COMPLEX *
			     q, LAPACK_INTEGER *ldq, LAPACK_REAL *d__, LAPACK_REAL *rho, LAPACK_INTEGER *cutpnt, LAPACK_REAL *z__, 
			     LAPACK_REAL *dlamda, LAPACK_COMPLEX *q2, LAPACK_INTEGER *ldq2, LAPACK_REAL *w, LAPACK_INTEGER *indxp, 
			     LAPACK_INTEGER *indx, LAPACK_INTEGER *indxq, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, 
			     LAPACK_INTEGER *givcol, LAPACK_REAL *givnum, LAPACK_INTEGER *info);
 
/* Subroutine */ int claein_(LAPACK_LOGICAL *rightv, LAPACK_LOGICAL *noinit, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_COMPLEX *w, LAPACK_COMPLEX *v, LAPACK_COMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_REAL *rwork, LAPACK_REAL *eps3, LAPACK_REAL *smlnum, LAPACK_INTEGER *info);
 
/* Subroutine */ int claesy_(LAPACK_COMPLEX *a, LAPACK_COMPLEX *b, LAPACK_COMPLEX *c__, LAPACK_COMPLEX *
			     rt1, LAPACK_COMPLEX *rt2, LAPACK_COMPLEX *evscal, LAPACK_COMPLEX *cs1, LAPACK_COMPLEX *sn1);
 
/* Subroutine */ int claev2_(LAPACK_COMPLEX *a, LAPACK_COMPLEX *b, LAPACK_COMPLEX *c__, LAPACK_REAL *rt1, 
			     LAPACK_REAL *rt2, LAPACK_REAL *cs1, LAPACK_COMPLEX *sn1);
 
/* Subroutine */ int clags2_(LAPACK_LOGICAL *upper, LAPACK_REAL *a1, LAPACK_COMPLEX *a2, LAPACK_REAL *a3, 
			     LAPACK_REAL *b1, LAPACK_COMPLEX *b2, LAPACK_REAL *b3, LAPACK_REAL *csu, LAPACK_COMPLEX *snu, LAPACK_REAL *csv, 
			     LAPACK_COMPLEX *snv, LAPACK_REAL *csq, LAPACK_COMPLEX *snq);
 
/* Subroutine */ int clagtm_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *
			     alpha, LAPACK_COMPLEX *dl, LAPACK_COMPLEX *d__, LAPACK_COMPLEX *du, LAPACK_COMPLEX *x, LAPACK_INTEGER *
			     ldx, LAPACK_REAL *beta, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int clahef_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_INTEGER *kb,
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *w, LAPACK_INTEGER *ldw, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int clahqr_(LAPACK_LOGICAL *wantt, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_COMPLEX *w, 
			     LAPACK_INTEGER *iloz, LAPACK_INTEGER *ihiz, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int clahrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *nb, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_COMPLEX *y, 
			     LAPACK_INTEGER *ldy);
 
/* Subroutine */ int claic1_(LAPACK_INTEGER *job, LAPACK_INTEGER *j, LAPACK_COMPLEX *x, LAPACK_REAL *sest,
			     LAPACK_COMPLEX *w, LAPACK_COMPLEX *gamma, LAPACK_REAL *sestpr, LAPACK_COMPLEX *s, LAPACK_COMPLEX *c__);
 
/* Subroutine */ int clals0_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *bx, 
			     LAPACK_INTEGER *ldbx, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, 
			     LAPACK_INTEGER *ldgcol, LAPACK_REAL *givnum, LAPACK_INTEGER *ldgnum, LAPACK_REAL *poles, LAPACK_REAL *
			     difl, LAPACK_REAL *difr, LAPACK_REAL *z__, LAPACK_INTEGER *k, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int clalsa_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *bx, LAPACK_INTEGER *ldbx, 
			     LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, LAPACK_INTEGER *k, LAPACK_REAL *difl, LAPACK_REAL *difr, 
			     LAPACK_REAL *z__, LAPACK_REAL *poles, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *
			     ldgcol, LAPACK_INTEGER *perm, LAPACK_REAL *givnum, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int clapll_(LAPACK_INTEGER *n, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx, LAPACK_COMPLEX *
			     y, LAPACK_INTEGER *incy, LAPACK_REAL *ssmin);
 
/* Subroutine */ int clapmt_(LAPACK_LOGICAL *forwrd, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX 
			     *x, LAPACK_INTEGER *ldx, LAPACK_INTEGER *k);
 
/* Subroutine */ int claqgb_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL 
			     *colcnd, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int claqge_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL *colcnd, LAPACK_REAL *amax, char *
			     equed);
 
/* Subroutine */ int claqhb_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int claqhe_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int claqhp_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_REAL *s, 
			     LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int claqp2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, LAPACK_COMPLEX 
			     *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_COMPLEX *tau, LAPACK_REAL *vn1, LAPACK_REAL *vn2, 
			     LAPACK_COMPLEX *work);
 
/* Subroutine */ int claqps_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, LAPACK_INTEGER 
			     *nb, LAPACK_INTEGER *kb, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_COMPLEX *
			     tau, LAPACK_REAL *vn1, LAPACK_REAL *vn2, LAPACK_COMPLEX *auxv, LAPACK_COMPLEX *f, LAPACK_INTEGER *ldf);
 
/* Subroutine */ int claqsb_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int claqsp_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_REAL *s, 
			     LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int claqsy_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int clar1v_(LAPACK_INTEGER *n, LAPACK_INTEGER *b1, LAPACK_INTEGER *bn, LAPACK_REAL *
			     sigma, LAPACK_REAL *d__, LAPACK_REAL *l, LAPACK_REAL *ld, LAPACK_REAL *lld, LAPACK_REAL *gersch, LAPACK_COMPLEX 
			     *z__, LAPACK_REAL *ztz, LAPACK_REAL *mingma, LAPACK_INTEGER *r__, LAPACK_INTEGER *isuppz, LAPACK_REAL *
			     work);
 
/* Subroutine */ int clar2v_(LAPACK_INTEGER *n, LAPACK_COMPLEX *x, LAPACK_COMPLEX *y, LAPACK_COMPLEX *z__,
			     LAPACK_INTEGER *incx, LAPACK_REAL *c__, LAPACK_COMPLEX *s, LAPACK_INTEGER *incc);
 
/* Subroutine */ int clarcm_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *rwork);
 
/* Subroutine */ int clarf_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *v, 
			    LAPACK_INTEGER *incv, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *
			    work);
 
/* Subroutine */ int clarfb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, 
			     LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, 
			     LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int clarfg_(LAPACK_INTEGER *n, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *x, LAPACK_INTEGER *
			     incx, LAPACK_COMPLEX *tau);
 
/* Subroutine */ int clarft_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt);
 
/* Subroutine */ int clarfx_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *v, 
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work);
 
/* Subroutine */ int clargv_(LAPACK_INTEGER *n, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx, LAPACK_COMPLEX *
			     y, LAPACK_INTEGER *incy, LAPACK_REAL *c__, LAPACK_INTEGER *incc);
 
/* Subroutine */ int clarnv_(LAPACK_INTEGER *idist, LAPACK_INTEGER *iseed, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *x);
 
/* Subroutine */ int clarrv_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *l, LAPACK_INTEGER *isplit, 
			     LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_INTEGER *iblock, LAPACK_REAL *gersch, LAPACK_REAL *tol, 
			     LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_REAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int clartg_(LAPACK_COMPLEX *f, LAPACK_COMPLEX *g, LAPACK_REAL *cs, LAPACK_COMPLEX *sn, 
			     LAPACK_COMPLEX *r__);
 
/* Subroutine */ int clartv_(LAPACK_INTEGER *n, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx, LAPACK_COMPLEX *
			     y, LAPACK_INTEGER *incy, LAPACK_REAL *c__, LAPACK_COMPLEX *s, LAPACK_INTEGER *incc);
 
/* Subroutine */ int clarz_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, 
			    LAPACK_COMPLEX *v, LAPACK_INTEGER *incv, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, 
			    LAPACK_COMPLEX *work);
 
/* Subroutine */ int clarzb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_COMPLEX *v, 
			     LAPACK_INTEGER *ldv, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int clarzt_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt);
 
/* Subroutine */ int clascl_(char *type__, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_REAL *
			     cfrom, LAPACK_REAL *cto, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int claset_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *
			     alpha, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int clasr_(char *side, char *pivot, char *direct, LAPACK_INTEGER *m,
			    LAPACK_INTEGER *n, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int classq_(LAPACK_INTEGER *n, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx, LAPACK_REAL *
			     scale, LAPACK_REAL *sumsq);
 
/* Subroutine */ int claswp_(LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *
			     k1, LAPACK_INTEGER *k2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *incx);
 
/* Subroutine */ int clasyf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_INTEGER *kb,
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *w, LAPACK_INTEGER *ldw, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int clatbs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *
			     x, LAPACK_REAL *scale, LAPACK_REAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int clatdf_(LAPACK_INTEGER *ijob, LAPACK_INTEGER *n, LAPACK_COMPLEX *z__, LAPACK_INTEGER 
			     *ldz, LAPACK_COMPLEX *rhs, LAPACK_REAL *rdsum, LAPACK_REAL *rdscal, LAPACK_INTEGER *ipiv, LAPACK_INTEGER 
			     *jpiv);
 
/* Subroutine */ int clatps_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *x, LAPACK_REAL *scale, LAPACK_REAL *cnorm,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int clatrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *e, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *w, LAPACK_INTEGER *ldw);
 
/* Subroutine */ int clatrs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *x, LAPACK_REAL *scale,
			     LAPACK_REAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int clatrz_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work);
 
/* Subroutine */ int clatzm_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *v, 
			     LAPACK_INTEGER *incv, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c1, LAPACK_COMPLEX *c2, LAPACK_INTEGER *ldc, 
			     LAPACK_COMPLEX *work);
 
/* Subroutine */ int clauu2_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int clauum_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *afb, LAPACK_INTEGER *ldafb, 
			     LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *
			     berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbstf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			    nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			    info);
 
/* Subroutine */ int cpbsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *afb, LAPACK_INTEGER *
			     ldafb, char *equed, LAPACK_REAL *s, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, 
			     LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbtf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbtrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab,
			     LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpbtrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cpocon_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpoequ_(LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_REAL *s, 
			     LAPACK_REAL *scond, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int cporfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb,
			     LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cposv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *a,
			    LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cposvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, char *
			     equed, LAPACK_REAL *s, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cpotf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cpotrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cpotri_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cpotrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cppcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_REAL *anorm,
			     LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cppequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_REAL *s, 
			     LAPACK_REAL *scond, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     ap, LAPACK_COMPLEX *afp, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cppsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			    ap, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cppsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *afp, char *equed, LAPACK_REAL *s, LAPACK_COMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL 
			     *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cpptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cpptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     ap, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cptcon_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_COMPLEX *e, LAPACK_REAL *anorm, 
			     LAPACK_REAL *rcond, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cptrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__,
			     LAPACK_COMPLEX *e, LAPACK_REAL *df, LAPACK_COMPLEX *ef, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX 
			     *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cptsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__, LAPACK_COMPLEX *e, 
			    LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cptsvx_(char *fact, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__,
			     LAPACK_COMPLEX *e, LAPACK_REAL *df, LAPACK_COMPLEX *ef, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX 
			     *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpttrf_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_COMPLEX *e, LAPACK_INTEGER *info);
 
/* Subroutine */ int cpttrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__,
			     LAPACK_COMPLEX *e, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cptts2_(LAPACK_INTEGER *iuplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *
			     d__, LAPACK_COMPLEX *e, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int crot_(LAPACK_INTEGER *n, LAPACK_COMPLEX *cx, LAPACK_INTEGER *incx, LAPACK_COMPLEX *
			   cy, LAPACK_INTEGER *incy, LAPACK_REAL *c__, LAPACK_COMPLEX *s);
 
/* Subroutine */ int cspcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cspmv_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *
			    ap, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *y, LAPACK_INTEGER *
			    incy);
 
/* Subroutine */ int cspr_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *x,
			   LAPACK_INTEGER *incx, LAPACK_COMPLEX *ap);
 
/* Subroutine */ int csprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     ap, LAPACK_COMPLEX *afp, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x,
			     LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cspsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			    ap, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int cspsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *afp, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *
			     ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int csptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int csptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int csptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     ap, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int csrot_(LAPACK_INTEGER *n, LAPACK_COMPLEX *cx, LAPACK_INTEGER *incx, LAPACK_COMPLEX *
			    cy, LAPACK_INTEGER *incy, LAPACK_REAL *c__, LAPACK_REAL *s);
 
/* Subroutine */ int csrscl_(LAPACK_INTEGER *n, LAPACK_REAL *sa, LAPACK_COMPLEX *sx, LAPACK_INTEGER *incx);
 
/* Subroutine */ int cstedc_(char *compz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *
			     rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cstein_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_INTEGER *m, LAPACK_REAL 
			     *w, LAPACK_INTEGER *iblock, LAPACK_INTEGER *isplit, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int csteqr_(char *compz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int csycon_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int csymv_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *
			    a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *x, LAPACK_INTEGER *incx, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *y,
			    LAPACK_INTEGER *incy);
 
/* Subroutine */ int csyr_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *x,
			   LAPACK_INTEGER *incx, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int csyrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *
			     b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int csysv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *a,
			    LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *work,
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int csysvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *
			     ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond,
			     LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int csytf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int csytrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int csytri_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int csytrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ctbcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctbrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, 
			     LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctbtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_COMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, 
			     LAPACK_INTEGER *m, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgex2_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *q, 
			     LAPACK_INTEGER *ldq, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *j1, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgexc_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *q, 
			     LAPACK_INTEGER *ldq, LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *ifst, LAPACK_INTEGER *
			     ilst, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgsen_(LAPACK_INTEGER *ijob, LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, 
			     LAPACK_LOGICAL *select, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_COMPLEX *alpha, LAPACK_COMPLEX *beta, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq,
			     LAPACK_COMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *m, LAPACK_REAL *pl, LAPACK_REAL *pr, LAPACK_REAL *
			     dif, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgsja_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_COMPLEX *a, LAPACK_INTEGER *
			     lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_REAL *tola, LAPACK_REAL *tolb, LAPACK_REAL *alpha, 
			     LAPACK_REAL *beta, LAPACK_COMPLEX *u, LAPACK_INTEGER *ldu, LAPACK_COMPLEX *v, LAPACK_INTEGER *ldv, 
			     LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_COMPLEX *work, LAPACK_INTEGER *ncycle, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ctgsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_REAL *s, LAPACK_REAL 
			     *dif, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER 
			     *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgsy2_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *d__, LAPACK_INTEGER *ldd, LAPACK_COMPLEX *e, LAPACK_INTEGER *lde, 
			     LAPACK_COMPLEX *f, LAPACK_INTEGER *ldf, LAPACK_REAL *scale, LAPACK_REAL *rdsum, LAPACK_REAL *rdscal, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctgsyl_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *d__, LAPACK_INTEGER *ldd, LAPACK_COMPLEX *e, LAPACK_INTEGER *lde, 
			     LAPACK_COMPLEX *f, LAPACK_INTEGER *ldf, LAPACK_REAL *scale, LAPACK_REAL *dif, LAPACK_COMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctpcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *ap, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctprfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_COMPLEX *x, 
			     LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctptri_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctptrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_REAL *rcond, LAPACK_COMPLEX *work, LAPACK_REAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, 
			     LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_COMPLEX *work, 
			     LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrexc_(char *compq, LAPACK_INTEGER *n, LAPACK_COMPLEX *t, LAPACK_INTEGER *
			     ldt, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *ifst, LAPACK_INTEGER *ilst, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ctrrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_COMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_COMPLEX *work, LAPACK_REAL 
			     *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrsen_(char *job, char *compq, LAPACK_LOGICAL *select, LAPACK_INTEGER 
			     *n, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_COMPLEX *w, 
			     LAPACK_INTEGER *m, LAPACK_REAL *s, LAPACK_REAL *sep, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_COMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_COMPLEX *vl, LAPACK_INTEGER *ldvl, 
			     LAPACK_COMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_REAL *s, LAPACK_REAL *sep, LAPACK_INTEGER *mm, LAPACK_INTEGER *
			     m, LAPACK_COMPLEX *work, LAPACK_INTEGER *ldwork, LAPACK_REAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrsyl_(char *trana, char *tranb, LAPACK_INTEGER *isgn, LAPACK_INTEGER 
			     *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrti2_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrtri_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctrtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ctzrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int ctzrzf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cung2l_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cung2r_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cungbr_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int cunghr_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *
			     a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER 
			     *info);
 
/* Subroutine */ int cungl2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunglq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cungql_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cungqr_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cungr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cungrq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cungtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda,
			     LAPACK_COMPLEX *tau, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunm2l_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunm2r_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmbr_(char *vect, char *side, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, 
			     LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cunmhr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, 
			     LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cunml2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmlq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmql_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmqr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmr2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmr3_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, 
			     LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmrq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cunmrz_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, 
			     LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int cunmtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_COMPLEX *a, LAPACK_INTEGER *lda, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_COMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int cupgtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *
			     tau, LAPACK_COMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int cupmtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_COMPLEX *ap, LAPACK_COMPLEX *tau, LAPACK_COMPLEX *c__, LAPACK_INTEGER *ldc, 
			     LAPACK_COMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dbdsdc_(char *uplo, char *compq, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *vt, 
			     LAPACK_INTEGER *ldvt, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *iq, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dbdsqr_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ncvt, LAPACK_INTEGER *
			     nru, LAPACK_INTEGER *ncc, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *vt, 
			     LAPACK_INTEGER *ldvt, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *
			     ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int ddisna_(char *job, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *sep, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbbrd_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *ncc,
			     LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *pt, 
			     LAPACK_INTEGER *ldpt, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbcon_(char *norm, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *rowcnd, LAPACK_DOUBLEREAL *colcnd, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dgbrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *afb, 
			     LAPACK_INTEGER *ldafb, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_INTEGER *
			    nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, 
			    LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl,
			     LAPACK_INTEGER *ku, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLEREAL *afb, LAPACK_INTEGER *ldafb, LAPACK_INTEGER *ipiv, char *equed, 
			     LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbtf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbtrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgbtrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgebak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *scale, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *
			     ldv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgebal_(char *job, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgebd2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tauq, LAPACK_DOUBLEREAL *
			     taup, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgebrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tauq, LAPACK_DOUBLEREAL *
			     taup, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgecon_(char *norm, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *rowcnd, LAPACK_DOUBLEREAL 
			     *colcnd, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgees_(char *jobvs, char *sort, LAPACK_L_FP select, LAPACK_INTEGER *n, 
			    LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, LAPACK_DOUBLEREAL *wr, 
			    LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *vs, LAPACK_INTEGER *ldvs, LAPACK_DOUBLEREAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeesx_(char *jobvs, char *sort, LAPACK_L_FP select, char *
			     sense, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, 
			     LAPACK_DOUBLEREAL *wr, LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *vs, LAPACK_INTEGER *ldvs, 
			     LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *rcondv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			    a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *wr, LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *vl, 
			    LAPACK_INTEGER *ldvl, LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLEREAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *wr, 
			     LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLEREAL *vr, 
			     LAPACK_INTEGER *ldvr, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *scale, 
			     LAPACK_DOUBLEREAL *abnrm, LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *rcondv, LAPACK_DOUBLEREAL 
			     *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgegs_(char *jobvsl, char *jobvsr, LAPACK_INTEGER *n, 
			    LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *
			    alphar, LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *vsl, 
			    LAPACK_INTEGER *ldvsl, LAPACK_DOUBLEREAL *vsr, LAPACK_INTEGER *ldvsr, LAPACK_DOUBLEREAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgegv_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			    a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alphar, 
			    LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, 
			    LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dgehd2_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgehrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgelq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgelqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgels_(char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *
			    nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgelsd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *
			     s, LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgelss_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *
			     s, LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgelsx_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     jpvt, LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dgelsy_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     jpvt, LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeql2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeqlf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeqp3_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *jpvt, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeqpf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *jpvt, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeqr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgeqrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgerfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *
			     ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgerq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgerqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgesc2_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLEREAL *rhs, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv, LAPACK_DOUBLEREAL *scale);
 
/* Subroutine */ int dgesdd_(char *jobz, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, 
			     LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *ldvt, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgesv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER 
			    *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *
			     ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *ldvt, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgesvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_INTEGER *ipiv, char *equed, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgetc2_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER 
			     *ipiv, LAPACK_INTEGER *jpiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgetf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgetrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgetri_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER 
			     *ipiv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgetrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dggbak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *lscale, LAPACK_DOUBLEREAL *rscale, LAPACK_INTEGER *m, 
			     LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *ldv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dggbal_(char *job, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLEREAL *lscale, LAPACK_DOUBLEREAL *rscale, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dgges_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			    delctg, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, 
			    LAPACK_INTEGER *ldb, LAPACK_INTEGER *sdim, LAPACK_DOUBLEREAL *alphar, LAPACK_DOUBLEREAL *alphai, 
			    LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *vsl, LAPACK_INTEGER *ldvsl, LAPACK_DOUBLEREAL *vsr, 
			    LAPACK_INTEGER *ldvsr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_LOGICAL *bwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dggesx_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			     delctg, char *sense, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *sdim, LAPACK_DOUBLEREAL *alphar, 
			     LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *vsl, LAPACK_INTEGER *ldvsl,
			     LAPACK_DOUBLEREAL *vsr, LAPACK_INTEGER *ldvsr, LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *
			     rcondv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     liwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dggev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			    a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alphar, 
			    LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, 
			    LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alphar, LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *
			     beta, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *lscale, LAPACK_DOUBLEREAL *rscale, 
			     LAPACK_DOUBLEREAL *abnrm, LAPACK_DOUBLEREAL *bbnrm, LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *
			     rcondv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_LOGICAL *
			     bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dggglm_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *x, LAPACK_DOUBLEREAL *y, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dgghrd_(char *compq, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *
			     ldz, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgglse_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *x, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dggqrf_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *taua, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *taub, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dggrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *taua, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *taub, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dggsvd_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLEREAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alpha, 
			     LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER 
			     *ldv, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dggsvp_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *tola, LAPACK_DOUBLEREAL *tolb, LAPACK_INTEGER *k, LAPACK_INTEGER 
			     *l, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *ldv, 
			     LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *iwork, LAPACK_DOUBLEREAL *tau, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgtcon_(char *norm, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *dl, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *du2, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgtrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *dl, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *dlf, 
			     LAPACK_DOUBLEREAL *df, LAPACK_DOUBLEREAL *duf, LAPACK_DOUBLEREAL *du2, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dgtsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *dl, 
			    LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER 
			    *info);
 
/* Subroutine */ int dgtsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *dl, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *
			     dlf, LAPACK_DOUBLEREAL *df, LAPACK_DOUBLEREAL *duf, LAPACK_DOUBLEREAL *du2, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgttrf_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *dl, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *du2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgttrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *dl, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *du2, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dgtts2_(LAPACK_INTEGER *itrans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *dl, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, LAPACK_DOUBLEREAL *du2, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int dhgeqz_(char *job, char *compq, char *compz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *
			     b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alphar, LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *
			     beta, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dhsein_(char *side, char *eigsrc, char *initv, LAPACK_LOGICAL *
			     select, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLEREAL *wr, 
			     LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLEREAL *vr, 
			     LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     ifaill, LAPACK_INTEGER *ifailr, LAPACK_INTEGER *info);
 
/* Subroutine */ int dhseqr_(char *job, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo,
			     LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLEREAL *wr, 
			     LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlabad_(LAPACK_DOUBLEREAL *small, LAPACK_DOUBLEREAL *large);
 
/* Subroutine */ int dlabrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tauq, 
			     LAPACK_DOUBLEREAL *taup, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *y, LAPACK_INTEGER 
			     *ldy);
 
/* Subroutine */ int dlacon_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *v, LAPACK_DOUBLEREAL *x, 
			     LAPACK_INTEGER *isgn, LAPACK_DOUBLEREAL *est, LAPACK_INTEGER *kase);
 
/* Subroutine */ int dlacpy_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int dladiv_(LAPACK_DOUBLEREAL *a, LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *p, LAPACK_DOUBLEREAL *q);
 
/* Subroutine */ int dlae2_(LAPACK_DOUBLEREAL *a, LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *c__, 
			    LAPACK_DOUBLEREAL *rt1, LAPACK_DOUBLEREAL *rt2);
 
/* Subroutine */ int dlaebz_(LAPACK_INTEGER *ijob, LAPACK_INTEGER *nitmax, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *mmax, LAPACK_INTEGER *minp, LAPACK_INTEGER *nbmin, LAPACK_DOUBLEREAL *abstol, 
			     LAPACK_DOUBLEREAL *reltol, LAPACK_DOUBLEREAL *pivmin, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *
			     e, LAPACK_DOUBLEREAL *e2, LAPACK_INTEGER *nval, LAPACK_DOUBLEREAL *ab, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_INTEGER *mout, LAPACK_INTEGER *nab, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed0_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *qsiz, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, 
			     LAPACK_DOUBLEREAL *qstore, LAPACK_INTEGER *ldqs, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed1_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *q, 
			     LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, LAPACK_DOUBLEREAL *rho, LAPACK_INTEGER *cutpnt, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed2_(LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER *n1, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, LAPACK_DOUBLEREAL *rho, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *dlamda, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *q2, 
			     LAPACK_INTEGER *indx, LAPACK_INTEGER *indxc, LAPACK_INTEGER *indxp, LAPACK_INTEGER *coltyp, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed3_(LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER *n1, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *rho, LAPACK_DOUBLEREAL *dlamda,
			     LAPACK_DOUBLEREAL *q2, LAPACK_INTEGER *indx, LAPACK_INTEGER *ctot, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLEREAL *s, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed4_(LAPACK_INTEGER *n, LAPACK_INTEGER *i__, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *delta, LAPACK_DOUBLEREAL *rho, LAPACK_DOUBLEREAL *dlam,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed5_(LAPACK_INTEGER *i__, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_DOUBLEREAL *delta, LAPACK_DOUBLEREAL *rho, LAPACK_DOUBLEREAL *dlam);
 
/* Subroutine */ int dlaed6_(LAPACK_INTEGER *kniter, LAPACK_LOGICAL *orgati, LAPACK_DOUBLEREAL *
			     rho, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *finit, LAPACK_DOUBLEREAL *
			     tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed7_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *n, LAPACK_INTEGER *qsiz, 
			     LAPACK_INTEGER *tlvls, LAPACK_INTEGER *curlvl, LAPACK_INTEGER *curpbm, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, LAPACK_DOUBLEREAL *rho, LAPACK_INTEGER 
			     *cutpnt, LAPACK_DOUBLEREAL *qstore, LAPACK_INTEGER *qptr, LAPACK_INTEGER *prmptr, LAPACK_INTEGER *
			     perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_DOUBLEREAL *givnum, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed8_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER 
			     *qsiz, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, 
			     LAPACK_DOUBLEREAL *rho, LAPACK_INTEGER *cutpnt, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *dlamda,
			     LAPACK_DOUBLEREAL *q2, LAPACK_INTEGER *ldq2, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *perm, LAPACK_INTEGER 
			     *givptr, LAPACK_INTEGER *givcol, LAPACK_DOUBLEREAL *givnum, LAPACK_INTEGER *indxp, LAPACK_INTEGER 
			     *indx, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaed9_(LAPACK_INTEGER *k, LAPACK_INTEGER *kstart, LAPACK_INTEGER *kstop, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *
			     rho, LAPACK_DOUBLEREAL *dlamda, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *s, LAPACK_INTEGER *lds, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaeda_(LAPACK_INTEGER *n, LAPACK_INTEGER *tlvls, LAPACK_INTEGER *curlvl, 
			     LAPACK_INTEGER *curpbm, LAPACK_INTEGER *prmptr, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, 
			     LAPACK_INTEGER *givcol, LAPACK_DOUBLEREAL *givnum, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *qptr, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *ztemp, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaein_(LAPACK_LOGICAL *rightv, LAPACK_LOGICAL *noinit, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLEREAL *wr, LAPACK_DOUBLEREAL *wi, 
			     LAPACK_DOUBLEREAL *vr, LAPACK_DOUBLEREAL *vi, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *work, LAPACK_DOUBLEREAL *eps3, LAPACK_DOUBLEREAL *smlnum, LAPACK_DOUBLEREAL *
			     bignum, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlaev2_(LAPACK_DOUBLEREAL *a, LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *rt1, LAPACK_DOUBLEREAL *rt2, LAPACK_DOUBLEREAL *cs1, LAPACK_DOUBLEREAL *sn1);
 
/* Subroutine */ int dlaexc_(LAPACK_LOGICAL *wantq, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *t, 
			     LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *j1, LAPACK_INTEGER *n1, 
			     LAPACK_INTEGER *n2, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlag2_(LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, 
			    LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *safmin, LAPACK_DOUBLEREAL *scale1, LAPACK_DOUBLEREAL *
			    scale2, LAPACK_DOUBLEREAL *wr1, LAPACK_DOUBLEREAL *wr2, LAPACK_DOUBLEREAL *wi);
 
/* Subroutine */ int dlags2_(LAPACK_LOGICAL *upper, LAPACK_DOUBLEREAL *a1, LAPACK_DOUBLEREAL *a2, 
			     LAPACK_DOUBLEREAL *a3, LAPACK_DOUBLEREAL *b1, LAPACK_DOUBLEREAL *b2, LAPACK_DOUBLEREAL *b3, 
			     LAPACK_DOUBLEREAL *csu, LAPACK_DOUBLEREAL *snu, LAPACK_DOUBLEREAL *csv, LAPACK_DOUBLEREAL *snv, 
			     LAPACK_DOUBLEREAL *csq, LAPACK_DOUBLEREAL *snq);
 
/* Subroutine */ int dlagtf_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_DOUBLEREAL *lambda, 
			     LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *tol, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_INTEGER *in, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlagtm_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *dl, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *du, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER 
			     *ldb);
 
/* Subroutine */ int dlagts_(LAPACK_INTEGER *job, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, 
			     LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *d__, LAPACK_INTEGER *in, 
			     LAPACK_DOUBLEREAL *y, LAPACK_DOUBLEREAL *tol, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlagv2_(LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alphar, LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *
			     beta, LAPACK_DOUBLEREAL *csl, LAPACK_DOUBLEREAL *snl, LAPACK_DOUBLEREAL *csr, LAPACK_DOUBLEREAL *
			     snr);
 
/* Subroutine */ int dlahqr_(LAPACK_LOGICAL *wantt, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLEREAL 
			     *wr, LAPACK_DOUBLEREAL *wi, LAPACK_INTEGER *iloz, LAPACK_INTEGER *ihiz, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlahrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *nb, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, 
			     LAPACK_DOUBLEREAL *y, LAPACK_INTEGER *ldy);
 
/* Subroutine */ int dlaic1_(LAPACK_INTEGER *job, LAPACK_INTEGER *j, LAPACK_DOUBLEREAL *x, 
			     LAPACK_DOUBLEREAL *sest, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *gamma, LAPACK_DOUBLEREAL *
			     sestpr, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *c__);
 
/* Subroutine */ int dlaln2_(LAPACK_LOGICAL *ltrans, LAPACK_INTEGER *na, LAPACK_INTEGER *nw, 
			     LAPACK_DOUBLEREAL *smin, LAPACK_DOUBLEREAL *ca, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLEREAL *d1, LAPACK_DOUBLEREAL *d2, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *wr, LAPACK_DOUBLEREAL *wi, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *xnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlals0_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL 
			     *bx, LAPACK_INTEGER *ldbx, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, 
			     LAPACK_INTEGER *ldgcol, LAPACK_DOUBLEREAL *givnum, LAPACK_INTEGER *ldgnum, LAPACK_DOUBLEREAL *
			     poles, LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *difr, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *
			     k, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlalsa_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *bx, LAPACK_INTEGER *
			     ldbx, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *difr, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *
			     poles, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_INTEGER *
			     perm, LAPACK_DOUBLEREAL *givnum, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlalsd_(char *uplo, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, LAPACK_INTEGER 
			     *nrhs, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlamc1_(LAPACK_INTEGER *beta, LAPACK_INTEGER *t, LAPACK_LOGICAL *rnd, LAPACK_LOGICAL 
			     *ieee1);
 
/* Subroutine */ int dlamc2_(LAPACK_INTEGER *beta, LAPACK_INTEGER *t, LAPACK_LOGICAL *rnd, 
			     LAPACK_DOUBLEREAL *eps, LAPACK_INTEGER *emin, LAPACK_DOUBLEREAL *rmin, LAPACK_INTEGER *emax, 
			     LAPACK_DOUBLEREAL *rmax);
 
/* Subroutine */ int dlamc4_(LAPACK_INTEGER *emin, LAPACK_DOUBLEREAL *start, LAPACK_INTEGER *base);
 
/* Subroutine */ int dlamc5_(LAPACK_INTEGER *beta, LAPACK_INTEGER *p, LAPACK_INTEGER *emin, 
			     LAPACK_LOGICAL *ieee, LAPACK_INTEGER *emax, LAPACK_DOUBLEREAL *rmax);
 
/* Subroutine */ int dlamrg_(LAPACK_INTEGER *n1, LAPACK_INTEGER *n2, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER 
			     *dtrd1, LAPACK_INTEGER *dtrd2, LAPACK_INTEGER *index);
 
/* Subroutine */ int dlanv2_(LAPACK_DOUBLEREAL *a, LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *rt1r, LAPACK_DOUBLEREAL *rt1i, LAPACK_DOUBLEREAL *rt2r,
			     LAPACK_DOUBLEREAL *rt2i, LAPACK_DOUBLEREAL *cs, LAPACK_DOUBLEREAL *sn);
 
/* Subroutine */ int dlapll_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLEREAL *y, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *ssmin);
 
/* Subroutine */ int dlapmt_(LAPACK_LOGICAL *forwrd, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_INTEGER *k);
 
/* Subroutine */ int dlaqgb_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *rowcnd, LAPACK_DOUBLEREAL *colcnd, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int dlaqge_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *rowcnd, LAPACK_DOUBLEREAL 
			     *colcnd, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int dlaqp2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_DOUBLEREAL *tau, 
			     LAPACK_DOUBLEREAL *vn1, LAPACK_DOUBLEREAL *vn2, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlaqps_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, LAPACK_INTEGER 
			     *nb, LAPACK_INTEGER *kb, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, 
			     LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *vn1, LAPACK_DOUBLEREAL *vn2, LAPACK_DOUBLEREAL *auxv, 
			     LAPACK_DOUBLEREAL *f, LAPACK_INTEGER *ldf);
 
/* Subroutine */ int dlaqsb_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *
			     ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax,
			     char *equed);
 
/* Subroutine */ int dlaqsp_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int dlaqsy_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int dlaqtr_(LAPACK_LOGICAL *ltran, LAPACK_LOGICAL *lLAPACK_REAL, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *b, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL 
			     *scale, LAPACK_DOUBLEREAL *x, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlar1v_(LAPACK_INTEGER *n, LAPACK_INTEGER *b1, LAPACK_INTEGER *bn, LAPACK_DOUBLEREAL 
			     *sigma, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *l, LAPACK_DOUBLEREAL *ld, LAPACK_DOUBLEREAL *
			     lld, LAPACK_DOUBLEREAL *gersch, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *ztz, LAPACK_DOUBLEREAL 
			     *mingma, LAPACK_INTEGER *r__, LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlar2v_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *x, LAPACK_DOUBLEREAL *y, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *incx, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, 
			     LAPACK_INTEGER *incc);
 
/* Subroutine */ int dlarf_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *v,
			    LAPACK_INTEGER *incv, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, 
			    LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *
			     ldv, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int dlarfg_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *x, 
			     LAPACK_INTEGER *incx, LAPACK_DOUBLEREAL *tau);
 
/* Subroutine */ int dlarft_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *t, 
			     LAPACK_INTEGER *ldt);
 
/* Subroutine */ int dlarfx_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     v, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlargv_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLEREAL *y, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *incc);
 
/* Subroutine */ int dlarnv_(LAPACK_INTEGER *idist, LAPACK_INTEGER *iseed, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *x);
 
/* Subroutine */ int dlarrb_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *l, 
			     LAPACK_DOUBLEREAL *ld, LAPACK_DOUBLEREAL *lld, LAPACK_INTEGER *ifirst, LAPACK_INTEGER *ilast, 
			     LAPACK_DOUBLEREAL *sigma, LAPACK_DOUBLEREAL *reltol, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *
			     wgap, LAPACK_DOUBLEREAL *werr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dlarre_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLEREAL *tol, LAPACK_INTEGER *nsplit, LAPACK_INTEGER *isplit, LAPACK_INTEGER *m, 
			     LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *woff, LAPACK_DOUBLEREAL *gersch, LAPACK_DOUBLEREAL *work,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dlarrf_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *l, 
			     LAPACK_DOUBLEREAL *ld, LAPACK_DOUBLEREAL *lld, LAPACK_INTEGER *ifirst, LAPACK_INTEGER *ilast, 
			     LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *dplus, LAPACK_DOUBLEREAL *lplus, LAPACK_DOUBLEREAL *work,
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlarrv_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *l, 
			     LAPACK_INTEGER *isplit, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *iblock, 
			     LAPACK_DOUBLEREAL *gersch, LAPACK_DOUBLEREAL *tol, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlartg_(LAPACK_DOUBLEREAL *f, LAPACK_DOUBLEREAL *g, LAPACK_DOUBLEREAL *cs, 
			     LAPACK_DOUBLEREAL *sn, LAPACK_DOUBLEREAL *r__);
 
/* Subroutine */ int dlartv_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLEREAL *y, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_INTEGER 
			     *incc);
 
/* Subroutine */ int dlaruv_(LAPACK_INTEGER *iseed, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *x);
 
/* Subroutine */ int dlarz_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, 
			    LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *incv, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *c__, 
			    LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlarzb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLEREAL *v,
			     LAPACK_INTEGER *ldv, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *
			     ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int dlarzt_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *t, 
			     LAPACK_INTEGER *ldt);
 
/* Subroutine */ int dlas2_(LAPACK_DOUBLEREAL *f, LAPACK_DOUBLEREAL *g, LAPACK_DOUBLEREAL *h__, 
			    LAPACK_DOUBLEREAL *ssmin, LAPACK_DOUBLEREAL *ssmax);
 
/* Subroutine */ int dlascl_(char *type__, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, 
			     LAPACK_DOUBLEREAL *cfrom, LAPACK_DOUBLEREAL *cto, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd0_(LAPACK_INTEGER *n, LAPACK_INTEGER *sqre, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *
			     ldvt, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *iwork, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dlasd1_(LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, LAPACK_INTEGER *sqre, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *u, 
			     LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *ldvt, LAPACK_INTEGER *idxq, LAPACK_INTEGER *
			     iwork, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd2_(LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, LAPACK_INTEGER *sqre, LAPACK_INTEGER 
			     *k, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *
			     beta, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *ldvt, 
			     LAPACK_DOUBLEREAL *dsigma, LAPACK_DOUBLEREAL *u2, LAPACK_INTEGER *ldu2, LAPACK_DOUBLEREAL *vt2, 
			     LAPACK_INTEGER *ldvt2, LAPACK_INTEGER *idxp, LAPACK_INTEGER *idx, LAPACK_INTEGER *idxc, LAPACK_INTEGER *
			     idxq, LAPACK_INTEGER *coltyp, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd3_(LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, LAPACK_INTEGER *sqre, LAPACK_INTEGER 
			     *k, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *dsigma, 
			     LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *u2, LAPACK_INTEGER *ldu2, 
			     LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *ldvt, LAPACK_DOUBLEREAL *vt2, LAPACK_INTEGER *ldvt2, 
			     LAPACK_INTEGER *idxc, LAPACK_INTEGER *ctot, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd4_(LAPACK_INTEGER *n, LAPACK_INTEGER *i__, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *delta, LAPACK_DOUBLEREAL *rho, LAPACK_DOUBLEREAL *
			     sigma, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd5_(LAPACK_INTEGER *i__, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_DOUBLEREAL *delta, LAPACK_DOUBLEREAL *rho, LAPACK_DOUBLEREAL *dsigma, LAPACK_DOUBLEREAL *
			     work);
 
/* Subroutine */ int dlasd6_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *vf, LAPACK_DOUBLEREAL *vl, 
			     LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *beta, LAPACK_INTEGER *idxq, LAPACK_INTEGER *perm, 
			     LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_DOUBLEREAL *givnum,
			     LAPACK_INTEGER *ldgnum, LAPACK_DOUBLEREAL *poles, LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *
			     difr, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd7_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_DOUBLEREAL *zw, LAPACK_DOUBLEREAL *vf, LAPACK_DOUBLEREAL *vfw, LAPACK_DOUBLEREAL *vl, 
			     LAPACK_DOUBLEREAL *vlw, LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *
			     dsigma, LAPACK_INTEGER *idx, LAPACK_INTEGER *idxp, LAPACK_INTEGER *idxq, LAPACK_INTEGER *perm, 
			     LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_DOUBLEREAL *givnum,
			     LAPACK_INTEGER *ldgnum, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd8_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *vf, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *difl, 
			     LAPACK_DOUBLEREAL *difr, LAPACK_INTEGER *lddifr, LAPACK_DOUBLEREAL *dsigma, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasd9_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *ldu, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *vf, LAPACK_DOUBLEREAL *vl, 
			     LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *difr, LAPACK_DOUBLEREAL *dsigma, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasda_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *sqre, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER 
			     *ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *difr, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *poles, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, 
			     LAPACK_INTEGER *ldgcol, LAPACK_INTEGER *perm, LAPACK_DOUBLEREAL *givnum, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasdq_(char *uplo, LAPACK_INTEGER *sqre, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     ncvt, LAPACK_INTEGER *nru, LAPACK_INTEGER *ncc, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *ldvt, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, 
			     LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasdt_(LAPACK_INTEGER *n, LAPACK_INTEGER *lvl, LAPACK_INTEGER *nd, LAPACK_INTEGER *
			     inode, LAPACK_INTEGER *ndiml, LAPACK_INTEGER *ndimr, LAPACK_INTEGER *msub);
 
/* Subroutine */ int dlaset_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     alpha, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int dlasq1_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasq2_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasq3_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *pp, LAPACK_DOUBLEREAL *dmin__, LAPACK_DOUBLEREAL *sigma, LAPACK_DOUBLEREAL *desig,
			     LAPACK_DOUBLEREAL *qmax, LAPACK_INTEGER *nfail, LAPACK_INTEGER *iter, LAPACK_INTEGER *ndiv, 
			     LAPACK_LOGICAL *ieee);
 
/* Subroutine */ int dlasq4_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *pp, LAPACK_INTEGER *n0in, LAPACK_DOUBLEREAL *dmin__, LAPACK_DOUBLEREAL *dmin1, 
			     LAPACK_DOUBLEREAL *dmin2, LAPACK_DOUBLEREAL *dn, LAPACK_DOUBLEREAL *dn1, LAPACK_DOUBLEREAL *dn2, 
			     LAPACK_DOUBLEREAL *tau, LAPACK_INTEGER *ttype);
 
/* Subroutine */ int dlasq5_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *pp, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *dmin__, LAPACK_DOUBLEREAL *dmin1, 
			     LAPACK_DOUBLEREAL *dmin2, LAPACK_DOUBLEREAL *dn, LAPACK_DOUBLEREAL *dnm1, LAPACK_DOUBLEREAL *dnm2,
			     LAPACK_LOGICAL *ieee);
 
/* Subroutine */ int dlasq6_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *pp, LAPACK_DOUBLEREAL *dmin__, LAPACK_DOUBLEREAL *dmin1, LAPACK_DOUBLEREAL *dmin2,
			     LAPACK_DOUBLEREAL *dn, LAPACK_DOUBLEREAL *dnm1, LAPACK_DOUBLEREAL *dnm2);
 
/* Subroutine */ int dlasr_(char *side, char *pivot, char *direct, LAPACK_INTEGER *m,
			    LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			    lda);
 
/* Subroutine */ int dlasrt_(char *id, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dlassq_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *sumsq);
 
/* Subroutine */ int dlasv2_(LAPACK_DOUBLEREAL *f, LAPACK_DOUBLEREAL *g, LAPACK_DOUBLEREAL *h__, 
			     LAPACK_DOUBLEREAL *ssmin, LAPACK_DOUBLEREAL *ssmax, LAPACK_DOUBLEREAL *snr, LAPACK_DOUBLEREAL *
			     csr, LAPACK_DOUBLEREAL *snl, LAPACK_DOUBLEREAL *csl);
 
/* Subroutine */ int dlaswp_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER 
			     *k1, LAPACK_INTEGER *k2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *incx);
 
/* Subroutine */ int dlasy2_(LAPACK_LOGICAL *ltranl, LAPACK_LOGICAL *ltranr, LAPACK_INTEGER *isgn, 
			     LAPACK_INTEGER *n1, LAPACK_INTEGER *n2, LAPACK_DOUBLEREAL *tl, LAPACK_INTEGER *ldtl, LAPACK_DOUBLEREAL *
			     tr, LAPACK_INTEGER *ldtr, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *scale, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *xnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlasyf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_INTEGER *kb,
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *
			     ldw, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlatbs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLEREAL *x, LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlatdf_(LAPACK_INTEGER *ijob, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *rhs, LAPACK_DOUBLEREAL *rdsum, LAPACK_DOUBLEREAL *rdscal, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv);
 
/* Subroutine */ int dlatps_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *x, LAPACK_DOUBLEREAL *scale, 
			     LAPACK_DOUBLEREAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlatrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *w, 
			     LAPACK_INTEGER *ldw);
 
/* Subroutine */ int dlatrs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *x, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlatrz_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlatzm_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     v, LAPACK_INTEGER *incv, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *c1, LAPACK_DOUBLEREAL *c2, 
			     LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int dlauu2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dlauum_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dopgtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, 
			     LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dopmtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER 
			     *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorg2l_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorg2r_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgbr_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorghr_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgl2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorglq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgql_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgqr_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgrq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dorgtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorm2l_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dorm2r_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormbr_(char *vect, char *side, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, 
			     LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dormhr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *
			     tau, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dorml2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormlq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormql_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormqr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormr2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormr3_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, 
			     LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormrq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dormrz_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, 
			     LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dormtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *
			     ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *
			     ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *afb, LAPACK_INTEGER *ldafb, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dpbstf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *
			     ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			    nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *afb, 
			     LAPACK_INTEGER *ldafb, char *equed, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr,
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbtf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *
			     ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbtrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *
			     ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpbtrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dpocon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpoequ_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int dporfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dposv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL 
			    *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dposvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *af, LAPACK_INTEGER *ldaf, 
			     char *equed, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *
			     x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *
			     berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpotf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpotrf_(const char *uplo, const LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, const LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpotri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpotrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dppcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, 
			     LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dppequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *afp, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dppsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL 
			    *ap, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dppsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *afp, char *equed, LAPACK_DOUBLEREAL *s, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dpptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dpptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dptcon_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dpteqr_(char *compz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dptrfs_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *df, LAPACK_DOUBLEREAL *ef, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER 
			     *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr,
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dptsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *d__, 
			    LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dptsvx_(char *fact, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *df, LAPACK_DOUBLEREAL *ef, 
			     LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *
			     rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dpttrf_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dpttrs_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dptts2_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int drscl_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *sa, LAPACK_DOUBLEREAL *sx, 
			    LAPACK_INTEGER *incx);
 
/* Subroutine */ int dsbev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			    LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			    LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *
			     ldq, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, 
			     LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbgst_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			     LAPACK_INTEGER *kb, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *bb, LAPACK_INTEGER *
			     ldbb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbgv_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			    LAPACK_INTEGER *kb, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *bb, LAPACK_INTEGER *
			    ldbb, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbgvd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			     LAPACK_INTEGER *kb, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *bb, LAPACK_INTEGER *
			     ldbb, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbgvx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ka, LAPACK_INTEGER *kb, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *
			     bb, LAPACK_INTEGER *ldbb, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *vl, 
			     LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER 
			     *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsbtrd_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER 
			     *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspev_(char *jobz, char *uplo, LAPACK_INTEGER *n, const LAPACK_DOUBLEREAL *
			    ap, const LAPACK_DOUBLEREAL *w, const LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dspevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     ap, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *
			     iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dspgst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *bp, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspgv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *bp, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			    LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspgvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *bp, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspgvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *bp, LAPACK_DOUBLEREAL *vl, 
			     LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER 
			     *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *afp, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL 
			    *ap, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dspsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *afp, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, 
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dsptrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *ap, LAPACK_INTEGER *
			     ipiv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dstebz_(char *range, char *order, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL 
			     *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_INTEGER *m, LAPACK_INTEGER *nsplit, 
			     LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *iblock, LAPACK_INTEGER *isplit, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dstedc_(char *compz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dstegr_(char *jobz, char *range, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, 
			     LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dstein_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *iblock, LAPACK_INTEGER *isplit, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsteqr_(char *compz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dsterf_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dstev_(char *jobz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			    LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dstevd_(char *jobz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dstevr_(char *jobz, char *range, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, 
			     LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dstevx_(char *jobz, char *range, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, 
			     LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsycon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsyev_(const char *jobz, const char *uplo, const LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a,
			    const LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *work, const LAPACK_INTEGER *lwork, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int dsyevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsyevr_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *
			     il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsyevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *
			     il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsygs2_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dsygst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dsygv_(const LAPACK_INTEGER *itype, const char *jobz, const char *uplo, const LAPACK_INTEGER *
			    n, LAPACK_DOUBLEREAL *a, const LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, const LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *work, const LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsygvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsygvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER 
			     *ldb, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, 
			     LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLEREAL *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsyrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *
			     ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dsysv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL 
			    *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsysvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *
			     ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsytd2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsytf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsytrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsytrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsytri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dsytrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtbcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *rcond, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtbrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL 
			     *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtbtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL 
			     *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER 
			     *mm, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgex2_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *
			     q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *j1, LAPACK_INTEGER *
			     n1, LAPACK_INTEGER *n2, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgexc_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *
			     q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *ifst, 
			     LAPACK_INTEGER *ilst, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgsen_(LAPACK_INTEGER *ijob, LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, 
			     LAPACK_LOGICAL *select, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *
			     b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alphar, LAPACK_DOUBLEREAL *alphai, LAPACK_DOUBLEREAL *
			     beta, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *pl, LAPACK_DOUBLEREAL *pr, LAPACK_DOUBLEREAL *dif, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgsja_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLEREAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *tola, 
			     LAPACK_DOUBLEREAL *tolb, LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLEREAL *u, 
			     LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *
			     ldq, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *ncycle, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *dif, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgsy2_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *d__, LAPACK_INTEGER *ldd, 
			     LAPACK_DOUBLEREAL *e, LAPACK_INTEGER *lde, LAPACK_DOUBLEREAL *f, LAPACK_INTEGER *ldf, LAPACK_DOUBLEREAL *
			     scale, LAPACK_DOUBLEREAL *rdsum, LAPACK_DOUBLEREAL *rdscal, LAPACK_INTEGER *iwork, LAPACK_INTEGER 
			     *pq, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtgsyl_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *d__, LAPACK_INTEGER *ldd, 
			     LAPACK_DOUBLEREAL *e, LAPACK_INTEGER *lde, LAPACK_DOUBLEREAL *f, LAPACK_INTEGER *ldf, LAPACK_DOUBLEREAL *
			     scale, LAPACK_DOUBLEREAL *dif, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtpcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int dtprfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtptri_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     ap, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtptrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *ap, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int dtrcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *
			     ldvl, LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrexc_(char *compq, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *
			     ldt, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *ifst, LAPACK_INTEGER *ilst, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_DOUBLEREAL *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrsen_(char *job, char *compq, LAPACK_LOGICAL *select, LAPACK_INTEGER 
			     *n, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *q, LAPACK_INTEGER *ldq, 
			     LAPACK_DOUBLEREAL *wr, LAPACK_DOUBLEREAL *wi, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL 
			     *sep, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLEREAL *vl, LAPACK_INTEGER *
			     ldvl, LAPACK_DOUBLEREAL *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *sep, 
			     LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *ldwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrsyl_(char *trana, char *tranb, LAPACK_INTEGER *isgn, LAPACK_INTEGER 
			     *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrti2_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrtri_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtrtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtzrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int dtzrzf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLEREAL *tau, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
LAPACK_INTEGER icmax1_(LAPACK_INTEGER *n, LAPACK_COMPLEX *cx, LAPACK_INTEGER *incx);
 
LAPACK_INTEGER ieeeck_(LAPACK_INTEGER *ispec, LAPACK_REAL *zero, LAPACK_REAL *one);
 
LAPACK_INTEGER ilaenv_(const LAPACK_INTEGER *ispec, char *name__, const char *opts, const LAPACK_INTEGER *n1, 
		       const LAPACK_INTEGER *n2, const LAPACK_INTEGER *n3, const LAPACK_INTEGER *n4, LAPACK_FTNLEN name_len, LAPACK_FTNLEN 
		       opts_len);
 
LAPACK_INTEGER izmax1_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *cx, LAPACK_INTEGER *incx);
 
/* Subroutine */ int sbdsdc_(char *uplo, char *compq, LAPACK_INTEGER *n, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, LAPACK_INTEGER *ldvt, LAPACK_REAL *q, 
			     LAPACK_INTEGER *iq, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sbdsqr_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ncvt, LAPACK_INTEGER *
			     nru, LAPACK_INTEGER *ncc, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *vt, LAPACK_INTEGER *ldvt, LAPACK_REAL *
			     u, LAPACK_INTEGER *ldu, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sdisna_(char *job, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *d__, 
			     LAPACK_REAL *sep, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbbrd_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *ncc,
			     LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *d__, LAPACK_REAL *
			     e, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *pt, LAPACK_INTEGER *ldpt, LAPACK_REAL *c__, LAPACK_INTEGER 
			     *ldc, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbcon_(char *norm, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, 
			     LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL *
			     colcnd, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *afb, LAPACK_INTEGER *ldafb,
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *
			     ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_INTEGER *
			    nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl,
			     LAPACK_INTEGER *ku, LAPACK_INTEGER *nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *afb, 
			     LAPACK_INTEGER *ldafb, LAPACK_INTEGER *ipiv, char *equed, LAPACK_REAL *r__, LAPACK_REAL *c__, 
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr,
			     LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbtf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbtrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgbtrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgebak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_REAL *scale, LAPACK_INTEGER *m, LAPACK_REAL *v, LAPACK_INTEGER *ldv, LAPACK_INTEGER 
			     *info);
 
/* Subroutine */ int sgebal_(char *job, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgebd2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *tauq, LAPACK_REAL *taup, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgebrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *tauq, LAPACK_REAL *taup, LAPACK_REAL *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgecon_(char *norm, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL *colcnd, LAPACK_REAL *amax, LAPACK_INTEGER 
			     *info);
 
/* Subroutine */ int sgees_(char *jobvs, char *sort, LAPACK_L_FP select, LAPACK_INTEGER *n, 
			    LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL *vs, 
			    LAPACK_INTEGER *ldvs, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *
			    info);
 
/* Subroutine */ int sgeesx_(char *jobvs, char *sort, LAPACK_L_FP select, char *
			     sense, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, LAPACK_REAL *wr, 
			     LAPACK_REAL *wi, LAPACK_REAL *vs, LAPACK_INTEGER *ldvs, LAPACK_REAL *rconde, LAPACK_REAL *rcondv, LAPACK_REAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_LOGICAL *bwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, 
			    LAPACK_INTEGER *ldvr, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL *
			     vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *ilo, LAPACK_INTEGER *
			     ihi, LAPACK_REAL *scale, LAPACK_REAL *abnrm, LAPACK_REAL *rconde, LAPACK_REAL *rcondv, LAPACK_REAL *work,
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgegs_(char *jobvsl, char *jobvsr, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL 
			    *beta, LAPACK_REAL *vsl, LAPACK_INTEGER *ldvsl, LAPACK_REAL *vsr, LAPACK_INTEGER *ldvsr, LAPACK_REAL *
			    work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgegv_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL 
			    *beta, LAPACK_REAL *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, LAPACK_INTEGER *ldvr, LAPACK_REAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgehd2_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgehrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgelq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgelqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgels_(char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *
			    nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgelsd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *s, LAPACK_REAL *rcond, LAPACK_INTEGER *
			     rank, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgelss_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *s, LAPACK_REAL *rcond, LAPACK_INTEGER *
			     rank, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgelsx_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *jpvt, LAPACK_REAL *rcond, 
			     LAPACK_INTEGER *rank, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgelsy_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *jpvt, LAPACK_REAL *rcond, 
			     LAPACK_INTEGER *rank, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeql2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeqlf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeqp3_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *jpvt, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeqpf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *jpvt, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeqr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgeqrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgerfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgerq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgerqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgesc2_(LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *rhs, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv, LAPACK_REAL *scale);
 
/* Subroutine */ int sgesdd_(char *jobz, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *s, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, LAPACK_INTEGER *ldvt,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgesv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			    LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgesvd_(char *jobu, char *jobvt, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *s, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, 
			     LAPACK_INTEGER *ldvt, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgesvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, 
			     char *equed, LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, 
			     LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgetc2_(LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv,
			     LAPACK_INTEGER *jpiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgetf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgetrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgetri_(LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgetrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggbak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_REAL *lscale, LAPACK_REAL *rscale, LAPACK_INTEGER *m, LAPACK_REAL *v, 
			     LAPACK_INTEGER *ldv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggbal_(char *job, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *lscale, LAPACK_REAL 
			     *rscale, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgges_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			    selctg, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, 
			    LAPACK_INTEGER *sdim, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL *beta, LAPACK_REAL *vsl, 
			    LAPACK_INTEGER *ldvsl, LAPACK_REAL *vsr, LAPACK_INTEGER *ldvsr, LAPACK_REAL *work, LAPACK_INTEGER *lwork,
			    LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggesx_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			     selctg, char *sense, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *sdim, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL *beta, 
			     LAPACK_REAL *vsl, LAPACK_INTEGER *ldvsl, LAPACK_REAL *vsr, LAPACK_INTEGER *ldvsr, LAPACK_REAL *rconde, 
			     LAPACK_REAL *rcondv, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     liwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL 
			    *beta, LAPACK_REAL *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, LAPACK_INTEGER *ldvr, LAPACK_REAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL 
			     *alphar, LAPACK_REAL *alphai, LAPACK_REAL *beta, LAPACK_REAL *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, 
			     LAPACK_INTEGER *ldvr, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *lscale, LAPACK_REAL *rscale,
			     LAPACK_REAL *abnrm, LAPACK_REAL *bbnrm, LAPACK_REAL *rconde, LAPACK_REAL *rcondv, LAPACK_REAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggglm_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *d__, LAPACK_REAL *x, LAPACK_REAL *y, 
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgghrd_(char *compq, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL 
			     *q, LAPACK_INTEGER *ldq, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgglse_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *c__, LAPACK_REAL *d__, LAPACK_REAL *x, 
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggqrf_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *taua, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *taub, LAPACK_REAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *taua, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *taub, LAPACK_REAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggsvd_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_REAL *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *alpha, LAPACK_REAL *beta, LAPACK_REAL *u, LAPACK_INTEGER *
			     ldu, LAPACK_REAL *v, LAPACK_INTEGER *ldv, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sggsvp_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_REAL *tola, LAPACK_REAL *tolb, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_REAL *u, LAPACK_INTEGER *ldu,
			     LAPACK_REAL *v, LAPACK_INTEGER *ldv, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *iwork, LAPACK_REAL *
			     tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgtcon_(char *norm, LAPACK_INTEGER *n, LAPACK_REAL *dl, LAPACK_REAL *d__, 
			     LAPACK_REAL *du, LAPACK_REAL *du2, LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_REAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgtrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *dl,
			     LAPACK_REAL *d__, LAPACK_REAL *du, LAPACK_REAL *dlf, LAPACK_REAL *df, LAPACK_REAL *duf, LAPACK_REAL *du2, 
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *
			     ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgtsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *dl, LAPACK_REAL *d__, 
			    LAPACK_REAL *du, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgtsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *dl, LAPACK_REAL *d__, LAPACK_REAL *du, LAPACK_REAL *dlf, LAPACK_REAL *df, LAPACK_REAL *duf, 
			     LAPACK_REAL *du2, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *
			     ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int sgttrf_(LAPACK_INTEGER *n, LAPACK_REAL *dl, LAPACK_REAL *d__, LAPACK_REAL *du, LAPACK_REAL *
			     du2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int sgttrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *dl,
			     LAPACK_REAL *d__, LAPACK_REAL *du, LAPACK_REAL *du2, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int sgtts2_(LAPACK_INTEGER *itrans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL 
			     *dl, LAPACK_REAL *d__, LAPACK_REAL *du, LAPACK_REAL *du2, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *
			     ldb);
 
/* Subroutine */ int shgeqz_(char *job, char *compq, char *compz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL *beta, LAPACK_REAL *q, LAPACK_INTEGER *ldq, 
			     LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int shsein_(char *side, char *eigsrc, char *initv, LAPACK_LOGICAL *
			     select, LAPACK_INTEGER *n, LAPACK_REAL *h__, LAPACK_INTEGER *ldh, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL 
			     *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, 
			     LAPACK_REAL *work, LAPACK_INTEGER *ifaill, LAPACK_INTEGER *ifailr, LAPACK_INTEGER *info);
 
/* Subroutine */ int shseqr_(char *job, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo,
			     LAPACK_INTEGER *ihi, LAPACK_REAL *h__, LAPACK_INTEGER *ldh, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL *z__,
			     LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slabad_(LAPACK_REAL *small, LAPACK_REAL *large);
 
/* Subroutine */ int slabrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *tauq, LAPACK_REAL *taup, LAPACK_REAL *x, 
			     LAPACK_INTEGER *ldx, LAPACK_REAL *y, LAPACK_INTEGER *ldy);
 
/* Subroutine */ int slacon_(LAPACK_INTEGER *n, LAPACK_REAL *v, LAPACK_REAL *x, LAPACK_INTEGER *isgn, 
			     LAPACK_REAL *est, LAPACK_INTEGER *kase);
 
/* Subroutine */ int slacpy_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int sladiv_(LAPACK_REAL *a, LAPACK_REAL *b, LAPACK_REAL *c__, LAPACK_REAL *d__, LAPACK_REAL *p, 
			     LAPACK_REAL *q);
 
/* Subroutine */ int slae2_(LAPACK_REAL *a, LAPACK_REAL *b, LAPACK_REAL *c__, LAPACK_REAL *rt1, LAPACK_REAL *rt2);
 
/* Subroutine */ int slaebz_(LAPACK_INTEGER *ijob, LAPACK_INTEGER *nitmax, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *mmax, LAPACK_INTEGER *minp, LAPACK_INTEGER *nbmin, LAPACK_REAL *abstol, LAPACK_REAL *
			     reltol, LAPACK_REAL *pivmin, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *e2, LAPACK_INTEGER *nval, 
			     LAPACK_REAL *ab, LAPACK_REAL *c__, LAPACK_INTEGER *mout, LAPACK_INTEGER *nab, LAPACK_REAL *work, LAPACK_INTEGER 
			     *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed0_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *qsiz, LAPACK_INTEGER *n, LAPACK_REAL 
			     *d__, LAPACK_REAL *e, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *qstore, LAPACK_INTEGER *ldqs, 
			     LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed1_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *q, LAPACK_INTEGER *ldq, 
			     LAPACK_INTEGER *indxq, LAPACK_REAL *rho, LAPACK_INTEGER *cutpnt, LAPACK_REAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed2_(LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER *n1, LAPACK_REAL *d__, 
			     LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, LAPACK_REAL *rho, LAPACK_REAL *z__, LAPACK_REAL *
			     dlamda, LAPACK_REAL *w, LAPACK_REAL *q2, LAPACK_INTEGER *indx, LAPACK_INTEGER *indxc, LAPACK_INTEGER *
			     indxp, LAPACK_INTEGER *coltyp, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed3_(LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER *n1, LAPACK_REAL *d__, 
			     LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *rho, LAPACK_REAL *dlamda, LAPACK_REAL *q2, LAPACK_INTEGER *
			     indx, LAPACK_INTEGER *ctot, LAPACK_REAL *w, LAPACK_REAL *s, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed4_(LAPACK_INTEGER *n, LAPACK_INTEGER *i__, LAPACK_REAL *d__, LAPACK_REAL *z__, 
			     LAPACK_REAL *delta, LAPACK_REAL *rho, LAPACK_REAL *dlam, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed5_(LAPACK_INTEGER *i__, LAPACK_REAL *d__, LAPACK_REAL *z__, LAPACK_REAL *delta, 
			     LAPACK_REAL *rho, LAPACK_REAL *dlam);
 
/* Subroutine */ int slaed6_(LAPACK_INTEGER *kniter, LAPACK_LOGICAL *orgati, LAPACK_REAL *rho, 
			     LAPACK_REAL *d__, LAPACK_REAL *z__, LAPACK_REAL *finit, LAPACK_REAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed7_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *n, LAPACK_INTEGER *qsiz, 
			     LAPACK_INTEGER *tlvls, LAPACK_INTEGER *curlvl, LAPACK_INTEGER *curpbm, LAPACK_REAL *d__, LAPACK_REAL *q, 
			     LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, LAPACK_REAL *rho, LAPACK_INTEGER *cutpnt, LAPACK_REAL *
			     qstore, LAPACK_INTEGER *qptr, LAPACK_INTEGER *prmptr, LAPACK_INTEGER *perm, LAPACK_INTEGER *
			     givptr, LAPACK_INTEGER *givcol, LAPACK_REAL *givnum, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed8_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER 
			     *qsiz, LAPACK_REAL *d__, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *indxq, LAPACK_REAL *rho, 
			     LAPACK_INTEGER *cutpnt, LAPACK_REAL *z__, LAPACK_REAL *dlamda, LAPACK_REAL *q2, LAPACK_INTEGER *ldq2, 
			     LAPACK_REAL *w, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_REAL *
			     givnum, LAPACK_INTEGER *indxp, LAPACK_INTEGER *indx, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaed9_(LAPACK_INTEGER *k, LAPACK_INTEGER *kstart, LAPACK_INTEGER *kstop, 
			     LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *rho, LAPACK_REAL *dlamda,
			     LAPACK_REAL *w, LAPACK_REAL *s, LAPACK_INTEGER *lds, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaeda_(LAPACK_INTEGER *n, LAPACK_INTEGER *tlvls, LAPACK_INTEGER *curlvl, 
			     LAPACK_INTEGER *curpbm, LAPACK_INTEGER *prmptr, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, 
			     LAPACK_INTEGER *givcol, LAPACK_REAL *givnum, LAPACK_REAL *q, LAPACK_INTEGER *qptr, LAPACK_REAL *z__, 
			     LAPACK_REAL *ztemp, LAPACK_INTEGER *info);
 
/* Subroutine */ int slaein_(LAPACK_LOGICAL *rightv, LAPACK_LOGICAL *noinit, LAPACK_INTEGER *n, 
			     LAPACK_REAL *h__, LAPACK_INTEGER *ldh, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL *vr, LAPACK_REAL *vi, LAPACK_REAL 
			     *b, LAPACK_INTEGER *ldb, LAPACK_REAL *work, LAPACK_REAL *eps3, LAPACK_REAL *smlnum, LAPACK_REAL *bignum, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slaev2_(LAPACK_REAL *a, LAPACK_REAL *b, LAPACK_REAL *c__, LAPACK_REAL *rt1, LAPACK_REAL *
			     rt2, LAPACK_REAL *cs1, LAPACK_REAL *sn1);
 
/* Subroutine */ int slaexc_(LAPACK_LOGICAL *wantq, LAPACK_INTEGER *n, LAPACK_REAL *t, LAPACK_INTEGER *
			     ldt, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *j1, LAPACK_INTEGER *n1, LAPACK_INTEGER *n2, 
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slag2_(LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, 
			    LAPACK_REAL *safmin, LAPACK_REAL *scale1, LAPACK_REAL *scale2, LAPACK_REAL *wr1, LAPACK_REAL *wr2, LAPACK_REAL *
			    wi);
 
/* Subroutine */ int slags2_(LAPACK_LOGICAL *upper, LAPACK_REAL *a1, LAPACK_REAL *a2, LAPACK_REAL *a3, 
			     LAPACK_REAL *b1, LAPACK_REAL *b2, LAPACK_REAL *b3, LAPACK_REAL *csu, LAPACK_REAL *snu, LAPACK_REAL *csv, LAPACK_REAL *
			     snv, LAPACK_REAL *csq, LAPACK_REAL *snq);
 
/* Subroutine */ int slagtf_(LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_REAL *lambda, LAPACK_REAL *b, LAPACK_REAL 
			     *c__, LAPACK_REAL *tol, LAPACK_REAL *d__, LAPACK_INTEGER *in, LAPACK_INTEGER *info);
 
/* Subroutine */ int slagtm_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *
			     alpha, LAPACK_REAL *dl, LAPACK_REAL *d__, LAPACK_REAL *du, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *
			     beta, LAPACK_REAL *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int slagts_(LAPACK_INTEGER *job, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_REAL *b, LAPACK_REAL 
			     *c__, LAPACK_REAL *d__, LAPACK_INTEGER *in, LAPACK_REAL *y, LAPACK_REAL *tol, LAPACK_INTEGER *info);
 
/* Subroutine */ int slagv2_(LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, 
			     LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL *beta, LAPACK_REAL *csl, LAPACK_REAL *snl, LAPACK_REAL *
			     csr, LAPACK_REAL *snr);
 
/* Subroutine */ int slahqr_(LAPACK_LOGICAL *wantt, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *h__, LAPACK_INTEGER *ldh, LAPACK_REAL *wr, LAPACK_REAL *
			     wi, LAPACK_INTEGER *iloz, LAPACK_INTEGER *ihiz, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int slahrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *nb, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *t, LAPACK_INTEGER *ldt, LAPACK_REAL *y, LAPACK_INTEGER *ldy);
 
/* Subroutine */ int slaic1_(LAPACK_INTEGER *job, LAPACK_INTEGER *j, LAPACK_REAL *x, LAPACK_REAL *sest, 
			     LAPACK_REAL *w, LAPACK_REAL *gamma, LAPACK_REAL *sestpr, LAPACK_REAL *s, LAPACK_REAL *c__);
 
/* Subroutine */ int slaln2_(LAPACK_LOGICAL *ltrans, LAPACK_INTEGER *na, LAPACK_INTEGER *nw, LAPACK_REAL *
			     smin, LAPACK_REAL *ca, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *d1, LAPACK_REAL *d2, LAPACK_REAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_REAL *wr, LAPACK_REAL *wi, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *scale, 
			     LAPACK_REAL *xnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int slals0_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_INTEGER *nrhs, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *bx, 
			     LAPACK_INTEGER *ldbx, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, 
			     LAPACK_INTEGER *ldgcol, LAPACK_REAL *givnum, LAPACK_INTEGER *ldgnum, LAPACK_REAL *poles, LAPACK_REAL *
			     difl, LAPACK_REAL *difr, LAPACK_REAL *z__, LAPACK_INTEGER *k, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slalsa_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *bx, LAPACK_INTEGER *ldbx, LAPACK_REAL *
			     u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, LAPACK_INTEGER *k, LAPACK_REAL *difl, LAPACK_REAL *difr, LAPACK_REAL *
			     z__, LAPACK_REAL *poles, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, 
			     LAPACK_INTEGER *perm, LAPACK_REAL *givnum, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slalsd_(char *uplo, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, LAPACK_INTEGER 
			     *nrhs, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *rcond, 
			     LAPACK_INTEGER *rank, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slamc1_(LAPACK_INTEGER *beta, LAPACK_INTEGER *t, LAPACK_LOGICAL *rnd, LAPACK_LOGICAL 
			     *ieee1);
 
/* Subroutine */ int slamc2_(LAPACK_INTEGER *beta, LAPACK_INTEGER *t, LAPACK_LOGICAL *rnd, LAPACK_REAL *
			     eps, LAPACK_INTEGER *emin, LAPACK_REAL *rmin, LAPACK_INTEGER *emax, LAPACK_REAL *rmax);
 
/* Subroutine */ int slamc4_(LAPACK_INTEGER *emin, LAPACK_REAL *start, LAPACK_INTEGER *base);
 
/* Subroutine */ int slamc5_(LAPACK_INTEGER *beta, LAPACK_INTEGER *p, LAPACK_INTEGER *emin, 
			     LAPACK_LOGICAL *ieee, LAPACK_INTEGER *emax, LAPACK_REAL *rmax);
 
/* Subroutine */ int slamrg_(LAPACK_INTEGER *n1, LAPACK_INTEGER *n2, LAPACK_REAL *a, LAPACK_INTEGER *
			     strd1, LAPACK_INTEGER *strd2, LAPACK_INTEGER *index);
 
/* Subroutine */ int slanv2_(LAPACK_REAL *a, LAPACK_REAL *b, LAPACK_REAL *c__, LAPACK_REAL *d__, LAPACK_REAL *
			     rt1r, LAPACK_REAL *rt1i, LAPACK_REAL *rt2r, LAPACK_REAL *rt2i, LAPACK_REAL *cs, LAPACK_REAL *sn);
 
/* Subroutine */ int slapll_(LAPACK_INTEGER *n, LAPACK_REAL *x, LAPACK_INTEGER *incx, LAPACK_REAL *y, 
			     LAPACK_INTEGER *incy, LAPACK_REAL *ssmin);
 
/* Subroutine */ int slapmt_(LAPACK_LOGICAL *forwrd, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *x,
			     LAPACK_INTEGER *ldx, LAPACK_INTEGER *k);
 
/* Subroutine */ int slaqgb_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL *
			     colcnd, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int slaqge_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *r__, LAPACK_REAL *c__, LAPACK_REAL *rowcnd, LAPACK_REAL *colcnd, LAPACK_REAL *amax, char *
			     equed);
 
/* Subroutine */ int slaqp2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, LAPACK_REAL *a,
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_REAL *tau, LAPACK_REAL *vn1, LAPACK_REAL *vn2, LAPACK_REAL *
			     work);
 
/* Subroutine */ int slaqps_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, LAPACK_INTEGER 
			     *nb, LAPACK_INTEGER *kb, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_REAL *tau, 
			     LAPACK_REAL *vn1, LAPACK_REAL *vn2, LAPACK_REAL *auxv, LAPACK_REAL *f, LAPACK_INTEGER *ldf);
 
/* Subroutine */ int slaqsb_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, 
			     LAPACK_INTEGER *ldab, LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int slaqsp_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *s, LAPACK_REAL *
			     scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int slaqsy_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, char *equed);
 
/* Subroutine */ int slaqtr_(LAPACK_LOGICAL *ltran, LAPACK_LOGICAL *lLAPACK_REAL, LAPACK_INTEGER *n, LAPACK_REAL 
			     *t, LAPACK_INTEGER *ldt, LAPACK_REAL *b, LAPACK_REAL *w, LAPACK_REAL *scale, LAPACK_REAL *x, LAPACK_REAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slar1v_(LAPACK_INTEGER *n, LAPACK_INTEGER *b1, LAPACK_INTEGER *bn, LAPACK_REAL *
			     sigma, LAPACK_REAL *d__, LAPACK_REAL *l, LAPACK_REAL *ld, LAPACK_REAL *lld, LAPACK_REAL *gersch, LAPACK_REAL *
			     z__, LAPACK_REAL *ztz, LAPACK_REAL *mingma, LAPACK_INTEGER *r__, LAPACK_INTEGER *isuppz, LAPACK_REAL *
			     work);
 
/* Subroutine */ int slar2v_(LAPACK_INTEGER *n, LAPACK_REAL *x, LAPACK_REAL *y, LAPACK_REAL *z__, LAPACK_INTEGER 
			     *incx, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_INTEGER *incc);
 
/* Subroutine */ int slarf_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *v, 
			    LAPACK_INTEGER *incv, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *work);
 
/* Subroutine */ int slarfb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *v, LAPACK_INTEGER *ldv, 
			     LAPACK_REAL *t, LAPACK_INTEGER *ldt, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *work, LAPACK_INTEGER *
			     ldwork);
 
/* Subroutine */ int slarfg_(LAPACK_INTEGER *n, LAPACK_REAL *alpha, LAPACK_REAL *x, LAPACK_INTEGER *incx, 
			     LAPACK_REAL *tau);
 
/* Subroutine */ int slarft_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_REAL *v, LAPACK_INTEGER *ldv, LAPACK_REAL *tau, LAPACK_REAL *t, LAPACK_INTEGER *ldt);
 
/* Subroutine */ int slarfx_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *v, 
			     LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *work);
 
/* Subroutine */ int slargv_(LAPACK_INTEGER *n, LAPACK_REAL *x, LAPACK_INTEGER *incx, LAPACK_REAL *y, 
			     LAPACK_INTEGER *incy, LAPACK_REAL *c__, LAPACK_INTEGER *incc);
 
/* Subroutine */ int slarnv_(LAPACK_INTEGER *idist, LAPACK_INTEGER *iseed, LAPACK_INTEGER *n, LAPACK_REAL 
			     *x);
 
/* Subroutine */ int slarrb_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *l, LAPACK_REAL *ld, LAPACK_REAL *
			     lld, LAPACK_INTEGER *ifirst, LAPACK_INTEGER *ilast, LAPACK_REAL *sigma, LAPACK_REAL *reltol, LAPACK_REAL 
			     *w, LAPACK_REAL *wgap, LAPACK_REAL *werr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slarre_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *tol, 
			     LAPACK_INTEGER *nsplit, LAPACK_INTEGER *isplit, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *woff, 
			     LAPACK_REAL *gersch, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slarrf_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *l, LAPACK_REAL *ld, LAPACK_REAL *
			     lld, LAPACK_INTEGER *ifirst, LAPACK_INTEGER *ilast, LAPACK_REAL *w, LAPACK_REAL *dplus, LAPACK_REAL *
			     lplus, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slarrv_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *l, LAPACK_INTEGER *isplit, 
			     LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_INTEGER *iblock, LAPACK_REAL *gersch, LAPACK_REAL *tol, LAPACK_REAL *
			     z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slartg_(LAPACK_REAL *f, LAPACK_REAL *g, LAPACK_REAL *cs, LAPACK_REAL *sn, LAPACK_REAL *r__);
 
/* Subroutine */ int slartv_(LAPACK_INTEGER *n, LAPACK_REAL *x, LAPACK_INTEGER *incx, LAPACK_REAL *y, 
			     LAPACK_INTEGER *incy, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_INTEGER *incc);
 
/* Subroutine */ int slaruv_(LAPACK_INTEGER *iseed, LAPACK_INTEGER *n, LAPACK_REAL *x);
 
/* Subroutine */ int slarz_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, 
			    LAPACK_REAL *v, LAPACK_INTEGER *incv, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *
			    work);
 
/* Subroutine */ int slarzb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_REAL *v, 
			     LAPACK_INTEGER *ldv, LAPACK_REAL *t, LAPACK_INTEGER *ldt, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *
			     work, LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int slarzt_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_REAL *v, LAPACK_INTEGER *ldv, LAPACK_REAL *tau, LAPACK_REAL *t, LAPACK_INTEGER *ldt);
 
/* Subroutine */ int slas2_(LAPACK_REAL *f, LAPACK_REAL *g, LAPACK_REAL *h__, LAPACK_REAL *ssmin, LAPACK_REAL *
			    ssmax);
 
/* Subroutine */ int slascl_(char *type__, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_REAL *
			     cfrom, LAPACK_REAL *cto, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd0_(LAPACK_INTEGER *n, LAPACK_INTEGER *sqre, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, LAPACK_INTEGER *ldvt, LAPACK_INTEGER *smlsiz, 
			     LAPACK_INTEGER *iwork, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd1_(LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, LAPACK_INTEGER *sqre, LAPACK_REAL *
			     d__, LAPACK_REAL *alpha, LAPACK_REAL *beta, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, 
			     LAPACK_INTEGER *ldvt, LAPACK_INTEGER *idxq, LAPACK_INTEGER *iwork, LAPACK_REAL *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int slasd2_(LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, LAPACK_INTEGER *sqre, LAPACK_INTEGER 
			     *k, LAPACK_REAL *d__, LAPACK_REAL *z__, LAPACK_REAL *alpha, LAPACK_REAL *beta, LAPACK_REAL *u, LAPACK_INTEGER *
			     ldu, LAPACK_REAL *vt, LAPACK_INTEGER *ldvt, LAPACK_REAL *dsigma, LAPACK_REAL *u2, LAPACK_INTEGER *ldu2, 
			     LAPACK_REAL *vt2, LAPACK_INTEGER *ldvt2, LAPACK_INTEGER *idxp, LAPACK_INTEGER *idx, LAPACK_INTEGER *idxc,
			     LAPACK_INTEGER *idxq, LAPACK_INTEGER *coltyp, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd3_(LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, LAPACK_INTEGER *sqre, LAPACK_INTEGER 
			     *k, LAPACK_REAL *d__, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *dsigma, LAPACK_REAL *u, LAPACK_INTEGER *
			     ldu, LAPACK_REAL *u2, LAPACK_INTEGER *ldu2, LAPACK_REAL *vt, LAPACK_INTEGER *ldvt, LAPACK_REAL *vt2, 
			     LAPACK_INTEGER *ldvt2, LAPACK_INTEGER *idxc, LAPACK_INTEGER *ctot, LAPACK_REAL *z__, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int slasd4_(LAPACK_INTEGER *n, LAPACK_INTEGER *i__, LAPACK_REAL *d__, LAPACK_REAL *z__, 
			     LAPACK_REAL *delta, LAPACK_REAL *rho, LAPACK_REAL *sigma, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd5_(LAPACK_INTEGER *i__, LAPACK_REAL *d__, LAPACK_REAL *z__, LAPACK_REAL *delta, 
			     LAPACK_REAL *rho, LAPACK_REAL *dsigma, LAPACK_REAL *work);
 
/* Subroutine */ int slasd6_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_REAL *d__, LAPACK_REAL *vf, LAPACK_REAL *vl, LAPACK_REAL *alpha, LAPACK_REAL *beta,
			     LAPACK_INTEGER *idxq, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, 
			     LAPACK_INTEGER *ldgcol, LAPACK_REAL *givnum, LAPACK_INTEGER *ldgnum, LAPACK_REAL *poles, LAPACK_REAL *
			     difl, LAPACK_REAL *difr, LAPACK_REAL *z__, LAPACK_INTEGER *k, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd7_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_INTEGER *k, LAPACK_REAL *d__, LAPACK_REAL *z__, LAPACK_REAL *zw, LAPACK_REAL *vf, 
			     LAPACK_REAL *vfw, LAPACK_REAL *vl, LAPACK_REAL *vlw, LAPACK_REAL *alpha, LAPACK_REAL *beta, LAPACK_REAL *dsigma,
			     LAPACK_INTEGER *idx, LAPACK_INTEGER *idxp, LAPACK_INTEGER *idxq, LAPACK_INTEGER *perm, LAPACK_INTEGER *
			     givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_REAL *givnum, LAPACK_INTEGER *
			     ldgnum, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd8_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *k, LAPACK_REAL *d__, LAPACK_REAL *
			     z__, LAPACK_REAL *vf, LAPACK_REAL *vl, LAPACK_REAL *difl, LAPACK_REAL *difr, LAPACK_INTEGER *lddifr, 
			     LAPACK_REAL *dsigma, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasd9_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *ldu, LAPACK_INTEGER *k, LAPACK_REAL *
			     d__, LAPACK_REAL *z__, LAPACK_REAL *vf, LAPACK_REAL *vl, LAPACK_REAL *difl, LAPACK_REAL *difr, LAPACK_REAL *
			     dsigma, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasda_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *sqre, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *vt, 
			     LAPACK_INTEGER *k, LAPACK_REAL *difl, LAPACK_REAL *difr, LAPACK_REAL *z__, LAPACK_REAL *poles, LAPACK_INTEGER *
			     givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_INTEGER *perm, LAPACK_REAL *givnum,
			     LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasdq_(char *uplo, LAPACK_INTEGER *sqre, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     ncvt, LAPACK_INTEGER *nru, LAPACK_INTEGER *ncc, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *vt, 
			     LAPACK_INTEGER *ldvt, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasdt_(LAPACK_INTEGER *n, LAPACK_INTEGER *lvl, LAPACK_INTEGER *nd, LAPACK_INTEGER *
			     inode, LAPACK_INTEGER *ndiml, LAPACK_INTEGER *ndimr, LAPACK_INTEGER *msub);
 
/* Subroutine */ int slaset_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *alpha, 
			     LAPACK_REAL *beta, LAPACK_REAL *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int slasq1_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slasq2_(LAPACK_INTEGER *n, LAPACK_REAL *z__, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasq3_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_REAL *z__, LAPACK_INTEGER *pp,
			     LAPACK_REAL *dmin__, LAPACK_REAL *sigma, LAPACK_REAL *desig, LAPACK_REAL *qmax, LAPACK_INTEGER *nfail, 
			     LAPACK_INTEGER *iter, LAPACK_INTEGER *ndiv, LAPACK_LOGICAL *ieee);
 
/* Subroutine */ int slasq4_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_REAL *z__, LAPACK_INTEGER *pp,
			     LAPACK_INTEGER *n0in, LAPACK_REAL *dmin__, LAPACK_REAL *dmin1, LAPACK_REAL *dmin2, LAPACK_REAL *dn, 
			     LAPACK_REAL *dn1, LAPACK_REAL *dn2, LAPACK_REAL *tau, LAPACK_INTEGER *ttype);
 
/* Subroutine */ int slasq5_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_REAL *z__, LAPACK_INTEGER *pp,
			     LAPACK_REAL *tau, LAPACK_REAL *dmin__, LAPACK_REAL *dmin1, LAPACK_REAL *dmin2, LAPACK_REAL *dn, LAPACK_REAL *
			     dnm1, LAPACK_REAL *dnm2, LAPACK_LOGICAL *ieee);
 
/* Subroutine */ int slasq6_(LAPACK_INTEGER *i0, LAPACK_INTEGER *n0, LAPACK_REAL *z__, LAPACK_INTEGER *pp,
			     LAPACK_REAL *dmin__, LAPACK_REAL *dmin1, LAPACK_REAL *dmin2, LAPACK_REAL *dn, LAPACK_REAL *dnm1, LAPACK_REAL *
			     dnm2);
 
/* Subroutine */ int slasr_(char *side, char *pivot, char *direct, LAPACK_INTEGER *m,
			    LAPACK_INTEGER *n, LAPACK_REAL *c__, LAPACK_REAL *s, LAPACK_REAL *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int slasrt_(char *id, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_INTEGER *info);
 
/* Subroutine */ int slassq_(LAPACK_INTEGER *n, LAPACK_REAL *x, LAPACK_INTEGER *incx, LAPACK_REAL *scale, 
			     LAPACK_REAL *sumsq);
 
/* Subroutine */ int slasv2_(LAPACK_REAL *f, LAPACK_REAL *g, LAPACK_REAL *h__, LAPACK_REAL *ssmin, LAPACK_REAL *
			     ssmax, LAPACK_REAL *snr, LAPACK_REAL *csr, LAPACK_REAL *snl, LAPACK_REAL *csl);
 
/* Subroutine */ int slaswp_(LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *k1, 
			     LAPACK_INTEGER *k2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *incx);
 
/* Subroutine */ int slasy2_(LAPACK_LOGICAL *ltranl, LAPACK_LOGICAL *ltranr, LAPACK_INTEGER *isgn, 
			     LAPACK_INTEGER *n1, LAPACK_INTEGER *n2, LAPACK_REAL *tl, LAPACK_INTEGER *ldtl, LAPACK_REAL *tr, LAPACK_INTEGER *
			     ldtr, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *scale, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL 
			     *xnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int slasyf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_INTEGER *kb,
			     LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_REAL *w, LAPACK_INTEGER *ldw, LAPACK_INTEGER 
			     *info);
 
/* Subroutine */ int slatbs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *x, 
			     LAPACK_REAL *scale, LAPACK_REAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int slatdf_(LAPACK_INTEGER *ijob, LAPACK_INTEGER *n, LAPACK_REAL *z__, LAPACK_INTEGER *
			     ldz, LAPACK_REAL *rhs, LAPACK_REAL *rdsum, LAPACK_REAL *rdscal, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *
			     jpiv);
 
/* Subroutine */ int slatps_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *x, LAPACK_REAL *scale, LAPACK_REAL *cnorm, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slatrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *e, LAPACK_REAL *tau, LAPACK_REAL *w, LAPACK_INTEGER *ldw);
 
/* Subroutine */ int slatrs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *x, LAPACK_REAL *scale, LAPACK_REAL 
			     *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int slatrz_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work);
 
/* Subroutine */ int slatzm_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *v, 
			     LAPACK_INTEGER *incv, LAPACK_REAL *tau, LAPACK_REAL *c1, LAPACK_REAL *c2, LAPACK_INTEGER *ldc, LAPACK_REAL *
			     work);
 
/* Subroutine */ int slauu2_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int slauum_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int sopgtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *tau, 
			     LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sopmtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc, LAPACK_REAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int sorg2l_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorg2r_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgbr_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER 
			     *info);
 
/* Subroutine */ int sorghr_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgl2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorglq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgql_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgqr_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgrq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorgtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorm2l_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorm2r_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormbr_(char *vect, char *side, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormhr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sorml2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormlq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormql_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormqr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormr2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormr3_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormrq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormrz_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sormtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *tau, LAPACK_REAL *c__, LAPACK_INTEGER *ldc,
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, 
			     LAPACK_INTEGER *ldab, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int spbequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, 
			     LAPACK_INTEGER *ldab, LAPACK_REAL *s, LAPACK_REAL *scond, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *afb, LAPACK_INTEGER *ldafb, LAPACK_REAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbstf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, 
			     LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			    nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_INTEGER *nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *afb, LAPACK_INTEGER *ldafb, 
			     char *equed, LAPACK_REAL *s, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, 
			     LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int spbtf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, 
			     LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbtrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_REAL *ab, 
			     LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int spbtrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int spocon_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int spoequ_(LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *s, LAPACK_REAL 
			     *scond, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int sporfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *af, LAPACK_INTEGER *ldaf, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x,
			     LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int sposv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sposvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *af, LAPACK_INTEGER *ldaf, char *equed, 
			     LAPACK_REAL *s, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, 
			     LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int spotf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int spotrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int spotri_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int spotrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sppcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *anorm, 
			     LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sppequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *s, LAPACK_REAL *
			     scond, LAPACK_REAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int spprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, 
			     LAPACK_REAL *afp, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, 
			     LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sppsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, 
			    LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sppsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *ap, LAPACK_REAL *afp, char *equed, LAPACK_REAL *s, LAPACK_REAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL 
			     *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int spptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_INTEGER *info);
 
/* Subroutine */ int spptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_INTEGER *info);
 
/* Subroutine */ int spptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, 
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sptcon_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *anorm, 
			     LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int spteqr_(char *compz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sptrfs_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *df, LAPACK_REAL *ef, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, 
			     LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sptsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__, LAPACK_REAL *e, 
			    LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sptsvx_(char *fact, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__,
			     LAPACK_REAL *e, LAPACK_REAL *df, LAPACK_REAL *ef, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER 
			     *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int spttrf_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_INTEGER *info);
 
/* Subroutine */ int spttrs_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sptts2_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int srscl_(LAPACK_INTEGER *n, LAPACK_REAL *sa, LAPACK_REAL *sx, LAPACK_INTEGER *incx);
 
/* Subroutine */ int ssbev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			    LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work,
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work,
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *vl,
			     LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *
			     w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbgst_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			     LAPACK_INTEGER *kb, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *bb, LAPACK_INTEGER *ldbb, LAPACK_REAL *
			     x, LAPACK_INTEGER *ldx, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbgv_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			    LAPACK_INTEGER *kb, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *bb, LAPACK_INTEGER *ldbb, LAPACK_REAL *
			    w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbgvd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			     LAPACK_INTEGER *kb, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *bb, LAPACK_INTEGER *ldbb, LAPACK_REAL *
			     w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbgvx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ka, LAPACK_INTEGER *kb, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *bb, LAPACK_INTEGER *
			     ldbb, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER 
			     *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL 
			     *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssbtrd_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *q, LAPACK_INTEGER *ldq, 
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_INTEGER *ipiv, 
			     LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, 
			    LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, 
			     LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER 
			     *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_REAL *ap, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, 
			     LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspgst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap,
			     LAPACK_REAL *bp, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspgv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_REAL *ap, LAPACK_REAL *bp, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int sspgvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_REAL *ap, LAPACK_REAL *bp, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspgvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *bp, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il,
			     LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *
			     ldz, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, 
			     LAPACK_REAL *afp, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *
			     ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int sspsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, 
			    LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sspsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *ap, LAPACK_REAL *afp, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL 
			     *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssptrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_REAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_INTEGER *ipiv, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ssptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *ap, LAPACK_INTEGER *ipiv, 
			     LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstebz_(char *range, char *order, LAPACK_INTEGER *n, LAPACK_REAL *vl, 
			     LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_INTEGER *m, LAPACK_INTEGER *nsplit, LAPACK_REAL *w, LAPACK_INTEGER *iblock, LAPACK_INTEGER *
			     isplit, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstedc_(char *compz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstegr_(char *jobz, char *range, LAPACK_INTEGER *n, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, 
			     LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_REAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstein_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_INTEGER *m, LAPACK_REAL 
			     *w, LAPACK_INTEGER *iblock, LAPACK_INTEGER *isplit, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssteqr_(char *compz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, 
			     LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssterf_(LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstev_(char *jobz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *
			    z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstevd_(char *jobz, LAPACK_INTEGER *n, LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL 
			     *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstevr_(char *jobz, char *range, LAPACK_INTEGER *n, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, 
			     LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_REAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int sstevx_(char *jobz, char *range, LAPACK_INTEGER *n, LAPACK_REAL *d__, 
			     LAPACK_REAL *e, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, 
			     LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssycon_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *anorm, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ssyev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_REAL *w, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssyevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *w, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssyevr_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, 
			     LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *
			     isuppz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ssyevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, 
			     LAPACK_REAL *abstol, LAPACK_INTEGER *m, LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssygs2_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssygst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssygv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *w, LAPACK_REAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssygvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *w, LAPACK_REAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssygvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *
			     vl, LAPACK_REAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_REAL *abstol, LAPACK_INTEGER *m, 
			     LAPACK_REAL *w, LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER 
			     *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssyrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_REAL *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, 
			     LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *
			     work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssysv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			    LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *work, 
			    LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssysvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *af, LAPACK_INTEGER *ldaf, LAPACK_INTEGER *ipiv, 
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *rcond, LAPACK_REAL *ferr,
			     LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ssytd2_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssytf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssytrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *d__, LAPACK_REAL *e, LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ssytrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssytri_(char *uplo, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int ssytrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int stbcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *rcond, LAPACK_REAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stbrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *b, LAPACK_INTEGER 
			     *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stbtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_REAL *ab, LAPACK_INTEGER *ldab, LAPACK_REAL *b, LAPACK_INTEGER 
			     *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int stgevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *vl, 
			     LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_REAL 
			     *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int stgex2_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, LAPACK_REAL 
			     *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *
			     z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *j1, LAPACK_INTEGER *n1, LAPACK_INTEGER *n2, LAPACK_REAL *work, 
			     LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stgexc_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, LAPACK_REAL 
			     *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *
			     z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *ifst, LAPACK_INTEGER *ilst, LAPACK_REAL *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stgsen_(LAPACK_INTEGER *ijob, LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, 
			     LAPACK_LOGICAL *select, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *
			     ldb, LAPACK_REAL *alphar, LAPACK_REAL *alphai, LAPACK_REAL *beta, LAPACK_REAL *q, LAPACK_INTEGER *ldq, 
			     LAPACK_REAL *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *m, LAPACK_REAL *pl, LAPACK_REAL *pr, LAPACK_REAL *dif, 
			     LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int stgsja_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_REAL *a, LAPACK_INTEGER *lda,
			     LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *tola, LAPACK_REAL *tolb, LAPACK_REAL *alpha, LAPACK_REAL *
			     beta, LAPACK_REAL *u, LAPACK_INTEGER *ldu, LAPACK_REAL *v, LAPACK_INTEGER *ldv, LAPACK_REAL *q, LAPACK_INTEGER *
			     ldq, LAPACK_REAL *work, LAPACK_INTEGER *ncycle, LAPACK_INTEGER *info);
 
/* Subroutine */ int stgsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *vl, 
			     LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, LAPACK_INTEGER *ldvr, LAPACK_REAL *s, LAPACK_REAL *dif, LAPACK_INTEGER *
			     mm, LAPACK_INTEGER *m, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int stgsy2_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *c__, LAPACK_INTEGER *
			     ldc, LAPACK_REAL *d__, LAPACK_INTEGER *ldd, LAPACK_REAL *e, LAPACK_INTEGER *lde, LAPACK_REAL *f, LAPACK_INTEGER 
			     *ldf, LAPACK_REAL *scale, LAPACK_REAL *rdsum, LAPACK_REAL *rdscal, LAPACK_INTEGER *iwork, LAPACK_INTEGER 
			     *pq, LAPACK_INTEGER *info);
 
/* Subroutine */ int stgsyl_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *c__, LAPACK_INTEGER *
			     ldc, LAPACK_REAL *d__, LAPACK_INTEGER *ldd, LAPACK_REAL *e, LAPACK_INTEGER *lde, LAPACK_REAL *f, LAPACK_INTEGER 
			     *ldf, LAPACK_REAL *scale, LAPACK_REAL *dif, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stpcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_REAL *ap, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stprfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, LAPACK_INTEGER *ldx,
			     LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int stptri_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_REAL *ap, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int stptrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_REAL *ap, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int strcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *rcond, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int strevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_REAL *t, LAPACK_INTEGER *ldt, LAPACK_REAL *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, 
			     LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_REAL *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int strexc_(char *compq, LAPACK_INTEGER *n, LAPACK_REAL *t, LAPACK_INTEGER *ldt, 
			     LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *ifst, LAPACK_INTEGER *ilst, LAPACK_REAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int strrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *x, 
			     LAPACK_INTEGER *ldx, LAPACK_REAL *ferr, LAPACK_REAL *berr, LAPACK_REAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int strsen_(char *job, char *compq, LAPACK_LOGICAL *select, LAPACK_INTEGER 
			     *n, LAPACK_REAL *t, LAPACK_INTEGER *ldt, LAPACK_REAL *q, LAPACK_INTEGER *ldq, LAPACK_REAL *wr, LAPACK_REAL *wi, 
			     LAPACK_INTEGER *m, LAPACK_REAL *s, LAPACK_REAL *sep, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int strsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_REAL *t, LAPACK_INTEGER *ldt, LAPACK_REAL *vl, LAPACK_INTEGER *ldvl, LAPACK_REAL *vr, 
			     LAPACK_INTEGER *ldvr, LAPACK_REAL *s, LAPACK_REAL *sep, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_REAL *
			     work, LAPACK_INTEGER *ldwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int strsyl_(char *trana, char *tranb, LAPACK_INTEGER *isgn, LAPACK_INTEGER 
			     *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_REAL *
			     c__, LAPACK_INTEGER *ldc, LAPACK_REAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int strti2_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int strtri_(char *uplo, char *diag, LAPACK_INTEGER *n, LAPACK_REAL *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int strtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_REAL *a, LAPACK_INTEGER *lda, LAPACK_REAL *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int stzrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int stzrzf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_REAL *a, LAPACK_INTEGER *lda, 
			     LAPACK_REAL *tau, LAPACK_REAL *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int xerbla_(char *srname, LAPACK_INTEGER *info);
 
/* Subroutine */ int zbdsqr_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ncvt, LAPACK_INTEGER *
			     nru, LAPACK_INTEGER *ncc, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *vt, 
			     LAPACK_INTEGER *ldvt, LAPACK_DOUBLECOMPLEX *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLECOMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zdrot_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *cx, LAPACK_INTEGER *incx, 
			    LAPACK_DOUBLECOMPLEX *cy, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s);
 
/* Subroutine */ int zdrscl_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *sa, LAPACK_DOUBLECOMPLEX *sx, 
			     LAPACK_INTEGER *incx);
 
/* Subroutine */ int zgbbrd_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *ncc,
			     LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, 
			     LAPACK_DOUBLECOMPLEX *pt, LAPACK_INTEGER *ldpt, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgbcon_(char *norm, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgbequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *rowcnd, LAPACK_DOUBLEREAL *colcnd, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgbrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *
			     afb, LAPACK_INTEGER *ldafb, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgbsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, LAPACK_INTEGER *
			    nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *
			    b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgbsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl,
			     LAPACK_INTEGER *ku, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLECOMPLEX *afb, LAPACK_INTEGER *ldafb, LAPACK_INTEGER *ipiv, char *equed, 
			     LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgbtf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgbtrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgbtrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *
			     ku, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgebak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *scale, LAPACK_INTEGER *m, LAPACK_DOUBLECOMPLEX *v, 
			     LAPACK_INTEGER *ldv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgebal_(char *job, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER 
			     *lda, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *scale, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgebd2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *tauq, 
			     LAPACK_DOUBLECOMPLEX *taup, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgebrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *tauq, 
			     LAPACK_DOUBLECOMPLEX *taup, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgecon_(char *norm, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeequ_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *rowcnd, 
			     LAPACK_DOUBLEREAL *colcnd, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgees_(char *jobvs, char *sort, LAPACK_L_FP select, LAPACK_INTEGER *n, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, LAPACK_DOUBLECOMPLEX *w, 
			    LAPACK_DOUBLECOMPLEX *vs, LAPACK_INTEGER *ldvs, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			    LAPACK_DOUBLEREAL *rwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeesx_(char *jobvs, char *sort, LAPACK_L_FP select, char *
			     sense, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *sdim, 
			     LAPACK_DOUBLECOMPLEX *w, LAPACK_DOUBLECOMPLEX *vs, LAPACK_INTEGER *ldvs, LAPACK_DOUBLEREAL *
			     rconde, LAPACK_DOUBLEREAL *rcondv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_DOUBLEREAL *rwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, 
			    const LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, const LAPACK_DOUBLECOMPLEX *w, const LAPACK_DOUBLECOMPLEX *vl, 
			    LAPACK_INTEGER *ldvl, const LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, const LAPACK_DOUBLECOMPLEX *work, 
			    LAPACK_INTEGER *lwork, const LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *w, 
			     LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *abnrm, 
			     LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *rcondv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgegs_(char *jobvsl, char *jobvsr, LAPACK_INTEGER *n, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *beta, LAPACK_DOUBLECOMPLEX *vsl, 
			    LAPACK_INTEGER *ldvsl, LAPACK_DOUBLECOMPLEX *vsr, LAPACK_INTEGER *ldvsr, LAPACK_DOUBLECOMPLEX *
			    work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgegv_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *beta, LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER 
			    *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER 
			    *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgehd2_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgehrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgelq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgelqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zgels_(char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *
			    nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgelsx_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *jpvt, LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgelsy_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *jpvt, LAPACK_DOUBLEREAL *rcond, LAPACK_INTEGER *rank, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeql2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeqlf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeqp3_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeqpf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeqr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgeqrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zgerfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work,
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgerq2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgerqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zgesc2_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLECOMPLEX *rhs, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv, LAPACK_DOUBLEREAL *scale);
 
/* Subroutine */ int zgesv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *a, 
			    LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *
			    info);
 
/* Subroutine */ int zgesvd_(char *jobu, char *jobvt, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLECOMPLEX *u, LAPACK_INTEGER *
			     ldu, LAPACK_DOUBLECOMPLEX *vt, LAPACK_INTEGER *ldvt, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgesvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *
			     ldaf, LAPACK_INTEGER *ipiv, char *equed, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgetc2_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgetf2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgetrf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgetri_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgetrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggbak_(char *job, char *side, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, 
			     LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *lscale, LAPACK_DOUBLEREAL *rscale, LAPACK_INTEGER *m, 
			     LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggbal_(char *job, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER 
			     *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLEREAL *lscale, LAPACK_DOUBLEREAL *rscale, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgges_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			    delctg, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			    LAPACK_INTEGER *ldb, LAPACK_INTEGER *sdim, LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *
			    beta, LAPACK_DOUBLECOMPLEX *vsl, LAPACK_INTEGER *ldvsl, LAPACK_DOUBLECOMPLEX *vsr, LAPACK_INTEGER 
			    *ldvsr, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, 
			    LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggesx_(char *jobvsl, char *jobvsr, char *sort, LAPACK_L_FP 
			     delctg, char *sense, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *sdim, LAPACK_DOUBLECOMPLEX *alpha, 
			     LAPACK_DOUBLECOMPLEX *beta, LAPACK_DOUBLECOMPLEX *vsl, LAPACK_INTEGER *ldvsl, 
			     LAPACK_DOUBLECOMPLEX *vsr, LAPACK_INTEGER *ldvsr, LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *
			     rcondv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, 
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggev_(char *jobvl, char *jobvr, LAPACK_INTEGER *n, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *beta, LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER 
			    *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER 
			    *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
			     sense, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *beta, 
			     LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLEREAL *lscale, LAPACK_DOUBLEREAL *rscale, 
			     LAPACK_DOUBLEREAL *abnrm, LAPACK_DOUBLEREAL *bbnrm, LAPACK_DOUBLEREAL *rconde, LAPACK_DOUBLEREAL *
			     rcondv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, 
			     LAPACK_INTEGER *iwork, LAPACK_LOGICAL *bwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggglm_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *x, LAPACK_DOUBLECOMPLEX *y, LAPACK_DOUBLECOMPLEX 
			     *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgghrd_(char *compq, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgglse_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *p, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggqrf_(LAPACK_INTEGER *n, LAPACK_INTEGER *m, LAPACK_INTEGER *p, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *taua, LAPACK_DOUBLECOMPLEX *b,
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *taub, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *p, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *taua, LAPACK_DOUBLECOMPLEX *b,
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *taub, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggsvd_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *p, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *alpha, 
			     LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLECOMPLEX *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLECOMPLEX *v, 
			     LAPACK_INTEGER *ldv, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zggsvp_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX 
			     *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *tola, LAPACK_DOUBLEREAL *tolb, LAPACK_INTEGER *k, 
			     LAPACK_INTEGER *l, LAPACK_DOUBLECOMPLEX *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER 
			     *ldv, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *iwork, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgtcon_(char *norm, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *dl, 
			     LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *du, LAPACK_DOUBLECOMPLEX *du2, LAPACK_INTEGER *
			     ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zgtrfs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *dl, LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *du, 
			     LAPACK_DOUBLECOMPLEX *dlf, LAPACK_DOUBLECOMPLEX *df, LAPACK_DOUBLECOMPLEX *duf, 
			     LAPACK_DOUBLECOMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zgtsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *dl, 
			    LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *du, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb,
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int zgtsvx_(char *fact, char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *dl, LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *du, 
			     LAPACK_DOUBLECOMPLEX *dlf, LAPACK_DOUBLECOMPLEX *df, LAPACK_DOUBLECOMPLEX *duf, 
			     LAPACK_DOUBLECOMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgttrf_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *dl, LAPACK_DOUBLECOMPLEX *
			     d__, LAPACK_DOUBLECOMPLEX *du, LAPACK_DOUBLECOMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zgttrs_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *dl, LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *du, 
			     LAPACK_DOUBLECOMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zgtts2_(LAPACK_INTEGER *itrans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *dl, LAPACK_DOUBLECOMPLEX *d__, LAPACK_DOUBLECOMPLEX *du, 
			     LAPACK_DOUBLECOMPLEX *du2, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int zhbev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			    LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, 
			    LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhbevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, 
			     LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhbevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *q, 
			     LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *
			     iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__,
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork,
			     LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhbgst_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			     LAPACK_INTEGER *kb, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *bb, 
			     LAPACK_INTEGER *ldbb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhbgv_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *ka, 
			    LAPACK_INTEGER *kb, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *bb, 
			    LAPACK_INTEGER *ldbb, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, 
			    LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhbgvx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ka, LAPACK_INTEGER *kb, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLECOMPLEX *bb, LAPACK_INTEGER *ldbb, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, 
			     LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *
			     abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhbtrd_(char *vect, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhecon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zheev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX 
			    *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, 
			    LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zheevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zheevr_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, 
			     LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *
			     w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_INTEGER *isuppz, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *
			     iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zheevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, 
			     LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *
			     w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zhegs2_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zhegst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zhegv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork,
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int zhegvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork,
			     LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhegvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *
			     iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__,
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork,
			     LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int zherfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work,
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhesv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, 
			    LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhesvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *
			     ldaf, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x,
			     LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhetf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhetrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhetrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zhetri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhetrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhgeqz_(char *job, char *compq, char *compz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *
			     beta, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *
			     ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zhpcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhpev_(char *jobz, char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX 
			    *ap, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *
			    work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhpevd_(char *jobz, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhpevx_(char *jobz, char *range, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLEREAL *vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, 
			     LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, 
			     LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhpgst_(LAPACK_INTEGER *itype, char *uplo, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *bp, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhpgv_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			    n, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *bp, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX 
			    *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			    info);
 
/* Subroutine */ int zhpgvd_(LAPACK_INTEGER *itype, char *jobz, char *uplo, LAPACK_INTEGER *
			     n, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *bp, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX 
			     *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zhpgvx_(LAPACK_INTEGER *itype, char *jobz, char *range, char *
			     uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *bp, LAPACK_DOUBLEREAL *
			     vl, LAPACK_DOUBLEREAL *vu, LAPACK_INTEGER *il, LAPACK_INTEGER *iu, LAPACK_DOUBLEREAL *abstol, 
			     LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *
			     ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *afp, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *
			     b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zhpsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			    LAPACK_DOUBLECOMPLEX *ap, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int zhpsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *afp, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhptrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zhsein_(char *side, char *eigsrc, char *initv, LAPACK_LOGICAL *
			     select, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLECOMPLEX *
			     w, LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr,
			     LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, 
			     LAPACK_INTEGER *ifaill, LAPACK_INTEGER *ifailr, LAPACK_INTEGER *info);
 
/* Subroutine */ int zhseqr_(char *job, char *compz, LAPACK_INTEGER *n, LAPACK_INTEGER *ilo,
			     LAPACK_INTEGER *ihi, LAPACK_DOUBLECOMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLECOMPLEX *w, 
			     LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zlabrd_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_DOUBLECOMPLEX *tauq, LAPACK_DOUBLECOMPLEX *taup, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *
			     ldx, LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *ldy);
 
/* Subroutine */ int zlacgv_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx);
 
/* Subroutine */ int zlacon_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *v, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_DOUBLEREAL *est, LAPACK_INTEGER *kase);
 
/* Subroutine */ int zlacp2_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *
			     a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int zlacpy_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int zlacrm_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *rwork);
 
/* Subroutine */ int zlacrt_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *cx, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLECOMPLEX *cy, LAPACK_INTEGER *incy, LAPACK_DOUBLECOMPLEX *c__, LAPACK_DOUBLECOMPLEX *
			     s);
 
/* Subroutine */ int zlaed0_(LAPACK_INTEGER *qsiz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *qstore, 
			     LAPACK_INTEGER *ldqs, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlaed7_(LAPACK_INTEGER *n, LAPACK_INTEGER *cutpnt, LAPACK_INTEGER *qsiz, 
			     LAPACK_INTEGER *tlvls, LAPACK_INTEGER *curlvl, LAPACK_INTEGER *curpbm, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *rho, LAPACK_INTEGER *indxq, 
			     LAPACK_DOUBLEREAL *qstore, LAPACK_INTEGER *qptr, LAPACK_INTEGER *prmptr, LAPACK_INTEGER *perm, 
			     LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_DOUBLEREAL *givnum, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlaed8_(LAPACK_INTEGER *k, LAPACK_INTEGER *n, LAPACK_INTEGER *qsiz, 
			     LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *rho, 
			     LAPACK_INTEGER *cutpnt, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *dlamda, LAPACK_DOUBLECOMPLEX *
			     q2, LAPACK_INTEGER *ldq2, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *indxp, LAPACK_INTEGER *indx, 
			     LAPACK_INTEGER *indxq, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, 
			     LAPACK_DOUBLEREAL *givnum, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlaein_(LAPACK_LOGICAL *rightv, LAPACK_LOGICAL *noinit, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *h__, LAPACK_INTEGER *ldh, LAPACK_DOUBLECOMPLEX *w, LAPACK_DOUBLECOMPLEX *v, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *rwork, LAPACK_DOUBLEREAL *eps3, 
			     LAPACK_DOUBLEREAL *smlnum, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlaesy_(LAPACK_DOUBLECOMPLEX *a, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_DOUBLECOMPLEX *rt1, LAPACK_DOUBLECOMPLEX *rt2, 
			     LAPACK_DOUBLECOMPLEX *evscal, LAPACK_DOUBLECOMPLEX *cs1, LAPACK_DOUBLECOMPLEX *sn1);
 
/* Subroutine */ int zlaev2_(LAPACK_DOUBLECOMPLEX *a, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_DOUBLEREAL *rt1, LAPACK_DOUBLEREAL *rt2, LAPACK_DOUBLEREAL *cs1,
			     LAPACK_DOUBLECOMPLEX *sn1);
 
/* Subroutine */ int zlags2_(LAPACK_LOGICAL *upper, LAPACK_DOUBLEREAL *a1, LAPACK_DOUBLECOMPLEX *
			     a2, LAPACK_DOUBLEREAL *a3, LAPACK_DOUBLEREAL *b1, LAPACK_DOUBLECOMPLEX *b2, LAPACK_DOUBLEREAL *b3,
			     LAPACK_DOUBLEREAL *csu, LAPACK_DOUBLECOMPLEX *snu, LAPACK_DOUBLEREAL *csv, LAPACK_DOUBLECOMPLEX *
			     snv, LAPACK_DOUBLEREAL *csq, LAPACK_DOUBLECOMPLEX *snq);
 
/* Subroutine */ int zlagtm_(char *trans, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLECOMPLEX *dl, LAPACK_DOUBLECOMPLEX *d__, 
			     LAPACK_DOUBLECOMPLEX *du, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *beta, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int zlahef_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_INTEGER *kb,
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *w, 
			     LAPACK_INTEGER *ldw, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlahqr_(LAPACK_LOGICAL *wantt, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLECOMPLEX *h__, LAPACK_INTEGER *ldh, 
			     LAPACK_DOUBLECOMPLEX *w, LAPACK_INTEGER *iloz, LAPACK_INTEGER *ihiz, LAPACK_DOUBLECOMPLEX *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlahrd_(LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *nb, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *t, 
			     LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *ldy);
 
/* Subroutine */ int zlaic1_(LAPACK_INTEGER *job, LAPACK_INTEGER *j, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_DOUBLEREAL *sest, LAPACK_DOUBLECOMPLEX *w, LAPACK_DOUBLECOMPLEX *gamma, LAPACK_DOUBLEREAL *
			     sestpr, LAPACK_DOUBLECOMPLEX *s, LAPACK_DOUBLECOMPLEX *c__);
 
/* Subroutine */ int zlals0_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *nl, LAPACK_INTEGER *nr, 
			     LAPACK_INTEGER *sqre, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *bx, LAPACK_INTEGER *ldbx, LAPACK_INTEGER *perm, LAPACK_INTEGER *givptr, 
			     LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_DOUBLEREAL *givnum, LAPACK_INTEGER *ldgnum,
			     LAPACK_DOUBLEREAL *poles, LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *difr, LAPACK_DOUBLEREAL *
			     z__, LAPACK_INTEGER *k, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zlalsa_(LAPACK_INTEGER *icompq, LAPACK_INTEGER *smlsiz, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *bx, 
			     LAPACK_INTEGER *ldbx, LAPACK_DOUBLEREAL *u, LAPACK_INTEGER *ldu, LAPACK_DOUBLEREAL *vt, LAPACK_INTEGER *
			     k, LAPACK_DOUBLEREAL *difl, LAPACK_DOUBLEREAL *difr, LAPACK_DOUBLEREAL *z__, LAPACK_DOUBLEREAL *
			     poles, LAPACK_INTEGER *givptr, LAPACK_INTEGER *givcol, LAPACK_INTEGER *ldgcol, LAPACK_INTEGER *
			     perm, LAPACK_DOUBLEREAL *givnum, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlapll_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *ssmin);
 
/* Subroutine */ int zlapmt_(LAPACK_LOGICAL *forwrd, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_INTEGER *k);
 
/* Subroutine */ int zlaqgb_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku,
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, 
			     LAPACK_DOUBLEREAL *rowcnd, LAPACK_DOUBLEREAL *colcnd, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int zlaqge_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *r__, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *rowcnd, 
			     LAPACK_DOUBLEREAL *colcnd, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int zlaqhb_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, 
			     LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int zlaqhe_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, 
			     char *equed);
 
/* Subroutine */ int zlaqhp_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int zlaqp2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLEREAL *vn1, LAPACK_DOUBLEREAL *vn2, LAPACK_DOUBLECOMPLEX *work);
 
/* Subroutine */ int zlaqps_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *offset, LAPACK_INTEGER 
			     *nb, LAPACK_INTEGER *kb, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *jpvt, 
			     LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLEREAL *vn1, LAPACK_DOUBLEREAL *vn2, LAPACK_DOUBLECOMPLEX *
			     auxv, LAPACK_DOUBLECOMPLEX *f, LAPACK_INTEGER *ldf);
 
/* Subroutine */ int zlaqsb_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, 
			     LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int zlaqsp_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, char *equed);
 
/* Subroutine */ int zlaqsy_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, 
			     char *equed);
 
/* Subroutine */ int zlar1v_(LAPACK_INTEGER *n, LAPACK_INTEGER *b1, LAPACK_INTEGER *bn, LAPACK_DOUBLEREAL 
			     *sigma, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *l, LAPACK_DOUBLEREAL *ld, LAPACK_DOUBLEREAL *
			     lld, LAPACK_DOUBLEREAL *gersch, LAPACK_DOUBLECOMPLEX *z__, LAPACK_DOUBLEREAL *ztz, 
			     LAPACK_DOUBLEREAL *mingma, LAPACK_INTEGER *r__, LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work);
 
/* Subroutine */ int zlar2v_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *x, LAPACK_DOUBLECOMPLEX *y, 
			     LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *incx, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLECOMPLEX *s, 
			     LAPACK_INTEGER *incc);
 
/* Subroutine */ int zlarcm_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *a, LAPACK_INTEGER *
			     lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc,
			     LAPACK_DOUBLEREAL *rwork);
 
/* Subroutine */ int zlarf_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX 
			    *v, LAPACK_INTEGER *incv, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *
			    ldc, LAPACK_DOUBLECOMPLEX *work);
 
/* Subroutine */ int zlarfb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER 
			     *ldv, LAPACK_DOUBLECOMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *
			     ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int zlarfg_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *
			     x, LAPACK_INTEGER *incx, LAPACK_DOUBLECOMPLEX *tau);
 
/* Subroutine */ int zlarft_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     t, LAPACK_INTEGER *ldt);
 
/* Subroutine */ int zlarfx_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *v, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *
			     ldc, LAPACK_DOUBLECOMPLEX *work);
 
/* Subroutine */ int zlargv_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *c__, LAPACK_INTEGER *incc);
 
/* Subroutine */ int zlarnv_(LAPACK_INTEGER *idist, LAPACK_INTEGER *iseed, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *x);
 
/* Subroutine */ int zlarrv_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *l, 
			     LAPACK_INTEGER *isplit, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *iblock, 
			     LAPACK_DOUBLEREAL *gersch, LAPACK_DOUBLEREAL *tol, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz,
			     LAPACK_INTEGER *isuppz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlartg_(LAPACK_DOUBLECOMPLEX *f, LAPACK_DOUBLECOMPLEX *g, LAPACK_DOUBLEREAL *
			     cs, LAPACK_DOUBLECOMPLEX *sn, LAPACK_DOUBLECOMPLEX *r__);
 
/* Subroutine */ int zlartv_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLECOMPLEX *s, 
			     LAPACK_INTEGER *incc);
 
/* Subroutine */ int zlarz_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, 
			    LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER *incv, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			    c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work);
 
/* Subroutine */ int zlarzb_(char *side, char *trans, char *direct, char *
			     storev, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLECOMPLEX 
			     *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLECOMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *c__, 
			     LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *ldwork);
 
/* Subroutine */ int zlarzt_(char *direct, char *storev, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     k, LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     t, LAPACK_INTEGER *ldt);
 
/* Subroutine */ int zlascl_(char *type__, LAPACK_INTEGER *kl, LAPACK_INTEGER *ku, 
			     LAPACK_DOUBLEREAL *cfrom, LAPACK_DOUBLEREAL *cto, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlaset_(char *uplo, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *beta, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *
			     lda);
 
/* Subroutine */ int zlasr_(char *side, char *pivot, char *direct, LAPACK_INTEGER *m,
			    LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLECOMPLEX *a, 
			    LAPACK_INTEGER *lda);
 
/* Subroutine */ int zlassq_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *sumsq);
 
/* Subroutine */ int zlaswp_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_INTEGER *k1, LAPACK_INTEGER *k2, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *incx);
 
/* Subroutine */ int zlasyf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, LAPACK_INTEGER *kb,
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *w, 
			     LAPACK_INTEGER *ldw, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlatbs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlatdf_(LAPACK_INTEGER *ijob, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *z__, 
			     LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *rhs, LAPACK_DOUBLEREAL *rdsum, LAPACK_DOUBLEREAL *
			     rdscal, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *jpiv);
 
/* Subroutine */ int zlatps_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *x, LAPACK_DOUBLEREAL *
			     scale, LAPACK_DOUBLEREAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlatrd_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nb, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *w, LAPACK_INTEGER *ldw);
 
/* Subroutine */ int zlatrs_(char *uplo, char *trans, char *diag, char *
			     normin, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *cnorm, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlatrz_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *l, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work);
 
/* Subroutine */ int zlatzm_(char *side, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER *incv, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     c1, LAPACK_DOUBLECOMPLEX *c2, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work);
 
/* Subroutine */ int zlauu2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int zlauum_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *
			     rcond, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, 
			     LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *afb, LAPACK_INTEGER *
			     ldafb, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx,
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbstf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			    nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *
			    ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *afb, 
			     LAPACK_INTEGER *ldafb, char *equed, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER 
			     *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *
			     ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbtf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbtrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, 
			     LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpbtrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *kd, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *
			     ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpocon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpoequ_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int zporfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zposv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int zposvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *
			     ldaf, char *equed, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zpotf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpotrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpotri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpotrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zppcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL 
			     *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zppequ_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *scond, LAPACK_DOUBLEREAL *amax, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *afp, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb,
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zppsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			    LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zppsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *afp, char *equed, LAPACK_DOUBLEREAL *
			     s, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zpptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zpptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zptcon_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLECOMPLEX *e, 
			     LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zptrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLECOMPLEX *e, LAPACK_DOUBLEREAL *df, LAPACK_DOUBLECOMPLEX *ef, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zptsv_(LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, LAPACK_DOUBLEREAL *d__, 
			    LAPACK_DOUBLECOMPLEX *e, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int zptsvx_(char *fact, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLECOMPLEX *e, LAPACK_DOUBLEREAL *df, LAPACK_DOUBLECOMPLEX *ef, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zpttrf_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLECOMPLEX *e, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zpttrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLECOMPLEX *e, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zptts2_(LAPACK_INTEGER *iuplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLECOMPLEX *e, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb);
 
/* Subroutine */ int zrot_(LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *cx, LAPACK_INTEGER *incx, 
			   LAPACK_DOUBLECOMPLEX *cy, LAPACK_INTEGER *incy, LAPACK_DOUBLEREAL *c__, LAPACK_DOUBLECOMPLEX *s);
 
/* Subroutine */ int zspcon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zspmv_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *alpha, 
			    LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, LAPACK_DOUBLECOMPLEX *
			    beta, LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *incy);
 
/* Subroutine */ int zspr_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *alpha, 
			   LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, LAPACK_DOUBLECOMPLEX *ap);
 
/* Subroutine */ int zsprfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *afp, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *
			     b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zspsv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			    LAPACK_DOUBLECOMPLEX *ap, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			    LAPACK_INTEGER *info);
 
/* Subroutine */ int zspsvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *afp, LAPACK_INTEGER *ipiv, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsptrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsptri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsptrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zstedc_(char *compz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *lrwork, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *liwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zstein_(LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, LAPACK_DOUBLEREAL *e, 
			     LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *w, LAPACK_INTEGER *iblock, LAPACK_INTEGER *isplit, 
			     LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, LAPACK_INTEGER *iwork, 
			     LAPACK_INTEGER *ifail, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsteqr_(char *compz, LAPACK_INTEGER *n, LAPACK_DOUBLEREAL *d__, 
			     LAPACK_DOUBLEREAL *e, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, LAPACK_DOUBLEREAL *work, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zsycon_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLEREAL *anorm, LAPACK_DOUBLEREAL *rcond, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsymv_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *alpha, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, 
			    LAPACK_DOUBLECOMPLEX *beta, LAPACK_DOUBLECOMPLEX *y, LAPACK_INTEGER *incy);
 
/* Subroutine */ int zsyr_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *alpha, 
			   LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *incx, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda);
 
/* Subroutine */ int zsyrfs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *ldaf, 
			     LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, 
			     LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work,
			     LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsysv_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			    LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, 
			    LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsysvx_(char *fact, char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *
			     nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *af, LAPACK_INTEGER *
			     ldaf, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x,
			     LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsytf2_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsytrf_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zsytri_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zsytrs_(char *uplo, LAPACK_INTEGER *n, LAPACK_INTEGER *nrhs, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *ipiv, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztbcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, LAPACK_DOUBLEREAL *rcond, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztbrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, 
			     LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *
			     rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztbtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *kd, LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ab, LAPACK_INTEGER *ldab, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER 
			     *ldb, LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *
			     ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgex2_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_INTEGER *j1, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgexc_(LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *ldz, 
			     LAPACK_INTEGER *ifst, LAPACK_INTEGER *ilst, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgsen_(LAPACK_INTEGER *ijob, LAPACK_LOGICAL *wantq, LAPACK_LOGICAL *wantz, 
			     LAPACK_LOGICAL *select, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *alpha, LAPACK_DOUBLECOMPLEX *
			     beta, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *z__, LAPACK_INTEGER *
			     ldz, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *pl, LAPACK_DOUBLEREAL *pr, LAPACK_DOUBLEREAL *dif, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *liwork, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgsja_(char *jobu, char *jobv, char *jobq, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *p, LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, LAPACK_DOUBLEREAL *tola, 
			     LAPACK_DOUBLEREAL *tolb, LAPACK_DOUBLEREAL *alpha, LAPACK_DOUBLEREAL *beta, LAPACK_DOUBLECOMPLEX *
			     u, LAPACK_INTEGER *ldu, LAPACK_DOUBLECOMPLEX *v, LAPACK_INTEGER *ldv, LAPACK_DOUBLECOMPLEX *q, 
			     LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *ncycle, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER 
			     *ldb, LAPACK_DOUBLECOMPLEX *vl, LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *
			     ldvr, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *dif, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztgsy2_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *d__, LAPACK_INTEGER *ldd, 
			     LAPACK_DOUBLECOMPLEX *e, LAPACK_INTEGER *lde, LAPACK_DOUBLECOMPLEX *f, LAPACK_INTEGER *ldf, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *rdsum, LAPACK_DOUBLEREAL *rdscal, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ztgsyl_(char *trans, LAPACK_INTEGER *ijob, LAPACK_INTEGER *m, LAPACK_INTEGER *
			     n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *d__, LAPACK_INTEGER *ldd, 
			     LAPACK_DOUBLECOMPLEX *e, LAPACK_INTEGER *lde, LAPACK_DOUBLECOMPLEX *f, LAPACK_INTEGER *ldf, 
			     LAPACK_DOUBLEREAL *scale, LAPACK_DOUBLEREAL *dif, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *iwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztpcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL 
			     *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztprfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, LAPACK_DOUBLEREAL *berr, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztptri_(char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *ap, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztptrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *b, LAPACK_INTEGER *ldb, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrcon_(char *norm, char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLEREAL *rcond, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrevc_(char *side, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *vl, 
			     LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_INTEGER *mm, LAPACK_INTEGER 
			     *m, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrexc_(char *compq, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *t, 
			     LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_INTEGER *ifst, LAPACK_INTEGER *
			     ilst, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrrfs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *x, LAPACK_INTEGER *ldx, LAPACK_DOUBLEREAL *ferr, 
			     LAPACK_DOUBLEREAL *berr, LAPACK_DOUBLECOMPLEX *work, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int ztrsen_(char *job, char *compq, LAPACK_LOGICAL *select, LAPACK_INTEGER 
			     *n, LAPACK_DOUBLECOMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, 
			     LAPACK_DOUBLECOMPLEX *w, LAPACK_INTEGER *m, LAPACK_DOUBLEREAL *s, LAPACK_DOUBLEREAL *sep, 
			     LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrsna_(char *job, char *howmny, LAPACK_LOGICAL *select, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *t, LAPACK_INTEGER *ldt, LAPACK_DOUBLECOMPLEX *vl, 
			     LAPACK_INTEGER *ldvl, LAPACK_DOUBLECOMPLEX *vr, LAPACK_INTEGER *ldvr, LAPACK_DOUBLEREAL *s, 
			     LAPACK_DOUBLEREAL *sep, LAPACK_INTEGER *mm, LAPACK_INTEGER *m, LAPACK_DOUBLECOMPLEX *work, 
			     LAPACK_INTEGER *ldwork, LAPACK_DOUBLEREAL *rwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrsyl_(char *trana, char *tranb, LAPACK_INTEGER *isgn, LAPACK_INTEGER 
			     *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLEREAL *scale, 
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrti2_(char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrtri_(char *uplo, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztrtrs_(char *uplo, char *trans, char *diag, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *nrhs, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *b, 
			     LAPACK_INTEGER *ldb, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztzrqf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_INTEGER *info);
 
/* Subroutine */ int ztzrzf_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zung2l_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zung2r_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungbr_(char *vect, LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunghr_(LAPACK_INTEGER *n, LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungl2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunglq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungql_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungqr_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungr2_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungrq_(LAPACK_INTEGER *m, LAPACK_INTEGER *n, LAPACK_INTEGER *k, 
			     LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zungtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, 
			     LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zunm2l_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunm2r_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmbr_(char *vect, char *side, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX 
			     *tau, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmhr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *ilo, LAPACK_INTEGER *ihi, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, 
			     LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunml2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmlq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmql_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmqr_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmr2_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmr3_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX 
			     *tau, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     info);
 
/* Subroutine */ int zunmrq_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmrz_(char *side, char *trans, LAPACK_INTEGER *m, LAPACK_INTEGER *n, 
			     LAPACK_INTEGER *k, LAPACK_INTEGER *l, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX 
			     *tau, LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *
			     lwork, LAPACK_INTEGER *info);
 
/* Subroutine */ int zunmtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *a, LAPACK_INTEGER *lda, LAPACK_DOUBLECOMPLEX *tau, 
			     LAPACK_DOUBLECOMPLEX *c__, LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *lwork,
			     LAPACK_INTEGER *info);
 
/* Subroutine */ int zupgtr_(char *uplo, LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, 
			     LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *q, LAPACK_INTEGER *ldq, LAPACK_DOUBLECOMPLEX *
			     work, LAPACK_INTEGER *info);
 
/* Subroutine */ int zupmtr_(char *side, char *uplo, char *trans, LAPACK_INTEGER *m, 
			     LAPACK_INTEGER *n, LAPACK_DOUBLECOMPLEX *ap, LAPACK_DOUBLECOMPLEX *tau, LAPACK_DOUBLECOMPLEX *c__,
			     LAPACK_INTEGER *ldc, LAPACK_DOUBLECOMPLEX *work, LAPACK_INTEGER *info);

#endif /* __CLAPACK_H */
