#ifndef __DEKLARACIJE__
#define __DEKLARACIJE__
extern double daxpy_(
    const int *n,
    const double *alpha,
    double *x,
    const int *incx,
    double *y,
    const int *incy);

extern double dscal_(
		const int *n,
		const double *da,
		double *dx,
		const int *incx);

extern void dlacpy_(
    const char *uplo,
    const int *m,
    const int *n,
    double *A,
    const int *lda,
    double *B,
    const int *ldb);

extern void dtrsm_(
    const char *side,
    const char *uplo,
    const char *transa,
    const char *diag,
    const int *m,
    const int *n,
    const double *alpha,
    double *A,
    const int *lda,
    double *B,
    const int *ldb);

extern double dlange_(
    const char *norm,
    const int *m,
    const int *n,
    double *A,
    const int *lda,
    double *work);

extern void dcopy_(
    const int *n,
    double *x,
    const int *incx,
    double *y,
    const int *incy);

extern void dgemv_(
    const char *trans,
    const int *m,
    const int *n,
    const double *alpha,
    double *A,
    const int *lda,
    double *x,
    const int *incx,
    const double *beta,
    double *y,
    const int *incy);

extern void dsymv_(
		const char *uplo,
		const int *n,
		const double *alpha,
		double *A,
		const int *lda,
		double *x,
		const int *incx,
		const double *beta,
		double *y,
		const int *incy);

extern double dnrm2_(
		const int *n,
		double *x,
		const int *incx);

extern double ddot_(
		const int *n,
		double *dx,
		const int *incx,
		double *dy,
		const int *incy);

extern void dlarnv_(
		const int *dist,
		int *seed,
		const int *n,
		double *x);

extern void dgeqrf_(
		const int *m, 
		const int *n, 
		double *A,
		const int *lda,
		double *tau,
		double *work,
		int *lwork,
		int *info);

extern void dormqr_(
		const char *side,
		const char *trans,
		const int *m,
		const int *n,
		const int *k,
		double *A,
		const int *lda,
		double *tau,
		double *C,
		const int *ldc,
		double *work,
		int *lwork,
		int *info);
#endif
