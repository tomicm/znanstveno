extern double daxpy_(
    const int *n,
    const double *alpha,
    double *x,
    const int *incx,
    double *y,
    const int *incy);

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
