#include "deklaracije.h"
#include "korisno.h"

void ic(int n, double *A) {
  int i, j, k;
  for (i = 0; i < n; ++i) {
    for (k = 0; k < i; ++k)
      A[i*n+i] -= sqr(A[i*n+k]);
    A[i*n+i] = sqrt(A[i*n+i]);
    for (j = i+1; j < n; ++j) {
      if (A[j*n+i] != 0.0) {
	for (k = 0; k < i; ++k)
	  A[j*n+i] -= A[i*n+k] * A[j*n+k];
	A[j*n+i] /= A[i*n+i];
      }
    }
  }
}

void pcg(int n, double *A, double *b, double *x, double tol) {
  int sz, k;
  double *d, *r, *p, *v;
  double *R;
  double alpha, beta;
  double norma_r, norma_b;
  double dot_rp, aa, bb;

  R = (double*)malloc(n * n * sizeof(double));

  sz = n * n;
  dcopy_(&sz, A, &INT1, R, &INT1);
  ic(n, R);
  
  d = (double*)malloc(n * sizeof(double));
  r = (double*)malloc(n * sizeof(double));
  p = (double*)malloc(n * sizeof(double));
  v = (double*)malloc(n * sizeof(double));

  dcopy_(&n, b, &INT1, r, &INT1);

  alpha = -1.0; beta = 1.0;
  dsymv_("U", &n, &alpha, A, &n, x, &INT1, &beta, r, &INT1);

  dcopy_(&n, r, &INT1, p, &INT1);
  dtrsv_("U", "T", "N", &n, R, &n, p, &INT1);
  dtrsv_("U", "N", "N", &n, R, &n, p, &INT1);
  
  dcopy_(&n, p, &INT1, d, &INT1);

  norma_r = dnrm2_(&n, r, &INT1);
  norma_b = dnrm2_(&n, b, &INT1);

  dot_rp = ddot_(&n, r, &INT1, p, &INT1);

  k = 0;
  
  while (norma_r / norma_b >= tol) {
    // v = A * d
    alpha = 1.0; beta = 0.0;
    dsymv_("U", &n, &alpha, A, &n, d, &INT1, &beta, v, &INT1);

    aa = dot_rp / ddot_(&n, d, &INT1, v, &INT1);

    alpha = aa;
    daxpy_(&n, &alpha, d, &INT1, x, &INT1);

    alpha = -aa;
    daxpy_(&n, &alpha, v, &INT1, r, &INT1);

    dcopy_(&n, r, &INT1, p, &INT1);
    dtrsv_("U", "T", "N", &n, R, &n, p, &INT1);
    dtrsv_("U", "N", "N", &n, R, &n, p, &INT1);

    bb = 1.0 / dot_rp;
    dot_rp = ddot_(&n, r, &INT1, p, &INT1);
    bb *= dot_rp;

    alpha = bb;
    dscal_(&n, &alpha, d, &INT1);

    alpha = 1.0;
    daxpy_(&n, &alpha, p, &INT1, d, &INT1);

    norma_r = dnrm2_(&n, r, &INT1);

    ++k;
  }

  printf("broj iteracija: %d\n", k);

  free(d);
  free(r);
  free(p);
  free(v);
  free(R);

}
