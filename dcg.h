#include "deklaracije.h"
#include "korisno.h"

void dcg(int n, double *A, double *b, double *x, double tol) {
  int sz, k, i;
  double *dg, *d, *r, *p, *v;
  double alpha, beta;
  double norma_r, norma_b;
  double dot_rp, aa, bb;

  d = (double*)malloc(n * sizeof(double));
  r = (double*)malloc(n * sizeof(double));
  p = (double*)malloc(n * sizeof(double));
  v = (double*)malloc(n * sizeof(double));
  dg = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; ++i) dg[i] = 1 / sqrt(A[i*n+i]);
  
  dcopy_(&n, b, &INT1, r, &INT1);

  alpha = -1.0; beta = 1.0;
  dsymv_("U", &n, &alpha, A, &n, x, &INT1, &beta, r, &INT1);

  dcopy_(&n, r, &INT1, p, &INT1);
  for (i = 0; i < n; ++i) p[i] *= dg[i];

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
    for (i = 0; i < n; ++i) p[i] *= dg[i];
    
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
  free(dg);

}
