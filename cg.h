#include "deklaracije.h"
#include "korisno.h"

void cg(int n, double *A, double *b, double *x, double tol) {
	double *d, *r, *v;
	double alpha, beta;
	double norma_r, norma_b;
	int k;
	
	d = (double*)malloc(n * sizeof(double));
	r = (double*)malloc(n * sizeof(double));
	v = (double*)malloc(n * sizeof(double));

	dcopy_(&n, b, &INT1, d, &INT1);
	
	alpha = -1.0; beta = 1.0;
	dsymv_("U", &n, &alpha, A, &n, x, &INT1, &beta, d, &INT1);
	dcopy_(&n, d, &INT1, r, &INT1);

	norma_r = dnrm2_(&n, r, &INT1);
	norma_b = dnrm2_(&n, b, &INT1);

	k = 0;

	while (norma_r / norma_b >= tol) {
		double r_sqsum, r2_sqsum, aa, bb;

		++k;

		r_sqsum = ddot_(&n, r, &INT1, r, &INT1);

		// v = A * d
		alpha = 1.0; beta = 0.0;
		dsymv_("U", &n, &alpha, A, &n, d, &INT1, &beta, v, &INT1);

		aa = r_sqsum / ddot_(&n, d, &INT1, v, &INT1);

		alpha = aa;
		daxpy_(&n, &alpha, d, &INT1, x, &INT1);

		alpha = -aa;
		daxpy_(&n, &alpha, v, &INT1, r, &INT1);

		r2_sqsum = ddot_(&n, r, &INT1, r, &INT1);
		bb = r2_sqsum / r_sqsum;

		alpha = bb;
		dscal_(&n, &alpha, d, &INT1);

		alpha = 1.0;
		daxpy_(&n, &alpha, r, &INT1, d, &INT1);

		norma_r = sqrt(r2_sqsum);
	}

	free(d);
	free(r);
	free(v);

	printf("broj iteracija: %d\n", k);
}
