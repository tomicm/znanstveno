#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "deklaracije.h"
#include "korisno.h"

#include "cg.h"

void build_matrix(int n, double *A) {
	double *T, *tau, *work;
	int seed[] = { time(NULL) % 4096, 0, 0, 1 };
	int info, size, lwork;
	int N, i, j;

	N = n; n = n * n;

	T = (double*)malloc(n * n * sizeof(double));
	work = (double*)malloc(n * n * sizeof(double));
	tau = (double*)malloc(n * sizeof(double));

	size = n * n;
	dlarnv_(&INT3, seed, &size, T);

	lwork = n; // ??
	dgeqrf_(&n, &n, T, &n, tau, work, &lwork, &info);

	int pos = 0;
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j, ++pos)
			A[pos*n+pos] = (double)(i+1);

	dormqr_("L", "N", &n, &n, &n, T, &n, tau, A, &n, work, &lwork, &info);
	dormqr_("R", "T", &n, &n, &n, T, &n, tau, A, &n, work, &lwork, &info);

	free(T);
	free(work);
	free(tau);
}

int main(void) {

	int N = 10;
	int n = N*N;
	double *A, *x, *x0, *b;
	double alpha, beta;
	double tol = 1e-8;
	
	A = (double*)malloc(n * n * sizeof(double));
	x = (double*)malloc(n * sizeof(double));
	x0 = (double*)malloc(n * sizeof(double));
	b = (double*)malloc(n * sizeof(double));

	build_matrix(N, A);
	// egzaktno rjesenje x je [1 ... 1]^T
	fill(n, 1, x, 1.0);

	// gradimo b
	alpha = 1.0; beta = 0.0;
	dsymv_("U", &n, &alpha, A, &n, x, &INT1, &beta, b, &INT1);

	// pocetna iteracija x0 je [0 ... 0]^T
	fill(n, 1, x0, 0.0);

	cg(n, A, b, x0, tol);

	printf("relativna greska: %g\n", distance(n, x, x0) / dnrm2_(&n, x, &INT1));

	return 0;

}
