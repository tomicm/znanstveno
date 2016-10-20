/* 
   Ulazni podaci:
   Program ucitava 100x100 matricu A iz datoteke stieltjes_matr.txt.
   Namjesta desnu stranu sustava (b) tako da je rjesenje [ 1 ... 1 ]^T
   i nakon toga rjesava sustav SOR metodom.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "deklaracije.h"
#include "korisno.h"

void sor_matrica(
    int n,
    double omega,
    double *A,
    double *T)
{
    int i, j;
    double *M;
    double alpha;

    M = (double*)malloc(n * n * sizeof(double));

    dlacpy_("L", &n, &n, A, &n, M, &n);
    for (i = 0; i < n; ++i)
	M[i*n+i] *= 1 / omega;

    fill(n, n, T, 0.0);
    dlacpy_("U", &n, &n, A, &n, T, &n);
    for (i = 0; i < n; ++i) {
	T[i*n+i] *= (1-omega)/omega;
	for (j = i+1; j < n; ++j)
	    T[j*n+i] *= -1;
    }

    alpha = 1.0;
    dtrsm_("L", "L", "N", "N", &n, &n, &alpha, M, &n, T, &n);

    free(M);
}

double sor_norma(
    int n,
    double omega,
    double *A,
    char norm)
{
    double result;
    double *T, *work;
    T = (double*)malloc(n * n * sizeof(double));
    work = (double*)malloc(n * sizeof(double));
    
    sor_matrica(n, omega, A, T);
    result = dlange_(&norm, &n, &n, T, &n, work);

    free(T);
    free(work);

    return result;
}

int sor_konvergencija(
    int n,
    double omega,
    double *A,
    char norm)
{
    return sor_norma(n, omega, A, norm) < 1;
}

void sor_iteracija(
    int n,
    double omega,
    double *A,
    double *b,
    double *x)
{
    int i, j;
    for (i = 0; i < n; ++i) {
	double pom;
	x[i] *= (1-omega);
	pom = b[i];
	for (j = 0; j < i; ++j)
	    pom -= A[j*n+i] * x[j];
	for (j = i+1; j < n; ++j)
	    pom -= A[j*n+i] * x[j];
	x[i] += pom * omega / A[i*n+i];
    }
}

void sor_rjesavac(
    int n,
    double omega,
    double *A,
    double *b,
    double *x,
    double epsilon,
    char norm)
{
    int i;
    double alpha, beta;
    double norm_b;
    double *res, *work;

    res = (double*)malloc(n * sizeof(double));
    work = (double*)malloc(n * sizeof(double));

    norm_b = dlange_(&norm, &n, &INT1, b, &n, work);
    
    while (1) {
	dcopy_(&n, b, &INT1, res, &INT1);
	alpha = 1.0; beta = -1.0;
	dgemv_("N", &n, &n, &alpha, A, &n, x, &INT1, &beta, res, &INT1);
	if (dlange_(&norm, &n, &INT1, res, &n, work) / norm_b < epsilon) break;
	sor_iteracija(n, omega, A, b, x);
    };

    free(res);
    free(work);
}

int main(void)
{
    int n, i, j;
    char norm;
    double omega, epsilon;
    double alpha, beta;
    double *A, *b, *x0, *x, *r, *work;
    FILE *f;
    
    n = 100;

    A = (double*)malloc(n * n * sizeof(double));
    b = (double*)malloc(n * sizeof(double));
    x0 = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    work = (double*)malloc(n * sizeof(double));

    f = fopen("stieltjes_matr.txt", "r");
    for (i = 0; i < n; ++i)
	for (j = 0; j < n; ++j)
	    fscanf(f, "%lf", &A[j*n+i]);
    fclose(f);

    // egzaktno rjesenje x = [ 1 .. 1 ]^T
    fill(n, 1, x, 1.0);
    // namjestamo b
    alpha = 1.0; beta = 0.0;
    dgemv_("N", &n, &n, &alpha, A, &n, x, &INT1, &beta, b, &INT1);
    
    // pocetna iteracija x0 = [ 0 .. 0 ]^T
    fill(n, 1, x0, 0.0);
    
    omega = 1.3;
    epsilon = 1e-5;
    norm = '1';

    printf("Rjesavamo sustav Ax = b SOR metodom uz omega = %g s pocetnom iteracijom x0.\n",
	   omega);
    printf("Trazimo rjesenje koje relativno odstupa od egzaktnog za %g u %c-normi.\n",
	   epsilon, norm);

    printf("\n|T|_%c = %.6lf, ", norm, sor_norma(n, omega, A, norm));
    printf("uvjet za konvergenciju u %c-normi %s ispunjen\n",
	   norm, sor_konvergencija(n, omega, A, norm) ? "je" : "nije");

    sor_rjesavac(n, omega, A, b, x0, epsilon, norm);

    r = (double*)malloc(n * sizeof(double));
    dcopy_(&n, x0, &INT1, r, &INT1);
    alpha = -1.0;
    daxpy_(&n, &alpha, x, &INT1, r, &INT1);

    printf("relativna pogreska = %g ...\n",
	   dlange_(&norm, &n, &INT1, r, &n, work) / dlange_(&norm, &n, &INT1, x, &n, work));

    free(A);
    free(b);
    free(x);
    free(x0);
    free(r);
    free(work);
    
    return 0;
}
