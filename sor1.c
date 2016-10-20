/* 
   Ulazni podaci:
   n       - dimenzija matrice
   A       - n x n matrica sustava
   b       - vektor iz R^n, desna strana sustava
   x0      - pocetna iteracija
   omega   - omega za SOR metodu
   epsilon - zeljena apsolutna preciznost
   norma   - '1' ili 'I' (norma koju cemo koristiti, 1-norma ili oo-norma)
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
    double *x0,
    double epsilon,
    char norm)
{
    int i;
    double alpha, norma_d, norma_T;
    double *x1, *work;

    work = (double*)malloc(n * sizeof(double));
    
    x1 = (double*)malloc(n * sizeof(double));
    dcopy_(&n, x0, &INT1, x1, &INT1);
    sor_iteracija(n, omega, A, b, x1);

    alpha = 1.0;
    daxpy_(&n, &alpha, x0, &INT1, x1, &INT1);
    norma_d = dlange_(&norm, &n, &INT1, x1, &n, work);
    norma_T = sor_norma(n, omega, A, norm);

    double k = ceil(log(epsilon*(1-norma_T) / norma_d) / log(norma_T));
    printf("dovoljno je %.0lf iteracija.\n", k);

    for (i = 0; i < k; ++i) sor_iteracija(n, omega, A, b, x0);
    
    free(x1);
    free(work);
}

int main(void)
{
    int n, i, j;
    char norm;
    double omega, epsilon;
    double *A, *T, *b, *x;
    
    scanf("%d", &n);

    A = (double*)malloc(n * n * sizeof(double));
    b = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    
    for (i = 0; i < n; ++i)
	for (j = 0; j < n; ++j)
	    scanf("%lf", &A[j*n+i]);

    for (i = 0; i < n; ++i)
	scanf("%lf", &b[i]);

    for (i = 0; i < n; ++i)
	scanf("%lf", &x[i]);

    scanf("%lf", &omega);
    scanf("%lf", &epsilon);
    scanf(" %c", &norm);

    printf("Rjesavamo sustav Ax = b SOR metodom uz omega = %g s pocetnom iteracijom x0.\n",
	   omega);
    print_matrix("A  = ", n, n, A);
    print_vector("b  = ", n, b);
    print_vector("x0 = ", n, x);
    printf("omega = %g\n", omega);
    printf("--------------------------------------\n");
    printf("Trazimo rjesenje koje odstupa od egzaktnog za %g u %c-normi.\n",
	   epsilon, norm);

    T = (double*)malloc(n * n * sizeof(double));
    sor_matrica(n, omega, A, T);

    print_matrix("matrica iteracija T = ", n, n, T);
    printf("\n|T|_%c = %.6lf, ", norm, sor_norma(n, omega, A, norm));
    printf("uvjet za konvergenciju u %c-normi %s ispunjen\n",
	   norm, sor_konvergencija ? "je" : "nije");

    sor_rjesavac(n, omega, A, b, x, epsilon, norm);

    print_vector("izracunato rjesenje = ", n, x);
    
    return 0;
}
