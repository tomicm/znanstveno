#ifndef __KORISNO__
#define __KORISNO__

#include <math.h>

const int INT1 = 1;
const int INT2 = 2;
const int INT3 = 3;

void fill(int m, int n, double *A, double v) {
    int i;
    for (i = 0; i < m*n; ++i)
	A[i] = v;
}

double sqr(double x) { return x*x; }

double distance(int n, double *x, double *y) {
	int i;
	double ans = 0.0;
	for (i = 0; i < n; ++i) ans += sqr(x[i]-y[i]);
	return sqrt(ans);
}

void print_matrix(char *msg, int m, int n, double *A) {
    int i, j;
    printf("%s\n", msg);
    for (i = 0; i < m; ++i) {
	for (j = 0; j < n; ++j)
	    printf("%10.6lf ", A[j*m+i]);
	printf("\n");
    }
}

void print_vector(char *msg, int n, double *x) {
    int i;
    printf("%s [ ", msg);
    for (i = 0; i < n; ++i)
	printf("%.6lf ", x[i]);
    printf("]^T\n");
}
#endif
