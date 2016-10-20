const int INT1 = 1;

void fill(int m, int n, double *A, double v) {
    int i;
    for (i = 0; i < m*n; ++i)
	A[i] = v;
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
