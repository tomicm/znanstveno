#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "deklaracije.h"
#include "korisno.h"

#include "cg.h"

int main(void) {
	// primjer s Wikipedije
	int n = 2;

	double A[] = {4, 1, 1, 3},
				 b[] = {1, 2},
				 x[] = {2, 1};

	printf("Rjesavamo Ax = b za: \n");
	print_matrix("A = ", n, n, A);
	print_vector("b = ", n, b);

	cg(n, A, b, x, 1e-5);

	printf("----------------------------\n");
	print_vector("rjesenje = ", n, x);

	return 0;

}
