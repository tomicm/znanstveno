#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "deklaracije.h"
#include "korisno.h"
#include "cg.h"

int main(void) {
  int n, i, j;
  FILE *f;
  double *A;

  n = 100;

  A = (double*)malloc(n * n * sizeof(double));

  f = fopen("stieltjes_matr.txt", "r");
  for (i = 0; i < n*n; ++i)
      fscanf(f, "%lf", A+i);
  fclose(f);
  
  double *x, *x0, *b;
  double alpha, beta;
  double tol = 1e-8;
	
  x = (double*)malloc(n * sizeof(double));
  x0 = (double*)malloc(n * sizeof(double));
  b = (double*)malloc(n * sizeof(double));
  
  // egzaktno rjesenje x je [1 ... 1]^T
  fill(n, 1, x, 1.0);
  
  // gradimo b
  alpha = 1.0; beta = 0.0;
  dsymv_("U", &n, &alpha, A, &n, x, &INT1, &beta, b, &INT1);

  // pocetna iteracija x0 je [0 ... 0]^T
  fill(n, 1, x0, 0.0);

  cg(n, A, b, x0, tol);
  
  print_vector("rjesenje = ", n, x0);
  printf("relativna greska: %g\n", distance(n, x, x0) / dnrm2_(&n, x, &INT1));
  
  return 0;
}
