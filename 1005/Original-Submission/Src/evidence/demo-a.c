#include <stdio.h>
void daxpy_(int*, double*, double*, int*, double*, int*);
void daxpy_rmd_(int*, double*, int*, int*, double*, double*);
int main() {
   double x[2] = {1., 2.};
   double y[2] = {3., 4.};
   double xa[4] = {1., 1.};
   double ya[4] = {2., 3.};
   double alpha = 2.5;
   int n = 2, one = 1;
   daxpy_(&n, &alpha, x, &one, y, &one);
   daxpy_rmd_(&n, &alpha, &one, &one, xa, ya);
   printf("xa = [%.2f, %.2f]\n", xa[0], xa[1]);
}
