#include <stdio.h>
double ddot_(int*, double*, int*, double*, int*);
void ddot_rmd_(int*, double*, int*, double*, int*, double*, double*, double*, char*);
int main() {
   double x[2] = {1., 2.};
   double y[2] = {3., 4.};
   double xa[2] = {1., 1.};
   double ya[2] = {2., 3.};
   double dot;
   int n = 2, one = 1;
   dot = ddot_(&n, x, &one, y, &one);
   ddot_rmd_(&n, x, &one, y, &one, &dot, xa, ya, "11");
   printf("xa = [%.2f, %.2f]\n", xa[0], xa[1]);
   printf("ya = [%.2f, %.2f]\n", ya[0], ya[1]);
}
