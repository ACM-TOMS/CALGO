#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "tswdft2d.h"

/*
  Swap pointers of two equal-sized complex arrays  

  INPUT
  a --- double complex ** --- Pointer to a pointer of a double complex 
  b --- double complex ** --- Pointer to a pointer of a double complex
*/
void swap(double complex **a, double complex **b) {
  double complex *temp = *a;
  *a = *b;
  *b = temp;
}
