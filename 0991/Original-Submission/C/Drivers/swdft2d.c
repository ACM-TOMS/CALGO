/* 
  This file implements a 2D SWDFT using a 2D discrete Fourier transform (DFT) in each 
  window position. This is straight from the 2D SWDFT definition
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "tswdft2d.h"

double complex * swdft2d(double complex *x, int n0, int n1, int N0, int N1) {
	int P0 = N0 - n0 + 1;
	int P1 = N1 - n1 + 1;

	double complex *a;	
	a = malloc(sizeof(double complex) * P0 * P1 * n0 * n1);
	
	int phat0, phat1;
	double complex tot, xval, tf0, tf1;

	int p0, p1, k0, k1, j0, j1;

	for (p0 = (n0 - 1); p0 < N0; p0++) {
		phat0 = p0 - n0 + 1;

		for (p1 = (n1 - 1); p1 < N1; p1++) {
			phat1 = p1 - n1 + 1;

			for (k0 = 0; k0 < n0; k0++) {
				for (k1 = 0; k1 < n1; k1++) {

					tot = 0;

					for (j0 = 0; j0 < n0; j0++) {
						for (j1 = 0; j1 < n1; j1++) {
							
							xval = x[X_LOOKUP(phat0 + j0, phat1 + j1)];
							tf0 = cpow(cexp((2 * PI * I) / n0), -(j0 * k0));
							tf1 = cpow(cexp((2 * PI * I) / n1), -(j1 * k1));
							tot = tot + (xval * tf0 * tf1);
						}
					}

					a[COEF_LOOKUP(phat0, phat1) + WINDOW_LOOKUP(k0, k1)] = tot;
				}
			}
		}
	}

	return a;
}
