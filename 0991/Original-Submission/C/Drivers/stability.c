/* 
  This programs verifies that the 2D Tree SWDFT algorithm gives the same coefficients as the 2D SWFFT. 
*/ 

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include "tswdft2d.h"

int main() {
	/* Generate a random complex-valued 2D array */  
	int N0 = 9;
	int N1 = 11;
	double complex x[N0][N1];
	srand(100);

	for (int i = 0; i < N0; i++) {
		for (int j = 0; j < N1; j++) {
			x[i][j] = (rand() % 50) + ((rand() % 50) * _Complex_I);
		}
	}	

	/* Set the window sizes	 */
	int n0 = (int) pow(2, 2);
	int n1 = (int) pow(2, 3);
	
	/* Call the three different algorithms on the same input data */
	double complex *a_t, *a_fft;
	a_t  = tswdft2d(*x, n0, n1, N0, N1);
	a_fft = swfft2d(*x, n0, n1, N0, N1);

	for (int a_print = 0; a_print < ((N0 - n0 + 1) * (N1 - n1 + 1) * n0 * n1); a_print++) {
		/* Compare the 2D Tree SWDFT with the 2D SWFFT algorithm */ 
		printf("a_t[%d]: %18.13f + i %18.13f ..... a_fft[%d]: %18.13f + i %18.13f\n", 
			a_print, creal(a_t[a_print]), cimag(a_t[a_print]), 
			a_print, creal(a_fft[a_print]), cimag(a_fft[a_print]));

		assert(creal(a_t[a_print]) - creal(a_fft[a_print]) == 0);
		assert(cimag(a_t[a_print]) - cimag(a_fft[a_print]) == 0);
	}

	printf("Results are Identical! \n");

	free(a_fft);
	free(a_t);
}
