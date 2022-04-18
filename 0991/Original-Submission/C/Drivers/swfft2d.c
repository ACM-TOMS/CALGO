/* 
  This file implements the 2D SWDFT using a 2D FFT algorithm in each window. The 2D FFT 
  works by taking a 1D FFT of each row, followed by a 1D FFT of the resulting columns. Therefore
  we use two extra functions for the 2D FFT and 1D FFTs
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "tswdft2d.h"

/* This function implements a 1D FFT */
void fft(double complex *x, int n, double complex *twf, double complex *a) {
	int shift, nodes;

	for (int level = 1; level <= log2(n); level++) {
		nodes = (int) pow(2, level);
		shift = floor((n / nodes)); 

		for (int tree = 0; tree < shift; tree++) {
			for (int node = 0; node < nodes; node++) {

				a[(tree * nodes) + node] = x[((tree * nodes) / 2) + IMOD(node, level - 1)] 
										 + (twf[node * shift] 
										 * x[(((tree + shift) * nodes) / 2) + IMOD(node, level - 1)]); 
			}
		}

		swap(&x, &a);
	}

	memmove(a, x, sizeof(double complex) * n);
}

/* This function implements a 2D FFT */
void fft2d(double complex *x, int n0, int n1) {
	// Generate the two twiddle factor vectors!
	double complex twiddle0[n0], twiddle1[n1];

	for (int init_tf0 = 0; init_tf0 < n0; init_tf0++) {
		twiddle0[init_tf0] = cpow(cexp((2 * PI * I) / n0), -init_tf0);
	}

	for (int init_tf1 = 0; init_tf1 < n1; init_tf1++) {
		twiddle1[init_tf1] = cpow(cexp((2 * PI * I) / n1), -init_tf1);
	}

	int size = max(n0, n1);
	double complex *fft_in, *fft_out;
	fft_in = malloc(sizeof(double complex) * size);
	fft_out = malloc(sizeof(double complex) * size);

	for (int row = 0; row < n0; row++) {
		for (int col = 0; col < n1; col++) {			
			fft_in[col] = x[WINDOW_LOOKUP(row, col)];
		}

		fft(fft_in, n1, twiddle1, fft_out);

		for (int col = 0; col < n1; col++) {
			x[WINDOW_LOOKUP(row, col)] = fft_out[col];
		}
	}

	for (int col = 0; col < n1; col++) {
		for (int row = 0; row < n0; row++) {
			fft_in[row] = x[WINDOW_LOOKUP(row, col)];
		}

		fft(fft_in, n0, twiddle0, fft_out);

		for (int row = 0; row < n0; row++) {
			x[WINDOW_LOOKUP(row, col)] = fft_out[row];
		}
	}

	free(fft_in);
	free(fft_out);
}

/* This function implements a 2D SWFT by calling a 2D FFT in each window */
double complex * swfft2d(double complex *x, int n0, int n1, int N0, int N1) {
	int P0 = N0 - n0 + 1;
	int P1 = N1 - n1 + 1;
	int phat0, phat1;

	double complex *a, *a_window;
	a = malloc(sizeof(double complex) * P0 * P1 * n0 * n1);
	a_window = malloc(sizeof(double complex) * n0 * n1);

	for (int p0 = (n0 - 1); p0 < N0; p0++) {
		phat0 = p0 - n0 + 1;

		for (int p1 = (n1 - 1); p1 < N1; p1++) {
			phat1 = p1 - n1 + 1;

			// Move the data into the temporary 2D array 
			for (int j0 = 0; j0 < n0; j0++) {
				for (int j1 = 0; j1 < n1; j1++) {
					a_window[WINDOW_LOOKUP(j0, j1)] = x[X_LOOKUP(phat0 + j0, phat1 + j1)];
				} 	
			}

			// Take the 2D FFT of the 2D temporary array 
			fft2d(a_window, n0, n1);

			// Store the correct portion of the 
			memcpy(&a[COEF_LOOKUP(phat0, phat1)], a_window, sizeof(double complex) * n0 * n1);
		}
	}

	free(a_window);

	return a;
}
