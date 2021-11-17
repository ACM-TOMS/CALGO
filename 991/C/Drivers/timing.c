#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include "tswdft2d.h"

int main() {
	// Timings for varying window sizes ---
	int N0 = 100;
	int N1 = 100;

	// Generate the random array 	
	double complex x[N0][N1];
	srand(100);

	for (int i = 0; i < N0; i++) {
		for (int j = 0; j < N1; j++) {
			x[i][j] = (rand() % 50) + ((rand() % 50) * _Complex_I);
		}
	}

	int window_range[5] = {1, 2, 3, 4, 5};
	for (int window_size = 0; window_size < 5; window_size++) {
		int m0 = window_range[window_size];
		int m1 = window_range[window_size];
		int n0 = (int) pow(2, m0);
		int n1 = (int) pow(2, m1);
		printf("N0: %d, N1: %d, n0: %d, n1: %d \n", N0, N1, n0, n1);

		/* 2D Sliding Window Discrete Fourier Transform */
		double complex *a_dft;
		clock_t start_swdft = clock(), diff_swdft;
		a_dft = swdft2d(*x, n0, n1, N0, N1);
		diff_swdft = clock() - start_swdft;
		int msec_swdft = (diff_swdft * 1000) / CLOCKS_PER_SEC;
		printf("2D SWDFT: %d seconds %d Milliseconds \n", msec_swdft / 1000, msec_swdft % 1000);
		free(a_dft);

		/* 2D Sliding Window Fast Fourier Transform */
		double complex *a_fft;
		clock_t start_swfft = clock(), diff_swfft;
		a_fft = swfft2d(*x, n0, n1, N0, N1);
		diff_swfft = clock() - start_swfft;
		int msec_swfft = (diff_swfft * 1000) / CLOCKS_PER_SEC;
		printf("2D SWFFT: %d seconds %d Milliseconds \n", msec_swfft / 1000, msec_swfft % 1000);
		free(a_fft);

		/* 2D Tree Sliding Window Discrete Fourier Transform */
		double complex *a_t;
		clock_t start_tswfft = clock(), diff_tswfft;
		a_t = tswdft2d(*x, n0, n1, N0, N1);
		diff_tswfft = clock() - start_tswfft;
		int msec_tswfft = (diff_tswfft * 1000) / CLOCKS_PER_SEC;
		printf("2D Tree SWDFT: %d seconds %d Milliseconds \n", msec_tswfft / 1000, msec_tswfft % 1000);
		free(a_t);
	}
}
