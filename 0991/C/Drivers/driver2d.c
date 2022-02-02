/*
  This driver demonstrates the tswdft2d function. First, we randomly generate a 2D array of 
  double complex numbers. Then we call the tswdft2d function on the 2D array. Finally, we print
  the output 2D SWDFT coefficients 
*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "tswdft2d.h"

int main() {
	/*
      Generate a N0 x N1 2D array called x. srand(100) seeds the random number 
	  generator. rand() generates a random number between 0 and RAND_MAX (defined in stdlib.h), 
	  so (rand % 50) randomly generates a number between 0 and 49, Since we are generating 
	  complex numbers, we assign a random number to both the real and imaginary parts. 
	*/
	int N0 = 11;
	int N1 = 9;
	double complex x[N0][N1];
	srand(100);
	for (int i = 0; i < N0; i++) {
		for (int j = 0; j < N1; j++) {
			x[i][j] = (rand() % 50) + ((rand() % 50) * _Complex_I);
			printf("x[%d, %d] = %f + i %f \n", i, j, creal(x[i][j]), cimag(x[i][j]));
		}
	}		

	/* Call the 2D Tree SWDFT algorithm on the randomly generated 2D array for window size n0 x n1 */
	int n0 = (int) pow(2, 3);
	int n1 = (int) pow(2, 2);

	double complex *a;
	a = tswdft2d(*x, n0, n1, N0, N1);

	/* Print the resulting 4D array using macros NODE and COEF_LOOKUP defined in the header fswft.h */		
	int P0 = N0 - n0 + 1;
	int P1 = N1 - n1 + 1;
	for (int p0 = 0; p0 < P0; p0++) {
		for (int p1 = 0; p1 < P1; p1++) {
			for (int node0 = 0; node0 < n0; node0++) {
				for (int node1 = 0; node1 < n1; node1++) {
					printf("a[%d, %d, %d, %d]: %f + i %f \n", p0, p1, node0, node1, 
														      creal(a[COEF_LOOKUP(p0, p1) + NODE(node0, node1)]), 
														      cimag(a[COEF_LOOKUP(p0, p1) + NODE(node0, node1)]));
				}
			}
		}
	}

	free(a);
}
