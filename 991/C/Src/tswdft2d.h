#include <math.h>
#include <complex.h>

#define PI 3.14159265358979323846
#define max(a, b) (((a) > (b)) ? (a) : (b)) 

/*  Gets the 1D node index based on the 2D binary tree position used by the tswdft2d function */ 
#define NODE(node0, node1) (((node0) * (n1)) + (node1))

/* Calculates the node index for a particular level of a tree */ 
#define IMOD(node, level) ((node) % (int) (pow((2), (level))))

/* Get index of tree inside the level_prev or level_cur arrays in the tswdft2d function */
#define TREE_2D(x, y) ((((x) * (N1)) + (y)) * ((n0) * (n1)))

/* Used to initialize level 0 of binary trees to the data  */
#define X_LOOKUP(x, y) (((x) * (N1)) + (y))

/* Used as a lookup for the output 4D array */
#define COEF_LOOKUP(x, y) ((((x) * (P1)) + (y)) * ((n0) * (n1)))

/* Lookups used in 2D SWFFT and 2D SWDFT sliding window functions  */
#define WINDOW_LOOKUP(x, y) (((x) * (n1)) + (y))

/* Declare all C programs in this package  */

/* The 2D Radix-2 Tree Sliding Window Discrete Fourier Transform */ 
void swap(double complex **a, double complex **b);
double complex * tswdft2d(double complex *x, int n0, int n1, int N0, int N1);

/* The 2D Sliding Window discrete Fourier transform */ 
double complex * swdft2d(double complex *x, int n0, int n1, int N0, int N1);

/* The 2D Sliding Window Fast Fourier Transform */
double complex * swfft2d(double complex *x, int n0, int n1, int N0, int N1);
