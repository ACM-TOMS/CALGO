/*
________________________________________________________________

Name:      Gauss_Pivot.c

Purpose:   This program implements Gaussian Elimination
           with pivoting. The program also demonstrates
		   the use of offsets so that array indices start
		   at 1 not 0. The program will read multiple sets
		   of problems and print the answers.
			
Author(s): R. Sureshkumar (10 January 1997)
           Gregory J. McRae (22 October 1997)

Address:   Department of Chemical Engineering
           Room 66-372
           Massachusetts Institute of Technology
           Cambridge, MA 02139
		   mcrae@mit.edu

Usage:     The program expects a data file 'gauss.dat' that is 
           structured as follows:
		   
		   Line 1:   A brief description (< 80 characters)
		   Line 2:   The dimension of the problem (n >= 1)
		   Line 3+:  Row 1 to n of the matrix A
		   Line n+3: Row 1 to n of the matrix b

Example format of "gauss.dat"
________________________________________________________________		     
Strang (Applied Mathematics, p.10)  Solution = x = {1, 0, 0, 4}  
4
 2.  1.  0.  0
 1.  2.  1.  0
 0.  1.  2.  1
 0.  0.  1.  2
2.  
1.  
4.  
8.
______________________________________________________________

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ICHAR 80  /*  Length of array holding description of
                      the problem  */

#define IDEBUG 1  /*  Flag to enable printing of intermediate
                      results of decomposition =1 yes, =0 no. */


/*   Function Prototypes  */

int  matrix_print_off (int nr, int nc, double **A);
int  vector_print_off (int nr, double *x);
void gauss(double **a, double *b, double *x, int n);

int main (void)
{
double  **a,*b,*x;
float  aij,bi;
char   desc[ICHAR];  
int    i,j,n;  
FILE   *finput;    

/* Open the file containing a description, n,A and b  */

finput = fopen("gauss.dat","r");  
if (finput == NULL) { 
						printf("Data file gauss.dat not found\n");
						return(-1);
					}

/*  
    Get a one line description of the matrix problem and
    then the dimension, n, of the system A[n][n] and b[n]  
*/
fgets(desc, ICHAR , finput);      
fscanf(finput, "%d",&n);
printf("%s", desc);  
printf("\nDimension of matrix = %d\n\n",n);

/* 
   Dynamic allocation of the arrays and vectors, 
   notice the --a; to offset the pointer so that 
   the array index goes from 1 to n rather that 
   from 0 to n-1 
*/ 

a = calloc(n, sizeof(double *)); --a;

for(i = 1; i<=n; ++i) { 
	a[i] =  calloc(n, sizeof(double)); --a[i];
}

b = calloc(n, sizeof(double)); --b;
x = calloc(n, sizeof(double)); --x;

/*   Read the elements of A */

for (i=1;i<=n;i++){
	for (j=1;j<=n;j++) {
		fscanf(finput,"%f ",&aij);
		a[i][j] = (double) aij;
	}
}

/*  Read the elements of b */

for (i=1;i<=n;i++){
   fscanf(finput,"%f ",&bi);
   b[i] = (double) bi;
}

fclose(finput); /*  Close the input file  */

printf("\nMatrices read from input file\n");

printf("\nMatrix A\n\n");
matrix_print_off (n,n,a);

printf("\nVector b\n\n");
vector_print_off (n,b);

/* Call the Gaussian elimination function */

	gauss(a,b,x,n);

printf("\nSolution x\n\n");
vector_print_off (n,x);

return(0);
}

void gauss(double **a, double *b, double *x, int n)
/*
________________________________________________________________

Name:      gauss

Purpose:   This program uses Gaussian Elimination with
           with pivoting to solve the problem A x =b.
		   Use is made of array offsets so that the
		   indices go from 1 to n.

Author(s): R. Sureshkumar (10 January 1997)
           Modified by: Gregory J. McRae (22 October 1997)

Address:   Department of Chemical Engineering
           Room 66-372
           Massachusetts Institute of Technology
           Cambridge, MA 02139
		   mcrae@mit.edu

Usage:     gauss(a,b,x,n)

		   a   - Matrix a[n][n]
		   b   - Right hand side vector b[n]
		   x   - Desired solution vector
		   n   - Matrix dimensions

           If IDEBUG is set to 1 the program will print
           out the intermediate decompositions as well
           as the number of row interchanges.
*/

{
int   i,j,k,m,rowx;
double xfac,temp,temp1,amax;

/*
_______________________________________

  Do the forward reduction step.
_______________________________________

*/

rowx = 0;   /* Keep count of the row interchanges */
for (k=1; k<=n-1; ++k) {


     amax = (double) fabs(a[k][k]) ;
     m = k;
     for (i=k+1; i<=n; i++){   /* Find the row with largest pivot */
               xfac = (double) fabs(a[i][k]);
               if(xfac > amax) {amax = xfac; m=i;}
     }
     if(m != k) {  /* Row interchanges */
                 rowx = rowx+1;
                 temp1 = b[k];
                 b[k]  = b[m];
                 b[m]  = temp1;
                 for(j=k; j<=n; j++) {
                       temp = a[k][j];
                       a[k][j] = a[m][j];
                       a[m][j] = temp;
                 }
      }
       for (i=k+1; i<=n; ++i) {
          xfac = a[i][k]/a[k][k];

               for (j=k+1; j<=n; ++j) {
                   a[i][j] = a[i][j]-xfac*a[k][j];
               }
          b[i] = b[i]-xfac*b[k];
       }

if(IDEBUG == 1) {printf("\n A after decomposition step %d\n\n",k);
			     matrix_print_off (n, n, a);}

}
/*
_______________________________________

  Do the back substitution step
_______________________________________

*/

for (j=1; j<=n; ++j) {
  k=n-j+1;
  x[k] = b[k];
       for(i=k+1; i<=n; ++i) {
         x[k] = x[k]-a[k][i]*x[i];
       }
  x[k] = x[k]/a[k][k];
}

if(IDEBUG == 1) printf("\nNumber of row exchanges = %d\n",rowx);

}

int matrix_print_off (int nr, int nc, double **A)

/*
________________________________________________________________

Name:     matrix_print_off.c

Purpose:  This function will print out a general two-dimensional
          matrix in a user defined format using offsets that assume
		  the indices are from 1 to nr and 1 to nc.

Author:   Gregory J. McRae (22 October 1997)

Address:  Department of Chemical Engineering
          Room 66-372
          Massachusetts Institute of Technology
          Cambridge, MA 02139
		  mcrae@mit.edu

Usage:    mat_print_off(nr, nc, A);


Input:    nr     - number of rows (must be >= 1)
		  nc     - number of columns (must be >= 1)
		  A      - Matrix A[nr][nc] to be printed

Return:   =  0 Matrix printed correctly
          = -1 Number of rows not >= 0
		  = -2 Number of columns not >= 0
________________________________________________________________
*/
{
int i,j;

if ( nr <= 0 ) return (-1);
if ( nc <= 0 ) return (-2);

for (i = 1; i <= nr; i++) {

 	for (j = 1; j <= nc; j++) {
		printf ("%9.4f  ", A[i][j]);
	}

	printf("\n"); /* Insert a new line at end of each row */
}
return (0);
}


int vector_print_off (int nr, double *x)

/*   
________________________________________________________________

Name:     vector_print_off.c

Purpose:  This function will print out a one-dimensional 
          vector in a user defined format, using offsets
		  that assume the indices are from 1 to nr.
			
Author:   Gregory J. McRae (22 October 1997)

Address:  Department of Chemical Engineering
          Room 66-372
          Massachusetts Institute of Technology
          Cambridge, MA 02139
		  mcrae@mit.edu

Usage:    vector_print_off(nr, x);


Input:    nr     - number of rows (must be >= 1)
          x      - Vector x[nr] to be printed

Return:   =  0 Vector printed correctly
          = -1 Number of rows not >= 0
________________________________________________________________
*/
{
int i;
  
if ( nr <= 0 ) return (-1);

for (i = 1; i <= nr; i++) {
	printf ("%9.4f  \n", x[i]);
}
printf("\n");  /* Insert a new line at the end  */
return (0);
}
