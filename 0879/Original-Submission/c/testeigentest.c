#include "eigentest.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/*
  This program tests the package Eigentest.  It runs 64 testcases
  probing various aspects of the package, as described below.  The
  numbers in the output should be within one or two orders of
  magnitude of the rounding unit.
*/

/*
  These definitions may be adjusted at compile time to change the
  tests, which involve computing (A - shift*I) op B,
  where op is *, '*, \, or '\ (' denoting transpose).
*/

# define N          10   /* Order of A.                        */
# define NCOLS       3   /* Number of columns in B.            */
# define SHIFT    -1.0   /* A shift.                           */
# define SPINRAND    0   /* Number of times to call the        */
                         /* random number generator initially. */

int main(int argc, char *argv[]){
  

  struct eigenmat A;
  int n=N, ncols = NCOLS, spinrand = SPINRAND;
  int i, j, i1, i2;
  int cs, ztype, ytype, atype;
  int zblk[4] = {1, 2, 3, 2};
 
  double B[N][N], C[N][N], D[N][N];
  double shift = SHIFT, norm1 = 0.0;
  double eigr[N], eigi[N], ev1[n];

  /* Spin the random number generator. */

  for (i=1; i<=spinrand; i++)
      rand();

  /* Loop on test case number. */

  for(cs=0; cs<64; cs++){
    printf("Test case %d\n", cs);

    /*
      Choose the type of the outer hsvd Y.  It can be an identity
      (ytype==0)  or a random hsvd (ytype==1).
    */

    ytype = cs%2;

    /*
      Choose the typep of the inner hsvd Z.    

         If ztype == 0, Z is an identity.
         If ztype == 1, Z has two blocks: (1:6)(7:10).
           The first block is identity.
         If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10),
           and the second block is an indentity.
         If ztype == 3, Z has two blocks: (1:4)(5::10),
           and the second block is identity.
    */

    ztype = (cs/2)%4;

    /*
      Chose the types of eigenvalues.
      
        If atype == 0, all eigenvalues are real.
        If atype == 1, all eigenvalues are complex.
        If atype == 2, eigenvalues 1,2,5,6,9,10 are complex.
        If atype == 3, eigenvalues 3,4,7,8 are complex.
    */

    atype = (cs/8)%8;

    /* Initialize A. */

    EigenmatAlloc(&A, n, zblk[ztype], ytype == 0, ztype == 0);

    /* Generate eigenvalues. */

    for(i = 0; i< n; i++){
      A.eig[i]  =  ((double) rand())/RAND_MAX  + 0.5;
    }

    switch(atype){

    case 0:
      for(i=0; i<n; i++) 
        A.type[i] = 1;
      break;

    case 1:
      for(i=0; i<n; i+=2){ 
        A.type[i] = 2;
        A.type[i+1] = 3;
      }
      break;

    case 2:
      A.type[0] = A.type[4] = A.type[8] = 2;
      A.type[1] = A.type[5] = A.type[9] = 3;
      A.type[2] = A.type[3] = A.type[6] = A.type[7] = 1;
      break;

    case 3:
      A.type[2] = A.type[6] = 2;
      A.type[3] = A.type[7] = 3;
      A.type[0] = A.type[1] = A.type[4] = A.type[5] = 1;
      A.type[8] = A.type[9] = 1;
      break;

    case 4:
      A.type[0] = -n;
      for(i=1; i<n; i++) 
        A.type[i] = -1;
      break;

    case 5:
      A.type[2] = -2;
      A.type[6] = -3;
      A.type[3] = A.type[8] = A.type[7] = -1;
      A.type[0] = A.type[1] = A.type[4] = A.type[5] = A.type[9] = 1;
      break;

    case 6:
      A.type[0] = A.type[7] = -3;
      A.type[1] = A.type[2] = -1;
      A.type[3] = A.type[4] = A.type[5] = A.type[6] = 1;
      A.type[8] = A.type[9] = -1;
      break;

    case 7:
      A.type[1] = -3;
      A.type[2] = A.type[3] = -1;
      A.type[0] = A.type[4] = A.type[9] = 1;
      A.type[5] = A.type[7] = 2;
      A.type[6] = A.type[8] = 3;
      break;

    }

    printf("A.type: ");
    for(j=0; j<n; j++){
      printf("%d ", A.type[j]);
    }
    printf("\n");

    
    /* Generate Y. */

    if (ytype == 1) {

      for(i=0; i<n; i++){
        A.Y.sig[i]  =  ((double) rand())/RAND_MAX  + 1.0;
        A.Y.u[i] = ((double) rand())/RAND_MAX  - 0.5;
        A.Y.v[i] = ((double) rand())/RAND_MAX  - 0.5;
      }
      hscal(n, &A.Y.u[0]);
      hscal(n, &A.Y.v[0]);
    }

    printf("Y.bs = %d %d\n", A.Y.bs[0], A.Y.bs[1]);

    /* generate Z. */

    if(ztype >0 ){

      for(i=0; i<n; i++){
        A.Z.sig[i]  =  ((double) rand())/RAND_MAX  + 1.0;
        A.Z.u[i] = ((double) rand())/RAND_MAX  - 0.5;
        A.Z.v[i] = ((double) rand())/RAND_MAX  - 0.5;
      }

      switch(ztype){

      case 1:
        /*
          If ztype == 1, Z has two blocks: (1:6)(7:10).
          The first block is identity.
        */

        A.Z.bs[1] = -6;
        break;

      case 2:
        /*
          If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10).
          The second block is indentity.
         */

        A.Z.bs[1] = 2;
        A.Z.bs[2] = -7;
        break;

      case 3:
        /*
          If ztype == 3, Z has two blocks: (1:4)(5::10).
          The second block is identity.
         */

        A.Z.bs[1] = 4;
        A.Z.bs[2] = -n;
        break;
      }
      
      for(j=0; j<zblk[ztype]; j++){
        i1 = abs(A.Z.bs[j]);
        i2 = A.Z.bs[j+1];
        if(i2>0){
	  hscal(i2-i1, &A.Z.u[i1]);
	  hscal(i2-i1, &A.Z.v[i1]);
        }
      }
    }


    printf("Z.bs = ");
    for(j=0; j<=zblk[ztype]; j++){
      printf("%d ", A.Z.bs[j]);
    }
    printf("\n");


    /*
      Generate B for A, Y, and Z to operate on.
    */

    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	B[i][j] = C[i][j] = ((double) rand())/RAND_MAX - 0.5;
      }
    }
    
    /* Test Z. */

    HsvdProd(&A.Z, ncols, (double *)C, n, "ab");
    HsvdProd(&A.Z, ncols, (double *)C, n, "aib");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Z\\Z*x - x|        = %e\n", norm1);

    HsvdProd(&A.Z, ncols, (double *)C, n, "aib");
    HsvdProd(&A.Z, ncols, (double *)C, n, "ab");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Z*Z\\x - x|        = %e\n", norm1);

    HsvdProd(&A.Z, ncols, (double *) C, n, "atb");
    HsvdProd(&A.Z, ncols, (double *) C, n, "aitb");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Z'\\Z'*x - x|      = %e\n", norm1);

    HsvdProd(&A.Z, ncols, (double *) C, n, "aitb");
    HsvdProd(&A.Z, ncols, (double *) C, n, "atb");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Z'*Z'\\x - x|      = %e\n", norm1);

    /* Test Y. */
    HsvdProd(&A.Y, ncols, (double *) C, n, "ab");
    HsvdProd(&A.Y, ncols, (double *) C, n, "aib");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Y\\Y*x - x|        = %e\n", norm1);

    HsvdProd(&A.Y, ncols, (double *) C, n, "aib");
    HsvdProd(&A.Y, ncols, (double *) C, n, "ab");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Y*Y\\x - x|        = %e\n", norm1);

    HsvdProd(&A.Y, ncols, (double *) C, n, "atb");
    HsvdProd(&A.Y, ncols, (double *) C, n, "aitb");

    norm1 = 0.0;
    for(j=0; j<ncols; j++){
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
	C[i][j] = B[i][j];
      }
    }
    printf("|Y'\\Y'*x - x|      = %e\n", norm1);

    HsvdProd(&A.Y, ncols, (double *) C, n, "aitb");
    HsvdProd(&A.Y, ncols, (double *) C, n, "atb");

    norm1 = 0.0;
    for(j=0; j<ncols; j++) {
      for(i=0; i<n; i++){
	norm1 += fabs(C[i][j]-B[i][j]);
      }
    }
    printf("|Y'*Y'\\x - x|      = %e\n", norm1);

    /* Test A. */
    EigenmatProd(&A, ncols, (double *) B, n, (double *) C, n, shift, "ab");
    EigenmatProd(&A, ncols, (double *) C, n, (double *) D, n, shift, "aib");
    
    norm1 = 0.0;
    for(j=0; j<ncols; j++) {
      for(i=0; i<n; i++){
	norm1 += fabs(D[i][j]-B[i][j]);
      }
    }
    printf("|A\\A*x - x|        = %e\n", norm1);
    
    EigenmatProd(&A, ncols, (double *) B, n, (double *) C, n, shift, "aib");
    EigenmatProd(&A, ncols, (double *) C, n, (double *) D, n, shift, "ab");
    
    norm1 = 0.0;
    for(j=0; j<ncols; j++) {
      for(i=0; i<n; i++){
	norm1 += fabs(D[i][j]-B[i][j]);
      }
    }
    printf("|A*A\\x - x|        = %e\n", norm1);

    EigenmatProd(&A, ncols, (double *) B, n, (double *) C, n, shift, "atb");
    EigenmatProd(&A, ncols, (double *) C, n, (double *) D, n, shift, "aitb");
    
    norm1 = 0.0;
    for(j=0; j<ncols; j++) {
      for(i=0; i<n; i++){
	norm1 += fabs(D[i][j]-B[i][j]);
      }
    }
    printf("|A'\\A'*x - x|      = %e\n", norm1);

    EigenmatProd(&A, ncols, (double *) B, n, (double *) C, n, shift, "aitb");
    EigenmatProd(&A, ncols, (double *) C, n, (double *) D, n, shift, "atb");
    
    norm1 = 0.0;
    for(j=0; j<ncols; j++) {
      for(i=0; i<n; i++){
	norm1 += fabs(D[i][j]-B[i][j]);
      }
    }    
    printf("|A'*A'\\x - x|      = %e\n", norm1);

    /*
       Test the eigenvector calculations by computing all
       left and right eigenvectors and their residuals.
    */

    /* Right eigensystem. */

    for (i=0; i<n;){
       EigenmatVecs(&A, i, &eigr[i], &eigi[i], &C[i][0], &D[i][0],
		    NULL, NULL, NULL, 'r');

       if(A.type[i] == 1 || A.type[i] < -1){
	  ev1[i] = 0;
          i++;
       } 
       else if (A.type[i] == -1){
	  ev1[i] = A.eig[i];
          i++;  
       }
       else  {

          for(j=0;j<n;j++) {
             C[i+1][j] =  C[i][j];
             D[i+1][j] = -D[i][j];
          }
	  eigr[i+1] =  eigr[i];
	  eigi[i+1] = -eigi[i];
	  ev1[i] = ev1[i+1] = 0;
          i+=2;
       }
    }

    norm1 = 0.0;

    /* Real part. */
    EigenmatProd(&A, 1, &C[0][0], 1, &B[0][0], 1, 0.0, "ab");
    /* Imaginary part. */
    EigenmatProd(&A, 1, &D[0][0], 1, &B[1][0], 1, 0.0, "ab");

    for (j=0;j<n;j++) {
      norm1 += fabs(eigr[0]*C[0][j] - eigi[0]*D[0][j] - B[0][j]);
      norm1 += fabs(eigr[0]*D[0][j] + eigi[0]*C[0][j] - B[1][j]);
    } 

    for (i=1; i<n; i++) {

      /* Real part. */
      EigenmatProd(&A, 1, &C[i][0], 1, &B[0][0], 1, 0.0, "ab");

       /* Imaginary part. */
      EigenmatProd(&A, 1, &D[i][0], 1, &B[1][0], 1, 0.0, "ab");

      for (j=0;j<n;j++) {
        norm1 += fabs(ev1[i]*C[i-1][j] + eigr[i]*C[i][j] - eigi[i]*D[i][j] - B[0][j]);
	norm1 += fabs(eigr[i]*D[i][j] + eigi[i]*C[i][j] - B[1][j]);
      } 
    }    

    printf("|A*REV - REV*EV|   = %e\n",norm1);

    /* Left eigensystem. */

    for (i=0; i<n;) {

      EigenmatVecs(&A, i, &eigr[i], &eigi[i],
		   NULL, NULL, &C[i][0], &D[i][0], NULL, 'l');

      if (A.type[i] == 1 || A.type[i] < -1) {
	ev1[i] = 0;
        i++;
      } 
      else if (A.type[i] == -1) {
	ev1[i] = A.eig[i];
	i++;
      } 
      else{

        for(j=0;j<n;j++) {
          C[i+1][j] =  C[i][j];
          D[i+1][j] = -D[i][j];
        }
	eigr[i+1]= eigr[i];
	eigi[i+1]= -eigi[i];
	ev1[i] = ev1[i+1] = 0;
        i+=2;
      } 
    }

    norm1 = 0.0;
    for (i=0; i<n-1; i++) {

      /* real part. */
      EigenmatProd(&A, 1, &C[i][0], 1, &B[0][0], 1, 0.0, "atb");

      /* complex part. */
      EigenmatProd(&A, 1, &D[i][0], 1, &B[1][0], 1, 0.0, "atb");

      for (j=0;j<n;j++) {
        norm1 += fabs(ev1[i+1]*C[i+1][j] + eigr[i]*C[i][j] + eigi[i]*D[i][j] - B[0][j]);
        norm1 += fabs(eigr[i]*D[i][j] - eigi[i]*C[i][j] - B[1][j]);
      } 
    }

    /* real part. */
    EigenmatProd(&A, 1, &C[n-1][0], 1, &B[0][0], 1, 0.0, "atb");
    /* complex part. */
    EigenmatProd(&A, 1, &D[n-1][0], 1, &B[1][0], 1, 0.0, "atb");
    for (j=1;j<n;j++) {
      norm1 += fabs(eigr[n-1]*C[n-1][0] + eigi[n-1]*D[n-1][0] - B[0][0]);
      norm1 += fabs(eigr[n-1]*D[n-1][0] - eigi[n-1]*C[n-1][0] - B[1][0]);
    } 

    
    printf("|A'*LEV - LEV*EV'| = %e\n\n", norm1);

    EigenmatFree(&A);
  }

  return 0;
}



