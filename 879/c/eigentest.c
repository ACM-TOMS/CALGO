/*
  Eigentest is a package that produces real test matrices with known
  eigensystems.  A test matrix, called an eigenmat, is generated in a
  factored form, in which the user can specify the eigenvalues
  (including complex conjugate eigenvalues and Jordan blocks) and has
  some control over the condition of the eigenvalues and eigenvectors.
  An eigenmat $A$ of order $n$ requires only $O(n)$ storage for its
  representation.  Auxiliary programs permit the computation of
 
   (A - s*I)*b, (A - s*I)'*b, inv(A - s*I)*b, and inv(A -s*I)'*b
 
  in $O(n)$ operations.  Thus eigenmats are suitable for testing
  algorithms based on Krylov sequences, as well as others based on
  matrix vector products.  A special routine computes specified
  eigenvectors of an eigenmat and the condition of its eigenvalue.
 
  For more details see

     1. Eigentest: A Test Matrix Generator for Large Scale
        Eigenproblem

     2. Eigentest: The C User's Guide

  Pdf files are included with this distribution.

  Functions;
 
     EigenmatAlloc
        Allocates storage for an eigenmat.
 
     EigenmatFree
        Frees the storage for an eigenmat.
 
     EigenmatProd
        Computes the above products.
 
     HsvdProd
        Computes products of special matrices in the factorization.
 
     hhp
        Computes products with a Householder transformation.
 
     EigenmatVecs
        Computes eigenvectors of an eigenmat.
 
     hscal
        Scales a nonzero vector so that its norms is sqrt(2).

  Coded by Roger Lee and Pete Stewart
*/ 

#include "eigentest.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

# define M_SQRT2        1.41421356237309504880  /* sqrt(2) */
# define M_SQRT1_2      0.70710678118654752440  /* 1/sqrt(2) */

void EigenmatAlloc(struct eigenmat *A, int n, int nblocks,
                   int yident, int zident){
/*
   EigenmatAlloc allocates zeroed memory for an eigenmat A.
   In addition it initializes A.n, A.Z.n, A.Z.nblocks, A.Z.bs[0]
   A.Z.bs[n], A.y.n, A.Y.nblocks, and A.Y.bs.

   A        Pointer to the eigenmat to be initialized.
   n        The order of A.
   nblocks  The number of blocks in the hsvdmat A.Z.
   yident   If yident!=0, A.Y is initialized as an identity matrix.
   zident   If zident!=0, A.Z is initialized as an identity matrix.
*/

   A->n = n;
   A->eig = (double *) calloc(n, sizeof(double));
   A->type = (int *) calloc(n, sizeof(int));

   A->Y.n = n;
   A->Y.nblocks = 1;
   A->Y.bs = (int *) calloc(2, sizeof(int));
   if (yident){
     A->Y.bs[0] = 0;
     A->Y.bs[1] = -n;
     A->Y.u = (double *) NULL;
     A->Y.v = (double *) NULL;
     A->Y.sig = (double *) NULL;
   } 
   else {
     A->Y.bs[0] = 0;
     A->Y.bs[1] = n;
     A->Y.u = (double *) calloc(n, sizeof(double));
     A->Y.v = (double *) calloc(n, sizeof(double));
     A->Y.sig = (double *) calloc(n, sizeof(double));
   }

   A->Z.n = n;
   if (zident){  /* Z is an identity. */
     A->Z.nblocks = 1;
     A->Z.u = (double *) NULL;
     A->Z.v = (double *) NULL;
     A->Z.sig = (double *) NULL;
     A->Z.bs = (int *) calloc(2, sizeof(int));
     A->Z.bs[0] = 0;
     A->Z.bs[1] = -n;
   }
   else{ /* Z is not an identity. */
     A->Z.nblocks = nblocks;
     A->Z.u = (double *) calloc(n, sizeof(double));
     A->Z.v = (double *) calloc(n, sizeof(double));
     A->Z.sig = (double *) calloc(n, sizeof(double));
     A->Z.bs = (int *) calloc(nblocks+1, sizeof(int));
     A->Z.bs[0] = 0;
     A->Z.bs[nblocks] = n;
   }
}


void EigenmatFree(struct eigenmat *A){

   /* EigenmatFree deallocates the storage of the  eigenmat A. */

   if (A->eig != (double*) NULL) 
     free(A->eig);

   if (A->type != (int*) NULL) 
     free(A->type);
   
   if (A->Y.u != (double*) NULL) {
     free(A->Y.u);
     free(A->Y.v);
     free(A->Y.sig);
   }

   if (A->Y.bs != (int*) NULL) 
     free(A->Y.bs);

   if (A->Z.u != (double*) NULL) {
     free(A->Z.u);
     free(A->Z.v);
     free(A->Z.sig);
   }

   if (A->Z.bs != (int*) NULL) 
     free(A->Z.bs);
}


void HsvdProd(struct hsvd *X, int ncols, double *B, int tdb, char *job){

/*
   HsvdProd computes the product of computes the product of an hsvdmat
   matrix X and a matrix B, overwriting B with the product.

   X         Pointer to the hsvdmat.
   ncols     The number of columns in B.
   B         The array B.
   tdb       The trailing dimension of B.
   job       A string specifying the operation to be performed.

             "ab"   B <- X*B
             "atb"  B <- X^t*B
             "aib"  B <- X^-1*B
             "aitb" B <- X^-T*B
*/
   int i, i1, i2, j, k, n, nblocks;
   double sig1;

   n = X->n;  /* Get the order of X. */

   nblocks = X->nblocks;  /* Get the number of blocks in X. */

   if (strncmp(job, "ab", 2) == 0){

      /* Compute B = X*B. */

      /* Loop over the blocks. */

      for (k=0; k<nblocks; k++){

         if (X->bs[k+1] > 0) {
            i1 = abs(X->bs[k]);
            i2 = abs(X->bs[k+1]);

            /* B = (I - vk*vk^T)*B. */

            hhp(i1, i2, ncols, X->v, B, tdb);
        
            /* B = S*B. */

            for (i=i1; i<i2; i++)
                for (j=0; j<ncols; j++)
                    BB(i,j) = X->sig[i]*BB(i,j);

            /* B = (I - uk*uk^T)*B. */

            hhp(i1, i2, ncols, X->u, B, tdb);

         }
      }
   }

   else if (strncmp(job, "atb", 3) == 0){

      /* Compute B = X^T*B. */

      /* Loop over the blocks. */

      for (k=0; k<nblocks; k++){

         if (X->bs[k+1] > 0) {
   
            i1 = abs(X->bs[k]);
            i2 = abs(X->bs[k+1]);

            /* B = (I - uk*uk^T)*B. */
  
            hhp(i1, i2, ncols, X->u, B, tdb);
         
            /* B = Sk*B. */

            for (i=i1; i<i2; i++)
               for (j=0; j<ncols; j++)
                  BB(i,j) = (X->sig[i])*BB(i,j);
        
            /* B = (I - vk*vk^T)*B. */

            hhp(i1, i2, ncols, X->v, B, tdb);
         }
      }
   }

   else if (strncmp(job, "aib", 3) == 0){

      /* Compute B = X^-1*B. */

      /* Loop over the blocks. */

      for (k=0; k<nblocks; k++){

         if (X->bs[k+1] > 0) {
            i1 = abs(X->bs[k]);
            i2   = abs(X->bs[k+1]); 

            /* B = (I - uk*uk^T)*B. */

            hhp(i1, i2, ncols, X->u, B, tdb);
        
            /* B = Sk^{-1}*B. */

            for (i=i1; i<i2; i++){
               sig1 = 1.0/X->sig[i];
               for (j=0; j<ncols; j++)
                  BB(i,j) = sig1*BB(i,j);
            }
        
            /* B = (I - vk*vk^T)*B. */

            hhp(i1, i2, ncols, X->v, B, tdb);
         }
      }
   }

   else if (strncmp(job, "aitb", 4) == 0){

      /* Compute B = X^-T*B. */

      /* Loop over the blocks. */

      for (k=0; k<nblocks; k++){

         if (X->bs[k+1] > 0) {
            i1 = abs(X->bs[k]);
            i2   = abs(X->bs[k+1]);

            /* B = (I - vk*vk^T)*B. */

            hhp(i1, i2, ncols, X->v, B, tdb);
        
            /* B = S^{-1}*B. */

            for (i=i1; i<i2; i++) {
               sig1 = 1.0/X->sig[i];
               for (j=0; j<ncols; j++)
                   BB(i,j) = sig1*BB(i,j);
            }
        
            /* B = (I - uk*uk^T)*B. */

            hhp(i1, i2, ncols, X->u, B, tdb);

         }
      }
   }
   else{
      fprintf(stderr, "Error in HsvdProd: Illegal operation.");
      exit(0);
   }
}

void hhp(int i1, int i2, int ncols, double *w, double *B, int tdb){

/* 
   hhp computes the product of a block Householder transformation
   in a hsvdmat with a matrix contained in B.  Specifically,

      B[i1:i2,1:ncols] = (I - w[i1:i2]*w[i1:i2]^T)*B[i1:i2,1:ncols].

   i1   (in)     Points to the beginning of the block.
   i2   (in)     Points to the end of the block.
   ncol (in)     The number of columns of B
   w    (in)     The vector from the hsvd.
   B    (inout)  The array containing the matrix B.
   tdb  (in)     The trailing dimension of B.
*/

  int i, j;
  double dot;

  for (j=0; j<ncols; j++){
    dot = 0.0;
    for (i=i1; i<i2; i++)
       dot = dot + w[i]*BB(i,j);
    for (i=i1; i<i2; i++)
       BB(i,j) = BB(i,j) - w[i]*dot;
  }

/*
  This variation, which is not ansi C, is row oriented and may
  improve performance.

  int i, j;
  double dot[ncols];

  for (j=0; j<ncols; j++)
    dot[j] = 0.0;

  for (i=i1; i<i2; i++)
    for (j=0; j<ncols; j++) 
      dot[j] = dot[j] + w[i]*BB(i,j);

  for (i=i1; i<i2; i++) 
    for (j=0; j<ncols; j++)
      BB(i,j) = BB(i,j) - w[i]*dot[j];
*/

}


void EigenmatProd(struct eigenmat *A, int ncols,
                  double *B, int tdb,
                  double *C, int tdc,
                  double shift, char *job){
/* 
   EigenmatProd computes the product, perhaps shifted and inverted, of
   the eigenmat A and a matrix B, placing the result in C.

   A        Pointer to the eigenmat.
   ncols    Number of columns in the matrix B.
   B        Pointer to an array containing the matrix B.
   tdb      The trailing dimension of the array B.
   C        Pointer to an array containing the matrix C.
   tdc      The trailing dimension of C.
   shift    A shift.
   job      A string specifying the operation to be performed.

            "ab"    C = (A - shift*I)*B
            "atb"   C = (A - shift*I)^T*B
            "aib"   C = (A - shift*I)^-1*B
            "aitb"  C = (A - shift*I)^-T*B
*/

   int i, j, k, n, endj;
   double lam, lam1, mu, mult, nu, temp, u11, u12, u22;

   n = A->n;  /* Get the order of the matrix. */

   /* Copy B into C. */

   for (i=0; i<n; i++)
      for (j=0; j<ncols; j++)
         CC(i,j) = BB(i,j);

   if (strcmp(job, "ab") == 0){

      /* Compute C = AB  = Y*Z*L*Z^-1*Y^-1. */


      /* Compute Z^-1*Y^-2*C. */

      HsvdProd(&(A->Y), ncols, C, tdc, "aib");
      HsvdProd(&(A->Z), ncols, C, tdc, "aib");

      /*
         Compute L*C.  The index i points to successive blocks
         of L.
      */

      i = 0;
      while (i < n){
         if (A->type[i] == 1){
                  
            /* 1x1 block. */
         
            for (j=0; j<ncols; j++)
               CC(i,j) = (A->eig[i] - shift)*CC(i,j);
            i++;
         }    
         else if (i == n-1){
            fprintf(stderr,
                    "Error in EigenmatProd: 2x2 block starts at bs(n).\n");
            exit(0);
         }
         else if (A->type[i]==2 && A->type[i+1]==3){
         
            /* A 2x2 block. */
         
            mu = A->eig[i] - shift;
            nu = A->eig[i+1];
            for (j=0; j<ncols; j++){
                temp = mu*CC(i,j) + nu*CC(i+1,j);
                CC(i+1,j) = -nu*CC(i,j) + mu*CC(i+1,j);
                CC(i,j) = temp;
            }
            i = i+2;
         }

	 else if (A->type[i]<-1){
         
	    /* A Jordan block. */

	   endj = i+abs(A->type[i])-1;
               
	   if (endj > A->n) {
	     fprintf(stderr,
		     "Error in EigenmatProd: Jordan block too large.\n");
	     exit(0);
	   }

	   mu = A->eig[i]-shift;
	   for (j=i; j < endj; j++){
	      if (A->type[j+1] != -1) {
		 fprintf(stderr,
		     "Error in EigenmatProd: Illegal type.\n");
		 exit(0);
	      }
	      for (k=0; k<ncols; k++) {
		CC(j,k) = mu*CC(j,k) + (A->eig[j+1])*CC(j+1,k);
	      }	      
	   }

	   for (k=0; k<ncols; k++) {
	     CC(endj,k) = mu*CC(endj,k);
	   }

	   i = endj + 1;

         }
         else{
            fprintf(stderr,
                   "Error in EigenmatProd: Illegal type.\n");
            exit(0);
         }
      }

         /* Compute Y*Z*C. */

      HsvdProd(&(A->Z), ncols, C, tdc, "ab");
      HsvdProd(&(A->Y), ncols, C, tdc, "ab");
   }
   else if (strcmp(job, "atb") == 0){

      /* Compute C = A^T*B = (Y^-T*Z^-T*L^T*Z^T*Y^T)*C. */

      /* compute Z^T*Y^T*C. */

      HsvdProd(&(A->Y), ncols, C, tdc, "atb");
      HsvdProd(&(A->Z), ncols, C, tdc, "atb");
      /*
         Compute L*C.  The index i points to successive blocks
         of L.
      */
      i = 0;
      while (i<n){
         if (A->type[i] == 1){

            /* 1x1 block. */

            lam = A->eig[i] - shift;
            for (j=0; j<ncols; j++)
               CC(i,j) = lam*CC(i,j);
            i++;
         }
         else if (i == n-1){
            fprintf(stderr,
                    "Error in EigenmatProd: 2x2 block starts at bs(n).\n");
            exit(0);
         }
         else if (A->type[i]==2 && A->type[i+1]==3){

            /* A 2x2 block. */

            mu = A->eig[i] - shift;
            nu = A->eig[i+1];

            for (j=0; j<ncols; j++){
               temp = mu*CC(i,j) - nu*CC(i+1,j);
               CC(i+1,j) = nu*CC(i,j) + mu*CC(i+1,j);
               CC(i,j) = temp;
            }

            i = i+2;

         } else if (A->type[i]<-1){
         
	    /* A Jordan block. */

	   endj = i+abs(A->type[i])-1;
               
	   if (endj > A->n) {
	     fprintf(stderr,
		     "Error in EigenmatProd: Jordan block too large.\n");
	     exit(0);
	   }

	   mu = A->eig[i]-shift;

	   for (j=endj; j > i; j--){
	      if (A->type[j] != -1) {
		 fprintf(stderr,
		     "Error in EigenmatProd: Illegal type.\n");
		 exit(0);
	      }
	      for (k=0; k<ncols; k++) {
		CC(j,k) = mu*CC(j,k) + (A->eig[j])*CC(j-1,k);
	      }	      
	   }

	   for (k=0; k<ncols; k++) {
	     CC(i,k) = mu*CC(i,k);
	   }

	   i = endj + 1;

         }

         else{
            fprintf(stderr,
                   "Error in EigenmatProd: Illegal type.\n");
            exit(0);
         }
      }

      /* Compute C = Y^T*Z^T*C. */

      HsvdProd(&(A->Z), ncols, C, tdc, "aitb");
      HsvdProd(&(A->Y), ncols, C, tdc, "aitb");
   }

   else if (strcmp(job, "aib") == 0){

      /* Compute C = A^-1*C = (Y*Z*L^-1*Y^-1*Z^-1)*C. */

      /* Compute Z^-1*Y^-2*C. */

      HsvdProd(&(A->Y), ncols, C, tdc, "aib");
      HsvdProd(&(A->Z), ncols, C, tdc, "aib");

      /*
         Compute L*C.  The index i points to successive blocks
         of L.
      */

      i = 0;
      while (i<n){
         if (A->type[i] == 1){

            /* 1x1 block. */
  
            lam1 = 1.0/(A->eig[i] - shift);
            for (j=0; j<ncols; j++)
               CC(i,j) = lam1*CC(i,j);
            i++;
         }
         else if (i == n-1){
            fprintf(stderr,
                    "Error in EigenmatProd: 2x2 block starts at bs(n).\n");
            exit(0);
         }
         else if (A->type[i]==2 && A->type[i+1]==3){

            /*
               A 2x2 block.  Compute The inverse of the block
               times C using Gaussian elimination with partial
               pivoting.
            */

            mu = A->eig[i] - shift;
            nu = A->eig[i+1];
            if (fabs(mu) >= fabs(nu)){

               /* No pivoting. */

               u11 = mu;
               u12 = nu;
               mult = -nu/mu;
               u22 = mu - mult*nu;
               for (j=0; j<ncols; j++){
                  CC(i+1,j) = (CC(i+1,j) - mult*CC(i,j))/u22;
                  CC(i,j) = (CC(i,j) - u12*CC(i+1,j))/u11;
               }
            }
            else{

               /* Pivot. */

               u11 = -nu;
               u12 = mu;
               mult = -mu/nu;
               u22 = nu - mult*mu;
               for (j=0; j<ncols; j++){
                  temp = CC(i+1,j);
                  CC(i+1,j) = (CC(i,j) - mult*temp)/u22;
                  CC(i,j) = (temp - u12*CC(i+1,j))/u11;
               }
            }
            i = i+2;

         } 
	 else if (A->type[i]<-1){
         
	    /* A Jordan block. */

	   endj = i+abs(A->type[i])-1;
               
	   if (endj > A->n) {
	     fprintf(stderr,
		     "Error in EigenmatProd: Jordan block too large.\n");
	     exit(0);
	   }

	   mu = A->eig[i]-shift;

	   for (k=0; k<ncols; k++) {
	     CC(endj,k) = CC(endj,k)/mu;
	   }

	   for (j=endj-1; j>=i; j--){
	      if (A->type[j+1] != -1) {
		 fprintf(stderr,
		     "Error in EigenmatProd: Illegal type.\n");
		 exit(0);
	      }
	      for (k=0; k<ncols; k++) {
		CC(j,k) = (CC(j,k) - (A->eig[j+1])*CC(j+1,k))/mu;
	      }	      
	   }

	   i = endj + 1;

         }
         else{
            fprintf(stderr,
                   "Error in EigenmatProd: Illegal type.\n");
            exit(0);
         }
      }

      /* Compute C = Y*Z*C. */

      HsvdProd(&(A->Z), ncols, C, tdc, "ab");
      HsvdProd(&(A->Y), ncols, C, tdc, "ab");
   }

   else if (strcmp(job, "aitb") == 0){

      /* Compute C = A^-T*C = (Y^-T*Z^-T*L^-T*Z^T*Y^T)*C. */

      /* Compute C = Z^T*Y^T*C. */

      HsvdProd(&(A->Y), ncols, C, tdc, "atb");
      HsvdProd(&(A->Z), ncols, C, tdc, "atb");

      /*
         Compute L*C.  The index i points to successive blocks
         of L.
      */

      i = 0;
      while (i<n){
         if (A->type[i] == 1){

            /* 1x1 block. */

            lam1 = 1.0/(A->eig[i] - shift);
            for (j=0; j<ncols; j++)
               CC(i,j) = lam1*CC(i,j);
            i++;
         }
         else if (i == n-1){
            fprintf(stderr,
                    "Error in EigenmatProd: 2x2 block starts at bs(n).\n");
            exit(0);
         }
         else if (A->type[i]==2 && A->type[i+1]==3){
   
            /*
               A 2x2 block.  Compute The inverse of the block
               times C using Gaussian elimination with partial
               pivoting.
            */

            mu = A->eig[i] - shift;
            nu = A->eig[i+1];
            if (fabs(mu) >= fabs(nu)){
      
               /* No pivoting. */

               u11 = mu;
               u12 = -nu;
               mult = nu/mu;
               u22 = mu + mult*nu;
               for (j=0; j<ncols; j++){
                  CC(i+1,j) = (CC(i+1,j) - mult*CC(i,j))/u22;
                  CC(i,j) = (CC(i,j) - u12*CC(i+1,j))/u11;
               }
            }
            else{

               /* Pivot. */

               u11 = nu;
               u12 = mu;
               mult = mu/nu;
               u22 = -nu - mult*mu;
               for (j=0; j<ncols; j++){
                  temp = CC(i+1,j);
                  CC(i+1,j) = (CC(i,j) - mult*temp)/u22;
                  CC(i,j) = (temp - u12*CC(i+1,j))/u11;
               }
            }
            i = i+2;
         }

	 else if (A->type[i]<-1){
         
	    /* A Jordan block. */

	   endj = i+abs(A->type[i])-1;
               
	   if (endj > A->n) {
	     fprintf(stderr,
		     "Error in EigenmatProd: Jordan block too large.\n");
	     exit(0);
	   }

	   mu = A->eig[i]-shift;

	   for (k=0; k<ncols; k++) {
	     CC(i,k) = CC(i,k)/mu;
	   }

	   for (j=i+1; j<=endj; j++){
	      if (A->type[j] != -1) {
		 fprintf(stderr,
		     "Error in EigenmatProd: Illegal type.\n");
		 exit(0);
	      }
	      for (k=0; k<ncols; k++) {
		CC(j,k) = (CC(j,k) - (A->eig[j])*CC(j-1,k))/mu;
	      }	      
	   }

	   i = endj + 1;

         }

         else{
            fprintf(stderr,
                   "Error in EigenmatProd: Illegal type.\n");
            exit(0);
         }
      }

      /* Compute Y^-T*Z^-T*C. */

      HsvdProd(&(A->Z), ncols, C, tdc, "aitb");
      HsvdProd(&(A->Y), ncols, C, tdc, "aitb");
   }
   else{
      fprintf(stderr,
              "Error in EigenmatProd: Illegal operation.\n");
      exit(0);
   }
}

void EigenmatVecs(struct eigenmat *A, int eignum,
                  double *eigr, double *eigi, 
                  double *xr, double *xi,
                  double *yr, double *yi,
                  double *cond, char job){

/*
      EigenmatVecs computes an eigenvalue and the corresponding left or
      right eigenvectors or principal vectors as specified job.  If
      both are computed, EigenmatVecs also returns the condition
      number of the eigenvalue.  The eigenvectors are scaled to have
      Euclidean norm 1.
 
      A        The eigenmat whose vectors are to be computed.
      eignum   The position in A.eig of the eigenvalue.
      eigr     The real part of the eigenvalue.
      eigi     The imaginary part of the eigenvalue.
      xr(:)    The real part of the right eigenvector or
               principal vector.
      xi(:)    The imaginary part of the right eigenvector or
               principal vector.
      yr(:)    The real part of the left eigenvector or
               principal vector.
      yi(:)    The imaginary part of the left eigenvector or
               principla vector.
      cond     The condition number of the eigenvalue
               (or -1, if the eigenvalue belongs to a
               Jordan block).
      job      A string specifying what to compute.
 
               "r"  Compute the right eigenvector or
                    principal vector.
               "l"  Compute the left eigenvector or
                    principal vector.
               "b"  Compute both and the condition number.
               (Note: For Jordan blocks, principal vectors
               are computed and -1 is returned for the
               condition number.)
*/

   int i, found;
   int j = eignum;
   double cr, ci, nrm2;

   if (eignum > A->n) {
      fprintf(stderr,
               "Error in EigenmatVecs: eignum > n.\n");
      exit(0);
   }

   if (job!='b' && job!='r' && job!='l'){
      fprintf(stderr,
              "Error in EigenmatVecs: Illegal option. \n");
      exit(0);
   }

   if (job == 'r' || job == 'b') {
     for (i=0; i < A->n; i++) {
       xr[i] = 0.0;
       xi[i] = 0.0;
     }
   } else if (job == 'l' || job == 'b') {
     for (i=0; i < A->n; i++) {
       yr[i] = 0.0;
       yi[i] = 0.0;  
     }   
   }
   
   if (job == 'b')
     *cond = 0.0;

   cr = 0.0;
   ci = 0.0;
   
   if (A->type[j] == 1 || A->type[j] <= -1) {

      /* Real eigenvalue. */
      
      *eigi = 0.0;

      if (A->type[j] == 1 || A->type[j] < -1) { 
	 *eigr = A->eig[j];
      }
      else if (A->type[j] == -1) {
         found = 0;
         for (i=j-1; i>=0; i--) {
	    if (A->type[i] < -1){
               *eigr = A->eig[i];
	       found = 1;
	       break;
	    }
	    else if (A->type[i] != -1) {
 	       fprintf(stderr,
		       "Error in EigenmatVecs: Illegal type in Jordan block.\n");
	       exit(0);
	    }
	 }
	 if (found == 0){
	    fprintf(stderr,
		    "Error in EigenmatVecs: Illegal type in Jordan block.\n");
	    exit(0);
	 }
      }
      else {
  	 fprintf(stderr,
		 "Error in EigenmatVecs: Illegal type in Jordan block.\n");
	 exit(0);	
      }

      if (job == 'r' || job == 'b') {

         /* Compute the right eigenvector.*/

         xr[j] = 1.0;
         HsvdProd(&(A->Z), 1, xr, 1, "ab");
         HsvdProd(&(A->Y), 1, xr, 1, "ab");
	 if (A->type[j] == 1) {
	   nrm2 = 0.0;
	   for (i=0; i<A->n; i++)
	     nrm2 = nrm2 + xr[i]*xr[i];
	   nrm2 = sqrt(nrm2);
	   for (i=0; i<A->n; i++)
	     xr[i] = xr[i]/nrm2;
	 }
      }

      if (job == 'l' || job == 'b') {

         /* Compute the left eigenvector.*/

         yr[j] = 1.0;
         HsvdProd(&(A->Z), 1, yr, 1, "aitb");
         HsvdProd(&(A->Y), 1, yr, 1, "aitb");

	 if (A->type[j] == 1) {
	   nrm2 = 0.0;
	   for (i=0; i<A->n; i++)
	     nrm2 = nrm2 + yr[i]*yr[i];
	   nrm2 = sqrt(nrm2);
	   for (i=0; i<A->n; i++)
            yr[i] = yr[i]/nrm2;
	 }
      }

      if (job == 'b') {

         /* Compute the condition number. */

 	 if (A->type[j] == 1) {
	    for (i=0; i<A->n; i++)
	      *cond = *cond + yr[i]*xr[i];
	    *cond = 1.0/(*cond);
	 } 
	 else {
	    *cond = -1.0;
	 }
      }

   } 
   else if (j == A->n) {
      fprintf(stderr,
              "Error in EigenmatVecs: Type error.\n");
      exit(0);
   }

   else if (A->type[j] == 2) {
      if (j >= A->n || A->type[j+1] != 3 ) {
         fprintf(stderr,
                 "Error in EigenmatVecs: Type error.\n");
         exit(0);
      }

      /* Complex eigenvalue. */

      *eigr = A->eig[j];
      *eigi = A->eig[j+1];

      /* Compute the right eigenvector. */

      if (job == 'r' || job == 'b') {
         xr[j] = M_SQRT1_2;
         xi[j+1] = M_SQRT1_2;
         HsvdProd(&(A->Z), 1, xr, 1, "ab");
         HsvdProd(&(A->Y), 1, xr, 1, "ab");
         HsvdProd(&(A->Z), 1, xi, 1, "ab");
         HsvdProd(&(A->Y), 1, xi, 1, "ab");
         nrm2 = 0.0;
         for (i=0; i<A->n; i++)
            nrm2 = nrm2 + xr[i]*xr[i] + xi[i]*xi[i];
         nrm2 = sqrt(nrm2);
         for (i=0; i<A->n; i++){
            xr[i] = xr[i]/nrm2;
            xi[i] = xi[i]/nrm2;
         }
      }

      if (job == 'l' || job == 'b') {

         /* Compute the right eigenvector. */

         yr[j] = M_SQRT1_2;
         yi[j+1] = M_SQRT1_2;
         HsvdProd(&(A->Z), 1, yr, 1, "aitb");
         HsvdProd(&(A->Y), 1, yr, 1, "aitb");
         HsvdProd(&(A->Z), 1, yi, 1, "aitb");
         HsvdProd(&(A->Y), 1, yi, 1, "aitb");
         nrm2 = 0.0;
         for (i=0; i<A->n; i++)
            nrm2 = nrm2 + yr[i]*yr[i] + yi[i]*yi[i];
         nrm2 = sqrt(nrm2);
         for (i=0; i<A->n; i++){
            yr[i] = yr[i]/nrm2;
            yi[i] = yi[i]/nrm2;
         }
      }

      /* Compute the condition number. */

      if (job == 'b') {
         for(i=0; i<A->n; i++){
            cr = cr + xr[i]*yr[i]+xi[i]*yi[i];
            ci = ci - xr[i]*yi[i]+xi[i]*yr[i];
         }
         *cond = 1.0/sqrt(cr*cr+ci*ci);
      }      

   }
   else {
      fprintf(stderr,
              "Error in EigenmatVecs: Type error %d.\n", A->type[j]);
      exit(0);
   }
}

void hscal(int n, double *u) {

/* 
   hscal scales vector u such that norm(u) = sqrt(2) 
   n is the length of vector u.
*/

  double max=0.0, tmp, nrm = 0.0;
  int i;

  /* Find the maximum magnitude element */

  max = fabs(u[0]);
  for(i=1; i<n; i++){
    tmp = fabs(u[i]);
    if (max<tmp) {
      max = tmp;
    }
  }

  /* Compute nrm = u'u/max^2 */

  for(i=0; i<n; i++){
    tmp = u[i]/max;
    nrm += tmp*tmp;
  }

  if (nrm == 0.0)  {
    fprintf(stderr,
	    "Error in hnorm: zero argument.\n");
    exit(0);
  }

  /* Scale u */

  nrm = M_SQRT2/sqrt(nrm)/max;

  for(i=0; i<n; i++){
    u[i] = nrm*u[i];
  }
}
