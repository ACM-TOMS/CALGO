#define NULL 0

/*
   Definitions to overcome the limitations of array passing in C.
   If B pooints to an n x tdb array, then BB(i,j) computes B[i][j].
   Similarly for CC.
*/

#define BB(I,J)  (*(B + tdb*(I) + J))
#define CC(I,J)  (*(C + tdc*(I) + J))

/*
  A hsvd matrix has the form

      X = diag(X1, X2, ..., Xnblocks),

   where

      Xk = (I - uk*uk^T)*Sigmak*(I - vk*vk^t).

   with Sk diagonal having positive diagonal entries.  The vectors,
   uk, vk, and the diagonals of Sk are stacked in order in the arrays
   X.u, X.v, and X.sig.  abs(X.bs[i]) points to the beginning of the
   ith block and abs(X.bs[X.nblocks-1]) = n.  If X.bs[i+1]<0, the
   ith block is an identity matrix.
*/

struct hsvd{
   int n;          /* The order of the matrix. */
   int nblocks;    /* The number of blocks in the hsvdmat. */
   int *bs;        /* abs(bs[i]) is the index of the start of the
                      block i+1.  abs(bs[nblock]) is equal to n.
                      If bs[i+1]<0, the i-th block is an identity. */
   double *u;      /* The vectors generating the left Householder
                      transformation. */
   double *v;      /* The vectors generating the right Householder
                      transformation. */
   double *sig;    /* The diagonals of the Sk. */
};

/*
    An eigenmat has the form

        A = Y*Z*(L - shift*I)*Z^-1*Y^-1.

    Here L is a block triangular with 1x1 blocks containing the
    real eigenvalues of A, 2x2 blocks conaining the complex
    eignevalues of A in the form

                 |  mu  nu |
                 |         |,
                 | -nu  mu |

    and Jordon blocks of order k having the form (for k=3)

                 | lambda     eta_1      0      |
                 |                              |
                 |   0       lambda     eta_2   |.
                 |                              |
                 |   0         0       lambda   |

   Y and Z are hsvd transformations, with Y having only one block.
*/

struct eigenmat{
   int n;            /* The order of the matrix. */
   double *eig;      /* Array containing the eigenvalues of A
                        or the superdiagonals of a Jordan block. */
   int    *type;     /* Types according to the following scheme.
                        Real eigenvalue
                           type[i]=1      eig[i]=the eigenvalue.
                        Complex conjugate eigenvalues
                           type[i]=2      eig[i]=mu
                           type[i]=3      eig[i]=nu
                        Jordan block of order k
                           type[i]=-k     eig[i]=lambda
                           type[i+j]=-1   eig[i+j]=eta_j */

   struct hsvd Y, Z; /* The outer and inner hsvd transformations */
};
void HsvdProd(struct hsvd *X, int ncols, double *B, int tdb, char *job);
void EigenmatAlloc(struct eigenmat *A, int n, int nblocks,
		  int yident, int zident);
void EigenmatProd(struct eigenmat *A, int ncols, double *B, int tdb,
                 double *C, int tdc, double shift, char *job);
void EigenmatFree(struct eigenmat *A);
void EigenmatVecs(struct eigenmat *A, int eignum, double *eigr, double *eigi, 
		 double *xr, double *xi, double *yr, double *yi, double *cond, 
		 char job);
void hhp(int i1, int i2, int ncols, double *w, double *B, int tdb);
void hscal(int n, double *u);
