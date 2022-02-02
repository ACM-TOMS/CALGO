
/*
 *  glsdat.c:  glsurf data set creation program.
 *
 *  This program creates a data set, suitable for input to
 *  glsurf, representing the graph of bivariate function f
 *  on [xmin,xmax] X [ymin,ymax].
 *
 *  R. Renka
 *  renka@cs.unt.edu
 *
 *  11/29/03
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 80   /* Number of vertices in the x direction */
#define NY 80   /* Number of vertices in the y direction */

double f(double x, double y, unsigned int nd, double *fx,
         double *fy);


/*
 *----------------------------------------------------------
 *
 *  f:  Bivariate function from which the data set is
 *      created.  Note that this should be consistent with
 *      the domain [xmin,xmax] X [ymin,ymax] defined below.
 *
 *  On input:  x,y = Evaluation point.
 *
 *             nd = Derivative flag:
 *                  nd == 0 if only a function value is to
 *                          be returned.
 *                  nd >= 1 if first partial derivatives
 *                          are also to be computed.
 *
 *  On output:  f = Function value f(x,y).
 *
 *              fx,fy = First partials at (x,y) if nd >= 1.
 *
 *----------------------------------------------------------
 */
double f(double x, double y, unsigned int nd, double *fx,
         double *fy)
{
   double t, tx, ty, z;
   t = sqrt( (80.0*x - 40.0)*(80.0*x - 40.0) +
             (90.0*y - 45.0)*(90.0*y - 45.0) );
   z = exp(-.04*t)*cos(.15*t);
   if (nd == 0) return z;
   if (t == 0.0) {
      *fx = -3.2;
      *fy = -3.6;
   } else {
      tx = (6400.0*x - 3200.0)/t;
      ty = (8100.0*y - 4050.0)/t;
      t = -(.15*exp(-.04*t)*sin(.15*t) + .04*z);
      *fx = t*tx;
      *fy = t*ty;
   }
   return z;
}


int main(int argc, char *argv[])
{
   FILE *fp;                   /* File pointer */
   char *filemsg = "\n\n"
     "  This program must be called with one argument --\n"
     "  a file name for the output data set.  If the file\n"
     "  already exists, it will be overwritten.\n";
   char *writerror = "  Error writing to file.\n";

   unsigned int i, j;          /* Row and column indexes */
   unsigned int k;             /* Vertex index */
   unsigned int nt = 2*(NX-1)*(NY-1);  /* Number of triangles */
   unsigned int nv = NX*NY;            /* Number of vertices */

   double x, y, z;             /* Vertex components */
   double dx, dy;              /* Domain width and height */
   double fx, fy;              /* First partial derivatives */

   double xmin = -1.0;         /* Domain extrema */
   double xmax = 2.0;
   double ymin = -1.0;
   double ymax = 2.0;
/*
 *  The following parameter controls the option of including
 *  vertex normal vectors in the data set.
 */
   unsigned int nd = 0;        /* Order of derivatives */
/*
 *  Create a file with user-specified file name and write nv.
 */
   if (argc < 2) {
      printf("%s", filemsg);
      exit(1);
   }
   fp = fopen(argv[1], "w");
   if (fp == NULL) {
      printf("Cannot open file %s.\n", argv[1]);
      exit(1);
   }
   if (fprintf(fp, "%u\n", nv) < 0) {
      printf("%s", writerror);
      exit(1);
   }
/*
 *  Write vertices to fp.
 */
   dx = (xmax-xmin)/(double) (NX-1);
   dy = (ymax-ymin)/(double) (NY-1);
   for (j = 0; j < NY; j++) {
      y = ymin + (double) j*dy;
      for (i = 0; i < NX; i++) {
         x = xmin + (double) i*dx;
         z = f(x,y,0,&fx,&fy);
         fprintf(fp, "%23.15e %23.15e %23.15e\n", x, y, z);
      }
   }
/*
 *  Write nt and the triangles (CCW-ordered triples of
 *  vertex indices) to fp.  For each of the (NX-1)*(NY-1)
 *  cells, k is the vertex index of the upper right corner,
 *  and the cell is partitioned by the diagonal with slope 1.
 */
   fprintf(fp, "%u\n", nt);
   for (j = 2; j <= NY; j++) {
      for (i = 2; i <= NX; i++) {
         k = NX*(j-1) + i;
         fprintf(fp, "%6.1u %6.1u %6.1u\n", k-NX-1, k-NX, k);
         fprintf(fp, "%6.1u %6.1u %6.1u\n", k-NX-1, k, k-1);
      }
   }
   if (nd == 0) return 0;
/*
 *  Write nv and vertex unit normal vectors (-fx,-fy,1)
 *  (normalized) to fp.
 */
   fprintf(fp, "%u\n", nv);
   for (j = 0; j < NY; j++) {
      for (i = 0; i < NX; i++) {
         y = ymin + (double) j*dy;
         x = xmin + (double) i*dx;
         z = f(x,y,nd,&fx,&fy);
         z = 1.0/sqrt(1.0 + fx*fx + fy*fy);
         x = -fx*z;
         y = -fy*z;
         fprintf(fp, "%23.15e %23.15e %23.15e\n", x, y, z);
      }
   }
   return 0;
}
