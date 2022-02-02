/*
 *  gts2srf.c:  Surface file conversion utility.
 *
 *  R. Renka
 *  renka@cs.unt.edu
 *  11/18/03
 *
 *  This program reads in a data set defining a surface
 *  in the format used by the Gnu Triangulated Surface
 *  (gts) program, and writes out a file suitable for
 *  input to glsurf.
 *
 */

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
   unsigned int i, i1, i2, i3, i11, i12, i21, i22, i31, i32,
                ie1, ie2, ie3, k;
   unsigned int ne, nt, nv;
   unsigned int *inde, *indt;
   double *vertices;
   FILE *fpin, *fpout;
   char *filemsg = "\n\n"
     "  This program must be called with two arguments:\n"
     "  a path to a gts-formatted file followed by a file\n"
     "  name for a glsurf-formatted file.\n\n";
   char *memerror = "  Unable to allocate sufficient memory.\n";
   char *readerror = "  Error reading input file.\n";
   char *writerror = "  Error writing output file.\n";
/*
 *  Print a message.
 */
   printf( "\n\nRunning surface file conversion utility.\n\n" );
/*
 *  Test for file names, and open files.
 */
   if (argc < 3) {
      printf("%s", filemsg);
      exit(1);
   }
   fpin = fopen(argv[1], "r");
   if (fpin == NULL) {
     printf("Cannot open input file %s.\n", argv[1]);
     exit(1);
   }
   fpout = fopen(argv[2], "w");
   if (fpout == NULL) {
      printf("Cannot open output file %s.\n", argv[2]);
      exit(1);
   }
/*
 *  Read counts nv, ne, and nt (vertices, edges, and
 *  triangles).
 */
   if (fscanf(fpin, "%u", &nv) != 1) {
      printf("%s", readerror);
      exit(1);
   }
   if (nv < 3) {
      printf("Invalid value:  nv = %d.\n", nv);
      exit(1);
   }
   if (fscanf(fpin, "%u", &ne) != 1) {
      printf("%s", readerror);
      exit(1);
   }
   if (ne < 3) {
      printf("Invalid value:  ne = %u.\n", ne);
      exit(1);
   }
   if (fscanf(fpin, "%u", &nt) != 1) {
      printf("%s", readerror);
      exit(1);
   }
   if (nt < 1) {
      printf("Invalid value:  nt = %u.\n", nt);
      exit(1);
   }
/*
 *  Allocate storage and read vertices.
 */
   vertices = (double *) calloc (3*nv, sizeof(double));
   if (vertices == NULL) {
      printf("%s", memerror);
      exit(1);
   }
   for (i = 0; i < nv; i++) {
      if (fscanf(fpin, "%lf%lf%lf", vertices+3*i,
                 vertices+3*i+1, vertices+3*i+2) != 3) {
         printf("%s", readerror);
         exit(1);
      }
   }
/*
 *  Allocate storage and read edge vertex indices inde.
 */
   inde = (unsigned int *) calloc(2*ne, sizeof(unsigned int));
   if (inde == NULL) {
      printf("%s", memerror);
      exit(1);
   }
   for (i = 0; i < ne; i++) {
      if (fscanf(fpin,"%u%u",inde+2*i,inde+2*i+1) != 2) {
         printf("%s", readerror);
         exit(1);
      }
      if (inde[2*i  ] < 1 || inde[2*i  ] > nv ||
          inde[2*i+1] < 1 || inde[2*i+1] > nv) {
         printf("The i-th edge index pair is invalid for "
                "i = %u.\n", i);
         exit(1);
      }
   }
/*
 *  Allocate storage and read triangle edge indices indt.
 */
   indt = (unsigned int *) calloc(3*nt, sizeof(unsigned int));
   if(indt == NULL) {
      printf("%s", memerror);
      exit(1);
   }
   for (k = 0; k < nt; k++) {
      if (fscanf(fpin,"%u%u%u",indt+3*k,indt+3*k+1,
                 indt+3*k+2) != 3) {
         printf("%s", readerror);
         exit(1);
      }
      if (indt[3*k  ] < 1 || indt[3*k  ] > ne ||
          indt[3*k+1] < 1 || indt[3*k+1] > ne ||
          indt[3*k+2] < 1 || indt[3*k+2] > ne) {
         printf("The k-th triangle triple is invalid for "
                "k = %u.\n", k);
         exit(1);
      }
   }
/*
 *  Convert indt to triangle vertex indices.
 */
   for (k = 0; k < nt; k++) {
      ie1 = indt[3*k]-1;
      ie2 = indt[3*k+1]-1;
      ie3 = indt[3*k+2]-1;
      i11 = inde[2*ie1];
      i12 = inde[2*ie1+1];
      i21 = inde[2*ie2];
      i22 = inde[2*ie2+1];
      i31 = inde[2*ie3];
      i32 = inde[2*ie3+1];
      if (i21 == i11) {
         i1 = i12;
         i2 = i11;
         i3 = i22;
      } else if (i21 == i12) {
         i1 = i11;
         i2 = i12;
         i3 = i22;
      } else if (i22 == i11) {
         i1 = i12;
         i2 = i11;
         i3 = i21;
      } else if (i22 == i12) {
         i1 = i11;
         i2 = i12;
         i3 = i21;
      } else {
         printf("Triangle %u is invalid.\n", k);
         exit(1);
      }
      if (i1 == i2 || i1 == i3 || i2 == i3) {
         printf("Triangle %u has duplicate vertices.\n", k);
         exit(1);
      }
      if ((i31 != i1  &&  i31 != i3) ||
          (i32 != i1  &&  i32 != i3)) {
         printf("Triangle %u is not valid.\n", k);
         exit(1);
      }
/*
 *  Overwrite the edge indices by vertex indices.
 */
      indt[3*k] = i1;
      indt[3*k+1] = i2;
      indt[3*k+2] = i3;
   }
/*
 *  Write the output file.
 */
   if (fprintf(fpout, "%u\n", nv) < 0) {
      printf("%s", writerror);
      exit(1);
   }
   for (i = 0; i < nv; i++) {
      if (fprintf(fpout, "% 24.16e % 24.16e % 24.16e\n",
          vertices[3*i], vertices[3*i+1], vertices[3*i+2]) < 0) {
         printf("%s", writerror);
         exit(1);
      }
   }
   if (fprintf(fpout, "%u\n", nt) < 0) {
      printf("%s", writerror);
      exit(1);
   }
   for (k = 0; k < nt; k++) {
      if (fprintf(fpout,"%u %u %u\n",indt[3*k],indt[3*k+1],
                  indt[3*k+2]) < 0) {
         printf("%s", writerror);
         exit(1);
      }
   }
   fclose(fpin);
   fclose(fpout);
   return 0;
}
