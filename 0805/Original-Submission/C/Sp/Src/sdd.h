/*----------------------------------------------------------------------
  
  SDDPACK: Software for the Semidiscrete Decomposition.
  Copyright (c) 1999 Tamara G. Kolda and Dianne P. O'Leary. 

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

  MODULE: sdd.o
  
  DESCRIPTION: Header file for sdd.c, routines for creating and
  updating the SDD decomposition and sparse matrices.
  
  MODULE DEPENDENCIES: sdd.c sdd.h 
  
  BUG REPORTS: Email Tamara.Kolda@na-net.ornl.gov.
  
  ----------------------------------------------------------------------*/

#ifndef _SDD_H
#define _SDD_H

/*----------------------------------------------------------------------
  Includes
  ----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------
  Definitions
  ----------------------------------------------------------------------*/

#define MAXLINE 1000		/* Max length of comment line in data file. */

/* WARNING: The IDXSHIFT definition only works on 32 or 64 bit
   architectures.  If you have a different architecture (e.g., 128
   bit), you should replace the next line with the appropriate value
   (e.g., #define IDXSHIFT 7). */

#ifndef IDXSHIFT
#define IDXSHIFT ((sizeof(ulong) == 4) ? 5 : 6)	/* log2(bits_per_word) */
#endif

#define IDXMASK ((sizeof(ulong) * 8) - 1) /* bits_per_word - 1 */

#define ONE ((ulong)(1))	/* one */

#define MAXMASK (ONE << IDXMASK) /* one in leftmost position */

/*----------------------------------------------------------------------
  Type Definitions
  ----------------------------------------------------------------------*/

/* Stores a sparse m x n matrix in compressed sparse column format.
   jc is an (n+1)-long array. jc[j] gives the starting index in val
   and ir for the j-th column, and jc[j+1] gives the ending index +
   1. */
typedef struct matrix_struct {
  int m;			/* number of rows */
  int n;			/* number of columns */
  int nnz;			/* number of nonzero entries */
  int *jc;			/* array of length (n+1), col starts */
  int *ir;			/* array of length nnz, row indices */
  float *val;			/* array of length nnz, values */
} matrix;

/* Used for packed bit vectors */
typedef unsigned long ulong;	

/* The type sddfloat is used as the final precision that the diagonal
   values of the matrix D are stored in.  The type sdddouble is used for
   all calculations.  Experience has shown us that float is sufficient
   for our calculations and requires less computation time.  */
typedef float sddfloat;		/* used for low-precision reals  */
typedef float sdddouble;	/* used for high-precision reals */

/* The following structure is used for sorting.  The function qsortopt
   in sdd.c is specially modified for this datatype and is used if
   QSORTOPT is defined. */
typedef struct sdddouble_plus_struct {
  sdddouble val;
  int idx;
} sdddouble_plus;

/* An s-vector is used to store an m-long array of s-values.  Each
   s-value is represented by a val bit and a sgn bit.
            0: val = 0, sgn = undefined
            1: val = 1, sgn = 0
           -1: val = 1, val = 1
   These bits are packed into ulong's to form the val and sgn
   arrays. */
typedef struct svector_struct {
  int m;			/* number of entries stored */
  ulong *val;			/* zero (0) or non-zero (1) */
  ulong *sgn;			/* plus (0) or minus (1) */
} svector;

/* An s-matrix is used to store up to an m by kmax matrix of s-values.
   The dimension of the stored matrix is m by k.  The col array stores
   kmax pointers to s-vectors which are the columns of the matrix.
   Specifically, col[i] is a pointer to the (i+1)st column of the
   matrix (0 <= i < k).  */
typedef struct smatrix_struct {
  int k, kmax;
  int m;
  svector **col;
} smatrix;

/* A dmatrix is used to store the diagonal entries of up to a kmax by
   kmax diagonal matrix.  The dimension of the stored matrix is given
   by k. d is an array of length kmax such that d[i] is the (i+1)st
   diagonal value (0 <= i < k).  */
typedef struct dmatrix_struct {
  int k, kmax;
  sddfloat *d;
} dmatrix;

/* An sdd stores a semi-discrete decomposition (SDD) to an m by n
   matrix up to kmax terms.  The current number of terms is given by
   k.  X and Y are pointers to the left and right smatrices
   respectively, and D is a pointer to the diagonal matrix.  */
typedef struct sdd_struct {
  int m, n, k, kmax;
  smatrix *X, *Y;
  dmatrix *D;
} sdd;

/*----------------------------------------------------------------------
  Function Declarations
  ----------------------------------------------------------------------*/

void free_sdd(sdd*);
sdd* read_sdd(char*, int);
void write_sdd(sdd*, char*, int, char*);

void free_matrix(matrix*);
matrix* read_matrix(char*, int);
void write_matrix(matrix*, char*, int);

sdd* compute_sdd (matrix*, sdd*, int, int, float, int, float, int);
sdd* update_sdd(matrix*, sdd*, int, int);

#endif





