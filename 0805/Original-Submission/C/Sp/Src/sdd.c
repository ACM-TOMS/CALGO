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
  
  DESCRIPTION: Routines for creating and updating SDD's and sparse
  matrices.
  
  EXTERNAL SUBROUTINES: No external subroutines
  
  MODULE DEPENDENCIES: sdd.c sdd.h
  
  BUG REPORTS: Email Tamara.Kolda@na-net.ornl.gov.
  
  ----------------------------------------------------------------------*/

#include "sdd.h"		/* header file */

/*----------------------------------------------------------------------
  Internal Definitions, Macros, and Structures
  ----------------------------------------------------------------------*/

/* Max length of comment line in data file. */
#define MAXLINE 1000

/* Information level. Higher numbers yield more info. */
#ifndef INFO			
#define INFO 10			
#endif

/* Sets every bit in the first y bytes pointed to by x to zero or one
   respectively. */
#define bzero(x, y) memset(x, 0, y) 
#define bone(x, y) memset(x, 0xff, y) 

/* Structure for reading in and sorting matrix entries from
   MatrixMarket formatted file. */
typedef struct entry_struct {  
  int i, j;			/* indicies */
  float val;			/* entry */
} entry;

/*----------------------------------------------------------------------
  Internal Function Declarations
  ----------------------------------------------------------------------*/

/* Internal svector functions */
int numwords(int);
void free_svector(svector*);
svector* create_svector(svector*, int);
void write_svector(svector*, FILE*, int);
svector* read_svector(int, FILE*, int);
int svcount(svector*);

/* Initialization functions */
sdddouble *init_threshold(svector*, matrix*, sdd*, sddfloat, int*);
void init_cycle(svector*, int);
void init_ones(svector*);
void init_pones(svector*);

/* Internal smatrix functions */
void free_smatrix(smatrix*);
smatrix* create_smatrix(smatrix*, int, int, int);
void write_smatrix(smatrix*, FILE*, int);
smatrix* read_smatrix(int, int, FILE*, int);
void expand_smatrix(smatrix*, svector*);
void smxv(smatrix*, int, sdddouble*, int, sdddouble*);

/* Internal dmatrix functions */
void free_dmatrix(dmatrix*);
dmatrix* create_dmatrix(dmatrix*, int, int);
void write_dmatrix(dmatrix*, FILE*, int);
dmatrix* read_dmatrix(int, FILE*, int);
void expand_dmatrix(dmatrix*, sddfloat);

/* Internal SDD functions */
sdd* create_sdd(sdd*, int, int, int, int);
void expand_sdd(sdd*, sddfloat, svector*, svector*);
void sddxsv(sdddouble*, sdd*, svector*, int);

/* Internal sparse matrix (matrix) functions */
double fnormsq(matrix*);

/* Internal functions the combine objects */
void matrixxsv(sdddouble*, matrix*, svector*, int);
sdddouble_plus subproblem(matrix*, sdd*, svector*, svector*, int, sdddouble*);

/* Internal utility functions */
int comparentry(const entry *a, const entry *b) {
  if (a->j == b->j) return (a->i > b->i ? 1 : -1);
  else return (a->j > b->j ? 1 : -1);
} 

#ifdef QSORTOPT			/* use optimized GNU qsort */
#include "qsortopt.h"
#else                           /* use system-provided qsort */
int compar(const sdddouble_plus *a, const sdddouble_plus *b) {
  return(a->val > b->val ? -1 : 1);
}
#endif

/*------------------------------------------------------------------------
  Functions for S-Vectors 
  ------------------------------------------------------------------------*/

int numwords(int x) 
{
  /* Returns the number of words needed to store x bits. On error,
     returns 0. */

  int i;			
  if (x < 0) return (0);	/* error */
  i = x >> IDXSHIFT;		/* whole words */
  if (x & IDXMASK) i ++;	/* partial word */
  return (i);

} /* numwords */

/*------------------------------------------------------------------------*/

void free_svector(svector *x)
{
  /* Frees the memory used by x. Does nothing if x is NULL. */

  if (x == NULL) {
    fprintf(stderr, "Warning: Trying to free null svector.\n");
    return;
  }
  free(x->val);
  free(x->sgn);
  free(x);
  return;

} /* free_svector */

/*------------------------------------------------------------------------*/

svector* create_svector(svector *x, int m) 
{
  /* Creates an svector that holds m svalues. If x is non-null, the
     first min(x->m, m) entries are preserved. On error, returns
     NULL. */

  int nwrds = numwords(m);

  if (x == NULL) { /* create new svector */

    if ((x = (svector*) calloc(1, sizeof(svector))) == NULL) {
      fprintf(stderr, "Error allocating space for svector.\n");
      return (NULL);
    }
    x->m = m;
    if (((x->val = (ulong*) calloc(nwrds, sizeof(ulong))) == NULL) ||
	((x->sgn = (ulong*) calloc(nwrds, sizeof(ulong))) == NULL)) {
      fprintf(stderr, "Error allocating space for svector components.\n");
      return (NULL);
    }

  } /* if - create new svector */

  else { /* enlarge or shrink previous svector */

    x->m = m;			
    if (((x->val = (ulong*) realloc(x->val, nwrds * sizeof(ulong))) == NULL) ||
	((x->sgn = (ulong*) realloc(x->sgn, nwrds * sizeof(ulong))) == NULL)) {
      fprintf(stderr, "Error shrinking/enlarging svector.\n");
      return (NULL);
    }

  } /* else - enlarge or shrink previous svector */

  return(x);

} /* create_svector */

/*------------------------------------------------------------------------*/

void write_svector(svector *x, FILE *fptr, int bflag) 
{  
  /* Writes x to the file pointed to by fptr in text (bflag=0) or
     binary (bflag=1) format. On error, does nothing. */

  int i;			/* counter */
  ulong *valptr, *sgnptr;	/* local pointers */
  ulong mask;			/* bit mask */
  ulong val, sgn;		/* current val and sgn words */

  if (x == NULL) {
    fprintf(stderr, "Warning: Trying to write NULL svector to file.\n");
    return;
  }

  if (bflag) { /* binary */

    i = numwords(x->m);		
    if (((fwrite(x->val, sizeof(ulong), i, fptr)) < i) ||
	((fwrite(x->sgn, sizeof(ulong), i, fptr)) < i)) {
      fprintf(stderr, "Error writing svector to binary file.\n");
      return;
    }

  } /* if - binary */

  else { /* text */

    valptr = x->val;		
    sgnptr = x->sgn;		
    val = *valptr;		
    sgn = *sgnptr;		
    mask = ONE;		

    /* Loop through each s-value in x. */
    for (i = 0; i < x->m; i ++) {
      
      /* Print out the appropriate value. */
      if (val & mask) {	
	if (sgn & mask)	
	  fprintf(fptr, " -1");
	else
	  fprintf(fptr, "  1");
      }
      else
	fprintf(fptr, "  0");

      /* Update mask, val, and sgn. */
      if (mask == MAXMASK) {
	mask = ONE;	
	val = *(++valptr);	
	sgn = *(++sgnptr);	
	fprintf(fptr, "\n");
      }
      else
	mask <<= 1;		
      
    } /* i-loop */

    if (mask != (ONE)) fprintf(fptr, "\n");

  } /* else - text */

  return;

} /* write_svector */

/*------------------------------------------------------------------------*/

svector* read_svector(int m, FILE *fptr, int bflag) 
{  
  /* Reads an svector from the file pointed to by fptr in text
     (bflag=0) or binary (bflag=1) format. Returns a pointer to the
     svector. On error, returns NULL. */

  int i;			/* counter */
  svector *x;			/* new svector */
  int tmp;			/* svalue read from file */
  int idx;			/* index */
  ulong mask;			/* bit mask */

  i = numwords(m);		

  if ((x = create_svector(NULL, m)) == NULL) 
    return (NULL);
  
  if (bflag) { /* binary */

    if (((fread(x->val, sizeof(ulong), i, fptr)) < i) ||
	((fread(x->sgn, sizeof(ulong), i, fptr)) < i)) {
      fprintf(stderr, "Error reading svector from binary file.\n");
      return (NULL);
    }

  } /* if - binary */

  else { /* text */
    
    bzero(x->val, numwords(m) * sizeof(ulong));   /* set val to all zeros */
    bzero(x->sgn, numwords(m) * sizeof(ulong));   /* set sgn to all zeros */

    for (i = 0; i < m; i ++) {
      
      if ((fscanf(fptr, "%d", &tmp)) < 1) { 
	fprintf(stderr, "Error reading svector from text file.\n");
	return (NULL);
      }

      if (tmp != 0) {
	idx = i >> IDXSHIFT;	/* word index */
	mask = ONE << (i & IDXMASK);  /* bit mask */
	x->val[idx] |= mask;	/* turn correct val bit to 1 */
	if (tmp == -1)
	  x->sgn[idx] |= mask;	/* turn correct sgn bit to 1 */
      }

    } /* i-loop */

  } /* else - text */

  return (x);

} /* read_svector */

/*------------------------------------------------------------------------*/

sdddouble *init_threshold(svector *y, matrix *A, sdd *B, sddfloat rho, int *idx)
{
  int i;			/* counter */
  int m, n;			/* size of A */
  int localidx;			/* initialization index */
  sdddouble *s;			/* return vector */
  sdddouble sqnorms = 0;	/* squared norm of s */

  if ((y == NULL) || (A == NULL) || (B == NULL)) { /* error */
    fprintf(stderr, "Error in initialization of y.\n");
    return (NULL);
  }

  /* Allocate space for s and init to zero. */
  m = (A->m > A->n) ? A->m : A->n;
  if ((s = (sdddouble*) calloc(m, sizeof(sdddouble))) == NULL) {
    fprintf(stderr, "Error allocating space for s.\n");
    return (NULL);
  }

  /* Set sizes. */
  m = A->m;
  n = A->n;

  while (sqnorms < (rho / n)) {

    /* Set y */
    bzero(y->val, numwords(n) * sizeof(ulong)); 
    bzero(y->sgn, numwords(n) * sizeof(ulong)); 
    localidx = (*idx) % n;
    y->val[localidx >> IDXSHIFT] |= ONE << (localidx & IDXMASK);

    /* Update iidx */
    (*idx) ++;

    /* Compute s */
    matrixxsv(s, A, y, 0);	/* s = A * y */
    sddxsv(s, B, y, 0);		/* s = s - B * y */

    /* Compute squared norm of s */
    sqnorms = 0;
    for (i = 0; i < m; i ++)
      sqnorms += s[i] * s[i];
    
  } /* while loop */

  return (s);
}

/*------------------------------------------------------------------------*/

void init_cycle(svector *y, int idx)
{
  int n;			/* size of A */
  int localidx;			/* initialization index */

  if (y == NULL) {
    fprintf(stderr, "Error trying to initialize NULL y-vector.\n");
    return;	
  }

  n = y->m;

  /* Set y */
  bzero(y->val, numwords(n) * sizeof(ulong)); 
  bzero(y->sgn, numwords(n) * sizeof(ulong)); 
  localidx = idx % n;
  y->val[localidx >> IDXSHIFT] |= ONE << (localidx & IDXMASK);

  return;
}

/*------------------------------------------------------------------------*/

void init_ones(svector *x) 
{
  /* Initializes x so that every 100th s-value is 1, and all other
     values are set to 0. Does nothing on error. */

  int m;			/* size of x */

  if (x == NULL) {
    fprintf(stderr, "Error trying to initialize NULL y-vector.\n");
    return;	
  }

  m = x->m;

  bone(x->val, numwords(m) * sizeof(ulong));   /* set x to all ones */
  bzero(x->sgn, numwords(m) * sizeof(ulong));   /* set x to all oness */

  return;

} /* init_pones */

/*------------------------------------------------------------------------*/

void init_pones(svector *x) 
{
  /* Initializes x so that every 100th s-value is 1, and all other
     values are set to 0. Does nothing on error. */

  int i;			/* counter */
  int m;			/* size of x */

  if (x == NULL) {
    fprintf(stderr, "Error trying to initialize NULL y-vector.\n");
    return;	
  }

  m = x->m;

  bzero(x->val, numwords(m) * sizeof(ulong));   /* set x to all zeros */
  bzero(x->sgn, numwords(m) * sizeof(ulong));   /* set x to all zeros */

  /* set every 100th element in x to 1 */
  for (i = 0 ; i < m ; i += 100) {
    x->val[i >> IDXSHIFT] |= ONE << (i & IDXMASK); 
  }

  return;

} /* init_pones */

/*------------------------------------------------------------------------*/

int svcount(svector *x) 
{
  /* Returns the number of non-zero values in x. On error, returns
     zero. */

  int i, j;			/* counters */
  int cnt;			/* counter */
  int nwrds;			/* number of whole words used in x */
  ulong val;			/* temporary variable */
  ulong mask;			/* mask for final word */
  static int *cntlookup = NULL;	/* static look-up table */

  if (x == NULL) return(0);	/* error */

  /* The array cntlookup is a lookup table of length 2^8 that, for
     each byte used as an index, returns the number of 1-bits in that
     byte. Allocate space for and compute the static look-up table
     cntlookup.  Note that since this variable is static, we only
     compute the look-up table the first time this function is
     called. We use this table for speed. For a given byte x,
     cntlookup[x] = number of one bits in x.  */

  if (cntlookup == NULL) {
    if ((cntlookup = (int*) calloc(256, sizeof(int))) == NULL) {
      fprintf(stderr, "Error allocating space for look-up table.\n");
      return(0);
    }
    for (i = 0; i < 256; i ++) {
      mask = i;
      cnt = 0;
      for (j = 0; j < 8; j ++) {
	if (mask & 1) cnt ++;
	mask >>= 1;
      }
      cntlookup[i] = cnt;
    }
  }

  nwrds = numwords(x->m);	/* number of words in x */
  
  /* compute mask for final partial word with (x->m & idxmask) 1-bits */
  mask = (ONE << (x->m & IDXMASK)) - 1;
  if (mask) x->val[nwrds - 1] &= mask;

  /* Count the number of 1 bits in the whole words of val. For
     efficiency, we look at the val array one byte at a time (rather
     than one bit at a time) using the cntlookup array.  The means we
     trade 8 right shifts, 8 AND's, and 8 adds for one right shift,
     one and, one add, and one array look-up. Generally, the cntlookup
     array should easily fit into cache memory. */

  cnt = 0;
  for (i = 0 ; i < nwrds; i ++) {
    val = x->val[i];
    for (j = 0 ; j < sizeof(ulong) ; j ++) {
      cnt += cntlookup[val & 0xff];
      val >>= 8;
    }
  }

  return(cnt);

} /* svcount */

/*------------------------------------------------------------------------
  Functions for S-Matrices
  ------------------------------------------------------------------------*/

smatrix* create_smatrix(smatrix *A, int kmin, int kmax, int m) 
{
  /* Returns a ptr to an smatrix of maximum dimension m by kmax, and
     of current dimension m by kmin.  If A is non-null, then the first
     min(A->m, m) s-values of the first kmin columns of A are
     preserved. On error, returns NULL. */

  int k;			/* counter */
  
  if ((kmax <= 0) || (m <= 0) || (kmin < 0)) {	/* error checking */
    fprintf(stderr, "Error in requested size of smatrix.\n");
    return (NULL);
  }

  if (kmin > kmax) {		/* error checking */
    fprintf(stderr, "Max k-dimension less than min dimension.\n");
    return (NULL);
  }

  if (A == NULL) { /* starting from scratch */

    if (kmin != 0) {		/* error checking */
      fprintf(stderr, "Invalid kmin dimension.\n");
      return (NULL);
    }

    if ((A = (smatrix*) calloc(1, sizeof(smatrix))) == NULL) {
      fprintf(stderr, "Error allocating memory for smatrix.\n");
      return (NULL);
    }

    A->k = 0;			/* current columns in A */
    A->kmax = kmax;		/* max size columns in A */
    A->m = m;			/* number of rows in A */

    /* allocate column pointers */
    if ((A->col = (svector**) calloc(kmax, sizeof(svector*))) == NULL) {
      fprintf(stderr, "Error allocating space for smatrix columns.\n");
      return (NULL);
    }

  } /* if - starting from scratch */

  else { /* enlarging or shrinking A */

    /* If nothing needs to be changed, just return A. */
    if ((A->k == kmin) && (A->kmax == kmax) && (A->m == m)) 
      return (A);

   /* If the current A has too many columns defined, delete some. */
    if (kmin < A->k) {      
      for (k = kmin; k < A->k; k ++) 
	free_svector(A->col[k]);
      A->k = kmin;
    }

    /* If the current A does not have the right number of rows,
       re-create each column vector, adjusting to the new size. */
    if (m != A->m) {        
      for (k = 0; k < A->k; k ++) 
	if ((A->col[k] = create_svector(A->col[k], m)) == NULL) {
	  fprintf(stderr, "Error recreating svectors.\n");
	  return (NULL);
	}
      A->m = m;
    }

    /* If A does not have the right maximum number of columns, adjust
       that. */
    if (kmax != A->kmax) {  
      if ((A->col = (svector**) realloc(A->col, kmax * sizeof(svector*))) 
	  == NULL) {
	fprintf(stderr, "Error enlarging number of columns in A matrix,\n");
	return (NULL);
      }
      A->kmax = kmax;
    }

  } /* else - enlarging / shrinking A */

  return (A);

} /* create_smatrix */

/*------------------------------------------------------------------------*/

void write_smatrix(smatrix *A, FILE *fptr, int bflag) 
{
  /* Writes the smatrix A to the file pointed to by fptr in text
     (flag=0) or binary (bflag=1) format. */

  int k;			/* counter */

  /* Simply write each svector to the file in sequence. */
  for (k = 0; k < A->k; k ++)
    write_svector(A->col[k], fptr, bflag);

  return;

} /* write_smatrix */

/*------------------------------------------------------------------------*/

smatrix* read_smatrix(int m, int k, FILE *fptr, int bflag) 
{
  /* Read in an smatrix of size m by k from the file pointed to by
     fptr in text (bflag=0) or binary (bflag=1) format.  Return a
     pointer to the new smatrix. On error, returns NULL. */

  smatrix *A;			/* smatrix to be read in */

  /* Create the matrix */
  if ((A = (smatrix*) calloc(1, sizeof(smatrix))) == NULL) {
    fprintf(stderr, "Error allocating memory for smatrix.\n");
    return (NULL);
  }

  A->m = m;
  A->k = k;
  A->kmax = k;

  if ((A->col = (svector**) calloc(k, sizeof(svector*))) == NULL) {
    fprintf(stderr, "Error allocating memory for smatrix columns.\n");
    return (NULL);
  }

  /* Read in the columns in sequence */
  for (k = 0; k < A->k; k ++) 
    if ((A->col[k] = read_svector(m, fptr, bflag)) == NULL) {
      fprintf(stderr, "Error reading smatrix from file.\n");
      return (NULL);
    }

  return (A);

} /* read_smatrix */

/*------------------------------------------------------------------------*/

void expand_smatrix(smatrix *A, svector *x) 
{  
  /* Add x as the next column in A. On error, prints message. */

  if (x->m != A->m) {
    fprintf(stderr, "Error: Size mismatch.\n");
    return;
  }

  if (A->k == A->kmax) {
    fprintf(stderr, "Error: Smatrix is full.\n");
    return;
  }

  A->col[A->k] = x;
  A->k ++;
  return;

} /* expand_smatrix */

/*------------------------------------------------------------------------*/

void free_smatrix(smatrix *A) 
{
  /* Frees the memory used by A. */

  int k;			/* counter */

  if (A == NULL) {
    fprintf(stderr, "Warning: Trying to free NULL smatrix.\n");
    return;
  }

  if (A->col != NULL) {
    for (k = 0; k < A->k; k ++) 
      free_svector(A->col[k]);
    free(A->col);
  }
  free(A);
  return;

} /* free_smatrix */


/*------------------------------------------------------------------------
  Functions for D-Matrices
  ------------------------------------------------------------------------*/

dmatrix* create_dmatrix(dmatrix *D, int kmin, int kmax) 
{
  /* Returns a ptr to a dmatrix of maximum size kmax.  If D is
     non-null then the first kmin elements of D are preserved.  On
     error, returns NULL. */

  if ((kmin < 0) || (kmax <= 0) || (kmin > kmax)) { /* error checking */
    fprintf(stderr, "Error in parameters for dmatrix.\n");
    return (NULL);
  }

  if (D == NULL) {		/* start from scratch */

    if ((D = (dmatrix*) calloc(1, sizeof(dmatrix))) == NULL) {
      fprintf(stderr, "Error allocating space for dmatrix.\n");
      return (NULL);
    }

    D->k = 0;
    D->kmax = kmax;

    if ((D->d = (sddfloat*) calloc(kmax, sizeof(sddfloat))) == NULL) {
      fprintf(stderr, "Error allocating space for dmatrix values,\n");
      return (NULL);
    }

  } /* if - start from scratch */

  else {			/* enlarging or shrinking D */

    /* Delete extra elements, if necessary. */

    if (kmin > D->k) {
      fprintf(stderr, "Error: kmin too big for dmatrix.\n");
      return (NULL);
    }
    
    D->k = kmin;

    /* Expand or shrink maximum number of elements, if necessary. */

    if (D->kmax != kmax) {
      D->kmax = kmax;
      if ((D->d = (sddfloat*) realloc(D->d, kmax * sizeof(sddfloat))) == NULL) {
	fprintf(stderr, "Error in reallocating memory for dmatrix.\n");
	return (NULL);
      }
    }

  } /* else - enlarge/shrink */

  return (D);

} /* create_dmatrix */
 
/*------------------------------------------------------------------------*/

void write_dmatrix(dmatrix *D, FILE *fptr, int bflag) 
{
  /* Writes D to the file pointed to by fptr in text (bflag=0) or
     binary (bflag=1) format. On error, prints message.*/

  int k;			/* counter */

  if (bflag) { /* binary */

    if ((fwrite(D->d, sizeof(sddfloat), D->k, fptr)) < D->k) {
      fprintf(stderr, "Error writing dmatrix to binary file.\n");
      return;
    }

  } /* if - binary */

  else { /* text */
    
    for (k = 0; k < D->k; k ++) 
      fprintf(fptr, "%30.25e\n", D->d[k]); 

  } /* else - text */

  return;

} /* write_dmatrix */

/*------------------------------------------------------------------------*/

dmatrix* read_dmatrix(int k, FILE *fptr, int bflag) 
{
  /* Read a dmatrix of size k from the file pointed to by fptr in text
     format (bflag=0) or binary format (bflag=1). Return a ptr to the
     dmatrix. On error, returns NULL. */

  int i;			/* counter */
  dmatrix *D;			/* matrix to be read in */

  if ((D = (dmatrix*) calloc(1, sizeof(dmatrix))) == NULL) {
    fprintf(stderr, "Error allocating space for dmatrix.\n");
    return (NULL);
  }

  D->k = k;			/* set current size */
  D->kmax = k;			/* set max size */

  if ((D->d = (sddfloat*) calloc(k, sizeof(sddfloat))) == NULL) {
    fprintf(stderr, "Error allocating space for dmatrix entries.\n");
    return (NULL);
  }
  
  if (bflag) { /* binary */

    if ((fread(D->d, sizeof(sddfloat), k, fptr)) < k) {
      fprintf(stderr, "Error reading D from file.\n");
      return (NULL);
    }

  } /* if - binary */

  else { /* text */

    for (i = 0; i < k; i ++) 
      if ((fscanf(fptr, "%e", D->d + i)) < 1) {
	fprintf(stderr, "Error reading D from file.\n");
	return (NULL);
      }

  } /* else - text */

  return (D);

} /* read_dmatrix */

/*------------------------------------------------------------------------*/

void expand_dmatrix(dmatrix *D, sddfloat d) 
{
  /* Add d to the next open spot in D. On error, prints message. */

  if (D == NULL) {
    fprintf(stderr, "Cannot expand NULL Dmatrix.\n");
    return;
  }

  if (D->k == D->kmax) {
    fprintf(stderr, "Error: Dmatrix is full.\n");
    return;
  }

  D->d[D->k] = d;
  D->k ++;

  return;

} /* expand_dmatrix */

/*------------------------------------------------------------------------*/

void free_dmatrix(dmatrix *D)
{
  /* Free the memory used by D. Does nothing if D is null. */
  
  if (D == NULL) {
    fprintf(stderr, "Warning: Trying to free NULL D-Matrix.\n");
    return;
  }

  free(D->d);
  free(D);
  return;
 
} /* free_dmatrix */

/*------------------------------------------------------------------------
  Functions for SDD's
  ------------------------------------------------------------------------*/

sdd* create_sdd(sdd *A, int m, int n, int kmin, int kmax) 
{  
  /* Returns a pointer to an SDD for an m by n matrix. The maximum
     allowable number of terms is kmax. If A is non-null, then the
     first kmin terms are preserved. On error, returns NULL. */

  /* Error checking */
  if ((m <= 0) || (n <= 0) || (kmin < 0) || (kmax <= 0) || (kmin > kmax)) {
    fprintf(stderr, "Error in size parameters for SDD.\n");
    return (NULL);
  }

  if (A == NULL) { /* start from scratch */

    if (kmin != 0) {
      fprintf(stderr, "Warning: Specified nonzero existing kmax for new SDD.\n");
      kmin = 0;
    }

    /* Allocate space for A, and initialize pointers to NULL. */
    if ((A = (sdd*) calloc(1, sizeof(sdd))) == NULL) {
      fprintf(stderr, "Error allocating space for SDD.\n");
      return (NULL);
    }

  } /* if start from scratch */

  else { /* appending A */

    if (kmin > A->k) {
      fprintf(stderr, "Error in kmin parameter for SDD.\n");
      return (NULL);
    }

  } /* appending A */

  /* Fill in A */
  A->m = m;
  A->n = n;
  A->k = kmin;
  A->kmax = kmax;
  if (((A->X = create_smatrix(A->X, kmin, kmax, m)) == NULL) ||
      ((A->Y = create_smatrix(A->Y, kmin, kmax, n)) == NULL) ||
      ((A->D = create_dmatrix(A->D, kmin, kmax)) == NULL)) {
    fprintf(stderr, "Error allocating memory for SDD components.\n");
    return (NULL);
  }

  return (A);

} /* create_sdd */

/*------------------------------------------------------------------------*/

void write_sdd(sdd *A, char *fname, int bflag, char *note)
{
  /* Write A to the file named fname in text (bflag=0) or binary
     (bflag=1) format. The extra stuff (matrixf, edecompf, kmin, tol)
     is written in the comments for a text file. On error, print
     message. */

  FILE *fptr;

  if (A == NULL) {
    fprintf(stderr, "Error trying to write NULL SDD to file,\n");
    return;
  }

  /* Write Header Information */
  if (bflag) { /* binary */
  
    if ((fptr = fopen(fname, "wb")) == NULL) {
      fprintf(stderr, "Error opening output file.\n");
      return;
    }
    
    fwrite(&A->k, sizeof(int), 1, fptr);
    fwrite(&A->m, sizeof(int), 1, fptr);
    fwrite(&A->n, sizeof(int), 1, fptr);

  } /* if - binary */

  else { /* text */
    
    if ((fptr = fopen(fname, "w")) == NULL) {
      fprintf(stderr, "Error opening output file.\n");
      return;
    }

    fprintf(fptr, "%%%% Semidiscrete Decomposition (SDD)\n");
    if (note != NULL)
      fprintf(fptr, "%%%% %s\n", note);
    fprintf(fptr, "%d %d %d\n", A->k, A->m, A->n);

  } /* else - text */

  /* Write out components of SDD */
  write_dmatrix(A->D, fptr, bflag);
  write_smatrix(A->X, fptr, bflag);
  write_smatrix(A->Y, fptr, bflag);

  fclose(fptr);
  return;

} /* write_sdd */

/*------------------------------------------------------------------------*/

sdd* read_sdd(char *fname, int bflag) 
{  
  /* Read an sdd from the file named fname in text (bflag=0) or binary
     (bflag=1) format. On error, return NULL.  */

  sdd *A;			/* sdd to be read in */
  char s[MAXLINE];		/* line of input file */
  FILE *fptr;			/* input file pointer */

  if ((A = (sdd*) calloc(1, sizeof(sdd))) == NULL) {
    fprintf(stderr, "Error allocating memory for SDD.\n");
    return (NULL);
  }

  /* Get header information */
  if (bflag) { /* binary */

    if ((fptr = fopen(fname, "rb")) == NULL) {
      fprintf(stderr, "Error in opening input file.\n");
      free(A);
      return (NULL);
    }
    
    if (((fread(&A->k, sizeof(int), 1, fptr)) < 1) ||
	((fread(&A->m, sizeof(int), 1, fptr)) < 1) ||
	((fread(&A->n, sizeof(int), 1, fptr)) < 1)) {
      fprintf(stderr, "Error reading input file.\n");
      return (NULL);
    }

  } /* if - binary */

  else { /* text */

    if ((fptr = fopen(fname, "r")) == NULL) {
      fprintf(stderr, "Error in opening input file.\n");
      return (NULL);
    }

    /* Scan through comment lines */
    do 
      if (fgets(s, MAXLINE, fptr) == NULL) {
	fprintf(stderr, "Error reading header of text SDD input file\n");
	return (NULL);
      }
    while(s[0] == '%');

    if ((sscanf(s, "%d %d %d", &A->k, &A->m, &A->n)) < 3) {
      fprintf(stderr, "Error reading dimensions from SDD file.\n");
      return (NULL);
    }  

  } /* else - text */


  /* Fill in remainder of SDD */

  A->kmax = A->k;

  if (((A->D = read_dmatrix(A->k, fptr, bflag)) == NULL) ||
      ((A->X = read_smatrix(A->m, A->k, fptr, bflag)) == NULL) ||
      ((A->Y = read_smatrix(A->n, A->k, fptr, bflag)) == NULL)) {
    fprintf(stderr, "Error reading SDD.\n");
    return (NULL);
  }

  fclose(fptr);
  return (A);

} /* read_sdd */

/*------------------------------------------------------------------------*/

void expand_sdd(sdd *A, sddfloat d, svector *x, svector *y) 
{  
  /* Adds d and the svectors x and y to the next spot in A, and
     updates the number of terms in A. On error, prints warning. */

  if (A->k == A->kmax) {
    fprintf(stderr, "Error: SDD is full.\n");
    return;
  }
  
  expand_dmatrix(A->D, d);
  expand_smatrix(A->X, x);
  expand_smatrix(A->Y, y);
  A->k ++;
  return;

} /* expand_sdd */

void free_sdd(sdd *A){

  /* Frees the memory used by A. Does nothing if A is NULL. */
  
  if (A == NULL) {
    fprintf(stderr, "Warning: Trying to free NULL matrix.\n");
    return;
  }

  free_smatrix(A->X);
  free_smatrix(A->Y);
  free_dmatrix(A->D);
  free(A);
  return;

} /* free_sdd */

/*------------------------------------------------------------------------
  Subroutines which use a variety of data types.
  ------------------------------------------------------------------------*/

void sddxsv(sdddouble *s, sdd *A, svector *y, int tflag) 
{
  /* SDD times Svector. Computes s = s - A y where s is a A->m long
     real vector (of type sdddouble), A is an SDD, and y is an A->n long
     svector.  A is transposed if tflag=1, and then we assume s is
     A->n long and y is A->m long.  The vector s is modified.  All the
     other vectors are unchanged. On error, prints message. */

  int i, j, k;			/* counters */
  int kmax;			/* number of terms in SDD A */
  int nwrds;			/* number of words in y */
  int cnt;			/* counter */
  sdddouble tmp;		/* tmp value */
  smatrix *X, *Y;		/* local pointers */
  sddfloat *D;			/* local pointer */
  ulong val, sgn;		/* used in inner product calculation */
  ulong mask;			/* used for last partial word */
  ulong *xvalptr, *xsgnptr;	/* pointers to walk through arrays */
  ulong *yvalptr, *ysgnptr;	/* pointers to walk through arrays */
  sdddouble *tmpvec = NULL;	/* temporary workspace holds A y */
  static int *cntlookup = NULL;	/* static look-up table */

  if (A->k == 0) return;   /* A is empty, do nothing */

  /* Allocate space for tmpvec. */

  if ((tmpvec = (sdddouble*) calloc(A->kmax, sizeof(sdddouble))) == NULL) {
    fprintf(stderr, "Error allocating work space for sdd times vec.\n");
    return;
  }

  /* Allocate space for and compute the static look-up table
     cntlookup.  Note that since this variable is static, we only
     compute the look-up table the first time this function is
     called. We use this table for speed. For a given byte x,
     cntlookup[x] = number of one bits in x. */

  if (cntlookup == NULL) {
    if ((cntlookup = (int*) calloc(256, sizeof(int))) == NULL) {
      fprintf(stderr, "Error allocating space for look-up table.\n");
      return;
    }
    for (i = 0; i < 256; i ++) {
      mask = i;
      cnt = 0;
      for (j = 0; j < 8; j ++) {
	if (mask & 1) cnt ++;
	mask >>= 1;
      }
      cntlookup[i] = cnt;
    }
  }

  /* Assign X, Y, and D to local variables.  Note that transposing the
     SDD A only has the effect of swapping the X and Y matrices, so
     that is the only thing we need to do if tflag=1. */

  kmax = A->k;

  if (tflag) { /* Transpose => Swap X and Y */
    X = A->Y;
    Y = A->X;
  }
  else {
    X = A->X;
    Y = A->Y;
  }

  D = (A->D)->d;

  /* COMPUTE tmpvec = Y'y */

  if (y->m != Y->m) {		/* error checking */
    fprintf(stderr, "Size mismatch.\n");
    return;
  }
  
  yvalptr = y->val;
  ysgnptr = y->sgn;
  
  nwrds = numwords(y->m); 
  mask = (ONE << (Y->m & IDXMASK)) - 1;
  if (mask) ysgnptr[nwrds - 1] &= mask;

  /* Loop through each column of Y, computing its inner product with y
     which will be added to s[k]. */
  
  for (k = 0; k < kmax; k ++) {
    
    /* Set xvalptr and xsgnptr to the beginnings of the (k+1)st column
       of Y's val and sgn arrays. */
    
    xvalptr = (Y->col[k])->val;
    xsgnptr = (Y->col[k])->sgn;

    /* Mask final partial word, if it exists */
    if (mask) {
      xvalptr[nwrds-1] &= mask;
      xsgnptr[nwrds-1] &= mask;
    }

    cnt = 0;     /* This will be the value of the innner product when
		    we are done. */

    /* Do whole words first. */

    for (i = 0; i < nwrds; i ++) {
      
      /* Rather than going through comparing y and the (k+1)st column
         of Y s-value by s-value, we will handle BITS_PER_WORD
         s-values at once.  Assuming there are 32 bits_per_word, we
         are trading a _MINIMUM_ of 32 shifts, 32 AND's, and 32
         compares for 2 AND's, 1 XOR, 8 adds, 8 array look-ups, 4
         multiplies, and 8 right shifts. Generally, the cntlookup
         array should easily fit into cache memory.
        
         Here is a brief mathematical explanation of what we are
         doing: suppose we have two s-values a and b.  If we AND the
         val bits, we get 1 if a*b = +1 or -1, and 0 if a*b = 0; call
         this bit c.  If a*b = -1, the XOR of the sgn bits is 1, and
         if a*b = 1, the XOR of the sgn bits is 0.  If a*b = 0, the
         XOR of the sgn bits could be anything.  In order to make sure
         it is zero if a*b = 0, we AND c with the result of the XOR of
         the sgn bits; Call this d.  Finally, a * b = c - 2 * d.  We
         are doing the same thing below, but rather than dealing with
         one pair of s-values at a time, we are dealing with 32! */

      val = xvalptr[i] & yvalptr[i];
      sgn = xsgnptr[i] ^ ysgnptr[i];
      sgn &= val;
      
      for (j = 0; j < sizeof(ulong); j ++) {
	cnt += cntlookup[val & 0xff];
	cnt -= 2 * cntlookup[sgn & 0xff];
	val >>= 8;
	sgn >>= 8;
      }
    }

    tmpvec[k] = cnt;
  
  } /* k-loop */


  /*------------------------------
    Compute tmpvec = D * tmpvec
    ------------------------------*/

  for (k = 0; k < kmax; k ++)
    tmpvec[k] *= D[k];

  /*------------------------------
    Compute s -= X * tmpvec
    ------------------------------*/

  /* Cycle through the columns of X. */

  for (k = 0; k < kmax; k ++) {
    
    tmp = tmpvec[k];		/* Copy tmpvec[k] to a local variable */
    xvalptr = (X->col[k])->val;	/* Set xvalptr to the beginning to val */
    xsgnptr = (X->col[k])->sgn; /* Set xsgnptr to the beginning to sgn */ 
    val = *xvalptr;		/* Dereference the xvalptr */
    sgn = *xsgnptr;		/* Dereference the xsgnptr */
    mask = ONE;			/* Initialize the mask */


    /* This loop will go through each s-value in X->m. */

    for (i = 0; i < X->m; i ++) {
      
      /* Mask picks off the current s-value. */

      if (val & mask) {		/* i-th s-value is nonzero */
	if (sgn & mask)		/* X[i,k] = -1 */
	  s[i] += tmp;
	else			/* X[i,k] = 1 */
	  s[i] -= tmp;
      }

      /* Update mask and, if necessary, val and sgn. */

      if (mask == MAXMASK) {
	mask = ONE;		/* Reinitialize the mask. */
	val = *(++xvalptr);	/* Increment the pointer and dereference. */
	sgn = *(++xsgnptr);	/* Increment the pointer and dereference. */
      }
      else
	mask <<= 1;		/* Update the mask. */
      
    } /* i-loop */

  } /* k-loop */

  free(tmpvec);
  return;

} /* smxsv */

/*------------------------------------------------------------------------*/

void smxv(smatrix *X, int kmax, sdddouble *y, int tflag, sdddouble *r) 
{  
  /* Smatrix times vector. Computes r = X * y using only the
     first kmax columns of X where X is an smatrix, y is an kmax-long
     real vector (of type sdddouble), and r is an X->m long real vector
     (of type double-type).  If tflag is TRUE, then the transpose of X
     is used, y should be of length X->m and s should be of length
     kmax. On error, prints message. */

  int j, k, m;			/* counters */
  ulong *valptr, *sgnptr;	/* tmp pointers */
  ulong mask;			/* mask out certain bits */
  ulong val, sgn;		/* tmp values */
  sdddouble tmp;		/* tmp value */

  m = X->m;			/* copy to a local variable */

  if (tflag) { /* tranpose */

    /* Cycle through each column of X. */

    for (k = 0; k < kmax; k ++) {

      tmp = 0;			/* This will be r[k]. */
      valptr = (X->col[k])->val; /* Set to beginning of the val array
				  * for the (k+1)st column of X. */
      sgnptr = (X->col[k])->sgn; /* Set to the beginning of the sgn
				  * array for the (k+1)st column of X. */
      val = *valptr;		/* Dereference into local variable. */
      sgn = *sgnptr;		/* Dereference into local variable. */
      mask = ONE;		/* Initialize the mask. */

      /* Cycle through each value in y. */

      for (j = 0; j < m; j ++) {

	/* Compute y(j) * X(j,k) and update tmp. Mask picks off the
           correct s-value. */

	if (val & mask)	{	/* X(j,k) is nonzero */
	  if (sgn & mask)	/* X(j,k) is -1 */
	    tmp -= y[j];
	  else			/* X(j,k) is 1 */
	    tmp += y[j];
	}

	/* Update the mask and, if necessary, val and sgn. */

	if (mask == MAXMASK) {
	  mask = ONE;	/* Reinitialize the mask. */
	  val = *(++valptr);	/* Increment the pointer and dereference. */
	  sgn = *(++sgnptr);	/* Increment the pointer and dereference. */
	}
	else
	  mask <<= 1;		/* Update the mask. */

      } /* j-loop */

      r[k] = tmp;		/* Copy tmp to r[k]. */

    } /* k-loop */

  } /* if - transpose */
	  
  else { /* transpose */

    bzero(r, m * sizeof(sdddouble)); /* Set all elements of r to zero. */
    
    /* Cycle through each value in y. */

    for (k = 0; k < kmax; k ++) {

      tmp = y[k];		/* Copy to a local variable. */
      valptr = (X->col[k])->val; /* Set to beginning of the val array
				    for the (k+1)st column of X. */
      sgnptr = (X->col[k])->sgn; /* Set to the beginning of the sgn
				    array for the (k+1)st column of X. */ 
      val = *valptr;		/* Dereference into local variable. */
      sgn = *sgnptr;		/* Dereference into local variable. */
      mask = ONE;		/* Initialize the mask. */

      for (j = 0; j < m; j ++) {
	
	if (val & mask) {	/* X(j,k) is nonzero */
	  if (sgn & mask)	/* X(j,k) is -1 */
	    r[j] -= tmp;
	  else			/* X(j,k) is 1 */
	    r[j] += tmp;
	}

	/* Update the mask and, if necessary, val and sgn. */

	if (mask == MAXMASK) {
	  mask = ONE;		/* Reinitialize the mask. */
	  val = *(++valptr);	/* Increment the pointer and dereference. */
	  sgn = *(++sgnptr);	/* Increment the pointer and dereference. */
	}
	else
	  mask <<=1;		/* Update the mask. */

      } /* j-loop */

    } /* k-loop */
  
  } /* else - no transpose*/

  return;

} /* smxv */

/*------------------------------------------------------------------------*/

void free_matrix(matrix *A)
{

  /* Free memory used by A. Does nothing if A is NULL. */

  if (A == NULL) {
    fprintf(stderr, "Warning: Trying to free NULL matrix.\n");
    return;
  }
  free(A->jc);
  free(A->ir);
  free(A->val);
  free(A);
  return;

} /* free_matrix */

void write_matrix(matrix *A, char *fname, int bflag)
{
  /* Writes A to file named fname in text (bflag=0) or binary
     (bflag=1) format. Text output is compatible with MatrixMarket
     format. Does nothing if A is null. On error, prints message. */

  int i, j;			/* counters */
  FILE *out;			/* pointer to output file */

  if (A == NULL) return;

  if (bflag) { /* binary */

    if ((out = fopen(fname, "wb")) == NULL) {
      fprintf(stderr, "Error opening binary sparse matrix output file.\n");
      return;
    }
    
    fwrite(&A->m, sizeof(int), 1, out);
    fwrite(&A->n, sizeof(int), 1, out);
    fwrite(&A->nnz, sizeof(int), 1, out);
    fwrite(A->jc, sizeof(int),  A->n + 1, out);
    fwrite(A->ir, sizeof(int), A->nnz, out);
    fwrite(A->val, sizeof(float), A->nnz, out);

    fclose(out);

  } /* binary */

  else { /* text */

    if ((out = fopen(fname, "w")) == NULL) {
      fprintf(stderr, "Error opening text sparse matrix output file.\n");
      return;
    }

    fprintf(out, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(out, "%% File generated by SDDPACK.\n");
    fprintf(out, "%d ", A->m);
    fprintf(out, "%d ", A->n);
    fprintf(out, "%d\n", A->nnz);

    for (j = 0; j < A->n; j++)
      for (i = A->jc[j]; i < A->jc[j+1]; i ++) 
	fprintf(out, "%4d %4d %15.6e\n", A->ir[i] + 1, j + 1, A->val[i]);
    
    fclose(out);
    
  } /* text */

  return;

} /* write_matrix */

/*------------------------------------------------------------------------*/

matrix *read_matrix(char *fname, int bflag)
{

  /* Reads a sparse matrix in text (bflag=0) or binary (bflag=1)
     format from file named fname. (The text matrix should be in
     general real coordinate MatrixMarket format. If you prefer to use
     binary format, the convertmtx utility will convert a general real
     coordinate MatrixMarket file to binary format.) Creates and
     allocates memory for a new matrix structure to contain the sparse
     matrix. Returns a pointer to the matrix. On error, returns
     NULL. */

  int i, j, curj, cnt;		/* counters */
  float val;			/* matrix entry */
  char s[MAXLINE];		/* line of input file */
  FILE *in;			/* pointer to input file */
  entry *a;			/* entries of the sparse matrix */
  matrix *A;			/* pointer to sparse matrix */

  if (bflag) { /* binary */

    if ((in = fopen(fname, "rb")) == NULL) {
      fprintf(stderr, "Error opening binary sparse matrix input file.\n");
      return (NULL);
    }
    
    if ((A = (matrix *) malloc(sizeof(matrix))) == NULL) {
      fprintf(stderr, "Memory allocation error.\n");
      return (NULL);
    }
    
    if (((fread(&A->m, sizeof(int), 1, in)) < 1) ||
	((fread(&A->n, sizeof(int), 1, in)) < 1) ||
	((fread(&A->nnz, sizeof(int), 1, in)) < 1)) {
      fprintf(stderr, "Error reading sparse matrix input file.\n");
      return (NULL);
    }
    
    if (((A->jc = (int *) calloc(A->n + 1, sizeof(int))) == NULL) ||
	((A->ir = (int *) calloc(A->nnz, sizeof(int))) == NULL) ||
	((A->val = (float *) calloc(A->nnz, sizeof(float))) == NULL)) {
      fprintf(stderr, "Memory allocation error.\n");
      return (NULL);
    }
    
    if (((fread(A->jc, sizeof(int), A->n + 1, in)) < A->n + 1) ||
	((fread(A->ir, sizeof(int), A->nnz, in)) < A->nnz) ||
	((fread(A->val, sizeof(float), A->nnz, in)) < A->nnz)) {
      fprintf(stderr, "Error reading sparse matrix input file.\n");
      return (NULL);
    }    
    
    fclose(in);

    return (A);

  } /* if - binary */

  else { /* text */

    if ((in = fopen(fname, "r")) == NULL) {
      fprintf(stderr, "Error opening text sparse matrix input file.\n");
      return (NULL);
    }
  
    /* Scan past comment lines. */
    do 
      if (fgets(s, MAXLINE, in) == NULL) {
	fprintf(stderr, "Error reading header of sparse matrix input file\n");
	return (NULL);
      }
    while(s[0] == '%');
    
    if ((A = (matrix*) malloc(sizeof(matrix))) == NULL) {
      fprintf(stderr, "Memory allocation error.\n");
      return (NULL);
    }

    if (sscanf(s, "%d %d %d", &A->m, &A->n, &A->nnz) != 3) {
      fprintf(stderr, "Error reading sparse matrix input file.\n");
      return (NULL);
    }

    if ((a = (entry*) calloc(A->nnz, sizeof(entry))) == NULL) {
      fprintf(stderr, "Memory allocation error.\n");
      return (NULL);
    }
    
    cnt = 0;
    while ((fscanf(in, "%d %d %f", &i, &j, &val)) == 3) {

      if ((i <= 0) || (i > A->m)) {
	fprintf(stderr, "Invalid row index i=%d, not in range 1 to %d\n", i, A->m);
	exit(1);
      }
      if ((j <= 0) || (j > A->n)) {
	fprintf(stderr, "Invalid row index j=%d, not in range 1 to %d\n", j, A->n);
	exit(1);
      }

      a[cnt].i = i - 1;
      a[cnt].j = j - 1;
      a[cnt ++].val = val;
    }
    
    fclose(in);

    if (cnt < A->nnz) {
      fprintf(stderr, "WARNING: Reported number of nonzeros (%d) ", A->nnz);
      fprintf(stderr, "is less than number read in (%d). Adjusting. \n", cnt);
      A->nnz = cnt;
    }
    
    qsort(a, cnt, sizeof(entry), (int(*)(const void*, const void*))comparentry);
    
    if (((A->jc = (int *) calloc(A->n + 1, sizeof(int))) == NULL) ||
	((A->ir = (int *) calloc(A->nnz, sizeof(int))) == NULL) ||
	((A->val = (float *) calloc(A->nnz, sizeof(float))) == NULL)) 
      {
	fprintf(stderr, "Memory allocating error.\n");
	return (NULL);
      }
    
    cnt = 0;
    curj = 0;
    A->jc[0] = 0;
    for (i = 0; i < A->nnz; i ++) {
      while (curj < a[i].j) A->jc[++ curj] = cnt;
      A->ir[cnt] = a[i].i;
      A->val[cnt ++] = a[i].val;
    }
    
    while (curj < A->n) A->jc[++ curj] = cnt;
    
    free(a);
    
    return (A);
    
  } /* else - text */

} /* read_matrix */

/*------------------------------------------------------------------------*/

void matrixxv(float *y, matrix *A, float *x, int tflag) 
{
  /* Computes y = A * x (tflag = 0) or y = A' * x (tflag = 1). Vectors
     x and y are assumed to be the appropriate sizes. The vector y is
     initialized to zero. On error, prints message. */

  int i, j;			/* counters */

  if ((y == NULL) || (A == NULL) || (x == NULL)) {
    fprintf(stderr, "Null input argument for matrix-vector multiply.\n");
    return;
  }

  if (tflag) { /* tranpose */
    for (j = 0; j < A->n; j ++) y[j] = 0;
    for (j = 0; j < A->n; j ++) 
      for (i = A->jc[j]; i < A->jc[j+1]; i ++)
	y[j] += A->val[i] * x[A->ir[i]];
  }
  else { /* no transpose */
    for (i = 0; i < A->m; i ++) y[i] = 0;
    for (j = 0; j < A->n; j ++) 
      for (i = A->jc[j]; i < A->jc[j+1]; i ++)
	y[A->ir[i]] += A->val[i] * x[j];
  }

  return;

} /* matrixxv */

/*------------------------------------------------------------------------*/

double fnormsq(matrix *A) 
{
  /* Returns the sum of all the squares of the elements of
     sparse matrix A. On error, returns zero. */

  int i;			/* counter */
  double f = 0;			/* sum */

  if (A == NULL) return (f);

  for (i = 0; i < A->nnz; i ++)
    f += A->val[i] * A->val[i];

  return (f);

} /* fnormsq */

/*------------------------------------------------------------------------*/

void matrixxsv(sdddouble *s, matrix *A, svector *y, int tflag) 
{
  /* Sparse matrix times svector. Computes s = A * y where A
     is transposed if tflag is TRUE. On error, prints message. */

  int i, j;			/* couters */
  int m, n;			/* size of A */
  int *ir, *jc;			/* tmp pointers */
  sddfloat *val;		/* tmp pointer */
  int row;			/* row number */
  int idx;			/* word index of y to get svalue */
  ulong mask;			/* mask to extract svalue */
  ulong *valptr, *sgnptr;	/* tmp pointers */

  m = A->m; 
  n = A->n; 
  ir = A->ir; 
  jc = A->jc; 
  val = A->val;
  valptr = y->val; 
  sgnptr = y->sgn;

  if (tflag) { /* tranpose A */
    
    bzero(s, n * sizeof(sdddouble)); /* set s to zero */

    for (j = 0; j < n; j ++) 

      for (i = jc[j]; i < jc[j+1]; i ++) {
	
	row = ir[i];
	idx = row >> IDXSHIFT;        /* element of y to look in */
	mask = ONE << (row & IDXMASK);  /* bit in y[idx] to pick off */
	
	if (valptr[idx] & mask) { /* nonzero */
	  if (sgnptr[idx] & mask) /* y[row] = -1 */
	    s[j] -= val[i];
	  else /* y[row] = 1 */
	    s[j] += val[i];
	}

      } /* i-loop */

  } /* if - transpose */

  else { /* no transpose */
    
    bzero(s, m * sizeof(sdddouble)); /* set s to zero */

    idx = 0;
    mask = 1;
    for (j = 0; j < n ; j ++) {

      if  (valptr[idx] & mask) { /* nonzero */
	if (sgnptr[idx] & mask) /* y[j] = -1 */
	  for (i = jc[j]; i < jc[j+1]; i ++) {
	    s[ir[i]] -= val[i];
	  }
	else /* y[j] = 1 */
	  for (i = jc[j]; i < jc[j+1]; i ++) {
	    s[ir[i]] += val[i];
	  }
      }

      if (mask == MAXMASK) {
	mask = 1;
	idx ++;
      }
      else
	mask <<= 1;

    } /* j-loop */

  } /* else - no transpose*/

} /* matrixxsv */

/*------------------------------------------------------------------------*/

sdddouble_plus subproblem(matrix *A, sdd *B, svector *y, svector *x, 
		       int tflag, sdddouble *s) 
{
  /* Solve the SDD subproblem.
    
     Returns the maximum value of 
    
     max (x' * s) / |x|_2^2 
     
     in 'val' and the number of nonzeros in x in 'idx'.  Here x is
     constrained to be an (A->m)-long svector and s equals (A - B)y.
     A is a sparse matrix, B is an SDD of the same size, and y is an
     (A->n)-long svector.  If tflag is TRUE, we use the transpose of A
     and B, x should be (A->n)-long, and y should be (A->m)-long.
     Only x is modified. On error, returns f with 0,0. */

  sdddouble_plus f;		/* final answer */
  int i, j, m;			/* counters */
  ulong *valptr, *sgnptr;	/* tmp pointer */
  ulong mask;			/* tmp mask */
  ulong onesmask;		/* word of all ones */
  sdddouble_plus *sorts = NULL;	/* workspace */
  int sflag;			/* compute s? */
  
  f.val = 0; f.idx = 0;		/* init f for error returns */
  
  if ((A == NULL) || (y == NULL) || (x == NULL)) {
    fprintf(stderr, "Error in subproblem input arguments.\n");
    return (f);
  }

  bone(&onesmask, sizeof(ulong)); /* word of all ones */

  /* allocate memory for the s and sorts */

  m = (A->m > A->n) ? A->m : A->n;

  if (s == NULL) {		/* s is not given */
    sflag = 1;
    if ((s = (sdddouble*) calloc(m, sizeof(sdddouble))) == NULL) {
      fprintf(stderr, "Error allocating work space.\n");
      return (f);
    }
  }
  else
    sflag = 0;			/* don't compute s */

  if ((sorts = (sdddouble_plus*) calloc(m, sizeof(sdddouble_plus))) == NULL) {
    fprintf(stderr, "Error allocating work space.\n");
    return (f);
  }
  
  m = x->m;			/* copy to a local variable */

  if (sflag) {
    matrixxsv(s, A, y, tflag);	/* s = A * y */
    sddxsv(s, B, y, tflag);	/* s = s - B * y */
  }

  /* Initialize of s-values in x to 1 */

  bone(x->val, numwords(m) * sizeof(ulong)); /* set every bit to 1 */
  bzero(x->sgn, numwords(m) * sizeof(ulong)); /* set every bit to 0 */

  /* Fill in the sorts array with abs(s), and change the i-th s-value
     in x to -1 if s[i] is negative. */

  valptr = x->val;		/* copy to local pointer */
  sgnptr = x->sgn;		/* copy to a local pointer */
  mask = ONE;			/* initialize */

  for (i = 0; i < m; i ++) {

    sorts[i].idx = i;		/* copy the index into array to be sorted */

    /* If s[i] is negative, set the i-th bit of x->sgn to 1 w/o
       affecting any other bits. Also, fill in sorts[i].val with
       abs(s[i]). */

    if (s[i] < 0) {		/* make x[i] = -1 */
      sorts[i].val = -s[i];	/* copy -s[i] into array to be sorted */
      *sgnptr |= mask;		/* swaps the correct bit of x->sgn to 1 */
    }
    else
      sorts[i].val = s[i];	/* copy s[i] into array to be sorted */

    /* Update the mask and, if necessary, the pointer. */

    if (mask == MAXMASK) {
      mask = ONE;
      sgnptr ++;
    }
    else
      mask <<= 1;
  }      

  /* sort sorts */
#ifdef QSORTOPT
  qsortopt((char *) sorts, m);
#else
  qsort(sorts, m, sizeof(sdddouble_plus), (int(*)(const void*, const void*)) compar);
#endif  

  /* compute partial sums */
  s[0] = sorts[0].val;
  for (i = 1; i < m; i ++)
    s[i] = s[i-1] + sorts[i].val;

  /* compute function values */
  for (i = 0; i < m; i ++)
    s[i] = (s[i] * s[i]) / (i + 1);

  /* find the maximum of the array s */
  f.val = s[0];
  f.idx = 0;
  for (i = 1; i < m; i ++)
    if (s[i] > f.val) {
      f.val = s[i];
      f.idx = i;
    }

  /* increment f.idx so that it is the number of nonzeros in x */
  f.idx ++;

  /* zero out certain elements of x with indices > f.idx in the sorted
     array */

  for (i = f.idx; i < m; i ++) {
    
    /* (j+1)st s-value in x should be zero */
    j = sorts[i].idx;

    /* causes bit j in x->val to be zero w/o affecting other bits */
    valptr[j >> IDXSHIFT] &= ((ONE << (j & IDXMASK)) ^ onesmask);

  }

  free(s); 
  free(sorts);

  return (f);

} /* subproblem */

/*------------------------------------------------------------------------*/

sdd* compute_sdd (matrix *A, sdd *B, int kmin, int kmax, float rhomin,
		  int lmax, float alphamin, int initflag) 
{

  /* Computes the SDD of A. If the sdd B is non-null, preserves the
     first kmin triplets, and then continues to expand. The SDD of A
     will be such that either the norm of the residual is less than
     rhomin or the number of terms is kmax. The inner iterations are
     controlled by alphamin (the improvement tolerance) and lmax (the
     maximum number of inner iterations). The initflag controls the
     method used to initialize y; the choices are ... (TGK - fill this
     in). */

  int k;			/* current number of terms in of SDD */
  int m, n;			/* size of matrix */
  sddfloat d;			/* current d-value */
  svector *x, *y;		/* current x & y svectors */
  sdddouble_plus xmax, ymax;	/* solutions to subproblems */
  int l, totall;		/* iteration count and total */
  sddfloat rho0 = 0;		/* square of initial residual norm */
  sddfloat rho = 0;		/* square of residual norm  */
  sddfloat alpha;		/* improvement in change in residual */
  sddfloat beta, betabar;	/* change in residual, and previous */
  sdddouble *s = NULL;		/* product A*y produced by initialization */
  int initidx = 0;		/* index used in initialization */

  if (A == NULL) {
    fprintf(stderr, "Error trying to compute SDD of NULL matrix.\n");
    return (NULL);
  }

  if (kmin == kmax) {
    fprintf(stderr, "No expansion is necessary.\n");
    return (NULL);
  }

  /* initialization of lengths */
  m = A->m;
  n = A->n;

  /* (Re)initialization of sdd decomposition (stored in B). */
  if ((B = create_sdd(B, m, n, kmin, kmax)) == NULL) {
    fprintf(stderr, "Error initializing SDD.\n");
    return (NULL);
  }

#if INFO > 1
  fprintf(stdout, "\n");
#endif

  /* Compute initial residual if kmin is zero. */
  rho = rho0 = fnormsq(A);
  
  /* Compute actual rho if expanding existing decomp. */
  if (kmin != 0) {
    fprintf(stdout, "rho0 = %e\n", rho);
    for (k = 0; k < kmin; k ++) {
      rho -= (B->D->d[k] * B->D->d[k]) * svcount(B->X->col[k]) * svcount(B->Y->col[k]);
      fprintf(stdout, "rho[%d] = %e\n", k+1, rho);
    }
  }

  totall = 0;

#if INFO > 1
  fprintf(stdout, "Iteration Residual Improvement Inner  Total\n");
  fprintf(stdout, "  Number  Squared    (beta)     Its  InnerIts\n");
  fprintf(stdout, "--------- -------- ----------- ----- --------\n");
#endif
  
  /* Outer iterations */
  for (k = kmin; k < kmax; k ++) {

#if INFO > 1
    fprintf(stdout, "  %3d    ", k+1);
#endif

    if (((x = create_svector(NULL, m)) == NULL) ||
	((y = create_svector(NULL, n)) == NULL)) {
      fprintf(stderr, "Error creating svectors.\n");
      return (NULL);
    }

    switch(initflag) {
    case 1: /* Threshold */
      s = init_threshold(y, A, B, rho, &initidx);
#if INFO > 1
    fprintf(stdout, "   %3d    ", initidx);
#endif
      break;
    case 2: /* Cycling */
      init_cycle(y, k);
      break;
    case 3: /* All Ones */
      init_ones(y);
      break;
    case 4: /* Periodic Ones */
      init_pones(y);
      break;
    default:
      fprintf(stderr, "Error in initflag.\n");
      return (NULL);
    }    

    beta = betabar = 0;

    for (l = 0; l < lmax; l ++) { /* inner iteration */
      
      xmax = subproblem(A, B, y, x, 0, s);
      s = NULL;
      ymax = subproblem(A, B, x, y, 1, NULL);

      beta = ymax.val / xmax.idx;

      if (l > 0) {
	alpha = (beta - betabar) / betabar;
	if (alpha < alphamin) {
	  l ++;
	  break;
	}
      }
      betabar = beta;

    } /* l-loop */

    d = (sqrt(ymax.val * ymax.idx)) / (xmax.idx * ymax.idx);
    expand_sdd(B, d, x, y);
    totall += l;
    rho = rho - beta;

#if INFO > 1
    fprintf(stdout, "%7.2e %10.5e  %2d    %4d\n", rho, beta, l, totall);
#endif

    if (rho <= rhomin) {
      if (rho < 0) rho = 0;
      k ++;
      break;
    }

  } /* end k-loop */


#if INFO > 0

  /* Output information about the decomposition. */
  fprintf(stdout, "\n");
  fprintf(stdout, "      -- SDD information --\n");
  fprintf(stdout, "final residual norm         : ");
  fprintf(stdout, "%10.4e\n", sqrt(rho));
  fprintf(stdout, "final relative residual norm: ");
  fprintf(stdout, "%5.3f\n", sqrt(rho / rho0));
  fprintf(stdout, "total outer iterations      : ");
  fprintf(stdout, "%d\n", (k - kmin));
  fprintf(stdout, "average inner iterations    : ");
  fprintf(stdout, "%5.3f\n", ((float) totall) / (k - kmin));
  if (initflag == 1) {
    fprintf(stdout, "average init iterations     : ");
    fprintf(stdout, "%5.3f\n", ((float) initidx) / (k - kmin));
  }
  fprintf(stdout, "\n");

#endif  

  return (B);

} /* compute_sdd */  

