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
  
  PROGRAM: decomp
  
  DESCRIPTION: Computes (or expands) the Semidiscrete Decomposition
  (SDD) of a sparse matrix.
  
  USAGE: decomp [options] infile outfile
  
  OPTIONS: -k terms       Desired number of terms in SDD (Default: 100)
	   -a accuarcy    Desired SDD accuracy (Default: 0)
           -t tolerance   Tolerance for inner iterations (Default: 0.01)
	   -i its         Max inner iterations (Default: 100)
	   -y choice      Initialization for inner loop...
	                    1 = Threshold (Default)
			    2 = Cycling
			    3 = Ones
			    4 = Periodic Ones
           -e filename    Expand existing SDD  
           -b             I/O in binary format
	   
	   
  FILE FORMATS: 
  
  EXTERNAL SUBROUTINES:   
  read_sparse, free_sparse, read_sdd, compute_sdd, write_sdd, free_sdd - sdd.c
  
  PROGRAM DEPENDENCIES: decomp.c, sdd.h, sdd.o
  
  BUG REPORTS: Email Tamara.Kolda@na-net.ornl.gov.
  
  ----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sdd.h"		/* semidiscrete decomposition */

#define TERMS   100		/* default number of terms */
#define ACCR    0.0		/* default desired accuarcy */
#define TOL     0.01		/* default tolerance */
#define ITS     100		/* default max inner its */
#define YINIT   1		/* default initialization choice */
#define BFLAG   0		/* default I/O in binary format? */

void usage(char *);

int main(int argc, char* argv[])
{
  int i;			/* counter */
  int terms = TERMS;		/* number of terms in SDD */
  float accr = ACCR;		/* desired accuary of SDD */
  float tol = TOL;		/* tolerance */
  int its = ITS;		/* max its for inner iteration */
  int yinit = YINIT;		/* choice for initialization */
  int bflag = BFLAG;		/* I/O in binary format? */
  char *matrixf;		/* matrix file name */
  matrix *A;			/* matrix */
  char *decompf;		/* sdd file name */
  sdd *S;			/* sdd */
  char *edecompf = NULL;	/* existing sdd file name */
  int  eterms = 0;		/* no. of terms in existing sdd */
  sdd *E = NULL;		/* existing sdd  */
  char note[MAXLINE];		/* note for sdd file */

  /* Parse command line */
  if (argc < 3) {
    usage(argv[0]);
    return(-1);
  }
  matrixf = argv[argc-2];
  decompf = argv[argc-1];

  i = 1;
  while (i < argc - 2) {
    if (argv[i][0] != '-') {
      usage(argv[0]);
      return(-1);
    }
    switch (argv[i][1]) {
    case 'e':
      edecompf = argv[++i];
      break;
    case 'k':
      terms = atoi(argv[++i]);
      break;
    case 'a':
      accr = atof(argv[++i]);
      break;
    case 't':
      tol = atof(argv[++i]);
      break;
    case 'i':
      its = atoi(argv[++i]);
      break;
    case 'y':
      yinit = atoi(argv[++i]);
      break;
    case 'b':
      bflag = 1;
      break;
    default:
      usage(argv[0]);
      return(-1);
    }
    i++;
  }

  /* Output values being used */
  fprintf(stdout, "*** output from %s ***\n", argv[0]);
  fprintf(stdout, "matrix file        : %s\n", matrixf);
  fprintf(stdout, "SDD file           : %s\n", decompf);
  fprintf(stdout, "existing SDD file  : ");
  if (edecompf) fprintf(stdout, "%s\n", edecompf);
  else fprintf(stdout, "<none>\n");
  fprintf(stdout, "terms              : %d\n", terms);
  fprintf(stdout, "accuracy           : %f\n", accr);
  fprintf(stdout, "tolerance          : %f\n", tol);
  fprintf(stdout, "max inner its      : %d\n", its);
  fprintf(stdout, "y init choice      : %d\n", yinit);
  fprintf(stdout, "input type         : ");
  if (bflag) fprintf(stdout, "binary\n");
  else fprintf(stdout, "text\n");

  /* Create note for SDD file */
  if (!bflag) {
    if (edecompf)
      sprintf(note, "Matrix: %s Old SDD: %s Terms: %d Accr: %.2e Tol: %.2e InnIts: %d Init: %d",
	      matrixf, edecompf, terms, accr, tol, its, yinit);
    else
      sprintf(note, "Matrix: %s Terms: %d Accr: %.2e Tol: %.2e InnIts: %d Init: %d",
	      matrixf, terms, accr, tol, its, yinit);
  }

  /* Read matrix matrix */
  if ((A = read_matrix(matrixf, bflag)) == NULL) {
    fprintf(stderr, "Error reading sparse matrix.\n");
    return (-1);
  }

  /* Read exisiting SDD if one exists */
  if (edecompf) {
    if ((E = read_sdd(edecompf, bflag)) == NULL) {
      fprintf(stderr, "Error reading existing SDD from file.\n");
      return (-1);
    }
    eterms = E->k;
  }

  /* Compute SDD */
  S = compute_sdd(A, E, eterms, terms, accr, its, tol, yinit);

  /* Save SDD */
  write_sdd(S, decompf, bflag, note);

  /* Free memory and exit */
  free_matrix(A);
  free_sdd(S);
  return (0);

} /* main */

/*----------------------------------------------------------------------*/

void usage(char *name) {

  fprintf(stderr, "\n");
  fprintf(stderr, "USAGE: decomp [options] infile outfile\n");
  fprintf(stderr, "Computes Semidiscrete Decomposition (SDD).\n");
  fprintf(stderr, "OPTIONS: -k terms       Terms in SDD (Default: 100)\n");
  fprintf(stderr, "         -a accuarcy    Desired residual norm (Default: 0)\n");
  fprintf(stderr, "         -t tol         Tolerance (Default: 0.01)\n");
  fprintf(stderr, "         -i its         Max inner iterations (Default: 100)\n");
  fprintf(stderr, "         -y (1|2|3|4)   Initialization Choice\n");
  fprintf(stderr, "                          1 = Threshold (Default)\n");
  fprintf(stderr, "                          2 = Cycling\n");
  fprintf(stderr, "                          3 = All Ones\n");
  fprintf(stderr, "                          4 = Periodic Ones\n");
  fprintf(stderr, "         -e filename    Expand exisiting SDD\n");
  fprintf(stderr, "         -b             I/O in binary format\n");
  fprintf(stderr, "\n");
  return;

} /* usage */










