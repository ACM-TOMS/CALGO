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
  
  PROGRAM: convertmtx
  
  DESCRIPTION: Converts a text matrix (in MatrixMarket coordinate,
  real, general format) to SDDPACK binary format and vice versa.
  
  USAGE: convertmtx [options] infile outfile
  
  OPTIONS: -b Convert from binary to text. 
           -t Convert from text to binary. (DEFAULT)
  
  FILE FORMATS: Text files should be in MatrixMarket coordinate,
  real, general format. The header lines (beginning with '%') are NOT
  checked and need not be present. The output files will be in binary
  format.
  
  EXTERNAL SUBROUTINES: read_matrix, write_matrix (from sdd.c)
  
  PROGRAM DEPENDENCIES: convertmtx.c sdd.o sdd.h
  
  BUG REPORTS: Email Tamara.Kolda@na-net.ornl.gov.
  
  ----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sdd.h"		/* sparse matrix */

#define INDATA 't'		/* default input data type ('t' = text) */

void usage(char*);

int main (int argc, char *argv[]) 
{
  int i;			/* counter */
  char *infile;			/* input file name */
  char *outfile;		/* output file name */
  char indata = INDATA;		/* input file type: [t]ext or [b]inary */
  matrix *A = NULL;		/* pointer to sparse matrix */

  /* Parse command line */
  if (argc < 3) {
    usage(argv[0]);
    exit(-1);
  }
  infile = argv[argc - 2];
  outfile = argv[argc - 1];
  i = 1;
  while (i < argc - 2) {
    if (argv[i][0] != '-') {
      usage(argv[0]);
      exit(-1);
    }
    switch (argv[i][1]) {
    case 'b':
      indata = 'b';
      break;
    case 't':
      indata = 't';
      break;
    default:
      usage(argv[0]);
      exit(-1);
      break;
    } /* switch */
    i++;
  } /* while */

  /* Echo selections to standard output */
  fprintf(stdout, "*** output from %s ***\n", argv[0]);
  fprintf(stdout, "Input File  : %s\n", infile);
  fprintf(stdout, "Output File : %s\n", outfile);

  /* Call appropriate conversion routines */
  switch(indata) {
      
  case 'b': /* Convert binary matrix to text */
    
    if ((A = read_matrix(infile, 1)) == NULL) {
      fprintf(stderr, "Error reading matrix.\n");
      exit(-1);
    }
    write_matrix(A, outfile, 0);
    break;
    
    case 't': default: /* Convert text matrix to binary */
    
    if ((A = read_matrix(infile, 0)) == NULL) {
      fprintf(stderr, "Error reading matrix.\n");
      exit(-1);
    }
    write_matrix(A, outfile, 1);
    break;
    
  } /* switch */

  /* Free memory and exit */
  free_matrix(A);
  return(0);
  
} /* main */

/*----------------------------------------------------------------------*/

void usage(char *name) 
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s [options] infile outfile\n", name);
  fprintf(stderr, "Converts a matrix text file to binary or vice versa.\n");
  fprintf(stderr, "Options: -b Convert from binary to text\n");
  fprintf(stderr, "         -t Convert from text to binary (DEFAULT)\n");
  fprintf(stderr, "\n");
  return;

} /* usage */

