// This code the work of Richard J. Hanson & peturbed by Krogh
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <float.h>
#include <complex.h>
#include <stdlib.h>
#include <omp.h>

#ifndef ENUM_
#define ENUM_ 1
enum MESSY_ENUMERATIONS
{ /*Used to flag groups of arguments in a variable argument list.*/
	m_rdat  =2013, /* Arbitrary, but avoids the value 0, ending the var arg list */
	m_rmat, m_zdat, m_zmat, m_idat, m_imat, m_ix, m_ptext
};
#endif
struct smessy_ty{
	/*Match 1-1 with derived type ?messy_ty
	found in ?messycall_m.F90. These values are not
	initialized in C so a call is made that gets a 
	copy with the Fortran default values.*/
	char ename[32];
	int fpprec;
	int kdf;
	int line_len;
	int munit;
	int eunit;
	int maxerr;
	int lstop;
	int lprint;
	int errcnt;
	int dblev;
	int cinit; /* Ttem used as a signature indicating a valid copy of cmessy_ty */
	int cthread; /*Thread index when using OpenMP; default is 0, no threading.*/
	int cstruct; /* Used in callmessy() to pick out the copy of messy_ty to use.*/
};

#define rk C_FLOAT
#define ck C_FLOAT_COMPLEX

/* Define the prototypes of the functions used for  interfacing to messy().*/
void get_smessy_defaults(struct smessy_ty_* e); /*Get default values for messy_ty into C struct cmessy_ty*/
void smessy(struct smessy_ty* e, char * const, ...); /* the ... designate var arg groups of [flags, sizes, pointers,] */

void callsmessy( /* The Fortran wrapper routine that calls messy() */
	struct smessy_ty* e, char s[], int slen,
	int, rk*,
	int, int, rk*,
	int, ck*,
	int, int, ck*,
	int, int*,
	int, int, int*,
	int, int*,
	char[],int);
/*
SUBROUTINE CALLSMESSY(e, STRING, slen, ND, D,&
MDA, NDA, DA,&
NZ, Z,&
MZA, NZA, ZA,&
NI, I,&
MIA, NIA, IA,&
NX, IX,&
PT, plen)&
BIND(C,name='callmessy') !The Fortran wrapper code ...
*/
void open_smessy_files( //.Open, close named files in Fortran for collecting output.
struct smessy_ty* e,
	int task, // Task:      1 for messages, 2 for errors, 3 for both
	char file_name[]); // Name of the file 

void close_smessy_files(
struct smessy_ty* e);// Fortran unit numbers are within structure e

void allocate_smessy_interface_(int maxcstructs, int maxcthreads);
//     The value of maxcstructs must be no smaller than the number
//     of different C structures smessy_ty  used in a C program.

//     The value of maxcthreads must be no smaller than the
//     number of parallel OpenMP threads used in a C program.

void deallocate_smessy_interface(void); // Resets maxcstructs, maxcthreads

