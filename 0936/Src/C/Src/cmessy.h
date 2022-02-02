// This code the work of Richard J. Hanson
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <float.h>
#include <complex.h>
#include <stdlib.h>
#include <omp.h>

enum MESSY_ENUMERATIONS
{ /*Used to flag groups of arguments in a variable argument list.*/
	m_rdat  =2013, /* Arbitrary, but avoids the value 0, ending the var arg list */
	m_rmat, m_zdat, m_zmat, m_idat, m_imat, m_ix, m_ptext
};

struct CMESSY_TY_{
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
	int cstruct; /* Index using internally in callmessy() to pick out the copy of messy_ty to use.*/
};

#define ck CTYPE_ _Complex

/* Define the prototypes of the functions used for  interfacing to messy().*/
void GET_CMESSY_DEFAULTS_(struct CMESSY_TY_* e); /*Get default values for messy_ty into C struct cmessy_ty*/
void MESSY_(struct CMESSY_TY_* e, char * const, ...); /* the ... designate var arg groups of [flags, sizes, pointers,] */

void CALLMESSY_( /* The Fortran wrapper routine that calls messy() */
	struct CMESSY_TY_* e, char s[], int slen,
	int, CTYPE_*,
	int, int, CTYPE_*,
	int, ck*,
	int, int, ck*,
	int, int*,
	int, int, int*,
	int, int*,
	char[],int);
/*
SUBROUTINE CALLMESSY(e, STRING, slen, ND, D,&
MDA, NDA, DA,&
NZ, Z,&
MZA, NZA, ZA,&
NI, I,&
MIA, NIA, IA,&
NX, IX,&
PT, plen)&
BIND(C,name='callmessy') !The Fortran wrapper code ...
*/
void OPEN_CMESSY_FILES_( //.Open, close named files in Fortran for collecting output.
struct CMESSY_TY_* e,
	int task, // Task:      1 for messages, 2 for errors, 3 for both
	char file_name[]); // Name of the file 

void CLOSE_CMESSY_FILES_(
struct CMESSY_TY_* e);// Fortran unit numbers are within structure e

void ALLOCATE_CMESSY_INTERFACE_(int maxcstructs, int maxcthreads);
//     The value of maxcstructs must be no smaller than the number
//     of different C structures cmessy_ty  used in a C program.

//     The value of maxcthreads must be no smaller than the
//     number of parallel OpenMP threads used in a C program.

void DEALLOCATE_CMESSY_INTERFACE_(void); // Needed for resetting maxcstructs, maxcthreads

