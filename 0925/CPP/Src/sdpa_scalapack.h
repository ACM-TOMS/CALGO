/* --------------------------------------------------
   prototype  collection of scalapack/BLACS routine
   from f2c.h, PBtools.h, PBblacs.h, Bdef.h PBpblas.h

   collected by Makoto Yamshita 2008.08.30
------------------------------------------------- */


/*--------------------------------------------------
  sdpa_scalapack.h
--------------------------------------------------*/
#ifndef __sdpa_scalapack_h__
#define __sdpa_scalapack_h__

extern "C" {

  // from PBtools.h  
#define DLEN1_ 9

  // from PBblacs.h
  void Cblacs_pinfo(int *, int *);
  void Cblacs_get(int, int, int *);
  void Cblacs_gridinit(int *, char *, int, int);
  void Cblacs_gridinfo(int, int *, int *, int *, int * );
  void Cblacs_barrier(int, char *);
  void Cblacs_gridexit(int );
  void Cblacs_abort(int, int);
  void Cblacs_exit(int );
  // Note : Cigebs2d, Cdgebs2d, Cigebr2d, Cdgebr2d, Cdgesum2d
  //        Cigesd2d, Cdgesd2d, Cigerv2d, Cdgerv2d
  // Address poineters are passed by (char*)
  void Cigebs2d(int, char *, char *, int, int, char *, int );
  void Cdgebs2d(int, char *, char *, int, int, char *, int );
  void Cigebr2d(int, char *, char *, int, int, char *,  int,
		int, int);
  void Cdgebr2d(int, char *, char *, int, int, char *, int,
		int, int);
  void Cdgsum2d(int, char *, char *, int, int, char*,
		int, int, int);
  void Cigesd2d(int, int, int, char*, int, int, int);
  void Cdgesd2d(int, int, int, char*, int, int, int);
  void Cigerv2d(int, int, int, char*, int, int, int);
  void Cdgerv2d(int, int, int, char*, int, int, int);
  // Bdef.h
  typedef char* F_CHAR_T;

  // from PBpblas.h
  void pdtradd_(F_CHAR_T, F_CHAR_T, int *, int *, double *,
		double *, int *, int *, int *, double *,
		double *, int *, int *, int * );
  void pdtrsv_(F_CHAR_T, F_CHAR_T, F_CHAR_T, int *,
	       double *, int *, int *, int *, double *, int *,
	       int *, int *, int *);
  void pdtrsm_(F_CHAR_T, F_CHAR_T, F_CHAR_T, F_CHAR_T,
	       int *, int *, double *, double *, int *,
	       int *, int *, double *, int *, int *, int *);
  void pdsyrk_(F_CHAR_T, F_CHAR_T, int *, int *,
	       double *, double *, int *, int *, int *,
	       double *, double *, int *, int *, int *);
  // #define PB_topget_ pb_topget__
  // #define PB_topset_ pb_topset__
#define PB_topget_ pb_topget_
#define PB_topset_ pb_topset_
  void PB_topget_(int *, F_CHAR_T, F_CHAR_T, F_CHAR_T);
  void PB_topset_(int *, F_CHAR_T, F_CHAR_T, F_CHAR_T);


  // from f2c.h
  typedef int integer;
  typedef unsigned long uinteger;
  typedef char *address;
  typedef short int shortint;
  typedef float real;
  typedef double doublereal;
  typedef struct { real r, i; } complex;
  typedef struct { doublereal r, i; } doublecomplex;
  typedef long int ftnlen;
  
  // from scalapack.h
  integer numroc_(integer *n, integer *nb, integer *iproc,
		  integer *isrcproc, integer *nprocs);
  int descinit_(integer *desc, integer *m, integer *n,
		integer *mb, integer *nb, integer *irsrc,
		integer *icsrc, integer *ictxt, integer *lld,
		integer *info);
  int pdlaprnt_(integer *m, integer *n, doublereal *a,
		integer *ia, integer *ja, integer *desca,
		integer *irprnt, integer *icprnt,
		char *cmatnm, integer *nout, doublereal *work,
		ftnlen cmatnm_len);
  int chk1mat_(integer *ma, integer *mapos0, integer *na,
	       integer *napos0, integer *ia, integer *ja,
	       integer *desca, integer *descapos0, integer *info);
  int pchk1mat_(integer *ma, integer *mapos0, integer *na,
		integer *napos0, integer *ia, integer *ja,
		integer *desca, integer *descapos0,
		integer *nextra, integer *ex, integer *expos,
		integer *info);
  int infog2l_(integer *grindx, integer *gcindx, integer *desc,
	       integer *nprow, integer *npcol, integer *myrow,
	       integer *mycol, integer *lrindx, integer *lcindx,
	       integer *rsrc, integer *csrc);
  integer iceil_(integer *inum, integer *idenom);

/*-------------------------------------------------
   scalapackmr.h

  This is header file to indicate
  Cp*gemr2d, Cp*trmr2d

  These functions distributs matrix from
  matirx with respect to one description array
  to another matrix with respect to
  different description.

   modified by Makoto Yamshita 2002.07.11 
-------------------------------------------------*/

  void Cpcgemr2d(int m, int n, complex* ptrmyblock,
		 int ia, int ja, int* desca,
		 complex* ptrmynewblock,
		 int  ib, int jb, int* descb,
		 int globcontext);
  void Cpdgemr2d(int m, int n, double* ptrmyblock,
		 int ia, int ja, int* desca,
		 double* ptrmynewblock,
		 int  ib, int jb, int* descb,
		 int globcontext);
  void Cpigemr2d(int m, int n, int* ptrmyblock,
		 int ia, int ja, int* desca,
		 int* ptrmynewblock,
		 int  ib, int jb, int* descb,
		 int globcontext);
  void Cpsgemr2d(int m, int n, float* ptrmyblock,
		 int ia, int ja, int* desca,
		 float* ptrmynewblock,
		 int ib, int jb, int* descb,
		 int globcontext);
  void Cpzgemr2d(int m, int n, doublecomplex* ptrmyblock,
		 int ia, int ja, int* desca,
		 doublecomplex* ptrmynewblock,
		 int ib, int jb, int* descb,
		 int globcontext);
  void Cpctrmr2d(char* uplo, char* diag,
		 int m, int n, complex* ptrmyblock,
		 int ia, int ja, int* desca,
		 complex* ptrmynewblock,
		 int  ib, int jb, int* descb,
		 int globcontext);
  void Cpdtrmr2d(char* uplo, char* diag,
		 int m, int n, double* ptrmyblock,
		 int ia, int ja, int* desca,
		 double* ptrmynewblock,
		 int  ib, int jb, int* descb,
		 int globcontext);
  void Cpitrmr2d(char* uplo, char* diag,
		 int m, int n, int* ptrmyblock,
		 int ia, int ja, int* desca,
		 int* ptrmynewblock,
		 int  ib, int jb, int* descb,
		 int globcontext);
  void Cpstrmr2d(char* uplo, char* diag,
		 int m, int n, float* ptrmyblock,
		 int ia, int ja, int* desca,
		 float* ptrmynewblock,
		 int ib, int jb, int* descb,
		 int globcontext);
  void Czstrmr2d(char* uplo, char* diag,
		 int m, int n, doublecomplex* ptrmyblock,
		 int ia, int ja, int* desca,
		 doublecomplex* ptrmynewblock,
		 int ib, int jb, int* descb,
		 int globcontext);


  // from Unknown
  // but I guess from scalapack

  void pxerbla_(int* ictxt, char* srname, int* info,
		ftnlen srname_len);
  
}; // end of extern "C"

#endif // __sdpa_scalapack_h__
