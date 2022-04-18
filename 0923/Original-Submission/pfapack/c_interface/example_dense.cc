/*****************************************************************
 Helper functions/macros and definitions to set up dense or banded
 matrices in Fortran format
******************************************************************/

#ifdef __cplusplus

inline int dense_fortran(int i, int j, int ldim)
{
  return i-1+(j-1)*ldim;
}

#include <complex>

typedef complex<double> doublecmplx
typedef complex<float> floatcmplx

#else

#define dense_fortran(i, j, ldim) (i-1+(j-1)*ldim)

typedef struct {
  double re, im;
} doublecmplx;

typedef struct {
  float re, im;
} floatcmplx;

#endif

void sskpfa(int, float *, int, float *, const char *, const char *, int *);
void dskpfa(int, double *, int, double *, const char *, const char *, int *);
void cskpfa(int, floatcmplx *, int, floatcmplx *,
	    const char *, const char *, int *);
void zskpfa(int, doublecmplx *, int, doublecmplx *,
	    const char *, const char *, int *);

void sskpfa10(int, float *, int, float *, const char *, const char *, int *);
void dskpfa10(int, double *, int, double *, const char *, const char *, int *);
void cskpfa10(int, floatcmplx *, int, floatcmplx *,
	      const char *, const char *, int *);
void zskpfa10(int, doublecmplx *, int, doublecmplx *,
	      const char *, const char *, int *);

void sskbpfa(int, int, float *, int, float *, const char *, int *);
void dskbpfa(int, int, double *, int, double *, const char *, int *);
void cskbpfa(int, int, floatcmplx *, int, floatcmplx *,
	     const char *, int *);
void zskbpfa(int, int, doublecmplx *, int, doublecmplx *,
	     const char *, int *);
