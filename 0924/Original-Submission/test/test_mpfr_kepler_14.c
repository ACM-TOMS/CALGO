/****************************************************************************
  This file is part of TIDES.
 
 Contributors:
 
 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"
#include "mpfr_kepler.h"

int main() {

	set_precision_digits(14);

	int i;
	int nvar = 4;
	int npar = 1;
	int nfun = 0;
	int nipt = 2;
	mpfr_t aux, error;
	mpfr_t v[nvar], x[nvar], p[npar], lt[nipt];
	mpfr_t tolrel, tolabs; 
	FILE   *fd;
	extern double   fac1,fac2,fac3;
	extern double   rmaxstep,rminstep; 
	extern int      nitermax, nordinc, minord, maxord; 
	extern int      defect_error_control;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	for(i=0; i<npar; i++) mpfrts_init(&p[i]);
	mpfrts_set_str(&p[0], "1."); 

/* --- INITIAL VALUES --- */
	for(i=0; i<nvar; i++) mpfrts_init(&v[i]);
	mpfrts_set_str(&v[0], "0.30000000000000000"); 
	mpfrts_set_str(&v[1], "0."); 
	mpfrts_set_str(&v[2], "0."); 
	mpfrts_set_str(&v[3], "2.3804761428476167"); 

/* ---     INTEGRATION POINTS    --- */
	for(i=0; i<nipt; i++) mpfrts_init(&lt[i]);
	mpfrts_set_str(&lt[0], "0." ); 
	mpfrts_set_str(&lt[1], "62.83185307179586" ); 

/* --- REQUIRED TOLERANCES --- */
	mpfrts_init(&tolrel); 
	mpfrts_init(&tolabs); 
	mpfrts_set_str(&tolrel, "1.e-14");
	mpfrts_set_str(&tolabs, "1.e-14");

/***********************************************************/
/***********************************************************/
/*        OUTPUT:                                           */
/***********************************************************/
/***********************************************************/

/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

	for (i=0; i<nvar; i++) {
			mpfrts_init(&x[i]);
			mpfrts_set(&x[i], v[i]);
		}
	
	mpfrts_set_info_taylor();

	mp_tides(mpfr_kepler, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, NULL, NULL);

  mpfrts_init(&aux);
  mpfrts_init(&error);

	mpfrts_set_str (&aux, "0."); mpfrts_set_str (&error, "-1.");
	for (i=0; i<nvar; i++) {
		mpfrts_sub (&aux, v[i],x[i]);
		mpfrts_abs (&aux, aux);
		/* mpfrts_div (&aux, aux, x[i]);*/
		if (mpfrts_greater (aux, error))
			mpfrts_set (&error, aux);
	}

	if (mpfrts_greater (aux, "1.e-10"))
		return 0;
	else {
    mpfrts_write ("Error", error);
		return 1;
	}

	return 0;
}


