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
#include "std_lorenz.h"

int main() {

	int nvar = 3;
	int npar = 3;
	int nfun = 0;
	int nipt = 2;
	int i;
	double aux, error;
	double v[nvar], p[npar], lt[nipt];
	double x[nvar];
	double tolrel, tolabs; 
	FILE   *fd;
	extern double   fac1,fac2,fac3;
	extern double   rmaxstep,rminstep; 
	extern int      nitermax, nordex, minord, maxord; 
	extern int      defect_error_control;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 10. ; 
	p[1] = 28. ; 
	p[2] = 2.666666666666667 ; 

/* --- INITIAL VALUES --- */
	v[0] = -13.7636106821342 ; 
	v[1] = -19.5787519424518 ; 
	v[2] = 27. ; 

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0 ; 
	lt[1] = 1.558652210716175 ; 

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-14 ;
	tolabs = 1.e-14 ;

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

	for (i=0; i<nvar; i++) x[i] = v[i];

  set_info_taylor();

	dp_tides(std_lorenz, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, NULL, NULL);

  aux = 0.;
  error = -1.;

	for (i=0; i<nvar; i++) {
		aux = fabs(v[i]-x[i]);
    if (aux > error) error = aux;
  }

	if (error < 1.e-12)
		return 0;
	else
		return 1;
}


