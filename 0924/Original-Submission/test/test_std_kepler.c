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
#include "std_kepler.h"

int main() {

	int nvar = 4;
	int npar = 1;
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
	p[0] = 1. ; 

/* --- INITIAL VALUES --- */
	v[0] = 0.30000000000000000 ; 
	v[1] = 0. ; 
	v[2] = 0. ; 
	v[3] = 2.3804761428476167 ; 

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0. ; 
	lt[1] = 62.83185307179586 ; 

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

	dp_tides(std_kepler, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, NULL, NULL);

  aux = 0.;
  error = -1.;

	for (i=0; i<nvar; i++) {
		aux = fabs(v[i]-x[i]);
    if (aux > error) error = aux;
  }


	if (error < 1.e-10)
		return 0;
	else {
    printf("%25.17e\n",error);
		return 1;
	}
}
