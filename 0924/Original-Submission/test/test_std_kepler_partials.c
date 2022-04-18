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
#include "std_kepler_partials.h"

int main() {

	int nvar = 4;
	int npar = 1;
	int nfun = 0;
	int nipt = 2;
	int i;
	double aux, error;
	double v[nvar], p[npar], lt[nipt];
	double x[25];
	double tolrel, tolabs; 
	double** data;
	FILE   *fd;


/************************************************************/
/************************************************************/
/*     INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 1. ; 

/* --- INITIAL VALUES --- */
	v[0] = 0.30000000000000004 ; 
	v[1] = 0 ; 
	v[2] = 0 ; 
	v[3] = 2.3804761428476167 ; 

/* ---     INTEGRATION POINTS    --- */
	lt[0] = 0 ; 
    lt[1] = 60.0;

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-14 ;
	tolabs = 1.e-14 ;

/***********************************************************/
/***********************************************************/
/*        OUTPUT:         data matrix                     */
/***********************************************************/
/***********************************************************/

	Array2DB_init(&data, nipt, std_kepler_partials_columns());

	fd = fopen ("out_std_kepler_partials.dat", "r");

  for (i=0; i<=std_kepler_partials_columns(); i++) fscanf (fd,"%lf", &x[i]);
  fclose(fd);  

/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

  set_info_taylor();

	dp_tides(std_kepler_partials, nvar, npar, nfun, v, p, 
			lt, nipt, tolrel, tolabs, data, NULL);


	aux = 0.;
  error = -1.;

	for (i=0; i<std_kepler_partials_columns(); i++) {
		aux = fabs(x[i]-data[nipt-1][i])/fabs(data[nipt-1][i]);
    if (aux > error) error = aux;
  }

	if (error < 1.e-08)
		return 0;
	else {
    printf("%25.17e\n",error);
		return 1;
	}
}


