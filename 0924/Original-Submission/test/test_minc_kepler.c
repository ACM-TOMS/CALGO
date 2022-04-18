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

#include "minc_tides.h"

int main() {

	int  i, VARS, PARS; 
	VARS = 4;
	PARS = 1;
	double aux, error;
	double tolrel, tolabs, tini, tend, dt; 
	double v[VARS], p[PARS]; 
	double x[VARS];
	extern double   fac1,fac2,fac3;
	extern double   rmaxstep,rminstep; 
	extern int      nitermax, nordinc, minord, maxord;
	extern int      dense_output, defect_error_control;
	extern int      accepted_steps, rejected_steps;
	extern FILE     *fd; 



/************************************************************/
/************************************************************/
/*      INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 1.e0 ; 

/* --- INITIAL VALUES --- */
	v[0] = 0.30000000000000000e0 ; 
	v[1] = 0e0 ; 
	v[2] = 0e0 ; 
	v[3] = 2.3804761428476167e0 ; 

/* --- INITIAL INTEGRATION POINT --- */
	tini = 0e0 ;

/* --- ENDPOINT OF INTEGRATION   --- */
	tend = 62.83185307179586e0 ;

/* --- DELTA t FOR DENSE OUTPUT  --- */
	dt   = 62.83185307179586e0 ;

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-14 ;
	tolabs = 1.e-14 ;

/***********************************************************/
/***********************************************************/
/*             DENSE OUTPUT (file, screen or none)          */
/***********************************************************/
/***********************************************************/

	dense_output = 0;

/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/

	for (i=0; i<VARS; i++) x[i] = v[i];

	minc_tides(v,VARS,p,PARS,tini,tend,dt,tolrel,tolabs);

  aux = 0.;
  error = -1.;

	for (i=0; i<VARS; i++) {
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




