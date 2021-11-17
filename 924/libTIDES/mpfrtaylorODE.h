/****************************************************************************
 libTIDES. Version 1.3.0.
 This file is part of TIDES.
 
 Contributors:
 
 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
 *****************************************************************************/

#ifndef _taylorODE_H
#define _taylorODE_H
#ifndef real_MP
#define real_MP
#endif
#include <stdio.h>
#include <stdlib.h>
#include "mpfr.h"
#include <math.h>

#include "mpfrNUMdef.h"
#include "mpfrNUM.h"


typedef long (*LinkedFunction)(mpfr_t t, mpfr_t v[], 
		mpfr_t p[], int orden, mpfr_t cvfd[][orden+1]);

void mpfrts_use_default_step_estimator ();
void mpfrts_set_info_taylor();
void mpfrts_unset_info_taylor();
void mpfrts_str_info_taylor();
void mpfrts_add_info_step(mpfr_t tstep);


int mpfrts_taylor_order(mpfr_t eps);

void mpfrts_norm_inf(mpfr_t *rop, int n, int k, mpfr_t coef[][k+1]);

void mpfrts_compute_step(mpfr_t *rop, mpfr_t tol, int n, int ord, 
		mpfr_t coef[][ord+1]);

void mpfrts_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, int n, int ord, mpfr_t coef[][ord+1]);

void mpfrts_taylor_horner(int n, int ord, mpfr_t coef[][ord+1],
		mpfr_t t, mpfr_t x[]);
void mpfrts_taylor_horner_der(int n, int ord, mpfr_t coef[][ord+1],
						  mpfr_t t, mpfr_t x[]);

void mpfrts_write_taylor_solution( int n, int j, mpfr_t tini, 
	mpfr_t x[],  mpfr_t** mat, FILE* fileout);
	

int mpfrts_valid_step (LinkedFunction fcn, mpfr_t *step, mpfr_t tip, 
		mpfr_t eps, int nvar, int ncol, int order,
		mpfr_t cvfd[][order+1], mpfr_t p[]);

void mp_tides(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	mpfr_t x[], mpfr_t p[],
	mpfr_t lt[], int ntes, 	
	mpfr_t tolrel, mpfr_t tolabs, 
	mpfr_t** mat, FILE* fileout);

void mp_tides_point(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	mpfr_t x[], mpfr_t p[],
	mpfr_t t0, mpfr_t tf, mpfr_t dt, 	
	mpfr_t tolrel, mpfr_t tolabs,   
	mpfr_t** mat, FILE* fileout);

void mp_tides_poaux(LinkedFunction fcn, 
	int nvar, int npar, 
	mpfr_t x[], mpfr_t p[],
	mpfr_t tini, mpfr_t tend,  	
	mpfr_t tolrel, mpfr_t tolabs,   
	mpfr_t *derini, mpfr_t *derend, mpfr_t *partials);

int getOrder ();
int getNsteps ();
#endif
