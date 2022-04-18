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
#ifndef real_Double
#define real_Double
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "doubNUMdef.h"
#include "doubNUM.h"




typedef long (*LinkedFunction)(double t, double v[], 
		double p[], int orden, double cvfd[][orden+1]);

void use_default_step_estimator ();
void set_info_taylor();
void unset_info_taylor();
void str_info_taylor();
void add_info_step(realNUM tstep);


int  taylor_order(double eps);

void norm_inf(double *rop, int n, int k, double coef[][k+1]);

void compute_step(double *rop, double tol, int n, int ord, double coef[][ord+1]);

void compute_tol (double *tol, double tolrel, double tolabs, int n, int ord, double coef[][ord+1]);

void taylor_horner(int n, int ord, double coef[][ord+1], double t, double x[]);
void taylor_horner_der (int n, int ord, double coef[][ord+1], double t, double x[]);

void write_taylor_solution(int n,  int j, double tini, 
	double x[], double **mat, FILE* fileout);
	
void dp_tides(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	double x[], double p[],
	double lt[], int ntes, 	
	double tolrel, double tolabs,   
	double **mat, FILE* fileout);

void dp_tides_point(LinkedFunction fcn, 
	int nvar, int npar, int nfun, 
	double x[], double p[],
	double t0, double tf, double dt, 	
	double tolrel, double tolabs,   
	double **mat, FILE* fileout);

void dp_tides_poaux(LinkedFunction fcn, 
	int nvar, int npar, 
	double *x, double *p,
	double tini, double tend, 	
	double tolrel, double tolabs,   
	double* derini, double *derfin, double *partials) ;


#endif
