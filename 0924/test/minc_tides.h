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

#ifndef minc_tides_HeadFile 
#define minc_tides_HeadFile 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

double mul_mc(double* u,   double* v,   int k);     
double div_mc(double* u,   double* v,   double* w,int k);   
double exp_mc(double* u,   double* v,   int k);     
double pow_mc_c(double* u,   double e,   double* w,   int k);   
double log_mc(double* u,   double* w,   int k);     
double sin_mc(double* u,   double* v,   double* w,   int k);   
double cos_mc(double* u,   double* v,   double* w,   int k);

void	horner_mc( double *v, double t, int ORDER) ;
void	hornerd_mc(double *v, double t, int ORDER) ; 
void	dense_output_mc(double t0, double step, int ORDER);

double	norm_inf_vec_mv();
double	norm_inf_mat_mc(int ord);
void	declare_matrix_coefs_mc();
void	tolerances_mc(double *tol, double *tolo, int *ORDER);
double	steps_mc(int ORDER, double tol);
void	steps_DEC_mc(double t0, double tol, int ORDER, double *step);


void	mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO); 
void	minc_tides(double *var, int nvar, double *par, int npar,  double tini, double tend, double dt,
                double tol_rel, double tol_abs); 


static inline int min_i(int a, int b) { 
	return a < b ? a: b; 
} 

static  inline int max_i(int a, int b) { 
	return a > b ? a: b; 
} 

static  inline double min_d(double a, double b) { 
	return a < b ? a: b; 
} 

static  inline double max_d(double a, double b) { 
	return a > b ? a: b; 
} 


#endif 



