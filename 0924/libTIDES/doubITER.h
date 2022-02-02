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

#ifndef _mpfrITER_H
#define _mpfrITER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "doubNUM.h"


extern	int		MAX_ORDER;
extern	long	NDER;
extern	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
extern	int		*PARTIAL_LIST, *FUNCTION_LIST;
extern	char 	**STR_DER;
extern	long	*POS_DER;
extern	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
extern	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;


		
void	varDB_init(double var[NVARS+1][NDER][MAX_ORDER+1], 
		double v[], double t);

void	parDB_init(double par[NPARS][NDER][MAX_ORDER+1], double p[]);

void	linkDB_init(double lk[NLINKS][NDER][MAX_ORDER+1]);

void	derDB_init(double var[NVARS+1][NDER][MAX_ORDER+1],
		double par[NPARS][NDER][MAX_ORDER+1], double v[]);

void	write_solution_DB(double cvf[][MAX_ORDER+1],
		double var[NVARS+1][NDER][MAX_ORDER+1], 
		double link[NLINKS][NDER][MAX_ORDER+1]);

void	double_htilde(double h[NDER][MAX_ORDER+1], long j, long v, long i, 
		double *ht, int ORDER_INDEX);

void	double_var_t(double f[NDER][MAX_ORDER+1], 
		double u[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_var_t_c(char* cs, double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_var_t_cc(double c, double w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_add_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_add_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_add_t_cc(double c, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_sub_t_cc(double c, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);


void	double_mul_t(double u[NDER][MAX_ORDER+1], double v[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_mul_t_c(char* cs, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_mul_t_cc(double c, double u[NDER][MAX_ORDER+1], 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_div_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_div_t_vc(double u[NDER][MAX_ORDER+1], double c, 
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_div_t_cv(double c, double u[NDER][MAX_ORDER+1],
						double w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_inv_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1],
		int ORDER_INDEX);
void	double_exp_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1],
		int ORDER_INDEX);
void	double_pow_t_c(double u[NDER][MAX_ORDER+1], char* cs, 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_pow_t_cc(double u[NDER][MAX_ORDER+1], double c, 
					   double w[NDER][MAX_ORDER+1], int ORDER_INDEX);

void	double_sct_0(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i);
void	double_sct_i(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	double_sin_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_cos_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    double_sin_cos_t (double f[NDER][MAX_ORDER+1], 
		double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);
void	double_sinh_t(double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_cosh_t(double c[NDER][MAX_ORDER+1], double s[NDER][MAX_ORDER+1], 
		double w[NDER][MAX_ORDER+1], int ORDER_INDEX);
void    double_sinh_cosh_t (double f[NDER][MAX_ORDER+1], 
		double s[NDER][MAX_ORDER+1], double c[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);

void	double_fgt_0(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], long i);
void	double_fgt_i(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], long i, int ORDER_INDEX);
void	double_log_t(double u[NDER][MAX_ORDER+1], double w[NDER][MAX_ORDER+1], 
		int ORDER_INDEX);
void	double_asin_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_acos_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_atan_t(double f[NDER][MAX_ORDER+1], double g[NDER][MAX_ORDER+1], 
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_asinh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_acosh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);
void	double_atanh_t(double f[NDER][MAX_ORDER+1],double g[NDER][MAX_ORDER+1],
		double h[NDER][MAX_ORDER+1], int ORDER_INDEX);


#endif

